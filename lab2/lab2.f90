module variants_lib
    implicit none
    
    integer, parameter :: vUNKNOWN = 0 !метод для отслеживания ошибок
    integer, parameter :: vSERIAL = 1 !последовательный алгоритм 
    integer, parameter :: vSPMD = 2 !SPMD
    integer, parameter :: vREDUCTION = 3 !автоматический цикл с редукцией
    
    contains
    
    !превращает константы vSERIAL, vSPMD,... в текстовое слово
    function var2str(var) result(str)
        integer, intent(in) :: var
        character(len=:), allocatable :: str
        select case (var)
        case (vSERIAL); str = "SERIAL"
        case (vSPMD); str = "SPMD"
        case (vREDUCTION); str = "REDUCTION"
        case default; str = "UNKNOWN"
        end select
    end function var2str
   !аналогично в обратную сторону
    function str2var(str) result(var)
        character(len=*), intent(in) :: str
        integer :: var
        select case (str)
        case ("SERIAL"); var = vSERIAL
        case ("SPMD"); var = vSPMD
        case ("REDUCTION"); var = vREDUCTION
        case default; var = vUNKNOWN
        end select
    end function str2var

end module variants_lib

!начало основной программы
program main
    use variants_lib 
    use omp_lib, only: omp_get_wtime, omp_get_max_threads 
    implicit none
    
    integer, parameter :: Action_RUN = 0, Action_ERROR = 1 !флаги проверки правильности программы
  
    type :: Options !структура настройки, внутри код выбранного метода (variant) и число разбиений интеграла (ndiv)
        integer :: variant = vUNKNOWN
        integer(8) :: ndiv = 0 !размер массива N
    end type Options

    type ParseResult !структура результат разбора: флаг действия (action), настройки (options),текстовые поля под описание ошибки (error) и аргумент, который эту ошибку вызвал (argument)
        integer :: action = Action_RUN
        type(Options) :: options
        character(len=:), allocatable :: error
        character(len=:), allocatable :: argument
    end type ParseResult

    type(ParseResult) :: parsed !итог проверки аргументов терминала
    character(len=256) :: prog !строка для хранения имени запущенной программы
    character(len=32) :: Sum_str, Time_str 
    integer(8) :: n, i
    real(8) :: total_sum, exec_time 
    
    ! динамический массив (вектор) 
    real(8), allocatable :: vector(:)
    
   !запуск всех вариантов вычислений
    interface
        function sum_sequential(vec, n) result(res)
            integer(8), intent(in) :: n
            real(8), intent(in) :: vec(n)
            real(8) :: res
        end function
        function sum_spmd(vec, n) result(res)
            integer(8), intent(in) :: n
            real(8), intent(in) :: vec(n)
            real(8) :: res
        end function
        function sum_reduction(vec, n) result(res)
            integer(8), intent(in) :: n
            real(8), intent(in) :: vec(n)
            real(8) :: res
        end function
    end interface

!вывод в консоль ошибок 
    call get_command_argument(0, prog) 
    parsed = parse_input_arguments() 
    
    if (parsed%action == Action_ERROR) then
        call report_error(trim(prog), parsed%error, parsed%argument)
        stop 1
    end if

    n = parsed%options%ndiv
    total_sum = 0.d0
    
    allocate(vector(1:n))
    
   !заполнение массива тестовыми данными (вынесено из замера времени)
    !$omp parallel do
    do i = 1, n
        vector(i) = 1.0d0 / real(i, 8)
    end do

    exec_time = omp_get_wtime() !запуск секундомера 
    select case (parsed%options%variant) !переключатель кода метода
    case (vSERIAL); total_sum = sum_sequential(vector, n)
    case (vSPMD); total_sum = sum_spmd(vector, n)
    case (vREDUCTION); total_sum = sum_reduction(vector, n)
    case default; continue
    end select
    exec_time = omp_get_wtime() - exec_time !остановка секундомера

!вывод результатов
    write(Sum_str, '(F20.9)') total_sum 
    write(Time_str, '(F20.9)') exec_time
    write (*, '(A,",",I0,",",I0,",",A,",",A)') & 
        var2str(parsed%options%variant), omp_get_max_threads(),&
        n, trim(adjustl(Sum_str)), trim(adjustl(Time_str))
    
    deallocate(vector)
    
    contains
   !защитный механизм-превращает аргумент терминала в число 
    subroutine parse_positive_int(token, val, ok) 
        character(len=*), intent(in) :: token
        integer(8), intent(inout) :: val 
        logical, intent(inout) :: ok 
        integer :: i, digit

        ok = .false.
        val = int(0, 8)
        if (len(token) == 0) return
        do i = 1, len_trim(token)
            if (token(i:i) < '0' .or. token(i:i) > '9') return; 
            digit = iachar(token(i:i)) - iachar('0') 

            val = val * int(10, 8) + int(digit, 8) !сбор из цифр числа 
        end do
        ok = (val > int(0, 8))
    end subroutine parse_positive_int

!запись информации об ошибке в структуру res (флаг Action_Error, сообщение message и проблемное слово argument)
    subroutine set_failure(res, message, argument)
       type(ParseResult), intent(inout) :: res       
       character(len=*), intent(in)     :: message   
       character(len=*), intent(in)     :: argument  

       res%action = Action_Error
       res%error = message
       res%argument = argument
    end subroutine set_failure

    ! разбор аргументов командной строки и возвращает структуру с результатом парсинга
    function parse_input_arguments() result(result_)
       type(ParseResult) :: result_
       integer :: i, argc, positional_count, arg_len
       character(len=:), allocatable :: token
       logical :: ok

    ! инициализация структуры: по умолчанию ошибок нет
       result_%action = Action_Run
       result_%error = ""
       result_%argument = ""
    
    ! получение общего количества переданных аргументов
       argc = command_argument_count()
       positional_count = 0 

    ! перебор всех аргументов командной строки в цикле
       do i = 1, argc
        ! динамическое выделение памяти под размер текущего аргумента
           if (allocated(token)) deallocate(token)
           call get_command_argument(i, length=arg_len)
           allocate(character(len=arg_len) :: token)
           call get_command_argument(i, value=token)

        ! первый аргумент - должен быть именем метода 
           if (positional_count == 0) then
               positional_count = positional_count + 1
               result_%options%variant = str2var(token)
               if (result_%options%variant == vUNKNOWN) then
                   call set_failure(result_, "Invalid directive name", token)
                   return
               end if
               cycle
            
        ! второй аргумент - должен быть числом разбиений 
           else if (positional_count == 1 .and. len(token) > 0) then
               positional_count = positional_count + 1
               call parse_positive_int(token, result_%options%ndiv, ok)
               if (.not. ok) then
                  call set_failure(result_, "Invalid array size", token)
                  return
               end if
               cycle
           end if

        ! все последующие аргументы - ошибка
           call set_failure(result_, "Unexpacted argument", token)
           return
       end do

    ! если введено меньше двух обязательных аргументов- ошибка
       if (positional_count < 2) then
           call set_failure(result_, "Missing required arguments. Expected: NAME N", "")
       end if
    end function parse_input_arguments


    ! выводит сообщение об ошибке, правила использования программы и список доступных директив
    subroutine report_error(prog_name, message, argument)
       character(len=*), intent(in) :: prog_name  
       character(len=*), intent(in) :: message    
       character(len=*), intent(in) :: argument   
       write (*, '(2A)', advance='no') "Error: ", message 
       if (len(argument) > 0) then
           write (*, '(3A)', advance='no') " -> '", argument, "'"
       end if
       write (*, '(A)') new_line('a')
       write (*, '(A)') "Usage: " // prog_name // " <DIRECTIVE NAME> <ARRAY SIZE>"
       write (*, '(A)') "Available names:"
       write (*, '(A)') "  SERIAL SPMD REDUCTION"
    end subroutine report_error

end program main

!последовательный (Sequential) обычный цикл на одном потоке без директив OpenMP.
function sum_sequential(vec, n) result(res)
    implicit none
    integer(8), intent(in) :: n
    real(8), intent(in) :: vec(n)
    real(8) :: res, s
    integer(8) :: i

    s = 0.0d0
    do i = 1, n 
        s = s + vec(i)
    end do
    res = s
end function sum_sequential

!SPMD (ручное циклическое распределение)
function sum_spmd(vec, n) result(res)
    use omp_lib
    implicit none
    integer(8), intent(in) :: n
    real(8), intent(in) :: vec(n)
    real(8) :: res, total_sum, local_sum
    integer(8) :: i
    integer :: tid, nt

    total_sum = 0.0d0
    !$omp parallel private(tid, nt, i, local_sum) shared(vec, n, total_sum)
    tid = omp_get_thread_num() 
    nt = omp_get_num_threads()
    local_sum = 0.0d0
    
    ! распределение витков цикла вручную с шагом nt
    do i = tid + 1, n, nt 
        local_sum = local_sum + vec(i)
    end do

    ! безопасное (атомик) сложение локальных результатов нитей в общую переменную
    !$omp atomic
    total_sum = total_sum + local_sum
    !$omp end parallel

    res = total_sum
end function sum_spmd

!Reduction (вариант с директивой !$omp do (автоматическое распределение))
function sum_reduction(vec, n) result(res)
    use omp_lib
    implicit none
    integer(8), intent(in) :: n
    real(8), intent(in) :: vec(n)
    real(8) :: res, s
    integer(8) :: i

    s = 0.0d0
    ! автоматическая нарезка цикла do между запущенными нитями
    ! reduction(+:s) изолирует промежуточные переменные и суммирует их в конце
    !$omp parallel do reduction(+:s)
    do i = 1, n
        s = s + vec(i)
    end do
    !$omp end parallel do
    
    res = s
end function sum_reduction
