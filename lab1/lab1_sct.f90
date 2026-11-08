module variants_lib
    implicit none
    
    integer, parameter :: vUNKNOWN = 0 !метод для отслеживания ошибок
    integer, parameter :: vSERIAL = 1 !последовательный алгоритм 
    integer, parameter :: vSPMD = 2 !SPMD
    integer, parameter :: vREDUCTION = 3 !автоматический цикл с редукцией
    integer, parameter :: vATOMIC = 4 !автоматический цикл с атомарными операциями 
    integer, parameter :: vCRITICAL = 5 !автоматический цикл с критической секцией 
    integer, parameter :: vFALSE_SHARING = 6 !вариант с общим массивом

    contains
    
    !превращает константы vSERIAL, vSPMD,... в текстовое слово
    function var2str(var) result(str)
        integer, intent(in) :: var
        character(len=:), allocatable :: str
        select case (var)
        case (vSERIAL); str = "SERIAL"
        case (vSPMD); str = "SPMD"
        case (vREDUCTION); str = "REDUCTION"
        case (vATOMIC); str = "ATOMIC"
        case (vCRITICAL); str = "CRITICAL"
        case (vFALSE_SHARING); str = "FALSE_SHARING"
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
        case ("ATOMIC"); var = vATOMIC
        case ("CRITICAL"); var = vCRITICAL
        case ("FALSE_SHARING"); var = vFALSE_SHARING
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
    real(8), parameter :: a = 0.0d0, b = 1.0d0 !границы интегрирования

    type :: Options !структура настройки, внутри код выбранного метода (variant) и число разбиений интеграла (ndiv)
        integer :: variant = vUNKNOWN
        integer(8) :: ndiv = 0
    end type Options

    type ParseResult !структура результат разбора: флаг действия (action), настройки (options),текстовые поля под описание ошибки (error) и аргумент, который эту ошибку вызвал (argument)
        integer :: action = Action_RUN
        type(Options) :: options
        character(len=:), allocatable :: error
        character(len=:), allocatable :: argument
    end type ParseResult

    type(ParseResult) :: parsed !итог проверки аргументов терминала
    character(len=256) :: prog !строка для хранения имени запущенной программы
    character(len=32) :: I_str, Time_str !для перевода значений интеграла и времени в текст для вывода на экран
    integer(8) :: n !итоговое число разбиений отрезка.
    real(8) :: I, exec_time !итоговое значение интеграла, время выполнения расчета
    
   !запуск всех вариантов вычислений
    interface
        function integrate_sequential(a, b, n) result(res)
            real(8), intent(in) :: a, b
            integer(8), intent(in) :: n
            real(8) :: res
        end function
        function integrate_spmd(a, b, n) result(res)
            real(8), intent(in) :: a, b
            integer(8), intent(in) :: n
            real(8) :: res
        end function
        function integrate_reduction(a, b, n) result(res)
            real(8), intent(in) :: a, b
            integer(8), intent(in) :: n
            real(8) :: res
        end function
        function integrate_atomic(a, b, n) result(res)
            real(8), intent(in) :: a, b
            integer(8), intent(in) :: n
            real(8) :: res
        end function
        function integrate_critical(a, b, n) result(res)
            real(8), intent(in) :: a, b
            integer(8), intent(in) :: n
            real(8) :: res
        end function
        function integrate_false_sharing(a, b, n) result(res)
            real(8), intent(in) :: a, b
            integer(8), intent(in) :: n
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
    I = 0.d0

    exec_time = omp_get_wtime() !запуск секундомера 
    select case (parsed%options%variant) !переключатель кода метода
    case (vSERIAL); I = integrate_sequential(a, b, n)
    case (vSPMD); I = integrate_spmd(a, b, n)
    case (vREDUCTION); I = integrate_reduction(a, b, n)
    case (vATOMIC); I = integrate_atomic(a, b, n)
    case (vCRITICAL); I = integrate_critical(a, b, n)
    case (vFALSE_SHARING); I = integrate_false_sharing(a, b, n)
    case default; continue
    end select
    exec_time = omp_get_wtime() - exec_time !остановка секундомера

!вывод результатов
    write(I_str, '(F20.9)') I 
    write(Time_str, '(F20.9)') exec_time
    write (*, '(A,",",I0,",",I0,",",A,",",A)') & 
        var2str(parsed%options%variant), omp_get_max_threads(),&
        n, trim(adjustl(I_str)), trim(adjustl(Time_str))

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
                  call set_failure(result_, "Invalid number of segment division", token)
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
       character(len=*), intent(in) :: prog_name  ! имя исполняемого файла программы
       character(len=*), intent(in) :: message    ! текст сообщения об ошибке
       character(len=*), intent(in) :: argument   ! некорректный аргумент, вызвавший ошибку

    ! вывод основное сообщение об ошибке 
       write (*, '(2A)', advance='no') "Error: ", message 
    
    ! если проблемный аргумент передан, дописываем его в той же строке
       if (len(argument) > 0) then
           write (*, '(3A)', advance='no') " -> '", argument, "'"
       end if
    
    
       write (*, '(A)') new_line('a')
    
    ! выводим инструкцию по правильному запуску программы и список доступных методов OpenMP
       write (*, '(A)') "Usage: " // prog_name // " <DIRECTIVE NAME> <NUMBER OF DIVISIONS>"
       write (*, '(A)') "Available names:"
       write (*, '(A)') "  SERIAL SPMD REDUCTION ATOMIC CRITICAL FALSE_SHARING"
    end subroutine report_error


end program main

!последовательный (Sequential) обычный цикл на одном потоке без директив OpenMP.
function integrate_sequential(a, b, n) result(res)
    implicit none
    real(8), intent(in) :: a, b
    integer(8), intent(in) :: n
    real(8) :: res, h, s, x
    integer(8) :: i

    h = (b - a) / real(n, 8) 
    s = 0.0d0
    do i = 1, n 
        x = a + (real(i, 8) - 0.5d0) * h 
        s = s + log(x) / (1.0d0 + x) 
    end do
    res = s * h
end function integrate_sequential

!SPMD (Single Program, Multiple Data)
function integrate_spmd(a, b, n) result(res)
    use omp_lib
    implicit none
    real(8), intent(in) :: a, b
    integer(8), intent(in) :: n
    real(8) :: res, h, total_sum, local_sum, x
    integer(8) :: i
    integer :: tid, nt
! вычисление шага интегрирования и сброс общей суммы
    h = (b - a) / real(n, 8)
    total_sum = 0.0d0
! открытие параллельной области: у каждого потока свои копии private-переменных
    !$omp parallel private(tid, nt, i, x, local_sum) shared(a, h, n, total_sum)
! получение номера потока, общего числа потоков и инициализация локальной суммы
    tid = omp_get_thread_num() 
    nt = omp_get_num_threads()
    local_sum = 0.0d0
!каждый поток считает площади только своих прямоугольников с шагом по циклу равным nt
    do i = tid + 1, n, nt 
        x = a + (real(i, 8) - 0.5d0) * h
        local_sum = local_sum + log(x) / (1.0d0 + x)
    end do

    !$omp atomic
    total_sum = total_sum + local_sum
    !$omp end parallel

    res = total_sum * h
end function integrate_spmd

!Reduction (параллельное вычисление интеграла с автоматическим распределением цикла и сборкой суммы через reduction)
function integrate_reduction(a, b, n) result(res)
    use omp_lib
    implicit none
    real(8), intent(in) :: a, b
    integer(8), intent(in) :: n
    real(8) :: res, h, s, x
    integer(8) :: i

    h = (b - a) / real(n, 8)
    s = 0.0d0
! автоматическое деление витков цикла между потоками; reduction(+:s) создает локальную копию s для каждого потока и в конце безопасно объединяет их (суммирует) в общую s
    !$omp parallel do private(x) reduction(+:s)
    do i = 1, n
        x = a + (real(i, 8) - 0.5d0) * h
        s = s + log(x) / (1.0d0 + x)
    end do
    !$omp end parallel do
    
    res = s * h
end function integrate_reduction

!Parallel For/Do + Atomic (параллельное вычисление интеграла с использованием атомарной операции для сборки промежуточных сумм)
function integrate_atomic(a, b, n) result(res)
    use omp_lib
    implicit none
    real(8), intent(in) :: a, b
    integer(8), intent(in) :: n
    real(8) :: res, h, s, p_s, x, val
    integer(8) :: i

    h = (b - a) / real(n, 8)
    s = 0.0d0
! открытие параллельной области с изоляцией локальных переменных потоков
    !$omp parallel private(x, val, p_s) shared(s, a, h, n)
    p_s = 0.0d0
! распределение витков цикла do между потоками 
    !$omp do
 
    do i = 1, n
        x = a + (real(i, 8) - 0.5d0) * h
        val = log(x) / (1.0d0 + x)
        p_s = p_s + val ! накопление суммы в быструю приватную переменную
    end do
    !$omp end do
!атомарное (безопасное) суммирование локальных итогов p_s в общую переменную s
    !$omp atomic
    s = s + p_s
    !$omp end parallel

    res = s * h
end function integrate_atomic

!Parallel For/Do + Critical (параллельное вычисление интеграла с использованием критической секции для финальной сборки данных)
function integrate_critical(a, b, n) result(res)
    use omp_lib
    implicit none
    real(8), intent(in) :: a, b
    integer(8), intent(in) :: n
    real(8) :: res, h, s, p_s, x, val
    integer(8) :: i

    h = (b - a) / real(n, 8)
    s = 0.0d0

    !$omp parallel private(x, val, p_s) shared(s, a, h, n)
    p_s = 0.0d0
    !$omp do
    do i = 1, n
        x = a + (real(i, 8) - 0.5d0) * h
        val = log(x) / (1.0d0 + x)
        p_s = p_s + val
    end do
    !$omp end do
! критическая секция: строго один поток в один момент времени прибавляет свой результат к s
    !$omp critical
    s = s + p_s
    !$omp end critical
    !$omp end parallel

    res = s * h
end function integrate_critical

!False Sharing (Ложное разделение)
function integrate_false_sharing(a, b, n) result(res)
    use omp_lib
    implicit none
    real(8), intent(in) :: a, b
    integer(8), intent(in) :: n
    real(8) :: res, h, x
    integer(8) :: i
    integer :: max_threads, tid
    real(8), allocatable :: partial_sums(:)

    h = (b - a) / real(n, 8)
! выделение массива под промежуточные суммы по количеству доступных потоков
    max_threads = omp_get_max_threads()
    allocate(partial_sums(0:max_threads-1))
    partial_sums = 0.0d0

    !$omp parallel private(x, tid) shared(partial_sums, a, h, n)
    tid = omp_get_thread_num()
 
    !$omp do
    do i = 1, n
        x = a + (real(i, 8) - 0.5d0) * h
        !потоки постоянно пишут в соседние элементы одного общего массива (это заставляет процессоры непрерывно перезаписывать кэш-линии друг друга)
        partial_sums(tid) = partial_sums(tid) + log(x) / (1.0d0 + x)
    end do
    !$omp end do
    !$omp end parallel

    res = sum(partial_sums) * h
    deallocate(partial_sums)
end function integrate_false_sharing
