module variants_lib
    implicit none

    integer, parameter :: vUNKNOWN    = 0 !метод отслеживания ошибок
    integer, parameter :: vSERIAL     = 1 !последовательный kji 
    integer, parameter :: vSERIAL_IJK = 2 !последовательный ijk 
    integer, parameter :: vSERIAL_JKI = 3 !последовательный jki 
    integer, parameter :: vOMP_OUTER  = 4 !!$omp do на внешнем цикле k + schedule(runtime) 
    integer, parameter :: vOMP_MID    = 5 !!$omp do на среднем цикле j 
    integer, parameter :: vOMP_INNER  = 6 !!$omp do на внутреннем цикле i 
    integer, parameter :: vCOLLAPSE2  = 7 !!$omp do collapse(2) schedule(runtime) 
    integer, parameter :: vCOLLAPSE3  = 8 !!$omp do collapse(3) schedule(runtime) 
    integer, parameter :: vSECTIONS   = 9 !!$omp sections по направлениям x, y и z

contains

    !превращает константы vSERIAL, vSECTIONS,... в текстовое слово
    function var2str(var) result(str)
        integer, intent(in) :: var
        character(len=:), allocatable :: str
        select case (var)
            case (vSERIAL);     str = "SERIAL"
            case (vSERIAL_IJK); str = "serial_ijk"
            case (vSERIAL_JKI); str = "serial_jki"
            case (vOMP_OUTER);  str = "omp_outer"
            case (vOMP_MID);    str = "omp_mid"
            case (vOMP_INNER);  str = "omp_inner"
            case (vCOLLAPSE2);  str = "collapse2"
            case (vCOLLAPSE3);  str = "collapse3"
            case (vSECTIONS);   str = "SECTIONS"
            case default;       str = "UNKNOWN"
        end select
    end function var2str

    !аналогично в обратную сторону
    function str2var(str) result(var)
        character(len=*), intent(in) :: str
        integer :: var
        select case (str)
            case ("SERIAL");     var = vSERIAL
            case ("serial_ijk"); var = vSERIAL_IJK
            case ("serial_jki"); var = vSERIAL_JKI
            case ("omp_outer");  var = vOMP_OUTER
            case ("omp_mid");    var = vOMP_MID
            case ("omp_inner");  var = vOMP_INNER
            case ("collapse2");  var = vCOLLAPSE2
            case ("collapse3");  var = vCOLLAPSE3
            case ("SECTIONS");   var = vSECTIONS
            case default;        var = vUNKNOWN
        end select
    end function str2var

end module variants_lib

!начало основной программы
program main
    use variants_lib
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none

    integer, parameter :: Action_RUN = 0, Action_ERROR = 1 !флаги проверки правильности программы

    type :: Options !структура настройки код выбранного метода (variant) и размер сетки N (ndiv)
        integer :: variant = vUNKNOWN
        integer(8) :: ndiv = 0 !число узлов сетки N по каждому направлению
    end type Options

    type ParseResult !структура результата разбора аргументов командной строки
        integer :: action = Action_RUN
        type(Options) :: options
        character(len=:), allocatable :: error
        character(len=:), allocatable :: argument
    end type ParseResult

    type(ParseResult) :: parsed !итог проверки аргументов терминала
    character(len=256) :: prog  !строка для хранения имени запущенной программы
    character(len=32) :: Sum_str, Time_str
    integer :: n, m, i, j, k, k_mid, unit_id
    real(8) :: total_sum, exec_time

    !динамические трехмерные массивы узлов сетки N x N x N и расстояний (N-1) x (N-1) x (N-1)
    real(8), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(8), allocatable :: dx(:,:,:), dy(:,:,:), dz(:,:,:), dxy(:,:,:), dxz(:,:,:), dyz(:,:,:)

    !интерфейсы всех вариантов вычислений
    interface
        subroutine calc_serial(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_serial_ijk(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_serial_jki(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_omp_outer(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_omp_mid(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_omp_inner(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_collapse2(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_collapse3(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
        subroutine calc_sections(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
            integer, intent(in) :: n, m
            real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
            real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
        end subroutine
    end interface

    !вывод в консоль ошибок разбора аргументов
    call get_command_argument(0, prog)
    parsed = parse_input_arguments()

    if (parsed%action == Action_ERROR) then
        call report_error(trim(prog), parsed%error, parsed%argument)
        stop 1
    end if

    n = int(parsed%options%ndiv)
    if (n < 2) then
        call report_error(trim(prog), "Grid size N must be >= 2", "")
        stop 1
    end if
    m = n - 1
    total_sum = 0.0d0

    allocate(x(n,n,n), y(n,n,n), z(n,n,n))
    allocate(dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m))

    !инициализация исходной декартовой сетки 
    !$omp parallel do collapse(2) private(i, j, k)
    do k = 1, n
        do j = 1, n
            do i = 1, n
                x(i, j, k) = (real(i - 1, 8) / real(m, 8)) ** 1.15d0
                y(i, j, k) = (real(j - 1, 8) / real(m, 8)) ** 1.10d0
                z(i, j, k) = (real(k - 1, 8) / real(m, 8)) ** 1.05d0
            end do
        end do
    end do
    !$omp end parallel do

    exec_time = omp_get_wtime() !запуск секундомера
    select case (parsed%options%variant) !переключатель кода метода
        case (vSERIAL);     call calc_serial(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vSERIAL_IJK); call calc_serial_ijk(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vSERIAL_JKI); call calc_serial_jki(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vOMP_OUTER);  call calc_omp_outer(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vOMP_MID);    call calc_omp_mid(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vOMP_INNER);  call calc_omp_inner(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vCOLLAPSE2);  call calc_collapse2(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vCOLLAPSE3);  call calc_collapse3(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case (vSECTIONS);   call calc_sections(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
        case default; continue
    end select
    exec_time = omp_get_wtime() - exec_time !остановка секундомера

    !контрольная сумма для проверки правильности расчета
    !$omp parallel do collapse(2) reduction(+:total_sum) private(i, j, k)
    do k = 1, m
        do j = 1, m
            do i = 1, m
                total_sum = total_sum + dx(i,j,k) + dy(i,j,k) + dz(i,j,k) &
                                      + dxy(i,j,k) + dxz(i,j,k) + dyz(i,j,k)
            end do
        end do
    end do
    !$omp end parallel do

    !вывод результатов
    write(Sum_str, '(F20.6)') total_sum
    write(Time_str, '(F20.9)') exec_time
    write (*, '(A,",",I0,",",I0,",",A,",",A)') &
        var2str(parsed%options%variant), omp_get_max_threads(), &
        n, trim(adjustl(Sum_str)), trim(adjustl(Time_str))

    !сохранение среза сетки и расстояний в файл для визуализации при N <= 50
    if (n <= 50) then
        k_mid = max(1, m / 2)
        open(newunit=unit_id, file="grid_slice.csv", status="replace")
        write(unit_id, '(A)') "i,j,x,y,dx,dy,dxy"
        do j = 1, m
            do i = 1, m
                write(unit_id, '(I0,",",I0,5(",",F12.6))') i, j, &
                    x(i,j,k_mid), y(i,j,k_mid), dx(i,j,k_mid), dy(i,j,k_mid), dxy(i,j,k_mid)
            end do
        end do
        close(unit_id)
    end if

    deallocate(x, y, z, dx, dy, dz, dxy, dxz, dyz)

contains

!защитный механизм который превращает аргумент терминала в число
    subroutine parse_positive_int(token, val, ok)
        character(len=*), intent(in) :: token
        integer(8), intent(inout) :: val
        logical, intent(inout) :: ok
        integer :: i, digit

        ok = .false.
        val = int(0, 8)
        if (len(token) == 0) return

        do i = 1, len_trim(token)
            if (token(i:i) < '0' .or. token(i:i) > '9') return
            digit = iachar(token(i:i)) - iachar('0')
            val = val * int(10, 8) + int(digit, 8) !сбор из цифр числа
        end do

        ok = (val > int(0, 8))
    end subroutine parse_positive_int

!запись информации об ошибке в res
    subroutine set_failure(res, message, argument)
        type(ParseResult), intent(inout) :: res
        character(len=*), intent(in) :: message
        character(len=*), intent(in) :: argument

        res%action = Action_Error
        res%error = message
        res%argument = argument
    end subroutine set_failure

!разбор аргументов командной строки и возврат структуры с результатом парсинга
    function parse_input_arguments() result(result_)
        type(ParseResult) :: result_
        integer :: i, argc, positional_count, arg_len
        character(len=:), allocatable :: token
        logical :: ok

        result_%action = Action_Run
        result_%error = ""
        result_%argument = ""

        argc = command_argument_count()
        positional_count = 0

        do i = 1, argc
            if (allocated(token)) deallocate(token)
            call get_command_argument(i, length=arg_len)
            allocate(character(len=arg_len) :: token)
            call get_command_argument(i, value=token)

            !первый аргумент имя метода
            if (positional_count == 0) then
                positional_count = positional_count + 1
                result_%options%variant = str2var(token)
                if (result_%options%variant == vUNKNOWN) then
                    call set_failure(result_, "Invalid directive name", token)
                    return
                end if
                cycle

            !второй аргумент размер сетки N
            else if (positional_count == 1 .and. len(token) > 0) then
                positional_count = positional_count + 1
                call parse_positive_int(token, result_%options%ndiv, ok)
                if (.not. ok) then
                    call set_failure(result_, "Invalid grid size", token)
                    return
                end if
                cycle
            end if

            call set_failure(result_, "Unexpected argument", token)
            return
        end do

        if (positional_count < 2) then
            call set_failure(result_, "Missing required arguments. Expected: NAME N", "")
        end if
    end function parse_input_arguments

!вывод сообщения об ошибке, правила использования и список директив
    subroutine report_error(prog_name, message, argument)
        character(len=*), intent(in) :: prog_name
        character(len=*), intent(in) :: message
        character(len=*), intent(in) :: argument

        write (*, '(2A)', advance='no') "Error: ", message
        if (len(argument) > 0) then
            write (*, '(3A)', advance='no') " -> '", argument, "'"
        end if
        write (*, '(A)') new_line('a')
        write (*, '(A)') "Usage: " // prog_name // " <DIRECTIVE NAME> <GRID SIZE>"
        write (*, '(A)') "Available names:"
        write (*, '(A)') "  SERIAL serial_ijk serial_jki omp_outer omp_mid omp_inner collapse2 collapse3 SECTIONS"
    end subroutine report_error

end program main

!последовательный kji (непрерывный доступ к памяти в fortran)
subroutine calc_serial(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    do k = 1, m
        do j = 1, m
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
    end do
end subroutine calc_serial

!последовательный ijk (скачки по памяти и промахи мимо кэша)
subroutine calc_serial_ijk(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    do i = 1, m
        do j = 1, m
            do k = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
    end do
end subroutine calc_serial_ijk

!последовательный jki
subroutine calc_serial_jki(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    do j = 1, m
        do k = 1, m
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
    end do
end subroutine calc_serial_jki

!!$omp do перед внешним циклом k (непрерывный блок памяти обрабатывает каждая нить)
subroutine calc_omp_outer(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    use omp_lib
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    !$omp parallel do schedule(runtime) private(i, j, k, lx, ly, lz)
    do k = 1, m
        do j = 1, m
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
    end do
    !$omp end parallel do
end subroutine calc_omp_outer


!!$omp do перед средним циклом j 
subroutine calc_omp_mid(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    use omp_lib
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    !$omp parallel private(i, j, k, lx, ly, lz)
    do k = 1, m
        !$omp do
        do j = 1, m
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
        !$omp end do
    end do
    !$omp end parallel
end subroutine calc_omp_mid


!!$omp do перед внутренним циклом i 
subroutine calc_omp_inner(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    use omp_lib
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    !$omp parallel private(i, j, k, lx, ly, lz)
    do k = 1, m
        do j = 1, m
            !$omp do
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
            !$omp end do
        end do
    end do
    !$omp end parallel
end subroutine calc_omp_inner
!!$omp do с collapse(2) и schedule(runtime) 
subroutine calc_collapse2(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    use omp_lib
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    !$omp parallel do collapse(2) schedule(runtime) private(i, j, k, lx, ly, lz)
    do k = 1, m
        do j = 1, m
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
    end do
    !$omp end parallel do
end subroutine calc_collapse2


!!$omp do с collapse(3) и schedule(runtime) 
subroutine calc_collapse3(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    use omp_lib
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    !$omp parallel do collapse(3) schedule(runtime) private(i, j, k, lx, ly, lz)
    do k = 1, m
        do j = 1, m
            do i = 1, m
                lx = x(i+1,j,k) - x(i,j,k)
                ly = y(i,j+1,k) - y(i,j,k)
                lz = z(i,j,k+1) - z(i,j,k)
                dx(i,j,k)  = lx
                dy(i,j,k)  = ly
                dz(i,j,k)  = lz
                dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
            end do
        end do
    end do
    !$omp end parallel do
end subroutine calc_collapse3


!!$omp sections три независимые секции по координатным направлениям x, y, z 
subroutine calc_sections(x, y, z, dx, dy, dz, dxy, dxz, dyz, n, m)
    use omp_lib
    implicit none
    integer, intent(in) :: n, m
    real(8), intent(in) :: x(n,n,n), y(n,n,n), z(n,n,n)
    real(8), intent(out) :: dx(m,m,m), dy(m,m,m), dz(m,m,m), dxy(m,m,m), dxz(m,m,m), dyz(m,m,m)
    integer :: i, j, k
    real(8) :: lx, ly, lz

    !$omp parallel sections private(i, j, k, lx, ly, lz)
        !$omp section
        do k = 1, m
            do j = 1, m
                do i = 1, m
                    lx = x(i+1,j,k) - x(i,j,k)
                    ly = y(i,j+1,k) - y(i,j,k)
                    dx(i,j,k)  = lx
                    dxy(i,j,k) = 0.5d0 * sqrt(lx*lx + ly*ly)
                end do
            end do
        end do

        !$omp section
        do k = 1, m
            do j = 1, m
                do i = 1, m
                    ly = y(i,j+1,k) - y(i,j,k)
                    lz = z(i,j,k+1) - z(i,j,k)
                    dy(i,j,k)  = ly
                    dyz(i,j,k) = 0.5d0 * sqrt(ly*ly + lz*lz)
                end do
            end do
        end do

        !$omp section
        do k = 1, m
            do j = 1, m
                do i = 1, m
                    lx = x(i+1,j,k) - x(i,j,k)
                    lz = z(i,j,k+1) - z(i,j,k)
                    dz(i,j,k)  = lz
                    dxz(i,j,k) = 0.5d0 * sqrt(lx*lx + lz*lz)
                end do
            end do
        end do
    !$omp end parallel sections
end subroutine calc_sections

