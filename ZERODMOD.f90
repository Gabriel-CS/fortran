module routines_p
    implicit none
    public
    integer :: n
    real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI
    real(kind=8), dimension(:), allocatable :: RESISTI, INDUCTI, COMPLII
    real(kind=8), dimension(:,:), allocatable :: A
    real(kind=8), dimension(:), allocatable :: b, Y
    real(kind=8) :: MEY, MI, RO, PI, dt
contains


    subroutine DADOSMOD()
        implicit none
        character(len=50) :: filename
        integer :: i, ios
        PI = 3.14d0

        print *, "Digite o nome do arquivo de dados (ex: dados.txt):"
        read(*,'(A)') filename

        open(unit=10, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Erro ao abrir o arquivo ", trim(filename)
            stop
        endif

        ! Lê número de elementos
        read(10,*) n

        allocate(COMPI(n), RAIOI(n), ESPHI(n))
        allocate(RESISTI(n), INDUCTI(n), COMPLII(n))

        ! Lê dados geométricos
        do i = 1, n
            read(10,*,iostat=ios) COMPI(i), RAIOI(i), ESPHI(i)
            if (ios /= 0) then
                print *, "Erro ao ler linha ", i
                stop
            endif
        end do

        ! Lê propriedades do material
        read(10,*) MEY, MI, RO
        close(10)

        call GEOMAT()
    end subroutine DADOSMOD


    subroutine GEOMAT()
        implicit none
        integer :: m
        do m = 1, n
            RESISTI(m) = (8d0 * MI * COMPI(m)) / (PI * (RAIOI(m)**4))
            INDUCTI(m) = (9d0 * RO * COMPI(m)) / (4d0 * PI * (RAIOI(m)**2))
            COMPLII(m) = (3d0 * PI * COMPI(m) * (RAIOI(m)**3)) / (2d0 * MEY * ESPHI(m))
        end do

        print *, "Cálculo das propriedades geométricas concluído."
    end subroutine GEOMAT
 

    subroutine FORMAB()
        implicit none
        integer :: i, j
        allocate(A(n,n))
        allocate(b(n))

        print *, "Leitura da matriz A (", n, "x", n, "):"
        do i = 1, n
            do j = 1, n
                read *, A(i,j)
            end do
        end do

        print *, "Leitura do vetor b:"
        do i = 1, n
            read *, b(i)
        end do
    end subroutine FORMAB


    subroutine MIVEULER()
        implicit none
        integer :: i, j, passos
        real(kind=8), allocatable :: dYdt(:)

        allocate(Y(n), dYdt(n))

        print *, "Digite o passo de tempo dt:"
        read *, dt
        print *, "Digite o número de passos:"
        read *, passos
        print *, "Digite o vetor inicial Y0:"
        do i = 1, n
            read *, Y(i)
        end do

        do j = 1, passos
            dYdt = matmul(A, Y) + b
            Y = Y + dYdt * dt
            print '(A,I3,A,*(F10.5,1X))', "Passo ", j, ": Y = ", Y
        end do

        deallocate(dYdt)
    end subroutine MIVEULER
 

end module routines_p


program ZERODMOD
    use routines_p
    implicit none

    call DADOSMOD()   ! Lê dados e calcula as propriedades
    call FORMAB()     ! Monta A e b
    call MIVEULER()   ! Executa o método de Euler

    print *, "Execução concluída com sucesso."
end program ZERODMOD
