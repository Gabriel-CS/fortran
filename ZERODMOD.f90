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
        integer :: i, ios
        character(len=100) :: linha
        PI = 3.14d0

        open(unit=10, file='DADOSMOD.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "ERRO: Nao foi possivel abrir o arquivo DADOSMOD.txt"
            stop
        endif

        ! Lê número de elementos
        read(10,*,iostat=ios) n
        if (ios /= 0) then
            print *, "ERRO: Nao foi possivel ler o numero de elementos"
            stop
        endif
        
        print *, "========================================="
        print *, "Numero de elementos: ", n
        print *, "========================================="

        allocate(COMPI(n), RAIOI(n), ESPHI(n))
        allocate(RESISTI(n), INDUCTI(n), COMPLII(n))

        ! Lê dados geométricos (comprimento, raio, espessura)
        print *, "Lendo dados geometricos..."
        do i = 1, n
            read(10,*,iostat=ios) COMPI(i), RAIOI(i), ESPHI(i)
            if (ios /= 0) then
                print *, "ERRO: Falha ao ler linha ", i, " dos dados geometricos"
                stop
            endif
            print '(A,I3,A,F8.3,A,F8.3,A,F8.3)', "  Elemento ", i, &
                  ": COMP=", COMPI(i), " RAIO=", RAIOI(i), " ESPH=", ESPHI(i)
        end do

        ! Lê propriedades do material (Módulo de Young, Permeabilidade magnética, Densidade)
        read(10,*,iostat=ios) MEY, MI, RO
        if (ios /= 0) then
            print *, "ERRO: Nao foi possivel ler as propriedades do material"
            stop
        endif
        
        close(10)

        print *, "========================================="
        print *, "Propriedades do Material:"
        print *, "  MEY:", MEY
        print *, "  MI:", MI
        print *, "  RO:", RO
        print *, "========================================="
        
        call GEOMAT()
    end subroutine DADOSMOD

    subroutine GEOMAT()
        implicit none
        integer :: m
        real(kind=8) :: raio_quad, raio_cub, raio_quarta
        
        print *, "Calculando propriedades dos elementos..."
        print *, ""
        
        do m = 1, n
            raio_quad = RAIOI(m)**2
            raio_cub = RAIOI(m)**3
            raio_quarta = RAIOI(m)**4
            
            ! Resistência
            RESISTI(m) = (8.0d0 * MI * COMPI(m)) / (PI * raio_quarta)
            
            ! Indutância
            INDUCTI(m) = (9.0d0 * RO * COMPI(m)) / (4.0d0 * PI * raio_quad)
            
            ! Compliância (flexibilidade)
            COMPLII(m) = (3.0d0 * PI * COMPI(m) * raio_cub) / (2.0d0 * MEY * ESPHI(m))
            
            print '(A,I3)', "  Elemento ", m
            print '(A,ES14.6)', "    Resistencia:  ", RESISTI(m)
            print '(A,ES14.6)', "    Indutancia:   ", INDUCTI(m)
            print '(A,ES14.6)', "    Compliancia:  ", COMPLII(m)
            print *, ""
        end do

        print *, "========================================="
        print *, "Calculo das propriedades concluido!"
        print *, "========================================="
    end subroutine GEOMAT
 
    subroutine FORMAB()
        implicit none
        integer :: i, j
        
        allocate(A(n,n))
        allocate(b(n))
        
        ! Inicializa matriz e vetor
        A = 0.0d0
        b = 0.0d0

        print *, ""
        print *, "========================================="
        print *, "ENTRADA DA MATRIZ A E VETOR b"
        print *, "========================================="
        print *, "Digite a matriz A (", n, "x", n, "):"
        print *, "Digite linha por linha, separando valores por espaco:"
        
        do i = 1, n
            print *, "Linha", i, ":"
            read *, (A(i,j), j = 1, n)
        end do
        
        print *, ""
        print *, "Matriz A lida:"
        do i = 1, n
            print '(*(ES12.4,1X))', (A(i,j), j = 1, n)
        end do

        print *, ""
        print *, "Digite o vetor b (", n, " elementos):"
        read *, (b(i), i = 1, n)
        
        print *, "Vetor b lido:"
        print '(*(ES12.4,1X))', b
        print *, "========================================="
    end subroutine FORMAB

    subroutine MIVEULER()
        implicit none
        integer :: i, j, passos
        real(kind=8), allocatable :: dYdt(:)

        allocate(Y(n), dYdt(n))

        print *, ""
        print *, "========================================="
        print *, "METODO DE EULER"
        print *, "========================================="
        print *, "Digite o passo de tempo dt:"
        read *, dt
        print *, "Digite o numero de passos:"
        read *, passos
        print *, "Digite o vetor inicial Y0 (", n, " elementos):"
        read *, (Y(i), i = 1, n)

        print *, ""
        print *, "Executando metodo de Euler..."
        print *, "----------------------------------------"
        do j = 1, passos
            dYdt = matmul(A, Y) + b
            Y = Y + dYdt * dt
            print '(A,I5,A,*(ES13.5,1X))', "Passo ", j, ": Y = ", Y
        end do
        print *, "----------------------------------------"
        print *, "Metodo de Euler concluido!"
        print *, "========================================="

        deallocate(dYdt)
    end subroutine MIVEULER
 
end module routines_p

program ZERODMOD
    use routines_p
    implicit none

    print *, ""
    print *, "========================================="
    print *, "  PROGRAMA ZERODMOD - INICIANDO"
    print *, "========================================="
    print *, ""

    call DADOSMOD()   ! Lê dados e calcula as propriedades
    call FORMAB()     ! Monta A e b
    call MIVEULER()   ! Executa o método de Euler

    ! Libera memória
    if (allocated(COMPI)) deallocate(COMPI, RAIOI, ESPHI)
    if (allocated(RESISTI)) deallocate(RESISTI, INDUCTI, COMPLII)
    if (allocated(A)) deallocate(A)
    if (allocated(b)) deallocate(b)
    if (allocated(Y)) deallocate(Y)

    print *, ""
    print *, "========================================="
    print *, "  EXECUCAO CONCLUIDA COM SUCESSO"
    print *, "========================================="
end program ZERODMOD
