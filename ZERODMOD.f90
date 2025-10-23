module routines_p
    implicit none
    public
    REAL, DIMENSION(:,:), ALLOCATABLE :: A
    REAL, DIMENSION(:), ALLOCATABLE :: b
    real, dimension(:), ALLOCATABLE :: Y0
    real :: t
    integer :: n

contains
    
    subroutine DADOSMOD !Função para leitura de dados de um arquivo
        
    end subroutine DADOSMOD

    
    subroutine FORMAB() !Fução para forma as matrizes A e b

        ALLOCATE(A(n, n))
        ALLOCATE(b(n))
        
        integer :: i, j

        print *, "Leitura da Matriz A:"
        do i = 1, n
            do j = 1, n
                read *, A(i, j)
            end do
        end do
        
        print *, "\nLeitura do vetor b"
        do i = 1, n
            b(i) = REAL(i)
        end do

    end subroutine FORMAB
    
    subroutine MIVEULER()
        real, dimension(n) :: Y, dYdt
        integer :: i

        print *, "\nLeitura do vetor Y"
        do i = 1, n
            Y(i) =  real(i)
        end do

        dy = matmul(A, Y**t) + b

    end subroutine MIVEULER

end module routines_p


program ZERODMOD
    USE routines_p
    implicit none

    call FORMAB(n)
    print *, "Olá mundo!"    
end program ZERODMOD