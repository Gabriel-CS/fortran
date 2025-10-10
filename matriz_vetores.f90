program matriz_vetores
    implicit none
    integer :: n, i, j
    real, allocatable :: M(:,:), v1(:), v2(:), resultado(:)

    ! Leitura do tamanho da matriz e dos vetores
    print *, "Digite o tamanho n da matriz (n x n) e dos vetores:"
    read *, n

    ! Alocação das variáveis
    allocate(M(n,n), v1(n), v2(n), resultado(n))

    ! Leitura da matriz M
    print *, "Digite os elementos da matriz M (linha por linha):"
    do i = 1, n
        read *, (M(i,j), j = 1, n)
    end do

    ! Leitura do vetor v1
    print *, "Digite os elementos do vetor v1:"
    read *, (v1(i), i = 1, n)

    ! Leitura do vetor v2
    print *, "Digite os elementos do vetor v2:"
    read *, (v2(i), i = 1, n)

    ! Cálculo de resultado = M * v1 + v2
    resultado = matmul(M, v1) + v2

    ! Impressão do resultado
    print *, "Resultado (M * v1 + v2):"
    do i = 1, n
        print *, resultado(i)
    end do

    ! Liberação da memória
    deallocate(M, v1, v2, resultado)

end program matriz_vetores