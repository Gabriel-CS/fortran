program PROBLEMA_ABSTRATO
    implicit none
    integer, parameter :: n = 122
    integer :: i, IT, NT
    real(kind=8) :: HT, TT, tempo
    real(kind=8) :: lambda
    real(kind=8) :: Y(n), Ynew(n)
    real(kind=8) :: yanaly, erro
    real(kind=8) :: A(n,n), b(n)

    
    lambda = -100.0d0
    TT = 2.0d0         ! tempo total 
    HT = 0.009d0        ! passo de tempo 
    NT = int(TT/HT)

    ! ==========================
    ! Inicializações
    ! ==========================
    Y = 0.0d0
    A = 0.0d0
    b = 1.0d0

    do i = 1, n
        A(i,i) = lambda
    end do

    ! ==========================
    ! Abrir arquivo de saída
    ! ==========================
    open(unit=10, file="parte2.csv", status="replace")
    write(10,'(A)') "tempo,Y_num,Y_analy,erro"

    ! ==========================
    ! Loop temporal (Euler explícito)
    ! ==========================
    do IT = 0, NT

        tempo = IT * HT

        ! solução analítica
        yanaly = (1.0d0/abs(lambda))*(1.0d0 - exp(lambda*tempo))

        erro = abs(Y(1) - yanaly)

        write(10,'(ES16.8,",",ES16.8,",",ES16.8,",",ES16.8)') tempo, Y(1), yanaly, erro

        !
        if (IT < NT) then
            Ynew = Y + HT*(matmul(A,Y) + b)
            Y = Ynew
        end if

    end do

    close(10)

    print *, "Simulação finalizada."
    print *, "Arquivo gerado: resultados_pa.csv"

end program PROBLEMA_ABSTRATO
