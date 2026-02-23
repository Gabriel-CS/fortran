module routines_p
    implicit none
    public :: DADOSMOD, GEOMAT, PULSOPQ, FORMAB, MIVEULER
    private

contains

    subroutine DADOSMOD(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO, PI
        real(kind=8), dimension(:), allocatable :: COMPI, RAIOI, ESPHI, MEY
        integer :: i, ios
        real(kind=8) :: t_MI, t_RO, t_PI

        open(unit=10, file='DADOSMOD.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop "ERRO: DADOSMOD.txt nao encontrado."

        read(10,*) n
        allocate(COMPI(n), RAIOI(n), ESPHI(n), MEY(n))
        do i = 1, n
            read(10,*) COMPI(i), RAIOI(i), ESPHI(i), MEY(i), t_MI, t_RO, t_PI
            if (i == 1) then
                MI = t_MI
                RO = t_RO
                PI = t_PI
            endif
        end do
        close(10)
    end subroutine DADOSMOD

    subroutine GEOMAT(n, COMPI, RAIOI, ESPHI, MEY, MI, RO, PI, RESISTI, INDUCTI, COMPLII)
        implicit none
        integer :: n
        real(kind=8) :: MI, RO, PI
        real(kind=8) :: COMPI(n), RAIOI(n), ESPHI(n), MEY(n)
        real(kind=8) :: RESISTI(n), INDUCTI(n), COMPLII(n)
        real(kind=8) :: raio_quad, raio_cub, raio_quarta
        integer :: m, ios

        do m = 1, n
            raio_quad = RAIOI(m)**2
            raio_cub = RAIOI(m)**3
            raio_quarta = RAIOI(m)**4

            RESISTI(m) = (8.0d0 * MI * COMPI(m)) / (PI * raio_quarta)
            INDUCTI(m) = (9.0d0 * RO * COMPI(m)) / (4.0d0 * PI * raio_quad)
            COMPLII(m) = (3.0d0 * PI * COMPI(m) * raio_cub) / (2.0d0 * MEY(m) * ESPHI(m))
        end do

        open(unit=20, file='GEOMAT.csv', status='replace', action='write', iostat=ios)
        write(20,'(A)') 'Trecho,RESISTI,INDUCTI,COMPLII'
        
        ! CORRECTED FORMAT: Using explicit commas in the format string
        do m = 1, n
            write(20,'(I0,",",ES15.6,",",ES15.6,",",ES15.6)') m, RESISTI(m), INDUCTI(m), COMPLII(m)
        end do
        
        close(20)
    end subroutine GEOMAT

    subroutine PULSOPQ(IT, HT, PI, PIN, QIN)
        implicit none
        integer :: IT
        real(kind=8) :: HT, PI
        real(kind=8) :: PIN, QIN
        real(kind=8) :: TS, TCC, PS, PD, QMAX, TI, A, B

        TS   = 0.3d0
        TCC  = 1.0d0
        PS   = 120.0d0
        PD   = 80.0d0
        QMAX = 300.0d0

        ! CORREÇÃO: TI é o tempo DENTRO do ciclo atual (reseta a cada 1 segundo)
        TI = mod(dble(IT) * HT, TCC)

        if (TI <= TS) then
            PIN = (PD + ((PS - PD) / TS) * TI) * 1333.22d0
            QIN = QMAX * (sin(PI * TI / TS))**2
        else
            A   = (PS - PD) / (TS - TCC)
            B   = (PD * TS - PS * TCC) / (TS - TCC)
            PIN = (A * TI + B) * 1333.22d0
            QIN = 0.0d0
        end if

    end subroutine PULSOPQ

    subroutine FORMAB(n, RESISTI, INDUCTI, COMPLII, HT, PI, IT, A, b)
        implicit none
        integer:: n, IT
        real(kind=8) :: HT, PI
        real(kind=8) :: RESISTI(n), INDUCTI(n), COMPLII(n)
        real(kind=8) :: A(122,122)
        real(kind=8) :: b(122)
        real(kind=8) :: PIN, QIN

        A = 0.0d0
        b = 0.0d0

                ! --- LINHA 1 ---
        A(1,4)  =  1.0d0/COMPLII(1)
        A(1,10) =  1.0d0/COMPLII(1)
        A(1,16) =  1.0d0/COMPLII(1)

        ! --- LINHA 2 ---
        A(2,1)  =  1.0d0/INDUCTI(1)
        A(2,4)  = -RESISTI(1)/INDUCTI(1)
        A(2,10) = -RESISTI(1)/INDUCTI(1)
        A(2,16) = -RESISTI(1)/INDUCTI(1)

        ! --- LINHA 3 --- (bifurcacao: Q_3=col6, Q_4=col8)
        A(3,4) =  1.0d0/COMPLII(2)
        A(3,6) = -1.0d0/COMPLII(2)
        A(3,8) = -1.0d0/COMPLII(2)

        ! --- LINHA 4 ---
        A(4,1) =  1.0d0/INDUCTI(2)
        A(4,3) = -1.0d0/INDUCTI(2)
        A(4,4) = -RESISTI(2)/INDUCTI(2)

        ! --- LINHA 5 --- (terminal)
        A(5,6) = 1.0d0/COMPLII(3)

        ! --- LINHA 6 --- (terminal: upstream P=col3)
        A(6,3) =  1.0d0/INDUCTI(3)
        A(6,5) = -1.0d0/INDUCTI(3)
        A(6,6) = -RESISTI(3)/INDUCTI(3)

        ! --- LINHA 7 --- (terminal)
        A(7,8) = 1.0d0/COMPLII(4)

        ! --- LINHA 8 --- (terminal: upstream P=col3)
        A(8,3) =  1.0d0/INDUCTI(4)
        A(8,7) = -1.0d0/INDUCTI(4)
        A(8,8) = -RESISTI(4)/INDUCTI(4)

        ! --- LINHA 9 --- (bifurcacao: Q_6=col12, Q_7=col14)
        A(9,10) =  1.0d0/COMPLII(5)
        A(9,12) = -1.0d0/COMPLII(5)
        A(9,14) = -1.0d0/COMPLII(5)

        ! --- LINHA 10 ---
        A(10,1)  = -1.0d0/INDUCTI(5)
        A(10,9)  = -1.0d0/INDUCTI(5)
        A(10,10) = -RESISTI(5)/INDUCTI(5)

        ! --- LINHA 11 --- (terminal)
        A(11,12) = 1.0d0/COMPLII(6)

        ! --- LINHA 12 --- (terminal: upstream P=col9)
        A(12,9)  =  1.0d0/INDUCTI(6)
        A(12,11) = -1.0d0/INDUCTI(6)
        A(12,12) = -RESISTI(6)/INDUCTI(6)

        ! --- LINHA 13 --- (terminal)
        A(13,14) = 1.0d0/COMPLII(7)

        ! --- LINHA 14 --- (terminal: upstream P=col9)
        A(14,9)  =  1.0d0/INDUCTI(7)
        A(14,13) = -1.0d0/INDUCTI(7)
        A(14,14) = -RESISTI(7)/INDUCTI(7)

        ! --- LINHA 15 ---
        A(15,16) =  1.0d0/COMPLII(8)
        A(15,18) = -1.0d0/COMPLII(8)
        A(15,20) = -1.0d0/COMPLII(8)

        ! --- LINHA 16 ---
        A(16,1)  =  1.0d0/INDUCTI(8)
        A(16,15) = -1.0d0/INDUCTI(8)
        A(16,18) = -RESISTI(8)/INDUCTI(8)
        A(16,20) = -RESISTI(8)/INDUCTI(8)

        ! --- LINHA 17 ---
        A(17,18) =  1.0d0/COMPLII(9)
        A(17,42) = -1.0d0/COMPLII(9)
        A(17,44) = -1.0d0/COMPLII(9)

        ! --- LINHA 18 ---
        A(18,15) =  1.0d0/INDUCTI(9)
        A(18,17) = -1.0d0/INDUCTI(9)
        A(18,42) = -RESISTI(9)/INDUCTI(9)
        A(18,44) = -RESISTI(9)/INDUCTI(9)

        ! --- LINHA 19 ---
        A(19,20) =  1.0d0/COMPLII(10)
        A(19,22) = -1.0d0/COMPLII(10)
        A(19,24) = -1.0d0/COMPLII(10)

        ! --- LINHA 20 ---
        A(20,15) =  1.0d0/INDUCTI(10)
        A(20,19) = -1.0d0/INDUCTI(10)
        A(20,22) = -RESISTI(10)/INDUCTI(10)
        A(20,24) = -RESISTI(10)/INDUCTI(10)

        ! --- LINHA 21 --- (bifurcacao: Q_13=col26, Q_14=col28)
        A(21,22) =  1.0d0/COMPLII(11)
        A(21,26) = -1.0d0/COMPLII(11)
        A(21,28) = -1.0d0/COMPLII(11)

        ! --- LINHA 22 ---
        A(22,19) =  1.0d0/INDUCTI(11)
        A(22,21) = -1.0d0/INDUCTI(11)
        A(22,26) = -RESISTI(11)/INDUCTI(11)
        A(22,28) = -RESISTI(11)/INDUCTI(11)

        ! --- LINHA 23 --- (bifurcacao: Q_19=col38, Q_20=col40)
        A(23,24) =  1.0d0/COMPLII(12)
        A(23,38) = -1.0d0/COMPLII(12)
        A(23,40) = -1.0d0/COMPLII(12)

        ! --- LINHA 24 ---
        A(24,19) =  1.0d0/INDUCTI(12)
        A(24,23) = -1.0d0/INDUCTI(12)
        A(24,38) = -RESISTI(12)/INDUCTI(12)
        A(24,40) = -RESISTI(12)/INDUCTI(12)

        ! --- LINHA 25 --- (terminal)
        A(25,26) = 1.0d0/COMPLII(13)

        ! --- LINHA 26 --- (terminal: upstream P=col21)
        A(26,21) =  1.0d0/INDUCTI(13)
        A(26,25) = -1.0d0/INDUCTI(13)
        A(26,26) = -RESISTI(13)/INDUCTI(13)

        ! --- LINHA 27 --- (bifurcacao: Q_15=col30, Q_16=col32)
        A(27,28) =  1.0d0/COMPLII(14)
        A(27,30) = -1.0d0/COMPLII(14)
        A(27,32) = -1.0d0/COMPLII(14)

        ! --- LINHA 28 ---
        A(28,21) =  1.0d0/INDUCTI(14)
        A(28,27) = -1.0d0/INDUCTI(14)
        A(28,30) = -RESISTI(14)/INDUCTI(14)
        A(28,32) = -RESISTI(14)/INDUCTI(14)

        ! --- LINHA 29 --- (terminal)
        A(29,30) = 1.0d0/COMPLII(15)

        ! --- LINHA 30 --- (terminal: upstream P=col27)
        A(30,27) =  1.0d0/INDUCTI(15)
        A(30,29) = -1.0d0/INDUCTI(15)
        A(30,30) = -RESISTI(15)/INDUCTI(15)

        ! --- LINHA 31 --- (bifurcacao: Q_17=col34, Q_18=col36)
        A(31,32) =  1.0d0/COMPLII(16)
        A(31,34) = -1.0d0/COMPLII(16)
        A(31,36) = -1.0d0/COMPLII(16)

        ! --- LINHA 32 ---
        A(32,27) =  1.0d0/INDUCTI(16)
        A(32,31) = -1.0d0/INDUCTI(16)
        A(32,34) = -RESISTI(16)/INDUCTI(16)
        A(32,36) = -RESISTI(16)/INDUCTI(16)

        ! --- LINHA 33 --- (terminal)
        A(33,34) = 1.0d0/COMPLII(17)

        ! --- LINHA 34 --- (terminal: upstream P=col31)
        A(34,31) =  1.0d0/INDUCTI(17)
        A(34,33) = -1.0d0/INDUCTI(17)
        A(34,34) = -RESISTI(17)/INDUCTI(17)

        ! --- LINHA 35 --- (terminal)
        A(35,36) = 1.0d0/COMPLII(18)

        ! --- LINHA 36 --- (terminal: upstream P=col31)
        A(36,31) =  1.0d0/INDUCTI(18)
        A(36,35) = -1.0d0/INDUCTI(18)
        A(36,36) = -RESISTI(18)/INDUCTI(18)

        ! --- LINHA 37 --- (terminal)
        A(37,38) = 1.0d0/COMPLII(19)

        ! --- LINHA 38 --- (terminal: upstream P=col23)
        A(38,23) =  1.0d0/INDUCTI(19)
        A(38,37) = -1.0d0/INDUCTI(19)
        A(38,38) = -RESISTI(19)/INDUCTI(19)

        ! --- LINHA 39 --- (terminal)
        A(39,40) = 1.0d0/COMPLII(20)

        ! --- LINHA 40 --- (terminal: upstream P=col23)
        A(40,23) =  1.0d0/INDUCTI(20)
        A(40,39) = -1.0d0/INDUCTI(20)
        A(40,40) = -RESISTI(20)/INDUCTI(20)

        ! --- LINHA 41 ---
        A(41,42) =  1.0d0/COMPLII(21)
        A(41,50) = -1.0d0/COMPLII(21)
        A(41,52) = -1.0d0/COMPLII(21)

        ! --- LINHA 42 ---
        A(42,17) =  1.0d0/INDUCTI(21)
        A(42,41) = -1.0d0/INDUCTI(21)
        A(42,50) = -RESISTI(21)/INDUCTI(21)
        A(42,52) = -RESISTI(21)/INDUCTI(21)

        ! --- LINHA 43 --- (bifurcacao: Q_23=col46, Q_24=col48)
        A(43,44) =  1.0d0/COMPLII(22)
        A(43,46) = -1.0d0/COMPLII(22)
        A(43,48) = -1.0d0/COMPLII(22)

        ! --- LINHA 44 ---
        A(44,17) =  1.0d0/INDUCTI(22)
        A(44,43) = -1.0d0/INDUCTI(22)
        A(44,46) = -RESISTI(22)/INDUCTI(22)
        A(44,48) = -RESISTI(22)/INDUCTI(22)

        ! --- LINHA 45 --- (terminal)
        A(45,46) = 1.0d0/COMPLII(23)

        ! --- LINHA 46 --- (terminal: upstream P=col43)
        A(46,43) =  1.0d0/INDUCTI(23)
        A(46,45) = -1.0d0/INDUCTI(23)
        A(46,46) = -RESISTI(23)/INDUCTI(23)

        ! --- LINHA 47 --- (terminal)
        A(47,48) = 1.0d0/COMPLII(24)

        ! --- LINHA 48 --- (terminal: upstream P=col43)
        A(48,43) =  1.0d0/INDUCTI(24)
        A(48,47) = -1.0d0/INDUCTI(24)
        A(48,48) = -RESISTI(24)/INDUCTI(24)

        ! --- LINHA 49 ---
        A(49,50) =  1.0d0/COMPLII(25)
        A(49,66) = -1.0d0/COMPLII(25)
        A(49,68) = -1.0d0/COMPLII(25)

        ! --- LINHA 50 ---
        A(50,41) =  1.0d0/INDUCTI(25)
        A(50,49) = -1.0d0/INDUCTI(25)
        A(50,66) = -RESISTI(25)/INDUCTI(25)
        A(50,68) = -RESISTI(25)/INDUCTI(25)

        ! --- LINHA 51 --- (bifurcacao: Q_27=col54, Q_28=col56)
        A(51,52) =  1.0d0/COMPLII(26)
        A(51,54) = -1.0d0/COMPLII(26)
        A(51,56) = -1.0d0/COMPLII(26)

        ! --- LINHA 52 ---
        A(52,41) =  1.0d0/INDUCTI(26)
        A(52,51) = -1.0d0/INDUCTI(26)
        A(52,54) = -RESISTI(26)/INDUCTI(26)
        A(52,56) = -RESISTI(26)/INDUCTI(26)

        ! --- LINHA 53 --- (terminal)
        A(53,54) = 1.0d0/COMPLII(27)

        ! --- LINHA 54 --- (terminal: upstream P=col51)
        A(54,51) =  1.0d0/INDUCTI(27)
        A(54,53) = -1.0d0/INDUCTI(27)
        A(54,54) = -RESISTI(27)/INDUCTI(27)

        ! --- LINHA 55 --- (bifurcacao: Q_29=col58, Q_30=col60)
        A(55,56) =  1.0d0/COMPLII(28)
        A(55,58) = -1.0d0/COMPLII(28)
        A(55,60) = -1.0d0/COMPLII(28)

        ! --- LINHA 56 ---
        A(56,51) =  1.0d0/INDUCTI(28)
        A(56,55) = -1.0d0/INDUCTI(28)
        A(56,58) = -RESISTI(28)/INDUCTI(28)
        A(56,60) = -RESISTI(28)/INDUCTI(28)

        ! --- LINHA 57 --- (terminal)
        A(57,58) = 1.0d0/COMPLII(29)

        ! --- LINHA 58 --- (terminal: upstream P=col55)
        A(58,55) =  1.0d0/INDUCTI(29)
        A(58,57) = -1.0d0/INDUCTI(29)
        A(58,58) = -RESISTI(29)/INDUCTI(29)

        ! --- LINHA 59 --- (bifurcacao: Q_31=col62, Q_32=col64)
        A(59,60) =  1.0d0/COMPLII(30)
        A(59,62) = -1.0d0/COMPLII(30)
        A(59,64) = -1.0d0/COMPLII(30)

        ! --- LINHA 60 ---
        A(60,55) =  1.0d0/INDUCTI(30)
        A(60,59) = -1.0d0/INDUCTI(30)
        A(60,62) = -RESISTI(30)/INDUCTI(30)
        A(60,64) = -RESISTI(30)/INDUCTI(30)

        ! --- LINHA 61 --- (terminal)
        A(61,62) = 1.0d0/COMPLII(31)

        ! --- LINHA 62 --- (terminal: upstream P=col59)
        A(62,59) =  1.0d0/INDUCTI(31)
        A(62,61) = -1.0d0/INDUCTI(31)
        A(62,62) = -RESISTI(31)/INDUCTI(31)

        ! --- LINHA 63 --- (terminal)
        A(63,64) = 1.0d0/COMPLII(32)

        ! --- LINHA 64 --- (terminal: upstream P=col59)
        A(64,59) =  1.0d0/INDUCTI(32)
        A(64,63) = -1.0d0/INDUCTI(32)
        A(64,64) = -RESISTI(32)/INDUCTI(32)

        ! --- LINHA 65 --- (terminal)
        A(65,66) = 1.0d0/COMPLII(33)

        ! --- LINHA 66 --- (terminal: upstream P=col49)
        A(66,49) =  1.0d0/INDUCTI(33)
        A(66,65) = -1.0d0/INDUCTI(33)
        A(66,66) = -RESISTI(33)/INDUCTI(33)

        ! --- LINHA 67 ---
        A(67,68) =  1.0d0/COMPLII(34)
        A(67,70) = -1.0d0/COMPLII(34)
        A(67,72) = -1.0d0/COMPLII(34)

        ! --- LINHA 68 ---
        A(68,49) =  1.0d0/INDUCTI(34)
        A(68,67) = -1.0d0/INDUCTI(34)
        A(68,70) = -RESISTI(34)/INDUCTI(34)
        A(68,72) = -RESISTI(34)/INDUCTI(34)

        ! --- LINHA 69 --- (bifurcacao: Q_40=col80, Q_41=col82)
        A(69,70) =  1.0d0/COMPLII(35)
        A(69,80) = -1.0d0/COMPLII(35)
        A(69,82) = -1.0d0/COMPLII(35)

        ! --- LINHA 70 ---
        A(70,67) =  1.0d0/INDUCTI(35)
        A(70,69) = -1.0d0/INDUCTI(35)
        A(70,80) = -RESISTI(35)/INDUCTI(35)
        A(70,82) = -RESISTI(35)/INDUCTI(35)

        ! --- LINHA 71 --- (bifurcacao: Q_37=col74, Q_38=col76, Q_39=col78)
        A(71,72) =  1.0d0/COMPLII(36)
        A(71,74) = -1.0d0/COMPLII(36)
        A(71,76) = -1.0d0/COMPLII(36)
        A(71,78) = -1.0d0/COMPLII(36)

        ! --- LINHA 72 ---
        A(72,67) =  1.0d0/INDUCTI(36)
        A(72,71) = -1.0d0/INDUCTI(36)
        A(72,74) = -RESISTI(36)/INDUCTI(36)
        A(72,76) = -RESISTI(36)/INDUCTI(36)
        A(72,78) = -RESISTI(36)/INDUCTI(36)

        ! --- LINHA 73 --- (terminal)
        A(73,74) = 1.0d0/COMPLII(37)

        ! --- LINHA 74 --- (terminal: upstream P=col71)
        A(74,71) =  1.0d0/INDUCTI(37)
        A(74,73) = -1.0d0/INDUCTI(37)
        A(74,74) = -RESISTI(37)/INDUCTI(37)

        ! --- LINHA 75 --- (terminal)
        A(75,76) = 1.0d0/COMPLII(38)

        ! --- LINHA 76 --- (terminal: upstream P=col71)
        A(76,71) =  1.0d0/INDUCTI(38)
        A(76,75) = -1.0d0/INDUCTI(38)
        A(76,76) = -RESISTI(38)/INDUCTI(38)

        ! --- LINHA 77 --- (terminal)
        A(77,78) = 1.0d0/COMPLII(39)

        ! --- LINHA 78 --- (terminal: upstream P=col71)
        A(78,71) =  1.0d0/INDUCTI(39)
        A(78,77) = -1.0d0/INDUCTI(39)
        A(78,78) = -RESISTI(39)/INDUCTI(39)

        ! --- LINHA 79 --- (terminal)
        A(79,80) = 1.0d0/COMPLII(40)

        ! --- LINHA 80 --- (terminal: upstream P=col69)
        A(80,69) =  1.0d0/INDUCTI(40)
        A(80,79) = -1.0d0/INDUCTI(40)
        A(80,80) = -RESISTI(40)/INDUCTI(40)

        ! --- LINHA 81 --- (bifurcacao: Q_42=col84, Q_43=col86)
        A(81,82) =  1.0d0/COMPLII(41)
        A(81,84) = -1.0d0/COMPLII(41)
        A(81,86) = -1.0d0/COMPLII(41)

        ! --- LINHA 82 ---
        A(82,69) =  1.0d0/INDUCTI(41)
        A(82,81) = -1.0d0/INDUCTI(41)
        A(82,84) = -RESISTI(41)/INDUCTI(41)
        A(82,86) = -RESISTI(41)/INDUCTI(41)

        ! --- LINHA 83 --- (terminal)
        A(83,84) = 1.0d0/COMPLII(42)

        ! --- LINHA 84 --- (terminal: upstream P=col81)
        A(84,81) =  1.0d0/INDUCTI(42)
        A(84,83) = -1.0d0/INDUCTI(42)
        A(84,84) = -RESISTI(42)/INDUCTI(42)

        ! --- LINHA 85 --- (bifurcacao: Q_44=col88, Q_45=col90)
        A(85,86) =  1.0d0/COMPLII(43)
        A(85,88) = -1.0d0/COMPLII(43)
        A(85,90) = -1.0d0/COMPLII(43)

        ! --- LINHA 86 ---
        A(86,81) =  1.0d0/INDUCTI(43)
        A(86,85) = -1.0d0/INDUCTI(43)
        A(86,88) = -RESISTI(43)/INDUCTI(43)
        A(86,90) = -RESISTI(43)/INDUCTI(43)

        ! --- LINHA 87 --- (terminal)
        A(87,88) = 1.0d0/COMPLII(44)

        ! --- LINHA 88 --- (terminal: upstream P=col85)
        A(88,85) =  1.0d0/INDUCTI(44)
        A(88,87) = -1.0d0/INDUCTI(44)
        A(88,88) = -RESISTI(44)/INDUCTI(44)

        ! --- LINHA 89 --- (bifurcacao: Q_46=col92, Q_47=col94)
        A(89,90) =  1.0d0/COMPLII(45)
        A(89,92) = -1.0d0/COMPLII(45)
        A(89,94) = -1.0d0/COMPLII(45)

        ! --- LINHA 90 ---
        A(90,85) =  1.0d0/INDUCTI(45)
        A(90,89) = -1.0d0/INDUCTI(45)
        A(90,92) = -RESISTI(45)/INDUCTI(45)
        A(90,94) = -RESISTI(45)/INDUCTI(45)

        ! --- LINHA 91 --- (terminal)
        A(91,92) = 1.0d0/COMPLII(46)

        ! --- LINHA 92 --- (terminal: upstream P=col89)
        A(92,89) =  1.0d0/INDUCTI(46)
        A(92,91) = -1.0d0/INDUCTI(46)
        A(92,92) = -RESISTI(46)/INDUCTI(46)

        ! --- LINHA 93 ---
        A(93,94) = 1.0d0/COMPLII(47)
        A(93,96) = -1.0d0/COMPLII(47)
        A(93,98) = -1.0d0/COMPLII(47)

        ! --- LINHA 94 ---
        A(94,89) = 1.0d0/INDUCTI(47)
        A(94,93) = -1.0d0/INDUCTI(47)
        A(94,96) = -RESISTI(47)/INDUCTI(47)
        A(94,98) = -RESISTI(47)/INDUCTI(47)

        ! --- LINHA 95 ---
        A(95,96)  =  1.0d0/COMPLII(48)
        A(95,100) = -1.0d0/COMPLII(48)
        A(95,102) = -1.0d0/COMPLII(48)

        ! --- LINHA 96 ---
        A(96,93)  =  1.0d0/INDUCTI(48)
        A(96,95)  = -1.0d0/INDUCTI(48)
        A(96,100) = -RESISTI(48)/INDUCTI(48)
        A(96,102) = -RESISTI(48)/INDUCTI(48)

        ! --- LINHA 97 ---
        A(97,98)  =  1.0d0/COMPLII(49)
        A(97,112) = -1.0d0/COMPLII(49)
        A(97,114) = -1.0d0/COMPLII(49)

        ! --- LINHA 98 ---
        A(98,93)  =  1.0d0/INDUCTI(49)
        A(98,97)  = -1.0d0/INDUCTI(49)
        A(98,112) = -RESISTI(49)/INDUCTI(49)
        A(98,114) = -RESISTI(49)/INDUCTI(49)

        ! --- LINHA 99 ---
        A(99,100) =  1.0d0/COMPLII(50)
        A(99,104) = -1.0d0/COMPLII(50)
        A(99,106) = -1.0d0/COMPLII(50)

        ! --- LINHA 100 ---
        A(100,95)  =  1.0d0/INDUCTI(50)
        A(100,99)  = -1.0d0/INDUCTI(50)
        A(100,104) = -RESISTI(50)/INDUCTI(50)
        A(100,106) = -RESISTI(50)/INDUCTI(50)

        ! --- LINHA 101 --- (terminal)
        A(101,102) = 1.0d0/COMPLII(51)

        ! --- LINHA 102 --- (terminal: upstream P=col95)
        A(102,95)  =  1.0d0/INDUCTI(51)
        A(102,101) = -1.0d0/INDUCTI(51)
        A(102,102) = -RESISTI(51)/INDUCTI(51)

        ! --- LINHA 103 --- (bifurcacao: Q_54=col108, Q_55=col110)
        A(103,104) =  1.0d0/COMPLII(52)
        A(103,108) = -1.0d0/COMPLII(52)
        A(103,110) = -1.0d0/COMPLII(52)

        ! --- LINHA 104 ---
        A(104,100)  =  1.0d0/INDUCTI(52)
        A(104,104) = -1.0d0/INDUCTI(52)
        A(104,108) = -RESISTI(52)/INDUCTI(52)
        A(104,110) = -RESISTI(52)/INDUCTI(52)

        ! --- LINHA 105 --- (terminal)
        A(105,106) = 1.0d0/COMPLII(53)

        ! --- LINHA 106 --- (terminal: upstream P=col99)
        A(106,99)  =  1.0d0/INDUCTI(53)
        A(106,105) = -1.0d0/INDUCTI(53)
        A(106,106) = -RESISTI(53)/INDUCTI(53)

        ! --- LINHA 107 --- (terminal)
        A(107,108) = 1.0d0/COMPLII(54)

        ! --- LINHA 108 --- (terminal: upstream P=col103)
        A(108,103) =  1.0d0/INDUCTI(54)
        A(108,107) = -1.0d0/INDUCTI(54)
        A(108,108) = -RESISTI(54)/INDUCTI(54)

        ! --- LINHA 109 --- (terminal)
        A(109,110) = 1.0d0/COMPLII(55)

        ! --- LINHA 110 --- (terminal: upstream P=col103)
        A(110,103) =  1.0d0/INDUCTI(55)
        A(110,109) = -1.0d0/INDUCTI(55)
        A(110,110) = -RESISTI(55)/INDUCTI(55)

        ! --- LINHA 111 --- (bifurcacao: Q_59=col118, Q_58=col116)
        A(111,112) =  1.0d0/COMPLII(56)
        A(111,116) = -1.0d0/COMPLII(56)
        A(111,118) = -1.0d0/COMPLII(56)

        ! --- LINHA 112 ---
        A(112,97)  =  1.0d0/INDUCTI(56)
        A(112,111) = -1.0d0/INDUCTI(56)
        A(112,116) = -RESISTI(56)/INDUCTI(56)
        A(112,118) = -RESISTI(56)/INDUCTI(56)

        ! --- LINHA 113 --- (terminal)
        A(113,114) = 1.0d0/COMPLII(57)

        ! --- LINHA 114 --- (terminal: upstream P=col97)
        A(114,97)  =  1.0d0/INDUCTI(57)
        A(114,113) = -1.0d0/INDUCTI(57)
        A(114,114) = -RESISTI(57)/INDUCTI(57)

        ! --- LINHA 115 --- (bifurcacao: Q_60=col120, Q_61=col122)
        A(115,116) =  1.0d0/COMPLII(58)
        A(115,120) = -1.0d0/COMPLII(58)
        A(115,122) = -1.0d0/COMPLII(58)

        ! --- LINHA 116 ---
        A(116,111) =  1.0d0/INDUCTI(58)
        A(116,115) = -1.0d0/INDUCTI(58)
        A(116,120) = -RESISTI(58)/INDUCTI(58)
        A(116,122) = -RESISTI(58)/INDUCTI(58)

        ! --- LINHA 117 --- (terminal)
        A(117,118) = 1.0d0/COMPLII(59)

        ! --- LINHA 118 --- (terminal: upstream P=col111)
        A(118,111) =  1.0d0/INDUCTI(59)
        A(118,117) = -1.0d0/INDUCTI(59)
        A(118,118) = -RESISTI(59)/INDUCTI(59)

        ! --- LINHA 119 --- (terminal)
        A(119,120) = 1.0d0/COMPLII(60)

        ! --- LINHA 120 --- (terminal: upstream P=col115)
        A(120,115) =  1.0d0/INDUCTI(60)
        A(120,119) = -1.0d0/INDUCTI(60)
        A(120,120) = -RESISTI(60)/INDUCTI(60)

        ! --- LINHA 121 --- (terminal)
        A(121,122) = 1.0d0/COMPLII(61)

        ! --- LINHA 122 --- (terminal: upstream P=col115)
        A(122,115) =  1.0d0/INDUCTI(61)
        A(122,121) = -1.0d0/INDUCTI(61)
        A(122,122) = -RESISTI(61)/INDUCTI(61)

        call PULSOPQ(IT, HT, PI, PIN, QIN)
        b(1) = PIN / COMPLII(1)
        b(2) = QIN / INDUCTI(1)

    end subroutine FORMAB

    subroutine MIVEULER(n, A, b, HT, x_old, x_new)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: HT
        real(kind=8), intent(in) :: A(:,:)
        real(kind=8), intent(in) :: b(:), x_old(:)
        real(kind=8), intent(out) :: x_new(:)
        real(kind=8) :: dx(122)
        integer :: i, j
        integer :: dim

        dim = 2*n
        if (dim > 122) stop "ERRO: dimensao excede tamanho do array"

        dx = 0.0d0
        do i = 1, dim
            do j = 1, dim
                dx(i) = dx(i) + A(i,j) * x_old(j)
            end do
            dx(i) = dx(i) + b(i)
        end do

        do i = 1, dim
            x_new(i) = x_old(i) + HT * dx(i)
        end do
    end subroutine MIVEULER

end module routines_p

program ZERODMOD
    use routines_p
    implicit none

    integer :: n, IT, NT
    integer, dimension(10) :: indices = [16, 19, 20, 23, 24, 30, 48, 49, 52, 58]
    real(kind=8) :: MI, RO, PI, HT, TT
    real(kind=8), allocatable :: COMPLI(:), RAIOI(:), ESPHI(:), MEY(:)
    real(kind=8), allocatable :: RESISTI(:), INDUCTI(:), COMPLII(:)
    real(kind=8), allocatable :: A(:,:), b(:), x_old(:), x_new(:), resultados(:,:)

    write(*,*) "Iniciando simulacao..."

    call DADOSMOD(n, COMPLI, RAIOI, ESPHI, MEY, MI, RO, PI)
    write(*,*) "Trechos lidos: ", n

    allocate(RESISTI(n), INDUCTI(n), COMPLII(n))
    call GEOMAT(n, COMPLI, RAIOI, ESPHI, MEY, MI, RO, PI, RESISTI, INDUCTI, COMPLII)
    write(*,*) "DAdos Geometricos Calculados."

    HT = 0.001d0
    TT = 1.0d0
    NT = int(TT / HT)
    write(*,*) "Integracao temporal (", NT, " passos)..."

    allocate(A(2*n, 2*n), b(2*n), x_old(2*n), x_new(2*n))
    allocate(resultados(NT+1, 1+10*2))

    x_old = 0.0d0
    do IT = 0, NT
        call FORMAB(n, RESISTI, INDUCTI, COMPLII, HT, PI, IT, A, b)
        call MIVEULER(n, A, b, HT, x_old, x_new)

        call store_results(IT+1, x_new, indices, resultados)

        x_old = x_new
    end do

    write(*,*) "Simulacao concluida."
    call save_results(resultados, NT)
    call save_info(HT, TT, NT, n, MI, RO, PI, indices)

contains

    subroutine store_results(idx, x, indices, res)
        integer, intent(in) :: idx
        real(kind=8), intent(in) :: x(122)
        integer, intent(in) :: indices(10)  ! CORRECTED: now declared as integer
        real(kind=8), intent(inout) :: res(:,:)
        integer :: k, j
        do k = 1, 10
            j = indices(k)
            res(idx, 2 + (k-1)*2) = x(2*j - 1)  ! Pressao
            res(idx, 3 + (k-1)*2) = x(2*j)      ! Vazao
        end do
    end subroutine store_results

    subroutine save_results(res, NT)
        real(kind=8), intent(in) :: res(:,:)
        integer, intent(in) :: NT
        integer :: IT, ios
        character(len=30), parameter :: filename = 'resultados_simulacao.csv'
        
        open(unit=30, file=filename, status='replace', iostat=ios)
        if (ios /= 0) stop "ERRO: Falha ao criar arquivo de resultados."
        
        write(30,'(A)') 'tempo,P_16,Q_16,P_19,Q_19,P_20,Q_20,P_23,Q_23,P_24,Q_24,P_30,Q_30,P_48,Q_48,P_49,Q_49,P_52,Q_52,P_58,Q_58'
        do IT = 1, NT+1
            write(30,'(*(ES15.6,","))') res(IT,:)
        end do
        close(30)
        write(*,*) "Resultados salvos: ", trim(filename)
    end subroutine save_results

    subroutine save_info(HT, TT, NT, n, MI, RO, PI, indices)
        real(kind=8), intent(in) :: HT, TT, MI, RO, PI
        integer, intent(in) :: NT, n, indices(10)
        integer :: ios
        character(len=20), parameter :: filename = 'info_simulacao.csv'
        
        open(unit=31, file=filename, status='replace', iostat=ios)
        if (ios /= 0) stop "ERRO: Falha ao criar arquivo de info."
        
        write(31,'(A)') 'Parametro,Valor'
        write(31,'(A,ES15.6)') 'HT,', HT
        write(31,'(A,ES15.6)') 'TT,', TT
        write(31,'(A,I0)') 'NT,', NT
        write(31,'(A,I0)') 'n,', n
        write(31,'(A,ES15.6)') 'MI,', MI
        write(31,'(A,ES15.6)') 'RO,', RO
        write(31,'(A,ES15.6)') 'PI,', PI
        write(31,'(A)') 'Indices_monitorados,16,19,20,23,24,30,48,49,52,58'
        close(31)
        write(*,*) "Info salva: ", trim(filename)
    end subroutine save_info

end program ZERODMOD