
  subroutine flavor(ID, IFL1, IFL2, IFL3, JSPIN, INDX)

! ======================================================================
!
!          THIS SUBROUTINE UNPACKS THE IDENT CODE ID=+/-IJKL
!          I.E., GIVEN THE PARTICLE "ID", LOOKUP CHARACTERISTICS AND
!          RETURN
!
!          MESONS--
!          I=0, J<=K, +/- IS SIGN FOR J
!          ID=110 FOR PI0, ID=220 FOR ETA, ETC.
!
!          BARYONS--
!          I<=J<=K IN GENERAL
!          J<I<K FOR SECOND STATE ANTISYMMETRI! IN (I,J), EG. L = 2130
!
!          OTHER--
!          ID=1,...,6 FOR QUARKS
!          ID=9 FOR GLUON
!          ID=10 FOR PHOTON
!          ID=11,...,16 FOR LEPTONS
!          ID=20 FOR KS, ID=-20 FOR KL
!
!          I=21...26 FOR SCALAR QUARKS
!          I=29 FOR GLUINO
!          I=30 FOR PHOTINO
!          I=31...36 FOR SCALAR LEPTONS
!          I=39 FOR WINO
!          I=40 FOR ZINO
!
!          ID=80 FOR W+
!          ID=81,...,89 FOR HIGGS MESONS
!          ID=90 FOR Z0
!
!          DIQUARKS--
!          ID=+/-IJ00, I<J FOR DIQUARK COMPOSED OF I,J.
!
!          INDX IS A SEQUENCE NUMBER USED INTERNALLY
!
! ======================================================================

    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32), intent(in   ) :: ID
    integer(int32), intent(  out) :: IFL1
    integer(int32), intent(  out) :: IFL2
    integer(int32), intent(  out) :: IFL3
    integer(int32), intent(  out) :: JSPIN
    integer(int32), intent(  out), optional :: INDX

    ! Random number generator calls
    real(real64) :: rndm

    ! Interim variables
    integer(int32) :: &
        & IDABS   ! Absolute value of the reaction ID
    integer(int32) :: I, J, K, switchType, switchType2

! ======================================================================

    ! AMELP not used
    real(real64) :: amlep
    integer(int32) :: nqlep, nmes, nbary
    COMMON/QLMASS/ AMLEP(52),NQLEP,NMES,NBARY
!@@@@@@@@@@   SIVOKL  @@@@@@@@@@
    integer(int32) :: nfla, nfl1, nfl2, nfl3, nspin, nndex
    COMMON/FLACOM/NFLA,NFL1,NFL2,NFL3,NSPIN,NNDEX
! ======================================================================

    ! Unpack characteristics associated to the reaction ID
    IDABS = IABS(ID)
    I = IDABS / 1000
    J = MOD(IDABS / 100, 10)
    K = MOD(IDABS /  10, 10)
    JSPIN = MOD(IDABS, 10)

    ! Evaluate which processing should be executed for the ID
    switchType = 0
    IF(ID /= 0 .AND. MOD(ID, 100) == 0) switchType = 300
    IF(J == 0 .AND. switchType == 0) switchType = 200
    IF(I == 0 .AND. switchType == 0) switchType = 100

    select case (switchType)
       case (100);
          ! Mesons
          IFL1 = 0
          IFL2 = ISIGN(J,  ID)
          IFL3 = ISIGN(K, -ID)
          INDX = J + K * (K - 1) / 2 + 36 * JSPIN + NQLEP + 11

          !@@@@@@@@@@@@@ SIVOKL - TONEEV @@@@@
          switchType2 = 0
          IF(ID == 110 .OR. ID == 111 .OR. ID == 221) switchType2 = 13
          IF((ID == 220 .OR. ID == 330) .AND. switchType2 == 0) switchType2 = 12

          select case (switchType2)
             case (12);
                IFL2 = 2 +INT(0.25 + RNDM(-1.))
                IF(IFL2 == 2) IFL2 = 1 + INT(0.5 + RNDM(-1.))
                IFL2 = ISIGN(IFL2,  ID)
                IFL3 = ISIGN(IFL2, -ID)

                IF(NFLA == -1) THEN
                   NFL1 = IFL1
                   NFL2 = IFL2
                   NFL3 = IFL3
                   NSPIN = JSPIN
                   NNDEX = INDX
                   NFLA = ID
                ENDIF

             case (13);
                IFL2 = 1 + INT(0.5 + RNDM(-1.))
                IFL2 = ISIGN(IFL2,ID)
                IFL3 = ISIGN(IFL2,-ID)

                IF(NFLA == -1) THEN
                   NFL1 = IFL1
                   NFL2 = IFL2
                   NFL3 = IFL3
                   NSPIN = JSPIN
                   NNDEX = INDX
                   NFLA = ID
! in the next line NFLB was not defined 11.16.94 V.T. !
!                 ELSE IF(NFLB == -1) THEN
!                    MFL1 = IFL1
!                    MFL2 = IFL2
!                    MFL3 = IFL3
!                    MSPIN = JSPIN
!                    MNDEX = INDX
!                    NFLB = ID
                ENDIF
          end select

       case (200);
          IFL1 = 0
          IFL2 = 0
          IFL3 = 0
          JSPIN = 0
          INDX = IDABS
          IF(IDABS < 20) RETURN

          ! DEFINE INDX=20 FOR KS, INDX=21 FOR KL
          INDX = IDABS + 1
          IF(ID == 20) INDX = 20
          ! INDX = NQLEP + 1, ... , NQLEP + 11 FOR W+, HIGGS, Z0
          IF(IDABS < 80) RETURN
          INDX = NQLEP + IDABS - 79

      case (300);
          IFL1 = ISIGN(I, ID)
          IFL2 = ISIGN(J, ID)
          IFL3 = 0
          JSPIN=0
          INDX=0

      case default;
          ! BARYONS
          !    ONLY X,Y BARYONS ARE QQX, QQY, Q=U,D,S.
          IFL1=ISIGN(I,ID)
          IFL2=ISIGN(J,ID)
          IFL3=ISIGN(K,ID)

          ! Determine INDX identifier
          INDX = MAX0(I - 1, J - 1)**2 + I + MAX0(I - J, 0)
          IF(K <= 6) then
             INDX = INDX + K * (K - 1) * (2 * K - 1) / 6
          else
             INDX = INDX + 9 * (K - 7) + 91
          end if
          INDX = INDX + 109 * JSPIN + 36 * NMES + NQLEP + 11

    end select

    return
  end subroutine flavor

