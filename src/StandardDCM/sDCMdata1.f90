
! ======================================================================
!
!     This subroutine sets up the coefficients ankj, bnkj, and ckj
!     for the main program.
!     ankj is related to angular distributions for elementary cross
!     sections in the INC.
!     Ankj used by COSTA;
!     Bnjk & Ckj used by PMOM; give momentum distributions of final-
!     state products in pion-production reactions.
!
!   The baryon channels are given in Barashenkov and Toneev's book
!   [in Russian] (1972).
!   For the photon channels, see
!   Barashenkov, et al., Nucl Phys. A231, 462 (1974).
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL T-2, February, 1996.
!    Edited by AJS, July-August, 1997.
!    Comments added by AJS, December, 1998.
!    Modified:    04-AUG-1998 BY NVMokhov
!    "Last" CHANGE: 12-AUG-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

     integer(int32) ::  k, n

! ======================================================================

     real(real64), public, protected, dimension(4, 4, 29) :: ankj
     real(real64), public, protected, dimension(4, 4,  8) :: bnkj
     real(real64), public, protected, dimension(3, 8)     :: ckj

! ======================================================================

!  ankj Angular distribution coefficients:
!  j = 1;  N + N elastic scattering; Tlab <= 2.8 GeV:
!          (n + n & p + p isotropic below 0.46 GeV.)
     data ((ankj(n,k, 1),n=1,4),k=1,4) / &
          & 2.7404d+0 , -9.6998d+0 ,  1.0400d+1 ,  2.3882d+0 , &
          & -7.5137d+0 ,  4.4096d+1 , -7.4379d+1 ,  4.6038d+1 , &
          & 7.5479d+0 , -3.9274d+1 ,  6.4835d+1 , -4.1609d+1 , &
          & -1.8369d+0 ,  8.6911d+0 , -1.3060d+1 ,  7.1880d+0 /

!  j = 2;  N + N elastic scattering; 2.8 < Tlab <= 10. GeV:
     data ((ankj(n,k, 2),n=1,4),k=1,4) / &
          & -3.0853d+1 ,  1.0624d+2 , -1.2939d+2 ,  5.4339d+1 , &
          & 1.9465d+1 , -6.8102d+1 ,  9.6358d+1 , -5.6827d+1 , &
          & -3.4831d+0 ,  1.2341d+1 , -1.8592d+1 ,  1.2024d+1 , &
          & 1.8941d-1 , -6.7880d-1 ,  1.0665d+0 , -7.2910d-1 /

!  j = 3;  n + p elastic scattering; Tlab <= 0.97 GeV:
     data ((ankj(n,k, 3),n=1,4),k=1,4) / &
          & 1.0258d-1 , -1.0542d+0 ,  1.1389d+1 , -1.6638d+1 , &
          & -4.9607d-1 ,  1.1800d+1 , -9.0857d+1 ,  1.6476d+2 , &
          & 1.5437d+0 , -3.3769d+1 ,  2.5192d+2 , -4.5071d+2 , &
          & -1.2021d+0 ,  2.5336d+1 , -1.8658d+2 ,  3.3254d+2 /

!  j = 4; pi+ p or pi- n elastic scattering; Tlab <= 0.080 GeV:
     data ((ankj(n,k, 4),n=1,4),k=1,4) / &
          & 1.5789d-1 ,  2.9671d+0 , -5.5251d+0 ,  6.8925d+0 , &
          & -7.0218d+0 , -2.0534d+2 ,  5.6951d+2 , -8.9858d+2 , &
          & 1.3496d+2 ,  4.8722d+3 , -1.4674d+4 ,  2.3924d+4 , &
          & -8.2116d+2 , -3.2586d+4 ,  1.0098d+5 , -1.6553d+5 /

!  j = 5; pi+ p or pi- n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
     data ((ankj(n,k, 5),n=1,4),k=1,4) / &
          & 3.1531d-1 , -7.4981d+0 ,  4.3295d+1 , -7.6360d+1 , &
          & -6.5373d+0 ,  1.9307d+2 , -1.0181d+3 ,  1.7426d+3 , &
          & 4.6864d+1 , -1.3030d+3 ,  6.7291d+3 , -1.1075d+4 , &
          & -9.5192d+1 ,  2.6373d+3 , -1.2857d+4 ,  2.0294d+4 /

!  j = 6; pi+ p or pi- n elastic scattering; 0.30 < Tlab <= 1.0 GeV:
     data ((ankj(n,k, 6),n=1,4),k=1,4) / &
          & -1.7953d+1 ,  1.0972d+2 , -2.3954d+2 ,  2.2826d+2 , &
          & 9.1968d+1 , -5.1963d+2 ,  1.1266d+3 , -1.0740d+3 , &
          & -1.3270d+2 ,  7.4112d+2 , -1.6000d+3 ,  1.5249d+3 , &
          & 5.8598d+1 , -3.1874d+2 ,  6.7751d+2 , -6.4011d+2 /

!  j = 7; pi+ p or pi- n elastic scattering; 1.0 < Tlab <= 2.4 GeV:
     data ((ankj(n,k, 7),n=1,4),k=1,4) / &
          & 4.2169d-1 ,  1.4705d+2 , -6.5335d+2 ,  9.1507d+2 , &
          & -3.5198d+0 , -2.6019d+2 ,  1.2250d+3 , -1.7481d+3 , &
          & 3.6373d+0 ,  1.5592d+2 , -7.5201d+2 ,  1.0796d+3 , &
          & -7.8041d-1 , -3.0563d+1 ,  1.4795d+2 , -2.1250d+2 /

!  j = 8; pi+ n or pi- p elastic scattering; Tlab <= 0.080 GeV:
     data ((ankj(n,k, 8),n=1,4),k=1,4) / &
          & -3.8288d-1 ,  3.7587d+0 , -6.5144d+0 ,  6.7740d+0 , &
          & 1.0381d+2 , -2.7282d+2 ,  4.7759d+2 , -5.1222d+2 , &
          & -1.7882d+3 ,  4.3052d+3 , -7.9314d+3 ,  9.3471d+3 , &
          & 7.1475d+3 , -3.3395d+3 , -4.1392d+3 , -4.4364d+3 /

!  j = 9; pi- p or pi+ n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
     data ((ankj(n,k, 9),n=1,4),k=1,4) / &
          & 2.4991d-1 ,  3.2028d+1 , -1.1882d+2 ,  1.5099d+2 , &
          & -2.6994d+0 , -4.6045d+2 ,  1.8959d+3 , -2.5190d+3 , &
          & 1.6268d+1 ,  2.1384d+3 , -9.1262d+3 ,  1.2431d+4 , &
          & -2.9654d+1 , -3.1823d+3 ,  1.3944d+4 , -1.9342d+4 /

!  j = 10; pi- p or pi+ n elastic or CX scattering;
!          0.30 < Tlab <= 1.0 GeV:
     data ((ankj(n,k,10),n=1,4),k=1,4) / &
          & 3.9025d+0 , -9.1126d+1 ,  3.2373d+2 , -4.0048d+2 , &
          & -2.0619d+1 ,  4.9170d+2 , -1.7155d+3 ,  2.1143d+3 , &
          & 3.3004d+1 , -7.6684d+2 ,  2.7003d+3 , -3.3525d+3 , &
          & -1.6367d+1 ,  3.7394d+2 , -1.3202d+3 ,  1.6423d+3 /

!  j = 11; pi- p or pi+ n elastic or CX scattering;
!          1.0 < Tlab <= 2.4 GeV:
     data ((ankj(n,k,11),n=1,4),k=1,4) / &
          & 1.9402d+1 , -2.2446d+2 ,  7.4733d+2 , -9.3570d+2 , &
          & -4.4180d+1 ,  4.7194d+2 , -1.4856d+3 ,  1.8055d+3 , &
          & 3.1567d+1 , -3.0176d+2 ,  9.0763d+2 , -1.0773d+3 , &
          & -6.8648d+0 ,  6.0476d+1 , -1.7520d+2 ,  2.0381d+2 /

!  j = 12;  gamma + N 'elastic' scattering, E-g <= 0.45 GeV.
!    gam + p --> p + pi0
     data ((ankj(n,k,12),n=1,4),k=1,4) / &
          & 4.0693d-1 , -4.1404d+0 ,  1.4044d+1 , -1.7265d+1 , &
          & -3.6799d+0 ,  5.9610d+1 , -1.6269d+2 ,  1.8873d+2 , &
          & 1.4556d+1 , -1.7550d+2 ,  4.5839d+2 , -5.3390d+2 , &
          & -1.2621d+1 ,  1.4964d+2 , -3.8118d+2 ,  4.5141d+2 /

!  j = 13;  gamma + N 'elastic' scattering, E-g > 0.45 GeV.
!    gam + p --> p + pi0
     data ((ankj(n,k,13),n=1,4),k=1,4) / &
          & -4.7554d-1 ,  2.2641d+0 , -1.2528d+1 ,  2.4647d+1 , &
          & 5.1620d+0 , -9.9236d+0 ,  5.5623d+1 , -1.0462d+2 , &
          & -8.1117d+0 ,  1.9315d+1 , -8.4255d+1 ,  1.3908d+2 , &
          & 3.5187d+0 , -9.1783d+0 ,  3.4950d+1 , -5.1243d+1 /

!  j = 14; gamma + p --> n + pi+; E-g <= 0.51 GeV:
     data ((ankj(n,k,14),n=1,4),k=1,4) / &
          & 4.8173d-1 ,  5.7726d+0 , -1.3745d+1 ,  2.7125d+1 , &
          & -4.4804d+0 , -3.8582d+1 ,  1.1159d+2 , -2.4305d+2 , &
          & 1.6306d+1 ,  1.1046d+2 , -3.3045d+2 ,  7.2270d+2 , &
          & -1.5968d+1 , -8.0140d+1 ,  2.4616d+2 , -6.0753d+2 /

!  j = 15; gamma + p --> n + pi+; 0.51 < E-g <= 1.0 GeV:
     data ((ankj(n,k,15),n=1,4),k=1,4) / &
          & -5.1646d+0 , -6.0776d+0 ,  7.8989d+1 , -1.0705d+2 , &
          & 2.1871d+1 ,  5.6915d+1 , -4.0159d+2 ,  5.1215d+2 , &
          & -2.7993d+1 , -9.4670d+1 ,  5.6928d+2 , -6.9621d+2 , &
          & 1.1587d+1 ,  4.5998d+1 , -2.4566d+2 ,  2.8452d+2 /

!  j = 16; gamma + p --> n + pi+; 1.0 < E-g <= 10 GeV:
     data ((ankj(n,k,16),n=1,4),k=1,4) / &
          & -5.3067d+1 ,  5.7612d+2 , -1.5438d+3 ,  1.6455d+5 , &
          & 1.4750d+2 , -1.6638d+3 ,  4.5923d+3 , -4.9949d+3 , &
          & -1.3436d+2 ,  1.5780d+3 , -4.4463d+3 ,  4.9022d+3 , &
          & 4.0253d+1 , -4.8860d+2 ,  1.4001d+3 , -1.5606d+3 /

!  j = 17; pi- + p --> pi0 n or pi+ + n --> pi0 + p Charge exchange
!          scattering; Tlab <= 0.08 GeV:
     data ((ankj(n,k,17),n=1,4),k=1,4) / &
          & 1.4988d-1 ,  2.8753d+0 , -5.3078d+0 ,  6.2233d+0 , &
          & -5.9558d+0 , -1.6203d+2 ,  4.3079d+2 , -6.2548d+2 , &
          & 1.2875d+2 ,  3.1402d+3 , -7.9189d+3 ,  1.0983d+4 , &
          & -8.5161d+2 , -1.8780d+4 ,  4.4607d+4 , -5.8790d+4 /

!  j = 18; pi- + p --> pi0 n or pi+ + n --> pi0 + p Charge exchange
!          scattering; 0.08 < Tlab <= 0.30 GeV:
     data ((ankj(n,k,18),n=1,4),k=1,4) / &
          & 5.3689d-1 , -1.3216d+1 ,  8.1011d+1 , -1.4285d+2 , &
          & -1.0550d+1 ,  2.9629d+2 , -1.6957d+3 ,  2.8935d+3 , &
          & 6.9621d+1 , -1.9245d+3 ,  1.0620d+4 , -1.7468d+4 , &
          & -1.3865d+2 ,  3.9281d+3 , -2.0293d+4 ,  3.2058d+4 /

!  j = 19; coefficients for gamma absorption on 2 nucleons; outgoing
!          nucleon ang. dist., Tlab <= 0.455
     data ((ankj(n,k,19),n=1,4),k=1,4) / &
          & 6.5288d-1 ,  3.8977d-1 ,  8.4078d-1 ,  1.8893d-1 , &
          & -4.3964d+0 ,  3.4309d+1 , -7.3692d+1 ,  8.4308d+1 , &
          & 1.4889d+1 , -1.4380d+2 ,  3.1227d+2 , -3.5014d+2 , &
          & -1.5658d+1 ,  1.7160d+2 , -3.7212d+2 ,  4.1299d+2 /

!  j = 20;  N + N --> N + N + pi; nucleon distributions:
     data ((ankj(n,k,20),n=1,4),k=1,4) / &
          & 8.5591d-2 ,  5.0390d+0 , -1.3782d+1 ,  1.4661d+1 , &
          & 5.4284d-2 , -9.2324d+0 ,  3.6397d+1 , -4.2962d+1 , &
          & -5.1111d-2 ,  4.6003d+0 , -2.0534d+1 ,  2.7731d+1 , &
          & 7.4514d-3 , -6.2529d-1 ,  2.9159d+0 , -4.1101d+0 /

!  j = 21;  N + N --> N + N + pi; pion distributions:
     data ((ankj(n,k,21),n=1,4),k=1,4) / &
          & 7.1622d-2 ,  3.0960d+0 , -1.1125d+1 ,  1.8130d+1 , &
          & 9.2581d-2 , -3.2186d+0 ,  2.0273d+1 , -3.3245d+1 , &
          & -5.1531d-2 ,  8.9886d-1 , -7.5084d+0 ,  1.3188d+1 , &
          & 5.8258d-3 , -1.7288d-3 ,  7.0224d-1 , -1.4854d+0 /

!  j = 22;  N + N --> N + N + n*pi, n > 1; nucleon distributions:
     data ((ankj(n,k,22),n=1,4),k=1,4) / &
          & 8.2300d-2 ,  1.5854d-1 ,  3.7716d+0 , -4.0562d+0 , &
          & 1.0802d-2 , -3.3688d-1 ,  1.1727d+0 , -6.7476d-1 , &
          & -2.1798d-3 ,  5.2166d-2 , -2.5816d-1 ,  3.2048d-1 , &
          & 6.5764d-5 , -1.4711d-3 ,  7.8209d-3 , -1.0580d-2 /

!  j = 23;  N + N --> N + N + n*pi, n > 1; pion distributions:
     data ((ankj(n,k,23),n=1,4),k=1,4) / &
          & 1.1138d-1 ,  6.0396d-1 ,  3.0174d+0 , -4.4190d+0 , &
          & -1.7709d-2 ,  2.3015d-1 , -1.8187d+0 ,  3.4518d+0 , &
          & 2.0977d-3 , -2.5458d-2 ,  2.1626d-1 , -4.0692d-1 , &
          & -5.4799d-5 ,  5.9111d-4 , -5.5552d-3 ,  1.0647d-2 /

!  j = 24;  pi + N --> pi + N + pi; nucleon distributions:
     data ((ankj(n,k,24),n=1,4),k=1,4) / &
          & 1.7288d-1 ,  7.1080d+0 , -1.7961d+1 ,  1.6403d+1 , &
          & -1.4504d-1 , -1.3032d+1 ,  4.1781d+1 , -4.0799d+1 , &
          & 4.5390d-2 ,  8.3515d+0 , -3.0260d+1 ,  3.2882d+1 , &
          & -4.7961d-3 , -1.4095d+0 ,  5.3505d+0 , -6.0946d+0 /

!  j = 25;  pi + N --> pi + N + pi; pion distributions:
     data ((ankj(n,k,25),n=1,4),k=1,4) / &
          & 3.7596d-2 ,  1.4331d+0 , -3.1350d+0 ,  6.4864d+0 , &
          & 2.3827d-1 ,  1.8253d+0 ,  1.7648d+0 , -1.6735d+1 , &
          & -1.5410d-1 , -1.5201d+0 , -1.5692d+0 ,  1.7185d+1 , &
          & 2.5037d-2 ,  3.0588d-1 ,  3.2520d-1 , -3.5277d+0 /

!  j = 26;  pi + N --> pi + N + n*pi, n > 1; nucleon distributions:
     data ((ankj(n,k,26),n=1,4),k=1,4) / &
          & 1.2489d-1 ,  1.3573d+0 ,  8.2338d-1 , -1.4595d+0 , &
          & -5.1577d-2 , -3.5778d-1 , -1.1690d+0 ,  1.8078d+0 , &
          & 7.4864d-3 ,  3.2888d-2 ,  2.3744d-1 , -3.9802d-1 , &
          & -2.9880d-4 , -7.5117d-4 , -1.1402d-2 ,  1.9505d-2 /

!  j = 27;  pi + N --> pi + N + n*pi, n > 1; pion distributions:
     data ((ankj(n,k,27),n=1,4),k=1,4) / &
          & 1.8470d-1 ,  1.9269d+0 , -3.2979d+0 ,  3.6843d+0 , &
          & -7.3932d-2 ,  2.7213d-1 ,  1.0600d+0 , -2.3354d+0 , &
          & 1.8907d-2 , -5.6473d-2 , -1.6487d-1 ,  3.8426d-1 , &
          & -9.2984d-4 ,  2.5506d-3 ,  7.3052d-3 , -1.7220d-2 /

!  j = 28; gamma + N --> N + delta + pi; pion distribution; T < 1.0
     data ((ankj(n,k,28),n=1,4),k=1,4) / &
          & -1.0306d+0 ,  3.2849d+1 , -7.5052d+1 ,  6.0255d+1 , &
          & 7.9586d+0 , -1.2572d+2 ,  2.5604d+2 , -1.6547d+2 , &
          & -1.4797d+1 ,  1.6590d+2 , -2.7991d+2 ,  1.1333d+2 , &
          & 8.2309d+0 , -6.7871d+1 ,  8.5762d+1 ,  5.9727d+0 /

!  j = 29; gamma + N --> N + delta + pi; pion distribution; T > 1.0
     data ((ankj(n,k,29),n=1,4),k=1,4) / &
          & -2.3722d+2 ,  9.6890d+2 , -1.6219d+3 ,  1.3637d+3 , &
          & 6.5800d+2 , -2.6941d+3 ,  4.5480d+3 , -3.8460d+3 , &
          & -6.0653d+2 ,  2.4983d+3 , -4.2498d+3 ,  3.6136d+3 , &
          & 1.8604d+2 , -7.6933d+2 ,  1.3166d+3 , -1.1242d+3 /

!    (bnkj coefficients)
!  j = 1 ; nucleon final state momentum distributions,
!          N + N --> 2N + pi
     data ((bnkj(n,k, 1),n=1,4),k=1,4) / &
          & 5.0278d-1 ,  3.1442d+0 , -7.8172d+0 ,  8.1667d+0 , &
          & 9.3482d-1 , -1.0590d+1 ,  2.9227d+1 , -3.4550d+1 , &
          & -9.6685d-2 ,  4.7335d+0 , -1.4298d+1 ,  1.7685d+1 , &
          & -2.5041d-2 , -6.2478d-1 ,  2.0282d+0 , -2.5895d+0 /

!  j = 2 ; pion final state momentum distributions,
!          N + N --> 2N + pi
     data ((bnkj(n,k, 2),n=1,4),k=1,4) / &
          & 1.1965d+0 , -8.2889d-1 ,  1.0426d+0 , -1.9090d+0 , &
          & 2.8703d-1 , -4.9065d+0 ,  1.6264d+1 , -1.9904d+1 , &
          & -2.4492d-1 ,  2.9191d+0 , -9.5776d+0 ,  1.1938d+1 , &
          & 3.7297d-2 , -4.2200d-1 ,  1.3883d+0 , -1.7476d+0 /

!  j = 3 ; nucleon final state momentum distributions,
!          N + N --> 2N + n*pi, n > 1
     data ((bnkj(n,k, 3),n=1,4),k=1,4) / &
          & 1.3508d+0 , -4.3139d+0 ,  1.2291d+1 , -1.5288d+1 , &
          & -2.0086d-1 ,  1.3641d+0 , -3.4030d+0 ,  3.8559d+0 , &
          & 1.2583d-2 , -8.3492d-2 ,  1.8600d-1 , -2.0043d-1 , &
          & -2.3628d-4 ,  1.3514d-3 , -2.4324d-3 ,  2.1906d-3 /

!  j = 4 ; pion final state momentum distributions,
!          N + N --> 2N + n*pi, n > 1
     data ((bnkj(n,k, 4),n=1,4),k=1,4) / &
          & 1.2419d+0 , -4.3633d+0 ,  1.3743d+1 , -1.8592d+1 , &
          & -2.4404d-1 ,  1.3158d+0 , -3.5691d+0 ,  4.3867d+0 , &
          & 1.5693d-2 , -8.2579d-2 ,  2.1427d-1 , -2.5846d-1 , &
          & -2.9386d-4 ,  1.4060d-3 , -3.3835d-3 ,  3.8664d-3 /

!  j = 5 ; nucleon final state momentum distributions,
!          pi + N --> N + 2*pi
     data ((bnkj(n,k, 5),n=1,4),k=1,4) / &
          & 6.3054d-1 , -3.7333d+0 ,  1.3464d+1 , -1.8594d+1 , &
          & 2.1801d+0 ,  1.5163d+0 , -1.6380d+1 ,  2.7944d+1 , &
          & -1.2886d+0 , -2.4570d+0 ,  1.5129d+1 , -2.3295d+1 , &
          & 2.0915d-1 ,  5.2279d-1 , -2.8687d+0 ,  4.2688d+0 /

!  j = 6 ; pion final state momentum distributions,
!          pi + N --> N + 2*pi
     data ((bnkj(n,k, 6),n=1,4),k=1,4) / &
          & 9.3363d-1 , -1.8181d+0 ,  5.5157d+0 , -8.5216d+0 , &
          & 1.7811d+0 , -8.2927d+0 ,  2.0607d+1 , -2.0827d+1 , &
          & -1.5264d+0 ,  6.8433d+0 , -1.6067d+1 ,  1.6845d+1 , &
          & 2.7128d-1 , -1.1944d+0 ,  2.7495d+0 , -2.9045d+0 /

!  j = 7 ; nucleon final state momentum distributions,
!          pi + N --> N + n*pi, n > 2
     data ((bnkj(n,k, 7),n=1,4),k=1,4) / &
          & 1.9439d+0 , -4.6268d+0 ,  9.7879d+0 , -9.6074d+0 , &
          & -3.4640d-1 ,  1.1093d+0 , -1.9313d+0 ,  1.7064d+0 , &
          & 2.7054d-2 , -1.1638d-1 ,  2.6969d-1 , -3.1853d-1 , &
          & -6.6092d-4 ,  5.0728d-3 , -1.4995d-2 ,  1.9605d-2 /

!  j = 8 ; pion final state momentum distributions,
!          pi + N --> N + n*pi, n > 2
     data ((bnkj(n,k, 8),n=1,4),k=1,4) / &
          & 1.8693d+0 , -5.5678d+0 ,  1.4795d+1 , -1.6903d+1 , &
          & -4.9965d-1 ,  1.7874d+0 , -4.1330d+0 ,  3.8393d+0 , &
          & 4.6194d-2 , -1.8536d-1 ,  4.5315d-1 , -4.6273d-1 , &
          & -1.3341d-3 ,  5.7710d-3 , -1.4554d-2 ,  1.5554d-2 /
!
!  Coefficients c_k(J) below are from Tab. 76, and J denotes the type of
!  interaction as following:
!
!    J=
!    1    N  + N --> 2N + pi          (nucleons)
!    2    N  + N --> 2N + pi          (pions)
!    3    N  + N --> 2N + n*pi  n>1   (nucleons)
!    4    N  + N --> 2N + n*pi  n>1   (pions)
!    5    pi + N -->  N + 2*pi        (nucleons)
!    6    pi + N -->  N + 2*pi        (pions)
!    7    pi + N -->  N + n*pi  n>2   (nucleons)
!    8    pi + N -->  N + n*pi  n>2   (pions)
!
     data (ckj(k, 1),k=1,3) / 1.4509d-1 , 4.6520d-1 , -3.3005d-2 /
     data (ckj(k, 2),k=1,3) / 1.5376d-1 , 2.7436d-1 , -1.4604d-2 /
     data (ckj(k, 3),k=1,3) / 6.2959d-1 , 1.7866d-1 , -2.6216d-3 /
     data (ckj(k, 4),k=1,3) / 8.3810d-1 , 8.6137d-3 ,  3.2946d-3 /
     data (ckj(k, 5),k=1,3) / 9.2852d-2 , 5.3886d-1 , -5.4493d-2 /
     data (ckj(k, 6),k=1,3) / 1.3032d-1 , 4.0709d-1 , -2.8782d-2 /
     data (ckj(k, 7),k=1,3) / 1.4909d-1 , 3.8502d-1 , -1.2775d-2 /
     data (ckj(k, 8),k=1,3) / 1.8024d-1 , 3.3022d-1 , -9.4491d-3 /
! Comments:
!
! (**)     At lower energies, angular distributions may be considered
!       approximatively as isotropic. Deviation from isotropy for angles
!       theta near 0 and 180, as a result of Coulomb effect, only
!       weakly affects results of cascade calculations.
!       It sould be noted that for a symmetric system p+p, the values
!       of a_nk in these tables correspond to the interval
!       0 < cos(Theta) < 1 only.
!       In Eq. (4.42) one should to use function w(z) instead of
!       w(2z-1), and cos(Theta) has to be calculated in this case not
!       under Eq. (4.44) but using the first formula on p. 296, in
!       Comments, after Tab. 72.
!
! (***)    Coefficients a_nk for this energy-region correspond to only
!       0.5 < cos(Theta) < 1. In this case, cos(Theta) has to be
!       calculated not using Eg. (4.44) but using the second formula
!       on p. 296, in Comments, after Tab. 72.
!
! (+)   For higher energies, the elastic scattering become almost
!       diffractive and the corresponding angular distributions may be
!       approximated with a good range of accuracy by the exponential
!       function (4.46).
!
! (++)  For higher energies, angular distributions may be taken the
!       same as for pp-scattering.
!
! (+++) For higher energies, angular distributions coincide with the
!       ones of pi-p elastic scattering.
!
! (1) For regions where |cos(Theta)| is approximatively equal to 1,
!     Eq. (4.44) may give rise to  |cos(Theta)| > 1. Is such cases
!     one should use |cos(Theta)|=1.
!
! (2) The coefficients for elastic n+n, pi(+,-)+n, pi0+n(p), and for CEX
!     pi0+n -> pi- + p scattering are defined using the charge symmetry
!     law.
!
! ======================================================================
