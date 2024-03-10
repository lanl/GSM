
  subroutine setupMDCMReaction()

! ====================================================================
!
!  LAQGSM Initial set module,  K. GUDIMA 04/27/07
!
!  This module is used for standalone LAQGSM calculations
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMData, only: mDCMDataInitialized, initializeModifiedDCMData

    implicit none
    integer(int32) :: ia1, ia2, k
    real(real64)   :: am, dlam, p0, temp1, temp2, temp3, temp4, &
         & tempA, tempZ, z0

! ====================================================================

    real(real64),   parameter, dimension(10) :: rms = &
         [ 0.85,2.095,1.976,1.671,2.50,2.57,2.45,2.519,2.45,2.42 ]

! ====================================================================

    ! Initial particles:
    real(real64) :: anucl1, anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    real(real64) :: stin, amin
    common /stin/ stin,amin

    ! For calculation (?):
    integer(int32) :: kobr
    real(real64) :: blab, glab
    common/bglab/blab,glab,kobr   ! Anti-lab system
    real(real64) :: xap
    integer(int32) :: inside, ivalon
    common/kappa/xap
    common/cinsid/inside
    common/cvalon/ivalon
    real(real64) :: cmali
    common/cmali/cmali
    real(real64) :: rcor
    common/rcor/ rcor
    real(real64) :: rint
    common/rint/ rint
    real(real64) :: pidabs
    common/pidabs/pidabs
    integer(int32) :: itmgo, min
    real(real64) :: timgo
    common/timgo/timgo,min,itmgo
    real(real64) :: mmes, mbar
    common/mmatur/mmes,mbar
    integer(int32) :: indint
    common/indint/indint
    integer(int32) :: intcc
    common/intcc/intcc
    real(real64) :: tlimit
    common/tlimit/tlimit
    integer(int32) :: iact
    common /iact/ iact
    integer(int32) :: inddec
    common/inddec/inddec
    ! isob3 looks like an option for mass computations, perhaps to signify
    ! baryon simulations?
    integer(int32) :: isob3
    common/isob3/isob3
    integer(int32) :: isob2
    common/isob2/isob2
    integer(int32) :: ifib0
    real(real64) :: xbmax
    common/xbmax/xbmax,ifib0
    integer(int32) :: ksyst
    common/ksyst/ksyst
    real(real64) :: bnn, gnn
    common/cmnn/bnn,gnn
    integer(int32) :: ncas, ncpri
    common/ncasca/ncas,ncpri
    integer(int32) :: iw10
    common /iw10/ iw10
    real(real64) :: ares1, ares2, collt
    common/ncollt/ ares1,ares2,collt(4)
    integer(int32) :: iran
    common/bran/ iran
    real(real64) :: rn, delta
    common/rndelt/ rn,delta
    integer(int32) :: ifspe, ncas1
    common/ifspe/ifspe,ncas1
    real(real64) :: sig1, sig2
    common/geocrs/ sig1,sig2
    integer(int32) :: isys
    real(real64) :: t0a, vla, gla, vev, gev, vcm, gcm
    common/vgsys/t0a,vla,gla,vev,gev,vcm,gcm,isys
! kkg 04/17/05
    real(real64) :: pud, sigma, alfa, beta
    common/comind/ pud,sigma,alfa,beta
    real(real64) :: ecmb
    common/comecb/ ecmb
    real(real64) :: enbou
    common /comenb/ enbou

    ! ------------------
    ! IN GSM AND MDCM:
    ! ------------------
    integer(int32) :: ibrems
    real(real64) :: tgmin, tgmax, teqv, sxabs
    common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems

! ====================================================================

    ! Read in data files and set various constants:
    if ( .not. mDCMDataInitialized ) then
       call initializeModifiedDCMData()
    end if

    iran=12345
!    CALL RDMINI


    ares1=0.
    ares2=0.
    do k=1,4
       collt(k)=0.
    enddo

    ! KKG  04/17/07
    ecmb  = 10.0
    enbou = 2.000
    pud=0.415
    sigma=0.51 ! 05.05.06

    tlimit=100.
    timgo =9000.
    xbmax =1.0
    rcor  =0.0
    rint  =0.0
    rn    =0.2
    delta =1.3
    pidabs=0.01
!    MMES  =0.1
    mmes  =0.7
    mbar  =0.0
    ifib0 =1

    ! Check for anti-lab boosting:
    ! KKG 07/13/04
    if(anucl1 > anucl2)  then
       ! Swap A1 with A2, and Z1 with Z2
       tempa = anucl1
       anucl1 = anucl2
       anucl2 = tempa
       tempz = znucl1
       znucl1 = znucl2
       znucl2 = tempz

       ! Boost results into antilab. system
       kobr  = 1_int32
       write(*,1200)
    else
       kobr  = 0_int32
    endif

    indint=1
    intcc =2
    ksyst =1 !  lab.  syst
!    KSYST =2 !  Eq.V. syst.
    iact  =2
    inddec=0
    isob2 =1
    isob3 =1

! For diagnostic printing (prints after NCAS > NCPRI)
    ncas = 0
    ncpri =999999999
!    NCPRI =5760
!    read(15,*)  NCPRI

    iw10  =0
    ifspe =1
    xap   =1.0
    cmali =1.1
    inside=0
    ivalon=0
    eps1  =0.007
    eps2  =0.007
    vpi   =0.025
    c1    =0.545
    c2    =0.545
    d1    =0.05
    d2    =0.05
    r0n1  =1.07
    r0n2  =1.07
    if(anucl1 > 2.1.and.anucl1 < 10.1)  then
       ia1=int(anucl1+0.1)
       r0n1=rms(ia1)
    endif
    if(anucl2 > 2.1.and.anucl2 < 10.1)  then
       ia2=int(anucl2+0.1)
       r0n2=rms(ia2)
    endif

    ! for D/p ratio
    if(anucl2 < 2.1)  d2=0.001

    ! Prints photon type to output file
    if(amin < 0.0001.and.ibrems == 0)  then
       ! Mono. Energetic Photon
       write(16,275) t0
    elseif(amin < 0.0001.and.ibrems == 1)  then
       ! Bremss. Photon
       write(16,276) tgmin,tgmax
    endif

    temp1 = t0+1.88
    ! call err_chk(1,'Laqi.f','206',1,temp1)
    temp2 = t0/temp1
    ! call err_chk(1,'Laqi.f','209',2,temp2)
    bnn=sqrt(temp2)

    temp1 = 1.+t0/1.88
    ! call err_chk(1,'Laqi.f','212',2,temp1)
    gnn=sqrt(temp1)

    ! Setup nuclei properties:
    rm1=0.
    rm2=0.
    if(anucl1 > 1.1)  call helpq(r0n1,anucl1,a1,c1,d1,tf01,rm1)
    if(anucl2 > 1.1)  call helpq(r0n2,anucl2,a2,c2,d2,tf02,rm2)
    if(anucl1 == 2.0)  rm1=2.158 ! kkg 26.11.01
    if(anucl2 == 2.0)  rm2=2.158 ! kkg 26.11.01
    am=0.940
    if(anucl1 <= 1.1)  am=amin
    temp1 = t0*(t0+2.*am)
    ! call err_chk(1,'Laqi.f','223',2,temp1)
    temp2 = 5.06*sqrt(temp1)
    ! call err_chk(1,'Laqi.f','226',1,temp2)

    dlam=1./temp2
    z0=rm1+rm2+dlam+delta
    sig1=31.459*(z0**2)*xbmax**2
    sig2=31.459*((rm1+rm2)**2)*xbmax**2

    temp1 = t0*(t0+2.*am)
    ! call err_chk(1,'Laqi.f','236',2,temp1)
    p0=sqrt(temp1)
    temp1 = t0+am
    ! call err_chk(1,'Laqi.f','241',1,temp1)
    blab=p0/temp1
    vla=blab

    if(abs(anucl1) < 0.001)  then
       temp1 = (t0+anucl2*0.940)
       ! call err_chk(1,'Laqi.f','244',1,temp1)
       vcm=p0/temp1
    else
       temp1 = (t0+(1.+anucl2/abs(anucl1))*0.94)
       ! call err_chk(1,'Laqi.f','248',1,temp1)
       vcm=p0/temp1
    endif

    !  kkg  12/10/04
    if(am <= 0.0001) then
       glab=0.0
       gla =0.0
    else
       temp1 = 1.-blab**2
       ! call err_chk(1,'Laqi.f','257',2,temp1)
       temp2 = sqrt(temp1)
       ! call err_chk(1,'Laqi.f','259',1,temp2)
       glab=1./temp2
       gla =glab
    endif
    t0a=t0

    temp1 = t0+1.88
    ! call err_chk(1,'Laqi.f','266',1,temp1)
    vev=p0/temp1

    temp1 = 1.-vcm**2
    ! call err_chk(1,'Laqi.f','270',2,temp1)
    temp2 = sqrt(temp1)
    ! call err_chk(1,'Laqi.f','272',1,temp2)
    gcm=1./temp2

    temp3 = 1.-vev**2
    ! call err_chk(1,'Laqi.f','276',2,temp3)
    temp4 = sqrt(temp3)
    ! call err_chk(1,'Laqi.f','278',1,temp4)
    gev=1./temp4

    ! call egsfill ! Not used in GSM calculations at all


    ! Establish photon information:
    if ( anucl1 < 0.1 .and. amin < 0.001 ) then
       call inigam0(t0)
       call inigamn(t0)
    end if

    return

! ====================================================================
275 format(/2x,'gamma + A interaction at E_g=',f7.3,' GeV'/)
276 format(/2x,'bremsstrahlung gamma(Egmin=',f7.3,', Egmax=',f7.3, &
         & ') + A'/)
1200 format(3x, 'Boosting to the anti-lab system.')
! ====================================================================
  end subroutine setupMDCMReaction
