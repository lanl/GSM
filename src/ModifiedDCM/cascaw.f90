
  subroutine cascaw(nel,rn,delta,iret)

! ====================================================================
!
!     Modified: 19-DEC-2003 BY NVM
!     Edited by KKG, May 2007 to include GDR region fot gamma + A
!     interaction below 20 MeV
!     Modified: 09-MAY-2007 BY NVM
!     Modified: 08/24/16 by CMJ to allow GSM use
!     Modified: 05/2018 by CMJ, XCP-3, to correct tallying of fragment
!          data
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMClass, only: results
    use standarddcmdata, only: photonEG

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    integer(int32), intent(  out) :: nel
    real(real64),   intent(in   ) :: rn
    real(real64),   intent(in   ) :: delta
    integer(int32), intent(  out) :: iret

    integer(int32), dimension(5) :: ip1, ip2, ip3
    real(real64),   dimension(3) :: v12, ptemp
    real(real64),   dimension(9) :: p1, p2, p3

    integer(int32) :: mv

! ====================================================================

    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre, &
         & vep(3),vet(3),gep,get
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/memorylaq/pmemo(9,5999),imemo(5,5999) &
         & /kypt/ kpt,kyp,kyt,kpt1,kyp1,kyt1,kyy,kyy1
    common/jgstar/ jg
    common/nrapi/nrapi
    common/nucoll/ nucoll,mvcoll(5999)
    common/actim/tint
    common/intcc/intcc
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/ipaul/ip
    common/timgo/timgo,minute,itmgo
    common/iprimc/iprimc
    common/idpme/idpme(5999)
    common/mvup/mvu
    common/bimp/bimp,bimx,bimy
    common/ncollt/ ares1,ares2,collt(4)
    common /stin/ stin,amin
    common/geocrs/ sig1,sig2
!  kkg 10/14/03
    common/inttyp/ ityp
    common /ncasca/ ncas,ncpri
    common/bglab/blab,glab,kobr
    common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems

! ====================================================================

    parameter (zro=0.0d0)

! ====================================================================

    nel=0
    iret=0

9   mvu=0
    nucoll=0
    an1=anucl1
    zn1=znucl1
    an2=anucl2
    zn2=znucl2
    enext1=zro
    enext2=zro
    pnucl1(1)=zro
    pnucl1(2)=zro
    pnucl1(3)=zro
    pnucl2(1)=zro
    pnucl2(2)=zro
    pnucl2(3)=zro
    amnuc1(1)=zro
    amnuc1(2)=zro
    amnuc1(3)=zro
    amnuc2(1)=zro
    amnuc2(2)=zro
    amnuc2(3)=zro
    tint=zro
    mv=0
    kpt=0
    kyp=0
    kyt=0
    kyy=0
    kpt1=0
    kyp1=0
    kyt1=0
    kyy1=0
    results%projExc%numTotal=0
    results%targExc%numTotal=zro
    results%projExc%numProtons=zro
    results%targExc%numProtons=zro
    results%projExc%numHoles=zro
    results%targExc%numHoles=zro
    do l=1,300
       mpa(l)=1
    end do

    call pinpnq(rn,r0x,r0y,r0z,mv,delta,nel)
    obr1=zro
    if(rm1 > 0.1d0) obr1=0.00144d0*zn1/rm1
    obr2=zro
    if(rm2 > 0.1d0) obr2=0.00144d0*zn2/rm2
    call  vinit(vpr,vta,anucl1,anucl2,t0)
    radp(1)=r0x
    radp(2)=r0y
    radp(3)=r0z
    radt(1)=zro
    radt(2)=zro
    radt(3)=zro

!  extension  to Egamma < 20 MeV
    if(amin < 0.0001d0.and.t0 < 0.020d0)  then
       t0mev = t0*1000.d0
       sgabs = photoneg%photocrosssection(t0mev, anucl2)
       pint = sgabs/sig1
       if(rndm(-1.0_real64) > pint)  then
          go  to  21
       else
          an2=anucl2
          zn2=znucl2
          enext2 = t0
          pnucl2(1)=zro
          pnucl2(2)=zro
          pnucl2(3)= t0
          amnuc2(1)=1.0d0
          amnuc2(2)=zro
          amnuc2(3)=zro
          results%targExc%numHoles = 1.0d0
          results%targExc%numTotal = 1.0d0
          results%targExc%numProtons = zro
          if(rndm(-1.0_real64) <= (znucl2/anucl2)) results%targExc%numProtons = 1.0d0
          if(ibrems == 1)  then
             teqv = teqv + t0
             sxabs = sxabs + sgabs
          endif
          return
       endif
    endif
!
    ip=-1
    ijpa=0
11  na1=an1+0.1d0
    na2=an2+0.1d0
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(na1 >= 1.and.na2 >= 1)  then    !
       do  l=1,na1                     !
          mpa(l)=1                      !
       enddo                           !
    endif
!   kkg 11/05/03                             !
    if(mv >= 1.and.nucoll.ne.0)   then ! 05.09.97
       do  l=1,mv                      !
          if(na1 >= 1)  myp(l)=1        !
          if(na2 >= 1)  myt(l)=1        !
          if(iabs(idpme(l)) <= 20) then !
             myp(l)=0        !
             myt(l)=0        !
             myy(l)=0        !
          endif                         !
       enddo                           !
    endif                              !
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iprimc=0
! NVM:      call  stopgu(timgo,istop)
! NVM:      if(istop.ne.0)   iret=1
! NVM:      if(istop.ne.0)   return
!      CALL  UPACOW(MV,MVU,0)
    call  upacow(mv,mvu,0)
    call  dissip(ip,t0)
!
! if(NCAS > NCPRI) then
!       write(16,*) ' NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2=',
!    &                NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2
!       write( *,*) ' NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2=',
!    &                NCAS,AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2
!       endif
!
    call pointn(m,p1,p2,ip1,ip2,v12,u12,t12,sig,sab,tint,nr,n1,n2, &
         & na1,na2,mv,delta,r0x,r0y,r0z)
    ntry=0
    if(m)   203,203,12
203 if(mv == 0)   go  to  205
!----> do  204  l=1,mv
    do   l=1,mv
       myp(l)=0
       myt(l)=0
       myy(l)=0
204    continue
    end do
205 continue
!  kkg 17.06.02
!     IF(MV == 0.and.AN2 == ANUCL2.and.ZN2 == ZNUCL2.and.
!    &   ENEXT2 < 0.001)  GO  TO  21
    if(enext2 < 1.d-7)  go  to  21
!      CALL  UPACOW(MV,MVU,1)
    call  upacow(mv,mvu,1)
    mv=mvu
    go  to  20
12  go  to  (13,14,13,112,15),nr
13  nu=2
    go to 15
14  nu=1
15  continue
    ntry=ntry+1
!  kkg  09/05/08
    call typnew(p1,ip1,p2,ip2,v12,u12,t12,sab,mv,np,nabs, &
         & p3,ip3,n3,nu,n2,na1,na2,iret)
    if(iret == 1)   go  to  11
!
    nrapi=nr
!
!
    ip=0
!  kkg 10/14/03
    if(ityp == 0)  go  to  11
!
    if(np)11,11,16
16  go to  (17,18,19,112,113),nr
17  continue
!
    call pauli1(p1,p2,ip1,ip2,n1,n2,v12,np,mv,tint,ip, results)
!
    kpt1=kpt1+ip
    kpt=kpt+1
    go  to  11
18  continue
!
    call pauli2(p1,p2,p3,ip1,ip2,ip3,n1,n2,n3,v12,np,mv,tint, &
         & ip,obr1, results)
!
    kyp1=kyp1+ip
    kyp=kyp+1
    go  to  11
19  continue
    if(anucl1 < 1.1d0.and.nucoll <= 1) iprimc=1
!
    call pauli3(p1,p2,p3,ip1,ip2,ip3,n1,n2,n3,v12,np,mv,tint, &
         & ip,obr2, results)
!
    kyt1=kyt1+ip
    kyt=kyt+1
    go  to  11
112 call  decren(n1,nur,mv,np)
!----> do 213 l=1,9
    do l=1,9
       p1(l)=pmemo(l,n1)
       if(l <= 5)  ip1(l)=imemo(l,n1)
213    continue
    end do
    call  paulid(mv,n1,nur,1,np,ip)
    go  to  11
113 continue
!
    call pauli4(p1,p2,ip1,ip2,n1,n2,v12,np,mv,tint,ip)
!
    kyy1=kyy1+ip
    kyy=kyy+1
!
    go  to  11
!  kkg 17.06.02
20  continue
!  kkg  18/04/03
    if(anucl1 >= 0.9d0)  then
       exmax2=anucl1*(t0+eps2)
    else
       exmax2=t0+amin
    endif
    exmin=0.0001d0
    if((an1 < zn1.or.zn1 < zro.or.enext1 < zro).or. &
         & (an2 < zn2.or.zn2 < zro.or.enext2 < exmin).or. &
         & (enext2 > exmax2)) then
!         write(*,*) 'AN1,ZN1,AN2,ZN2,ENEXT1,ENEXT2=',
!     &               AN1,ZN1,AN2,ZN2,ENEXT1,ENEXT2
       iret=1
       return
    else
       go  to  22
    endif
!  kkg  17.06.02
21  nel=nel+1
    go to 9
!----> do 23 l=1,3
22  do l=1,3
       amnuc1(l)=amnuc1(l)*5.06d0
23     amnuc2(l)=amnuc2(l)*5.06d0
    end do
!   kkg 12/10/04
    if(amin < 0.0001d0.and.ibrems == 1)  then   ! bremss. gamma
       teqv = teqv + t0
    endif
    if(anucl1-2.1d0)  24,26,26
26  am1=.94d0*an1
    call cinema(pnucl1,vpr,ptemp,ct,st,cf,sf,ttemp,am1)
    pnucl1(1)=ptemp(1)
    pnucl1(2)=ptemp(2)
    pnucl1(3)=ptemp(3)
!  kkg 07/13/04
    if(kobr == 1)  then
       pnucl1(3)=-glab*(pnucl1(3)-blab*(ttemp+am1))
    endif
24  continue
    if(anucl2 < 2.1d0)   go  to  27
    am2=.94d0*an2
    call cinema(pnucl2,vta,ptemp,ct,st,cf,sf,ttemp,am2)
    pnucl2(1)=ptemp(1)
    pnucl2(2)=ptemp(2)
    pnucl2(3)=ptemp(3)
!  kkg 07/13/04
    if(kobr == 1)  then
       pnucl2(3)=-glab*(pnucl2(3)-blab*(ttemp+am2))
    endif
27  continue
    collt(1)=collt(1)+kpt1
    collt(2)=collt(2)+kyp1
    collt(3)=collt(3)+kyt1
    collt(4)=collt(4)+kyy1
    ares1=ares1+an1
    ares2=ares2+an2


    if(ncas > ncpri) then
       write(*,*) ' cascaw: ncas,an1,an2,zn1,zn2,enext1,enext2=', &
            & ncas,an1,an2,zn1,zn2,enext1,enext2
       return
    endif

    return
! ====================================================================
  end subroutine cascaw
