
  subroutine rapid1(na1,na2,delta,p1,ip1,p2,ip2,n1,n2,tau,dl1)

! =====================================================================
!
!  Determination of interaction  pair (N1,P1,IP1),(N2,P2,IP2)
!  where nucleon N1 is from  projetile and
!        nucleon N2 is from  target; the time TAU (fm/c) is calculated
!  in observer's system (lab. or equal velocity system depending of
!  projectile VPR and target VTA nuclei velocities)
!  Calls: PARTNQ, TFERMIQ, KINEMQ
!  Input: NA1, NA2 - mass numbers of projectile and target nuclei
!         DL1 - interaction parameter
!  Output: P1,IP1,P2,IP2,N1,N2,TAU
!
! =====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    integer(int32), intent(in   ) :: na1
    integer(int32), intent(in   ) :: na2
    real(real64),   intent(in   ) :: delta
    real(real64),   intent(inout) :: p1(9)
    integer(int32), intent(inout) :: ip1(5)
    real(real64),   intent(inout) :: p2(9)
    integer(int32), intent(inout) :: ip2(5)
    integer(int32), intent(  out) :: n1
    integer(int32), intent(  out) :: n2
    real(real64),   intent(  out) :: tau
    real(real64),   intent(  out) :: dl1

    real(real64), dimension(3) :: ri, rj, rij, vij, pi, pj, r1, r2, pli, plj

! =====================================================================

    common/rcor/ rcor
    common/ncasca/ncas,ncpri
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common /holpt/ ph1(3),ph2(3),rh1(3),rh2(3),ehol1,ehol2,tfp,tft
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get

! =====================================================================

    temp1 = gre*gre-1.
!      call err_chk(1,'LAQ1.f','477',2,temp1)
    temp2 = 5.06*0.940*sqrt(temp1)
!      call err_chk(1,'LAQ1.f','480',1,temp2)

    dl1=1./temp2+delta

    temp3 = abs(1.-vpr(1)**2-vpr(2)**2-vpr(3)**2)
    temp4 = abs(1.-vta(1)**2-vta(2)**2-vta(3)**2)
!      call err_chk(1,'LAQ1.f','488',2,temp3)
!      call err_chk(1,'LAQ1.f','489',2,temp4)


    gpr=1./sqrt(temp3)
    gta=1./sqrt(temp4)
!----> do   8  k=1,3
    do    k=1,3
       pi(k)=0.940*gpr*vpr(k)
       pj(k)=0.940*gta*vta(k)
       vij(k)=vpr(k)-vta(k)
8      continue
    end do
    vij2=vij(1)**2+vij(2)**2+vij(3)**2
    ei=0.940*gpr
    ej=0.940*gta
    sij=(ei+ej)**2-(pi(1)+pj(1))**2-(pi(2)+pj(2))**2- &
         & (pi(3)+pj(3))**2
    tau=-0.1

    temp1 = gpr+1.
    temp2 = gta+1.
!      call err_chk(1,'LAQ1.f','510, 511, 512',1,temp1)
!      call err_chk(1,'LAQ1.f','516, 517, 518',1,temp2)

!----> do  13  i=1,na1
    do   i=1,na1
       if(mpa(i) == 0)   go  to  13
       ni=0
       vri=xc(1,i)*vpr(1)+yc(1,i)*vpr(2)+zc(1,i)*vpr(3)
       ri(1)=xc(1,i)-vpr(1)*vri*gpr/(temp1)+radp(1)
       ri(2)=yc(1,i)-vpr(2)*vri*gpr/(temp1)+radp(2)
       ri(3)=zc(1,i)-vpr(3)*vri*gpr/(temp1)+radp(3)
!----> do  12  j=1,na2
       do   j=1,na2
          vrj=xc(2,j)*vta(1)+yc(2,j)*vta(2)+zc(2,j)*vta(3)
          rj(1)=xc(2,j)-vta(1)*vrj*gta/(temp2)+radt(1)
          rj(2)=yc(2,j)-vta(2)*vrj*gta/(temp2)+radt(2)
          rj(3)=zc(2,j)-vta(3)*vrj*gta/(temp2)+radt(3)
!----> do  9  k=1,3
          do   k=1,3
             rij(k)=ri(k)-rj(k)
9            continue
          end do

          temp1 = vij2
!      call err_chk(1,'LAQ1.f','526',1,temp1)

          tij=-(rij(1)*vij(1)+rij(2)*vij(2)+rij(3)*vij(3))/temp1
          if(tij <= 0.)      go  to  12
          sp=rij(1)*pi(1)+rij(2)*pi(2)+rij(3)*pi(3)
          st=rij(1)*pj(1)+rij(2)*pj(2)+rij(3)*pj(3)
          rij2=rij(1)**2+rij(2)**2+rij(3)**2

          temp1 = sij-4.*0.940**2
          temp2 = sij
!      call err_chk(1,'LAQ1.f','537',1,temp1)
!      call err_chk(1,'LAQ1.f','537',1,temp2)

          b2=rij2+4.*(sij*sp*st-(0.940*(sp+st))**2)/temp2/temp1
          if(b2 > (dl1**2))          go  to  12
          if(ipo == 1.and.nrst == 1.and.i == ist.and.j == jst)  go  to 12
          ni=ni+1
          if(tau < 0.)        go  to  10
          if(tau < tij)       go  to  12
10        tau=tij
          n1=i
          n2=j
!----> do  11  k=1,3
          do   k=1,3
             r1(k)=ri(k)+tij*vpr(k)
             r2(k)=rj(k)+tij*vta(k)
11           continue
          end do
12        continue
       end do
       if(ni == 0)   mpa(i)=0
13     continue
    end do
    if(tau <= 0)  return
    call  partnq(1,n1,p1,ip1)
    tfp=tfermiq(p1(1),p1(2),p1(3),1)
    ehol1=p1(8)
    call  partnq(2,n2,p2,ip2)
    tft=tfermiq(p2(1),p2(2),p2(3),2)
    ehol2=p2(8)
!----> do  14  k=1,3
    do   k=1,3
       pi(k)=p1(3+k)
       pj(k)=p2(3+k)
14     continue
    end do
    call  kinemq(pi,vpr,pli,cti,sti,cfi,sfi,tli,p1(9))
    call  kinemq(pj,vta,plj,ctj,stj,cfj,sfj,tlj,p2(9))
!----> do  15  k=1,3
    do   k=1,3
       rh1(k)=p1(k)
       rh2(k)=p2(k)
       p1(k) =r1(k)
       p2(k) =r2(k)
       ph1(k)=-p1(3+k)
       ph2(k)=-p2(3+k)
       p1(3+k)=pli(k)
       p2(3+k)=plj(k)
15     continue
    end do
    p1(7)=0.
    p2(7)=0.
    p1(8)=sqrt(pli(1)**2+pli(2)**2+pli(3)**2+0.940**2)-0.940
    p2(8)=sqrt(plj(1)**2+plj(2)**2+plj(3)**2+0.940**2)-0.940
    p1(9)=0.940
    p2(9)=0.940
    return
! =====================================================================
  end subroutine rapid1


  subroutine rapid2(mv,na1,delta,p1,ip1,n1,n2,tau,dl)

! =====================================================================
!
!  Determination of interaction  pair (N1,P1,IP1),N2,
!  where nucleon N1 is number of cascade particle in memory array and
!  N2 is number of spectator nucleon from projectile nucleus;
!  the time interval TAU (fm/c) is calculated in observer's system
!  (lab. or equal velocity system in dependence of
!  projectile VPR and target VTA nuclei velocities)
!  Calls: CENUM1, CINEMA
!  Input: MV - total number of cascade particles,
!         NA1 -mass number of projectile nucleus
!         DL - interaction parameter
!  Output: P1,IP1,N1,N2,TAU
!
! =====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)(a-h,o-z), integer(int32) (i-n)
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(in   ) :: na1
    real(real64),   intent(in   ) :: delta
    real(real64),   intent(inout) :: p1(9)
    integer(int32), intent(inout) :: ip1(5)
    integer(int32), intent(  out) :: n1
    integer(int32), intent(  out) :: n2
    real(real64),   intent(  out) :: tau
    real(real64),   intent(  out) :: dl

    integer(int32), dimension(5) :: iin
    real(real64),   dimension(3) :: v0, cs, c
    real(real64),   dimension(9) :: pin

! =====================================================================

    common/rcor/ rcor
    common/ncasca/ncas,ncpri
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/cenpar/nukc(100)
    common/taue/tpte,type,tyte
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/cslid/clider(5999)
    common/iact/ iact/cvalon/ ivalon

! =====================================================================

    k = 1
    tau = -.1
    taul=-.1
    v0(1)=-vpr(1)
    v0(2)=-vpr(2)
    v0(3)=-vpr(3)

    temp1 = 1.-vpr(1)**2-vpr(2)**2-vpr(3)**2
!      call err_chk(1,'LAQ1.f','626',2,temp1)
    temp2 = sqrt(temp1)
!      call err_chk(1,'LAQ1.f','629',1,temp2)

    gpr=1./temp2
    gg=gpr*gpr/(gpr+1.)
9   if(myp(k)) 26,26,10
10  vrk=(pmemo(1,k)-radp(1))*vpr(1)+(pmemo(2,k)-radp(2))*vpr(2)+ &
         & (pmemo(3,k)-radp(3))*vpr(3)
    xk0=pmemo(1,k)-radp(1)+vpr(1)*vrk*gg
    yk0=pmemo(2,k)-radp(2)+vpr(2)*vrk*gg
    zk0=pmemo(3,k)-radp(3)+vpr(3)*vrk*gg
    cs(1)=pmemo(4,k)
    cs(2)=pmemo(5,k)
    cs(3)=pmemo(6,k)
    call cinema(cs,v0,c,ct,st,cf,sf,tl,pmemo(9,k))
    pin(4)=c(1) ! x momentum
    pin(5)=c(2) ! y momentum
    pin(6)=c(3) ! z momentum
    pin(7)=pmemo(7,k) ! cos(phi)
    pin(9)=pmemo(9,k) ! rest mass (gev)
    iin(1)=imemo(1,k) ! particle charge
    iin(2)=imemo(2,k) ! =1 for gammas, 0 otherwise
    iin(3)=imemo(3,k) ! strangeness number for particle
    iin(4)=imemo(4,k) ! baryon number
    iin(5)=imemo(5,k) ! zone number of nucleus
    clid=clider(k)
    pin(1)=xk0 ! x position at reaction start for particle k
    pin(2)=yk0 ! y position at reaction start for particle k
    pin(3)=zk0 ! z position at reaction start for particle k

    temp1 = tl*(tl+2.*pin(9))
    temp2 = tl+pin(9)
!      call err_chk(1,'LAQ1.f','661',2,temp1)
!      call err_chk(1,'LAQ1.f','661',1,temp2)

    vli=sqrt(temp1)/(temp2)
    pin(8)=tl ! kinetic energy of particle

    temp1 = pin(8)*(pin(8)+2.*pin(9))
!      call err_chk(1,'LAQ1.f','666',2,temp1)
    temp2 = 5.06*sqrt(temp1)
!      call err_chk(1,'LAQ1.f','669',1,temp2)

    dlk=1./(temp2)+delta
    ice=0
    if(ipo == 0)   go  to  12
    k1=int2(k)
    go  to  16
12  if(ice.ne.0)   go  to  13
    call  cenum1(na1,pin,iin(4),dlk,nc,k1,k2,0,1)
    ice=1
13  if(nc.ne.0)   go  to  15
    myp(k)=0
    go  to  26
15  k1=nukc(ice)
    ice=ice+1
    nc=nc-1
16  continue
    dr=((xc(1,k1)-pin(1))*cf+(yc(1,k1)-pin(2))*sf)*st &
         & +(zc(1,k1)-pin(3))*ct
    if(dr <= 0.)    go to  12
    if(ipo == 1.and.nrst == 2.and.k == ist.and.k1 == jst) go to 12
    pin(1)=pin(1)+dr*st*cf
    pin(2)=pin(2)+dr*st*sf
    pin(3)=pin(3)+dr*ct
    dx=pin(1)-xk0
    dy=pin(2)-yk0
    dz=pin(3)-zk0

    temp1 = dx**2+dy**2+dz**2
!      call err_chk(1,'LAQ1.f','697',2,temp1)
    drk=sqrt(temp1)

    temp2 = vli
!      call err_chk(1,'LAQ1.f','701',1,temp2)
    tauk=drk/temp2
    tauk1=tauk*gep*(1.-vli*(vep(1)*st*cf+vep(2)*st*sf+vep(3)+ct))
    taukl=tauk*gpr*(1.+vli*(vpr(1)*st*cf+vpr(2)*st*sf+vpr(3)*ct))
! TAUK  is time interval in projectile rest system
! TAUKL is time interval in observer's rest system
! TAUK1 is time interval in equal velocity  system
!    !!!!!
    if(taukl < pin(7).and.clid < 0.3)   go  to  12
    if(taukl < pin(7).and.ivalon == 0)   go  to  12
!    !!!!!
    int2(k)=k1
    if(tau)  23,23,22
! c22 IF(TAU-TAUK)  26,26,23
22  if(taul-taukl)  26,26,23
23  tau=tauk
    type=tauk1
    taul=taukl
!----> do 24 l=1,9
    do l=1,9
       p1(l)=pin(l)
24     continue
    end do
!----> do 25 l=1,5
    do l=1,5
       ip1(l)=iin(l)
25     continue
    end do
    dl=dlk
    n1=k
    n2=k1
26  if(k-mv) 27,28,28
27  k=k+1
    go to 9
28  continue
    tau=taul
    return
  end subroutine rapid2
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine rapid3(mv,na2,delta,p1,ip1,n1,n2,tau,dl)

! =====================================================================
!
!  Determination of interaction  pair (N1,P1,IP1),N2,
!  where nucleon N1 is number of cascade particle in memory array and
!  N2 is number of spectator nucleon from target nucleus;
!  the time interval TAU (fm/c) is calculated in observer's system
!  (lab. or equal velocity system in dependence of
!  projectile VPR and target VTA nuclei velocities)
!  Calls: CENUM1, CINEMA
!  Input: MV - total number of cascade particles,
!         NA2 -mass number of target nucleus,
!         DL - interaction parameter
!  Output: P1,IP1,N1,N2,TAU
!
! =====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)(a-h,o-z), integer(int32) (i-n)
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(in   ) :: na2
    real(real64),   intent(in   ) :: delta
    real(real64),   intent(inout) :: p1(9)
    integer(int32), intent(inout) :: ip1(5)
    integer(int32), intent(  out) :: n1
    integer(int32), intent(  out) :: n2
    real(real64),   intent(  out) :: tau
    real(real64),   intent(  out) :: dl

    integer(int32), dimension(5) :: iin
    real(real64),   dimension(3) :: v0, cs, c
    real(real64),   dimension(9) :: pin

! =====================================================================

    common/rcor/ rcor
    common/ncasca/ncas,ncpri
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/cenpar/nukc(100)
    common/taue/tpte,type,tyte
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/cslid/clider(5999)
    common/iact/ iact/cvalon/ ivalon

! =====================================================================

    v0(1)=-vta(1)
    v0(2)=-vta(2)
    v0(3)=-vta(3)

    temp1 = 1.-vta(1)**2-vta(2)**2-vta(3)**2
!      call err_chk(1,"LAQ1.f","773",2,temp1)
    temp2 = sqrt(temp1)
!      call err_chk(1,"LAQ1.f","775",1,temp2)
    gta=1./temp2

    temp1 = gta+1.
!      call err_chk(1,"LAQ1.f","779",1,temp1)
    gg=gta*gta/temp1
    k = 1
    tau = -.1
    taul=-.1
9   if(myt(k)) 26,26,10
10  vrk=(pmemo(1,k)-radt(1))*vta(1)+(pmemo(2,k)-radt(2))*vta(2)+ &
         & (pmemo(3,k)-radt(3))*vta(3)
    xk0=pmemo(1,k)-radt(1)+vta(1)*vrk*gg
    yk0=pmemo(2,k)-radt(2)+vta(2)*vrk*gg
    zk0=pmemo(3,k)-radt(3)+vta(3)*vrk*gg
    cs(1)=pmemo(4,k)
    cs(2)=pmemo(5,k)
    cs(3)=pmemo(6,k)
    call  cinema(cs,v0,c,ct,st,cf,sf,tl,pmemo(9,k))
    pin(4)=c(1)
    pin(5)=c(2)
    pin(6)=c(3)
    pin(7)=pmemo(7,k)
    pin(9)=pmemo(9,k)
    pin(8)=tl
    iin(1)=imemo(1,k)
    iin(2)=imemo(2,k)
    iin(3)=imemo(3,k)
    iin(4)=imemo(4,k)
    iin(5)=imemo(5,k)
    clid=clider(k)
    pin(1)=xk0
    pin(2)=yk0
    pin(3)=zk0

    temp1 = tl+pin(9)
!      call err_chk(1,"LAQ1.f","811",1,temp1)
    temp2 = tl*(tl+2.*pin(9))
!      call err_chk(1,"LAQ1.f","814",2,temp2)

    vli=sqrt(temp2)/temp1
    pin(8)=tl

    temp1 = pin(8)*(pin(8)+2.*pin(9))
!      call err_chk(1,"LAQ1.f","819",2,temp1)
    temp2 = 5.06*sqrt(temp1)
!      call err_chk(1,"LAQ1.f","821",1,temp2)
    dlk=1./(temp2)+delta
    ice=0
    if(ipo == 0)   go  to  12
    k1=int3(k)
    go  to  16
12  if(ice.ne.0)   go  to  13
    call  cenum1(na2,pin,iin(4),dlk,nc,k1,k2,0,2)
    ice=1
13  if(nc.ne.0)   go  to  15
    myt(k)=0
    go  to  26
15  k1=nukc(ice)
    ice=ice+1
    nc=nc-1
16  continue
    dr=((xc(2,k1)-pin(1))*cf+(yc(2,k1)-pin(2))*sf)*st+ &
         & (zc(2,k1)-pin(3))*ct
    if(dr <= 0.)    go to  12
    if(ipo == 1.and.nrst == 3.and.k == ist.and.k1 == jst) go to 12
    pin(1)=pin(1)+dr*st*cf
    pin(2)=pin(2)+dr*st*sf
    pin(3)=pin(3)+dr*ct
    dx=pin(1)-xk0
    dy=pin(2)-yk0
    dz=pin(3)-zk0

    temp1 = dx**2+dy**2+dz**2
!      call err_chk(1,"LAQ1.f","849",2,temp1)
    drk=sqrt(temp1)
    temp2 = vli
!      call err_chk(1,"LAQ1.f","852",1,temp2)
    tauk=drk/temp2
    tauk1=tauk*get*(1.-vli*(vet(1)*st*cf+vet(2)*st*sf+vet(3)*ct))
    taukl=tauk*gta*(1.+vli*(vta(1)*st*cf+vta(2)*st*sf+vta(3)*ct))
! TAUK  is time interval in target     rest system
! TAUKL is time interval in observer's rest system
! TAUK1 is time interval in equal velocity  system
!    !!!!!
    if(taukl < pin(7).and.clid < 0.3)   go  to  12
    if(taukl < pin(7).and.ivalon == 0)   go  to  12
!    !!!!!
    int3(k)=k1
    if(tau)  23,23,22
! c22 IF(TAU-TAUK) 26,26,23
22  if(taul-taukl) 26,26,23
23  tau=tauk
    tyte=tauk1
    taul=taukl
!----> do 24 l=1,9
    do l=1,9
       p1(l)=pin(l)
24     continue
    end do
!----> do 25 l=1,5
    do l=1,5
       ip1(l)=iin(l)
25     continue
    end do
    dl=dlk
    n1=k
    n2=k1
26  if(k-mv) 27,28,28
27  k=k+1
    go to 9
28  continue
    tau=taul
    return
  end subroutine rapid3
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  rapidd(mv,tau,md)

! =====================================================================
!
!     Find number MD of nearest in decay time of resonance particle
!     from MV cascade particles
!
! =====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in   ) :: mv
    real(real64),   intent(  out) :: tau
    integer(int32)                :: md

    integer(int32) :: m
    real(real64)   :: taum

! =====================================================================

    real(real64) :: pme
    integer(int32) :: ime
    common /memorylaq/ pme(9,5999),ime(5,5999)
    real(real64) :: tlimit, tint
    common/tlimit/tlimit
    common/actim/tint
    integer(int32) :: idpme
    common/idpme/ idpme(5999)

! =====================================================================

    tau=-.1
    if(mv <= 0)   return
    do m=1,mv
       if(ime(5,m) == 0) cycle
       taum=dble(ime(5,m))/1000.
       if(tau < 0.)  go  to  10
       if(taum > tau) cycle
10     tau=taum
       md=m
    end do
    return
! =====================================================================
  end subroutine rapidd



  subroutine  rapid4(mv,n1,n2,tau)

! =====================================================================
!
!  Determination of  interaction pair (N1,N2) from MV cascade particles
!  the time interval TAU (fm/c) is calculated
!  in observer's system (lab. or equal velocity system depending of
!  projectile VPR and target VTA nuclei velocities)
!  Calls: B2IJ, SGIJ
!  Input: MV- total number of cascade particles
!  Output: N1,N2,TAU
!  For first time all possible numbers of interacting pairs are
!  memorized in INT4 with respective increasing ordered
!  time intervals in TAUIJ4
!
! =====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(  out) :: n1
    integer(int32), intent(  out) :: n2
    real(real64),   intent(  out) :: tau

    real(real64),   dimension(   3) :: rij, vij
    integer(int32), dimension(5999) :: nyy

! =====================================================================

    common/ncasca/ncas,ncpri
    common /memorylaq/ pme(9,5999),ime(5,5999)
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/tauij/ tauk,tpts,typs,tyts,tyys,tij4(100000)
    common/nucoll/ nucoll,mvcoll(5999)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/tlimit/tlimit
    common/actim/tint
    common/ipaul/ip
    common/intcc/intcc
    common/iact/ iact/cvalon/ ivalon
    common/cslid/clider(5999)

! =====================================================================

    tau=-.1
    if(mv <= 1)                                           return
    if(ip == 0.or.ipo == 1)                               go  to  15
!----> do  8  i=1,mv
    do   i=1,mv
8      nyy(i)=myy(i)
    end do
!----> do  13 i=2,mv
    do  i=2,mv
       if(ime(2,i).ne.0)                                     go  to  13
       clidi=clider(i)
       ei=pme(8,i)+pme(9,i)
       inu=0
       if( ime(2,i) == 0.and.ime(3,i) == 0.and.ime(4,i) == 1 &
            & .and.ime(5,i) == 0)  inu=1
       jm=i-1
!----> do  12  j=1,jm
       do   j=1,jm
          if(ime(2,j).ne.0)                                     go  to  12
          if(mvcoll(i) == mvcoll(j))                            go  to  12
          if(myy(i) == 1)                                       go  to  9
          if(myy(j) == 1)                                       go  to  9
          go  to  12
9         jnu=0
          if( ime(2,j) == 0.and.ime(3,j) == 0.and.ime(4,j) == 1 &
               & .and.ime(5,j) == 0)  jnu=1
          if((inu+jnu) < 1.and.intcc == 1)                     go  to  12
          clidj=clider(j)
          ej=pme(8,j)+pme(9,j)
          temp1 = ei
!      call err_chk(1,"LAQ1.f","963",1,temp1)
          temp2 = ej
!      call err_chk(1,"LAQ1.f","968",1,temp2)
!----> do  10  k=1,3
          do   k=1,3
             rij(k)=pme(k,i)-pme(k,j)
             vij(k)=pme(k+3,i)/temp1-pme(k+3,j)/temp2
10           continue
          end do
          rvij=rij(1)*vij(1)+rij(2)*vij(2)+rij(3)*vij(3)
          if(rvij >= 0.)                                        go  to  12
          vij2=vij(1)**2+vij(2)**2+vij(3)**2
          temp1 = vij2
!      call err_chk(1,"LAQ1.f","975",1,temp1)
          tij=-rvij/temp1
!                ! 14.02.05
          if(tij < 0.001.or.(tint+tij) > tlimit)              go  to  12
          if(nrst == 5.and.((ist == i.and.jst == j).or. &
               & (jst == i.and.ist == j)).and. &
               & abs(tij-tauk) <= 0.001)           go  to  12       !!14.02.05
!   !!!!!
          if(tij < pme(7,i).and.clidi < 0.3)                  go  to  12
          if(tij < pme(7,i).and.ivalon == 0)                   go  to  12
          if(tij < pme(7,j).and.clidj < 0.3)                  go  to  12
          if(tij < pme(7,j).and.ivalon == 0)                   go  to  12
!   !!!!!
          call  b2ij(i,j,b2)
          call  sgij(i,j,sig)
          sigv=sig
          if(ivalon.ne.0)    sigv=clidi*clidj*sig
          if(b2 > (sigv/31.41592))                             go  to  12
          if(ijpa < 100000)                                     go  to  11
          write(16,100)
          write( *,100)
100       format(1x,'ijpa>100000')
          go  to  12
11        ijpa=ijpa+1
          tij4(ijpa)=tij
          int4(ijpa)=10000*i+j
!
!     IF(IME(4,I).EQ.0.AND.IME(4,J).EQ.0)
!    *CALL NAMIJ(IP,IPO,INT4(IJPA),I,J,IJPA,TIJ)
!
          nyy(j)=nyy(j)+1
          nyy(i)=nyy(i)+1
12        continue
       end do
13     continue
    end do
!----> do  14  i=1,mv
    do   i=1,mv
       if(nyy(i) == myy(i).and.myy(i) <= 1)    nyy(i)=0
       if(nyy(i) == 1.and.myy(i) == 0)         nyy(i)=2
       myy(i)=nyy(i)
14     continue
    end do
15  continue
    if(ijpa <= 0)                                          return
    ijm=1
    ijpa0=ijpa
!----> do  19  ijp=1,ijpa0
    do   ijp=1,ijpa0
16     if(ijp > ijpa)                                      go  to  20
       tij=tij4(ijp)
       if(tij < 0.000001)                                  go  to  17
       ij=int4(ijp)
       i=(ij+1)/10000
       j=ij-10000*i
       if(i == j)                                           go  to  17
       if(i < 1.or.i > 5999.or.j < 1.or.j > 5999) &
            & write(16,101) ij,i,j
101    format(' ij,i,j=',3i10)
       if(ipo == 1.and.nrst == 5.and.(i+j) == (ist+jst) &
            & .and.(i*j) == (ist*jst))    go  to  17
       if((tint+tij) > tlimit)                             go  to  17
!   !!!!!
       if(tij < pme(7,i).and.clider(i) < 0.3)             go  to  17
       if(tij < pme(7,i).and.ivalon == 0)                  go  to  17
       if(tij < pme(7,j).and.clider(j) < 0.3)             go  to  17
       if(tij < pme(7,j).and.ivalon == 0)                  go  to  17
!   !!!!!
       if(tau < 0.)                                        go  to  18
       if(tij > tau)                                       go  to  19
       go  to  18
17     tij4(ijp)=tij4(ijpa)
       int4(ijp)=int4(ijpa)
       ijpa=ijpa-1
       if(ijpa <= 0)                                        go  to  21
       go  to  16
18     tau=tij
       n1=i
       n2=j
       ijm=ijp
19     continue
    end do
20  continue
!      IJ=INT4(IJM)
!      INT4(IJM)=INT4(IJPA)
!      TIJ4(IJM)=TIJ4(IJPA)
!      IJPA=IJPA-1
21  if(tau < 0.)                                          return
!
!     IF(IME(4,N1).EQ.0.AND.IME(4,N2).EQ.0)
!    *CALL NAMIJ(IP,IPO,IJ,N1,N2,IJPA,TAU)
!
    if(intcc >= 2)                                         return
    if((ime(4,n1)+ime(4,n2)) == 2.or.ime(4,n2) == 1)       return
    m=n2
    n2=n1
    n1=m
    return
! =====================================================================
  end subroutine rapid4
