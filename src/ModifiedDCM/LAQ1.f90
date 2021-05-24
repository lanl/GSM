
!          For  LAQGSM-MARS interface      KKG  04/09/07
! aa2k7g.f,  gengamn7.f, gadd7.f,    ->     laqgsm2007_1.f
! qgsm7.f                            ->     laqgsm2007_2.f
! coales07.f, gemaa07.f, preco07g.f, spegem7g.f, fermi07.f,
! gambetm7.f                         ->     laqgsm2007_3.f
!
!
!      =================================================================
!          laqgsm2007_1.f   subroutines :
!      =================================================================
!
!  Dimensions of arrays UP1 and UP2 are 66000 and 70000
!  Dimensions of arrays XC,YC,ZC,IZ,MPA,INT1 are 300,kkg 23.11.04
! Fermi energies are changed via TF=TF*(AN1/ANUCL1)**(2./3.)
! see PANTNQ, POTENQ, TFERMIQ;                              20.06.95
! anti-baryon potential is added                         05.11.95
! introduced /hadr1/                                     16.03.96
! interactions in RIJ<DELTA are included in RAPID4       01.10.97
! pi+Y=>AK+N channel is added                            08.02.99
! scattering width in ELASLE is added for Deltas         27.04.99
! corr. B+B==>B+Y+K (DDNYK,DDDYK,DNNYK,DNDYK,NNNYK,NNDYK 14.05.01
! on unit 10 3 numbers are writen instead of IS9         21.05.01
! last corrections in PAULI2,PAULI3             26.06.02,13.12.04
! new angular distributions for n+p,p+p at T<2.0GeV(ELEX)22.10.03
! gamma as projectile is included                        28.10.03
! extended to Egamma up to 10 GeV, Dec. 2004
! corrections in RAPID4                                  14.02.05
! last changes by KKG for interface with MARS code       23.03.07
!  /isecgpn/ ==> /ixsgpn/, KINEMA ==> KINEMQ             28.03.07
!
  subroutine pointn(m,p1,p2,ip1,ip2,v1,u1,tr2,sig,sab,tint,nr,n1, &
       & n2,na1,na2,mv,delta,r0x,r0y,r0z)
! ======================================================================
!
!  Determination of interaction hadron pair (N1,P1,IP1),(N2,P2,IP2)
!  and type NR of collision:
!  NR =1,      collision of nucleon N1  from  projetile and
!                           nucleon N2  from  target;
!  NR =2,      collision of cascade particle N1 and
!                           nucleon N2  from  projetile
!  NR =3,      collision of cascade particle N1 and
!                           nucleon N2  from  target;
!  NR =5,      collision of cascade particle N1 and
!                           cascade particle N2
!  NR =4,      decay of resonance N1.
!  Interaction time TINT is replaced by TINT + TAU, where
!  TAU is min(TAU_NR), NR=1,2,3,4,5.
!  Calls:  RAPID1,RAPID2,RAPID3,RAPIDD and RAPID4 - to calculate
!          time of next interaction TAU_NR
!

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
    common/indint/indint
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/taue/tpte,type,tyte
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/dtint/dtint
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t00,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/ncasca/ncas,ncpri
    common/pidabs/pidabs
    common/intcc/ intcc
    common/tauij/ tauk,tpts,typs,tyts,tyys,tij4(100000)
    common/tlimit/tlimit
    common/cslid/clider(5999)
    common/iact/ iact/cvalon/ ivalon
    common /holpt/ ph1(3),ph2(3),rh1(3),rh2(3),ehol1,ehol2,tfp,tft
    common /idpme/ idpme(5999)
    common /idn12/ id1,id2
    common /barpot/ pot
    common /parinc/ pinc(3),einc
    dimension p1(9),p2(9),ip1(5),ip2(5),v1(3) &
         & ,pt1(9),pt2(9),ipt1(5),ipt2(5),yp1(9),iyp1(5),yt1(9),iyt1(5)
    tfro1=tint

8   continue
    r02 = r0x**2+r0y**2+r0z**2
    ! call err_chk(1,'laq1.f','84',2,r02)
    r0=sqrt(r02)*1.01
    ipo=0
9   tpt=-0.1
    typ=-0.1
    tyt=-0.1
    taud=-0.1
    tyy=-0.1
    tau=-0.1
    ty=-.1
    tfro2=tint
!      write( *,1991) NCAS,TINT
1991 format('+',50x,i6,1x,e10.3)
    temp1 = 0
!----> do 9999 i=1,3
    do i=1,3
       temp1 = temp1 + (radp(i)-radt(i))**2
9999   continue
    end do
      ! temp1 = (RADP(1)-RADT(1))**2+(RADP(2)-RADT(2))**2+(RADP(3)-RADT(3))**2
    ! call err_chk(1,'laq1.f','98',2,temp1)
    rt=sqrt(temp1)
    if(rt-r0)10,10,12
10  if(na1 < 1.or.na2 < 1)   go  to  12
    if(ipo.ne.1)               go  to  11
    if(tpts <= 0.)             go  to  12
    if((tpts-tauk) <= 0.0001)  go  to  11
    tpt=tpts-tauk
    go  to  12
11  continue
    call  rapid1(na1,na2,delta,pt1,ipt1,pt2,ipt2,npt1,npt2,tpt,dlpt)
12  if(mv) 124,124,100
100 if(na1 <= 0)                             go  to  101
    if(ipo.ne.1)                             go  to  13
    if(typs <= 0.)                           go  to  101
    if((typs-tauk) <= 0.0001)                go  to  13
    typ=typs-tauk
    if(typ < pmemo(7,nyp1).and.clider(nyp1) < 0.3) &
         & go  to  13
    if(typ < pmemo(7,nyp1).and.ivalon == 0) &
         & go  to  13
    go  to  101
13  continue
    call rapid2(mv,na1,delta,yp1,iyp1,nyp1,nyp2,typ,dlyp)
101 if(na2 < 1)   go  to  14
102 continue
    if(ipo.ne.1)                             go  to  103
    if(tyts <= 0.)                           go  to  14
    if((tyts-tauk) <= 0.0001)                go  to  103
    tyt=tyts-tauk
    if(tyt < pmemo(7,nyt1).and.clider(nyt1) < 0.3) &
         & go  to  103
    if(tyt < pmemo(7,nyt1).and.ivalon == 0) &
         & go  to  103
    go  to  14
103 continue
    call rapid3(mv,na2,delta,yt1,iyt1,nyt1,nyt2,tyt,dlyt)
14  continue
    call  rapidd(mv,taud,nd)
1992 continue
    if(mv < 2.or.intcc == 0.or.tint > tlimit) go  to 124
    if(ipo.ne.1)                               go  to  114
    if(tyys <= 0.)                             go  to  124
    if((tyy-tauk) <= 0.0001)                   go  to  114
    tyy=tyys-tauk
    if(tyy < pmemo(7,nyy1).and.clider(nyy1) < 0.3) &
         & go  to 114
    if(tyy < pmemo(7,nyy2).and.clider(nyy2) < 0.3) &
         & go  to 114
    if(tyy < pmemo(7,nyy1).and.ivalon == 0) &
         & go  to 114
    if(tyy < pmemo(7,nyy2).and.ivalon == 0) &
         & go  to 114
    if(nyy1 == nyy2)                            go  to 114
    go  to  124
114 continue
    call  rapid4(mv,nyy1,nyy2,tyy)
124 continue
    if(typ)19,19,15
15  if(tyt)16,16,17
16  ty=typ
    nr=2
    clid1=clider(nyp1)
    clid2=1.
!
    go  to  20
17  if(typ-tyt) 16,16,18
18  ty=tyt
    nr=3
    clid1=clider(nyt1)
    clid2=1.
    go to 20
19  if(tyt)20,20,18
20  if(tpt)24,24,21
21  if(ty)23,23,22
22  if(tpt-ty)  23,23,25
23  tau=tpt
    nr=1
    clid1=1.
    clid2=1.
!
    go  to   26
24  if(ty)26,26,25
25  tau=ty
!
    go  to   26
26  if(taud < 0.)                go  to  126
    if(taud > tau.and.tau >= 0.) go  to  126
    tau=taud
    n1=nd
    nr=4
126 if(tyy < 0.)                 go  to  226
    if(tyy > tau.and.tau >= 0.)  go  to  226
    tau=tyy
    n1=nyy1
    n2=nyy2
    nr=5
    clid1=clider(nyy1)
    clid2=clider(nyy2)
226 if(tau < 0.)                 go  to  41
!
    if(ncas >= ncpri)  then
!      write(16,6000) MV,IPO,NR,TPT,TYP,TYT,TYY,TAUD,TAU
6000   format(1x,'mv,ipo,nr,tpt,typ,tyt,tyy,taud,tau=',i4,2i2,6(1x,f6.3))
    endif
!
    tauk=tau
    if(mv)29,29,27
!----> do 28 k=1,mv
27  do k=1,mv
       ek=pmemo(8,k)+pmemo(9,k)

       temp1 = ek
      ! TAUK value is SMALL, skipping error protection.
       ! call err_chk(1,'laq1.f','212, 213, 214',1,temp1)

       pmemo(1,k)=pmemo(1,k)+pmemo(4,k)/temp1*tauk
       pmemo(2,k)=pmemo(2,k)+pmemo(5,k)/temp1*tauk
       pmemo(3,k)=pmemo(3,k)+pmemo(6,k)/temp1*tauk
       pmemo(7,k)=pmemo(7,k)-tauk
       if(pmemo(7,k) > 0.)          go  to  128
       pmemo(7,k)=0.
       clider(k)=1.
128    continue
       if(imemo(5,k) == 0)           go  to  28
       imemo(5,k)=imemo(5,k)-intg(tauk*1000.)
       if(imemo(5,k) <= 0)  imemo(5,k)=1
28     continue
    end do
29  tint=tint+tau
    if(ijpa <= 0)                 go  to  291
!----> do  290  ij=1,ijpa
    do   ij=1,ijpa
290    tij4(ij)=tij4(ij)-tauk
    end do
291 continue
!----> do  30  k=1,3
    do   k=1,3
       radt(k)=radt(k)+tau*vta(k)
       radp(k)=radp(k)+tau*vpr(k)
30     continue
    end do
    go  to (33,35,37,130,230),nr
!  ====> RAPIDD
130 continue
    id1=idpme(n1)
    m=1
    go  to  135
!  ====> RAPID4
!----> do 231  k=1,9
230 do  k=1,9
       p1(k)=pmemo(k,n1)
       p2(k)=pmemo(k,n2)
       if(k <= 5)   ip1(k)=imemo(k,n1)
       if(k <= 5)   ip2(k)=imemo(k,n2)
231    continue
    end do
    id1=idpme(n1)
    id2=idpme(n2)
    go  to  32
!  ====> RAPID1
!----> do  34  k=1,9
33  do   k=1,9
       p1(k)=pt1(k)
       p2(k)=pt2(k)
       if(k <= 5)   ip1(k)=ipt1(k)
       if(k <= 5)   ip2(k)=ipt2(k)
34     continue
    end do
    n1=npt1
    n2=npt2
    dlk=dlpt
    id1=1120
    if(ip1(1) == 0) id1=1220
    id2=1120
    if(ip2(1) == 0) id2=1220
    go  to   32
!  ====> RAPID2
!----> do  36  k=1,9
35  do   k=1,9
       p1(k)=yp1(k)
       if(k <= 5)   ip1(k)=iyp1(k)
36     continue
    end do
    pinc(1)=yp1(4)
    pinc(2)=yp1(5)
    pinc(3)=yp1(6)
    temp1 = pinc(1)**2+pinc(2)**2+pinc(3)**2+yp1(9)**2
    ! call err_chk(1,'laq1.f','274',2,temp1)
    einc=sqrt(temp1)
    n1=nyp1
    n2=nyp2
    dlk=dlyp
    go  to  39
!  ====> RAPID3
!----> do  38  k=1,9
37  do   k=1,9
       p1(k)=yt1(k)
       if(k <= 5)   ip1(k)=iyt1(k)
38     continue
    end do
    pinc(1)=yt1(4)
    pinc(2)=yt1(5)
    pinc(3)=yt1(6)
    temp = pinc(1)**2+pinc(2)**2+pinc(3)**2+yt1(9)**2
    ! call err_chk(1,'laq1.f','288',2,temp)
    einc=sqrt(temp)
    n1=nyt1
    n2=nyt2
    dlk=dlyt
    call  partnq(2,n2,p2,ip2)
    tft=tfermiq(p2(1),p2(2),p2(3),2)
    ehol2=p2(8)
    p2(9)=0.94
    pot=0.
    if(ip1(3) == 0.and.ip1(5) == 0.and.ip1(4) == 1)  pot=tft+eps2
!  kkg 28.10.03
!      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0.AND.IP1(4).EQ.0)  POT=VPI
    if(ip1(3) == 0.and.ip1(5) == 0.and.ip1(4) == 0 &
         & .and.abs(p1(9)-0.140) < 0.1)                  pot=vpi
! ANTI-NUCLEON
    if(ip1(3) == 0.and.ip1(5) == 0.and.ip1(4) == -1) &
         & pot=potbar(2,p2(1),p2(2),p2(3))
    p1(8)=p1(8)+pot

    temp1 = p1(8)-pot
    temp2 = p1(8)-pot+2.*p1(9)
    ! call err_chk(1,'laq1.f','312',1,temp1)
    ! call err_chk(1,'laq1.f','312',1,temp2)
    temp3 = p1(8)/(temp1) * (p1(8)+2.*p1(9))/(temp2)
    ! call err_chk(1,'laq1.f','315',2,temp3)

    c=sqrt( temp3 )
!----> do  339 k=1,3
    do  k=1,3
       p1(k+3)=p1(k+3)*c
       ph2(k)=-p2(3+k)
339    rh2(k)=p2(k)
    end do
    id1=idpme(n1)
    id2=1120
    if(ip2(1) == 0) id2=1220
    go  to  32
39  continue
    call  partnq(1,n2,p2,ip2)
    tfp=tfermiq(p2(1),p2(2),p2(3),1)
    ehol1=p2(8)
    p2(9)=0.94
    pot=0.
    if(ip1(3) == 0.and.ip1(5) == 0.and.ip1(4) == 1)  pot=tfp+eps1
    if(ip1(3) == 0.and.ip1(5) == 0.and.ip1(4) == 0)  pot=vpi
! ANTI-NUCLEON
    if(ip1(3) == 0.and.ip1(5) == 0.and.ip1(4) == -1) &
         & pot=potbar(1,p2(1),p2(2),p2(3))
    p1(8)=p1(8)+pot

    temp1 = p1(8)-pot
    temp2 = p1(8)-pot+2.*p1(9)
    ! call err_chk(1,'laq1.f','341',1,temp1)
    ! call err_chk(1,'laq1.f','341',1,temp2)
    temp3 = p1(8)/(temp1) * (p1(8)+2.*p1(9))/(temp2)
    ! call err_chk(1,'laq1.f','344',2,temp3)

    c=sqrt( temp3 )
!----> do  239 k=1,3
    do  k=1,3
       p1(k+3)=p1(k+3)*c
       ph1(k)=-p2(3+k)
239    rh1(k)=p2(k)
    end do
    id1=idpme(n1)
    id2=1120
    if(ip2(1) == 0) id2=1220
    go  to  32
32  continue
!
    if((id1 == 120.or.id1 == -120.or.id1 == 110.or.id1 == 1120.or. &
!  kkg  28.10.03     ! including gamma(ID1=10)
!     &    ID1.EQ.1220).AND.(ID2.EQ.1120.OR.ID2.EQ.1220)) THEN
         & id1 == 1220.or.id1 == 10) &
         & .and.(id2 == 1120.or.id2 == 1220)) then
!  kkg  28.10.03
       call  slqekq(l,ms,mq,ksin,me,ip1,ip2)
       call  tinvuq(p1,p2,u1,v1,tr2)
       if(id1 == 10)  then
          sito= csgntot(ip2(1),tr2,p2(9))/1000.   ! tot. g+n cr.sec.
       else
          sito=croseg(l,ms,mq,ksin,0,tr2,p1(9),ip1(5))
       endif
    else
       px1=p1(4)
       py1=p1(5)
       pz1=p1(6)
       am1=p1(9)
       px2=p2(4)
       py2=p2(5)
       pz2=p2(6)
       am2=p2(9)
       call crosec(1,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
       call  slqekq(l,ms,mq,ksin,me,ip1,ip2)
       call  tinvuq(p1,p2,u1,v1,tr2)
    endif
    sig=sito
    sab=0.
    go  to  (135,131,132,135,135),nr
131 if(an1 < 2.1.or.zn1 < 1.1)          go  to  135
    go  to  133
132 if(an2 < 2.1.or.zn2 < 1.1)          go  to  135
133 continue
!  kkg  28.10.03     ! including gamma(ID1=10)
    if(id1 == 120.or.id1 == -120.or.id1 == 110.or.id1 == 10) &
         & go  to  134
    go  to  135
134 continue
!  kkg  28.10.03     ! including gamma(ID1=10)
!   !!!  ONLY FOR PION(+-0) and gamma  !!!
    sab=croseg(l,ms,mq,ksin,3,p1(8),p1(9),ip1(5))
    if(id1.ne.10)  sab=sab*pidabs
!  kkg  28.10.03
135 continue
    ipo=1
    ist=n1
    jst=n2
    nrst=nr
    tpts=tpt
    typs=typ
    tyts=tyt
    tyys=tyy
    if(nr == 4)        go  to  42
    sigv=sig
!  kkg  28.10.03     ! not include  gamma(ID1=10)
    if(ivalon.ne.0.and.id1.ne.10)    sigv=clid1*clid2*sig
!  kkg  04/02/04
    if(id1 == 10)      sigv=sig+sab
    if(indint == 2)    go  to  139
    if(nr == 5)        go  to  42

    temp1 = 31.41592*(dlk**2)
    ! call err_chk(1,'laq1.f','419',1,temp1)

    temp1=1.-sigv/temp1
    drnd=rndm(-1.0_real64)
    if(drnd < temp1)  go  to  9
    go  to  42
139 temp1=sigv/31.41592
    bi2=bim2(p1,p2)
    if(bi2 > temp1)   go  to  9
    go  to  42
41  m=0
    tsh2=0.
!      IF(NCAS.GE.NCPRI)  write( *,*) ' M=0'
    return
42  continue
    m=1
    dtint=tint-tfro1
    if(ncas < ncpri)  return
    write(16,665) ijpa,nr,tauk,tint
665 format(1x,'ijpa,nr,tauk,tint=',i5,i2,2(1x,f7.3))
    if(nr >= 2.and.nr <= 4) &
         & write(16,666) nr,n1,(pmemo(k,n1),k=1,9),(imemo(k,n1),k=1,5), &
         & n2,(p2(k),k=1,9),(ip2(k),k=1,5)
666 format(' nr,n1=',i1,i5,9(1x,f10.3),4i2,i8/ &
         & '   +n2=',1x,i5,9(1x,f10.3),4i2,i8)
    if(nr == 1.or.nr == 5) &
         & write(16,667) nr,n1,(p1(k),k=1,9),(ip1(k),k=1,5), &
         & n2,(p2(k),k=1,9),(ip2(k),k=1,5)
667 format(' nr,n1=',i1,i5,9(1x,f10.3),4i2,i8/ &
         & '   +n2=',1x,i5,9(1x,f10.3),4i2,i8)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! ==============================================================================
! Moved RAPID routines to new file, rapid.f
! ==============================================================================

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine eraij(n)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Erase the particle N from INT4, TIJ4 arrays
!
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/tauij/ tauk,tpts,typs,tyts,tyys,tij4(100000)
    if(ijpa <= 0)          return
    ijpa0=ijpa
!----> do  12  ijp=1,ijpa0
    do   ijp=1,ijpa0
10     if(ijp > ijpa)        return
       ij=int4(ijp)
       i=(ij+1)/10000
       j=ij-10000*i
       if(n == i.or.n == j)  go  to  11
       go  to  12
11     tij4(ijp)=tij4(ijpa)
       int4(ijp)=int4(ijpa)
       ijpa=ijpa-1
       go  to  10
12     continue
    end do
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine repij(n,m)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Replace the particle M to particle in INT4 array
!
    common/intcen/ipo,int1(300),int2(5999),int3(5999), &
         & int4(100000),ijpa,ist,jst,nrst
    common/tauij/ tauk,tpts,typs,tyts,tyys,tij4(100000)
    if(ijpa <= 0)          return
!----> do  11  ijp=1,ijpa
    do   ijp=1,ijpa
       ij=int4(ijp)
       i=(ij+1)/10000
       j=ij-10000*i
       if(m == i.or.m == j)  go  to  10
       go  to  11
10     if(m == i)  i=n
       if(m == j)  j=n
       ij=10000*i+j
       int4(ijp)=ij
11     continue
    end do
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine namij(ip,ipo,ij,i,j,ijpa,tau)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Determines the name of interacting particles I,J
!     from memorylaq array
!
    character(len=8) :: pni,pnj
    common /memorylaq/ pme(9,5999),ime(5,5999)
    common /idpme/ idpme(5999)
    dimension pi(9),ipi(5),pj(9),ipj(5)
!----> do  10  k=1,9
    do   k=1,9
       pi(k)=pme(k,i)
       pj(k)=pme(k,j)
10     continue
    end do
!----> do  11  k=1,5
    do   k=1,5
       ipi(k)=ime(k,i)
       ipj(k)=ime(k,j)
11     continue
    end do
    idi=idpme(i)
    call panuid(idi,iki,pni)
    idj=idpme(j)
    call panuid(idj,ikj,pnj)
!     IF(IPI(4).EQ.0.AND.IPJ(4).EQ.0)
!    *write(16,100) IP,IPO,PNI,PNJ,IJ,I,J,IJPA,TAU
    if(ipi(4) == 0.and.ipj(4) == 0) &
         & write( *,100) ip,ipo,pni,pnj,ij,i,j,ijpa,tau
100 format(11x,2i3,1x,2a8,i7,2i5,i6,f8.3)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  b2ij(i,j,b2)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculates the impact parameter B2=B**2 of cascade particles I,J
!
    common /memorylaq/ pme(9,5999),ime(5,5999)
    dimension r(3),p(3),rs(3),ps(3)
    ei=pme(8,i)+pme(9,i)
    ej=pme(8,j)+pme(9,j)
    e=ei+ej
    pr=0.
    ppi=0.
    ppj=0.
!----> do  10  k=1,3
    do   k=1,3
       r(k)=pme(k,i)-pme(k,j)
       p(k)=pme(k+3,i)+pme(k+3,j)
       pr=pr+p(k)*r(k)
       ppi=ppi+p(k)*pme(k+3,i)
       ppj=ppj+p(k)*pme(k+3,j)
10     continue
    end do
    temp1 = e**2-p(1)**2-p(2)**2-p(3)**2
    ! call err_chk(1,"laq1.f","1167",2,temp1)
    u=sqrt(temp1)
    temp2 = 2.*u
    ! call err_chk(1,"laq1.f","1170",1,temp2)
    eis=(u**2+pme(9,i)**2-pme(9,j)**2)/(temp2)
    ejs=u-eis
    prs=0.
    temp1 = u
    temp2 = e + u
    ! call err_chk(1,"laq1.f","1179, 1180, 1181",1,temp1)
    ! call err_chk(1,"laq1.f","1179, 1180, 1181",1,temp2)
!----> do  11  k=1,3
    do   k=1,3
       ps(k)=(pme(k+3,i)+p(k)/temp1*(ppi/temp2-ei)- &
            & pme(k+3,j)-p(k)/temp1*(ppj/temp2-ej))/2.
       rs(k)=r(k)+p(k)*pr/temp2/temp1
       prs=prs+ps(k)*rs(k)
11     continue
    end do
    ps2=ps(1)**2+ps(2)**2+ps(3)**2
    rs2=rs(1)**2+rs(2)**2+rs(3)**2
    temp1 = ps2
    ! call err_chk(1,"laq1.f","1188",1,temp1)
    b2=rs2-(prs**2)/temp1
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  bim2(p1,p2)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculates the impact parameter BIM2=B**2 of particles P1,P2
!     Definition of P1(P2):
!                          P1(1) - x coordinate in observer's system
!                          P1(2) - y coordinate in observer's system
!                          P1(3) - z coordinate in observer's system
!                          P1(4) - x component of momentum (GeV/c)
!                          P1(5) - y component of momentum (GeV/c)
!                          P1(6) - z component of momentum (GeV/c)
!                          P1(7) - the maturity (formation) time (fm/c)
!                          P1(8) - kinetic energy (GeV)
!                          P1(9) - mass (GeV)
!
    dimension r(3),p(3),rs(3),ps(3),p1(9),p2(9)
    ei=p1(8)+p1(9)
    ej=p2(8)+p2(9)
    e=ei+ej
    pr=0.
    ppi=0.
    ppj=0.
!----> do  10  k=1,3
    do   k=1,3
       r(k)=p1(k)-p2(k)
       p(k)=p1(k+3)+p2(k+3)
       pr=pr+p(k)*r(k)
       ppi=ppi+p(k)*p1(k+3)
       ppj=ppj+p(k)*p2(k+3)
10     continue
    end do
    temp1 = e**2-p(1)**2-p(2)**2-p(3)**2
    ! call err_chk(1,"laq1.f","1223",2,temp1)
    u=sqrt(temp1)
    temp2 = 2.*u
    ! call err_chk(1,"laq1.f","1226",1,temp2)
    eis=(u**2+p1(9)**2-p2(9)**2)/temp2
    ejs=u-eis
    prs=0.
    temp1 = u
    temp2 = e + u
    ! call err_chk(1,"laq1.f","1234, 1235, 1236, 1237",1,temp1)
    ! call err_chk(1,"laq1.f","1234, 1235, 1236, 1237",1,temp2)
!----> do  11  k=1,3
    do   k=1,3
       ps(k)=(p1(k+3)+p(k)/temp1*(ppi/temp2-ei)- &
            & p2(k+3)-p(k)/temp1*(ppj/temp2-ej))/2.
       ps(k)=(p1(k+3)*(ej+ejs)-p2(k+3)*(ei+eis))/temp2
       rs(k)=r(k)+p(k)*pr/temp2/temp1
       prs=prs+ps(k)*rs(k)
11     continue
    end do
    ps2=ps(1)**2+ps(2)**2+ps(3)**2
    rs2=rs(1)**2+rs(2)**2+rs(3)**2
    temp1 = ps2
    ! call err_chk(1,"laq1.f","1244",1,temp1)
    bim2=rs2-(prs**2)/temp1
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  sgij(i,j,sg)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   Calculates the total interaction cross section SG of particles I,J
!
    character(len=8) :: pna1,pna2
    common/memorylaq/pme(9,5999),ime(5,5999)
    common /idpme/ idpme(5999)
    dimension p1(9),p2(9),ip1(5),ip2(5)
!----> do  10  k=1,9
    do   k=1,9
       p1(k)=pme(k,i)
       p2(k)=pme(k,j)
       if(k <= 5)  ip1(k)=ime(k,i)
       if(k <= 5)  ip2(k)=ime(k,j)
10     continue
    end do
    id1=idpme(i)
    call panuid(id1,ik1,pna1)
    id2=idpme(j)
    call panuid(id2,ik2,pna2)
    px1=p1(4)
    py1=p1(5)
    pz1=p1(6)
    am1=p1(9)
    px2=p2(4)
    py2=p2(5)
    pz2=p2(6)
    am2=p2(9)
    call crosec(1,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sg=sito
    return
  end

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine rpts(x,y,z,x1,y1,z1,nu)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Transforms the observer system's coordinates X,Y,Z into
!     nucleus NU rest system
!
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    dimension v(3),r(3)
!----> do  11  k=1,3
    do   k=1,3
       if(nu == 2)  go  to  10
       v(k)=vpr(k)
       r(k)=radp(k)
       go  to  11
10     v(k)=vta(k)
       r(k)=radt(k)
11     continue
    end do
    vr=(x-r(1))*v(1)+(y-r(2))*v(2)+(z-r(3))*v(3)
    g=1./sqrt(1.-v(1)**2-v(2)**2-v(3)**2)
    gg=g*g/(g+1.)
    x1=x-r(1)+v(1)*vr*gg
    y1=y-r(2)+v(2)*vr*gg
    z1=z-r(3)+v(3)*vr*gg
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  tfermiq(x,y,z,nu)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculates the Fermi energy of nucleus NU
!
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    tfermiq=0.
    if(nu == 1.and.anucl1 < 2.1)  return
    if(nu == 2.and.anucl2 < 2.1)  return
    if(nu == 1)  then
       an=anucl1
       a=a1
       c=c1
       d=d1
       tf0=tf01
       rm=rm1
       ares=an1
    else
       an=anucl2
       a=a2
       c=c2
       d=d2
       tf0=tf02
       rm=rm2
       ares=an2
    endif
    r=sqrt(x**2+y**2+z**2)
    if((r/rm) > 1.5)  return
    if(an <= 10.)  then
       tfermiq=tf0*exp(-(2./3.)*(r/a)**2) *(ares/an)**0.6666667
    else
       tfermiq=tf0*( (1.+exp(-a/c))/(1.+exp((r-a)/c)) &
            & *(ares/an) )**0.6666667            ! 09.04.95
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  dissip(ip,t0)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Change the projectile and target momenta due of Coulomb interaction
!     Calls: COTRAN, CINEMA
!
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    dimension ppl(3),ptl(3)
    if(ip == -1)  go to  14
    if(t0 < 5.0.and.an1 > 1.1) call  cotran
    if(ip)  14,14,9
9   if(an1 < 0.5)   go  to  11
    amp=0.940*an1
    call  cinema(pnucl1,vpr,ppl,ct,st,cf,sf,tlp,amp)
    ep=tlp+amp
    temp = ep
    ! call err_chk(1,'laq1.f','1463, 1464, 1465',1,temp)
    vpr(1)=ppl(1)/temp
    vpr(2)=ppl(2)/temp
    vpr(3)=ppl(3)/temp
!----> do  10  k=1,3
    do   k=1,3
10     pnucl1(k)=0.0
    end do
11  if(an2 < 0.5)  go  to  14
    amt=0.940*an2
    call  cinema(pnucl2,vta,ptl,ct,st,cf,sf,tlt,amt)
    et=tlt+amt
    temp = et
    ! call err_chk(1,'laq1.f','1474, 1475, 1476',1,temp)
    vta(1)=ptl(1)/temp
    vta(2)=ptl(2)/temp
    vta(3)=ptl(3)/temp
!----> do  12  k=1,3
    do   k=1,3
12     pnucl2(k)=0.0
    end do
14  temp1 = 1.-vpr(1)**2-vpr(2)**2-vpr(3)**2
    ! call err_chk(1,'laq1.f','1481',2,temp1)
    temp2 = sqrt(temp1)
    ! call err_chk(1,'laq1.f','1482',1,temp2)
    gp=1./temp2

    temp1 = 1.-vta(1)**2-vta(2)**2-vta(3)**2
    ! call err_chk(1,'laq1.f','1487',2,temp1)
    temp2 = sqrt(temp1)
    ! call err_chk(1,'laq1.f','1489',1,temp1)
    gt=1./temp2
    vtp=vpr(1)*vta(1)+vpr(2)*vta(2)+vpr(3)*vta(3)

    temp1 = gt+1.
    temp2 = 1.-vtp
    temp3 = gt+gp
    temp4 = gt
    ! call err_chk(1,'laq1.f','1505',1,temp1)
    ! call err_chk(1,'laq1.f','1505',1,temp2)
    ! call err_chk(1,'laq1.f','1506',1,temp3)
    ! call err_chk(1,'laq1.f','1505',1,temp4)
!----> do  15  k=1,3
    do   k=1,3
       vre(k)=(vpr(k)+gt*vta(k)*(vtp*gt/temp1-1.))/(temp2)/temp4
15     vev(k)=(vpr(k)*gp+vta(k)*gt)/(temp3)
    end do

    temp1 = 1.-vev(1)**2-vev(2)**2-vev(3)**2
    temp2 = 1.-vre(1)**2-vre(2)**2-vre(3)**2
    ! call err_chk(1,'laq1.f','1508',2,temp1)
    ! call err_chk(1,'laq1.f','1509',2,temp2)
    temp1 = sqrt(temp1)
    temp2 = sqrt(temp2)
    ! call err_chk(1,'laq1.f','1512',1,temp1)
    ! call err_chk(1,'laq1.f','1513',1,temp2)
    gev=1./temp1
    gre=1./temp2
    set=vev(1)*vta(1)+vev(2)*vta(2)+vev(3)*vta(3)
    sep=vev(1)*vpr(1)+vev(2)*vpr(2)+vev(3)*vpr(3)

    temp3 = gt+1.
    temp4 = 1.-set
    temp7 = gt
    ! call err_chk(1,'laq1.f','1531',1,temp3)
    ! call err_chk(1,'laq1.f','1531',1,temp4)
    ! call err_chk(1,'laq1.f','1531',1,temp7)
    temp5 = gp+1.
    temp6 = 1.-sep
    temp8 = gp
    ! call err_chk(1,'laq1.f','1532',1,temp5)
    ! call err_chk(1,'laq1.f','1532',1,temp6)
    ! call err_chk(1,'laq1.f','1532',1,temp8)

!----> do 115 k=1,3
    do k=1,3
       vet(k)=(vev(k)+gt*vta(k)*(set*gt/(temp3)-1.))/(temp4)/temp7
       vep(k)=(vev(k)+gp*vpr(k)*(sep*gp/(temp5)-1.))/(temp6)/temp8
115    continue
    end do
    get=gev*gt*(1.-set)
    gep=gev*gp*(1.-sep)
16  continue
200 format(5x,'vpr,vta,gp,gt',2(3(f7.4,2x),5x),2(f8.4,2x))
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! ==============================================================================
!
! Moved cenum1 to its own file - this routine is a 'bottleneck' routine
!
! ==============================================================================

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  abelq(pin,v,u,p1,p2,ct,fi,cm1,cm2)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculates the final momenta of elastic interacting particles
!     with masses CM1,CM2 of total cms energy U, velocity V
!     Calls: CINEMA, ROTORQ;
!     scattering angles CT and FI are as input, P1,P2 - output
!
    dimension  pin(9),v(3),b(3),p1(3),p2(3),pinl(3),pins(3)
    b(1)=-v(1)
    b(2)=-v(2)
    b(3)=-v(3)
    pinl(1)=pin(4)
    pinl(2)=pin(5)
    pinl(3)=pin(6)
    call  cinema(pinl,b,pins,cts,sts,cfs,sfs,ts,pin(9))
    st=sqrt(1.-ct*ct)
    cf=cos(fi)
    sf=sin(fi)
    e1=(u*u+cm1**2-cm2**2)/(2.*u)
    t1=abs(e1-cm1)
    p1m=sqrt(t1*(t1+2.*cm1))
    p2(1)=p1m*st*cf
    p2(2)=p1m*st*sf
    p2(3)=p1m*ct
    call  rotorq(pins,v,p2,p1)
    p2(1)=-p1(1)
    p2(2)=-p1(2)
    p2(3)=-p1(3)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  tinvuq(p1,p2,u,v,tin1)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculates the total cms energy, velocity of the system of two
!     particles (P1,P2), and the kinetic energy TIN1 of particle "1" in
!     rest system of particle "2"
!
    dimension  v(3),p1(9),p2(9)
! kkg 10/14/03
!      E1=P1(8)+P1(9)
!      E2=P2(8)+P2(9)
    e1=sqrt(p1(4)**2+p1(5)**2+p1(6)**2+p1(9)**2)
    e2=sqrt(p2(4)**2+p2(5)**2+p2(6)**2+p2(9)**2)
    v(1)=(p1(4)+p2(4))/(e1+e2)
    v(2)=(p1(5)+p2(5))/(e1+e2)
    v(3)=(p1(6)+p2(6))/(e1+e2)
    u=(e1+e2)*sqrt(1.-v(1)**2-v(2)**2-v(3)**2)
    tin1=(u**2-(p1(9)+p2(9))**2)/(2.*p2(9))
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine absorpq(partin,ipatin,partne,par1,ne,ie,mv,np,v,u)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Block of calculation of outgoing characteristics
!     in two nucleons (PARTNE,PAR1) absorption of incidente
!     particle (pion) (PARTIN,IPATIN);
!     Calls: TINVUQ, ABELQ
!
    real(real64) ::  masn
    common/memorylaq/pmemo,imemo
    common/tefabs/ tefabs,ehol3,tfr3
    dimension partin(9),ipatin(5),partne(9),par1(9), &
         & paf(9),v(3),pist(3),pnst(3),pmemo(9,5999),imemo(5,5999)
    masn=0.940
    if(ipatin(5).ne.0.and.ipatin(4) == 1)  go  to  13
    paf(4)=partne(4)+par1(4)
    paf(5)=partne(5)+par1(5)
    paf(6)=partne(6)+par1(6)
    paf(9)=2.*par1(9)
    paf(8)=sqrt(paf(4)**2+paf(5)**2+paf(6)**2+paf(9)**2)-paf(9)
    tefabs=paf(8)
    call  tinvuq(partin,paf,u,v,tin1)
!  kkg 03/19/04
13  ctst = 1.-2.*rndm(-1.0_real64)
    fist = 6.283185*rndm(-1.0_real64)
16  call abelq(partin,v,u,pist,pnst,ctst,fist,masn,masn)
    if(mv-5997)17,17,18
18  np = 0
    write(16,19)
19  format (25x,'memorylaq is exceeded in cascad')
    return
17  continue
    pmemo(1,mv+3) = 0.
    pmemo(2,mv+3) = 0.
    pmemo(3,mv+3) = 0.
    pmemo(4,mv+3) = pist(1)
    pmemo(5,mv+3) = pist(2)
    pmemo(6,mv+3) = pist(3)
    pmemo(7,mv+3) = 0.
    pmemo(9,mv+3) = 0.940
    imemo(1,mv+3) = ie
    imemo(2,mv+3) = 0
    imemo(3,mv+3) = 0
    imemo(4,mv+3) = 1
    imemo(5,mv+3) = 0
    pmemo(1,mv+1) = 0.
    pmemo(2,mv+1) = 0.
    pmemo(3,mv+1) = 0.
    pmemo(4,mv+1) = -pist(1)
    pmemo(5,mv+1) = -pist(2)
    pmemo(6,mv+1) = -pist(3)
    pmemo(7,mv+1) = 0.
    pmemo(9,mv+1) = 0.940
    imemo(1,mv+1) = ne
    imemo(2,mv+1) = 0
    imemo(3,mv+1) = 0
    imemo(4,mv+1) = 1
    imemo(5,mv+1) = 0
    np = 2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine partnq(nu,i,p,ip)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Partner selection in the projectile (NU=1) or target (NU=2)
!     nucleus.
!     P and IP refer to the partner nuclear particle which is
!     potentially interacting with the cascade particle.
!
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300) &
         & /hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    dimension p(9),ip(5)

    p(1)=xc(nu,i)
    p(2)=yc(nu,i)
    p(3)=zc(nu,i)
    r = sqrt(p(1)**2+p(2)**2+p(3)**2)
!     R=RPOTEN(R)
    if(nu-1) 10,10,11
10  if(anucl1-10.) 120,120,110
110 tfr=tf01*(((1.+exp(-a1/c1))/(1.+exp((r-a1)/c1)))**0.6666667 )
    go to 12
120 tfr=tf01*exp(-(2./3.)*(r/a1)**2)
    go to 12
11  if(anucl2-10.) 121,121,111
111 tfr = tf02*(((1.+exp(-a2/c2))/(1.+exp((r-a2)/c2)))**.6666667)
    go to 12
121 tfr=tf02*exp(-(2./3.)*(r/a2)**2)
    go to 12
12  continue
!                    !!! 20.06.1995
    if(nu == 1) tfr=tfr*(an1/anucl1)**(2./3.)
    if(nu == 2) tfr=tfr*(an2/anucl2)**(2./3.)
!
    tn = tfr*(rndm(-1.0_real64)**(2./3.))
    pn=sqrt(tn*(tn+1.88))
    ct=1.-2.*rndm(-1.0_real64)
    st=sqrt(1.-ct*ct)
    fi=6.283185*rndm(-1.0_real64)
    cf=cos(fi)
    sf=sin(fi)
    p(4)=pn*st*cf
    p(5)=pn*st*sf
    p(6)=pn*ct
    p(7)=0.
    p(8)=tn
    p(9)=0.94
    ip(2)=0
    ip(3)=0
    ip(4)=1
    ip(1)=iz(nu,i)
    ip(5)=0
    if(ip(1) > 1)  then
       write(16,*) ' partnq: nu,i,ip(1)=',nu,i,ip(1)
       ip(1)=0
    end if
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine directq (v,tin1,mq,mv,np,partin,kp,ith)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   DETERMINING OF DIRECTION OF SECONDARY PARTICLES MOTION.
!   Calls: COSTAQ JTYPAQ ROTORQ, CINEMA
!
    dimension pmemo(9,5999),imemo(5,5999),v(3),pakv(3),partin(9), &
         & plst(3),pl(3),pakst(3),pin(3),b(3),pil(3)
    common/memorylaq/pmemo,imemo
    nd = 0
    kp = 0
    if (mq-1) 12,12,10
10  if (rndm(-1.0_real64) -0.5) 15,11,11
11  if (rndm(-1.0_real64) -0.5) 14,16,16
12  if (rndm(-1.0_real64) -0.5) 14,13,13
13  if (rndm(-1.0_real64) -0.5) 15,16,16
14  m1 = 2
    m2 = 3
    go to 17
15  m1 = 1
    m2 = 3
    go to 17
16  m1 = 1
    m2 = 2
    go to 17
17  lambda=1
    pakv(1)=0.
    pakv(2)=0.
    pakv(3)=0.
    m1temp = mv+m1
    m2temp = mv+m2
18  ltemp = mv+lambda
    if (lambda-m1) 19,21,19
19  if (lambda-m2) 20,21,20
20  ja = jtypaq(ith,mq,lambda)
    ctl = costaq(ja,tin1)
    fl = 6.283185*rndm(-1.0_real64)
    stl = sqrt (1.-ctl**2)
    temp1=cos(fl)
    temp2=sin(fl)
    temp3=pmemo(8,ltemp)
    pmemo(1,ltemp) = 0.
    pmemo(2,ltemp) = 0.
    pmemo(3,ltemp) = 0.
    pmemo(4,ltemp) = temp3*stl*temp1
    pmemo(5,ltemp) = temp3*stl*temp2
    pmemo(6,ltemp) = temp3*ctl
    pmemo(7,ltemp) = 0.
    pakv(1) = pakv(1)+pmemo(4,ltemp)
    pakv(2) = pakv(2)+pmemo(5,ltemp)
    pakv(3) = pakv(3)+pmemo(6,ltemp)
21  if (lambda-np) 22,23,23
22  lambda = lambda+1
    go to 18
23  pakvm = sqrt (pakv(1)**2+pakv(2)**2+pakv(3)**2)
    if (np-3) 25,24,25
24  pil(1)=partin(4)
    pil(2)=partin(5)
    pil(3)=partin(6)
    b(1)=-v(1)
    b(2)=-v(2)
    b(3)=-v(3)
    call  cinema(pil,b,pin,cts,sts,cfs,sfs,tin,partin(9))
    lambda = 1
    go to 27
25  if (pmemo(8,m1temp)-pakvm-pmemo(8,m2temp)) 26,32,32
26  if (pmemo(8,m1temp)-abs(pakvm-pmemo(8,m2temp))) 32,32,24
27  ltemp = mv+lambda
    if (lambda-m1) 29,34,29
28  lambda = lambda+1
    go to 27
29  if (lambda-m2) 30,34,30
30  pl(1) = pmemo(4,ltemp)
    pl(2) = pmemo(5,ltemp)
    pl(3) = pmemo(6,ltemp)
    call rotorq (pin,v,pl,plst)
    pmemo(1,ltemp) = 0.
    pmemo(2,ltemp) = 0.
    pmemo(3,ltemp) = 0.
    pmemo(4,ltemp) = plst(1)
    pmemo(5,ltemp) = plst(2)
    pmemo(6,ltemp) = plst(3)
    pmemo(7,ltemp) = 0.
34  if (lambda-np) 28,31,31
31  call rotorq (pin,v,pakv,pakst)
    ctm1 = (pmemo(8,m2temp)**2-pmemo(8,m1temp)**2-pakvm**2)/ &
         & (2.*pakvm*pmemo(8,m1temp))
    ctm2 = (pmemo(8,m1temp)**2-pmemo(8,m2temp)**2-pakvm**2)/ &
          & (2.*pakvm*pmemo(8,m2temp))
    fm1 = 6.283185*rndm(-1.0_real64)
    fm2 = 3.141592+fm1
    stm1 = sqrt (1.-ctm1**2)
    stm2 = sqrt (1.-ctm2**2)
    cfm1 = cos(fm1)
    sfm1 = sin(fm1)
    cfm2 = cos(fm2)
    sfm2 = sin(fm2)
    pl(1) = pmemo(8,m1temp)*stm1*cfm1
    pl(2) = pmemo(8,m1temp)*stm1*sfm1
    pl(3) = pmemo(8,m1temp)*ctm1
    call rotorq (pakst,v,pl,plst)
    pmemo(1,m1temp) = 0.
    pmemo(2,m1temp) = 0.
    pmemo(3,m1temp) = 0.
    pmemo(4,m1temp) = plst(1)
    pmemo(5,m1temp) = plst(2)
    pmemo(6,m1temp) = plst(3)
    pmemo(7,m1temp) = 0.
    pl(1) = pmemo(8,m2temp)*stm2*cfm2
    pl(2) = pmemo(8,m2temp)*stm2*sfm2
    pl(3) = pmemo(8,m2temp)*ctm2
    call rotorq (pakst,v,pl,plst)
    pmemo(1,m2temp) = 0.
    pmemo(2,m2temp) = 0.
    pmemo(3,m2temp) = 0.
    pmemo(4,m2temp) = plst(1)
    pmemo(5,m2temp) = plst(2)
    pmemo(6,m2temp) = plst(3)
    pmemo(7,m2temp) = 0.
    return
32  nd = nd+1
    if (nd-100) 17,33,33
33  kp = 2
    return
  end
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine pinpnq(rn,r0x,r0y,r0z,mv,delta,nel)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculation of entry point of particle into nucleus and
!     positions of target nucleons;
!     in the case of A+A interactions simulates positions of
!     nucleons in projectile nucleus (see /CENTER/ )
!     in the case of bremsstrahlung gamma + A interaction
!     energy of gamma is samplied from Schiff's spectra (1/E)
!     Calls: RXYZ,PANUN, IDPANUN, ROTNUC
!
!     --Edited by--   --Date--   --Edit Summary--
!     CMJ             3/10/17    Added Error Protection to routine
!
!
!
!
    character(len=8) :: pna1
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2 &
         & /center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
!
    common/xbmax/xbmax,ifib0
    common /stin/ stin,amin
    common/cslid/clider(5999)
    common /idpme/ idpme(5999)
    common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems
    dimension partin(9),ipatin(5)
    t=(2.*rn)**2
    if(anucl1-1.1) 98,98,99
98  rm1=0.
99  if(anucl1 == 2.)  rm1=2.158
    temp1 = t0*(t0+1.88)
    ! call err_chk(1,'laq1.f','1968',2,temp1)
    temp2 = 5.06*sqrt(temp1)
    ! call err_chk(1,'laq1.f','1970',1,temp2)
    r12=rm1+rm2+delta+1./temp2
    call  rxyz(r12,r0x,r0y,r0z)
    if(anucl1-1.) 100,100,101
100 pmemo(1,1)=r0x
    pmemo(2,1)=r0y
    pmemo(3,1)=r0z
    imemo(1,1)=znucl1
    imemo(2,1)=0
!  kkg 10.12.04
    if(amin <= 0.00001)   then
       imemo(2,1)=1   ! gamma
       if(ibrems == 1)  then
!      sampling the energy of gamma from Schiff's spectra (1/E)
          rdm = rndm(-1.0_real64)
          temp1 = tgmin
          ! call err_chk(1,'laq1.f','1986',1,temp1)
          t0 = tgmin*(tgmax/temp1)**rdm
       endif
    endif
!  kkg 10.12.03
    imemo(3,1)=intg(stin)
    imemo(4,1)=intg(anucl1)
    if(anucl1 > 0.5)  imemo(4,1)=1
    imemo(5,1)=0
    pmemo(9,1)=amin
    pmemo(4,1)=0.
    pmemo(5,1)=0.
    temp1 = t0*(t0+2.*pmemo(9,1))
    ! call err_chk(1,'laq1.f','1998',2,temp1)
    pmemo(6,1)=sqrt(temp1)
    pmemo(7,1)=0.
    pmemo(8,1)=t0
    myp(1)=0
    myt(1)=1
    myy(1)=1
    clider(1)=1.
!----> do  1 k=1,9
    do  k=1,9
       partin(k)=pmemo(k,1)
       if(k <= 5)  ipatin(k)=imemo(k,1)
1      continue
    end do
!  kkg 29.10.03
    if(imemo(2,1).ne.0)  then
       id1=10                 ! gamma
    else
       call panun(partin,ipatin,ik1)
       call idpanu(id1,ik1,pna1)
    endif
!  kkg 29.10.03
    idpme(1)=id1
    mv=1
    an1=0.
    zn1=0.
    go  to  102
101 if(anucl1 > 2.1)   go  to  104
    g=1.+t0/0.940
    ct1=1.-2.*rndm(-1.0_real64)
    fi1=6.283185*rndm(-1.0_real64)
    temp1 = 1.-ct1*ct1
    ! call err_chk(1,'laq1.f','2029',2,temp1)
    st1=sqrt(temp1)
    temp1 = g
    ! call err_chk(1,'laq1.f','2032',1,temp1)
    z1=rm1*ct1/temp1
    y1=rm1*st1*sin(fi1)
    x1=rm1*st1*cos(fi1)
    pmemo(1,1)=r0x+x1
    pmemo(2,1)=r0y+y1
    pmemo(3,1)=r0z+z1
    pmemo(1,2)=r0x-x1
    pmemo(2,2)=r0y-y1
    pmemo(3,2)=r0z-z1
103 al=10.*rndm(-1.0_real64)
    br=rndm(-1.0_real64)
    fbr=4.*al*al/(1.+al*al)**2
    if(br > fbr)   go  to  103
    epm=0.940*0.00223
    temp1 = epm
    ! call err_chk(1,'laq1.f','2048',2,temp1)
    pm1=sqrt(temp1)*al
    ct1=1.-rndm(-1.0_real64)
    temp1 = 1.-ct1*ct1
    ! call err_chk(1,'laq1.f','2051',2,temp1)
    st1=sqrt(temp1)
    fi1=6.283185*rndm(-1.0_real64)
    px1=pm1*st1*cos(fi1)
    py1=pm1*st1*sin(fi1)
    pz1=pm1*ct1
    pmemo(4,1)=px1
    pmemo(4,2)=-px1
    pmemo(5,1)=py1
    pmemo(5,2)=-py1
    temp1 = pm1*pm1+0.94*0.94
    temp2 = g*g-1.
    ! call err_chk(1,'laq1.f','2063',2,temp1*temp2)
    pmemo(6,1)=g*pz1+sqrt(temp1*temp2)
    pmemo(6,2)=pmemo(6,1)-2.*g*pz1
    pmemo(7,1)=0.
    pmemo(7,2)=0.
    ! call err_chk(1,'laq1.f','2070',2,temp1)
    ! call err_chk(1,'laq1.f','2070, 2071',2,temp2)
    pmemo(8,1)=g*sqrt(temp1)+pz1*sqrt(temp2)-0.940
    pmemo(8,2)=pmemo(8,1)-2.*pz1*sqrt(temp2)
    pmemo(9,1)=0.940
    pmemo(9,2)=0.940
    imemo(1,1)=0
    idpme(1)=1220
    imemo(2,1)=0
    imemo(3,1)=0
    imemo(4,1)=1
    imemo(1,2)=1
    idpme(2)=1120
    imemo(2,2)=0
    imemo(3,2)=0
    imemo(4,2)=1
    imemo(5,1)=0
    imemo(5,2)=0
    myp(1)=0
    myp(2)=0
    myt(1)=1
    myy(1)=1
    myt(2)=1
    myy(2)=1
    clider(1)=1.
    clider(2)=1.
    mv=2
    an1=0
    zn1=0
    go  to  102
104 continue
    if(nel.ne.0)  call  rotnuc(1,anucl1)
    if(nel.ne.0)  go to 102
    im=intg(anucl1)
    nz1=intg(znucl1)
!----> do 21 i=1,im
    do i=1,im
11     b1=rndm(-1.0_real64)
       ri=rm1*(b1**(1./3.))
       if(anucl1-10.) 12,12,13
12     temp1 = a1**2
       ! call err_chk(1,'laq1.f','2109',1,temp1)
       fb=exp(-(ri**2)/temp1)
       go to 14
13     temp1 = c1
       ! call err_chk(1,'laq1.f','2113',1,temp1)  ! no need to error check denominator below
       fb=(1.+exp(-a1/temp1))/(1.+exp((ri-a1)/temp1))
14     drnd=rndm(-1.0_real64)
       if(drnd-fb) 15,15,11
15     ct=1.-2.*rndm(-1.0_real64)
       fi=6.283185*rndm(-1.0_real64)
       temp1 = 1.-ct**2
       ! call err_chk(1,'laq1.f','2120',2,temp1)
       st=sqrt(temp1)
       rs=ri*st
       xc(1,i)=rs*cos(fi)
       yc(1,i)=rs*sin(fi)
       zc(1,i)=ri*ct
       if(i-nz1) 16,16,17
16     iz(1,i)=1
       go to 18
17     iz(1,i)=0
18     if(i-1) 221,221,19
19     km=i-1
!----> do 20 k=1,km
       do k=1,km
          if((xc(1,i)-xc(1,k))**2+(yc(1,i)-yc(1,k))**2+(zc(1,i)-zc(1,k))**2 &
               & -t) 11,20,20
20        continue
       end do
221    continue
21     continue
    end do
102 if(anucl2 > 1.1)   go  to  114
    pmemo(1,1)=0.
    pmemo(2,1)=0.
    pmemo(3,1)=1.e-20
    pmemo(4,1)=0.
    pmemo(5,1)=0.
    pmemo(6,1)=1.e-20
    pmemo(7,1)=0.
    pmemo(8,1)=0.
    pmemo(9,1)=0.940
    if(anucl2 < 0.5)  pmemo(9,1)=0.140
    imemo(1,1)=znucl2
    idpme(2)=1120
    if(znucl2 <= 0.5) idpme(2)=1220
    imemo(2,1)=0
    imemo(3,1)=0
    imemo(4,1)=anucl2
    imemo(5,1)=0
    clider(1)=1.
    an2=0.
    zn2=0.
    mv=1
    myp(1)=1
    myt(1)=0
    myy(1)=1
    return
114 continue
    if(nel.ne.0)  go to 105
    im=intg(anucl2)
    nz2=intg(znucl2)
!----> do 33 i=1,im
    do i=1,im
22     b1=rndm(-1.0_real64)
       ri=rm2*(b1**(1./3.))
       if(anucl2-10.) 23,23,24
23     temp1 = a2**2
       ! call err_chk(1,'laq1.f','2173',1,temp1)
       fb=exp(-(ri**2)/temp1)
       go to 25
24     temp1 = c2
       ! call err_chk(1,'laq1.f','2177',1,temp1)
       fb=(1.+exp(-a2/temp1))/(1.+exp((ri-a2)/temp1))
25     drnd=rndm(-1.0_real64)
       if(drnd-fb) 26,26,22
26     ct=1.-2.*rndm(-1.0_real64)
       fi=6.283185*rndm(-1.0_real64)
       temp1 = 1.-ct**2
       ! call err_chk(1,'laq1.f','2184',2,temp1)
       st=sqrt(temp1)
       rs=ri*st
       xc(2,i)=rs*cos(fi)
       yc(2,i)=rs*sin(fi)
       zc(2,i)=ri*ct
       if (i-nz2) 27,27,28
27     iz(2,i)=1
       go to 29
28     iz(2,i)=0
29     if(i-1) 33,33,30
30     km=i-1
!----> do 32 k=1,km
       do k=1,km
          if((xc(2,i)-xc(2,k))**2+(yc(2,i)-yc(2,k))**2+(zc(2,i)-zc(2,k))**2 &
               & -t) 22,32,32
32        continue
       end do
33     continue
    end do
    return
105 call  rotnuc(2,anucl2)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine rotnuc(n,an)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     random rotation of projectile (N=1) or target nucleus (N=2)
!     changing the positions of nucleons
!     Calls: ROTORQ
!
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    dimension a(3),b(3),r(3),rr(3)
    if(an < 2.1)   return
    ct=1.-2.*rndm(-1.0_real64)
    st=sqrt(1.-ct*ct)
    fi=6.283185*rndm(-1.0_real64)
    a(1)=st*cos(fi)
    a(2)=st*sin(fi)
    a(3)=ct
    b(1)=0.
    b(2)=0.
    b(3)=1.
    ia=an+0.1
!----> do 10 i=1,ia
    do i=1,ia
       r(1)=xc(n,i)
       r(2)=yc(n,i)
       r(3)=zc(n,i)
       call rotorq(a,b,r,rr)
       xc(n,i)=rr(1)
       yc(n,i)=rr(2)
       zc(n,i)=rr(3)
10     continue
    end do
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!  KKG 09/05/08 : parameter IRET was added
  subroutine   typnew(partin,ipatin,partne,ipatne,v,u,tin1, &
       & sabs,mv,np,nabs,par1,ipa1,n3,nu,n2,na1,na2,iret)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!
!     Determine interaction type and calculate
!     secondary particles' characteristics.
!   Called by: CASCAW
!   Calls: CINEMA, PANUID, SLQEKQ, PARTNQ, CENUM1, CROSEG, ABSORPQ,
!          TFERMIQ, ELEXQ, PINETA, PINSEX, BBSEX, PIPIKK,AKANNI,
!          PIYSEX, BINELQ, HEINEN
!
!   Edited by CMJ (3/10/17): Added minor error protection
!
    character(len=4) :: atyp(5)
    character(len=8) :: dtyp(11),tdiag,pna1,pna2,pnak
    logical yesela
    common/yesela/ yesela
    common/ncasca/ncas,ncpri
    common/inttyp/ ityp
    common/ithea/ithea(11)
    common/lowmis/lowmis
    common/iact/ iact
    common /memorylaq/ pme(9,5999),ime(5,5999)
    common /hele/ ihele
    common /idn12/ id1,id2
    common /idn120/ id10,id20
    common /idpme/ idpme(5999)
    common/comelx/ sel
    common/comcro/ sto
    common /data2/ pud,ps1,sigma,cx2
    common /nrapi/ nrapi
    common /sori/ sori(5999),ssor
    common /hadr1/hadr1(4,5999),hadr2(4,5999),hadi1(4),hadi2(4)
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/tefabs/ tefabs,ehol3,tfr3

    common /comenb/ enbou
    common/uold/ uold
    common/strexc/strexc

    dimension partin(9),ipatin(5),partne(9),ipatne(5),v(3), &
         & par1(9),ipa1(5),pa(9),ipa(5),ps(3),pl(3)
    data atyp/'abs:','ele:','bin:','hei:','sex:'/
    data dtyp/'diftri','planar','uncyli','elast  ','annih','difsma', &
         & 'chains','binar ','??????','regtri ','doubdi'/
!
!
!  KKG 09/05/08
    iret=0
    ntry=0
    strexc = 0.0                 !   kkg 19.02.07
    am1=partin(9)
    ssor=u
!  for a 'projectile' hadron
    do j=1,3
       ps(j)=partin(3+j)
    enddo
    if(nrapi == 2)         then
       call  cinema(ps,vpr,pl,ctl,stl,cfl,sfl,tl,partin(9))
    elseif(nrapi == 3)     then
       call  cinema(ps,vta,pl,ctl,stl,cfl,sfl,tl,partin(9))
    else
       do j=1,3
          pl(j)=partin(3+j)
       enddo
    endif
    hadi1(1)=pl(1)
    hadi1(2)=pl(2)
    hadi1(3)=pl(3)
    hadi1(4)=partin(9)
!  for a 'target-partner' hadron
    do j=1,3
       ps(j)=partne(3+j)
    enddo
    if(nrapi == 2)         then
       call  cinema(ps,vpr,pl,ctl,stl,cfl,sfl,tl,partne(9))
    elseif(nrapi == 3)     then
       call  cinema(ps,vta,pl,ctl,stl,cfl,sfl,tl,partne(9))
    else
       do j=1,3
          pl(j)=partne(3+j)
       enddo
    endif
    hadi2(1)=pl(1)
    hadi2(2)=pl(2)
    hadi2(3)=pl(3)
    hadi2(4)=partne(9)
11  continue
!  KKG 09/05/08
    ntry=ntry + 1
    if(ntry > 100)  then
       iret=1
       return
    end if
!
    enbous=enbou
    uolds =uold
    if((ipatin(4)+ipatne(4)) == 1)  then
       enbou=10.0
       uold =1.3
    else
       enbou=2.0
       uold =2.0
    endif
!
    id10=id1
    id20=id2
    lowmis=0
    ityp=0
    nabs=0
    yesela=.true.
    call panuid(id1,ik1,pna1)
    call panuid(id2,ik2,pna2)
    if(ncas >= ncpri) then
       write(16,599) pna1,id1,pna2,id2,tin1,sabs
       write( *,599) pna1,id1,pna2,id2,tin1,sabs
599    format(1x,a4,'(',i5,')','+',a4,'(',i5,')','  tin1=',f7.3, &
            & ' sabs=',1pe9.3)
    endif
    if(tin1 <= 0) then
!  kkg 10/14/03
!       write(16,*) 'PI=', (PARTIN(KK),KK=4,9)
!       write(16,*) 'PN=', (PARTNE(KK),KK=4,9)
    endif
!       for gamma + N interaction, KKG, 12/10/04
98  continue
    if(id1 == 10)  then
       sto = csgntot(ipatne(1),tin1,partne(9))/1000.
       temp1 = sto+sabs
       ! call err_chk(1,'laq1.f','2370',1,temp1)
       betabs = sabs/(temp1)
       if(rndm(-1.0_real64) <= betabs) then
          go  to  101
       else
          go  to  13
       endif
    endif
!
112 continue
    if(ik2 < 37.or.ik2 > 38)           go  to  18
!  kkg 10/28/03
    if(ik1 == 1.or.ik1 == 2.or.ik1 == 7.or.id1 == 10) &
         & go  to  99
    if(ik1 < 37.or.ik1 > 38)           go  to  18
    go  to  13
99  continue
    call slqekq (l,ms,mq,ksi,me,ipatin,ipatne)
    sto=croseg(l,ms,mq,ksi,0,tin1,am1,ipatin(5))
    temp1 = sto+sabs
    ! call err_chk(1,'laq1.f','2390',1,temp1)
    betabs=sabs/(temp1)
    drnd=rndm(-1.0_real64)
    if(drnd-betabs) 101,101,13
101 continue
    dl=0.
    if(nu == 1)  call  cenum1(na1,partin,ipatin(4),dl,nc,n2,n3,n2,1)
    if(nu == 2)  call  cenum1(na2,partin,ipatin(4),dl,nc,n2,n3,n2,2)
    if(n3 == 0)   then
       write(*,*) ' typnew: n3=0,id1=',id1
       go  to  13
    endif
    call  partnq(nu,n3,par1,ipa1)
!  kkg 03/19/04   only gamma +(pn) !
!     if(ID1 == 10.and.(IPATNE(1)+IPA1(1)).ne.1)  go to 13
    tfr3=tfermiq(par1(1),par1(2),par1(3),nu)
    ehol3=par1(8)
    par1(9) = partne(9)
    ie2=ipatne(1)
    ie3=ipa1(1)
    ie1=ipatin(1)
104 if(ie1) 108,103,105
103 ne1=ie2
    ne2=ie3
    go to 14
105 if(ie2+ie3-1) 107,106,13
106 ne1=1
    ne2=1
    go to 14
107 ne1=1
    ne2=0
    go to 14
108 if(ie2+ie3-1) 13,109,110
109 ne1=0
    ne2=0
    go to 14
110 ne1=0
    ne2=1
14  continue
    call absorpq(partin,ipatin,partne,par1,ne1,ne2,mv,np,v,u)
    nabs = 1
    ityp=1
    go  to  27
13  continue
!      for gamma + N ==> hadrons, KKG, 12/10/04
    if(id1 == 10)  then
       if(ncas >= ncpri) write(*,*) ' to gntoh: tin1=',tin1
       call  gntoh(v,u,tin1,partin,ipatin,partne,ipatne,mv,np)
       if(ncas >= ncpri) write(*,*) ' from gntoh: np=',np
       if(np < 2)  go to 101
       ityp=5
       nin=0
       go  to  27
    endif
!
!    ****** TIN1=4.5  IS INTRODUCED ******* Kostya
!      IF(IACT.GE.2.AND.TIN1.GT.4.5)    GO  TO  18
    if(iact >= 2.and.tin1 > 0.8)    go  to  18   ! kkg 21.11.07
!      IF(IACT.GE.2.AND.U.GT.Uold)        GO  TO  18 ! KKG 04.03.07

    call slqekq (l,ms,mq,ksi,me,ipatin,ipatne)
    sto=croseg(l,ms,mq,ksi,0,tin1,am1,ipatin(5))
    sel=croseg(l,ms,mq,ksi,1,tin1,am1,ipatin(5))
    sex=croseg(l,ms,mq,ksi,2,tin1,am1,ipatin(5))
    temp1 = sto
    ! call err_chk(1,'laq1.f','2455',1,temp1)
    betael=(sel+sex)/temp1
    drnd=rndm(-1.0_real64)
    if (drnd-betael) 16,16,15
16  continue
    call elexq(v,u,tin1,partin,ipatin,ipatne,mv,np,l,ms, &
         & mq,ksi,me)
    ityp=2
    go  to  27
15  continue
!      IF(TIN1.LE.5.0d0)  GO  TO  20
    if(tin1 <= 0.8d0)  go  to  20      ! kkg 21.11.07

    yesela=.false.
    go  to  18
20  continue
!
    yesela=.false.

!     kkg  11/21/07
    ieta1=0
    ieta2=0
    isex1=0
    isex2=0
    ibbkak=0
    ipibkak=0

    call pineta(id1,id2,v,u,partin,ipatin,partne,ipatne, &
         & mv,np,ieta1)
    if((id1 == 1120.and.id2 == 1120).or. &
         & (id1 == 1220.and.id2 == 1220).or. &
         & (id1 == 1120.and.id2 == 1220).or. &
         & (id1 == 1220.and.id2 == 1120).and.u <= uold) & ! kkg 04.03.07 &
         & call pneta(id1,id2,v,u,tin1,partin,partne,mv,np,ieta2)

    if(ieta1 == 0) &
         & call pinsex(id1,id2,v,u,partin,ipatin,partne,ipatne, &
         & mv,np,isex1)
    if(ipatin(4) == 1.and.ipatin(3) == 0.and. &
         & ipatne(4) == 1.and.ipatne(3) == 0.and.u <= uold.and. &
         & ieta2 == 0) call bbsex(id1,id2,v,u,partin,partne,mv,np,isex2)
    if(ieta2 == 0.and.isex2 == 0.and. &
         & ipatin(4) == 1.and.ipatne(4) == 1) &
         & call  bbkak(id1,id2,v,u,partin,partne,ipatin,ipatne,mv,np, &
         & ibbkak)
    if(ibbkak.ne.0)   then
       ityp=5
       nin=0
       go to  27
    endif
    if(isex1 == 0.and.ieta1 == 0) &
         & call  pibkak(id1,id2,v,u,partin,partne,ipatin,ipatne,mv,np, &
         & ipibkak)
    if(ipibkak.ne.0)   then
       ityp=5
       nin=0
       go to  27
    endif
!
    if(isex1.ne.0.or.isex2.ne.0.or.ieta1.ne.0.or.ieta2.ne.0) &
         & then
       ityp=5
       nin=0
    else
       call binelq(partin,ipatin,ipatne,l,ms,mq,ksi,me,v,u,tin1, &
            & mv,np,nin)
       ityp=3
       if (nin.ne.0)   go  to  11
    endif
    go  to   27
18  continue
!
!     kkg  07/25/06
    ikak=0
    ieta1=0
    ieta2=0
    isex1=0
    isex2=0
    isex3=0
    isex4=0
    ietan=0
    ibbkak=0
    ipibkak=0

    call pipikk(id1,id2,v,u,partin,ipatin,partne,ipatne, &
         & mv,np,ikak)
    call pineta(id1,id2,v,u,partin,ipatin,partne,ipatne, &
         & mv,np,ieta1)
    if(ieta1 == 0) &
         & call pinsex(id1,id2,v,u,partin,ipatin,partne,ipatne, &
         & mv,np,isex1)
    if((id1 == 1120.and.id2 == 1120).or. &
         & (id1 == 1220.and.id2 == 1220).or. &
         & (id1 == 1120.and.id2 == 1220).or. &
         & (id1 == 1220.and.id2 == 1120).and.u <= uold) & ! kkg 04.03.07 &
         & call pneta(id1,id2,v,u,tin1,partin,partne,mv,np,ieta2)
    if(ieta1 == 0.and.isex1 == 0) &
         & call piysex(id1,id2,v,u,partin,partne,mv,np,isex2)
    if(ipatin(4) == 1.and.ipatin(3) == 0.and. &
!    &   IPATNE(4).EQ.1.AND.IPATNE(3).EQ.0.and.U <= Uold.
         & ipatne(4) == 1.and.ipatne(3) == 0.and.u <= 2.0d0 &
         & .and.ieta2 == 0) &
         & call bbsex(id1,id2,v,u,partin,partne,mv,np,isex3)
!
    call akanni(id1,id2,v,u,partin,mv,np,isex4)
!   kkg  20.09.06
    if(id1 == 220) call etan(id1,id2,v,u,partin,ipatin, &
         & partne,ipatne,mv,np,ietan)
    if(ieta2 == 0.and.isex3 == 0.and. &
         & ipatin(4) == 1.and.ipatne(4) == 1) &
         & call  bbkak(id1,id2,v,u,partin,partne,ipatin,ipatne,mv,np, &
         & ibbkak)
    if(ibbkak.ne.0)   then
       ityp=5
       nin=0
       go to  27
    endif
    if(isex1 == 0.and.ieta1 == 0) &
         & call  pibkak(id1,id2,v,u,partin,partne,ipatin,ipatne,mv,np, &
         & ipibkak)
    if(ipibkak.ne.0)   then
       ityp=5
       nin=0
       go to  27
    endif
    if(ikak.ne.0.or.ieta1.ne.0.or.ieta2.ne.0.or. &
         & isex1.ne.0.or.isex2.ne.0.or.isex3.ne.0.or. &
         & isex4.ne.0.or.ietan.ne.0) then

       ityp=5
       nin=0
    else
!                            12.03.94
       ps1=0.750
!
       call  heinen(partin,ipatin,partne,ipatne,mv,np,nin)
       ityp=4
!
!----> do  19  i=1,11
       do   i=1,11
19        if(ithea(i).ne.0) tdiag=dtyp(i)
       end do
       if(nin.ne.0)  go  to  11
    endif
27  continue
!
    enbou=enbous
    uold =uolds
!
    if(ityp == 0) go  to  603
    if(ncas >= ncpri.or.lowmis.ne.0)  then
!      IF(ITYP.NE.4)
!     *write(16,601) ATYP(ITYP),      PARTIN,IPATIN,PARTNE,IPATNE
! IF(ITYP.EQ.1) WRITE(16,600)    PAR1,IPA1
!      IF(ITYP.EQ.4.AND.IHELE.EQ.1)
!     *write(16,611) ATYP(ITYP),TDIAG,PARTIN,IPATIN,PARTNE,IPATNE
!      IF(ITYP.EQ.4.AND.IHELE.EQ.2)
!     *write(16,612) ATYP(ITYP),TDIAG,PARTIN,IPATIN,PARTNE,IPATNE
       do    k=1,np
          m=mv+k
          if(np == 2.and.k == 2.and.(iabs(ipatin(4))+iabs(ipatne(4))) > 0) &
               & m=m+1
!      write(16,602) M,(PME(I,M),I=4,6),PME(9,M),(IME(J,M),J=1,5)
       enddo
    endif
!
    if(ityp == 4)  go  to  3
!   kkg 12/13/04
    if(id1 == 10.and.ityp == 5)  go  to 3
!----> do  2  mm=1,np
    do   mm=1,np
       m=mm
       if(np == 2.and.m == 2.and.(iabs(ipatin(4))+iabs(ipatne(4))) > 0) &
            & m=3
!----> do  1 k=1,9
       do  k=1,9
          pa(k)=pme(k,mv+m)
          if(k <= 5)  ipa(k)=ime(k,mv+m)
1         continue
       end do
       call panun(pa,ipa,ikk)
       call idpanu(idk,ikk,pnak)
       idpme(mv+m)=idk
2      continue
    end do
3   continue
    return
600 format(1x,4x,9(1x,f10.3),4i2,i10)
601 format(1x,a4,9(1x,f10.3),4i2,i10/ &
         & 5x,9(1x,f10.3),4i2,i10)
611 format(1x,a4,1x,a8,'(le)'/ &
         & 5x,9(1x,f10.3),4i2,i10/ &
         & 5x,9(1x,f10.3),4i2,i10)
612 format(1x,a4,1x,a8,'(he)'/ &
         & 5x,9(1x,f10.3),4i2,i10/ &
         & 5x,9(1x,f10.3),4i2,i10)
602 format(1x,i3,4(1x,f10.3),4i2,i10)
603 write(16,604)
604 format(2x,'typnew: ityp=0 !!!!!'/)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  akanni(id1,id2,v,u,pin,mv,np,isex)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     charge exchange and annihilation of K- or AK0
!     Calls: ABELQ
!
    real(real64) ::  mpi,mn,mk,ml,ms,m1,m2
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension v(3),pin(9),p1(3),p2(3)
    data  mpi/0.140/,mn/0.940/,mk/0.494/,ml/1.116/,ms/1.189/
    isex=0
    ir=0
    if((id1 == -130.and.id2 == 1120).or.(id1 == 1120.and.id2 == -130)) &
         & ir=1  ! k- + p
    if((id1 == -130.and.id2 == 1220).or.(id1 == 1220.and.id2 == -130)) &
         & ir=2  ! k- + n
    if((id1 == -230.and.id2 == 1120).or.(id1 == 1120.and.id2 == -230)) &
         & ir=3  ! ak0 + p
    if((id1 == -230.and.id2 == 1220).or.(id1 == 1220.and.id2 == -230)) &
         & ir=4  ! ak0 + n
    if(ir == 0)        return
!
    if(u <= (mk+mn))   return
    ek0=(u**2-mn**2-mk**2)/(2.*mn)
    pk0=sqrt(ek0**2-mk**2)
!   Calculation of K- + p=>L + pi0 cross section
!   param. of G.Q.Li et al NP A625(1997)342
    if(pk0 <= 0.6)        then
       sla=1.205*pk0**(-1.428)
    elseif(pk0 <= 1.0)    then
       sla=3.5*pk0**0.659
    else
       sla=3.5*pk0**(-3.97)
    endif
!   Calculation of K- + p=>S0 + pi0 cross section
!   param. of G.Q.Li et al NP A625(1997)342
    if(pk0 <= 0.345)      then
       ssi=0.624*pk0**(-1.830)
    elseif(pk0 <= 0.425)  then
       ssi=0.0138/((pk0-0.385)**2+0.0017)
    else
       ssi=0.7*pk0**(-2.09)
    endif
!   Calculation of K- + p total,elastic,
!   and charge exchange cross sections
!   param. of G.Q.Li et al NP A625(1997)342
    if(pk0 <= 0.35)      then
       sto=23.5*pk0**(-1.04)
    elseif(pk0 <= 0.46)  then
       sto=0.504/((pk0-0.39)**2+0.0056)
    elseif(pk0 <= 1.05)  then
       sto=181.9*(pk0-0.75)**2+34.0
    else
       sto=55.2*pk0**(-1.85)
    endif
    if(pk0 <= 0.70)      then
       sel=11.2*pk0**(-0.986)
    else
       sel=5.0/((pk0-0.95)**2+0.25)
    endif
    if(pk0 <= 0.35)      then
       sex=1.813*pk0**(-1.14)
    elseif(pk0 <= 0.43)  then
       sex=0.0192/((pk0-0.39)**2+0.0016)
    else
       sex=15.9/((pk0-0.9)**2+2.65)
    endif
    san=sla+ssi
    ban=(san+sex)/sto
    if(rndm(-1.0_real64) > ban)  return
    bex=sex/(sex+san)
    if(rndm(-1.0_real64) <= bex)  then
!  charge exchange K- p=>AK0 n or AK0 n=> K- p
       is1=-1
       is2=0
       m1=mk
       m2=mn
       if(ir == 1)      then
! K- p=> AK0 n
          ie1=0
          ie2=0
       elseif(ir == 4)  then
! AK0 n=> K- p
          ie1=-1
!  kkg 07/26/06
!         IE2=0
          ie2=1
       else
          return
       endif
    else
!  annihilation: AK+N==>PI+Y
       is1=0
       is2=-1
       m1=mpi
       bla=sla/(sla+ssi)
       if(rndm(-1.0_real64) <= bla)  then
          m2=ml
          ie2=0
          if(ir == 1)      then
! K- p=> pi0 L
             ie1=0
          elseif(ir == 2)  then
! K- n=> pi- L
             ie1=-1
          elseif(ir == 3)  then
! AK0 p=> pi+ L
             ie1=1
          else
! AK0 n=> pi0 L
             ie1=0
          endif
       else
          m2=ms
          br=rndm(-1.0_real64)
          if(ir == 1)      then
             if(br <= 0.333333)     then
! K- p=> pi0 S0
                ie1=0
                ie2=0
             elseif(br >= 0.666667) then
! K- p=> pi- S+
                ie1=-1
                ie2= 1
             else
! K- p=> pi+ S-
                ie1= 1
                ie2=-1
             endif
          elseif(ir == 2)  then
             if(br <= 0.5)          then
! K- n=> pi- S0
                ie1=-1
                ie2= 0
             else
! K- n=> pi0 S-
                ie1= 0
                ie2=-1
             endif
          elseif(ir == 3)  then
             if(br <= 0.5)          then
! AK0 p=> pi0 S+
                ie1= 0
                ie2= 1
             else
! AK0 p=> pi+ S0
                ie1= 0
                ie2= 1
             endif
          else
             if(br <= 0.5)          then
! AK0 n=> pi0 S0
                ie1= 0
                ie2= 0
             else
! AK0 n=> pi- S+
                ie1=-1
                ie2= 1
             endif
          endif
       endif
    endif

!
    ct=2.*rndm(-1.0_real64)-1.
    fi=6.283185*rndm(-1.0_real64)
    call  abelq(pin,v,u,p1,p2,ct,fi,m1,m2)
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p1(1)
    pme(5,mv+3)=p1(2)
    pme(6,mv+3)=p1(3)
    pme(7,mv+3)=0.
    pme(9,mv+3)=m1
    ime(1,mv+3)=ie1
    ime(2,mv+3)=0
    ime(3,mv+3)=is1
    ime(4,mv+3)=0
    ime(5,mv+3)=0
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=-p1(1)
    pme(5,mv+1)=-p1(2)
    pme(6,mv+1)=-p1(3)
    pme(7,mv+1)=0.
    pme(9,mv+1)=m2
    ime(1,mv+1)=ie2
    ime(2,mv+1)=0
    ime(3,mv+1)=is2
    ime(4,mv+1)=1
    ime(5,mv+1)=0
    isex=1
    np = 2
!       write(*,*)  ' akanni: IR=',IR,' ID1,ID2=',ID1,ID2
!       write(*,*)  'IE1,IS1,IE2,IS2=',IE1,IS1,IE2,IS2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  etan(id1,id2,v,u,pin,iin,pn,ipn,mv,np,ietan)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    real(real64) ::  v(3),u,pin(9),pn(9),p1(3),p2(3),mpi,mn,meta,m3
    logical yesela
    common/yesela/ yesela
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension  iin(5),ipn(5)
    data  mpi/0.140/,mn/0.940/,meta/0.549/
    ietan=0
    if(id1 == 220.and.(id2 == 1120.or.id2 == 1220)) go  to  10
    return
10  continue
    ksi=1
    if(ipn(1) == 0) ksi=2
    tetn=(u**2-(pin(9)+pn(9))**2)/(2.*pn(9))
    stot=sigeta(ksi,0,tetn)
    sela=sigeta(ksi,1,tetn)
    rns=rndm(-1.0_real64)
    if (rns <= sela/stot)  then
       m3=meta              ! elastic scattering
       ie3=0
       ie1=ipn(1)
       iel=1
    else
       m3=mpi
       ie1=1
       if(rndm(-1.0_real64) > 0.5)  ie1=0
       ie3=ipn(1)-ie1
       iel=0
    endif
    ct=1.0 - 2.0*rndm(-1.0_real64)
    fi=6.283185*rndm(-1.0_real64)
    call  abelq(pin,v,u,p1,p2,ct,fi,m3,mn)
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p1(1)
    pme(5,mv+3)=p1(2)
    pme(6,mv+3)=p1(3)
    pme(7,mv+3)=0.
    pme(9,mv+3)=m3
    ime(1,mv+3)=ie3
    ime(2,mv+3)=0
    ime(3,mv+3)=0
    ime(4,mv+3)=0
    if(iel == 1)  then
       ime(5,mv+3)=idint(1000.0d0 *taun(8))
       if(ime(5,mv+3) < 1) ime(5,mv+3)=1
    else
       ime(5,mv+3)=0
    endif
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=-p1(1)
    pme(5,mv+1)=-p1(2)
    pme(6,mv+1)=-p1(3)
    pme(7,mv+1)=0.
    pme(9,mv+1)=mn
    ime(1,mv+1)=ie1
    ime(2,mv+1)=0
    ime(3,mv+1)=0
    ime(4,mv+1)=1
    ime(5,mv+1)=0
    ietan=2-iel
    np = 2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  bbkak(ida,idb,v,u,pin,pn,iin,ipn,mv,np,ibbkak)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    real(real64) ::  v(3),u,pin(9),pn(9),mk,mn
    logical yesela
    common/yesela/ yesela
    common/memorylaq/pme(9,5999),ime(5,5999)
    common /idpme/ idpme(5999)
    common/ncasca/ncas,ncpri
    dimension iin(5),ipn(5),amas4(4),ps4(5,4)
    data  mk/0.494d0/,mn/0.939d0/,zro/0.0d0/, &
         & onth/0.33333333d0/,twth/0.66666667d0/, &
         & onfo/0.25d0/,ontw/0.5d0/,thfo/0.75d0/
    ibbkak = 0
    umin=2.0*mn + 2.0*mk
    if(iin(4) == 1.and.iin(3) == 0.and.iin(2) == 0.and. &
         & ipn(4) == 1.and.ipn(3) == 0.and.ipn(2) == 0.and. &
         & u > umin) then
       go  to  1
    else
       return
    endif
1   continue
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sino=sito
    if(yesela)         go to 2
    call  crosec(0,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siel,0)
    call  crosec(2,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siex,0)
    sino=sito-siel-siex
2   if(sino <= 0.)                     return
    ru = (umin/u)**2
    skk = 1.5*((1.0d0 - ru)**3.17)*ru**1.96
    if(rndm(-1.0_real64) > skk/sino)  return
    if(iin(5) == 0.and.ipn(5) == 0)  then
!  N + N
       if(iin(1) == 1.and.ipn(1) == 1)  then
!  p + p
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       elseif((iin(1)+ipn(1)) == 1)  then
!  p + n
          rnd=rndm(-1.0_real64)
          if(rnd <= onfo)          then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= ontw)      then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          elseif(rnd <= thfo)      then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          else
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       else
!  n + n
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          endif
       endif
    elseif((iin(5).ne.0.and.ipn(5) == 0).or. &
            & (iin(5) == 0.and.ipn(5).ne.0))  then
!  D + N
       if((iin(1) == 2.and.ipn(1) == 1).or. &
            & (iin(1) == 1.and.ipn(1) == 2))  then
! D++ + p
          ime(1,mv+1) = 1    ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1    ! p
          idpme(mv+2) = 1120
          ime(1,mv+3) = 1    ! k+
          idpme(mv+3) = 130
          ime(1,mv+4) = 0    ! ak0
          idpme(mv+4) =-230
       elseif((iin(1) == 2.and.ipn(1) == 0).or. &
               & (iin(1) == 0.and.ipn(1) == 2))  then
! D++ + n
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       elseif(iin(1) == 1.and.ipn(1) == 1)  then
! D+ + p
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       elseif((iin(1)+ipn(1)) == 1)  then
! D+  + n or D0 + p
          rnd=rndm(-1.0_real64)
          if(rnd <= onfo)          then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= ontw)      then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          elseif(rnd <= thfo)      then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          else
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       elseif((iin(1)+ipn(1)) == 0)  then
! D0  + n or D- + p
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          endif
       elseif((iin(1) == -1.and.ipn(1) == 0).or. &
               & (iin(1) == 0.and.ipn(1) == -1))  then
! D- + n
          ime(1,mv+1) = 0    ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0    ! n
          idpme(mv+2) = 1220
          ime(1,mv+3) = 0    ! k0
          idpme(mv+3) = 230
          ime(1,mv+4) =-1    ! k-
          idpme(mv+4) =-130
       else
          write( *,*) ' bbkak1 d+n : ida,idb=',ida,idb
          write(16,*) ' bbkak1 d+n : ida,idb=',ida,idb
          return
       endif
    elseif(iin(5).ne.0.and.ipn(5).ne.0)  then
!  D + D
       if((iin(1)+ipn(1)) == 3)    then
!  D++  + D+
          ime(1,mv+1) = 1    ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1    ! p
          idpme(mv+2) = 1120
          ime(1,mv+3) = 1    ! k+
          idpme(mv+3) = 130
          ime(1,mv+4) = 0    ! ak0
          idpme(mv+4) =-230
       elseif((iin(1)+ipn(1)) == 2)  then
!  D++  + D0  or  D+  +  D+
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       elseif((iin(1)+ipn(1)) == 1)  then
!  D+  + D0  or  D-  +  D++
          rnd=rndm(-1.0_real64)
          if(rnd <= onfo)          then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= ontw)      then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          elseif(rnd <= thfo)      then
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 1    ! p
             idpme(mv+2) = 1120
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          else
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          endif
       elseif((iin(1)+ipn(1)) == 0)  then
!  D-  + D+  or  D0  +  D0
          rnd=rndm(-1.0_real64)
          if(rnd <= onth)  then
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 1    ! k+
             idpme(mv+3) = 130
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          elseif(rnd <= twth)  then
             ime(1,mv+1) = 0    ! n
             idpme(mv+1) = 1220
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) = 0    ! ak0
             idpme(mv+4) =-230
          else
             ime(1,mv+1) = 1    ! p
             idpme(mv+1) = 1120
             ime(1,mv+2) = 0    ! n
             idpme(mv+2) = 1220
             ime(1,mv+3) = 0    ! k0
             idpme(mv+3) = 230
             ime(1,mv+4) =-1    ! k-
             idpme(mv+4) =-130
          endif
       elseif((iin(1)+ipn(1)) == -1)  then
!  D-  + D0
          ime(1,mv+1) = 0    ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0    ! n
          idpme(mv+2) = 1220
          ime(1,mv+3) = 0    ! k0
          idpme(mv+3) = 230
          ime(1,mv+4) =-1    ! k-
          idpme(mv+4) =-130
       else
          write( *,*) ' bbkak2 d+n : ida,idb=',ida,idb
          write(16,*) ' bbkak2 d+n : ida,idb=',ida,idb
          return
       endif
    else
       write( *,*) ' bbkak3 d+n : ida,idb=',ida,idb
       write(16,*) ' bbkak3 d+n : ida,idb=',ida,idb
       return
    endif
    np=4
    amas4(1) = mn
    amas4(2) = mn
    amas4(3) = mk
    amas4(4) = mk
    call  genbodl(np,u,amas4,ps4,w4)
    do  k=1,np
       pme(4,mv+k) = ps4(1,k)
       pme(5,mv+k) = ps4(2,k)
       pme(6,mv+k) = ps4(3,k)
       pme(7,mv+k) = zro
       pme(9,mv+k) = amas4(k)
       ime(2,mv+k) = 0
       if(k <= 2) ime(3,mv+k) = 0  !  p or n
       if(k == 3) ime(3,mv+k) = 1  !  k+ or k0
       if(k == 4) ime(3,mv+k) =-1  !  k- or ak0
       if(k <= 2) ime(4,mv+k) = 1  !  p or n
       if(k > 2) ime(4,mv+k) = 0  !  k+ or k0,k- or ak0
       ime(5,mv+k) = 0
    enddo
    ibbkak = 1
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  pibkak(ida,idb,v,u,pin,pn,iin,ipn,mv,np,ipibkak)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    real(real64) ::  v(3),u,pin(9),pn(9),mk,mn
    logical yesela
    common/yesela/ yesela
    common/memorylaq/pme(9,5999),ime(5,5999)
    common /idpme/ idpme(5999)
    common/ncasca/ncas,ncpri
    dimension iin(5),ipn(5),amas3(3),ps3(5,3)
    data  mk/0.494d0/,mn/0.939d0/,zro/0.0d0/, &
         & onth/0.33333333d0/,twth/0.66666667d0/, &
         & onfo/0.25d0/,ontw/0.5d0/,thfo/0.75d0/, &
         & onfi/0.20d0/,thfi/0.6d0/,fofi/0.80d0/
    ipibkak = 0
    umin=mn + 2.0*mk
    if((ida == 120.or.ida == -120.or.ida == 110).and. &
         & (idb == 1120.or.idb == 1220.or. &
         & idb == 1111.or.idb == 1121.or.idb == 2221.or. &
         & idb == 1221).and.u > umin) then
       go  to  1
    else
       return
    endif
1   continue
    ru = (umin/u)**2
    sig0 = 1.121*((1.0d0 - ru)**1.86)*ru**2
    if(ida == -120.and.idb == 1120)      then
       skk=5.0*sig0                         ! pi-  +  p
    elseif(ida == -120.and.idb == 1220)  then
       skk=1.0*sig0                         ! pi-  +  n
    elseif(ida ==  120.and.idb == 1120)  then
       skk=1.0*sig0                         ! pi+  +  p
    elseif(ida ==  120.and.idb == 1220)  then
       skk=5.0*sig0                         ! pi+  +  n
    elseif(ida ==  110.and.idb == 1120)  then
       skk=3.0*sig0                         ! pi0  +  p
    elseif(ida ==  110.and.idb == 1220)  then
       skk=3.0*sig0                         ! pi0  +  n
    elseif(ida ==  120.and.idb == 1121)  then
       skk=2.0*sig0                         ! pi+  +  d+
    elseif(ida ==  110.and.idb == 1111)  then
       skk=3.0*sig0                         ! pi0  +  d++
    elseif(ida ==  120.and.idb == 2221)  then
       skk=4.0*sig0                         ! pi+  +  d0
    elseif(ida ==  110.and.idb == 1121)  then
       skk=5.0*sig0                         ! pi0  +  d+
    elseif(ida == -120.and.idb == 1111)  then
       skk=6.0*sig0                         ! pi-  +  d++
    elseif(ida ==  120.and.idb == 1221)  then
       skk=6.0*sig0                         ! pi+  +  d-
    elseif(ida ==  110.and.idb == 2221)  then
       skk=5.0*sig0                         ! pi0  +  d0
    elseif(ida == -120.and.idb == 1121)  then
       skk=4.0*sig0                         ! pi-  +  d+
    elseif(ida ==  110.and.idb == 1221)  then
       skk=3.0*sig0                         ! pi0  +  d-
    elseif(ida == -120.and.idb == 2221)  then
       skk=2.0*sig0                         ! pi-  +  d0
    else
       return
    endif
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sino=sito
    if(yesela)         go to 2
    call  crosec(0,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siel,0)
    call  crosec(2,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siex,0)
    sino=sito-siel-siex
2   if(sino <= 0.)  return
    if(rndm(-1.0_real64) > (skk/sino))  return
    if(ida == -120.and.idb == 1120)  then
!  pi-   +  p
       rnd=rndm(-1.0_real64)
       if(rnd <= onfi)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       elseif(rnd <= thfi)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida ==  120.and.idb == 1120)  then
!  pi+   +  p
       ime(1,mv+1) = 1              ! p
       idpme(mv+1) = 1120
       ime(1,mv+2) = 1              ! k+
       idpme(mv+2) = 130
       ime(1,mv+3) = 0              ! ak0
       idpme(mv+3) =-230
    elseif(ida == -120.and.idb == 1220)  then
!  pi-   +  n
       ime(1,mv+1) = 0              ! n
       idpme(mv+1) = 1220
       ime(1,mv+2) = 0              ! k0
       idpme(mv+2) = 230
       ime(1,mv+3) =-1              ! k-
       idpme(mv+3) =-130
    elseif(ida ==  120.and.idb == 1220)  then
!  pi+   +  n
       rnd=rndm(-1.0_real64)
       if(rnd <= onfi)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       elseif(rnd <= thfi)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida ==  110.and.idb == 1120)  then
!  pi0   +  p
       rnd=rndm(-1.0_real64)
       if(rnd <= twth)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       elseif(rnd <= 5./6.)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida ==  110.and.idb == 1220)  then
!  pi0   +  n
       rnd=rndm(-1.0_real64)
       if(rnd <= twth)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       elseif(rnd <= 5./6.)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida ==  120.and.idb == 1121)  then
!  pi+  +  D+
       ime(1,mv+1) = 1              ! p
       idpme(mv+1) = 1120
       ime(1,mv+2) = 1              ! k+
       idpme(mv+2) = 130
       ime(1,mv+3) = 0              ! ak0
       idpme(mv+3) =-230
    elseif(ida ==  110.and.idb == 1111)  then
!  pi0  +  D++
       ime(1,mv+1) = 1              ! p
       idpme(mv+1) = 1120
       ime(1,mv+2) = 1              ! k+
       idpme(mv+2) = 130
       ime(1,mv+3) = 0              ! ak0
       idpme(mv+3) =-230
    elseif(ida ==  120.and.idb == 2221)  then
!  pi+  +  D0
       rnd=rndm(-1.0_real64)
       if(rnd <= onfo)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       elseif(rnd <= ontw)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida ==  110.and.idb == 1121)  then
!  pi0  +  D+
       rnd=rndm(-1.0_real64)
       if(rnd <= onfi)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       elseif(rnd <= fofi)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida == -120.and.idb == 1111)  then
!  pi-  +  D++
       rnd=rndm(-1.0_real64)
       if(rnd <= ontw)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       else
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       endif
    elseif(ida ==  120.and.idb == 1221)  then
!  pi+  +  D-
       rnd=rndm(-1.0_real64)
       if(rnd <= ontw)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       else
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       endif
    elseif(ida ==  110.and.idb == 2221)  then
!  pi0  +  D0
       rnd=rndm(-1.0_real64)
       if(rnd <= onfi)         then
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       elseif(rnd <= fofi)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       else
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       endif
    elseif(ida == -120.and.idb == 1121)  then
!  pi-  +  D+
       rnd=rndm(-1.0_real64)
       if(rnd <= onfo)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 1              ! k+
          idpme(mv+2) = 130
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       elseif(rnd <= ontw)         then
          ime(1,mv+1) = 0              ! n
          idpme(mv+1) = 1220
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) = 0              ! ak0
          idpme(mv+3) =-230
       else
          ime(1,mv+1) = 1              ! p
          idpme(mv+1) = 1120
          ime(1,mv+2) = 0              ! k0
          idpme(mv+2) = 230
          ime(1,mv+3) =-1              ! k-
          idpme(mv+3) =-130
       endif
    elseif(ida ==  110.and.idb == 1221)  then
!  pi0  +  D-
       ime(1,mv+1) = 0              ! n
       idpme(mv+1) = 1220
       ime(1,mv+2) = 0              ! k0
       idpme(mv+2) = 230
       ime(1,mv+3) =-1              ! k-
       idpme(mv+3) =-130
    elseif(ida == -120.and.idb == 2221)  then
!  pi-  +  D0
       ime(1,mv+1) = 0              ! n
       idpme(mv+1) = 1220
       ime(1,mv+2) = 0              ! k0
       idpme(mv+2) = 230
       ime(1,mv+3) =-1              ! k-
       idpme(mv+3) =-130
    else
       write( *,*) ' pibkak3 pi+b : ida,idb=',ida,idb
       return
    endif
    np=3
    amas3(1) = mn
    amas3(2) = mk
    amas3(3) = mk
    call  genbodl(np,u,amas3,ps3,w3)
    do  k=1,np
       pme(4,mv+k) = ps3(1,k)
       pme(5,mv+k) = ps3(2,k)
       pme(6,mv+k) = ps3(3,k)
       pme(7,mv+k) = zro
       pme(9,mv+k) = amas3(k)
       ime(2,mv+k) = 0
       if(k == 1) ime(3,mv+k) = 0  !  p or n
       if(k == 2) ime(3,mv+k) = 1  !  k+ or k0
       if(k == 3) ime(3,mv+k) =-1  !  k- or ak0
       if(k == 1) ime(4,mv+k) = 1  !  p or n
       if(k > 1) ime(4,mv+k) = 0  !  k+ or k0,k- or ak0
       ime(5,mv+k) = 0
    enddo
    ipibkak = 1
!
    ie0=iin(1)+ipn(1)
    is0=iin(3)+ipn(3)
    ib0=iin(4)+ipn(4)
    ies=0
    iss=0
    ibs=0
    do  k=1,np
       ies=ies+ime(1,mv+k)
       iss=iss+ime(3,mv+k)
       ibs=ibs+ime(4,mv+k)
    enddo
!     if(ies.ne.ie0.or.iss.ne.is0.or.ibs.ne.ib0)  then
!       write( *,100)  ies,iss,ibs,ie0,is0,ib0
! 100   format(6I5)
!     endif
    return
  end
  subroutine  pneta(ida,idb,v,u,tin1,pin,pn,mv,np,ieta)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    real(real64) ::  v(3),u,tin1,pin(9),pn(9),mn,meta,m1,m2,m3
    logical yesela
    common/yesela/ yesela
    common/memorylaq/pme(9,5999),ime(5,5999)
    common/ncasca/ncas,ncpri
    common/psigeta/ seta
    dimension  p0(3),p0s(3),p1c(3),p2c(3),p3c(3), &
         & p1s(3),p2s(3),p3s(3),psum(3)
    data  mn/0.939/,meta/0.549/
    ieta=0
    seta = 0.0d0
    if(ida == 1120.and.idb == 1120)  then
       ksi=1
       ic1=1                 ! p
       ic2=1                 ! p
    elseif(ida == 1220.and.idb == 1220)  then
       ksi=1
       ic1=0                 ! n
       ic2=0                 ! n
    elseif(ida == 1120.and.idb == 1220)  then
       ksi=2
       ic1=1                 ! p
       ic2=0                 ! n
    elseif(ida == 1220.and.idb == 1120)  then
       ksi=2
       ic1=0                 ! n
       ic2=1                 ! p
    else
       return
    endif
    u0=meta+mn+mn
    if(u <= u0)                                return
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sino=sito
    if(yesela)         go to 1
    call  crosec(0,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siel,0)
    call  crosec(2,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siex,0)
    sino=sito-siel-siex
1   if(sino <= 0.)                     return
    seta=sppeta(tin1,ksi)
    reta=seta/sino
    if(rndm(-1.0_real64) > reta)         return
!
!     write(*,*) ' TIN1,SINO,seta,reta=',TIN1,SINO,seta,reta
!
    id1=ida
    id2=idb
    id3=220
    m1=mn
    m2=mn
    m3=meta
    ic3=0                  ! eta
!       SAMPLING ACCORDING TO 3-BODY PHASE SPACE VOLUME
    if(u <= (m1+m2+m3))  return
    em1=(u**2+m1**2-(m2+m3)**2)/(2.*u)
    em2=(u**2+m2**2-(m1+m3)**2)/(2.*u)
2   continue
    e1=m1+(em1-m1)*rndm(-1.0_real64)
    e2=m2+(em2-m2)*rndm(-1.0_real64)
    e3=u-e1-e2
    if(e3 <= m3)                                  go  to  2
    fnorm=27.*e1*e2*e3/u**3
    if(rndm(-1.0_real64) > fnorm)                        go  to  2
    p1=dsqrt(e1**2-m1**2)
    p2=dsqrt(e2**2-m2**2)
    p3=dsqrt(e3**2-m3**2)
    if(((p1+p2-p3)*(p1-p2+p3)*(p2+p3-p1)) <= 0.)  go  to  2
    ct3=1.-2.*rndm(-1.0_real64)
    fi3=6.283185*rndm(-1.0_real64)
    p0(1)=px1
    p0(2)=py1
    p0(3)=pz1
    call  kinemq(p0,v,p0s,ct0,st0,cf0,sf0,t0,am1)
    st3=dsqrt(1.-ct3**2)
    p3c(1)=p3*st3*dcos(fi3)
    p3c(2)=p3*st3*dsin(fi3)
    p3c(3)=p3*ct3
    call rotorq(p0s,v,p3c,p3s)
    ct1=-(p3**2+p1**2-p2**2)/(2.*p3*p1)
    ct2=-(p3**2+p2**2-p1**2)/(2.*p3*p2)
    st1=dsqrt(1.-ct1**2)
    st2=dsqrt(1.-ct2**2)
    fi1=6.283185*rndm(-1.0_real64)
    fi2=3.141592+fi1
    p1c(1)=p1*st1*dcos(fi1)
    p1c(2)=p1*st1*dsin(fi1)
    p1c(3)=p1*ct1
    p2c(1)=p2*st2*dcos(fi2)
    p2c(2)=p2*st2*dsin(fi2)
    p2c(3)=p2*ct2
    call rotorq(p3s,v,p1c,p1s)
    call rotorq(p3s,v,p2c,p2s)
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=p1s(1)
    pme(5,mv+1)=p1s(2)
    pme(6,mv+1)=p1s(3)
    pme(7,mv+1)=0.
    pme(8,mv+1)=p1
    pme(9,mv+1)=m1
    pme(1,mv+2)=0.
    pme(2,mv+2)=0.
    pme(3,mv+2)=0.
    pme(4,mv+2)=p3s(1)
    pme(5,mv+2)=p3s(2)
    pme(6,mv+2)=p3s(3)
    pme(7,mv+2)=0.
    pme(8,mv+2)=p3
    pme(9,mv+2)=m3
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p2s(1)
    pme(5,mv+3)=p2s(2)
    pme(6,mv+3)=p2s(3)
    pme(7,mv+3)=0.
    pme(8,mv+3)=p2
    pme(9,mv+3)=m2
    ime(1,mv+1)=ic1
    ime(2,mv+1)=0
    ime(3,mv+1)=0
    ime(4,mv+1)=1
    ime(5,mv+1)=0
    ime(1,mv+2)=ic3
    ime(2,mv+2)=0
    ime(3,mv+2)=0
    ime(4,mv+2)=0
    ime(5,mv+2)=idint(1000.0d0 *taun(8))
    if(ime(5,mv+2) < 1) ime(5,mv+2)=1
    ime(1,mv+3)=ic2
    ime(2,mv+3)=0
    ime(3,mv+3)=0
    ime(4,mv+3)=1
    ime(5,mv+3)=0
    ieta=1
    np = 3
    do  k=1,3
       psum(k)=p1s(k)+p2s(k)+p3s(k)
    enddo
!     IF(DABS(PSUM(1)).GT.1.D-3.OR.DABS(PSUM(2)).GT.1.D-3.OR.
!    &   DABS(PSUM(3)).GT.1.D-3) THEN
!    WRITE(16,*) 'PNETA: PSUM=',PSUM
!    WRITE( *,*) 'PNETA: PSUM=',PSUM
!     ENDIF
    if(ncas >= ncpri)  then
       pm1=dsqrt(p1s(1)**2+p1s(2)+p1s(3)**2)
       pm2=dsqrt(p2s(1)**2+p2s(2)+p2s(3)**2)
       pm3=dsqrt(p3s(1)**2+p3s(2)+p3s(3)**2)
       en1=dsqrt(pm1**2+m1**2)
       en2=dsqrt(pm2**2+m2**2)
       en3=dsqrt(pm3**2+m3**2)
       write(16,100) u,v
100    format(28x,'pneta: u,v=',1pe11.4,3(1pe11.4))
       write(16,101) id1,p1c,p1s,e1,en1,p1,pm1,m1,ic1
       write(16,101) id2,p2c,p2s,e2,en2,p2,pm2,m2,ic2
       write(16,101) id3,p3c,p3s,e3,en3,p3,pm3,m3,ic3
101    format(1x,i5,11(1pe11.4),i3)
    endif
!
!     write(*,*)  ' PNETA: TIN1,IC1,IC2,IC3',TIN1,IC1,IC2,IC3
!
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  piysex(id1,id2,v,u,pin,pn,mv,np,isex)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channel pi + Y ==> AKA + N
!     Calls: CROSEC, ABELQ
!
    real(real64) ::  mpi,mn,mk,ml,ms
!
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension v(3),pin(9),pn(9),p1(3),p2(3)
    data  mpi/0.140/,mn/0.940/,mk/0.494/,ml/1.116/,ms/1.189/
    isex=0
    ir=0
    if((id1 == 110.and.id2 == 2130).or.(id1 == 2130.and.id2 ==  110)) &
         & ir=1  ! pi0 + l
    if((id1 ==  120.and.id2 == 2130).or.(id1 == 2130.and.id2 ==  120)) &
         & ir=2  ! pi+ + l
    if((id1 == -120.and.id2 == 2130).or.(id1 == 2130.and.id2 == -120)) &
         & ir=3  ! pi- + l
    if((id1 ==  110.and.id2 == 1230).or.(id1 == 1230.and.id2 ==  110)) &
         & ir=4  ! pi0 + s0
    if((id1 ==  120.and.id2 == 1230).or.(id1 == 1230.and.id2 ==  120)) &
         & ir=5  ! pi+ + s0
    if((id1 == -120.and.id2 == 1230).or.(id1 == 1230.and.id2 == -120)) &
         & ir=6  ! pi- + s0
    if((id1 ==  110.and.id2 == 1130).or.(id1 == 1130.and.id2 ==  110)) &
         & ir=7  ! pi0 + s+
    if((id1 == -120.and.id2 == 1130).or.(id1 == 1130.and.id2 == -120)) &
         & ir=8  ! pi- + s+
    if((id1 ==  110.and.id2 == 2230).or.(id1 == 2230.and.id2 ==  110)) &
         & ir=9  ! pi0 + s-
    if((id1 ==  120.and.id2 == 2230).or.(id1 == 2230.and.id2 ==  120)) &
         & ir=10 ! pi+ + s-
    if(ir == 0)        return
!
    if(u <= (mn+mk))   return
    ek0=(u**2-mn**2-mk**2)/(2.*mn)
    pk0=sqrt(ek0**2-mk**2)
    eks=(u**2+mk**2-mn**2)/(2.*u)
    pks=sqrt(eks**2-mk**2)
    if(ir <= 3)  then
!   Calculation of K- + p=>L + pi0 cross section
!   param. of G.Q.Li et al NP A625(1997)342
       if(pk0 <= 0.6)        then
          skmp=1.205*pk0**(-1.428)
       elseif(pk0 <= 1.0)    then
          skmp=3.5*pk0**0.659
       else
          skmp=3.5*pk0**(-3.97)
       endif
       epis=(u**2+mpi**2-ml**2)/(2.*u)
       pis=sqrt(epis**2-mpi**2)
!   inverse cross section
       spiy=2.*((pks/pis)**2)*skmp
    else
!   Calculation of K- + p=>S0 + pi0 cross section
!   param. of G.Q.Li et al NP A625(1997)342
       if(pk0 <= 0.345)      then
          skmp=0.624*pk0**(-1.830)
       elseif(pk0 <= 0.425)  then
          skmp=0.0138/((pk0-0.385)**2+0.0017)
       else
          skmp=0.7*pk0**(-2.09)
       endif
       epis=(u**2+mpi**2-ms**2)/(2.*u)
       pis=sqrt(epis**2-mpi**2)
!   inverse cross section
       spiy=4./3.*((pks/pis)**2)*skmp
    endif
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    if(rndm(-1.0_real64) > (spiy/(sito+spiy)))  return
    go  to  (1,2,3,4,5,6,7,8,9,10),ir
1   if(rndm(-1.0_real64) <= 0.5)  then     ! pi0 + l  ==> k-   + p
       iek=-1
       ien= 1
    else                           ! pi0 + l  ==> ak0  + n
       iek= 0
       ien= 0
    endif
    go  to  11
2   continue                       ! pi+ + l  ==> ak0  + p
    iek= 0
    ien= 1
    go  to  11
3   continue                       ! pi- + l  ==> k-   + n
    iek=-1
    ien= 0
    go  to  11
4   if(rndm(-1.0_real64) <= 0.5)  then     ! pi0 + s0 ==> k-   + p
       iek=-1
       ien= 1
    else                           ! pi0 + s0 ==>ak0   + n
       iek= 0
       ien= 0
    endif
    go  to  11
5   continue                       ! pi+ + s0 ==>ak0   + p
    iek= 0
    ien= 1
    go  to  11
6   continue                       ! pi- + s0 ==> k-   + n
    iek=-1
    ien= 0
    go  to  11
7   continue                       ! pi0 + s+ ==>ak0   + p
    iek= 0
    ien= 1
    go  to  11
8   if(rndm(-1.0_real64) <= 0.5)  then     ! pi- + s+ ==> k-   + p
       iek=-1
       ien= 1
    else                           ! pi- + s+ ==>ak0   + n
       iek= 0
       ien= 0
    endif
    go  to  11
9   continue                       ! pi0 + s- ==> k-   + n
    iek=-1
    ien= 0
    go  to  11
10  if(rndm(-1.0_real64) <= 0.5)  then     ! pi+ + s- ==> k-   + p
       iek=-1
       ien= 1
    else                           ! pi+ + s- ==>ak0   + n
       iek= 0
       ien= 0
    endif
11  continue
    ct=2.*rndm(-1.0_real64)-1.
    fi=6.283185*rndm(-1.0_real64)
    call  abelq(pin,v,u,p1,p2,ct,fi,mk,mn)
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p1(1)
    pme(5,mv+3)=p1(2)
    pme(6,mv+3)=p1(3)
    pme(7,mv+3)=0.
    pme(9,mv+3)=mk
    ime(1,mv+3)=iek
    ime(2,mv+3)=0
    ime(3,mv+3)=-1
    ime(4,mv+3)=0
    ime(5,mv+3)=0
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=-p1(1)
    pme(5,mv+1)=-p1(2)
    pme(6,mv+1)=-p1(3)
    pme(7,mv+1)=0.
    pme(9,mv+1)=mn
    ime(1,mv+1)=ien
    ime(2,mv+1)=0
    ime(3,mv+1)=0
    ime(4,mv+1)=1
    ime(5,mv+1)=0
    isex=1
    np = 2
! write(*,*)  'piysex: IR=',IR
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  pinsex(id1,id2,v,u,pin,iin,pn,ipn,mv,np,isex)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of strange production by pion on nucleon
!     Calls: CROSEC, ABELQ
!
    real(real64) ::  mpi,mn
    logical yesela
    common/yesela/ yesela
    common/strexc/strexc
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension v(3),pin(9),pn(9),p1(3),p2(3)
    dimension  iin(5),ipn(5)
    data  mpi/0.140/,mn/0.940/
    isex=0
    strexc=0.
    ie1=iin(1)
    is1=iin(3)
    ib1=iin(4)
    ie2=ipn(1)
    is2=ipn(3)
    ib2=ipn(4)
    if((ib1+ib2).ne.1.or.is1.ne.0.or.is2.ne.0) return
    if((ie1+ie2) > 2.or.(ie1+ie2) < -1)      return
!      IF(PIN(9) > 1.d0.or.PN(9) > 1.d0)        RETURN
!      IF(IB1.EQ.0.AND.IS1.EQ.0.AND.(ID2.EQ.1120.OR.ID2.EQ.1220))
!     *   GO  TO  10
!      RETURN
!   10 CONTINUE
!      IF(ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110) THEN
!        TPIN=(U**2-(PIN(9)+PN(9))**2)/(2.*PN(9))
!      ELSE
    tpin=(u**2-(mpi+mn)**2)/2./mn
!      ENDIF
    if(tpin <= 0.759)           return
    ssp=sst(0,1,ie1,ie2,4,tpin,0)
    if((ie1+ie2) == 2.and.(ie1 == 2.or.ie2 == 2)) &
         & ssp=sst(0,1,  1,  1,4,tpin,0)
    ssm=sst(0,1,ie1,ie2,5,tpin,0)
    if((ie1+ie2) == -1.and.(ie1 == -1.or.ie2 == -1)) &
         & ssm=sst(0,1, -1,  0,5,tpin,0)
    ss0=sst(0,1,ie1,ie2,6,tpin,0)
    sl =sst(0,1,ie1,ie2,3,tpin,0)
    ssex=ssp+ssm+ss0+sl
    if(ssex <= 1.d-10)            return
    strexc=ssex
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sino=sito
    if(yesela)         go to 1
    call  crosec(0,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,siel,0)
    call  crosec(2,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,siex,0)
    sino=sito-siel-siex
1   if(sino <= 0.)                return
!
!      WRITE(16,100) TPIN,SL,SS0,SSM,SSP,SSEX,SINO
!  100 FORMAT(2X,'TPIN,SL,SS0,SSM,SSP,SSEX,SINO='/1X,7(1PE10.3,1X))
!
    rns=rndm(-1.0_real64)
!      WRITE(16,'(5x,''U,ID1,ID2,RNS='',5x,F8.3,1x,2I6,1x,F8.3)')
!     &  U,ID1,ID2,RNS
    if(rns > (ssex/sino))  return
    rn1=rndm(-1.0_real64)
    if(rn1 <= (sl/ssex))    then
!     K L
       amy=1.116
       iey= 0
    else if(rn1 <= ((sl+ss0)/ssex)) then
!     K S0
       amy=1.189
       iey= 0
    else if(rn1 <= ((sl+ss0+ssm)/ssex)) then
!     K S-
       amy=1.189
       iey=-1
    else
!     K S+
       amy=1.189
       iey=1
    endif
    amk=0.494
    if(u <= (amk+amy+0.001))       return
    if(pn(9) > 1..or.pin(9) > 1.)   amk=0.896
    if(u <= (amk+amy+0.001))       amk=0.494
    iek=ie1+ie2-iey
    if(iek < 0)  then
       write(16,*) '  pinsex: k- ???'
       write( *,*) '  pinsex: k- ???'
    endif

!      WRITE(16,'(1x,''U,AMK,AMY,ID1,ID2='',3(F8.3,1x),2I6))')
!    & U,AMK,AMY,ID1,ID2
    tmax=3.5344*tpin*(tpin+mpi)/(1.88*(tpin+mpi)+mpi*mpi+0.8836)
    bet=bhn(tpin,mpi)
    exm=bet*tmax
    if(exm > 30.)    exm=30.
    ct=1.+(2.*log(1.+rndm(-1.0_real64)*(exp(-exm)-1.)))/exm
    fi=6.283185*rndm(-1.0_real64)
    call  abelq(pin,v,u,p1,p2,ct,fi,amk,amy)
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p1(1)
    pme(5,mv+3)=p1(2)
    pme(6,mv+3)=p1(3)
    pme(7,mv+3)=0.
    pme(9,mv+3)=amk
    ime(1,mv+3)=iek
    ime(2,mv+3)=0
    ime(3,mv+3)=1
    ime(4,mv+3)=0
    ime(5,mv+3)=0
    if(amk > 0.50)  ime(5,mv+3)=intg(1000. *taun(12))
!      if(AMK > 0.50)  IME(5,MV+3)=INTG(1000. *TAUN(14)) !for K*0
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=-p1(1)
    pme(5,mv+1)=-p1(2)
    pme(6,mv+1)=-p1(3)
    pme(7,mv+1)=0.
    pme(9,mv+1)=amy
    ime(1,mv+1)=iey
    ime(2,mv+1)=0
    ime(3,mv+1)=-1
    ime(4,mv+1)=1
    ime(5,mv+1)=0
    isex=1
    np = 2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  pineta(id1,id2,v,u,pin,iin,pn,ipn,mv,np,ieta)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channel pi + N ==> eta + N
!     Calls: CROSEC, ABELQ
!
    real(real64) ::  mpi,mn,meta
    logical yesela
!
    common/yesela/ yesela
    common/spinec/spinec
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension v(3),pin(9),pn(9),p1(3),p2(3)
    dimension  iin(5),ipn(5)
    data  mpi/0.140/,mn/0.940/,meta/0.549/
    ieta=0
    spinec=0.
    if((id1 == 120.or.id1 == -120.or.id1 == 110).and. &
         & (id2 == 1120.or.id2 == 1220)) go  to  10
    return
10  continue
    ie1=iin(1)
    ie2=ipn(1)
    if((ie1+ie2) == 2.or.(ie1+ie2) == -1)     return
    tpin=(u**2-(pin(9)+pn(9))**2)/(2.*pn(9))
    spinec=spinet(tpin)
    if(spinec <= 1.d-10)                      return
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sino=sito
    if(yesela)         go to 1
    call  crosec(0,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,siel,0)
    call  crosec(2,id1,id2,px1,py1,pz1,am1,px2,py2,pz2,am2,siex,0)
    sino=sito-siel-siex
1   if(sino <= 0.)                return
!
!      WRITE(16,100) TPIN,SPINEC,SINO
!  100 FORMAT(2X,'TPIN,SPINEC,SINO='/1X,3(1PE10.3,1X))
!
    rns=rndm(-1.0_real64)
    if(rns > (spinec/sino))      return
    tmax=3.5344*tpin*(tpin+mpi)/(1.88*(tpin+mpi)+mpi*mpi+0.8836)
    bet=bhn(tpin,mpi)
    exm=bet*tmax
    if(exm > 30.)    exm=30.
    ct=1.+(2.*log(1.+rndm(-1.0_real64)*(exp(-exm)-1.)))/exm
    fi=6.283185*rndm(-1.0_real64)
    call  abelq(pin,v,u,p1,p2,ct,fi,meta,mn)
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p1(1)
    pme(5,mv+3)=p1(2)
    pme(6,mv+3)=p1(3)
    pme(7,mv+3)=0.
    pme(9,mv+3)=meta
    ime(1,mv+3)=0
    ime(2,mv+3)=0
    ime(3,mv+3)=0
    ime(4,mv+3)=0
    ime(5,mv+3)=intg(1000. *taun(8))
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=-p1(1)
    pme(5,mv+1)=-p1(2)
    pme(6,mv+1)=-p1(3)
    pme(7,mv+1)=0.
    pme(9,mv+1)=mn
    ime(1,mv+1)=ie1+ie2
    ime(2,mv+1)=0
    ime(3,mv+1)=0
    ime(4,mv+1)=1
    ime(5,mv+1)=0
    ieta=1
    np = 2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  pipikk(id1,id2,v,u,pin,iin,pn,ipn,mv,np,ikak)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!      calculation of channel pi + pi ==> K + AKA
!     Calls: PIPICS, ABELQ
!
    real(real64) ::  mk
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension v(3),pin(9),pn(9),p1(3),p2(3)
    dimension  iin(5),ipn(5)
    data  mk/0.495/
    data a1/8.6335481/,a2/0.554395/,a3/2.750998/,a4/0.893608/
    ikak=0
    if((id1 ==  120.and.id2 == -120).or. &
         & (id1 ==  120.and.id2 ==  110).or. &
         & (id1 == -120.and.id2 ==  120).or. &
         & (id1 == -120.and.id2 ==  110).or. &
         & (id1 ==  110.and.id2 ==  120).or. &
         & (id1 ==  110.and.id2 ==  110).or. &
         & (id1 ==  110.and.id2 == -120))  go  to 10
    return
10  continue
    if(u <= (2.*mk))  return

!  kkg  24.03.06
! x=(U/(2.*MK))**2
! SKK=a1*((x-1.)**a2)/(a3**2+x)**a4
    stk=1.6    ! b.tomasik, e.kolomeitsev,nucl-th/0512088
    if((id1 ==  120.and.id2 == -120).or. &
         & (id1 == -120.and.id2 ==  120))   skk=stk
    if( id1 ==  110.and.id2 ==  110)    skk=2.0/5.0*stk
    if((id1 ==  120.and.id2 ==  110).or. &
         & (id1 == -120.and.id2 ==  110))   skk=6.0/5.0*stk

    stot=pipics(id1,id2,u,1)
    if(rndm(-1.0_real64) > (skk/(skk+stot)))  return
!  write(16,*) 'PIPIKK:SKK,STOT=',SKK,STOT
    ct=1.-2.*rndm(-1.0_real64)
    fi=6.283185*rndm(-1.0_real64)
    call  abelq(pin,v,u,p1,p2,ct,fi,mk,mk)
    ie1=iin(1)
    ie2=ipn(1)
    if((ie1+ie2) == 0)         then
       if(rndm(-1.0_real64) < 0.5)  then
          ime(1,mv+1)= 1
          ime(1,mv+2)=-1
       else
          ime(1,mv+1)= 0
          ime(1,mv+2)= 0
       endif
    elseif((ie1+ie2) == 1)  then
       ime(1,mv+1)= 1
       ime(1,mv+2)= 0
    elseif((ie1+ie2) == -1) then
       ime(1,mv+1)= 0
       ime(1,mv+2)=-1
    else
       write(16,*) 'pipikk:id1,id2,ie1,ie2=',id1,id2,ie1,ie2
       return
    endif
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=p1(1)
    pme(5,mv+1)=p1(2)
    pme(6,mv+1)=p1(3)
    pme(7,mv+1)=0.
    pme(9,mv+1)=mk
    ime(2,mv+1)=0
    ime(3,mv+1)=1
    ime(4,mv+1)=0
    ime(5,mv+1)=0
    pme(1,mv+2)=0.
    pme(2,mv+2)=0.
    pme(3,mv+2)=0.
    pme(4,mv+2)=-p1(1)
    pme(5,mv+2)=-p1(2)
    pme(6,mv+2)=-p1(3)
    pme(7,mv+2)=0.
    pme(9,mv+2)=mk
    ime(2,mv+2)=0
    ime(3,mv+2)=-1
    ime(4,mv+2)=0
    ime(5,mv+2)=0
    ikak=1
    np = 2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  bbsex(ida,idb,v,u,pin,pn,mv,np,isex)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels B + B ==> Y + B + K,
!     where B=N, Delta, Y=Lambda,SIGMA, K=K+,K0
!     Calls: CROSEC, BBBYK, ROTORQ
!
    real(real64) ::  mn,mk,ml,ms,md,m1,m2,m3
!
    logical yesela
    common/yesela/ yesela
    common/memorylaq/pme(9,5999),ime(5,5999)
    common/ncasca/ncas,ncpri
    dimension v(3),pin(9),pn(9)
    dimension  p0(3),p0s(3),p1c(3),p2c(3),p3c(3), &
         & p1s(3),p2s(3),p3s(3),psum(3)
    data  mn/0.939/,mk/0.494/,ml/1.115/,ms/1.192/,md/1.236/
    isex=0
    u0=mk+ml+mn
    if(u <= u0)                                return
    px1=pin(4)
    py1=pin(5)
    pz1=pin(6)
    am1=pin(9)
    px2=pn(4)
    py2=pn(5)
    pz2=pn(6)
    am2=pn(9)
    call  crosec(1,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,sito,0)
    sino=sito
    if(yesela)         go to 1
    call  crosec(0,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siel,0)
    call  crosec(2,ida,idb,px1,py1,pz1,am1,px2,py2,pz2,am2,siex,0)
    sino=sito-siel-siex
1   if(sino <= 0.)                     return
    call bbbyk(ida,idb,id1,id2,id3,u,sigsex)
    if(sigsex < 1.d-10)               return
    if(rndm(-1.0_real64) > sigsex/sino)       return
    if(id1 == 1120.or.id1 == 1220)   then
       m1=mn
       if(id1 == 1120)   ic1= 1
       if(id1 == 1220)   ic1= 0
    elseif(id1 == 1111.or.id1 == 1121.or.id1 == 2221.or.id1 == 1221) &
            & then
       m1=md
       if(id1 == 1111)   ic1= 2
       if(id1 == 1121)   ic1= 1
       if(id1 == 2221)   ic1=-1
       if(id1 == 1221)   ic1= 0
       call  wmd(md,tau,fmd)
    else
       write(16,*) 'bbsex: id1=',id1
       write( *,*) 'bbsex: id1=',id1
    endif
    if(id2 == 2130)  then
       m2=ml
       ic2=0
    elseif(id2 == 1130.or.id2 == 2230.or.id2 == 1230)  then
       m2=ms
       if(id2 == 1130)   ic2= 1
       if(id2 == 2230)   ic2=-1
       if(id2 == 1230)   ic2= 0
    else
       write(16,*) 'bbsex: id2=',id2
       write( *,*) 'bbsex: id2=',id2
    endif
    if(id3 == 130.or.id3 == 230)  then
       m3=mk
       if(id3 == 230)   ic3= 0
       if(id3 == 130)   ic3= 1
    else
       write(16,*) 'bbsex: id3=',id3
       write( *,*) 'bbsex: id3=',id3
    endif
!       SAMPLING ACCORDING TO 3-BODY PHASE SPACE VOLUME
    if(u <= (m1+m2+m3))  return
    em1=(u**2+m1**2-(m2+m3)**2)/(2.*u)
    em2=(u**2+m2**2-(m1+m3)**2)/(2.*u)
2   continue
    e1=m1+(em1-m1)*rndm(-1.0_real64)
    e2=m2+(em2-m2)*rndm(-1.0_real64)
    e3=u-e1-e2
    if(e3 <= m3)                                  go  to  2
    fnorm=27.*e1*e2*e3/u**3
    if(rndm(-1.0_real64) > fnorm)                        go  to  2
    p1=sqrt(e1**2-m1**2)
    p2=sqrt(e2**2-m2**2)
    p3=sqrt(e3**2-m3**2)
    if(((p1+p2-p3)*(p1-p2+p3)*(p2+p3-p1)) <= 0.)  go  to  2
    ct3=1.-2.*rndm(-1.0_real64)
    fi3=6.283185*rndm(-1.0_real64)
    p0(1)=px1
    p0(2)=py1
    p0(3)=pz1
    call  kinemq(p0,v,p0s,ct0,st0,cf0,sf0,t0,am1)
    st3=sqrt(1.-ct3**2)
    p3c(1)=p3*st3*cos(fi3)
    p3c(2)=p3*st3*sin(fi3)
    p3c(3)=p3*ct3
    call rotorq (p0s,v,p3c,p3s)
    ct1=-(p3**2+p1**2-p2**2)/(2.*p3*p1)
    ct2=-(p3**2+p2**2-p1**2)/(2.*p3*p2)
    st1=sqrt(1.-ct1**2)
    st2=sqrt(1.-ct2**2)
    fi1=6.283185*rndm(-1.0_real64)
    fi2=3.141592+fi1
    p1c(1)=p1*st1*cos(fi1)
    p1c(2)=p1*st1*sin(fi1)
    p1c(3)=p1*ct1
    p2c(1)=p2*st2*cos(fi2)
    p2c(2)=p2*st2*sin(fi2)
    p2c(3)=p2*ct2
    call rotorq (p3s,v,p1c,p1s)
    call rotorq (p3s,v,p2c,p2s)
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=p1s(1)
    pme(5,mv+1)=p1s(2)
    pme(6,mv+1)=p1s(3)
    pme(7,mv+1)=0.
    pme(8,mv+1)=p1
    pme(9,mv+1)=m1
    pme(1,mv+2)=0.
    pme(2,mv+2)=0.
    pme(3,mv+2)=0.
    pme(4,mv+2)=p3s(1)
    pme(5,mv+2)=p3s(2)
    pme(6,mv+2)=p3s(3)
    pme(7,mv+2)=0.
    pme(8,mv+2)=p3
    pme(9,mv+2)=m3
    pme(1,mv+3)=0.
    pme(2,mv+3)=0.
    pme(3,mv+3)=0.
    pme(4,mv+3)=p2s(1)
    pme(5,mv+3)=p2s(2)
    pme(6,mv+3)=p2s(3)
    pme(7,mv+3)=0.
    pme(8,mv+3)=p2
    pme(9,mv+3)=m2
    ime(1,mv+1)=ic1
    ime(2,mv+1)=0
    ime(3,mv+1)=0
    ime(4,mv+1)=1
    if(id1 == 1120.or.id1 == 1220)   then
       ime(5,mv+1)=0
    else
       ime(5,mv+1)=intg(1000.*tau)
       if(ime(5,mv+1) < 1) ime(5,mv+1)=1
    endif
    ime(1,mv+2)=ic3
    ime(2,mv+2)=0
    ime(3,mv+2)=1
    ime(4,mv+2)=0
    ime(5,mv+2)=0
    ime(1,mv+3)=ic2
    ime(2,mv+3)=0
    ime(3,mv+3)=-1
    ime(4,mv+3)=1
    ime(5,mv+3)=0
    isex=1
    np = 3
    do  k=1,3
       psum(k)=p1s(k)+p2s(k)+p3s(k)
    enddo
    if(abs(psum(1)) > 1.d-6.or.abs(psum(2)) > 1.d-6.or. &
         & abs(psum(3)) > 1.d-6) then
!    WRITE(16,*) 'BBSEX: PSUM=',PSUM
!    WRITE( *,*) 'BBSEX: PSUM=',PSUM
    endif
    if(ncas >= ncpri)  then
       write(16,100) u,v
100    format(28x,'bbsex: u,v=',1pe11.4,3(1pe11.4))
       write(16,101) id1,p1c,p1s,e1,p1,m1,ic1
       write(16,101) id2,p2c,p2s,e2,p2,m2,ic2
       write(16,101) id3,p3c,p3s,e3,p3,m3,ic3
101    format(1x,i5,9(1pe11.4),i3)
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine bbbyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation probabilities of channels B + B ==> Y + B + K,
!     where B=N, Delta, Y=Lambda,SIGMA, K=K+,K0
!     Calls: DDNYK,DDDYK,DNNYK, DNDYK,NNNYK,NNDYK
!
    sig=0.
    if((ida == 1111.or.ida == 1121.or.ida == 2221.or.ida == 1221).and. &
         & (idb == 1111.or.idb == 1121.or.idb == 2221.or.idb == 1221)) &
         & then
!     Delta+Delta
       call  ddnyk(ida,idb,idn1,idn2,idn3,ss,sign)
       call  dddyk(ida,idb,idd1,idd2,idd3,ss,sigd)

!   kkg  04/12/07
       sig=sign+sigd
       if(sig < 1.d-10)     return
       if(rndm(-1.0_real64) <= sign/sig)  then
          id1=idn1
          id2=idn2
          id3=idn3
       else
          id1=idd1
          id2=idd2
          id3=idd3
       endif
    elseif((ida == 1120.or.ida == 1220).and. &
            & (idb == 1111.or.idb == 1121.or.idb == 2221.or.idb == 1221))  then
!  N+Delta
       call  dnnyk(idb,ida,idn1,idn2,idn3,ss,sign)
       call  dndyk(idb,ida,idd1,idd2,idd3,ss,sigd)
       sig=sign+sigd
       if(sig < 1.d-10)     return
       if(rndm(-1.0_real64) <= sign/sig)  then
          id1=idn1
          id2=idn2
          id3=idn3
       else
          id1=idd1
          id2=idd2
          id3=idd3
       endif
    elseif((idb == 1120.or.idb == 1220).and. &
            & (ida == 1111.or.ida == 1121.or.ida == 2221.or.ida == 1221)) then
!  Delta+N
       call  dnnyk(ida,idb,idn1,idn2,idn3,ss,sign)
       call  dndyk(ida,idb,idd1,idd2,idd3,ss,sigd)
       sig=sign+sigd
       if(sig < 1.d-10)     return
       if(rndm(-1.0_real64) <= sign/sig)  then
          id1=idn1
          id2=idn2
          id3=idn3
       else
          id1=idd1
          id2=idd2
          id3=idd3
       endif
    elseif((ida == 1120.or.ida == 1220).and. &
            & (idb == 1120.or.idb == 1220))         then
!  N+N
       call  nnnyk(idb,ida,idn1,idn2,idn3,ss,sign)
       call  nndyk(idb,ida,idd1,idd2,idd3,ss,sigd)
       sig=sign+sigd
       if(sig < 1.d-10)     return
       if(rndm(-1.0_real64) <= sign/sig)  then
          id1=idn1
          id2=idn2
          id3=idn3
       else
          id1=idd1
          id2=idd2
          id3=idd3
       endif
    else
       sig=0.
       write(16,*) 'bbbyk: ida,idb=',ida,idb
       write( *,*) 'bbbyk: ida,idb=',ida,idb
       return
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  ddnyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels Delta + Delta ==> Y + N + K,
!     K=K+,K0
!
    real(real64) ::  mn,mk,ms,ml
    data  mn/0.940/,mk/0.494/,ml/1.115/,ms/1.192/
! Parameters for isospin-averaged NN==>NLK cross section
    data al/2.330/,bl/2.140/,cl/5.024/
! Parameters for isospin-averaged NN==>NSK cross section
    data as/15.49/,bs/2.768/,cs/7.222/
    sig=0.
    if((ida == 1111.and.idb == 1111).or. &
         & (ida == 2221.and.idb == 2221))    return
    ss0l=mn+ml+mk
    if(ss <= ss0l)                  return
    ss0s=mn+ms+mk
    sig1=1./4.*sibir(al,bl,cl,ss,ss0l)
    sig2=1./4.*sibir(as,bs,cs,ss,ss0s)
    if((ida == 1111.and.idb == 1121).or. &
         & (ida == 1121.and.idb == 1111))    then
! D++ D+
       sig=1./4.*sig2
       if(sig < 1.d-10)               return
! p S+ K+
       id1=1120           ! p
       id2=1130           ! s+
       id3=130            ! k+
    elseif((ida == 1111.and.idb == 1221).or. &
            & (ida == 1221.and.idb == 1111).or. &
            & (ida == 1121.and.idb == 1121))    then
! D++ D0 or D+ D+
       sigl=1./4.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=3./4.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
! p L K+
          id1=1120           ! p
          id2=2130           ! l
          id3=130            ! k+
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 1./3.)  then
! p S+ K0
             id1=1120           ! p
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd < 2./3.)  then
! p S0 K+
             id1=1120           ! p
             id2=1230           ! s0
             id3=130            ! k+
          else
! n S+ K+
             id1=1220           ! n
             id2=1130           ! s+
             id3=130            ! k+
          endif
       endif
    elseif((ida == 1121.and.idb == 1221).or. &
            & (ida == 1221.and.idb == 1121).or. &
            & (ida == 1111.and.idb == 2221).or. &
            & (ida == 2221.and.idb == 1111))    then
! D+ D0 or D++ D-
       sigl=2./4.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=4./4.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
          if(rndm(-1.0_real64) < 0.5)  then
! p L K0
             id1=1120           ! p
             id2=2130           ! l
             id3=230            ! k0
          else
! n L K+
             id1=1220           ! n
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 0.25)  then
! p S0 K0
             id1=1120           ! p
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd < 0.50)  then
! p S+ K0
             id1=1120           ! p
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd < 0.75)  then
! n S0 K+
             id1=1220           ! n
             id2=1230           ! s0
             id3=130            ! k+
          else
! p S- K+
             id1=1220           ! n
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    elseif((ida == 1121.and.idb == 2221).or. &
            & (ida == 2221.and.idb == 1121).or. &
            & (ida == 1221.and.idb == 1221))    then
! D+ D- or D0 D0
       sigl=1./4.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=3./4.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
! n L K0
          id1=1220           ! n
          id2=2130           ! l
          id3=230            ! k0
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 1./3.)  then
! n S0 K0
             id1=1220           ! n
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd < 2./3.)  then
! n S- K+
             id1=1220           ! n
             id2=2230           ! s-
             id3=130            ! k+
          else
! p S- K0
             id1=1120           ! p
             id2=2230           ! s-
             id3=230            ! k0
          endif
       endif
    elseif((ida == 1221.and.idb == 2221).or. &
            & (ida == 2221.and.idb == 1221))    then
! D0 D-
       sig=1./4.*sig2
       if(sig < 1.d-10)               return
! n S- K0
       id1=1220           ! n
       id2=2230           ! s-
       id3=230            ! k0
    else
       write(16,*) 'ddnyk: ida,idb=',ida,idb
       write( *,*) 'ddnyk: ida,idb=',ida,idb
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  dddyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels Delta + Delta ==> Y + Delta + K,
!     K=K+,K0
!
    real(real64) ::  md,mk,ml,ms
    data  md/1.236/,mk/0.494/,ml/1.115/,ms/1.192/
! Parameters for isospin-averaged DD==>DLK cross section
    data al/0.337/,bl/2.149/,cl/7.967/
! Parameters for isospin-averaged DD==>DSK cross section
    data as/5.140/,bs/2.952/,cs/12.05/
    sig=0.
    ss0l=md+ml+mk
    if(ss <= ss0l)                  return
    ss0s=md+ms+mk
    sig1=1./2.*sibir(al,bl,cl,ss,ss0l)
    sig2=1./2.*sibir(as,bs,cs,ss,ss0s)
    if(ida == 1111.and.idb == 1111)        then
! D++ D++  ==>  DYK
       sig=sig2
       if(sig < 1.d-10)               return
! D++ S+ K+
       id1=1111           ! d++
       id2=1130           ! s+
       id3=130            ! k+
    elseif((ida == 1111.and.idb == 1121).or. &
            & (ida == 1121.and.idb == 1111))  then
! D++ D+  ==>  DYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=3.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
! D++ L K+
          id1=1111           ! d++
          id2=2130           ! l
          id3=130            ! k+
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 0.33333)  then
! D++ S+ K0
             id1=1111           ! d++
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd < 0.66667)  then
! D++ S0 K+
             id1=1111           ! d++
             id2=1230           ! s0
             id3=130            ! k+
          else
! D+ S+ K+
             id1=1121           ! d+
             id2=1130           ! s+
             id3=130            ! k+
          endif
       endif
    elseif((ida == 1111.and.idb == 1221).or. &
            & (ida == 1221.and.idb == 1111).or. &
            & (ida == 1121.and.idb == 1121))    then
! D++ D+ or D+ D+ ==>  DYK
       sigl=2.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=5.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
          if(rndm(-1.0_real64) < 0.5)  then
! D++ L K0
             id1=1111           ! d++
             id2=2130           ! l    kkg 01/19/09
             id3=230            ! k0
          else
! D+ L K+
             id1=1121             ! d+
             id2=2130             ! l   kkg 01/19/09
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 0.2)  then
! D++ S0 K0
             id1=1111           ! d++
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd < 0.4)  then
! D++ S- K+
             id1=1111           ! d++
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd < 0.6)  then
! D+ S0 K+
             id1=1121           ! d+
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd < 0.8)  then
! D+ S+ K0
             id1=1121           ! d+
             id2=1130           ! s+
             id3=230            ! k0
          else
! D0 S+ K+
             id1=1221           ! d0
             id2=1130           ! s+
             id3=130            ! k+
          endif
       endif
    elseif((ida == 1111.and.idb == 2221).or. &
            & (ida == 2221.and.idb == 1111).or. &
            & (ida == 1121.and.idb == 1221).or. &
            & (ida == 1221.and.idb == 1121))    then
! D++ D- or D+ D0
       sigl=2.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=6.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
          if(rndm(-1.0_real64) < 0.5)  then
! D+ L K0
             id1=1121           ! d+
             id2=2130           ! l
             id3=230            ! k0
          else
! D0 L K+
             id1=1221           ! d0
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 0.16667)  then
! D++ S- K0
             id1=1111           ! d++
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd < 0.33333)  then
! D+ S0 K0
             id1=1121           ! d+
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd < 0.5)  then
! D+ S- K+
             id1=1121           ! d+
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd < 0.66667)  then
! D0 S+ K0
             id1=1221           ! d0
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd < 0.83333)  then
! D0 S0 K+
             id1=1221           ! d0
             id2=1230           ! s0
             id3=130            ! k+
          else
! D- S+ K+
             id1=2221           ! d-
             id2=1130           ! s+
             id3=130            ! k+
          endif
       endif
    elseif((ida == 1121.and.idb == 2221).or. &
            & (ida == 2221.and.idb == 1121).or. &
            & (ida == 1221.and.idb == 1221))    then
! D+ D- or D0 D0 ==>  DYK
       sigl=2.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=5.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
          if(rndm(-1.0_real64) < 0.5)  then
! D0 L K0
             id1=1221           ! d0
             id2=2130           ! l
             id3=230            ! k0
          else
! D- L K+
             id1=2221           ! d-
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 0.2)  then
! D+ S- K0
             id1=1121           ! d+
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd < 0.4)  then
! D0 S0 K0
             id1=1221           ! d0
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd < 0.6)  then
! D0 S- K+
             id1=1221           ! d0
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd < 0.8)  then
! D- S+ K0
             id1=2221           ! d-
             id2=1130           ! s+
             id3=230            ! k0
          else
! D- S0 K+
             id1=2221           ! d-
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    elseif((ida == 2221.and.idb == 1221).or. &
            & (ida == 1221.and.idb == 2221))  then
! D- D0  ==>  DYK
       sigl=1.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=3.*sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)               return
       if(rndm(-1.0_real64) < sigl/sig)  then
! D- L K0
          id1=2221           ! d-
          id2=2130           ! l
          id3=230            ! k0
       else
          rnd=rndm(-1.0_real64)
          if(rnd < 0.33333)  then
! D- S0 K0
             id1=2221           ! d-
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd < 0.66667)  then
! D- S- K+
             id1=2221           ! d-
             id2=2230           ! s-
             id3=130            ! k+
          else
! D0 S- K0
             id1=1221           ! d0
             id2=2230           ! s-
             id3=230            ! k0
          endif
       endif
    elseif(ida == 2221.and.idb == 2221)        then
! D- D-  ==>  DYK
       sig=1.*sig2
       if(sig < 1.d-10)               return
! D- S- K0
       id1=2221           ! d-
       id2=2230           ! s-
       id3=230            ! k0
    else
       write(16,*) 'dddyk: ida,idb=',ida,idb
       write( *,*) 'dddyk: ida,idb=',ida,idb
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  dndyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels Delta + N ==> Y + Delta + K,
!     K=K+,K0
!
    real(real64) ::  md,mk,ml,ms
    data  md/1.236/,mk/0.494/,ml/1.115/,ms/1.192/
! Parameters for D++ p ==>D++ L  K+ cross section
    data  a1/2.704/, b1/2.303/,c1/5.551/
! Parameters for D++ p ==>D++ S0 K+ cross section
    data  a2/10.30/, b2/2.748/,c2/9.321/
! Parameters for D++ n ==>D++ S- K+ cross section
    data  a7/10.33/, b7/2.743/,c7/8.915/
! Parameters for D+  p ==>D+  L  K+ cross section
    data a15/2.917/,b15/2.350/,c15/6.557/
! Parameters for D+  p ==>D+  S0 K+ cross section
    data a17/10.62/,b17/2.759/,c17/10.20/
! Parameters for D+  p ==>D0  S+ K+ cross section
    data a18/0.647/,b18/2.830/,c18/3.862/
! Parameters for D0  p ==>D+  S- K+ cross section
    data a21/2.128/,b21/2.843/,c21/5.986/
! Parameters for D+  n ==>D+  S- K+ cross section
    data a29/10.57/,b29/2.757/,c29/10.11/
! Parameters for D+  n ==>D0  L  K+ cross section
    data a30/0.312/,b30/2.110/,c30/2.165/
! Parameters for D+  n ==>D0  S0 K+ cross section
    data a31/1.112/,b31/2.846/,c31/5.943/
    sig=0.
    ss0l=md+ml+mk
    if(ss <= ss0l)                  return
    ss0s=md+ms+mk
    sig1 =sibir(a1,b1,c1,ss,ss0l)
    sig2 =sibir(a2,b2,c2,ss,ss0s)
    sig7 =sibir(a7,b7,c7,ss,ss0s)
    sig11=0.
    sig15=sibir(a15,b15,c15,ss,ss0l)
    sig17=sibir(a17,b17,c17,ss,ss0s)
    sig18=sibir(a18,b18,c18,ss,ss0s)
    sig21=sibir(a21,b21,c21,ss,ss0s)
    sig29=sibir(a29,b29,c29,ss,ss0s)
    sig30=sibir(a30,b30,c30,ss,ss0l)
    sig31=sibir(a31,b31,c31,ss,ss0s)
    if(ida == 1111.and.idb == 1120) then
! D++ p ==>DYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig2+sig7+3./4.*sig18
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  D++ L K+
          id1=1111           ! d++
          id2=2130           ! l
          id3=130            ! k+
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig2/sigs)  then
!  D++ S0 K+
             id1=1111           ! d++
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= (sig2+sig7)/sigs)  then
!  D++ S+ K0
             id1=1111           ! d++
             id2=1130           ! s+
             id3=230            ! k0
          else
!  D+ S+ K+
             id1=1121           ! d+
             id2=1130           ! s+
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1111.and.idb == 1220) then
! D++ n ==>DYK
       sigl=sig1+3./4.*sig30
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig2+sig7+3./4.*sig21+3./4.*sig31+sig11
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= sig1/sigl)  then
!  D++ L K0
             id1=1111           ! d++
             id2=2130           ! l
             id3=230            ! k0
          else
!  D+ L K+
             id1=1121           ! d+
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig2/sigs)  then
!  D++ S0 K0
             id1=1111           ! d++
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (sig2+sig7)/sigs)  then
!  D++ S- K+
             id1=1111           ! d++
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= (sig2+sig7+3./4.*sig21)/sigs)  then
!  D+ S+ K0
             id1=1121           ! d+
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= (sig2+sig7+3./4.*sig21+sig11)/sigs)  then
!  D0 S+ K+
             id1=1221           ! d0
             id2=1130           ! s+
             id3=130            ! k+
          else
!  D+ S0 K+
             id1=1121           ! d+
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1121.and.idb == 1120) then
! D+ p ==>DYK
       sigl=3./4.*sig30+sig15
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=3./4.*sig31+3./4.*sig21+sig15+sig21+sig17+sig18
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 3./4.*sig30/sigl)  then
!  D++ L K0
             id1=1111           ! d++
             id2=2130           ! l
             id3=230            ! k0
          else
!  D+ L K+
             id1=1121           ! d+
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 3./4.*sig31/sigs)  then
!  D++ S0 K0
             id1=1111           ! d++
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (3./4.*(sig31+sig21))/sigs)  then
!  D++ S- K+
             id1=1111           ! d++
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= (3./4.*(sig31+sig21)+sig21)/sigs)  then
!  D+ S+ K0
             id1=1121           ! d+
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= (3./4.*(sig31+sig21)+sig21+sig18)/sigs) then
!  D0 S+ K+
             id1=1221           ! d0
             id2=1130           ! s+
             id3=130            ! k+
          else
!  D+ S0 K+
             id1=1121           ! d+
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1221.and.idb == 1120) then
! D0 p ==>DYK
       sigl=sig30+sig15
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig31+sig21+sig17+sig29+3./4.*sig18+3./2.*sig11
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= sig30/sigl)  then
!  D+ L K0
             id1=1121           ! d+
             id2=2130           ! l
             id3=230            ! k0
          else
!  D0 L K+
             id1=1221           ! d0
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig31/sigs)  then
!  D+ S0 K0
             id1=1121           ! d+
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (sig31+sig21)/sigs)  then
!  D+ S- K+
             id1=1121           ! d+
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= (sig31+sig21+sig17)/sigs)   then
!  D0 S0 K+
             id1=1221           ! d0
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= (sig31+sig21+sig17+sig29))  then
!  D0 S+ K0
             id1=1221           ! d0
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= (sig31+sig21+sig17+sig29+3./4.*sig18))  then
!  D- S+ K+
             id1=2221           ! d-
             id2=1130           ! s+
             id3=130            ! k+
          else
!  D++ S- K0
             id1=1111           ! d++
             id2=2230           ! s-
             id3=230            ! k0
          endif
       endif
    elseif(ida == 1121.and.idb == 1220) then
! D+ n ==>DYK
       sigl=sig15+sig30
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig17+sig29+sig31+sig21+sig11+3./4.*sig18
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= sig15/sigl)  then
!  D+ L K0
             id1=1121           ! d+
             id2=2130           ! l
             id3=230            ! k0
          else
!  D0 L K+
             id1=1221           ! d0
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig17/sigs)  then
!  D+ S0 K0
             id1=1121           ! d+
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (sig17+sig29)/sigs)  then
!  D+ S- K+
             id1=1121           ! d+
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= (sig17+sig29+sig31)/sigs)  then
!  D0 S0 K+
             id1=1221           ! d0
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= (sig17+sig29+sig31+sig21)/sigs)  then
!  D0 S+ K0
             id1=1221           ! d0
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= (sig17+sig29+sig31+sig21+sig11)/sigs)  then
!  D- S+ K+
             id1=2221           ! d-
             id2=1130           ! s+
             id3=130            ! k+
          else
!  D++ S- K0
             id1=1111           ! d++
             id2=2230           ! s-
             id3=230            ! k0
          endif
       endif
    elseif(ida == 2221.and.idb == 1120) then
! D- P ==>DYK
       sigl=3./4.*sig30+sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig11+3./4.*sig31+3./4.*sig21+sig7+sig2
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 3./4.*sig30/sigl)  then
!  D0 L K0
             id1=1221           ! d0
             id2=2130           ! l
             id3=230            ! k0
          else
!  D- L K+
             id1=2221           ! d-
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig11/sigs)  then
!  D+ S- K0
             id1=1121           ! d+
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd <= (sig11+3./4.*sig31)/sigs)  then
!  D0 S0 K0
             id1=1221           ! d0
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (sig11+3./4.*sig31+3./4.*sig21)/sigs)  then
!  D0 S- K+
             id1=1221           ! d0
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= (sig11+3./4.*sig31+3./4.*sig21+sig7)/sigs) then
!  D- S+ K0
             id1=2221           ! d-
             id2=1130           ! s+
             id3=230            ! k0
          else
!  D- S0 K+
             id1=2221           ! d-
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1221.and.idb == 1220) then
! D0 n ==>DYK
       sigl=sig17+3./4.*sig30
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig11+sig15+sig29+3./4.*sig21+3./4.*sig31
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= sig17/sigl)  then
!  D0 L K0
             id1=1221           ! d0
             id2=2130           ! l
             id3=230            ! k0
          else
!  D- L K+
             id1=2221           ! d-
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig11/sigs)  then
!  D+ S- K0
             id1=1121           ! d+
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd <= (sig11+sig15)/sigs)  then
!  D0 S0 K0
             id1=1221           ! d0
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (sig11+sig15+sig29)/sigs)  then
!  D0 S- K+
             id1=1221           ! d0
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= (sig11+sig15+sig29+3./4.*sig21)/sigs)  then
!  D- S+ K0
             id1=2221           ! d-
             id2=1130           ! s+
             id3=230            ! k0
          else
!  D- S0 K+
             id1=2221           ! d-
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    elseif(ida == 2221.and.idb == 1220) then
! D- N ==>DYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig18+sig2+sig7
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  D- L K0
          id1=2221           ! d-
          id2=2130           ! l
          id3=230            ! k0
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig18/sigs)  then
!  D0 S- K0
             id1=1221           ! d0
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd <= (sig18+sig2)/sigs)  then
!  D- S0 K0
             id1=2221           ! d-
             id2=1230           ! s0
             id3=230            ! k0
          else
!  D- S- K+
             id1=2221           ! d-
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    else
       write(16,*) 'dndyk: ida,idb=',ida,idb
       write( *,*) 'dndyk: ida,idb=',ida,idb
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  dnnyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels Delta + N ==> Y + N + K,
!     K=K+,K0
!
    real(real64) ::  mn,mk,ml,ms
    data  mn/0.939/,mk/0.494/,ml/1.115/,ms/1.192/
! Parameters for D++ n ==>p  L  K+ cross section
    data  a1/8.337/, b1/2.227/,c1/2.511/
! Parameters for D-  p ==>n  S- K+ cross section
    data a22/52.72/, b22/2.799/,c22/6.303/
    sig=0.
    ss0l=mn+ml+mk
    if(ss <= ss0l)                  return
    ss0s=mn+ms+mk
    sig1 =sibir(a1,b1,c1,ss,ss0l)
    sig22=sibir(a22,b22,c22,ss,ss0s)
    if(ida == 1111.and.idb == 1220) then
! D++ n ==>NYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=(0.5+0.5+1.)*sig22
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  p L K+
          id1=1120           ! p
          id2=2130           ! l
          id3=130            ! k+
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.25)  then
!  n S+ K+
             id1=1220           ! n
             id2=1130           ! s+
             id3=130            ! k+
          elseif(rnd <= 0.5)  then
!  p S0 K+
             id1=1120           ! p
             id2=1230           ! s0
             id3=130            ! k+
          else
!  p S+ K0
             id1=1120           ! p
             id2=1130           ! s+
             id3=230            ! k0
          endif
       endif
    elseif(ida == 1111.and.idb == 1120) then
! D++ p ==>NYK
       sig=0.5*sig22
       if(sig < 1.d-10)         return
!  p S+ K+
       id1=1120           ! p
       id2=1130           ! s+
       id3=130            ! k+
    elseif(ida == 1121.and.idb == 1220) then
! D+ n ==>NYK
       sigl=2./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig22
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.5)  then
!  p L K0
             id1=1120           ! p
             id2=2130           ! l
             id3=230            ! k0
          else
!  n L K+
             id1=1220           ! n
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 1./6.)  then
!  p S0 K0
             id1=1120           ! p
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= 2./6.)  then
!  n S0 K+
             id1=1220           ! n
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= 2./3.)  then
!  n S+ K0
             id1=1220           ! n
             id2=1130           ! s+
             id3=230            ! k0
          else
!  p S- K+
             id1=1120           ! p
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1121.and.idb == 1120) then
! D+ p ==>NYK
       sigl=1./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=1./3.*sig22
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  p L K+
          id1=1120           ! p
          id2=2130           ! l
          id3=130            ! k+
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.333333)  then
!  p S0 K+
             id1=1120           ! p
             id2=1230           ! s0
             id3=130            ! k+
          else
!  p S+ K0
             id1=1120           ! p
             id2=1130           ! s+
             id3=230            ! k0
          endif
       endif
    elseif(ida == 1221.and.idb == 1220) then
! D0 n ==>NYK
       sigl=1./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig22
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  n L K0
          id1=1220           ! n
          id2=2130           ! l
          id3=230            ! k0
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 1./6.)  then
!  n S0 K0
             id1=1220           ! n
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= 2./3.)  then
!  p S- K0
             id1=1120           ! p
             id2=2230           ! s-
             id3=230            ! k0
          else
!  n S- K+
             id1=1220           ! n
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    elseif(ida == 2221.and.idb == 1120) then
! D- p ==>NYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=2.*sig22
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  n L K0
          id1=1220           ! n
          id2=2130           ! l
          id3=230            ! k0
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.25)  then
!  n S0 K0
             id1=1220           ! n
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= 0.5)  then
!  p S- K0
             id1=1120           ! p
             id2=2230           ! s-
             id3=230            ! k0
          else
!  n S- K+
             id1=1220           ! n
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    elseif(ida == 2221.and.idb == 1220) then
! D- n ==>NYK
       sig=0.5*sig22
       if(sig < 1.d-10)         return
!  n S- K0
       id1=1220           ! n
       id2=2230           ! s-
       id3=230            ! k0
    elseif(ida == 1221.and.idb == 1120) then
! D0 p ==>NYK
       sigl=2./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig22
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.5)  then
!  n L K+
             id1=1220           ! n
             id2=2130           ! l
             id3=130            ! k+
          else
!  p L K0
             id1=1120           ! p
             id2=2130           ! l
             id3=230            ! k0
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 1./6.)  then
!  n S0 K+
             id1=1220           ! n
             id2=1230           ! s0
             id3=130            ! k+4
          elseif(rnd <= 1./3.)  then
!  p S0 K0
             id1=1120           ! p
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= 2./3.)  then
!  n S+ K0
             id1=1220           ! n
             id2=1130           ! s+
             id3=230            ! k0
          else
!  p S- K+
             id1=1120           ! p
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    else
       write(16,*) 'dnnyk: ida,idb=',ida,idb
       write( *,*) 'dnnyk: ida,idb=',ida,idb
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  nndyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels N + N ==> Y + Delta + K,
!     K=K+,K0
!
    real(real64) ::  md,mk,ml,ms
    data  md/1.236/,mk/0.494/,ml/1.115/,ms/1.192/
! Parameters for p p ==>D++ L  K0 cross section
    data  a1/6.166/, b1/2.842/,c1/1.960/
! Parameters for p p ==>D++ S- K+ cross section
    data a10/10.00/, b10/2.874/,c10/2.543/
    sig=0.
    ss0l=md+ml+mk
    if(ss <= ss0l)                  return
    ss0s=md+ms+mk
    sig1 =sibir(a1,b1,c1,ss,ss0l)
    sig10=sibir(a10,b10,c10,ss,ss0s)
    if(ida == 1120.and.idb == 1120) then
! p  p ==>DYK
       sigl=4./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=(1./3.+1./6.+1./3.+1.+1./2.)*sig10
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.75)  then
!  D++ L K0
             id1=1111           ! d++
             id2=2130           ! l
             id3=230            ! k0
          else
!  D+ L K+
             id1=1121           ! d+
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 2./14.)  then
!  D+ S+ K0
             id1=1121           ! d+
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= 3./14.)  then
!  D+ S0 K+
             id1=1121           ! d+
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= 5./14.)  then
!  D0 S+ K+
             id1=1221           ! d0
             id2=1130           ! s+
             id3=130            ! k+
          elseif(rnd <= 11./14.)  then
!  D++ S- K+
             id1=1111           ! d++
             id2=2230           ! s-
             id3=130            ! k+
          else
!  D++ S0 K0
             id1=1111           ! d++
             id2=1230           ! s0
             id3=230            ! k0
          endif
       endif
    elseif(ida == 1120.and.idb == 1220.or. &
            & ida == 1220.and.idb == 1120) then
! n p(p n) ==>DYK
       sigl=2./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=(1./3.+1./6.+1./3.+1./6.+1./2.+1./3.)*sig10
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.5)  then
!  D+ L K0
             id1=1121           ! d+
             id2=2130           ! l
             id3=230            ! k0
          else
!  D0 L K+
             id1=1221           ! d0
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 2./11.)  then
!  D++ S- K0
             id1=1111           ! d++
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd <= 3./11.)  then
!  D+ S0 K0
             id1=1121           ! d+
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= 5./11.)  then
!  D0 S+ K0
             id1=1221           ! d0
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= 6./11.)  then
!  D0 S0 K+
             id1=1221           ! d0
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= 9./11.)  then
!  D- S+ K+
             id1=2221           ! d-
             id2=1130           ! s+
             id3=130            ! k+
          else
!  D+ S- K+
             id1=1121           ! d+
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1220.and.idb == 1220) then
! n  n ==>DYK
       sigl=4./3.*sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=(1./6.+1./3.+1./6.+1.+0.5)*sig10
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.75)  then
!  D- L K+
             id1=2221           ! d-
             id2=2130           ! l
             id3=130            ! k+
          else
!  D0 L K0
             id1=1221           ! d0
             id2=2130           ! l
             id3=230            ! k0
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= 1./13.)  then
!  D+ S- K0
             id1=1121           ! d+
             id2=2230           ! s-
             id3=230            ! k0
          elseif(rnd <= 3./13.)  then
!  D0 S- K+
             id1=1221           ! d0
             id2=2230           ! s-
             id3=130            ! k+
          elseif(rnd <= 4./13.)  then
!  D0 S0 K0
             id1=1221           ! d0
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= 10./13.)  then
!  D- S+ K0
             id1=2221           ! d-
             id2=1130           ! s+
             id3=230            ! k0
          else
!  D- S0 K+
             id1=2221           ! d-
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    else
       write(16,*) 'nndyk: ida,idb=',ida,idb
       write( *,*) 'nndyk: ida,idb=',ida,idb
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  nnnyk(ida,idb,id1,id2,id3,ss,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     calculation of channels N + N ==> Y + N + K,
!     K=K+,K0
!
    real(real64) ::  mn,mk,ml,ms
    data  mn/0.939/,mk/0.494/,ml/1.115/,ms/1.192/
! Parameters for p p ==>p L  K+ cross section
    data  a1/1.879/, b1/2.176/,c1/5.264/
! Parameters for p n ==>n L  K+ cross section
    data  a2/2.812/, b2/2.121/,c2/4.893/
! Parameters for p p ==>n S+ K+ cross section
    data  a5/1.466/, b5/2.743/,c5/3.271/
! Parameters for p p ==>p S+ K0 cross section
    data  a6/7.079/, b6/2.760/,c6/8.164/
! Parameters for p p ==>p S0 K+ cross section
    data  a7/5.321/, b7/2.753/,c7/8.510/
! Parameters for p n ==>p S- K+ cross section
    data  a8/11.02/, b8/2.782/,c8/7.674/
! Parameters for p n ==>n S0 K+ cross section
    data  a9/6.310/, b9/2.773/,c9/7.820/
    sig=0.
    ss0l=mn+ml+mk
    if(ss <= ss0l)                  return
    ss0s=mn+ms+mk
    sig1 =sibir(a1,b1,c1,ss,ss0l)
    sig2 =sibir(a2,b2,c2,ss,ss0l)
    sig5 =sibir(a5,b5,c5,ss,ss0s)
    sig6 =sibir(a6,b6,c6,ss,ss0s)
    sig7 =sibir(a7,b7,c7,ss,ss0s)
    sig8 =sibir(a8,b8,c8,ss,ss0s)
    sig9 =sibir(a9,b9,c9,ss,ss0s)
    if(ida == 1120.and.idb == 1120) then
! p  p ==>NYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig5+sig6+sig7
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  p L K+
          id1=1120           ! p
          id2=2130           ! l
          id3=130            ! k+
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig5/sigs)  then
!  n S+ K+
             id1=1220           ! n
             id2=1130           ! s+
             id3=130            ! k+
          elseif(rnd <= (sig5+sig6)/sigs)  then
!  p S+ K0
             id1=1120           ! p
             id2=1130           ! s+
             id3=230            ! k0
          else
!  p S0 K+
             id1=1120           ! p
             id2=1230           ! s0
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1120.and.idb == 1220.or. &
            & ida == 1220.and.idb == 1120) then
! n p(p n) ==>NYK
       sigl=2.*sig2
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig8+sig9+sig9+sig8
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
          rnd=rndm(-1.0_real64)
          if(rnd <= 0.5)  then
!  p L K0
             id1=1120           ! p
             id2=2130           ! l
             id3=230            ! k0
          else
!  n L K+
             id1=1220           ! n
             id2=2130           ! l
             id3=130            ! k+
          endif
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig8/sigs)  then
!  n S+ K0
             id1=1220           ! n
             id2=1130           ! s+
             id3=230            ! k0
          elseif(rnd <= (sig8+sig9)/sigs)  then
!  n S0 K+
             id1=1220           ! n
             id2=1230           ! s0
             id3=130            ! k+
          elseif(rnd <= (sig8+sig9+sig9)/sigs)  then
!  p S0 K0
             id1=1120           ! p
             id2=1230           ! s0
             id3=230            ! k0
          else
!  p S- K+
             id1=1120           ! p
             id2=2230           ! s-
             id3=130            ! k+
          endif
       endif
    elseif(ida == 1220.and.idb == 1220) then
! n  n ==>NYK
       sigl=sig1
       if(ss <= ss0s)  then
          sigs=0.
       else
          sigs=sig7+sig6+sig5
       endif
       sig=sigl+sigs
       if(sig < 1.d-10)         return
       if(rndm(-1.0_real64) <= sigl/sig)  then
!  n L K0
          id1=1220           ! n
          id2=2130           ! l
          id3=230            ! k0
       else
          rnd=rndm(-1.0_real64)
          if(rnd <= sig7/sigs)  then
!  n S0 K0
             id1=1220           ! n
             id2=1230           ! s0
             id3=230            ! k0
          elseif(rnd <= (sig7+sig6)/sigs)  then
!  n S- K+
             id1=1220           ! n
             id2=2230           ! s-
             id3=130            ! k+
          else
!  p S- K0
             id1=1120           ! p
             id2=2230           ! s-
             id3=230            ! k0
          endif
       endif
    else
       write(16,*) 'nnnyk: ida,idb=',ida,idb
       write( *,*) 'nnnyk: ida,idb=',ida,idb
    endif
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  sibir(a,b,c,ss,ss0)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculation of Sibirtsev's functional form of
!     strange proguction cross sections in B+B==>Y+B+K
!
    sibir=0.0
    if(ss <= ss0)   return
    s=ss*ss
    s0=ss0*ss0
    sibir=a*((s/s0-1.)**b)*(s0/s)**c
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine elexq(v,u,tin1,partin,ipatin,ipatne, &
       & mv,np,l,ms,mq,ksi,me)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     BLOCK OF CALCULATION OF PARTICLE CHARACTERISTICS IN ELASTIC AND
!     CHARGE EXCHANGE SCATTERING.
!     Calls: PINDEL,cduarteq, COSELQ, STREX, ABELQ
!
    dimension pmemo(9,5999),imemo(5,5999),v(3),partin(9),ipatin(5), &
         & ipatne(5),pist(3),pnst(3)
    common/memorylaq/pmemo,imemo
    common/isob2/isob2
    common/ncasca/ncas,ncpri
    is=ipatin(3)
    ns=ipatne(3)
    if (ipatin(2)) 10,11,10
10  cmi = 0.14
    cmn = 0.94
    go to 12
11  cmi = partin(9)
    cmn = 0.94
    go to 12
12  sex=croseg(l,ms,mq,ksi,2,tin1,partin(9),ipatin(5))
    sel=croseg(l,ms,mq,ksi,1,tin1,partin(9),ipatin(5))
    betaex=sex/(sex+sel)
    drnd=rndm(-1.0_real64)
    if(drnd-betaex)14,13,13
13  ie = ipatin(1)
    ne = ipatne(1)
!  kkg 10/28/03
    if(isob2 == 0.or.l.ne.0)  go  to  113
    iso=0
    if(ipatin(4) == 0.and.ipatin(3) == 0.and.ipatin(5) == 0) &
         & call  pindel(tin1,u,mv,ie,ne,np,iso)
    if(iso == 1)  return
!  kkg 10/22/03
113 continue
    if(mq == 2.and.tin1 < 2.0)  then
       ctsti=cduarteq(tin1, ipatin(1), ipatne(1))
    else
       ctsti = coselq(l,mq,ksi,tin1,partin(9))
    endif
    go to 18
!  kkg 10/22/03
14  if(ipatin(5).ne.0.and.ipatin(4) == 1)    go  to  117
    if(ipatin(3) == 0)  go  to  114
    if(ipatin(3).ne.0) &
         & call strex(ipatin(3),ipatin(1),ipatne(1),cmi,cmn,is,ns,ie,ne)
    go  to  17
114 if (ipatin(1)) 15,16,15
15  ie = 0
    ne = me
    go to 17
16  ne = 1-ipatne(1)
    ie = me-ne
    go to 17
17  ctsti = cosexq (l,tin1,partin(9))
    go to 18
117 ne=1-ipatne(1)
    ie=me-ne
    go  to  113
18  fisti = 0.
    call abelq(partin,v,u,pist,pnst,ctsti,fisti,cmi,cmn)
    if(mv-5997) 19,19,20
20  np = 0
    write(16,21)
21  format (45x,32h memorylaq is exceeded in cascad)
    return
19  continue
    pmemo(1,mv+3)=0.
    pmemo(2,mv+3)=0.
    pmemo(3,mv+3)=0.
    pmemo(4,mv+3) = pist(1)
    pmemo(5,mv+3) = pist(2)
    pmemo(6,mv+3) = pist(3)
    pmemo(7,mv+3)=0.
    pmemo(9,mv+3)=cmi
    imemo(1,mv+3)=ie
    imemo(2,mv+3)=0
    imemo(3,mv+3) = is
    imemo(4,mv+3) = ipatin(4)
    imemo(5,mv+3) = ipatin(5)
    pmemo(1,mv+1)=0.
    pmemo(2,mv+1)=0.
    pmemo(3,mv+1)=0.
    pmemo(4,mv+1) = pnst(1)
    pmemo(5,mv+1) = pnst(2)
    pmemo(6,mv+1) = pnst(3)
    pmemo(7,mv+1)=0.
    pmemo(9,mv+1)=cmn
    imemo(1,mv+1)=ne
    imemo(2,mv+1)=0
    imemo(3,mv+1) = ns
    imemo(4,mv+1) = 1
    imemo(5,mv+1) = ipatne(5)
    np = 2
!
    if(ncas > ncpri)  then
       write(16,100) mv+1, tin1, pist,cmn,(imemo(k,mv+1),k=1,5), &
            & mv+3, ctsti,pnst,cmi,(imemo(k,mv+3),k=1,5)
       write( *,100) mv+1, tin1, pist,cmn,(imemo(k,mv+1),k=1,5), &
            & mv+3, ctsti,pnst,cmi,(imemo(k,mv+3),k=1,5)
100    format(' elexq:',i5,5e11.4,4i3,i10/ &
            & 6x,i5,5e11.4,4i3,i10)
    endif
!
    return
  end
!
! ********************************************************************
!
  double precision function cduarteq(tin1,ie1,ie2)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!     angular distribution simulation for n+p and p+p
!     for energies < 2GeV using Duarte's approximations:
!     n+p : ds/dom ~ exp(B*t)+a*exp(B*u)+c*exp(alc*u)
!     P+p:  ds/dom ~ exp(B*t)+a*exp(B1*t)
!     with Mandelstam's variables t,u
!     H.Duarte http://www.fjfi.cvut.cz/con_adtt99/papers/Mo-o-c17.pdf
    real(real64) ::  mn
    data mn/0.939/
    data  pi /3.1415926536d0/
!
    elab=tin1*1000.
    ps2=tin1*mn/2.
    if(ie1 == ie2)   then      ! p +p or n + n
       if(elab <= 300.)      then
          b = 0.
          a = 0.2
          b1= 0.
       elseif(elab <= 670.)  then
          b = 9.87e-8*(elab-300.)**3
          a = 0.2
          b1= 0.
       elseif(elab <= 1100.)  then
          b = 4.56e-3*(elab-670.)+4.76
          a = 97.02e+3*exp(-2.0e-2*elab)+0.053
          b1= 9.72e-8*exp(-5.0e-3*(elab-670.))*(elab-670.)**3
       else
          b = 7.4/(1.+3.0e+5/(elab-300.)**2.23)
          a = 0.28*exp(-1.5e-3*elab)
          b1= 1.94*exp(-7.0e-4*elab)
       endif
       if(b <= 0.0001)  then
          cn1 = 4.*pi
       else
          cn1 = 2.*pi*  (1.-exp(-2.*ps2*b ))/(ps2*b)
       endif
       if(b1 <= 0.0001) then
          cn2 = 4.*pi*a
       else
          cn2=  2.*pi*a*(1.-exp(-2.*ps2*b1))/(ps2*b1)
       endif
       cn = cn1+cn2
       c1 = cn1/cn
       c2 = cn2/cn
       r1=rndm(-1.0_real64)
       tm=-2.0*ps2
       if(r1 <= c1)   then
          r2=rndm(-1.0_real64)
          if(b <= 0.0001)  then
             t= tm*r2
          else
             t=log(1.0-r2*(1.0-exp(tm*b)))/b
          endif
       else
          r2=rndm(-1.0_real64)
          if(b1 <= 0.0001)  then
             t= tm*r2
          else
             t=log(1.0-r2*(1.0-exp(tm*b1)))/b1
          endif
       endif
       cts=1.0-t/tm
       if(rndm(-1.0_real64) <= 0.5) then
          cduarteq = cts
       else
          cduarteq =-cts
       endif
!
    elseif((ie1 == 0.and.ie2 == 1) &           ! n + p &
       & .or.(ie1 == 1.and.ie2 == 0))   then     ! p + n
       b = 25.0*exp(-6.0e-3*elab)+2.0e-3*elab+2.8
       a = (3.0e-3*elab+0.1)*exp(1.55-4.9e-3*elab)+1.0e-4*elab
       c = 6.0e-5*elab*elab*exp(-6.0e-3*elab)+15.0e-5*elab
       if(elab < 100.)  alc=330.-elab
       if(elab >= 100.)  alc=80.+15.0e+3/elab
       cn1 = pi/(ps2*b)*(1.-exp(-4.*ps2*b))
       cn2 = cn1*a
       cn3 = pi/(ps2*alc)*c*(1.-exp(-4.*ps2*alc))
       cn = cn1+cn2+cn3
       c1=cn1/cn
       c2=cn2/cn
       c3=cn3/cn
       tm=-4.0*ps2
       um=tm
       r1 = rndm(-1.0_real64)
       if(r1 <= c1)         then
          r2=rndm(-1.0_real64)
          if(b <= 0.0001)     then
             t= tm*r2
          else
             t=log(1.0-r2*(1.0-exp(tm*b)))/b
          endif
          cts = 1.0-2.0*t/tm
       elseif(r1 <= (c1+c2)) then
          r2=rndm(-1.0_real64)
          if(b <= 0.0001)     then
             u= um*r2
          else
             u=log(1.0-r2*(1.0-exp(um*b)))/b
          endif
          cts =-1.0+2.0*u/um
       else
          r2=rndm(-1.0_real64)
          if(alc <= 0.0001)     then
             u= um*r2
          else
             u=log(1.0-r2*(1.0-exp(um*alc)))/alc
          endif
          cts =-1.0+2.0*u/um
       endif
       if(ie1 == 1.and.ie2 == 0) then
!  parametrization is for n + p case
          cduarteq=-cts
       else
          cduarteq= cts
       endif

    else
       write(*,*) ' cduarteq: ie1 and ie2 are wrong, ie1,ie2=',ie1,ie2
       stop
    endif
    if(abs(cduarteq) > 1.0)         then
       cduarteq=sign(1.d0,cduarteq)
    elseif(abs(cduarteq) < 1.0d-10) then
       cduarteq=0.0
    endif
    return
  end
!
! ********************************************************************
!

  subroutine strex(is1,ie1,ne1,cmi,cmn,is,ns,ie,ne)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculation of  channels a) charge exchange K + N ==> K' + N'
!     b) strange exchange AKA + N ==> pi + Y,
!     where K=(K+,K0), AKA=(K-,AK0)
!
    if(is1 < 0)  go  to  10
    ie=1-ie1
9   is=is1
    ns=0
    cmi=0.492
    cmn=0.940
    go  to  23
10  if((ie1+ne1) == 0)  go  to  17
    rnd4=4.*rndm(-1.0_real64)
    nbr=intg(rnd4)+1
    go  to  (11,12,13,14,14),nbr
11  ie=1
    go  to  15
12  ie=0
    go  to  16
13  ie=1
    go  to  16
14  ie=-1
    go  to   9
15  is=0
    ns=-1
    cmi=0.140
    cmn=1.116
    go  to  23
16  is=0
    ns=-1
    cmi=0.140
    cmn=1.189
    go  to  23
17  rnd5=5.*rndm(-1.0_real64)
    nbr=intg(rnd5)+1
    go  to  (18,19,20,21,22,22),nbr
18  ie=0
    go  to  15
19  ie=0
    go  to  16
20  ie=-1
    go  to  16
21  ie=1
    go  to  16
22  if(ie1 ==  0)  ie=-1
    if(ie1 == -1)  ie= 0
    go  to  9
23  ne=ie1+ne1-ie
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine pindel(t1,dm,mv,ie,ne,np,iso)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)(a-h,o-z), integer(int32) (i-n)
!
!   forms the Delta from pi+N system
!
    common/memorylaq/pme(9,5999),ime(5,5999)
    iso=0
    if(t1 > 0.5)                  return
    if(dm <= 1.082)                return
    call  wmd(dm,tau0,fmd)
    drnd=rndm(-1.0_real64)
    if(drnd > fmd) return
    pme(1,mv+1)=0.
    pme(2,mv+1)=0.
    pme(3,mv+1)=0.
    pme(4,mv+1)=0.
    pme(5,mv+1)=0.
    pme(6,mv+1)=0.
    pme(7,mv+1)=0.
    pme(9,mv+1)=dm
    ime(1,mv+1)=ie+ne
    ime(2,mv+1)=0
    ime(3,mv+1)=0
    ime(4,mv+1)=1
    ime(5,mv+1)=intg(1000.*tau0)
    if(ime(5,mv+1) == 0)  ime(5,mv+1)=1
    np=1
    iso=1
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  tmatur(np,v,m,p1,p2,mv,ps,tl,tmat)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!  calculation of formation (maturity) time(fm/c) of produced particle
!
    real(real64) ::   mmat,mmes,mbar,mlid
    common/mmatur/mmes,mbar
    common/inttyp/ ityp
    common/ithea/ithea(11)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/iact/ iact
    common/cslid/clider(5999)
    common/cmali/cmali
    common/actim/tint
    dimension v(3),p1(9),p2(9),ps(3),rl(3)
!
!     ITYP=1(ABSORP),2(ELEX),3(BINEL),4(HEINEN)
!
    tmat=0.
    if(np <= 2.or.ityp <= 3.or. &
         & (ityp == 4.and.ithea(4) == 1))  then
!   formation time TMAT=0 for  2-body channel
!                              absorption
!                              elastic scattering
!                              low energy inelastic collision
!                              high energy elastic collision
       clider(m)=1.
    else
       if(iact <= 2) then
          mlid=mbar
          if(imemo(4,m) == 0)              mmat=mmes
          if(imemo(4,m).ne.0)              mmat=mbar
!  FOR ANTI-BARYON
          if(imemo(4,m) < 0)              mmat=mmes
!
!           IF(CLIDER(M).GT.0.)              MMAT=MLID !!12.10.97/12.03.98
!
          tau0=mmat/(5.06*pmemo(9,m))  !! 05.06.96,  04/11/07
!         TAU0=MMAT/(5.06*0.14)        !! 14.09.97
          gl=1.+pmemo(8,m)/pmemo(9,m)
          v2=v(1)**2+v(2)**2+v(3)**2
          u=(p1(8)+p1(9)+p2(8)+p2(9))*sqrt(1.-v2)
          x=abs(ps(3))/(u/2.)
!         FX=FXTMAT(X)
          fx=1.0d0
          tau0l=gl*tau0*fx
          brnd=rndm(-1.0_real64)
          af=-log(brnd)               !! 05.06.96
!           AF=1.                        !! 27.02.97
          tmat=tau0l*af
       else
!  For IACT=3
          rtau=tl-tint
          tmat=tl
          if(clider(m) > cmali)   then
             tmat=0.
             clider(m)=1.
          endif
       endif
       if(tmat < 0.) then
          write(16,100) ityp,m,imemo(5,m),(pmemo(k,m),k=8,9),x,tmat
100       format(2x,'tmatur:ityp=',i2,2x,'m=',i5,2x, &
               & 'im(5)=',i15,2x,'t=',f8.3,2x,'mass=',f8.3,2x, &
               & 'x=',e13.6,2x,'tmat=',e13.6)
          tmat=0.
       endif
    endif
!     CALL  HTFORL(M,TMAT)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine kinemr(v,m,rl,tl)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   calculation of proper and observer's system time of
!   produced particle M
!
    common/memorylaq/pme(9,5999),ime(5,5999)
    dimension v(3),rl(3)
! PROPER TIME
    tp2=pme(7,m)**2-pme(1,m)**2-pme(2,m)**2-pme(3,m)**2
    tp=sqrt(abs(tp2))
! TIME IN OBSERVER SYSTEM
    g=(pme(8,m)+pme(9,m))/pme(9,m)
    tl=g*tp
    rl(1)=0.
    rl(2)=0.
    rl(3)=0.
!     CALL  HTFORP(M)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine vmnspq (partin,ipatin,u,mv,np,ith,mq,tin1,lp)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculate secondary particle number and determine absolute
!     values of momenta in inelastic (pion production) interaction.
!
!   Called by: BINELQ
!
!   Calls: JTYPBQ PMOMQ
!
!
    dimension pmemo(9,5999),imemo(5,5999),partin(9),ipatin(5)
    common/memorylaq/pmemo,imemo
    kh = 0
    lp = 0
    am3=partin(9)
    if(ipatin(5).ne.0.and.ipatin(4) == 1)  am3=0.940
10  u1 = u
    lambda = 1
11  ltemp = mv+lambda
    if(ltemp-5999) 37,37,38
38  np = 0
    write(16,39)
39  format (25x,'memorylaq is exceeded in cascan')
    return
37  continue
    if (lambda-1) 12,12,13
12  pmemo(9,mv+1)=0.94
    imemo(2,mv+1)=0
    imemo(3,mv+1)=0
    imemo(4,mv+1) = 1
    imemo(5,mv+1) = 0
    go to 16
13  if (lambda-3) 15,14,15
14  pmemo(9,mv+3) = am3
    imemo(2,mv+3) = 0
    imemo(3,mv+3) = ipatin(3)
    imemo(4,mv+3) = ipatin(4)
    imemo(5,mv+3) = 0
    go to 16
15  pmemo(9,ltemp) = 0.14
    imemo(2,ltemp) = 0
    imemo(3,ltemp) = 0
    imemo(4,ltemp) = 0
    imemo(5,ltemp) = 0
    go to 16
16  jb = jtypbq(ith,mq,lambda)
    pmemo(8,ltemp) = pmomq(jb,tin1)
    el = sqrt (pmemo(8,ltemp)**2+pmemo(9,ltemp)**2)
    deltu = u1-el
    if (lambda-2) 23,17,23
17  if (deltu-am3)       35,35,18
18  if (ith) 19,22,19
19  pmemo(8,mv+3) = sqrt (deltu**2-am3**2)
    pmemo(9,mv+3) = am3
    imemo(2,mv+3) = 0
    imemo(3,mv+3) = ipatin(3)
    imemo(4,mv+3) = ipatin(4)
    imemo(5,mv+3) = 0
    if (pmemo(8,mv+1)-pmemo(8,mv+2)-pmemo(8,mv+3)) 20,20,35
20  if (pmemo(8,mv+1)-abs(pmemo(8,mv+2)-pmemo(8,mv+3))) 35,35,21
21  np = 3
    return
22  u1 = deltu
    lambda = lambda+1
    go to 11
23  if (deltu-0.14) 24,24,22
24  if (lambda-1) 35,35,25
25  if (lambda-3) 26,35,26
26  el=deltu+el
    pmemo(8,ltemp) = sqrt (el**2-pmemo(9,ltemp)**2)
    np = lambda
    i = 1
28  itemp = mv+i
    c = pmemo(8,itemp)
29  if (np-i) 34,34,30
30  if (c-pmemo(8,itemp+1)) 32,31,31
31  i=i+1
    itemp=itemp+1
    go  to  29
32  i = i+1
    go to 28
34  pmax = c
    sigma = 0.
!----> do 33 i=1,np
    do i=1,np
       itemp = mv+i
       sigma = sigma+pmemo(8,itemp)
33     continue
    end do
    if (2.*pmax-sigma) 27,35,35
27  return
35  kh = kh+1
    if (kh-100) 10,36,36
36  lp = 2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine chinelq (ipatin,l,ms,mq,ksi,np,mv,tin1,me,ipatne,amin)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Determinine secondary particles' charges in inelastic (pion prod.)
!     reaction.
!
!   Called by: BINELQ
!
!   Calls: CROSEG
!
    dimension ipatin(5),ipatne(5),pmemo(9,5999),imemo(5,5999)
    common/memorylaq/pmemo,imemo
    if(ipatin(5).ne.0.and.ipatin(4) == 1)  go  to  41
    if (np-3) 21,10,21
10  spi0 = croseg (l,ms,mq,ksi,4,tin1,amin,ipatin(5))
    sth = croseg (l,ms,mq,ksi,7,tin1,amin,ipatin(5))
    bpi0 = spi0/sth
    bpiex = (spi0+croseg(l,ms,mq,ksi,5,tin1,amin,ipatin(5)))/sth
    temp1 = rndm(-1.0_real64)
    if (temp1-bpi0) 19,11,11
11  if (temp1-bpiex) 20,12,12
12  imemo(1,mv+1)=ipatne(1)-(ipatin(4)-1)*ipatin(1)
    imemo(1,mv+3)=(ipatin(4)-1)*ipatin(1)**2-ipatin(4)*ipatin(1)+1
    go to 18
18  imemo(1,mv+2)=me-imemo(1,mv+1)-imemo(1,mv+3)
    return
19  imemo(1,mv+1) = ipatne(1)
    imemo(1,mv+2) = 0
    imemo(1,mv+3) = ipatin(1)
    return
20  imemo(1,mv+1) = 1-ipatne(1)
    imemo(1,mv+3) = ipatin(1)
    imemo(1,mv+2)=me-imemo(1,mv+1)-imemo(1,mv+3)
    return
21  if (rndm(-1.0_real64) -0.5) 22,22,23
22  imemo(1,mv+1) = 1
    go to 24
23  imemo(1,mv+1) = 0
    go to 24
24  if (mq-1) 28,28, 25
25  if (rndm(-1.0_real64) -0.5) 26,26,27
26  imemo(1,mv+3) = 1
    go to 28
27  imemo(1,mv+3) = 0
    go to 28
28  lambda = 2
29  if (mq-1) 31,31,30
30  if (lambda-3) 31,36,31
31  temp2 = rndm(-1.0_real64)
    mtemp = mv+lambda
    if (temp2-1./3.) 33,32,32
32  if (temp2-2./3.) 34,35,35
33  imemo(1,mtemp) = 1
    go to 36
34  imemo(1,mtemp) = 0
    go to 36
35  imemo(1,mtemp) = -1
    go to 36
36  if (lambda-np) 37,38,38
37  lambda = lambda+1
    go to 29
38  sigq = 0.
!----> do 39 i=1,np
    do i=1,np
       itemp = mv+i
       sigq = sigq+imemo(1,itemp)
39     continue
    end do
    if (me-sigq) 21,40,21
40  return
41  ie1=1
    if(rndm(-1.0_real64) > 0.5)  ie1=0
    br=rndm(-1.0_real64)
    if(br <= 0.33333333)                      ie2=-1
    if(br > 0.33333333.and.br <= 0.66666667) ie2= 0
    if(br > 0.66666667)                      ie2= 1
    ie3=me-ie1-ie2
    if(ie3 < 0.or.ie3 > 1)  go  to  41
    imemo(1,mv+1)=ie1
    imemo(1,mv+2)=ie2
    imemo(1,mv+3)=ie3
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  potbar(nu,x,y,z)
!
!  calculated as U=U0*rho(r)/rho(0) with U0=200 MeV according to the paper
!  by Ye.S.Golubeva, A.S.Iljinov at al. Nucl.Phys. A483 (1988) 539.
!

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    data pbar0 /0.200/
    potbar=0.
    if(nu == 1) then
       a=a1
       c=c1
       d=d1
       an=anucl1
       rmax=rm1
    else
       a=a2
       c=c2
       d=d2
       an=anucl2
       rmax=rm2
    endif
    r=sqrt(x**2+y**2+z**2)
    if(r > 1.5*rmax)  return
    if(an <= 10.) then
       roro0=exp(-(r/a)**2)
    else
       roro0=(1.+exp(-a/c)) / (1.+exp((r-a)/c))
    endif
    potbar=pbar0*roro0
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  potenq(p,ip,a,c,d,tf0,vpion,eps)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Calculation of particle nuclear potential.
!
!
    dimension   ip(5),p(9)
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    if(ip(5).ne.0)  go  to  10
    if(ip(3))  10,11,10
10  potenq=0.
    return
11  if(ip(4))   13,12,13
12  potenq = vpion
    return
13  if(abs(a-a1)-.001) 100,101,101
100 an=anucl1
    rmax=rm1
    nu=1
    go to 105
101 an=anucl2
    rmax=rm2
    nu=2
105 continue
    r=sqrt(p(1)**2+p(2)**2+p(3)**2)/rmax
!     R=RPOTEN(R)
    if(r-1.5)  14,10,10
14  r=r*rmax
    if(an-10.) 106,106,107
106 tf=tf0*exp(-(2./3.)*(r/a)**2)
    go to 108
107 tf=tf0*(((1.+exp(-a/c))/(1.+exp((r-a)/c)))**0.6666667)
108 continue
!     !!! 20.06.1995
    if(nu == 1) tf=tf*(an1/anucl1)**(2./3.)
    if(nu == 2) tf=tf*(an2/anucl2)**(2./3.)
!
    potenq=tf+eps
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine recul(nu,px,py,pz,x,y,z)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   Change of linear and angular momenta of projectile (NU=1)
!   or target (NU=2) by recul px,py,pz
!
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    if(nu-1) 10,10,11
10  pnucl1(1)=pnucl1(1)+px
    pnucl1(2)=pnucl1(2)+py
    pnucl1(3)=pnucl1(3)+pz
    amnuc1(1)=amnuc1(1)+z*py-y*pz
    amnuc1(2)=amnuc1(2)+x*pz-z*px
    amnuc1(3)=amnuc1(3)+y*px-x*py
    go to 12
11  pnucl2(1)=pnucl2(1)+px
    pnucl2(2)=pnucl2(2)+py
    pnucl2(3)=pnucl2(3)+pz
    amnuc2(1)=amnuc2(1)+z*py-y*pz
    amnuc2(2)=amnuc2(2)+x*pz-z*px
    amnuc2(3)=amnuc2(3)+y*px-x*py
12  continue
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!      DOUBLE PRECISION FUNCTION  RPOTEN(R)
!      REAL*8 R
!
!      RPOTEN=R
!      RETURN
!      END
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! COLDEV was here until removed by CMJ (XCP-3, LANL) to remove all
!    DEAD CODE.  COLDEV was NOT called by any part of LAQGSM (or GSM)
!    Date: 09/07/2017
! =====================================================================

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine helpq(r0n,anucl,a,c,d,tf0,rm)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    Determines nucleon density, Fermi momentum, Fermi energy
!
    dimension   w(8),fiks(8)
    data w/0.1012285363,0.2223810345,0.3137066459,0.3626837834, &
         & 0.3626837834,0.3137066459,0.2223810345,0.1012285363/, &
         & fiks/0.9602898565,0.7966664774,0.5255324099,0.1834346425, &
         & -0.1834346425,-0.5255324099,-0.7966664774,-0.9602898565/
    pi=3.141593
    if(anucl <= 10.) then
       a=r0n
       rm=a*sqrt(-log(d))
       ro0=anucl/((sqrt(pi)*a)**3)
    else
       s=0.
       a=r0n*anucl**0.333333333
       rm=a+c*log((1.-d)/d)
       do  k=1,8
          sk=(((fiks(k)+1.)**2)/(exp(rm*(fiks(k)+1.)/(2.*c))+ &
               & exp(a/c)))*w(k)
          s=s+sk
       enddo
       s=(rm**3)*s/8.
       ro0=anucl/(12.566370*(exp(a/c)+1.)*s)
    endif
    tf0=0.1985*(ro0/2.)**0.666666667
    rh=a
    if(anucl <= 10.)  rh=a*sqrt(log(2.d0))
!      write(16,100) ANUCL,RH,RO0,TF0
!      write( *,100) ANUCL,RH,RO0,TF0
!  100 FORMAT(1X,'ANUCL = ',F6.0,'  R(1/2)=  ',F6.3,
!     &'  RO0(NUCL)=',F6.4,'  TF0=',F6.4)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function cosexq (l,t,cm)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     BLOCK OF COSINUS CALCULATION FOR CHARGE EXCHANGE SCATTERING.
!  kkg 10/28/03        includes gamma +N --> pi+- + N1
!
    if(l.ne.0)              then
!  gamma + N => pi+-  + N1
       if(t <= 0.51)     then
!          COSEXQ = COSTAQ(24,T)
          cosexq = cosgaml(24)
       elseif(t <= 1.0)  then
!          COSEXQ = COSTAQ(25,T)
          cosexq = cosgaml(25)
       else
!          COSEXQ = COSTAQ(26,T)
          cosexq = cosgaml(26)
       endif
       return
    endif
11  if (t-0.08) 15,15,16
15  cosexq = costaq(12,t)
    return
16  if (t-0.3) 17,17,18
17  cosexq = costaq(13,t)
    return
18  if (t-1.0) 19,19,20
19  cosexq = costaq(10,t)
    return
20  if (t-2.4) 21,21,22
21  cosexq = costaq(11,t)
    return
22  tmax=3.5344*t*(t+cm)/(1.88*(t+cm)+cm*cm+0.8836)
    bet=bhn(t,cm)
    exm=bet*tmax
    if(exm > 30.)    exm=30.
    cosexq=1.+(2.*log(1.+rndm(-1.0_real64)*(exp(-exm)-1.)))/exm
    return
  end
!
! ********************************************************************
!
  double precision function cosgaml (j0)

! ======================================================================
!
!     Cosine calculation for elastic and
!     charge-exchange gamma+ N reactions using linear interpolation
!   energy of gamma is fixed in inigam
!   Called by: COSELQ COSEXQ
!
!   written by K. K. Gudima, Mar. 2004
! ======================================================================


    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit double precision (a-h,o-z), integer(int32) (i-n)
    common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)
    data degrad/0.0174532925199432958d0/
    if(j0 == 22.or.j0 == 23)  then
       jg=2
    elseif(j0 == 24.or.j0 == 25.or.j0 == 26)  then
       jg=1
    else
       write(*,*)  ' cosgaml: j0=',j0
       stop
    endif
    rr = rndm(-1.0_real64)
    do  ir=1,181
       if(abs(rr-ri(jg,ir)) < 0.00001)  then
          cth=cos(thetai(ir)*degrad)
          go  to  2
       endif
    enddo
    do  ir=2,181
       if(rr < ri(jg,ir))  then
          ir1=ir-1
          ir2=ir
          go  to  1
       endif
    enddo
    ir1=180
    ir2=181
1   continue
    x1=ri(jg,ir1)
    x2=ri(jg,ir2)
    y1=thetai(ir1)
    y2=thetai(ir2)
    th=y1+(y2-y1)*((rr-x1)/(x2-x1))
    cth=cos(th*degrad)
2   temp1 = abs(cth)
    if (temp1 <= 1.0) then
       cosgaml= cth
    else
       cosgaml= sign(1.d0, cth)
    endif
    return

! ======================================================================
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function bhn(t,cm)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!  Determines slope parameters of distribution exp(-b*t) for
!  high energy elastic scattering
!
    if(cm-0.9)  10,10,11
10  bhn=11.
    if(t < 10.)  bhn=7.5
    return
11  bhn=8.3+0.56*log(1.88*(t+1.88))
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function coselq (l,mq,ksi,t,cm)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     BLOCK OF COSINUS CALCULATION FOR ELASTIC SCATTERING.
!  kkg 10/28/03        includes gamma +N --> pi0 + N
!
13  if(mq-2) 24,14,78
!  N-He4 scattering
78  coselq=cosal(t)
    return
!  N-N scattering
14  go  to  (17,15,17,15,15,17),ksi
!   n + p:
15  if (t-0.97) 16,16,19
16  coselq = costaq(3,t)
    return
!   n + n or p + p:
17  if (t-0.46) 18,18,19
18  coselq=1.-2.*rndm(-1.0_real64)
    return
!   All types:
19  if (t-2.8) 20,20,21
20  coselq = (1.+costaq(1,t))/2.
    return
21  if (t-10.) 22,22,23
22  coselq = (3.+costaq(2,t))/4.
    return
23  tmax=3.5344*t*(t+cm)/(1.88*(t+cm)+cm*cm+0.8836)
    bet=bhn(t,cm)
    exm=bet*tmax
    if(exm > 30.)    exm=30.
    coselq=1.+(2.*log(1.+rndm(-1.0_real64)*(exp(-exm)-1.)))/exm
    return
24  if(l.ne.0)  go  to  42
!  pi-N scattering
    if (ksi-2) 25,34,33
!  pi+ p or pi- n scattering:
25  if (t-0.08) 26,26,27
26  coselq = costaq(4,t)
    return
27  if (t-0.3) 28,28,29
28  coselq = costaq(5,t)
    return
29  if (t-1.) 30,30,31
30  coselq = costaq(6,t)
    return
31  if (t-2.4) 32,32,23
32  coselq = costaq(7,t)
    return
!  pi0 p or pi0 n scattering:
33  if(rndm(-1.0_real64)-0.5) 25,25,34
!  pi+ n or pi- p scattering:
34  if (t-0.08) 35,35,36
35  coselq = costaq(8,t)
    return
36  if (t-0.3) 37,37,38
37  coselq = costaq(9,t)
    return
38  if (t-1.0) 39,39,40
39  coselq = costaq(10,t)
    return
40  if (t-2.4) 41,41,23
41  coselq = costaq(11,t)
    return
42  continue
!  gamma + N => pi0 + N
    if(t <= 0.45)  then
!        COSELQ = COSTAQ (22,T)
       coselq = cosgaml(22)
    else
!        COSELQ = COSTAQ (23,T)
       coselq = cosgaml(23)
    endif
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     * * * * * * * * * * * * * * * * *
  block data sigar

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     This subprogram puts the tables of cross sections
!     and energies into the appropriate common blocks.
!     EDITED by KKG, October, 2003, April 2007
!
    common /tabele/sigma,argus/typecsq/icst,nsicst
    dimension sigma(30,46) , argus(30,9) , icst(46) , nsicst(37) &
         & ,a1(30,3),a2(30,3),a3(30,3),a4(30,3), &
         & a5(30,3),a6(30,3),a7(30,3),a8(30,3),a9(30,2), &
         & b1(30,2),b2(30,2),b3(30,2),b4(30,2),ex(30,3), &
         & ar1(30,3),ar2(30,2),ar3(30,1), &
         & arg(30,3),ag1(30,2),ag2(30,3),ag3(30,3),ag4(30,1)
    equivalence (sigma(1, 1),a1(1,1)),(sigma(1, 4),a2(1,1)), &
         & (sigma(1, 7),a3(1,1)),(sigma(1,10),a4(1,1)), &
         & (sigma(1,13),a5(1,1)),(sigma(1,16),a6(1,1)), &
         & (sigma(1,19),a7(1,1)),(sigma(1,22),a8(1,1)), &
         & (sigma(1,25),a9(1,1)),(sigma(1,27),b1(1,1)), &
         & (sigma(1,29),b2(1,1)),(sigma(1,31),b3(1,1)), &
         & (sigma(1,33),b4(1,1)),(sigma(1,35),ex(1,1)), &
         & (argus(1, 1),ar1(1,1)),(argus(1,4),ar2(1,1)), &
         & (argus(1, 6),ar3(1,1)), &
         & (argus(1,7),arg(1,1)),(sigma(1,38),ag1(1,1)), &
         & (sigma(1,40),ag2(1,1)),(sigma(1,43),ag3(1,1)), &
         & (sigma(1,46),ag4(1,1))
    data  a1/ &
!   (1) N+N , P+P total
         & 17613.0,   302.911,   150.083,    95.861,    69.160, &
         & 54.056,    38.395,    29.639,    23.563,    22.800, &
         & 22.600,    23.100,    24.045,    25.957,    30.747, &
         & 42.160,    48.114,    47.430,    47.465,    47.307, &
         & 46.681,    45.082,    42.476,    41.236,    40.866, &
         & 40.143,    40.184,    38.352,    38.375,    38.406, &
!   (2) N+N,  P+P elastic
         & 17613.0,   302.911,   150.083,    95.861,    69.160, &
         & 54.056,    38.395,    29.639,    23.563,    22.800, &
         & 22.600,    23.097,    23.399,    23.824,    24.623, &
         & 25.060,    24.306,    23.388,    22.543,    21.101, &
         & 19.701,    17.264,    14.477,    12.848,    11.745, &
         & 10.691,     9.640,     8.583,     8.296,     7.956, &
!   (3) N+P  total
         & 20357.0,   912.566,   499.490,   288.179,   208.121, &
         & 161.621,   106.227,    71.137,    50.588,    42.483, &
         & 37.982,    35.381,    33.820,    33.299,    34.219, &
         & 37.566,    38.801,    39.255,    40.368,    41.570, &
         & 41.624,    42.144,    42.584,    42.091,    41.847, &
         & 41.26029,  39.48050,  39.40114,  39.37191,  39.33528/
    data  a2/ &
!   (4) N+P  elastic
         & 20357.0,   912.566,   499.490,   288.179,   208.121, &
         & 161.621,   106.227,    71.137,    50.588,    42.483, &
         & 37.982,    35.375,    33.000,    30.700,    27.800, &
         & 24.700,    21.100,    20.000,    18.750,    17.700, &
         & 17.000,    15.750,    14.200,    13.300,    12.700, &
         & 11.950,    11.100,     9.97127,   9.60346,   9.16243, &
!   (5) pi+ + N, pi- + P total
         & 5.200,     6.400,     9.300,    17.000,    37.150, &
         & 58.940,    72.000,    59.050,    39.880,    26.460, &
         & 27.430,    32.010,    46.730,    35.900,    45.000, &
         & 59.550,    43.570,    36.440,    36.430,    36.530, &
         & 34.780,    34.370,    36.020,    36.300,    34.600, &
         & 32.950,    32.200,    28.900,    26.000,    25.47546, &
!   (6) pi+ + N, pi- + P elastic
         & 1.593,     1.784,     2.281,     4.117,    10.467, &
         & 19.657,    23.245,    23.589,    15.349,     9.931, &
         & 10.087,    14.277,    20.290,    14.807,    18.139, &
         & 25.873,    18.879,    14.746,    13.019,    12.345, &
         & 9.897,    10.019,     9.553,     8.942,     8.448, &
         & 8.034,     7.188,     5.325,     4.258,     3.991/
    data  a3/ &
!   (7) pi(+0-) + (N,P)  single charge exchange
         & 3.607,     5.303,     8.281,    13.170,    26.134, &
         & 42.727,    48.927,    38.530,    25.550,    14.713, &
         & 10.465,     8.442,     7.448,     4.610,     6.693, &
         & 6.584,     2.828,     2.465,     2.438,     2.360, &
         & 2.026,     1.459,     1.052,     0.792,     0.610, &
         & 0.380,     0.255,     0.0903,    0.0381,    0.0252, &
!   (8) pi- + N,  pi+ + P total
         & 2.015,     5.866,    12.637,    32.909,    97.023, &
         & 172.340,   209.160,   167.843,   103.625,    52.363, &
         & 30.821,    18.862,    14.257,    15.979,    22.619, &
         & 26.113,    26.892,    31.413,    37.517,    41.354, &
         & 37.945,    31.563,    29.571,    29.083,    30.211, &
         & 31.127,    28.865,    26.390,    23.357,    23.348, &
!   (9) pi- + N,  pi+ + P elastic
         & 2.015,     5.866,    12.637,    32.909,    97.023, &
         & 172.340,   209.160,   167.843,   103.625,    52.363, &
         & 30.566,    17.281,    11.473,     8.214,     9.272, &
         & 12.004,    13.186,    14.052,    15.522,    18.586, &
         & 15.715,    12.785,    10.808,     9.548,     8.649, &
         & 7.936,     7.051,     5.668,     4.216,     3.951/
    data  a4/ &
!   (10) pion + np absorption cross section
!   (Must be multiplied by a function [currently = 4.0] for use)
         & 3.900,     4.350,     6.090,     8.540,    11.050, &
         & 11.750,    10.190,     6.540,     3.400,     1.700, &
         & 0.825,     0.380,     0.200,     0.130,     0.089, &
         & 0.063,     0.049,     0.040,     0.031,     0.026, &
         & 0.018,     0.0125,    0.0076,    0.0040,    0.0021, &
         & 0.00085,   0.00030,   0.00002,   0.000003,  0.0000003, &
!   (11) p + p --> p + p + pi0
         & 0.000,     0.000,     0.00295,   0.01969,   0.07956, &
         & 0.24419,   0.5399,    1.086,     1.953,     3.029, &
         & 3.522,     3.730,     3.845,     4.181,     4.526, &
         & 4.619,     4.571,     4.483,     4.282,     4.164, &
         & 4.099,     4.023,     3.960,     3.800,     3.374, &
         & 3.089,     2.603,     2.889,     1.532,     1.320, &
!   (12) p + p --> p + n + pi+
         & 0.000,     0.000,     0.000,     0.35747,   0.79723, &
         & 1.582,     2.053,     4.773,     7.542,    10.113, &
         & 13.214,    15.564,    16.767,    17.219,    17.589, &
         & 17.851,    18.190,    18.279,    18.330,    17.926, &
         & 17.546,    17.068,    16.420,    16.067,    13.130, &
         & 10.534,     8.688,     7.278,     4.480,     3.686/
    data  a5/ &
!  (13) p + n --> p + n + pi0
         & 0.000,     0.000,     0.05749,   0.28462,   0.5271, &
         & 0.9429,    1.8173,    3.142,     4.761,     6.508, &
         & 6.778,     7.169,     7.400,     7.550,     7.700, &
         & 7.806,     7.982,     8.195,     8.424,     8.640, &
         & 8.920,     9.201,     9.055,     8.600,     7.700, &
         & 6.550,     5.129,     3.985,     2.448,     1.300, &
!  (14) p + n --> p + p + pi-
         & 0.000,     0.000,     0.000,     0.000,     0.06818, &
         & 0.15784,   0.29075,   0.6300,    1.0910,    1.3824, &
         & 1.6719,    1.9102,    2.1401,    2.350,     2.526, &
         & 2.675,     2.806,     2.918,     3.009,     3.097, &
         & 3.168,     3.227,     3.199,     3.121,     2.875, &
         & 2.51422,   1.85326,   1.46821,   1.22737,   1.07373, &
!  (15) pi+ + p --> pi+ + p + pi0
         & 0.00227,   0.02183,   0.08394,   0.22867,   0.38106, &
         & 0.63783,   1.0998,    2.0452,    3.1995,    4.574, &
         & 6.276,     9.381,    10.310,    10.343,     9.484, &
         & 8.784,     8.277,    10.559,    12.108,    11.563, &
         & 10.076,     7.845,     6.545,     5.452,     3.589, &
         & 2.6616,    1.77222,   1.15771,   0.79248,   0.60762/
    data  a6/ &
!  (16) pi+ + p --> pi+ + n + pi+
         & 0.00147,   0.03056,   0.05862,   0.09654,   0.13958, &
         & 0.21361,   0.34633,   0.51922,   0.68048,   0.93144, &
         & 1.1971,    1.480,     1.781,     1.977,     2.083, &
         & 2.097,     2.117,     2.499,     3.107,     3.474, &
         & 3.648,     3.631,     3.099,     2.309,     2.467, &
         & 2.102,     1.2741,    0.85565,   0.70307,   0.53185, &
!  (17) pi- + p --> pi0 + n + pi0
         & 0.01433,   0.21686,   0.61527,   1.25405,   1.52076, &
         & 1.7479,    1.9039,    2.0563,    2.1526,    2.250, &
         & 2.402,     2.807,     3.294,     3.623,     3.085, &
         & 2.251,     1.841,     1.445,     1.183,     0.9801, &
         & 0.80905,   0.59528,   0.62035,   0.70180,   0.40717, &
         & 0.25280,   0.12009,   0.08365,   0.015573,  0.013499, &
!  (18) pi- + p --> pi- + p + pi0
         & 0.00193,   0.02508,   0.0919,    0.2654,    0.6346, &
         & 1.180,     2.020,     3.237,     5.008,     4.861, &
         & 4.529,     5.011,     5.687,     6.188,     6.426, &
         & 5.896,     5.286,     4.383,     4.632,     4.700, &
         & 4.737,     5.115,     4.128,     3.363,     3.395, &
         & 1.876,     2.0066,    1.1843,    0.68066,   0.60043/
    data  a7/ &
!  (19) pi- + p --> pi- + n + pi+
         & 0.00908,   0.12116,   0.64364,   1.5597,    3.153, &
         & 3.971,     5.004,     5.860,     6.306,     6.474, &
         & 7.194,     8.300,     9.978,    11.775,    11.789, &
         & 9.938,     7.522,     7.512,     7.930,     7.781, &
         & 7.089,     6.944,     6.127,     5.332,     4.000, &
         & 3.385,     2.967,     2.10879,   1.20557,   0.87941, &
         & 310.00,   270.00,   240.00,   213.00,   172.00, &
         & 156.00,   132.00,   116.00,   108.00,   102.00, &
         & 98.00,    96.00,    96.00,    98.00,   102.00, &
         & 109.00,   116.00,   122.00,   128.00,   132.00, &
         & 136.00,   138.00,   140.00,   140.00,   140.00, &
         & 141.00,   140.00,   138.00,   134.00,   132.00, &
         & 204.00,   165.00,   138.00,   117.00,    87.00, &
         & 75.00,    55.00,    40.00,    30.00,    24.00, &
         & 22.00,    20.00,    20.00,    20.00,    21.00, &
         & 23.00,    26.00,    28.00,    31.00,    33.00, &
         & 34.00,    36.00,    36.00,    36.00,    36.00, &
         & 36.00,     36.00,     36.00,     34.00,     32.00/
    data  a8/ &
         & 266.00,   228.00,   200.00,   178.00,   148.00, &
         & 136.00,   115.00,   102.00,    94.00,    86.00, &
         & 79.00,    76.00,    76.00,    78.00,    81.00, &
         & 85.00,    91.00,   100.00,   107.00,   112.00, &
         & 115.00,   118.00,   118.00,   118.00,   118.00, &
         & 118.00,   116.00,   112.00,   110.00,   106.00, &
         & 133.00,   105.00,    84.00,    68.00,    46.00, &
         & 38.00,    28.00,    22.00,   19.00,    15.00, &
         & 11.00,   10.00,     9.50,    10.00,   11.00, &
         & 12.00,    15.00,    17.00,    19.00,    21.00, &
         & 23.00,    24.00,    24.00,    24.00,    23.00, &
         & 22.00,    22.00,    21.00,   20.00,    20.00, &
         & 0.00,     0.00,     0.00,     0.20,     0.80, &
         & 1.60,     3.20,     6.20,    33.00,    40.30, &
         & 41.60,    42.00,    42.00,    42.00,    41.90, &
         & 41.50,    40.70,    40.00,    39.30,    38.40, &
         & 37.60,    36.60,    35.40,    34.20,    32.80, &
         & 31.20,     29.40,     27.70,     26.20,     16.00/
    data  a9/ &
         & 000000.00,     0.00,     0.00,     0.30,     1.20, &
         & 2.20,     3.40,     4.90,     9.00,    11.20, &
         & 13.40,    15.60,    18.10,    20.30,    22.30, &
         & 25.40,    27.10,    27.90,    28.30,    28.60, &
         & 28.80,    29.00,    29.10,    29.20,    29.30, &
         & 29.30,    29.30,    29.30,    29.30,    29.30, &
         & 0.00,     1.00,    30.00,    41.00,    43.20, &
         & 43.60,    43.40,    43.00,    41.00,    39.00, &
         & 32.00,    23.50,    16.00,    12.00,     9.40, &
         & 5.60,     3.80,     2.90,     2.20,    1.80, &
         & 1.60,     1.30,     1.00,     0.80,     0.70, &
         & 0.50,      0.40,      0.40,      0.30,      0.00/
    data  b1/ &
         & 10.0, 11.0, 11.5, 11.9, 12.0, 12.1, 12.2, 12.2, 12.8, 14.5, &
         & 16.3, 18.9, 18.0, 17.6, 17.2, 17.1, 17.0, 17.0, 17.1, 17.2, &
         & 17.5, 17.8, 18.0, 19.0, 19.4, 20.0, 20.7, 21.4, 22.0, 24.0, &
         & 10.0, 11.0, 11.5, 11.9, 12.0, 12.1, 12.2, 12.1, 12.0, 11.9, &
         & 11.6, 10.5,  9.0,  6.5,  5.0,  4.1,  3.5,  3.2,  3.2,  3.1, &
         & 2.9,  2.5,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4/
    data  b2/ &
         & 6.0,  7.0,  7.5, 10.0, 12.0, 13.0, 13.0, 13.5, 15.0, 16.5, &
         & 18.5, 20.5, 19.0, 18.5, 17.8, 17.3, 17.2, 17.2, 17.2, 17.5, &
         & 17.7, 18.0, 18.2, 18.5, 19.0, 19.7, 20.0, 20.3, 20.6, 21.0, &
         & 3.0,  4.0,  4.5,  6.0,  8.0,  5.5,  5.5,  5.7,  5.7,  5.7, &
         & 5.7,  5.5,  5.2,  5.0,  4.8,  4.5,  4.2,  3.8,  3.5,  3.1, &
         & 2.3,  2.5,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4/
    data  b3/ &
         & 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 20.0, 26.0, 29.0, 31.0, &
         & 31.0, 29.0, 28.0, 23.0, 22.0, 21.0, 20.7, 20.5, 20.3, 20.2, &
         & 20.0, 20.0, 20.0, 20.0, 20.1, 20.5, 21.0, 22.0, 23.0, 24.0, &
         & 5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  7.0,  9.0, 15.0, 18.0, &
         & 17.0, 10.0,  7.5,  5.0,  3.4,  3.1,  3.0,  2.9,  2.8,  2.7, &
         & 2.65, 2.60, 2.50, 2.51, 2.52, 2.55, 2.56, 2.57, 2.58, 2.58/
    data  b4/ &
         & 450.0,323.0,125.0, 82.0, 65.0, 38.0,30.75, 29.0, 32.5, 34.0, &
         & 50.0, 32.0, 33.0, 30.0, 28.0, 27.0, 25.0, 24.0, 23.0, 22.0, &
         & 21.0, 20.5, 20.0, 20.0, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, &
         & 150.0, 98.0, 60.0, 42.0, 33.0, 23.0, 18.0, 16.0, 17.0, 20.0, &
         & 22.0, 15.0,  9.0,  8.0,  5.5,  4.5,  4.0,  3.6,  3.2,  2.9, &
         & 2.8,  2.6,  2.5, 2.51, 2.52, 2.55, 2.56, 2.57, 2.58, 2.58/
    data  ex/ &
         & 3.0,  3.0,  3.0,  4.0,  4.0,  7.5,  7.5,  7.8,  7.8,  7.8, &
         & 7.8,  4.5,  3.0,  1.8,  0.8, 0.45, 0.30, 0.15,0.075,0.035, &
         & 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, &
         & 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  9.0,  8.0,  7.0, &
         & 6.5,  4.0, 2.75, 1.50, 0.65, 0.30,  0.0,  0.0,  0.0,  0.0, &
         & 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, &
         & 300.0,225.0, 65.0, 40.0, 32.0, 15.0,12.75, 12.0, 15.5, 14.0, &
         & 9.75,  5.0,  4.0,  2.0, 1.25,  0.5,  0.0,  0.0,  0.0,  0.0, &
         & 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
!
!
    data ag1/ &
!   (38) gamma + p --> p + pi0
!   New constants (August, 1998) (Valid to 10 GeV [use log for low E])
         & 0.00021,   0.00437,   0.02075,   0.09167,   0.2423, &
         & 0.2796,    0.1371,    0.08599,   0.05499,   0.04352, &
         & 0.03260,   0.03310,   0.03990,   0.04336,   0.03950, &
         & 0.03307,   0.02750,   0.02847,   0.02950,   0.02337, &
         & 0.01960,   0.01637,   0.01320,   0.00863,   0.00574, &
         & 0.00346,   0.00200,   0.00076,   0.000212,  0.000065, &
!   (39) gamma + p --> n + pi+:
!   New constants (August, 1998) (Valid to 10 GeV [use log for low E])
         & 0.000001,  0.08169,   0.10202,   0.15975,   0.2343, &
         & 0.2031,    0.1230,    0.10029,   0.09024,   0.08760, &
         & 0.08789,   0.09241,   0.10347,   0.09100,   0.06329, &
         & 0.04951,   0.04944,   0.05184,   0.05612,   0.05019, &
         & 0.03622,   0.02501,   0.01922,   0.01350,   0.00774, &
         & 0.00516,   0.00250,   0.00075,   0.000232,  0.000083/
    data ag2/ &
!   (40) gamma + [NN] absorption:
!   New constants (August, 1998) (Valid to 5 GeV)
!   (Must be multiplied by a function [currently = 5.0] for use)
         & 0.000,     1.14301,   1.8797,    2.5348,    1.3234, &
         & 0.57795,   0.3391,    0.2289,    0.1259,    0.08875, &
         & 0.06946,   0.06128,   0.05130,   0.05410,   0.05460, &
         & 0.05805,   0.06317,   0.05787,   0.04589,   0.03161, &
         & 0.02021,   0.01355,   0.01090,   0.00830,   0.00600, &
         & 0.00450,   0.00340,   0.00270,   0.00100,   0.00050, &
!   (41) gamma + p --> pi- + p + pi+
!   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
         & 1.82e-6,   0.000163,  0.00312,   0.00762,   0.02099, &
         & 0.04421,   0.06208,   0.07277,   0.07452,   0.07485, &
         & 0.07468,   0.07448,   0.07434,   0.07355,   0.07194, &
         & 0.07058,   0.06931,   0.06781,   0.06471,   0.06091, &
         & 0.05760,   0.05302,   0.04648,   0.03943,   0.03271, &
         & 0.02806,   0.02357,   0.01731,   0.01368,   0.01219, &
!   (42) gamma + p --> pi0 + p + pi0
!   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
!   These coefficients are found by multiplying those of # 41 by
!   0.142, found by optimally fitting the limited available data
!   in the energy range 464 to 779 MeV gamma lab energy to the #23
!   cross section linearly scaled.
         & 2.58e-7,   0.0000231, 0.000443,  0.00108,   0.00298, &
         & 0.00628,   0.00882,   0.01033,   0.01058,   0.01063, &
         & 0.01060,   0.01058,   0.01056,   0.01044,   0.01022, &
         & 0.01022,   0.00984,   0.00963,   0.00919,   0.00865, &
         & 0.00818,   0.00753,   0.00660,   0.00560,   0.00464, &
         & 0.00398,   0.00335,   0.00246,   0.00194,   0.00173/
    data ag3/ &
!   (43) gamma + p --> pi+ + n + pi0
!   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
!   These coefficients are found by fitting the 12 data points from
!   464 to 779 MeV gamma lab energy, and fixing the coefficients
!   from 800 to 20000 MeV to 2/3 * (#41 coefficients).  There is
!   an unphysical bump (with negligible cross section) at energy below
!   400 MeV due to limitations in the parabolic interpolation method.
         & 2.1e-9,    1.1e-6,    0.00061,   0.00124,   0.00610, &
         & 0.01491,   0.02359,   0.03218,   0.04855,   0.05341, &
         & 0.04979,   0.04965,   0.04956,   0.04903,   0.04796, &
         & 0.04705,   0.04621,   0.04520,   0.04314,   0.04061, &
         & 0.03849,   0.03532,   0.03099,   0.02629,   0.02181, &
         & 0.01871,   0.01570,   0.011557,  0.009120,  0.008126, &
!   (44) gamma + p --> delta + pi
!   New constants (August, 1998) (Valid to 10 GeV [use log for low E])
!   There is an unphysical zero at energy below 400 MeV due to
!   limitations in the parabolic interpolation method (the cross
!   section is negligible anyway).
         & 1.8e-11,   0.00002,   0.00082,   0.00762,   0.02099, &
         & 0.04421,   0.06208,   0.06935,   0.06387,   0.05537, &
         & 0.04570,   0.04594,   0.05130,   0.05416,   0.05146, &
         & 0.04503,   0.03727,   0.02564,   0.02613,   0.01998, &
         & 0.01725,   0.01511,   0.01071,   0.00634,   0.00414, &
         & 0.00289,   0.00174,   0.000821,  0.000257,  0.000124, &
!   (45)  gamma + p --> 2 pi + N (#41 + #42 + #43)
!   New constants (August, 1998) (Valid to 20 GeV [use log for low E])
         & 2.1e-6,    0.000179,  0.00417,   0.00994,   0.03007, &
         & 0.06540,   0.09449,   0.11528,   0.13365,   0.13889, &
         & 0.13507,   0.13471,   0.13446,   0.13302,   0.13012, &
         & 0.12785,   0.12536,   0.12264,   0.11704,   0.11017, &
         & 0.10427,   0.09587,   0.08407,   0.07132,   0.05916, &
         & 0.05075,   0.04262,   0.03133,   0.02474,   0.02205/
    data ag4/ &
!   (46)  gamma + p --> X (TOTAL)
!   New constants (August, 1998) (Valid to 100 GeV [use log for low E])
         & 0.00010,   0.12102,   0.24837,   0.4994,    0.4709, &
         & 0.2727,    0.1877,    0.1699,    0.1880,    0.2146, &
         & 0.2339,    0.2659,    0.2772,    0.2316,    0.2093, &
         & 0.2102,    0.2013,    0.2142,    0.2120,    0.1842, &
         & 0.1693,    0.1539,    0.1509,    0.1524,    0.1446, &
         & 0.1394,    0.1283,    0.1289,    0.11739,   0.11714/
    data  icst/ &
         & 210 ,   211 ,   220 ,   221 ,   120 ,   121 ,   122 , &
         & 110 ,   111 ,   123 ,   214 ,   215 ,   224 ,   225 , &
         & 114 ,   115 ,   126 ,   124 ,   125 , &
         & 00510 , 00511 , 00410 , 00411 , 00514 , 00515 , 00516, &
         & 1110 ,  1111 ,  1120 ,  1121 , -1110 , -1111 , -1120, &
         & -1121 ,  1122 , -1112 , -1122 , &
         & 10111,  10112,  10113,  10115,  10114,  10116,  10118, &
         & 10117,  10110/
    data  nsicst/ &
         & 112 ,   113 ,   116 ,   117 ,   127 ,   130 ,   131 , &
         & 132 ,   133 ,   134 ,   135 ,   136 ,   137 ,   212 , &
         & 213 ,   216 ,   217 ,   222 ,   223 ,   226 ,   227 , &
         & 230 ,   231 ,   232 ,   233 , &
         & 240 ,   241 ,   242 ,   243 , &
         & 250 ,   251 ,   252 ,   253 , &
         & 260 ,   261 ,   262 ,   263 /
    data  ar1/ &
!   Lab kinetic energies for nucleon-nucleon reactions (#1-#4):
         & 0.000,     0.010,     0.020,     0.030,     0.040, &
         & 0.050,     0.070,     0.100,     0.150,     0.200, &
         & 0.250,     0.300,     0.350,     0.400,     0.500, &
         & 0.650,     0.850,     0.950,     1.100,     1.300, &
         & 1.500,     2.000,     3.000,     4.000,     5.000, &
         & 7.000,    10.000,    16.000,    22.000,    30.000, &
!  Lab kinetic energies for pion induced reactions(2)  (#5-#10):
         & 0.000,     0.030,     0.050,     0.080,     0.120, &
         & 0.150,     0.175,     0.210,     0.250,     0.320, &
         & 0.400,     0.500,     0.600,     0.700,     0.800, &
         & 0.900,     1.000,     1.100,     1.200,     1.300, &
         & 1.450,     1.600,     1.800,     2.000,     2.250, &
         & 2.500,     3.000,     5.000,    12.000,    20.000, &
!   Lab kinetic energies for pion production reactions (3) (#11-19):
         & 0.200,     0.250,     0.300,     0.350,     0.400, &
         & 0.450,     0.500,     0.550,     0.600,     0.650, &
         & 0.700,     0.750,     0.800,     0.850,     0.900, &
         & 0.950,     1.000,     1.100,     1.200,     1.300, &
         & 1.400,     1.600,     1.800,     2.000,     2.500, &
         & 3.000,     4.000,     5.000,     7.000,    10.000/
    data  ar2/ &
!   Lab kinetic energies for N + He-4 and N + He-3(H-3) reactions
         & 0.05  ,   0.06  ,   0.07  ,   0.08  ,   0.10  , &
         & 0.11  ,   0.13  ,   0.15  ,   0.17  ,   0.20  , &
         & 0.25  ,   0.30  ,   0.35  ,   0.40  ,   0.45  , &
         & 0.50  ,   0.55  ,   0.60  ,   0.65  ,   0.70  , &
         & 0.75  ,   0.80  ,   0.85  ,   0.90  ,    1.00  , &
         & 2.00  ,   3.00  ,   5.00  ,   10.0  ,   30.0  , &
!   Lab kinetic energies for N + He-4 ==> N + N + (A=3)
         & 0.0185,   0.0200,   0.0224,   0.0250,   0.0282, &
         & 0.0316,   0.0355,   0.0400,   0.0500,   0.0560, &
         & 0.0630,   0.0710,   0.0800,   0.0900,   0.1000, &
         & 0.1260,   0.1570,   0.2000,   0.2500,   0.3160, &
         & 0.4000,   0.5000,   0.6300,   0.8000,   1.0000, &
         & 1.2600,   1.5700,   2.0000,   2.5000,   30.000/
    data  ar3/ &
!   Lab kinetic energies for K + N reactions
         & 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, &
         & 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0,15.0, &
         & 20.0,30.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0, &
         & 500.0, 1000.0/
    data arg/ &
!   Lab Energies for gamma-induced pion production cross sections
!   (new version, #'s 38,39,46)
         & 0.150,     0.165,     0.200,     0.250,     0.300, &
         & 0.350,     0.400,     0.450,     0.500,     0.550, &
         & 0.600,     0.650,     0.700,     0.750,     0.800, &
         & 0.850,     0.900,     0.950,     1.000,     1.050, &
         & 1.100,     1.150,     1.200,     1.350,     1.500, &
         & 1.750,     2.250,     3.500,     6.000,    10.000, &
!   Lab energies for gamma absorption on 2 nucleons (5) (#40):
         & 0.002226,  0.0025,    0.003,     0.005,     0.010, &
         & 0.020,     0.030,     0.040,     0.060,     0.080, &
         & 0.100,     0.120,     0.150,     0.170,     0.200, &
         & 0.220,     0.250,     0.300,     0.330,     0.360, &
         & 0.400,     0.440,     0.500,     0.600,     0.700, &
         & 0.800,     0.900,     1.000,     2.000,     3.000, &
!   Lab energies for gamma induced 2-pion production (6) (#41-45):
         & 0.320,     0.350,     0.400,     0.450,     0.500, &
         & 0.550,     0.600,     0.650,     0.700,     0.750, &
         & 0.800,     0.850,     0.900,     0.950,     1.000, &
         & 1.050,     1.100,     1.150,     1.200,     1.250, &
         & 1.300,     1.500,     1.700,     2.000,     2.500, &
         & 3.000,     4.000,     6.000,    10.000,    14.000/
!
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  block data coefabq

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     SUBROUTINE WHICH PUT THE COEFFICIENTS ANKJ,BNKJ AND CKJ
!     IN THE MAIN PROGRAM
!  kkg  10/28/03
    common /coefaq/ ankj  /coefbcq/ bnkj,ckj
    dimension ankj(4,4,29),bnkj(4,4,8),ckj(3,8) &
         & ,a1(4,4,4),a2(4,4,4),a3(4,4,4),a4(4,4,4), &
         & a5(4,4,4),a6(4,4,1),a7(4,4,7),a8(4,4,1), &
         & b1(4,4,4),b2(4,4,4)
    equivalence (ankj(1,1, 1),a1(1,1,1)), &
         & (ankj(1,1, 5),a2(1,1,1)), &
         & (ankj(1,1, 9),a3(1,1,1)), &
         & (ankj(1,1,13),a4(1,1,1)), &
         & (ankj(1,1,17),a5(1,1,1)), &
         & (ankj(1,1,21),a6(1,1,1)), &
         & (ankj(1,1,22),a7(1,1,1)), &
         & (ankj(1,1,29),a8(1,1,1)), &
         & (bnkj(1,1, 1),b1(1,1,1)), &
         & (bnkj(1,1, 5),b2(1,1,1))
    data  a1/ &
!  ankj Angular distribution coefficients:
!  j = 1;  N + N elastic scattering; Tlab <= 2.8 GeV:
!          (n + n & p + p isotropic below 0.46 GeV.)
         & 2.7404e+00 , -9.6998e+00 ,  1.0400e+01 ,  2.3882e+00 , &
         & -7.5137e+00 ,  4.4096e+01 , -7.4379e+01 ,  4.6038e+01 , &
         & 7.5479e+00 , -3.9274e+01 ,  6.4835e+01 , -4.1609e+01 , &
         & -1.8369e+00 ,  8.6911e+00 , -1.3060e+01 ,  7.1880e+00 , &
!  j = 2;  N + N elastic scattering; 2.8 < Tlab <= 10. GeV:
         & -3.0853e+01 ,  1.0624e+02 , -1.2939e+02 ,  5.4339e+01 , &
         & 1.9465e+01 , -6.8102e+01 ,  9.6358e+01 , -5.6827e+01 , &
         & -3.4831e+00 ,  1.2341e+01 , -1.8592e+01 ,  1.2024e+01 , &
         & 1.8941e-01 , -6.7880e-01 ,  1.0665e+00 , -7.2910e-01 , &
!  j = 3;  n + p elastic scattering; Tlab <= 0.97 GeV:
         & 1.0258e-01 , -1.0542e+00 ,  1.1389e+01 , -1.6638e+01 , &
         & -4.9607e-01 ,  1.1800e+01 , -9.0857e+01 ,  1.6476e+02 , &
         & 1.5437e+00 , -3.3769e+01 ,  2.5192e+02 , -4.5071e+02 , &
         & -1.2021e+00 ,  2.5336e+01 , -1.8658e+02 ,  3.3254e+02 , &
!  j = 4; pi+ p or pi- n elastic scattering; Tlab <= 0.080 GeV:
         & 1.5789e-01 ,  2.9671e+00 , -5.5251e+00 ,  6.8925e+00 , &
         & -7.0218e+00 , -2.0534e+02 ,  5.6951e+02 , -8.9858e+02 , &
         & 1.3496e+02 ,  4.8722e+03 , -1.4674e+04 ,  2.3924e+04 , &
         & -8.2116e+02 , -3.2586e+04 ,  1.0098e+05 , -1.6553e+05 /
    data  a2/ &
!  j = 5; pi+ p or pi- n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
         & 3.1531e-01 , -7.4981e+00 ,  4.3295e+01 , -7.6360e+01 , &
         & -6.5373e+00 ,  1.9307e+02 , -1.0181e+03 ,  1.7426e+03 , &
         & 4.6864e+01 , -1.3030e+03 ,  6.7291e+03 , -1.1075e+04 , &
         & -9.5192e+01 ,  2.6373e+03 , -1.2857e+04 ,  2.0294e+04 , &
!  j = 6; pi+ p or pi- n elastic scattering; 0.30 < Tlab <= 1.0 GeV:
         & -1.7953e+01 ,  1.0972e+02 , -2.3954e+02 ,  2.2826e+02 , &
         & 9.1968e+01 , -5.1963e+02 ,  1.1266e+03 , -1.0740e+03 , &
         & -1.3270e+02 ,  7.4112e+02 , -1.6000e+03 ,  1.5249e+03 , &
         & 5.8598e+01 , -3.1874e+02 ,  6.7751e+02 , -6.4011e+02 , &
!  j = 7; pi+ p or pi- n elastic scattering; 1.0 < Tlab <= 2.4 GeV:
         & 4.2169e-01 ,  1.4705e+02 , -6.5335e+02 ,  9.1507e+02 , &
         & -3.5198e+00 , -2.6019e+02 ,  1.2250e+03 , -1.7481e+03 , &
         & 3.6373e+00 ,  1.5592e+02 , -7.5201e+02 ,  1.0796e+03 , &
         & -7.8041e-01 , -3.0563e+01 ,  1.4795e+02 , -2.1250e+02 , &
!  j = 8; pi+ n or pi- p elastic scattering; Tlab <= 0.080 GeV:
         & -3.8288e-01 ,  3.7587e+00 , -6.5144e+00 ,  6.7740e+00 , &
         & 1.0381e+02 , -2.7282e+02 ,  4.7759e+02 , -5.1222e+02 , &
         & -1.7882e+03 ,  4.3052e+03 , -7.9314e+03 ,  9.3471e+03 , &
         & 7.1475e+03 , -3.3395e+03 , -4.1392e+03 , -4.4364e+03 /
    data  a3/ &
!  j = 9; pi- p or pi+ n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
         & 2.4991e-01 ,  3.2028e+01 , -1.1882e+02 ,  1.5099e+02 , &
         & -2.6994e+00 , -4.6045e+02 ,  1.8959e+03 , -2.5190e+03 , &
         & 1.6268e+01 ,  2.1384e+03 , -9.1262e+03 ,  1.2431e+04 , &
         & -2.9654e+01 , -3.1823e+03 ,  1.3944e+04 , -1.9342e+04 , &
!  j = 10; pi- p or pi+ n elastic or CX scattering;
!          0.30 < Tlab <= 1.0 GeV:
         & 3.9025e+00 , -9.1126e+01 ,  3.2373e+02 , -4.0048e+02 , &
         & -2.0619e+01 ,  4.9170e+02 , -1.7155e+03 ,  2.1143e+03 , &
         & 3.3004e+01 , -7.6684e+02 ,  2.7003e+03 , -3.3525e+03 , &
         & -1.6367e+01 ,  3.7394e+02 , -1.3202e+03 ,  1.6423e+03 , &
!  j = 11; pi- p or pi+ n elastic or CX scattering;
!          1.0 < Tlab <= 2.4 GeV:
         & 1.9402e+01 , -2.2446e+02 ,  7.4733e+02 , -9.3570e+02 , &
         & -4.4180e+01 ,  4.7194e+02 , -1.4856e+03 ,  1.8055e+03 , &
         & 3.1567e+01 , -3.0176e+02 ,  9.0763e+02 , -1.0773e+03 , &
         & -6.8648e+00 ,  6.0476e+01 , -1.7520e+02 ,  2.0381e+02 , &
!  j = 12; pi- + p --> pi0 n of pi+ + n --> pi0 + p Charge exchange
!          scattering; Tlab <= 0.08 GeV:
         & 1.4988e-01 ,  2.8753e+00 , -5.3078e+00 ,  6.2233e+00 , &
         & -5.9558e+00 , -1.6203e+02 ,  4.3079e+02 , -6.2548e+02 , &
         & 1.2875e+02 ,  3.1402e+03 , -7.9189e+03 ,  1.0983e+04 , &
         & -8.5161e+02 , -1.8780e+04 ,  4.4607e+04 , -5.8790e+04 /
    data  a4/ &
!  j = 13; pi- + p --> pi0 n of pi+ + n --> pi0 + p Charge exchange
!          scattering; 0.08 < Tlab <= 0.30 GeV:
         & 5.3689e-01 , -1.3216e+01 ,  8.1011e+01 , -1.4285e+02 , &
         & -1.0550e+01 ,  2.9629e+02 , -1.6957e+03 ,  2.8935e+03 , &
         & 6.9621e+01 , -1.9245e+03 ,  1.0620e+04 , -1.7468e+04 , &
         & -1.3865e+02 ,  3.9281e+03 , -2.0293e+04 ,  3.2058e+04 , &
!  j = 14;  N + N --> N + N + pi; nucleon distributions:
         & 8.5591e-02 ,  5.0390e+00 , -1.3782e+01 ,  1.4661e+01 , &
         & 5.4284e-02 , -9.2324e+00 ,  3.6397e+01 , -4.2962e+01 , &
         & -5.1111e-02 ,  4.6003e+00 , -2.0534e+01 ,  2.7731e+01 , &
         & 7.4514e-03 , -6.2529e-01 ,  2.9159e+00 , -4.1101e+00 , &
!  j = 15;  N + N --> N + N + pi; pion distributions:
         & 7.1622e-02 ,  3.0960e+00 , -1.1125e+01 ,  1.8130e+01 , &
         & 9.2581e-02 , -3.2186e+00 ,  2.0273e+01 , -3.3245e+01 , &
         & -5.1531e-02 ,  8.9886e-01 , -7.5084e+00 ,  1.3188e+01 , &
         & 5.8258e-03 , -1.7288e-03 ,  7.0224e-01 , -1.4854e+00 , &
!  j = 16;  N + N --> N + N + n*pi, n > 1; nucleon distributions:
         & 8.2300e-02 ,  1.5854e-01 ,  3.7716e+00 , -4.0562e+00 , &
         & 1.0802e-02 , -3.3688e-01 ,  1.1727e+00 , -6.7476e-01 , &
         & -2.1798e-03 ,  5.2166e-02 , -2.5816e-01 ,  3.2048e-01 , &
         & 6.5764e-05 , -1.4711e-03 ,  7.8209e-03 , -1.0580e-02 /
    data  a5/ &
!  j = 17;  N + N --> N + N + n*pi, n > 1; pion distributions:
         & 1.1138e-01 ,  6.0396e-01 ,  3.0174e+00 , -4.4190e+00 , &
         & -1.7709e-02 ,  2.3015e-01 , -1.8187e+00 ,  3.4518e+00 , &
         & 2.0977e-03 , -2.5458e-02 ,  2.1626e-01 , -4.0692e-01 , &
         & -5.4799e-05 ,  5.9111e-04 , -5.5552e-03 ,  1.0647e-02 , &
!  j = 18;  pi + N --> pi + N + pi; nucleon distributions:
         & 1.7288e-01 ,  7.1080e+00 , -1.7961e+01 ,  1.6403e+01 , &
         & -1.4504e-01 , -1.3032e+01 ,  4.1781e+01 , -4.0799e+01 , &
         & 4.5390e-02 ,  8.3515e+00 , -3.0260e+01 ,  3.2882e+01 , &
         & -4.7961e-03 , -1.4095e+00 ,  5.3505e+00 , -6.0946e+00 , &
!  j = 19;  pi + N --> pi + N + pi; pion distributions:
         & 3.7596e-02 ,  1.4331e+00 , -3.1350e+00 ,  6.4864e+00 , &
         & 2.3827e-01 ,  1.8253e+00 ,  1.7648e+00 , -1.6735e+01 , &
         & -1.5410e-01 , -1.5201e+00 , -1.5692e+00 ,  1.7185e+01 , &
         & 2.5037e-02 ,  3.0588e-01 ,  3.2520e-01 , -3.5277e+00 , &
!  j = 20;  pi + N --> pi + N + n*pi, n > 1; nucleon distributions:
         & 1.2489e-01 ,  1.3573e+00 ,  8.2338e-01 , -1.4595e+00 , &
         & -5.1577e-02 , -3.5778e-01 , -1.1690e+00 ,  1.8078e+00 , &
         & 7.4864e-03 ,  3.2888e-02 ,  2.3744e-01 , -3.9802e-01 , &
         & -2.9880e-04 , -7.5117e-04 , -1.1402e-02 ,  1.9505e-02 /
    data  a6/ &
!  j = 21;  pi + N --> pi + N + n*pi, n > 1; pion distributions:
         & 1.8470e-01 ,  1.9269e+00 , -3.2979e+00 ,  3.6843e+00 , &
         & -7.3932e-02 ,  2.7213e-01 ,  1.0600e+00 , -2.3354e+00 , &
         & 1.8907e-02 , -5.6473e-02 , -1.6487e-01 ,  3.8426e-01 , &
         & -9.2984e-04 ,  2.5506e-03 ,  7.3052e-03 , -1.7220e-02 /
    data  a7/ &
!  j = 22;  gamma + N --> pi0 + N, E-g <= 0.45 GeV.
         & 4.0693d-01 , -4.1404d+00 ,  1.4044d+01 , -1.7265d+01 , &
         & -3.6799d+00 ,  5.9610d+01 , -1.6269d+02 ,  1.8873d+02 , &
         & 1.4556d+01 , -1.7550d+02 ,  4.5839d+02 , -5.3390d+02 , &
         & -1.2621d+01 ,  1.4964d+02 , -3.8118d+02 ,  4.5141d+02 , &
!  j = 23;  gamma + N --> pi0 + N, E-g > 0.45 GeV.
         & -4.7554d-01 ,  2.2641d+00 , -1.2528d+01 ,  2.4647d+01 , &
         & 5.1620d+00 , -9.9236d+00 ,  5.5623d+01 , -1.0462d+02 , &
         & -8.1117d+00 ,  1.9315d+01 , -8.4255d+01 ,  1.3908d+02 , &
         & 3.5187d+00 , -9.1783d+00 ,  3.4950d+01 , -5.1243d+01 , &
!  j = 24; gamma + p --> n + pi+; E-g <= 0.51 GeV:
         & 4.8173d-01 ,  5.7726d+00 , -1.3745d+01 ,  2.7125d+01 , &
         & -4.4804d+00 , -3.8582d+01 ,  1.1159d+02 , -2.4305d+02 , &
         & 1.6306d+01 ,  1.1046d+02 , -3.3045d+02 ,  7.2270d+02 , &
         & -1.5968d+01 , -8.0140d+01 ,  2.4616d+02 , -6.0753d+02 , &
!  j = 25; gamma + p --> n + pi+; 0.51 < E-g <= 1.0 GeV:
         & -5.1646d+00 , -6.0776d+00 ,  7.8989d+01 , -1.0705d+02 , &
         & 2.1871d+01 ,  5.6915d+01 , -4.0159d+02 ,  5.1215d+02 , &
         & -2.7993d+01 , -9.4670d+01 ,  5.6928d+02 , -6.9621d+02 , &
         & 1.1587d+01 ,  4.5998d+01 , -2.4566d+02 ,  2.8452d+02 , &
!  j = 26; gamma + p --> n + pi+; 1.0 < E-g <= 10 GeV:
         & -5.3067d+01 ,  5.7612d+02 , -1.5438d+03 ,  1.6455d+05 , &
         & 1.4750d+02 , -1.6638d+03 ,  4.5923d+03 , -4.9949d+03 , &
         & -1.3436d+02 ,  1.5780d+03 , -4.4463d+03 ,  4.9022d+03 , &
         & 4.0253d+01 , -4.8860d+02 ,  1.4001d+03 , -1.5606d+03 , &
!  j = 27; gamma + N --> delta + pi; pion distribution; T < 1.0
         & -1.0306d+00 ,  3.2849d+01 , -7.5052d+01 ,  6.0255d+01 , &
         & 7.9586d+00 , -1.2572d+02 ,  2.5604d+02 , -1.6547d+02 , &
         & -1.4797d+01 ,  1.6590d+02 , -2.7991d+02 ,  1.1333d+02 , &
         & 8.2309d+00 , -6.7871d+01 ,  8.5762d+01 ,  5.9727d+00 , &
!  j = 28; gamma + N --> delta + pi; pion distribution; T > 1.0
         & -2.3722d+02 ,  9.6890d+02 , -1.6219d+03 ,  1.3637d+03 , &
         & 6.5800d+02 , -2.6941d+03 ,  4.5480d+03 , -3.8460d+03 , &
         & -6.0653d+02 ,  2.4983d+03 , -4.2498d+03 ,  3.6136d+03 , &
         & 1.8604d+02 , -7.6933d+02 ,  1.3166d+03 , -1.1242d+03 /
    data  a8/ &
!  j = 29; coefficients for K absorption ang. dist., Tlab <= 0.455
         & 6.5288d-01 ,  3.8977d-01 ,  8.4078d-01 ,  1.8893d-01 , &
         & -4.3964d+00 ,  3.4309d+01 , -7.3692d+01 ,  8.4308d+01 , &
         & 1.4889d+01 , -1.4380d+02 ,  3.1227d+02 , -3.5014d+02 , &
         & -1.5658d+01 ,  1.7160d+02 , -3.7212d+02 ,  4.1299d+02 /
    data  b1/ &
         & 5.0278e-01 ,  3.1442e+00 , -7.8172e+00 ,  8.1667e+00 , &
         & 9.3482e-01 , -1.0590e+01 ,  2.9227e+01 , -3.4550e+01 , &
         & -9.6685e-02 ,  4.7335e+00 , -1.4298e+01 ,  1.7685e+01 , &
         & -2.5041e-02 , -6.2478e-01 ,  2.0282e+00 , -2.5895e+00 , &
         & 1.1965e+00 , -8.2889e-01 ,  1.0426e+00 , -1.9090e+00 , &
         & 2.8703e-01 , -4.9065e+00 ,  1.6264e+01 , -1.9904e+01 , &
         & -2.4492e-01 ,  2.9191e+00 , -9.5776e+00 ,  1.1938e+01 , &
         & 3.7297e-02 , -4.2200e-01 ,  1.3883e+00 , -1.7476e+00 , &
         & 1.3508e+00 , -4.3139e+00 ,  1.2291e+01 , -1.5288e+01 , &
         & -2.0086e-01 ,  1.3641e+00 , -3.4030e+00 ,  3.8559e+00 , &
         & 1.2583e-02 , -8.3492e-02 ,  1.8600e-01 , -2.0043e-01 , &
         & -2.3628e-04 ,  1.3514e-03 , -2.4324e-03 ,  2.1906e-03 , &
         & 1.2419e+00 , -4.3633e+00 ,  1.3743e+01 , -1.8592e+01 , &
         & -2.4404e-01 ,  1.3158e+00 , -3.5691e+00 ,  4.3867e+00 , &
         & 1.5693e-02 , -8.2579e-02 ,  2.1427e-01 , -2.5846e-01 , &
         & -2.9386e-04 ,  1.4060e-03 , -3.3835e-03 ,  3.8664e-03 /
    data  b2/ &
         & 6.3054e-01 , -3.7333e+00 ,  1.3464e+01 , -1.8594e+01 , &
         & 2.1801e+00  ,  1.5163e+00 , -1.6380e+01 ,  2.7944e+01 , &
         & -1.2886e+00 , -2.4570e+00 ,  1.5129e+01 , -2.3295e+01 , &
         & 2.0915e-01 ,  5.2279e-01 , -2.8687e+00 ,  4.2688e+00 , &
         & 9.3363e-01 , -1.8181e+00 ,  5.5157e+00 , -8.5216e+00 , &
         & 1.7811e+00 , -8.2927e+00 ,  2.0607e+01 , -2.0827e+01 , &
         & -1.5264e+00 ,  6.8433e+00 , -1.6067e+01 ,  1.6845e+01 , &
         & 2.7128e-01 , -1.1944e+00 ,  2.7495e+00 , -2.9045e+00 , &
         & 1.9439e+00 , -4.6268e+00 ,  9.7879e+00 , -9.6074e+00 , &
         & -3.4640e-01 ,  1.1093e+00 , -1.9313e+00 ,  1.7064e+00 , &
         & 2.7054e-02 , -1.1638e-01 ,  2.6969e-01 , -3.1853e-01 , &
         & -6.6092e-04 ,  5.0728e-03 , -1.4995e-02 ,  1.9605e-02 , &
         & 1.8693e+00 , -5.5678e+00 ,  1.4795e+01 , -1.6903e+01 , &
         & -4.9965e-01 ,  1.7874e+00 , -4.1330e+00 ,  3.8393e+00 , &
         & 4.6194e-02 , -1.8536e-01 ,  4.5315e-01 , -4.6273e-01 , &
         & -1.3341e-03 ,  5.7710e-03 , -1.4554e-02 ,  1.5554e-02 /
    data  ckj/ &
         & 1.4509e-01 ,  4.6520e-01 , -3.3005e-02 ,  1.5376e-01 , &
         & 2.7436e-01 , -1.4604e-02 ,  6.2959e-01 ,  1.7866e-01 , &
         & -2.6216e-03 ,  8.3810e-01 ,  8.6137e-03 ,  3.2946e-03 , &
         & 9.2852e-02 ,  5.3886e-01 , -5.4493e-02 ,  1.3032e-01 , &
         & 4.0709e-01 , -2.8782e-02 ,  1.4909e-01 ,  3.8502e-01 , &
         & -1.2775e-02 ,  1.8024e-01 ,  3.3022e-01 , -9.4491e-03 /
  end
!     * * * * * * * * * * * * * * * * *
  double precision function sigmag (l,ms,mq,ksi,iks,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     SUBPROGRAM OF CHOOSING CROSS SECTION TYPE AND
!     CALCULATION CROSS SECTION VALUE FOR GIVEN ENERGY.
!     EDITED by KKG, October, 2003
!
    common /typecsq/icst,nsicst
    dimension icst(46),nsicst(37)
!     data ICST/
!  Translation:total:|elas:   |total: |elas:  |total:  |elas:   |pi+ n or
!              p+p or|p+p or  |       |       |pi- p or|pi- p or|pi- p  |
!              n+n   |n+n     |p+n    |n+p    |pi+ n   |pi+ n   |SCX    |
!    &         210,    211,    220,    221,     120,     121,     122,
!                           | Absorp: |
!  Translation:total:|elas: |pi+ np:pp|p + p->|p + p->|p + n->|p + n->|
!            pi+ p or|pi+ p |pi- pp:np|p p pi0|p n pi+|p + n +|p + p +|
!            pi- n   |pi- n |pi- pn:nn|n + n->|n + n->|pi0    |pi-    |
!                           |pi+ nn:np|n n pi0|p n pi-
!    &         110,    111,    123,    214,     215,     224,     225,
!
!  Transl:   pi+ p->|pi+ p->|pi- p ->|pi- + p->|pi- + p->|gam + p|gam +p|
!         pi0 pi+ p |n +2pi+|2pi0 + n|p + pi- +|n pi- pi+|-> p + |-> n +|
!            pi- n->|pi- n->|pi+ n ->|pi0      |pi+ + n->| pi0   |pi+
!         pi0 pi- n |p +2pi-|2pi0 + p|         |p pi- pi+| (?)   | (?)
!    &         114,    115,    126,     124,      125,    10111,  10112,
!            Absorp:
!  Transl:   gam +  |gam + p |gam + p |gam + p |gam + p |gam + p|gam + p|
!            2N ->  |-> pi+ p|->pi0 p |->n pi+ |-> delta|-> N + |total  |
!            2N     |+ pi-   |+ pi0   |+ pi0   |++ + pi-|2 pi   |       |
!    &       10113,  10115,  10114,   10116,    10118,   10117,   10110/
!
!  Transl:   total(elastic) for K + N
!            K+  + p or| K0  + p or| K-  + p or   |AK0 + p or  |
!            K0  + n   | K+  + n   |AK0  + n      | K- + n     |
!            1110(1111)| 1120(1121)| -1120(-1121) |-1110(-1111)|
!
!  Transl:   charge exchange for K + N
!            K+  + p or| K0  + p<=>| K-  + p <==> |AK0 + p or  |
!            K0  + n   | K+  + n   |AK0  + n      | K- + n     |
!              sig=0   |           |              |   sig=0.   |
!              1112    |    1122   |     -1122    |  -1112     |
!
!     data nsicst/   Absorp:  |  =0! |
!   Key:    pi+ p or|pi+ pn:pp|pi+ p:|pi+ p:|pi+ n:  |pi0 n |pi0 n  |
!            pi- n  |   or    |n 2pi0|N+ 2pi|N + 2pi | or   | or    |
!            SCX=0! |pi- pn:nn|pi- n:|pi- n:|pi- p:  |pi0 p |pi0 p  |
!                             |p 2pi0|N+ 2pi|N + 2pi |total |elastic|
!    &         112,    113,    116,    117,    127,    130,    131,
!                    Absorp:
!   Key:     pi0+p: |pi0+pn:pn|pi0+p:|pi0+p:  |pi0+p:  |pi0+p: |p + p: |
!            pi+ + n|pi0+pp:pp|p 2pi0|p pi+pi-|n pi+pi0|N + 2pi| SCX   |
!            pi0+n: |pi0+np:np|pi0+n:|pi0+n:  |pi0+n:  |pi0+n: |n + n: |
!            pi- + p|pi0+nn:nn|n 2pi0|n pi+pi-|p pi-pi0|N + 2pi| = 0!  |
!    &         132,    133,    134,    135,     136,     137,    212,
!
!   Key:     p + p: |p+p;n+n |p + p->|n + p: |n + p:  |n + p: |n + p->|
!      pi absorption|chg exch|2N + pi| SCX   | Pion   |n p pi0|N N pi |
!            n + n: |+ pi0   |n + n->|       |absorp. |Same as|       |
!             = 0!  |  = 0!  |2N + pi| = 0!  | = 0!   |  224! |       |
!    &         213,    216,    217,    222,     223,     226,    227/
!
!
!  Transl:      total(elastic) for Delta + N
!           D-  + p or| D0  + p or| D+  + p or| D++ + p or|
!           D++ + n   | D+  + n   | D0  + n   | D-  + n   |
!            230(231) |  240(241) |  250(251) |  260(261) |
!
!  Transl:     charge exchange for Delta + N
!           D-  + p =>| D0  + p <=>| D+  + p <=>| D++ + p or|
!           D0  + n or| D+  + n    | D++ + n or | D-  + n   |
!           D++ + n =>|            | D0  + n => |  sig=0.   |
!           D+  + p   |            | D-  +p     |           |
!            232      |  242       |     252    |  262      |
!
!  Transl:     absorption  (Delta + N ==> N + N)
!           D-  + p =>| D0  + p  =>| D+  + p =>| D++ + p or|
!           n   + n or| p   + n or | p   + p or| D-  + n   |
!           D++ + n =>| D+  + n  =>| D0  + n =>|  sig=0.   |
!           p   + p   | p   + n    | n   + n   | |         |
!            233      |  243       |     253   |  263      |
!
! ======================================================================
    data absgam/5.0/,abspi/4.0/
!   Start subroutine:

    ics = 10000*l+1000*iabs(ms)+100*mq+10*ksi+iks
    if(ms < 0)  ics=-ics
    js = 1
17  if (ics-icst(js)) 11,10,11
10  sigmag = qintg(t,js)
!     if(JS == 10)  SIGMAG = abspi*SIGMAG
    if(js == 40)  sigmag = absgam*sigmag
    return
11  if(js-46)12,13,13
12  js = js + 1
    go to 17
13  nsjs = 1
16  if (ics-nsicst(nsjs)) 15,14,15
15  if(nsjs-21)100,39,39
100 nsjs=nsjs+1
    go to 16
14  kns = nsjs
    go to (18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,&
         & 37,38)kns
18  sigmag = 0.
    return
19  sigmag = qintg(t,10)
    return
20  sigmag = 0.
    return
21  sigmag = qintg(t,15)+qintg(t,16)
    return
22  sigmag = qintg(t,18)+qintg(t,19)+qintg(t,17)
    return
23  sigmag = (qintg(t,8)+qintg(t,5))/2.
    return
24  sigmag = (qintg(t,9)+qintg(t,6)-qintg(t,7))/2.
    return
25  sigmag = qintg(t,7)
    return
26  sigmag = qintg(t,10)*0.5
    return
27  sigmag = (qintg(t,15)+qintg(t,18))/2.
    return
28  sigmag = (qintg(t,16)+qintg(t,19))/2.
    return
29  sigmag = (qintg(t,17))/2.
    return
30  sigmag = (qintg(t,15)+qintg(t,16)+qintg(t,18)+ &
         & qintg(t,19)+qintg(t,17))/2.
    return
31  sigmag = 0.
    return
32  sigmag = 0.
    return
33  sigmag = 0.
    return
34  sigmag = qintg(t,11)+qintg(t,12)
    return
35  sigmag = 0.
    return
36  sigmag = 0.
    return
37  sigmag = qintg(t,14)
    return
38  sigmag = qintg(t,13)+2.*qintg(t,14)
    return
!----> do  40  isd=1,16
39  do   isd=1,16
       ksd=isd
       if(ics == nsicst(21+isd))  go  to  41
40     continue
    end do
    go  to  35
41  continue
    w2=1.88*t+(1.232+0.940)**2
    tn=((sqrt(w2)+0.94-1.232)**2)/1.88-1.88
    pnn2=w2/4.-0.940**2
    pdn2=(w2+0.940**2-1.232**2)**2/(4.*w2)-0.940**2
    r=0.5*pnn2/pdn2
    hex=0.5
!  * * * * * * *  PRELIMINARY VERSION  * * * * * * * * * *
    go  to  (230,231,232,233, &
         & 240,241,242,243, &
         & 250,251,252,253, &
         & 260,261,262,263),ksd
230 sigmag=qintg(tn,3)
    return
240 go  to  230
250 sigmag=qintg(tn,1)
    return
260 go  to  250
231 sigmag=qintg(tn,4)*(1.-hex)
    return
241 go  to  231
251 sigmag=qintg(tn,2)*(1.-hex)
    return
261 go  to  251
232 sigmag=qintg(tn,4)*hex
    return
242 go  to  232
252 sigmag=qintg(tn,2)*hex
    return
262 go  to  263
233 sigmag=qintg(tn,12)/2.*r
    return
243 sigmag=(qintg(tn,12)/2.+qintg(tn,11))*r
    return
253 sigmag=(qintg(tn,14)/2.+qintg(tn,13)/2.)*r
    return
263 sigmag=0.
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  croseg(l,ms,mq,ksi,iks,t,am,ir)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Choose cross section type and calculate cross section for
!     a given struck-nucleon rest frame kinetic energy.
!     l always 1 for photons; ms always 0; mb = baryon number.
!     mb may be 1 or 2.
!     ksi = 1 for n - n, p - p, pi+ - p & pi- - n.
!     ksi = 2 for n - p, pi+ - n & pi- - p.
!     ksi = 3 for pi0 - n, pi0 - p
!     iks = 0:   total cross section
!     iks = 1:   elastic cross section
!     iks = 2:   pion charge exchange cross section
!     iks = 3:   pion or gamma absorption cross section
!     iks = 4:   neutral pion production cross section
!     iks = 5:   charged pion production cross section
!     iks = 6:   target & projectile isospin change & neutral pion
!                production (9/11/97)
!     iks = 7:   total one-pion production cross section (3/15/99)
!     iks = 8:   delta production (by gamma) cross section
!
    croseg=0.
    x=t
    w=sqrt(1.88*t+(1.88+am)**2)
    if(ir == 0.and.ms == 0)      go  to  10
    if(ir.ne.0.and.ms == 0)      go  to  11
    if(ir == 0.and.ms.ne.0)      go  to  12
    if(ir.ne.0.and.ms.ne.0)      go  to  14
    return
10  croseg=sigmag(l,ms,mq,ksi,iks,x)
    go  to  15
11  if(mq == 2)                  go  to  10
    if(iks >= 3)                 return
    t0=((w+0.140-am)**2-(0.940+0.140)**2)/1.88
    x=t0
    go  to  10
12  if(mq.ne.1)                  go  to  13
    p0=sqrt(t*(t+2.*am))
    x=p0
    go  to  10
13  if(iks > 1)                 return
    tn=((w+0.940-am)**2-(0.940+0.940)**2)/1.88
    tp=((w+0.140-am)**2-(0.940+0.140)**2)/1.88
    tk=((w+0.492-am)**2-(0.940+0.492)**2)/1.88
    pk=sqrt(tk*(tk+2.*0.492))
    spp  =qintg(tn, 1+iks)
    spn  =qintg(tn, 3+iks)
    spipp=qintg(tp, 8+iks)
    spimp=qintg(tp, 5+iks)
    skmp =qintg(pk,33+iks)
    skmn =qintg(pk,31+iks)
    slp=spn+skmp-spimp
    sln=spp+skmn-spipp
    if(ms == -1.and.ksi == 1)  croseg=   slp+spp-spn
    if(ms == -1.and.ksi == 2)  croseg=   slp+spn-spp
    if(ms == -1.and.ksi == 3)  croseg=   slp
    if(ms == -1.and.ksi == 4)  croseg=   sln+spn-spp
    if(ms == -1.and.ksi == 5)  croseg=   sln+spp-spn
    if(ms == -1.and.ksi == 6)  croseg=   sln
    if(ms == -2.and.ksi == 1)  croseg=2.*slp-spp
    if(ms == -2.and.ksi == 2)  croseg=2.*slp-spn
    if(ms == -2.and.ksi == 3)  croseg=2.*sln-spn
    if(ms == -2.and.ksi == 4)  croseg=2.*sln-spp
    if(ms == -3.and.ksi == 1)  croseg=3.*slp-spp-spn
    if(ms == -3.and.ksi == 2)  croseg=3.*sln-spp-spn
    go  to  15
14  if(mq == 2)           go  to  13
    t0=((w+0.492-am)**2-(0.940+0.492)**2)/1.88
    p0=sqrt(t0*(t0+2.*0.492))
    x=p0
    go  to  10
15  if(croseg < 0.)  croseg=0.
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine rotorq (ar,br,pstar,pr)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    BLOCK OF ROTATION.
    dimension ar(3),br(3),pstar(3),pr(3),an(3)
    sp = 0.
!----> do 31 ir=1,3
    do ir=1,3
       sp = sp+ar(ir)*br(ir)
31     continue
    end do
    amod = sqrt (ar(1)**2+ar(2)**2+ar(3)**2)
    if(amod < 1.e-30)  go  to  10
    alpha1 = sp/amod
    bmod2 = br(1)**2+br(2)**2+br(3)**2
    temp=bmod2-alpha1**2
    if(temp <= 0.d0)  go  to  10
    alpha2 = sqrt (temp)
    if(alpha2 < 1.e-30)  go  to  10
    an(1) = ar(2)*br(3)-ar(3)*br(2)
    an(2) = ar(3)*br(1)-ar(1)*br(3)
    an(3) = ar(1)*br(2)-ar(2)*br(1)
    pr(1)=pstar(1)*br(1)/alpha2+(pstar(3)-alpha1*pstar(1)/alpha2) &
         & *ar(1)/amod+(pstar(2)*an(1))/(alpha2*amod)
    pr(2)=pstar(1)*br(2)/alpha2+(pstar(3)-alpha1*pstar(1)/alpha2) & 
         & *ar(2)/amod+(pstar(2)*an(2))/(alpha2*amod)
    pr(3)=pstar(1)*br(3)/alpha2+(pstar(3)-alpha1*pstar(1)/alpha2) & 
         & *ar(3)/amod+(pstar(2)*an(3))/(alpha2*amod)
    return
!----> do  11  k=1,3
10  do   k=1,3
11     pr(k)=pstar(k)
    end do
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  subroutine binelq(partin,ipatin,ipatne,l,ms,mq,ksi,me,v,u, &
       & tin1,mv,np,nin)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     BLOCK OF INELASTIC SCATTERING at low energies TIN1 < 5-20 GeV
!
    common /i3act/ sig3,sigin,ind,ith
    common/ncasca/ncas,ncpri
    dimension partin(9),ipatin(5),ipatne(5),v(3)
    data uthr/1.220/
    nin = 0
    ik = 0
!  kkg 03/15/04
    if(ipatin(2).ne.0)  then
       if(u < uthr)  then
          nin = 2
          return
       endif
!                    gamma + N => Delta+pi or  => N+2pi
       betais = qintg(tin1,44)/qintg(tin1,45)
       if(rndm(-1.0_real64) <= betais)  then
!                    gamma + N => Delta+pi
          call isobal(u,v,tin1,partin,ipatne,mv,np)
          return
       else
!                    gamma + N  => N+2pi
          call statl(u,v,partin,ipatne,mv,np)
          return
       endif
    endif

!  kkg 03/15/04
    am3=partin(9)
    if(ipatin(5) == 0)  go  to  11
    if(ipatin(4) == 1)  am3=0.940
    go  to  14
11  if(tin1-4.)12,24,24
24  betath=0.
    go to 13
12  sig3=  croseg(l,ms,mq,ksi,7,tin1,partin(9),ipatin(5))
    sigin= croseg(l,ms,mq,ksi,0,tin1,partin(9),ipatin(5))- &
         & croseg(l,ms,mq,ksi,1,tin1,partin(9),ipatin(5))- &
         & croseg(l,ms,mq,ksi,2,tin1,partin(9),ipatin(5))
    betath=sig3/sigin
13  drnd=rndm(-1.0_real64)
    if(drnd-betath)14,15,15
14  ith=1
    th=1.
    go to 16
15  ith=0
    th=2.
16  if(u- am3     -0.14*th-0.96)19,19,17
17  continue
    if(ncas >= ncpri) write( *,101)
    call vmnspq (partin,ipatin,u,mv,np,ith,mq,tin1,lp)
    if (np) 26,26,27
26  return
27  continue
    if (lp) 19,22,19
22  continue
    if(ncas >= ncpri) write( *,102)
    call directq (v,tin1,mq,mv,np,partin,kp,ith)
    if (kp) 23,18,23
18  continue
    if(ncas >= ncpri) write( *,103)
    call chinelq (ipatin,l,ms,mq,ksi,np,mv,tin1,me,ipatne,partin(9))
    if(ipatin(5) == 0.and.ith == 1.and.ncas >= ncpri) write( *,104)
    if(ipatin(5) == 0.and.ith == 1) &
         & call  disob(mv,u,ind,np,mq)
    return
23  ik = ik+1
    if(ik < 50)   go  to  17
19  nin=2
101 format(' ====> vmnspq')
102 format(' ====> directq')
103 format(' ====> chinel')
104 format(' ====> disob')
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  subroutine statl(u,v,partin,ipatne,mv,np)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     determining of secondary particles characteristics for
!     gamma-n interaction with statistical model.
! last modification: 15 Mar. 2004 by KKG
!
!
    dimension v(3),pv1(3),pv2(3),pv3(3),ps1(3),ps2(3),ps3(3),pin(3), &
         & pinst(3),pmemo(9,5999),imemo(5,5999),partin(9),ipatne(5)
    common /memorylaq/ pmemo,imemo
    data emnucm, emnucg, emnucb / 938.919, 0.938919, 0.9315014/
    data emneut, emprot, empich, empi0 / 0.9395656, 0.9382723, &
         & 0.139568, 0.134973/
    data  pi /3.1415926536d0/
!
!  determine randomly the resulting particle types
    thrd = 1.d0/3.d0
    twthrd = 2.d0/3.d0
    temp4 = rndm(-1.0_real64)
    imemo(2,mv+1) = 0
    imemo(3,mv+1) = 0
    imemo(4,mv+1) = 0
    imemo(5,mv+1) = 0
    imemo(2,mv+2) = 0
    imemo(3,mv+2) = 0
    imemo(4,mv+2) = 0
    imemo(5,mv+2) = 0
    imemo(2,mv+3) = 0
    imemo(3,mv+3) = 0
    imemo(4,mv+3) = 1
    imemo(5,mv+3) = 0
    if (ipatne(1) <= 0) then
       if (temp4 <= thrd) then
          imemo(1,mv+1) = 0
          pmemo(9,mv+1) = empi0
          imemo(1,mv+2) = 0
          pmemo(9,mv+2) = empi0
          imemo(1,mv+3) = 0
          pmemo(9,mv+3) = emneut
       elseif (temp4 < twthrd) then
          imemo(1,mv+1) = 0
          pmemo(9,mv+1) = empi0
          imemo(1,mv+2) = -1
          pmemo(9,mv+2) = empich
          imemo(1,mv+3) = 1
          pmemo(9,mv+3) = emprot
       else
          imemo(1,mv+1) = -1
          pmemo(9,mv+1) = empich
          imemo(1,mv+2) = 1
          pmemo(9,mv+2) = empich
          imemo(1,mv+3) = 0
          pmemo(9,mv+3) = emneut
       endif
    else
       if (temp4 <= thrd) then
          imemo(1,mv+1) = 0
          pmemo(9,mv+1) = empi0
          imemo(1,mv+2) = 0
          pmemo(9,mv+2) = empi0
          imemo(1,mv+3) = 1
          pmemo(9,mv+3) = emprot
       elseif (temp4 <= twthrd) then
          imemo(1,mv+1) = 0
          pmemo(9,mv+1) = empi0
          imemo(1,mv+2) = 1
          pmemo(9,mv+2) = empich
          imemo(1,mv+3) = 0
          pmemo(9,mv+3) = emneut
       else
          imemo(1,mv+1) = 1
          pmemo(9,mv+1) = empich
          imemo(1,mv+2) = -1
          pmemo(9,mv+2) = empich
          imemo(1,mv+3) = 1
          pmemo(9,mv+3) = emprot
       endif
    endif
    np = 3
!
!  other combinations might be possible!?
    twopi = 2.0d0*pi
    empi = max(pmemo(9,mv+1),pmemo(9,mv+2))
    empi2 = empi**2
    emnu = pmemo(9,mv+3)
    emnupi = emnu + min(pmemo(9,mv+1),pmemo(9,mv+2))
!
    epim = (u**2+empi2-emnupi**2)/(2.d0*u)
    tpim = epim-empi
10  t1 = rndm(-1.0_real64)*tpim
    t2 = rndm(-1.0_real64)*tpim
    e1 = t1+pmemo(9,mv+1)
    e2 = t2+pmemo(9,mv+2)
    f1 = 27.d0*e1*e2*(u-e1-e2)/(u**3)
    b1 = rndm(-1.0_real64)
    if (b1 < f1)then
       e3 = u-e1-e2
       t3 = e3-pmemo(9,mv+3)
    else
       goto 10
    endif
    if (t3 <= 0.0) goto 10
    p1 = sqrt(t1*(t1+2.0d0*pmemo(9,mv+1)))
    p2 = sqrt(t2*(t2+2.0d0*pmemo(9,mv+2)))
    p3 = sqrt(t3*(t3+2.0*pmemo(9,mv+3)))
    temp1 = (p1+p2-p3)*(p1-p2+p3)*(p2+p3-p1)
    if (temp1 <= 0.0) goto 10
    ct3 = 1.d0-2.d0*rndm(-1.0_real64)
    fi3 = twopi*rndm(-1.0_real64)
    temp2 = sqrt(1.d0-ct3**2)
    pv3(1) = p3*temp2*cos(fi3)
    pv3(2) = p3*temp2*sin(fi3)
    pv3(3) = p3*ct3
    temp3=sqrt(partin(8)*(partin(8)+2.d0*partin(9)))
    pin(1)=temp3*partin(4)*partin(7)
    pin(2)=temp3*partin(4)*partin(6)
    pin(3) = temp3*partin(5)
    call cmsq(pin,v,pinst,partin(8),partin(9))
    call rotorq (pinst,v,pv3,ps3)
    ct1 = (p2**2-p1**2-p3**2)/(2.d0*p3*p1)
    ct2 = (p1**2-p2**2-p3**2)/(2.d0*p3*p2)
    fi1 = twopi*rndm(-1.0_real64)
    fi2 = pi+fi1
    st1 = sqrt(1.d0-ct1**2)
    st2 = sqrt(1.d0-ct2**2)
    pv1(1) = p1*st1*cos(fi1)
    pv1(2) = p1*st1*sin(fi1)
    pv1(3) = p1*ct1
    call rotorq (ps3,v,pv1,ps1)
    pv2(1) = p2*st2*cos(fi2)
    pv2(2) = p2*st2*sin(fi2)
    pv2(3) = p2*ct2
    call rotorq (ps3,v,pv2,ps2)
    pmemo(4,mv+1) = ps1(1)
    pmemo(5,mv+1) = ps1(2)
    pmemo(6,mv+1) = ps1(3)
    pmemo(7,mv+1) = 0.
    pmemo(4,mv+2) = ps2(1)
    pmemo(5,mv+2) = ps2(2)
    pmemo(6,mv+2) = ps2(3)
    pmemo(7,mv+2) = 0.
    pmemo(4,mv+3) = ps3(1)
    pmemo(5,mv+3) = ps3(2)
    pmemo(6,mv+3) = ps3(3)
    pmemo(7,mv+3) = 0.
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  subroutine isobal (u,v,tin1,partin,ipatne,mv,np)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     determining of secondary particles characteristics for gamma-n
!     interaction with (3/2,3/2) isobar production.
!
!   Last change: 15-MAR-2004 by KKG
!
    dimension vt(3),ppim(3),ppt(3),ppit(3),ppi(3),pp(3),pin(3), &
         & pinst(3),ppimst(3),ppist(3),ppst(3), &
         & pmemo(9,5999),imemo(5,5999), &
         & partin(9),ipatne(5),v(3)
    common /memorylaq/ pmemo,imemo
    data emnucm, emnucg, emnucb / 938.919, 0.938919, 0.9315014/
    data emneut, emprot, empich, empi0 / 0.9395656, 0.9382723, &
         & 0.139568, 0.134973/
    data  pi /3.1415926536d0/
!
    imemo(2,mv+1) = 0
    imemo(3,mv+1) = 0
    imemo(4,mv+1) = 0
    imemo(5,mv+1) = 0
    imemo(2,mv+2) = 0
    imemo(3,mv+2) = 0
    imemo(4,mv+2) = 0
    imemo(5,mv+2) = 0
    imemo(2,mv+3) = 0
    imemo(3,mv+3) = 0
    imemo(4,mv+3) = 1
    imemo(5,mv+3) = 0
    if (ipatne(1) > 0) then
       imemo(1,mv+1) = -1
       pmemo(9,mv+1) = empich
       imemo(1,mv+2) = 1
       pmemo(9,mv+2) = empich
       imemo(1,mv+3) = 1
       pmemo(9,mv+3) = emprot
    else
       imemo(1,mv+1) = 1
       pmemo(9,mv+1) = empich
       imemo(1,mv+2) = -1
       pmemo(9,mv+2) = empich
       imemo(1,mv+3) = 0
       pmemo(9,mv+3) = emneut
    endif
    np = 3
!
    twopi = 2.0d0*pi
    empi2 = empich**2
    emnu = pmemo(9,mv+3)
    emnupi = emnu + empich
    emnu2pi = emnupi + empich
    a1 = (u**2+empi2-emnupi**2)/(2.d0*u)
    a2 = sqrt(a1**2-empi2)
    a3 = u-a1
    f1 = a1*a2*a3/u
    alpha = 200.d0*f1
10  bms = rndm(-1.0_real64)*(u-emnu2pi)+emnupi
    epim = (u**2+empi2-bms**2)/(2.d0*u)
    edn = u-epim
    pim = sqrt(epim**2-empi2)
    f = (pim*epim*edn)/u
    ts = (bms**2-emnupi**2)/emnu/2.d0
    p = f*qintg(ts,9)/alpha
    b1 = rndm(-1.0_real64)
    if (p > b1) then
       tpi = epim-empich
    else
       goto 10
    endif
    if (tin1 < 1.d0) then
       ctpi = costaq(27,tin1)
    else
       ctpi = costaq(28,tin1)
    endif
    fipi = twopi*rndm(-1.0_real64)
    epit = (bms**2+empi2-emnu**2)/(2.d0*bms)
    ent = bms-epit
    temp1 = sqrt(1.d0-ctpi**2)
    ppim(1) = pim*temp1*cos(fipi)
    ppim(2) = pim*temp1*sin(fipi)
    ppim(3) = pim*ctpi
    temp2 = tpi+empich-u
    vt(1) = ppim(1)/temp2
    vt(2) = ppim(2)/temp2
    vt(3) = ppim(3)/temp2
    ctilpi=1.d0-2.d0*rndm(-1.0_real64)
    ftilpi=twopi*rndm(-1.0_real64)
    temp1 = sqrt(1.d0-ctilpi**2)
    pmt = sqrt(epit**2-empi2)
    ppit(1) = pmt*temp1*cos(ftilpi)
    ppit(2) = pmt*temp1*sin(ftilpi)
    ppit(3) = pmt*ctilpi
    tpi = epit-empich
    vt(1)=-vt(1)
    vt(2)=-vt(2)
    vt(3)=-vt(3)
    call cmsq (ppit,vt,ppi,tpi,empich)
    ppt(1) = -ppit(1)
    ppt(2) = -ppit(2)
    ppt(3) = -ppit(3)
    tn = ent-emnu
    call cmsq (ppt,vt,pp,tn,emnu)
    temp2 = sqrt(partin(8)*(partin(8)+2.d0*partin(9)))
    pin(1)=temp2*partin(4)*partin(7)
    pin(2)=temp2*partin(4)*partin(6)
    pin(3) = temp2*partin(5)
    call cmsq(pin,v,pinst,partin(8),partin(9))
    call rotorq(pinst,v,ppim,ppimst)
    pmemo(4,mv+1) = ppimst(1)
    pmemo(5,mv+1) = ppimst(2)
    pmemo(6,mv+1) = ppimst(3)
    pmemo(7,mv+1) = 0.
    call rotorq(pinst,v,ppi,ppist)
    pmemo(4,mv+2) = ppist(1)
    pmemo(5,mv+2) = ppist(2)
    pmemo(6,mv+2) = ppist(3)
    pmemo(7,mv+2) = 0.
    call rotorq(pinst,v,pp,ppst)
    pmemo(4,mv+3) = ppst(1)
    pmemo(5,mv+3) = ppst(2)
    pmemo(6,mv+3) = ppst(3)
    pmemo(7,mv+3) = 0.
!
    return
  end

!
!     ********************************************************************
!

  subroutine cmsq (p, v, pstar, t, cm)

! ======================================================================
!
!     Momentum calculation in system which has a relative velocity
!     v to given one.
!     Lorentz transformation; see Jackson, 2nd ed., p541.
!
!   Called by: ISOBAR STAT
!
!   Edited by AJS, August, 1997.
!
!   Last change: 15-MAR-2004 by KKG
! ======================================================================


    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit double precision (a-h,o-z), integer(int32) (i-n)

    dimension p(3), v(3), pstar(3)

! ======================================================================

    v2 = v(1)**2 + v(2)**2 + v(3)**2
    temp1 = sqrt(1.0 - v2)
    spv = p(1)*v(1) + p(2)*v(2) + p(3)*v(3)
    temp2 = spv/v2*(1.0/temp1 - 1.0)
    sv = (t + cm)/temp1

! should bellow be "+" instead of "-" ?, SGM, 05/25/03

    pstar(1) = p(1) + v(1)*temp2 - v(1)*sv
    pstar(2) = p(2) + v(2)*temp2 - v(2)*sv
    pstar(3) = p(3) + v(3)*temp2 - v(3)*sv
    return

! ======================================================================
  end

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  cosal(t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    sampling cos(theta) for elastic N+He-4 scattering
!
    if(t-0.147) 1,1,2
1   a1=25.2*t**(-0.843)
    go to 3
2   a1=130.*t**0.0145
3   a2=11.3*t**0.432
    if(t-0.055) 4,4,5
4   a3=0.22*t**(-1.35)
    go to 6
5   a3=0.000043*t**(-4.32)
6   a4=130.*t**1.33
    temp1=1.-exp(-a2*3.141592)
    temp3=2.-temp1
    temp2=1.-exp(-a4*3.141592)
    temp4=2.-temp2
    w=a3*temp4/((1.+a4**2)*(a1*temp3/(1.+a2**2)+a3*temp4/(1.+a4**2)))
    drnd=rndm(-1.0_real64)
    if(drnd-w) 7,7,9
7   tau=-log(1.-rndm(-1.0_real64)*temp2)/a4
    drnd=rndm(-1.0_real64)
    if(drnd-sin(tau))8,8,7
8   teta=3.141592-tau
    go to 10
9   teta=-log(1.-rndm(-1.0_real64)*temp1)/a2
    drnd=rndm(-1.0_real64)
    if(drnd-sin(teta)) 10,10,9
10  cosal= cos(teta)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function costaq(j,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Cosine calculation for elastic and charge-exchange reactions.
!
    common/coefaq/ankj
!   kkg  10/28/03
    dimension ankj(4,4,29),ank(4,4)
    if(t < 0.001)  go  to  14
!----> do 10 k=1,4
    do k=1,4
!----> do 10 n=1,4
       do n=1,4
          ank(n,k) = ankj(n,k,j)
10        continue
       end do
    end do
    s1 = 0.
    r1 = rndm(-1.0_real64)
    if(r1 < 1.e-10) r1 = rndm(-1.0_real64)
    s2 = 0.
!----> do 11 n=1,4
    do n=1,4
!----> do 11 k=1,4
       do k=1,4
          s1 = s1+ank(n,k)*(t**(k-1))*(r1**(n-1))
11        continue
       end do
    end do
!----> do 12 n=1,4
    do n=1,4
!----> do 12 k=1,4
       do k=1,4
          s2 = s2+ank(n,k)*t**(k-1)
12        continue
       end do
    end do
    cta = 2.*sqrt(r1)*(s1+(1.-s2)*r1**4)-1.
    temp1 = abs(cta)
    if (temp1-1.) 13,13,14
13  costaq = cta
    return
14  costaq=1.-2.*rndm(-1.0_real64)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function qintg (x,lq)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
    logical ilog
    common /tabele/sigma,argus
    dimension sigma(30,46),argus(30,9)
!      Parabolic interpolation of tabulated cross sections.
!      For lq = 1-4 and 15-19, use a parabola in log(sigma) vs. log(E)
!      for energies (x) less than the 2nd point.
!   Modified for new cross section tables using AJS's sugestion from
!   cem2k, October, 2003, KKG
! ======================================================================

    lpmax = 30
    ilog = .false.
    if (lq <= 4) then
!  1-4
       i = 1
!  kkg 24.03.06
       if(x > argus(30,i))  then
          qintg = hexsecg(x,lq)
          return
       endif
       ilog = .true.
    elseif (lq <= 10) then
!  5-10
       i = 2
!  kkg 24.03.06
       if(lq >= 5.and.lq <= 9.and.x > argus(30,i))  then
          qintg = hexsecg(x,lq)
          return
       endif
    elseif (lq <= 14) then
!  11-14
       i = 3
    elseif (lq <= 19) then
!  15-19
       i = 3
       ilog = .true.
    elseif(lq <= 23)  then
!  20-23
       i = 4
    elseif(lq <= 26)  then
!  24-26
       i = 5
    elseif(lq <= 37) then
       i = 6
!  27-37
    elseif(lq <= 39.or.lq == 46)  then
!  38,39,46
       i = 7
       ilog = .true.
    elseif(lq == 40)  then
!  40
       i = 8
    elseif(lq <= 45)  then
!  41-45
       i = 9
       ilog = .true.
    else
       write(*,1000)  lq
       stop
    endif
! ---------------------------------
    lpha = 1
10  continue
    if (x == argus(lpha,i)) then
!  If x = table value, set cross section to tabulated one.
       qintg = sigma(lpha,lq)
       return
    elseif (x > argus(lpha,i)) then
!  Increment lpha until x is <= tabulated x.
       if (lpha >= lpmax-1) then
          lpha = lpmax - 1
          go to 20
       else
          lpha = lpha + 1
          go to 10
       endif
    elseif (x < argus(lpha,i)) then
       if (lpha <= 1 .and. .not.ilog .and. i == 3) then
!  If x < first table value and not using parabolic log interpolation;
!  below threshold; set cross section to 0.
          qintg = 0.0
          return
       elseif (lpha <= 1 .and. i == 2) then
!  For pi+N reactions with energy < 0; use-energy cross section:
          qintg = sigma(1,lq)
          return
       else
          lpha = max(2,lpha)
          lpha = min(lpmax-1,lpha)
       endif
    endif
20  phi1 = sigma(lpha-1,lq)
    psi1 = argus(lpha-1,i)
    phi2 = sigma(lpha,lq)
    psi2 = argus(lpha,i)
    phi3 = sigma(lpha+1,lq)
    psi3 = argus(lpha+1,i)
    if (ilog .and. lpha == 2) then
       x = log(x)
       phi1 = log(sigma(lpha-1,lq))
       psi1 = log(max(argus(lpha-1,i), 1.e-6))
       phi2 = log(sigma(lpha,lq))
       psi2 = log(argus(lpha,i))
       phi3 = log(sigma(lpha+1,lq))
       psi3 = log(argus(lpha+1,i))
    endif
    a = psi2 - psi3
    b = psi3 - psi1
    c = psi1 - psi2
    delta = a*psi1**2 + b*psi2**2 + c*psi3**2
    deltaa = phi1*a + phi2*b + phi3*c
    deltab = (phi2 - phi3)*psi1**2 + (phi3 - phi1)*psi2**2 + &
         & (phi1 - phi2)*psi3**2
    deltac = (psi2*phi3 - psi3*phi2)*psi1**2 + &
         & (psi3*phi1 - psi1*phi3)*psi2**2 + &
         & (psi1*phi2 - psi2*phi1)*psi3**2
    fact = 0.0
    if (delta.ne.0.0) fact = 1.0/delta
    a = deltaa*fact
    b = deltab*fact
    c = deltac*fact
    qintg = a*x**2 + b*x + c
    if (ilog .and. lpha == 2) then
       x = exp(x)
       qintg = exp(qintg)
    else
       qintg = max(qintg, 0.0)
    endif
    return

! ======================================================================

1000 format (1x,'qintg called with improper value of lq = ',i4)

! ======================================================================
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function hexsecg(t,i)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!  kkg  10.03.06, 04/12/07
!  high energy approximation of total and elastic cross sections
!  pi+p, pi-p, pp,np
    real(real64) ::  mpi,mn
    dimension p1(8),p2(8),p3(8),p4(8)
    data mpi/0.139/,mn/0.939/
    data s1/1.00/,ss0/5.380/,et1/0.4580/,et2/0.5450/
    data p1/19.912,-0.5032,21.122,-1.4187,35.530,5.3828,34.425,35.80/
    data p2/5.4529,8.85970,4.7596,11.2230,4.3527,3.9741,9.6702,6.336/
    data p3/3.6838,8.66200,1.2521,11.5740,0.6752,0.4921,9.8362,5.477/
    data p4/0.6230,0.48470,0.5390,0.50790,0.5555,0.2944,0.6545,0.308/
    if(i <= 4)                then
       k=i+4
       s=2.0*t*mn+(mn+mn)**2
    elseif(i == 5.or.i == 6)  then
       k=i-2
       s=2.0*t*mn+(mpi+mn)**2
    elseif(i >= 8)            then
       k=i-7
       s=2.0*t*mn+(mpi+mn)**2
    else
    endif
    s0=ss0**2
    z = p1(k)
    y1= p2(k)**2
    y2= p3(k)**2
    b = p4(k)**2
    hexsecg=z+b*log(s/s0)**2+y1*(s1/s)**et1-y2*(s1/s)**et2
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine slqekq(l,ms,mq,ksi,me,iin,ipn)
    use, intrinsic:: iso_fortran_env, only: int32, real64
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     Determine cross section type.
!     Edited by KKG, October, 2003
    dimension iin(5),ipn(5)
!   |  i   |     IIN(i) or IPN(i) translation
!   |  1   |             charge
!   |  2   |           lepton number
!   |  3   |             strangeness
!   |  4   |           baryon number
!   |  5   |zero for stable particles, INT(1000*t_life) for resonances
!
    mein=iin(1)
!  Total strangenes in system:
    ms = iin(3)+ipn(3)
!  Total lepton number in system:
    l =  iin(2)+ipn(2)
!  Total baryon number in system:
    mq = iin(4)+ipn(4)
!  Total charge in system:
    me = iin(1)+ipn(1)
!    condition for Delta's
    if(iin(5).ne.0.and.iin(4) == 1.and.iin(3) == 0)  go  to  30
9   continue
    if (ms) 10,11,10
10  if(iin(4).ne.0)    go  to  110
!    strange mesons + p
    if(ms ==  1.and.iin(1) ==  1.and.ipn(1) == 1)    ksi=1
    if(ms ==  1.and.iin(1) ==  0.and.ipn(1) == 1)    ksi=2
    if(ms == -1.and.iin(1) == -1.and.ipn(1) == 1)    ksi=2
    if(ms == -1.and.iin(1) ==  0.and.ipn(1) == 1)    ksi=1
!    strange mesons + n
    if(ms ==  1.and.iin(1) ==  1.and.ipn(1) == 0)    ksi=2
    if(ms ==  1.and.iin(1) ==  0.and.ipn(1) == 0)    ksi=1
    if(ms == -1.and.iin(1) == -1.and.ipn(1) == 0)    ksi=1
    if(ms == -1.and.iin(1) ==  0.and.ipn(1) == 0)    ksi=2
    return
!    strange baryon + p or n
110 if(ms == -1.and.iin(1) ==  1.and.ipn(1) == 1)    ksi=1
    if(ms == -1.and.iin(1) ==  1.and.ipn(1) == 0)    ksi=4
    if(ms == -1.and.iin(1) == -1.and.ipn(1) == 1)    ksi=2
    if(ms == -1.and.iin(1) == -1.and.ipn(1) == 0)    ksi=5
    if(ms == -1.and.iin(1) ==  0.and.ipn(1) == 1)    ksi=3
    if(ms == -1.and.iin(1) ==  0.and.ipn(1) == 0)    ksi=6
    if(ms == -2.and.iin(1) == -1.and.ipn(1) == 1)    ksi=1
    if(ms == -2.and.iin(1) == -1.and.ipn(1) == 0)    ksi=3
    if(ms == -2.and.iin(1) ==  0.and.ipn(1) == 1)    ksi=2
    if(ms == -2.and.iin(1) ==  0.and.ipn(1) == 0)    ksi=4
    if(ms == -3.and.iin(1) == -1.and.ipn(1) == 1)    ksi=1
    if(ms == -3.and.iin(1) == -1.and.ipn(1) == 0)    ksi=2
    return
11  if (l) 13,13,12
!   gamma + p or n
12  ksi = 1
    return
13  if(mq-1) 15,15,100
!   baryon on baryon
100 if(mq-2) 14,14,16
14  if (me-1)16,17,16
!   2 neutrons or 2 protons
16  ksi = 1
    return
!   1 neutron and 1 proton
17  ksi = 2
    return
!   meson interacting with baryon
15  if (me-2) 19,18,19
!   pi+ on proton
18  ksi = 1
    return
19  if (me+1) 21,20,21
!   pi- on neutron
20  ksi = 1
    return
21  if (me) 22,23,22
22  if (mein-1) 27,26,27
!   pi+ incident on neutron
26  ksi = 2
    return
!   pi0 incident on proton
27  ksi = 3
    return
23  if (mein+1) 25,24,25
!   pi- incident on proton
24  ksi = 2
    return
!   pi0 incident on neutron
25  ksi = 3
    return
!   Delta + proton or neutron
!                    p             n
!     D-             3             6
!     D0             4             5
!     D+             5             4
!     D++            6             3
!
30  ksi=5+iin(1)*(ipn(1)-1)+ipn(1)*(iin(1)-1)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine cinema(pstar,v,p,ct,st,cfi,sfi,t,cm)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     KINEMATIC BLOCK.
    dimension pstar(3),v(3),p(3)
    spv = pstar(1)*v(1)+pstar(2)*v(2)+pstar(3)*v(3)
    temp1 = 1.d0-v(1)**2-v(2)**2-v(3)**2
    ! call err_chk(1,'laq1.f','8924',2,temp1)
    temp2 = sqrt(temp1)
    ! call err_chk(1,'laq1.f','8926',1,temp2)
    g=1./temp2
    temp3 = pstar(1)**2+pstar(2)**2+pstar(3)**2+cm**2
    ! call err_chk(1,'laq1.f','8929',2,temp3)
    estar=sqrt(temp3)
    temp4 = g+1.
    ! call err_chk(1,'laq1.f','8933',1,temp4)
    do  k=1,3
       p(k)=pstar(k)+g*v(k)*(spv*g/(temp4)+estar)
    enddo
    temp1 = p(1)**2+p(2)**2+p(3)**2
    ! call err_chk(1,'laq1.f','8937',2,temp1)
    pm = sqrt(temp1)
    if(pm < 1.e-30)   then
       ct=1.d0
       cfi=1.d0
       sfi=0.d0
       st=0.d0
    else
       ct = p(3)/pm
       temp4 = 1.d0-ct**2
       if(temp4 < 1.e-30) then
          ct=1.d0
          cfi=1.d0
          sfi=0.d0
          st=0.d0
       else
          st = sqrt(temp4)
          temp3 = pm*st
          cfi=p(1)/temp3
          sfi=p(2)/temp3
       endif
    endif
    temp1 = pm**2+cm**2
    ! call err_chk(1,'laq1.f','8960',2,temp1)
    t=sqrt(temp1)-cm
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function pmomq (j,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    real(real64) ::  mn,mpi
!     BLOCK OF CALCULATION OF SECONDARY PARTICLES MOMENTUM.
    common/coefbcq/bnkj,ckj
    dimension bnkj(4,4,8),ckj(3,8),bnk(4,4)
    data mn/0.939/,mpi/0.139/
!
!----> do 10 k=1,4
    do k=1,4
!----> do 10 n=1,4
       do n=1,4
          bnk(n,k) = bnkj(n,k,j)
10        continue
       end do
    end do
    s1 = 0.
    r1 = rndm(-1.0_real64)
    s2 = 0.
    pmax = 0.
!----> do 11 n=1,4
    do n=1,4
!----> do 11 k=1,4
       do k=1,4
          s1 = s1+bnk(n,k)*(t**(k-1))*(r1**(n-1))
11        continue
       end do
    end do
!----> do 12 n=1,4
    do n=1,4
!----> do 12 k=1,4
       do k=1,4
          s2 = s2+bnk(n,k)*t**(k-1)
12        continue
       end do
    end do
!  calculation of maximal momentum
!  back to old version, KKG, Dec. 2004
    do k=1,3
       pmax = pmax+ckj(k,j)*t**(k-1)
    enddo
!
    pmomq = pmax*sqrt(r1)*(s1+(1.-s2)*r1**4)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  integer function jtypaq (ith,mq,lamb)
    use, intrinsic:: iso_fortran_env, only: int32, real64
    implicit real(real64)(a-h,o-z), integer(int32) (i-n)
!     DETERMINING OF TYPE OF COEFFICIENTS A(N,K).
    if (ith) 10,18,10
10  if (mq-1) 11,11,14
11  if (lamb-1) 13,13,12
12  jtypaq = 19
    return
13  jtypaq = 18
    return
14  if (lamb-1) 16,16,15
15  if (lamb-3) 17,16,17
16  jtypaq = 14
    return
17  jtypaq = 15
    return
18  if (mq-1) 19,19,22
19  if (lamb-1) 20,20,21
20  jtypaq = 20
    return
21  jtypaq = 21
    return
22  if (lamb-1) 24,24,23
23  if (lamb-3) 25,24,25
24  jtypaq = 16
    return
25  jtypaq = 17
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  integer function jtypbq (ith,mq,lamb)
    use, intrinsic:: iso_fortran_env, only: int32, real64
    implicit real(real64)(a-h,o-z), integer(int32) (i-n)
!     DETERMINING OF TYPE OF COEFFICIENTS B(N,K).
    if (ith) 10,18,10
10  if (mq-1) 11,11,14
11  if (lamb-1) 13,13,12
12  jtypbq = 6
    return
13  jtypbq = 5
    return
14  if (lamb-1) 16,16,15
15  if (lamb-3) 17,16,17
16  jtypbq = 1
    return
17  jtypbq = 2
    return
18  if (mq-1) 19,19,22
19  if (lamb-1) 20,20,21
20  jtypbq = 7
    return
21  jtypbq = 8
    return
22  if (lamb-1) 24,24,23
23  if (lamb-3) 25,24,25
24  jtypbq = 3
    return
25  jtypbq = 4
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine cotran

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   Changes the projectile an target Coulomb trajectoty
!
    common/dtint/dtau
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    dimension vcm(3),p(3),p1(3),p2(3),r(3),b(3),rn(3)
    if((an1+an2) == 0.)   return
    amr=an1*an2/(an1+an2)*0.940
    if(amr <= 0.)  return
    if(dtau < 0.) return
    amp=an1*.940
    amt=an2*.940
    gp=1./sqrt(1.-vpr(1)**2-vpr(2)**2-vpr(3)**2)
    gt=1./sqrt(1.-vta(1)**2-vta(2)**2-vta(3)**2)
    g1=gp*an1/(gp*an1+gt*an2)
    g2=gt*an2/(gp*an1+gt*an2)
!----> do 10  k=1,3
    do  k=1,3
       rn(k)=radp(k)-radt(k)
       r(k)=rn(k)-(vpr(k)-vta(k))*dtau
       vcm(k)=g1*vpr(k)+g2*vta(k)
       b(k)=-vcm(k)
10     p1(k)=gp*amp*vpr(k)
    end do
    call cinema(p1,b,p,ct,st,cf,sf,t,amp)
    p0=sqrt(p(1)**2+p(2)**2+p(3)**2)
    r0=sqrt(r(1)**2+r(2)**2+r(3)**2)
    r1=sqrt(rn(1)**2+rn(2)**2+rn(3)**2)
    prn=p(1)*rn(1)+p(2)*rn(2)+p(3)*rn(3)
    pl2=(prn/r1)**2-2.*amr*(potnu(r1)-potnu(r0))
    spl=1.
    if(prn < 0.)   spl=-1.
    if(pl2 < 0.)   pl2=0.
    temp=(spl*sqrt(pl2)-prn/r1)/r1
!----> do  11  k=1,3
    do   k=1,3
11     p(k)=p(k)+temp*rn(k)
    end do
    call  cinema(p,vcm,p1,ct1,st1,cf1,sf1,t1,amp)
    p(1)=-p(1)
    p(2)=-p(2)
    p(3)=-p(3)
    call  cinema(p,vcm,p2,ct2,st2,cf2,sf2,t2,amt)
    ep=amp+t1
    et=amt+t2
!----> do  12  k=1,3
    do   k=1,3
       vpr(k)=p1(k)/ep
12     vta(k)=p2(k)/et
    end do
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  forcen (r1,r2)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    Calculates Coulomb force
!
    if(abs(r1-r2) == 0.)  go  to 10
    forcen=-(potnu(r2)-potnu(r1))/(r2-r1)
    return
10  forcen=0.
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  double precision function  potnu(r)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    Calculates Coulomb potential acting between projectile and target
!
!
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    r1=0.
    if(an1 > 0.1)  s1=r0n1*an1**(1./3.)
    if(an1 > 0.1)  r1=s1*(1.+(c1*3.141592/s1)**2)**(1./3.)
    r2=0.
    if(an2 > 0.1)  s2=r0n2*an2**(1./3.)
    if(an2 > 0.1)   r2=s2*(1.+(c2*3.141592/s2)**2)**(1./3.)
    if((r1+r2) <= 0.)  go  to  14
    rp=r1
    if(r1 > r2)  rp=r2
    rt=r2
    if(r1 >= r2)  rt=r1
    if(rp <= 0.0)  go  to  14
    x=r/(rp+rt)
    b=0.00144*zn1*zn2/(rp+rt)
    a=rt/rp
    bm=1.+1./a
    c=3./5./a/a
    d=(a+2.+1./a)/4.
    x0=(a-1.)/(a+1.)
    if(x-x0)  10,10,11
10  potnu=b*(3.-c-(bm*x)**2)*bm/2.
    return
11  if(x-1.) 12,12,13
12  potnu=b*(1.-3.*d*d*((1.-x)**4)*(1.-2./15.*d*(1.-x)*(5.+x)))/x
    return
13  potnu=b/x
    return
14  potnu=0.
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  kinemq(ps,v,p,ct,st,cf,sf,t,cm)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   KINEMATIC  BLOCK
!
    dimension  ps(3),v(3),p(3)
    psv=ps(1)*v(1)+ps(2)*v(2)+ps(3)*v(3)
    es=sqrt(ps(1)**2+ps(2)**2+ps(3)**2+cm**2)
    g=1./sqrt(1.-v(1)**2-v(2)**2-v(3)**2)
!----> do  10  k=1,3
    do   k=1,3
10     p(k)=ps(k)+g*v(k)*(psv*g/(g+1.)+es)
    end do
    e=g*(es+psv)
    t=e-cm
    pm=sqrt(p(1)**2+p(2)**2+p(3)**2)
    ct=p(3)/pm
    st2=1.-ct*ct
    if(st2 <= 0.)   go  to  11
    st=sqrt(st2)
    cf=p(1)/pm/st
    sf=p(2)/pm/st
    return
11  st=0.
    cf=1.
    sf=0.
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine rxyz(r12,r0x,r0y,r0z)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    Samples impact parameter and initial positions of
!    projectile and target
!
    common /xbmax/ xbmax,ifib0
    common /bimp/ b00,bx,by
    common /hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    dimension vpr(3),vta(3)
    if(ifib0 == 0)  go to 10
    b1=rndm(-1.0_real64)
    b2=6.283185*rndm(-1.0_real64)
    b3=r12*sqrt(b1)*xbmax
    r0x=b3*cos(b2)
    r0y=b3*sin(b2)
    go  to  11
10  r0x=r12*xbmax
    r0y=0.
    b3=r0x
11  r0z=-sqrt(r12**2-r0x**2-r0y**2)
    b00=b3
    bx=r0x
    by=r0y
    call  vinit(vpr,vta,anucl1,anucl2,t0)
    gpr=1./sqrt(1.-vpr(3)**2)
    gta=1./sqrt(1.-vta(3)**2)
    r0z=-rm1/gpr-rm2/gta
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine vinit(vpr,vta,a1,a2,t0)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!    Calculates initial velocities of projectile and target
!
! CMJ 3/10/17: Added error proection
!
    dimension vpr(3),vta(3)
    common /ksyst/ksyst
    vpr(1)=0.
    vpr(2)=0.
    vta(1)=0.
    vta(2)=0.
    if(a1 < 2.1)  go to 10
    go  to  (10,11,12), ksyst
10  temp1 = t0*(t0+1.88)
    temp2 = t0+0.94
    ! call err_chk(1,'laq1.f','9239',2,temp1)
    ! call err_chk(1,'laq1.f','9239',1,temp2)
    vpr(3)=sqrt(temp1)/(temp2)
    vta(3)=1.e-6
    return
11  temp1 = t0*(t0+1.88)
    temp2 = t0+1.88
    ! call err_chk(1,'laq1.f','9246',2,temp1)
    ! call err_chk(1,'laq1.f','9246',1,temp2)
    vpr(3)=sqrt(temp1)/(temp2)
    vta(3)=-vpr(3)
    return
12  temp1 = t0*(t0+1.88)
    temp2 = t0+0.94*(1.+a1/a2)
    ! call err_chk(1,'laq1.f','9253, 9254',2,temp1)
    ! call err_chk(1,'laq1.f','9253, 9254',1,temp2)
    vpr(3)= sqrt(temp1)/(temp2)
    vta(3)=-sqrt(temp1)/(temp2)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  disob(mv,u,ind,np,mq)
! new version 11-16-95 09:59am
!     Formation of Delta from pion and nucleon
!

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    real(real64) ::   md,m2
    common /primp/ pp(3)
    common /taud3/ tau0
    common /isob3/ isob3
    common /memorylaq/ pme(9,5999),ime(5,5999)
    ind=0
    if(isob3 == 0)  return
    iatt=0
    rnd=rndm(-1.0_real64)
    if(mq == 2) then
       n1=mv+1
       n2=mv+3
    else
       n1=mv+2
       n2=mv+3
    endif
1   continue
    iatt=iatt+1
    if(iatt > 2)        return
    if((rnd <= 0.5.and.iatt == 1).or.iatt == 2)  then
       do  k=1,9
          temp=pme(k,n1)
          pme(k,n1)=pme(k,n2)
          pme(k,n2)=temp
       enddo
       do  k=1,5
          itemp=ime(k,n1)
          ime(k,n1)=ime(k,n2)
          ime(k,n2)=itemp
       enddo
    endif
    l2=mv+3-mq
    ld=mv+mq
    if((ime(4,ld)+ime(4,mv+3)).ne.1)  go  to  1
    p2=sqrt(pme(4,l2)**2+pme(5,l2)**2+pme(6,l2)**2)
    pp(1)=pme(4,mv+2)
    pp(2)=pme(5,mv+2)
    pp(3)=pme(6,mv+2)
    m2=pme(9,l2)
    e2=sqrt(p2**2+m2**2)
    md=sqrt(u*(u-2.*e2)+m2**2)
    if(md < 1.082)                   go  to  1
    call  wmd(md,t0,fmd)
    if(isob3 == 2)                    go  to  2
    drnd=rndm(-1.0_real64)
    if(drnd > fmd)                   go  to  1
2   continue
    pme(4,ld)  =-pme(4,l2)
    pme(5,ld)  =-pme(5,l2)
    pme(6,ld)  =-pme(6,l2)
    pme(9,ld)  =md
    ime(1,ld)  = ime(1,ld)  +ime(1,mv+3)
    ime(2,ld)  =0
    ime(3,ld)  =0
    ime(4,ld)  =1
    ime(5,ld)  =intg(1000.*t0)
    if(ime(5,ld) < 1)  ime(5,ld)=1
    np=np-1
    ind=mv+2*mq-1
!----> do  11 k=1,9
    do  k=1,9
       pme(k,mv+3)=pme(k,mv+2)
       if(k <= 5) ime(k,mv+3)=ime(k,mv+2)
11     continue
    end do
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  wmd(md,tau0,fmd)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!     sampling the life time of Delta resonance with mass MD
!
    real(real64) ::  md,m0
    data   m0/1.232d0/
    en=(md**2+0.940**2-0.140**2)/(2.*md)
    q=sqrt(en**2-0.940**2)
    r=q/0.140
    g=0.47*q*r*r/(1.+0.6*r*r)
    drnd=rndm(-1.0_real64)
    if(drnd < 1.d-10)  drnd=1.d-10
    tau0=-1./(5.06*g)*log(drnd)
    if(tau0 <= 0.001)  tau0=0.0011
    a=g*g/4.
    fmd=a/((md-m0)**2+a)
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine decren(lr,nur,mv,np)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    character(len=8) :: pnar,pnaj
!
!   decays of instable particle LR
!   Calls: PANUID, DECAYQ
!
    common/partcl/pptcl(9,499),nptcl,iorig(499),ident(499) &
         & ,idcay(499)
    common /memorylaq/pme(9,5999),ime(5,5999)
    common/ncasca/ncas,ncpri
    common/idpme/ idpme(5999)
    common/tlimit/tlimit
    common/actim/tint
    common/porig/iori(3,5999)
!       COMMON/IDpart/IDp1,IDp2
    dimension pin(9),iin(5)
!
    ndecay=1
    if(tint > tlimit)  ndecay=5
!
    np=0
    do   k=1,9
       pin(k)=pme(k,lr)
       if(k <= 5)  iin(k)=ime(k,lr)
    enddo
    idr=idpme(lr)
    idp1=iori(1,lr)
    idp2=iori(2,lr)
    call  panuid(idr,nur,pnar)
    if(ncas >= ncpri) then
       write(16,18) pnar,(pme(i,lr),i=4,9),(ime(k,lr),k=1,5)
       write( *,18) pnar,(pme(i,lr),i=4,9),(ime(k,lr),k=1,5)
18     format(1x,'decay:',a5,6(1x,f7.3),4i2,i10)
    endif
    nptcl=1
    pptcl(1,1)=pin(4)
    pptcl(2,1)=pin(5)
    pptcl(3,1)=pin(6)
    pptcl(4,1)=pin(8)+pin(9)
    pptcl(5,1)=pin(9)
    ident(1)=idr
    iorig(1)=0
    idcay(1)=0
!
    do  nd=1,ndecay
       if(nptcl <= 0)   go  to  116
       ndec=0
       nptcl1=nptcl
!----> do 114  i=1,nptcl1
       do  i=1,nptcl1
          idi=ident(i)
!       IF(IDI.EQ.110.OR.IDI.EQ.220.OR.IDI.EQ.230.OR.IDI.EQ.-230)
          if(idi == 110.or.idi == 230.or.idi == -230) &
               & go  to  114
          if(idi == 1230.or.idi == -1230) &
               & go  to  114
          call decayq(i,ipoint)
          if(ipoint < 0)  go  to  114
          ndec=ndec+1
          do  j=1,9
             pptcl(j,i)=pptcl(j,nptcl)
          enddo
          ident(i)=ident(nptcl)
          iorig(i)=iorig(nptcl)
          idcay(i)=idcay(nptcl)
          nptcl=nptcl-1
114       continue
       end do
       if(ndec == 0)    go  to  116
    enddo
116 continue
    if(nptcl <= 1)  then
       write(16,21) nptcl,idi,ipoint,pin(9)
       write( *,21) nptcl,idi,ipoint,pin(9)
21     format(5x,'decren: nptcl=',i2,2x,'idi=',i6,2x,'ipoint=',i3, &
            & 1x,'m=',1pe11.4)
       return
    endif
    if((mv+nptcl) > 5999)  then
       write(16,'(5x,''decren : mv+nptcl>5999'')')
       return
    endif
    do  j=1,nptcl
       m=mv+j
       idj=ident(j)
       call  panuid(idj,jp,pnaj)
       ime(1,m)=idint(1.001*charge(idj))
       ime(2,m)=0
!  !!!
       if(jp >= 55.and.jp <= 65) ime(2,m)=idj
       ime(3,m)=is(idj)
       ime(4,m)=ib(idj)
       ime(5,m)=0
!  kkg 01.02.06
       if((jp >= 8.and.jp <= 18).or.(jp >= 27.and.jp <= 36).or. &
            & (jp >= 45.and.jp <= 53)) then
          tau0= taun(jp)
          temp1 = pptcl(5,j)
          ! call err_chk(1,'laq1.f','9447',1,temp1)
          taul=tau0*pptcl(4,j)/temp1
          ime(5,m)=intg(taul*1000.)
          if(ime(5,m) == 0)  ime(5,m)=1
       endif
!
       pme(4,m)=pptcl(1,j)
       pme(5,m)=pptcl(2,j)
       pme(6,m)=pptcl(3,j)
       pme(8,m)=pptcl(4,j)-pptcl(5,j)
       pme(9,m)=pptcl(5,j)
       idpme(m)=idj
       if(ncas >= ncpri) &
            & write(16,'(1x,''====>:'',a5,6(1x,f11.4),4i3,i10)') &
            & pnaj,(pme(i,m),i=4,9),(ime(k,m),k=1,5)
    enddo
    np=nptcl
    return
  end
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! =====================================================================
! DECDVM was here until removed by CMJ (XCP-3, LANL) to remove all
!    DEAD CODE.  DECDVM was NOT called by any part of LAQGSM (or GSM)
!    Date: 09/07/2017
! Purpose: Decay of Deltas with variable mass "MD"
! =====================================================================

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! =====================================================================
! PIIZO was here until removed by CMJ (XCP-3, LANL) to remove all
!    DEAD CODE.  PIIZO was NOT called by any part of LAQGSM (or GSM)
!    Date: 09/07/2017
! Purpose: simulate Delta production in N+N==>Delta + N
! =====================================================================
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! =====================================================================
! SARNDT was here until removed by CMJ (XCP-3, LANL) to remove all
!    DEAD CODE.  SARNDT was NOT called by any part of LAQGSM (or GSM)
!    Date: 09/07/2017
! Purpose: Calculates Arndt's N+N inelastic cross sections
! =====================================================================
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  szwer(jp,u)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!    B+B==>B+Y+K   CROSS SECTION  (ZWERMANN APPROXIMATION)
!
    dimension am(8)
    data  am/2*0.494,1.115,3*1.192,0.494,0.939/
    amj=am(jp)
    amk=0.494
    szwer=0.
    if(jp <= 2) amx=1.115+0.939
    if(jp > 2) amx=0.494+0.939
    if(u <= (amx+amj))  return
    pmax=sqrt((u**2-(amx+amj)**2)*(u**2-(amx-amj)**2))/(2.*u)
    szwer=0.049*(pmax/amk)**4
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  sst(ib1,ib2,ie1,ie2,jp,t,is)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
! TOTAL PRODUCTION CROSS SECTION OF strange particle JP
!
    sst=0.
    go  to  (101,102,103,104,105,106,107),jp
!   for K+      ***********************************
101 if((ib1+ib2) == 1)  go  to  51
    if((ie1+ie2) == 2)  go  to  31
    if((ie1+ie2) == 1)  go  to  32
    if((ie1+ie2) == 0)  sst=sig(t,13)
    go  to  200
31  s1=sig(t,1)
    s2=sig(t,5)
    s3=sig(t,7)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3)
       sst=s1
       if(rs > s1.and.rs <= (s1+s2))  sst=s2
       if(rs > (s1+s2))               sst=s3
    else
! !!!
       sst=s1+s2+s3
! !!!
    endif
    go  to  200
32  s1=sig(t,2)
    s2=sig(t,9)
    s3=sig(t,11)
    if(is == 1) then

       rs=rndm(-1.0_real64)*(s1+s2+s3)
       sst=s1
       if(rs > s1.and.rs <= (s1+s2))  sst=s2
       if(rs > (s1+s2))               sst=s3
    else
! !!!
       sst=s1+s2+s3
! !!!
    endif
    go  to  200
51  if(ie1 == 1.and.ie2 == 1)   sst=sig(t,15)
    if(ie1 == -1.and.ie2 == 1)  sst=sig(t,21)
    if(ie1 == 0.and.ie2 == 1)   go  to  21
    if(ie1 == 1.and.ie2 == 0)   go  to  22
    if(ie1 == 0.and.ie2 == 0)   sst=sig(t,27)
    go  to  200
21  s1=sig(t,23)
    s2=sig(t,25)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
22  s1=sig(t,17)
    s2=sig(t,18)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
!  for K0   **************************************
102 if((ib1+ib2) == 1)   go  to  52
    if((ie1+ie2) == 2)   sst=sig(t,6)
    if((ie1+ie2) == 1)   go  to  33
    if((ie1+ie2) == 0)   go  to  34
    go  to  200
33  s1=sig(t,3)
    s2=sig(t,8)
    s3=sig(t,10)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3)
       sst=s1
       if(rs > s1.and.rs <= (s1+s2))  sst=s2
       if(rs > (s1+s2))               sst=s3
    else
! !!!
       sst=s1+s2+s3
! !!!
    endif
    go  to  200
34  s1=sig(t,4)
    s2=sig(t,14)
    s3=sig(t,12)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3)
       sst=s1
       if(rs > s1.and.rs <= (s1+s2))  sst=s2
       if(rs > (s1+s2))               sst=s3
    else
! !!!
       sst=s1+s2+s3
! !!!
    endif
    go  to  200
52  if(ie1 == -1.and.ie2 == 1)  go  to  23
    if(ie1 == 0.and.ie2 == 1)   sst=sig(t,24)
    if(ie1 == 1.and.ie2 == 0)   sst=sig(t,16)
    if(ie1 == -1.and.ie2 == 0)  sst=sig(t,22)
    if(ie1 == 0.and.ie2 == 0)   go  to  24
    go  to  200
23  s1=sig(t,19)
    s2=sig(t,20)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
24  s1=sig(t,28)
    s2=sig(t,26)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
!  for Lambda  **********************
103 if((ib1+ib2) == 1)   go  to  53
    if((ie1+ie2) == 2)   sst=sig(t,1)
    if((ie1+ie2) == 1)   go  to  25
    if((ie1+ie2) == 0)   sst=sig(t,4)
    go  to  200
25  s1=sig(t,2)
    s2=sig(t,3)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
53  if(ie1 == -1.and.ie2 == 1)  sst=sig(t,19)
    if(ie1 == 0.and.ie2 == 1)   sst=sig(t,23)
    if(ie1 == 1.and.ie2 == 0)   sst=sig(t,18)
    if(ie1 == 0.and.ie2 == 0)   sst=sig(t,26)
    go  to  200
!   for S+   ******************************************
104 if((ib1+ib2) == 1)   go  to  54
    if((ie1+ie2) == 2)   go  to  26
    if((ie1+ie2) == 1)   sst=sig(t,8)
    go  to  200
26  s1=sig(t,5)
    s2=sig(t,6)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
54  if(ie1 == 1.and.ie2 == 1)   sst=sig(t,15)
    if(ie1 == 0.and.ie2 == 1)   sst=sig(t,24)
    if(ie1 == 1.and.ie2 == 0)   sst=sig(t,16)
    go  to  200
!  for S-   *****************************************
105 if((ib1+ib2) == 1)   go  to  55
    if((ie1+ie2) == 1)   sst=sig(t,11)
    if((ie1+ie2) == 0)   go  to  27
    go  to  200
27  s1=sig(t,13)
    s2=sig(t,14)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
55  if(ie1 == -1.and.ie2 == 1)  sst=sig(t,21)
    if(ie1 == -1.and.ie2 == 0)  sst=sig(t,22)
    if(ie1 == 0.and.ie2 == 0)   sst=sig(t,27)
    go  to  200
!    for    S0   *************************************
106 if((ib1+ib2) == 1)   go  to  56
    if((ie1+ie2) == 2)   sst=sig(t,7)
    if((ie1+ie2) == 1)   go  to  28
    if((ie1+ie2) == 0)   sst=sig(t,12)
    go  to  200
28  s1=sig(t,9)
    s2=sig(t,10)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       sst=s1
       if(rs > s1)  sst=s2
    else
! !!!
       sst=s1+s2
! !!!
    endif
    go  to  200
56  if(ie1 == -1.and.ie2 == 1)  sst=sig(t,20)
    if(ie1 == 0.and.ie2 == 1)   sst=sig(t,25)
    if(ie1 == 1.and.ie2 == 0)   sst=sig(t,17)
    if(ie1 == 0.and.ie2 == 0)   sst=sig(t,28)
    go  to  200
107 continue
200 continue
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function ssl(ib1,ib2,ie1,ie2,jp,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
! TOTAL PRODUCTION CROSS SECTION OF JP+L(AS A PARTENER)
!
    ssl=0.
    go  to  (201,202,103,300,300,300,300),jp
103 ssl=sst(ib1,ib2,ie1,ie2,jp,t,0)
    go  to  300
201 if((ib1+ib2) == 1)  go  to  61
    if((ie1+ie2) == 2)   ssl=sig(t,1)
    if((ie1+ie2) == 1)   ssl=sig(t,2)
    go  to  300
61  if(ie1 == 0.and.ie2 == 1)   ssl=sig(t,23)
    if(ie1 == 1.and.ie2 == 0)   ssl=sig(t,18)
    go  to  300
202 if((ib1+ib2) == 1)   go  to  62
    if((ie1+ie2) == 1)   ssl=sig(t,3)
    if((ie1+ie2) == 0)   ssl=sig(t,4)
    go  to  300
62  if(ie1 == -1.and.ie2 == 1)  ssl=sig(t,19)
    if(ie1 == 0.and.ie2 == 0)   ssl=sig(t,26)
300 continue
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  sig(t,n)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!                        N                  REACTION
! 		 1*             p+p==>p + L  + K+
! 		 2              p+n==>n + L  + K+
! 		 3              p+n==>p + L  + K0
! 		 4              n+n==>n + L  + K0
! 		 5*             p+p==>n + S+ + K+
! 		 6*             p+p==>p + S+ + K0
! 		 7*             p+p==>p + S0 + K+
! 		 8              p+n==>n + S+ + K0
! 		 9              p+n==>n + S0 + K+
! 		10              p+n==>p + S0 + K0
! 		11              p+n==>p + S- + K+
! 		12              n+n==>n + S0 + K0
! 		13              n+n==>n + S- + K+
! 		14              n+n==>p + S- + K0
! 		15*             (pi+)+p==>S+ + K+
! 		16              (pi+)+n==>S+ + K0
! 		17              (pi+)+n==>S0 + K+
! 		18              (pi+)+n==>L  + K+
! 		19*             (pi-)+p==>L  + K0
! 		20*             (pi-)+p==>S0 + K0
! 		21*             (pi-)+p==>S- + K+
! 		22              (pi-)+n==>S- + K0
! 		23              (pi0)+p==>L  + K+
! 		24              (pi0)+p==>S+ + K0
! 		25              (pi0)+p==>S0 + K+
! 		26              (pi0)+n==>L  + K0
! 		27              (pi0)+n==>S- + K+
! 		28              (pi0)+n==>S0 + K0
!                        *    -          basic channel
!   PP   --> L K+ P
    sn1(x)=far(0.132d0,0.709d0,3.432d0,-0.564d0,x)
!   PP   --> S+ K+ N
    sn5(x)=far(0.144d0,0.709d0,3.432d0,-0.564d0,x)
!   PP   --> S+ K0 P
    sn6(x)=far(0.062d0,0.709d0,3.432d0,-0.564d0,x)
!   PP   --> S0 K+ P
    sn7(x)=far(0.058d0,0.709d0,3.432d0,-0.564d0,x)
!   PI+P --> K+ S+
    s15(x)=far(0.565d0,1.129d0,0.662d0,-1.539d0,x)
!   PI-P --> K0 L
    s19(x)=far(0.183d0,1.375d0,0.130d0,-1.213d0,x)
!   PI-P --> K0 S0
    s20(x)=far(0.098d0,1.221d0,0.150d0,-1.073d0,x)
!   PI-P --> K+ S-
    s21(x)=far(0.112d0,0.873d0,0.457d0,-1.724d0,x)
!
    tl2=t-0.759
    ts2=t-0.888
    tl3=t-1.579
    ts3=t-1.780

!      TL2=T-0.760
!      TS2=T-0.900
!      TL3=T-1.570
!      TS3=T-1.800
    go  to  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, &
         & 21,22,23,24,25,26,27,28),n
1   sig=sn1(tl3)
    return
2   sig=5.*sn1(tl3)
    return
3   go  to  2
4   sig=1.*sn1(tl3)
    return
5   sig=sn5(ts3)
    return
6   sig=sn6(ts3)
    return
7   sig=sn7(ts3)
    return
8   sig=1./3.*sn5(ts3)-2.*sn7(ts3)+4.*sn6(ts3)
    return
9   sig=4./3.*sn5(ts3)-7.*sn7(ts3)+6.*sn6(ts3)
    return
10  go  to  9
11  sig=9.*sn5(ts3)-8.*sn7(ts3)+9.*sn6(ts3)
    return
12  go  to  7
13  sig=2./3.*sn5(ts3)-4.*sn7(ts3)+3.*sn6(ts3)
14  go  to  5
15  sig=s15(ts2)
    return
16  go  to  21
17  go  to  20
18  go  to  19
19  sig=s19(tl2)
    return
20  sig=s20(ts2)
    return
21  sig=s21(ts2)
    return
22  go  to  15
23  sig=s19(tl2)/2.
    return
24  sig=2./3.*(s21(ts2)+s15(ts2)/3.-s20(ts2)/2.)
    return
25  sig=1./2.*(s21(ts2)+s15(ts2)-s20(ts2))
    return
26  go  to  23
27  go  to  20
28  sig=1./2.*(s21(ts2)+s15(ts2)-s20(ts2))
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  far(a1,a2,a3,a4,x)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!   Function form of approximated strange particle
!   production cross sections
    far=0.0
    if(x <= 0.)   return
    s=a3**2+x**2
    far=a1*(x**a2)*(s**a4)
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function spinet(t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
    real(real64) ::   mpi,mn,meta
!      (pi-)+p=>eta+n   CROSS SECTIONS  [mb]
!     T - lab. kinetic energy of pion
!
    data mpi/0.140/,mn/0.940/,meta/0.549/
    ss0=mn+meta
    e=t+mpi
    p=sqrt(t*(t+2.*mpi))
    ss=sqrt((e+mn)**2-p**2)
    spinet=0.
    if(ss <= ss0) return
    x=ss-ss0
    spinet=far(0.485d0,0.643d0,0.039d0,-0.649d0,x)
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function sppeta0(t,ksi)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
    real(real64) ::   mn,meta
!      p+p=>eta+p+p   CROSS SECTIONS  [mb]
!  De Paoli et al. Phys.Lett. B219(1989),194 [A=0.167,B=0.253]
!  New parametrization (Gudima,1996)
!     T - lab. kinetic energy of proton
    data mn/0.939/,meta/0.549/,a/0.723/,b/0.411/
    ss0=2.*mn+meta
    e=t+mn
    p=sqrt(t*(t+2.*mn))
    ss=sqrt((e+mn)**2-p**2)
    sppeta0=0.
    if(ss <= ss0) return
    x=ss-ss0
    sppeta0=a*x/(b+x**2)
!  !!!   s(p+n=p+n+eta)=2*s(p+p=p+p+eta)
    if(ksi == 2)  sppeta0=2.*sppeta0
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function sppeta(t,ksi)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
    real(real64) ::   mn,meta
!      p+p=>eta+p+p   CROSS SECTIONS  [mb]
!  Parametrization of M. Thomere et al nucl-th/0702004v1
!  sig=a*x**b, sig(mkb), x=ss-ss0 (GeV)
!     T - lab. kinetic energy of proton
    data mn/0.939/,meta/0.549/,a/0.723/,b/0.411/
    ss0=2.*mn+meta
    e=t+mn
    p=sqrt(t*(t+2.*mn))
    ss=sqrt((e+mn)**2-p**2)
    sppeta=0.
    if(ss <= ss0) return
    x=ss-ss0
    xm=x*1000.
    if(ksi == 1)  then
       if(xm < 283.0)  then
          aw=1213.8
          bw=1.50
       elseif(xm < 651.)  then
          aw=162.1
          bw=-0.08
       else
          aw=99.6
          bw=-1.24
       endif
       sppeta=aw*(x**bw)/1000. ! in mb
    else
       if(xm < 200.0)  then
          aw=25623.0
          bw=2.03
       elseif(xm < 651.)  then
          aw=324.3
          bw=-0.08
       else
          aw=199.0
          bw=-1.24
       endif
       sppeta=aw*(x**bw)/1000. ! in mb
    endif
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function spidpp(tpi)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
    real(real64) ::  mpi,md
!  cross section for pi+d=p+p; Tpi is pion lab.kin.ener.
!  Parametrization B.G.Ritche,Phys.Rev. C28(1983)926
    data a/-1.2/,b/3.5/,c/7.4/,d/5600./,er/2136./
    data  mpi/139.0/,md/1878./
    tpim=tpi*1000.
    e=sqrt((mpi+md)**2+2.*tpim*md)
    spidpp=a+b/sqrt(tpim)+c*10000./((e-er)**2+d)
    if(spidpp <= 0.)  spidpp=0.
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function ssdt(ib1,ib2,ie1,ie2,jp,t,is)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!  PRODUCTION CROSS SECTION N+N=>DELTA+JP+X
!
    ssdt=0.
    if((ib1+ib2).ne.2)   return
    if(ie1+ie2-1)  100,101,102
100 go  to  (11,12,13,14,15,16,17),jp
11  s1=sigd(t,5)
    s2=sigd(t,17)
    s3=sigd(t,20)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3)
       ssdt=s1
       if(rs > s1.and.rs <= (s1+s2))  ssdt=s2
       if(rs > (s1+s2))               ssdt=s3
    else
! !!!
       ssdt=s1+s2+s3
! !!!
    endif
    return
12  s1=sigd(t,6)
    s2=sigd(t,16)
    s3=sigd(t,18)
    s4=sigd(t,19)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3+s4)
       ssdt=s1
       if(rs > s1.and.rs <= (s1+s2))          ssdt=s2
       if(rs > (s1+s2).and.rs <= (s1+s2+s3))  ssdt=s3
       if(rs > (s1+s2+s3))                    ssdt=s4
    else
! !!!
       ssdt=s1+s2+s3+s4
! !!!
    endif
    return
13  s1=sigd(t,5)
    s2=sigd(t,6)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
14  ssdt=sigd(t,19)
    return
15  s1=sigd(t,16)
    s2=sigd(t,17)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
16  s1=sigd(t,18)
    s2=sigd(t,20)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
17  continue
    return
101 go  to  (21,22,23,24,25,26,27),jp
21  s1=sigd(t,4)
    s2=sigd(t,14)
    s3=sigd(t,15)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3)
       ssdt=s1
       if(rs > s1.and.rs <= (s1+s2))  ssdt=s2
       if(rs > (s1+s2))               ssdt=s3
    else
! !!!
       ssdt=s1+s2+s3
! !!!
    endif
    return
22  s1=sigd(t,11)
    s2=sigd(t,12)
    s3=sigd(t,13)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3)
       ssdt=s1
       if(rs > s1.and.rs <= (s1+s2))  ssdt=s2
       if(rs > (s1+s2))               ssdt=s3
    else
! !!!
       ssdt=s1+s2+s3
! !!!
    endif
    return
23  s1=sigd(t,3)
    s2=sigd(t,4)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
24  s1=sigd(t,13)
    s2=sigd(t,14)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
25  ssdt=sigd(t,11)
    return
26  s1=sigd(t,12)
    s2=sigd(t,15)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
27  return
102 go  to  (31,32,33,34,35,36,37),jp
31  s1=sigd(t,2)
    s2=sigd(t,8)
    s3=sigd(t,9)
    s4=sigd(t,10)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2+s3+s4)
       ssdt=s1
       if(rs > s1.and.rs <= (s1+s2))          ssdt=s2
       if(rs > (s1+s2).and.rs <= (s1+s2+s3))  ssdt=s3
       if(rs > (s1+s2+s3))                    ssdt=s4
    else
! !!!
       ssdt=s1+s2+s3+s4
! !!!
    endif
    return
32  s1=sigd(t,1)
    s2=sigd(t,7)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
33  s1=sigd(t,1)
    s2=sigd(t,2)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
34  s1=sigd(t,7)
    s2=sigd(t,9)
    if(is == 1) then
       rs=rndm(-1.0_real64)*(s1+s2)
       ssdt=s1
       if(rs > s1)  ssdt=s2
    else
! !!!
       ssdt=s1+s2
! !!!
    endif
    return
35  ssdt=sigd(t,10)
    return
36  ssdt=sigd(t,8)
37  return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function ssdl(ib1,ib2,ie1,ie2,jp,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!  PRODUCTION CROSS SECTION N+N=>L+DELTA+K
!
    ssdl=0.
    if((ib1+ib2).ne.2)  return
    if(ie1+ie2-1)  100,101,102
100 go  to  (11,12,13,14,14,14,14),jp
11  ssdl=sigd(t,5)
    return
12  ssdl=sigd(t,6)
    return
13  ssdl=ssdt(ib1,ib2,ie1,ie2,jp,t,0)
14  return
101 go  to  (21,22,13,14,14,14,14),jp
21  ssdl=sigd(t,4)
    return
22  ssdl=sigd(t,3)
    return
102 go  to  (31,32,13,14,14,14,14),jp
31  ssdl=sigd(t,2)
    return
32  ssdl=sigd(t,1)
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  sigd(t,n)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!
!   PP  --->  D+ L K+
    sd2(x)=far(0.132d0,0.709d0,3.432d0,-0.564d0,x)
!   PP  --->  D+ S+ K0
    sd7(x)=far(0.062d0,0.709d0,3.432d0,-0.564d0,x)
!   PP  --->  D+ S0 K+
    sd8(x)=far(0.058d0,0.709d0,3.432d0,-0.564d0,x)
!   PP  --->  D0 S+ K+
    sd9(x)=0.25*far(0.144d0,0.709d0,3.432d0,-0.564d0,x)
!
    tdl=t-2.432
    tds=t-2.870
    go  to  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n
1   sigd=3.*sd2(tdl)
    return
2   sigd=sd2(tdl)
    return
3   sigd=2.*sd2(tdl)
    return
4   sigd=2.*sd2(tdl)
    return
5   sigd=3.*sd2(tdl)
    return
6   sigd=sd2(tdl)
    return
7   sigd=sd7(tds)
    return
8   sigd=sd8(tds)
    return
9   sigd=sd9(tds)
    return
10  sigd=1./2.*(9.*sd7(tds)-6.*sd8(tds)+2.*sd9(tds))
    return
11  sigd=3.*sd9(tds)
    return
12  sigd=1./6.*(9.*sd7(tds)-6.*sd8(tds)+8.*sd9(tds))
    return
13  sigd=1./6.*(15.*sd7(tds)-6.*sd8(tds)+2.*sd9(tds))
    return
14  go  to  11
15  go  to  12
16  go  to  9
17  sigd=1./3.*(9.*sd7(tds)-12.*sd8(tds)+8.*sd9(tds))
    return
18  go  to  8
19  go  to  10
20  sigd=1./2.*(9.*sd7(tds)-12.*sd8(tds)+8.*sd9(tds))
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function  sthe(u,jp,st4)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!   S19    =  TOTAL CROSS SECTION AT T0=19.0 GeV  JP+N
!   FOR  JP=     K+  K0  L  S+  S-  S0  K-
!
    dimension  s19(7)
    data  s19/1.8,1.8,0.9,0.9,0.9,0.9,0.6/
    u4=3.325
    u19=6.13
    sthe=st4+(s19(jp)-st4)*log(u/u4)/log(u19/u4)
    return
  end
!
! * * * * * * * * * *  G: 20.05.93  * * * * * * * * * *
!
  double precision function  sikmi(pmax)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!   PRODUCTION CROSS SECTION IN  NN ---> N N K+ K-
!     SIKMI=PMAX/40.   (mb)
!
    sikmi=pmax/40.
    return
  end
! ********************************************************

! =====================================================================
! TPIABS was here until removed by CMJ (XCP-3, LANL) to remove all
!    DEAD CODE.  TPIABS was NOT called by any part of LAQGSM (or GSM)
!    Date: 09/07/2017
! Purpose: test subroutine for absorption pi+(NN) cross section
! =====================================================================

!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  double precision function crosss(jp1,iks,a,z,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!     calculates isospin averaged cross section
!
    dimension k1(7),k2(7)
!     JP    1   2   3   4   5   6  7
!          PI- PI0 PI+  P   N   d eta
    data k1/2,2,1,1,2,1,1/
    data k2/1,3,3,2,1,2,2/
    jp=jp1-1
    if(rndm(-1.0_real64) <= (z/a))  then
       ksi=k1(jp)
    else
       ksi=k2(jp)
    endif
    if(jp <= 3) then
! pi+N cross section
       st=sigmag(0,0,1,ksi,iks,t)
    elseif(jp == 4.or.jp == 5) then
! N+N cross section
       st=sigmag(0,0,2,ksi,iks,t)
    elseif(jp == 6)  then
! d+N cross section
       st=sigmad(ksi,iks,t)
    elseif(jp == 7)  then
! eta+N cross section
       st=sigeta(ksi,iks,t)
    else
       st=0.
    endif
    crosss=st
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!

! =====================================================
! SIGPI removed by CMJ (XCP-3, LANL) on 09/08/2017
!    SIGPI unused anywhere in code.
! Purpose: Calculates pion absorption cross section on pair (pn)
! =====================================================

! ********************************************************
  double precision function sigeta(ksi,iks,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!   eta + N cross sections
!
    real(real64) ::  mn,mpi,met,ms,k
    data  mn/0.939/,mpi/0.139/,met/0.549/,ms/1.535/, &
         & g0/0.150/,sm/74./,bpi/0.55/
    p0=sqrt(t*(t+2.*met))
    e0=t+met
    ss=sqrt((e0+mn)**2-p0**2)
    s=ss**2
    eet=(s+met**2-mn**2)/(2.*ss)
    epi=(s+mpi**2-mn**2)/(2.*ss)
    q=sqrt(eet**2-met**2)
    eet0=(ms**2+met**2-mn**2)/(2.*ms)
    q0=sqrt(eet0**2-met**2)
    g=g0*(q/q0)
    k=sqrt(epi**2-mpi**2)
    t0pi=(s-(mpi+mn)**2)/(2.*mn)
    sinv=spinet(t0pi)
    sdir=3.*sinv*(k/q)**2
    sres=sm*((q0/q)**2)*(g/2.)**2/((ss-ms)**2+(g/2.)**2)
    if(iks == 0)      then
! total cr.sec.
!       SIGETA=sres
       sigeta=sdir/bpi
    elseif(iks == 1)  then
! elastic cr.sec.
!       SIGETA=0.5*sres
       sigeta=(1.-bpi)/bpi*sdir
    elseif(iks == 2)  then
! eta+N-->pi+N
       sigeta=sdir
    elseif(iks == 3)  then
! abs. cr. sec.
!       SIGETA=0.5*sres
       sigeta=sdir
    else
       sigeta=sres
    endif
    if(ksi == 2)  sigeta=sigeta*2.
    return
  end
! ********************************************************
  double precision function sigmad(ksi,iks,t)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!     Approximated d+N cross sections
!
    tn=t/2.
    if(iks == 0)      then
! total Nd cross section
       sigmad=qindcr(tn,0)
       return
    elseif(iks == 1)  then
! elastic Nd cross section
       sigmad=qindcr(tn,0)-qindcr(tn,1)
       return
    elseif(iks == 3)  then
! absorption Nd cross section=inelastic cross section
       sigmad=qindcr(tn,1)
       return
    else
       sigmad=0.
       write(16,*) 'sigmad:iks=',iks
       write( *,*) 'sigmad:iks=',iks
       return
    endif
!     RETURN
  end
! ********************************************************
  double precision function qindcr(tn,is)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64)  (a-h,o-z), integer(int32) (i-n)
!
!   Interpolation of d+N cross sections
!
    dimension tnto(30),tnin(30),sndto(30),sndin(30),xx(30),yy(30)
    data tnto/ &
         & .0000,     .0010,     .0020,     .0030,     .0040, &
         & .0050,     .0060,     .0070,     .0080,     .0090, &
         & .0100,     .0200,     .0300,     .0400,     .0500, &
         & .0600,     .0700,     .0800,     .0900,     .1000, &
         & .2000,     .3000,     .4000,     .5000,     .6000, &
         & .7000,     .8000,     .9000,    1.0000,    2.0000/
    data tnin/ &
         & .0003339,  .0034,     .0035,     .00375,    .0040, &
         & .0050,     .0060,     .0070,     .0080,     .0090, &
         & .0100,     .0200,     .0300,     .0400,     .0500, &
         & .0600,     .0700,     .0800,     .0900,     .1000, &
         & .2000,     .3000,     .4000,     .5000,     .6000, &
         & .7000,     .8000,     .9000,    1.0000,    2.0000/
    data sndin/ &
         & .0   ,    1.2   ,    2.9   ,    8.0   ,   13.5   , &
         & 37.0   ,   60.2   ,   82.0   ,  102.7   ,  121.0   , &
         & 137.0   ,  205.3   ,  175.0   ,  140.0   ,  105.0   , &
         & 85.0   ,   72.0   ,   61.0   ,   54.0   ,   50.0   , &
         & 30.0   ,   26.0   ,   24.0   ,   23.0   ,   22.0   , &
         & 21.5   ,   21.2   ,   21.0   ,   20.5   ,   18.0   /
    data sndto/ &
         & 3154.0000, 2889.6817, 2572.7316, 2246.4534, 1890.2873, &
         & 1688.6910, 1522.4142, 1357.7684, 1241.8348, 1158.5485, &
         & 1067.9506,  596.6951,  379.1020,  288.4043,  219.2028, &
         & 177.6805,  150.7435,  126.9888,  113.7065,  105.0809, &
         & 64.7585,   56.1204,   57.3817,   64.8824,   70.7762, &
         & 73.8892,   76.2111,   79.0000,   80.0000,   80.0000/
!     IS=0  for total Nd cross section
!     IS=1  for inelastic(n,2n) Nd cross section
    mq=30
    x=tn
    if(is == 0)     then
       do  k=1,30
          xx(k)=tnto(k)
          yy(k)=sndto(k)
       enddo
    elseif(is == 1) then
       do  k=1,30
          xx(k)=tnin(k)
          yy(k)=sndin(k)
       enddo
    else
       qindcr=0.
       write(16,*) 'qindcr: is=',is
       write( *,*) 'qindcr: is=',is
       return
    endif
    if(x < xx(1))   then
       qindcr=yy(1)
       return
    endif
    do  j=1,mq
       k=j
       if(abs(x-xx(j)) < 1.d-10)  go  to  16
    enddo
    k=1
10  if(x-xx(k))  11,16,17
11  if(k-1)      12,12,13
12  y1=yy(1)
    x1=xx(1)
    y2=yy(2)
    x2=xx(2)
    y3=yy(3)
    x3=xx(3)
    go  to  18
13  if(k-(mq-1)) 15,14,14
14  y1=yy(mq-2)
    x1=xx(mq-2)
    y2=yy(mq-1)
    x2=xx(mq-1)
    y3=yy(mq)
    x3=xx(mq)
    go  to  18
15  y1=yy(k-1)
    x1=xx(k-1)
    y2=yy(k)
    x2=xx(k)
    y3=yy(k+1)
    x3=xx(k+1)
    go  to  18
16  qindcr= yy(k)
    return
17  k = k + 1
    if(k > mq)  go  to  14
    go  to  10
18  dt =(x2-x3)*x1**2+(x3-x1)*x2**2+(x1-x2)*x3**2
    da=(x2-x3)*y1   +(x3-x1)*y2+(x1-x2)*y3
    db=(y2-y3)*x1**2+(y3-y1)*x2**2+(y1-y2)*x3**2
    dg=(x2*y3-x3*y2)*x1**2+(x3*y1-x1*y3)*x2**2+(x1*y2-x2*y1)*x3**2
    al=da/dt
    be=db/dt
    ga=dg/dt
    qindcr = al*x**2 + be*x + ga
    if(qindcr < 0.)  qindcr=0.
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  block  data  dgtab
    implicit real*8  (a-h,o-z), integer (i-n)
    common/sabcde/eg(8),st(8),a(8),b(8),c(8),d(8),e(8)
!
! Photodesintegration cr.section and coeff. for diff.cr.sec
! from F.Partovi, Ann.Physics 27,79(1964)
!
    data  eg/ &
         & 10.,   20.,    40.,  60.,   80.,   100.,   120.,   140./
    data  st/ &
         & 1387., 558.2, 224.2, 126.6,  87.4,  66.25,  53.41,  44.53/
    data   a/ &
         & 4.623, 5.387, 6.236, 6.101, 5.651,  5.114,  4.618,  4.156/
    data   b/ &
         & 162.0, 65.50, 19.14, 7.209, 3.009,  1.025, 0.0653,-0.4096/
    data   c/ &
         & 0.3744,0.7314,0.9983,1.0010,0.9420, 0.8808, 0.8230, 0.7737/
    data   d/ &
         & 29.96, 18.25, 7.914, 4.452, 2.774,  1.812,  1.237, 0.8757/
    data   e/ &
         & -4.230,-4.223,-2.164,-1.564,-1.317,-0.9853,-0.7715,-0.6367/
  end
! ********************************************************
  block  data  wxw
    implicit real*8  (a-h,o-z), integer (i-n)
    common/wx/w(8),xw(8)
!
!     Data for Gauss 8-point integration
!
    data w/ &
         & 0.1012285352, 0.2223809958, 0.3137066364, 0.3626837730, &
         & 0.3626837730, 0.3137066364, 0.2223809958, 0.1012285352/
    data xw/ &
         & -0.9602898359,-0.7966664433,-0.5255323648,-0.1834346056, &
         & 0.1834346056, 0.5255323648, 0.7966664433, 0.9602898359/
  end

! **********************************************************************

! =====================================================================
! gengamn was here until removed by CMJ (XCP-3, LANL) to remove all
!    DEAD CODE.  gengamn was NOT called by any part of LAQGSM (or GSM)
!    Date: 09/07/2017
! Purpose: This subroutine generates a event gamma + N ==> hadrons
! =====================================================================

!  **************************************************************
  subroutine gnappr(eg,sig)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!    The function presents a calculation of approximation gamma-
!    nucleon partial (gamma + N ==> m(pi) + N, m=2-8) cross-sections
!    from experiments: S.I.Alekhin et al., CERN-HERA 87-01,
!    received from Igor Pshenichnov, Modified by KKG for LAQGSM
!
    dimension sig(7),al(7,4),tth(7), wr(7), nal(7)
    dimension pl(0:3)
    data &
         & (al(1,j),j=1,4) / 0.3317900,-0.0821800, 0.0997600,0.0001300/, &
         & (al(2,j),j=1,4) / 0.0537900, 0.1350400,-0.0317300,0.0000000/, &
         & (al(3,j),j=1,4) /-0.0009800, 0.0983250,-0.0475448,0.0000000/, &
         & (al(4,j),j=1,4) / 0.2056800,-0.0629900, 0.0000000,0.0000000/, &
         & (al(5,j),j=1,4) / 0.0619000,-0.0192100, 0.0000000,0.0000000/, &
         & (al(6,j),j=1,4) / 0.1113700,-0.0409400, 0.0000000,0.0000000/, &
         & (al(7,j),j=1,4) / 0.0336780,-0.0130250, 0.0000000,0.0000000/
    data tth/ 0.321, 0.506 , 0.727, 0.952, 1.215, 1.481, 1.788/
    data wr/ 0.7, 0.75, 0.2667, 0.4381, 0.125, 0.2755, 0.05614/
    data nal/4,3,3,2,2,2,2/
    data alfa/2.0d0/
    do  nr=1,7
       if(eg <= tth(nr))  then
          sig(nr) = 0.0d0
       else
          x = log(eg/tth(nr))
          nl = nal(nr)-1
          call  plaguer(alfa,x,pl,nl)
          f = 0.0d0
          do  m=0,nl
             f = f + al(nr,m+1)*pl(m)
          enddo
          f=f*x/exp(x/2.0d0)
          sig(nr) = f*f/wr(nr)*1000.0d0   !  [mkb]
       endif
    enddo
    return
  end
!  **************************************************************
  subroutine plaguer(al,x,pl,n)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
! The recurrent calculation of LAGERR's polinoms up to order=3
! (See for details: Handbook of mathematical functions. Ed. by
! M.Abramowitz and I.A.Stegun), received from Igor Pshenichnov
! Modified by KKG for LAQGSM
!
    dimension pl(0:3)
    data one, two/1.0d0,2.0d0/
    pl(0) = one
    if(n == 0)  then
       return
    else
       pl(1) = al + one - x
       if(n == 1)  return
       do  m=1,n-1
          rm = float(m)
          pl(m+1) = ((two*rm + al + one - x)*pl(m) - &
               & (rm + al)*pl(m-1))/(rm+one)
       enddo
    endif
    return
  end
!  **************************************************************
  function csgntot(ip,eg,emn)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!     calculate total gamma + N cross section at energy eg
!
    common /csgn3_9/ s2_9(8)
    dimension s3_9(7)
    s2_9(1) = csgn2(ip,eg,emn)
    csgntot = s2_9(1)
    call  gnappr(eg,s3_9)
    do ir=1,7
       s2_9(ir+1) = s3_9(ir)
       csgntot = csgntot + s3_9(ir)
    enddo
    return
!
  end
!  **************************************************************
  function csgn2(ip,eg,emn)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!     calculate two body gamma + N cross section at energy eg
!
    common /nchapn/ nchap(13),nchan(13)
!
    csgn2 = 0.0d0
    do ir=1,13
       if(ip == 1) then
          ich = nchap(ir)   ! gamma + p
       else
          ich = nchan(ir)   ! gamma + n
       endif
       csgn2 = csgn2 + csgnh(ich,eg,emn)
!      for channel 20 isotopicaly symmetric channal is added:
!    sig0 + K+ = sig+ + K0 (for gamma + p) or
!    sig0 + K0 = sig- + K+ (for gamma + n)
       if(ich == 20) csgn2 = csgn2 + csgnh(ich,eg,emn)
    enddo
    return
!
  end
!  **************************************************************
  function csgnh(ich,eg,emn)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!     calculate n body gamma + N cross section at energy eg
!
    common /csgnt/ wr(251),sigr(22,251)
    data zro, two /0.d0, 2.d0/
    w = sqrt((two*eg+emn)*emn)
    csgnh = 0.
    if(w < wr(1))  then
       csgnh = zro
       return
    elseif(w > wr(251))  then
       w1=wr(250)
       w2=wr(251)
       s1=sigr(ich,250)
       s2=sigr(ich,251)
       csgnh = s1 + (s2-s1)/(w2-w1)*(w-w1)
       if(csgnh < zro) csgnh = zro
       return
    endif
    do  iw=1,251
       if(abs(w-wr(iw)) < 1.0d-6)  then
          csgnh = sigr(ich,iw)
          return
       elseif(w < wr(iw).or.(w > wr(iw).and.iw == 251))  then
          w1=wr(iw-1)
          w2=wr(iw)
          s1=sigr(ich,iw-1)
          s2=sigr(ich,iw)
          csgnh = s1 + (s2-s1)/(w2-w1)*(w-w1)
          return
       endif
    enddo
    return
  end
!   ****************************************************************
  subroutine gntoh(v,u,tin1,pin,iin,pn,ipn,mv,np)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!      main subroutine to calculate gamma + N ==> hadrons
!
    common/ncasca/ncas,ncpri
    common /csgn3_9/ s2_9(8)
    dimension v(3),pin(9),iin(5),pn(9),ipn(5)
    np = 0
    rnd = rndm(-1.0_real64)
    emn = pn(9)
!      separate two body channels
    stot = csgntot(ipn(1),tin1,emn)
    sig2  = s2_9(1)
    beta2 = sig2/stot
    if(rnd <= beta2)  then
       if(ncas >= ncpri) write(*,*) ' to gnto2h: tin1,stot,sig2=', &
            & tin1,stot,sig2
       call gnto2h(v,u,tin1,sig2,pin,emn,ipn,mv,np)
       if(ncas >= ncpri) write(*,*) ' from gnto2h: np=',np
       if(np < 2)  then
!     write(*,*) ' from gnto2h(np<2): eg,stot,sig2,np=',
!    &     tin1,stot,sig2,np
       endif
       return
    else
!         channels gamma + N ==> n*pi + N, n=2,3,4,5,6,7,8
       if(ncas >= ncpri) write(*,*) ' to gntonh: tin1,stot,sig2=', &
            & tin1,stot,sig2
       call gntonh(u,tin1,ipn,mv,np)
       if(ncas >= ncpri)  write(*,*) ' from gntonh: np=',np
       if(np < 3)  then
          write(*,*) ' from gntonh(np<3): eg,stot,s2_9,np=', &
               & tin1,stot,s2_9,np
       endif
       return
    endif
    return
  end
!   ****************************************************************
  subroutine gnto2h(v,u,tin1,sig2,pin,emn,ipn,mv,np)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!      main subroutine to calculate gamma + N ==> 2 hadrons
!
    common/memorylaq/pme(9,5999),ime(5,5999)
    common /idpme/ idpme(5999)
    common /nchapn/ nchap(13),nchan(13)
    common /gbrems/ tgmin, tgmax, teqv, sxabs, ibrems
    dimension v(3),pin(9),ipn(5),sig(13),pist(3),pnst(3)
    data zro  /0.0d0/, eps /0.0001d0/
    nrep = 0
9   ssum = zro
    rnd = rndm(-1.0_real64)
    do  ic=1,13
       if(ipn(1) == 1)  then
          ich = nchap(ic)
       else
          ich = nchan(ic)
       endif
       sig(ic) =  csgnh(ich,tin1,emn)
       ssum = ssum + sig(ic)
    enddo
    ss = zro
    do  ic=1,13
       ss = ss + sig(ic)
       if(rnd <= ss/ssum)  then
          nch = ic
          go  to  10
       endif
    enddo
    write(*,*) ' gnto2h: not found 2-body channel:', &
         & 'tin1,u,emn,sig2,sig,', &
         & 'ssum=',tin1,u,emn,sig2,sig,ssum
! stop
!     go  to  9
    return
10  continue
    if(ipn(1) == 1)  then
!      target is proton
       if(nch == 1)        then  ! gamma + p ==> pi+ + n
          id1 = 120
          id2 = 1220
       elseif(nch == 2)    then  ! gamma + p ==> pi0 + p
          id1 = 110
          id2 = 1120
       elseif(nch == 3)    then  ! gamma + p ==> pi- + d++
          id1 =-120
          id2 = 1111
       elseif(nch == 4)    then  ! gamma + p ==> pi0 + d+
          id1 = 110
          id2 = 1121
       elseif(nch == 5)    then  ! gamma + p ==> pi+ + d0
          id1 = 120
          id2 = 1221
       elseif(nch == 6)    then  ! gamma + p ==> rho0 + p
          id1 = 111
          id2 = 1120
       elseif(nch == 7)    then  ! gamma + p ==> rho+ + n
          id1 = 121
          id2 = 1220
       elseif(nch == 8)    then  ! gamma + p ==> eta + p
          id1 = 220
          id2 = 1120
       elseif(nch == 9)    then  ! gamma + p ==> omeg + p
          id1 = 221
          id2 = 1120
       elseif(nch == 10)   then  ! gamma + p ==> k+ + l
          id1 = 130
          id2 = 2130
       elseif(nch == 11)   then
          if(rndm(-1.0_real64) <= 0.5d0) then ! gamma + p ==> k+ + s0
             id1 = 130
             id2 = 1230
          else                        ! gamma + p ==> k0 + s+
             id1 = 230
             id2 = 1130
          endif
       elseif(nch == 12)   then  ! gamma + p ==> etap + p
          id1 = 330
          id2 = 1120
       elseif(nch == 13)   then  ! gamma + p ==> phi + p
          id1 = 331
          id2 = 1120
       endif
       ich = nchap(nch)
    else
!      target is neutron
       if(nch == 1)        then  ! gamma + n ==> pi- + p
          id1 =-120
          id2 = 1120
       elseif(nch == 2)    then  ! gamma + n ==> pi0 + n
          id1 = 110
          id2 = 1220
       elseif(nch == 3)    then  ! gamma + n ==> pi- + d+
          id1 =-120
          id2 = 1121
       elseif(nch == 4)    then  ! gamma + n ==> pi0 + d0
          id1 = 110
          id2 = 1221
       elseif(nch == 5)    then  ! gamma + n ==> pi+ + d-
          id1 = 120
          id2 = 2221
       elseif(nch == 6)    then  ! gamma + n ==> rho- + p
          id1 =-121
          id2 = 1120
       elseif(nch == 7)    then  ! gamma + n ==> rho0 + n
          id1 = 111
          id2 = 1220
       elseif(nch == 8)    then  ! gamma + n ==> eta + n
          id1 = 220
          id2 = 1220
       elseif(nch == 9)    then  ! gamma + n ==> omeg + n
          id1 = 221
          id2 = 1220
       elseif(nch == 10)   then  ! gamma + n ==> k0 + l
          id1 = 230
          id2 = 2130
       elseif(nch == 11)   then
          if(rndm(-1.0_real64) <= 0.5d0) then  ! gamma + n ==> k0 + s0
             id1 = 230
             id2 = 1230
          else                         ! gamma + n ==> k+ + s-
             id1 = 130
             id2 = 2230
          endif
       elseif(nch == 12)   then  ! gamma + n ==> etap + n
          id1 = 330
          id2 = 1220
       elseif(nch == 13)   then  ! gamma + n ==> phi + n
          id1 = 331
          id2 = 1220
       endif
       ich = nchan(nch)
    endif
!  compute mass, strangeness, electic and baryon charge of hadrons
    am1 = amassf(id1)
    am2 = amassf(id2)
    is1 = is(id1)
    is2 = is(id2)
    ic1 = charge(id1)
    ic2 = charge(id2)
    ib1 = 0
    ib2 = 1
!   check the threshold
    if(u <= (am1+am2))  then
       if(id2 == 1111.or.id2 == 1121.or.id2 == 1221.or.id2 == 2221) &
            & then  ! change mass of delta
          am2 = 1.080d0 + rndm(-1.0_real64)*(u-am1-1.080d0)
          if(u <= (am1+am2))  go  to  9
          go  to 11
       elseif(id1 == 121.or.id1 == -121.or.id1 == 111) &
               & then  ! change mass of rho
          am1 = 0.281d0 + rndm(-1.0_real64)*(u-am2-0.281d0)
          if(u <= (am1+am2))  go  to  9
          go  to  11
       elseif(id1 == 221)       then  ! change mass of omega
          am1 = 0.660d0 + rndm(-1.0_real64)*(u-am2-0.660d0)
          if(u <= (am1+am2))  go  to  9
          go  to  11
       else
!        write(*,*) ' gnto2h:  u < am1 + am2, u, am1,am2,id1,id2,',
!    &  'nch,sig(nch),ssum,sig2=',u, am1,am2,id1,id2,nch,sig(nch),sig2
          nrep = nrep + 1
          if(nrep < 100)  then
             go  to  9
          else
             np = 0
             return
          endif
       endif
    endif
!  compute the resonance life time
11  ir1 = 0
    ir2 = 0
    if(id1 == 121)        then      !  rho+
       tau = taun(10)
       ir1  = 1
    elseif(id1 == -121)   then      !  rho-
       tau = taun(11)
       ir1 = 1
    elseif(id1 ==  111)   then      !  rho0
       tau = taun(16)
       ir1 = 1
    elseif(id1 ==  221)   then      !  omeg
       tau = taun(17)
       ir1 = 1
    elseif(id1 ==  331)   then      !  phi
       tau = taun(18)
       ir1 = 1
    elseif(id1 ==  220)   then      !  eta
       tau = taun(8)
       ir1 = 1
    elseif(id1 ==  330)   then      !  etap
       tau = taun(9)
       ir1 = 1
    endif
    if(id2 == 1111.or.id2 == 1121.or.id2 == 1221.or.id2 == 2221) then
       call wmd(am2,tau,fmd)
       ir2 = 1
    endif
    r1 = rndm(-1.0_real64)
!  select the scattering angle in c.m.s.
!       for bremstrahlung gamma call inigamn for given tin1
    if(ibrems == 1) call inigamn (tin1)
!
    cts=  cosgamm(ich,r1)
    phis= zro
!  compute the momenta of produced hadrons
    call abelq(pin,v,u,pist,pnst,cts,phis,am1,am2)
    e1 = sqrt(pist(1)**2+pist(2)**2+pist(3)**2+am1**2)
    e2 = sqrt(pnst(1)**2+pnst(2)**2+pnst(3)**2+am2**2)
    pme(1,mv+3) = zro
    pme(2,mv+3) = zro
    pme(3,mv+3) = zro
    pme(4,mv+3) = pist(1)
    pme(5,mv+3) = pist(2)
    pme(6,mv+3) = pist(3)
    pme(7,mv+3) = zro
    pme(9,mv+3) = am1
    ime(1,mv+3) = ic1
    ime(2,mv+3) = 0
    ime(3,mv+3) = is1
    ime(4,mv+3) = ib1
    if(ir1.ne.0)  then
       ime(5,mv+3) = nint(1000. *tau)
    else
       ime(5,mv+3) = 0
    endif
    pme(8,mv+3) = e1 - am1
    pme(1,mv+1) = zro
    pme(2,mv+1) = zro
    pme(3,mv+1) = zro
    pme(4,mv+1) = pnst(1)
    pme(5,mv+1) = pnst(2)
    pme(6,mv+1) = pnst(3)
    pme(7,mv+1) = zro
    pme(9,mv+1) = am2
    ime(1,mv+1) = ic2
    ime(2,mv+1) = 0
    ime(3,mv+1) = is2
    ime(4,mv+1) = ib2
    if(ir2.ne.0)  then
       ime(5,mv+1) = nint(1000. *tau)
    else
       ime(5,mv+1) = 0
    endif
    pme(8,mv+1) = e2 - am2
    idpme(mv+3) = id1
    idpme(mv+1) = id2
    np = 2
    psx = pme(4,mv+1) + pme(4,mv+3)
    psy = pme(5,mv+1) + pme(5,mv+3)
    psz = pme(6,mv+1) + pme(6,mv+3)
    es  = e1 + e2
    ies = ime(1,mv+1) + ime(1,mv+3)
    iss = ime(3,mv+1) + ime(3,mv+3)
    ibs = ime(4,mv+1) + ime(4,mv+3)
!      check the conservation law
    if(abs(psx) > eps.or.abs(psy) > eps.or.abs(psz) > eps &
         & .or.abs(es-u) > eps.or.ies.ne.ipn(1).or.iss.ne.0.or.ibs.ne.1) &
         & write(*,*) ' gnto2h: psx,psy,psz,es,u,am1,am2,ies,iss,ibs=', &
         & psx,psy,psz,es,u,am1,am2,ies,iss,ibs
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  subroutine gntonh(u,tin1,ipatne,mv,np)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
!
!     determining of secondary particles characteristics for
!     gamma + N ==> npi + N, n=2-8, interaction using GENBODL
!    last modification: 6 Dec. 2004 by KKG
!
    common /su_2wc/ su_2w3(3),su_2w4(4), &
         & su_2w5(5),su_2w6(6),su_2w7(7),su_2w8(8),su_2w9(9), &
         & isu_2p3(3,3),isu_2p4(4,4), &
         & isu_2p5(5,5),isu_2p6(6,6),isu_2p7(7,7),isu_2p8(8,8),isu_2p9(9,9), &
         & isu_2n3(3,3),isu_2n4(4,4), &
         & isu_2n5(5,5),isu_2n6(6,6),isu_2n7(7,7),isu_2n8(8,8),isu_2n9(9,9)
    common /memorylaq/ pmemo(9,5999),imemo(5,5999)
    common /csgn3_9/ s2_9(8)
    common /idpme/ idpme(5999)
    dimension ipatne(5), &
         & ama3(3),ama4(4),ama5(5),ama6(6),ama7(7),ama8(8),ama9(9), &
         & ps3(5,3),ps4(5,4),ps5(5,5),ps6(5,6),ps7(5,7),ps8(5,8),ps9(5,9), &
!
! SGM, 10/27/2011 correction sugested by KKG
!     & sig(7),sigk(7),
         & sig(7),sigk(7)
!     & w3(3),w4(4),w5(5),w6(6),w7(7),w8(8),w9(9)
!
    data emneut, emprot, empich, empi0 /0.9395656d0,0.9382723d0, &
         & 0.139568d0, 0.134973d0/
    data zro /0.0d0/, eps/0.0001d0/
!
!  determine  the resulting pion multiplicity
    s28 = zro
    do  k=1,7
       sig(k) = s2_9(k+1)
       s28 = s28 + sig(k)
       sigk(k) = s28
    enddo
    if(s28 <= zro)  then
       np = 0
       return
    endif
    rnd = rndm(-1.0_real64)
    do  k=1,7
       npi = k + 1
       if(rnd  <=  sigk(k)/s28) go  to  10
    enddo
    write(*,*) ' gntonh: npi, rnd, sigk=', npi, rnd, sigk
10  continue
    np = npi + 1
!  determine  the resulting particle types
    if(np  ==  3)  then      !  gamma + n ==> 2pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w3(n)
          if(rnd  <=  s) go  to  30
       enddo
30     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p3(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n3(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama3(k) = emprot           ! proton
             else
                ama3(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama3(k) = empi0            ! pi0
             else
                ama3(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama3,ps3,w3)
       do k=1,np
          pmemo(4,mv+k) = ps3(1,k)
          pmemo(5,mv+k) = ps3(2,k)
          pmemo(6,mv+k) = ps3(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama3(k)
       enddo
    elseif(np  ==  4)  then      !  gamma + n ==> 3pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w4(n)
          if(rnd  <=  s) go  to  40
       enddo
40     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p4(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n4(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama4(k) = emprot           ! proton
             else
                ama4(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama4(k) = empi0            ! pi0
             else
                ama4(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama4,ps4,w4)
       do k=1,np
          pmemo(4,mv+k) = ps4(1,k)
          pmemo(5,mv+k) = ps4(2,k)
          pmemo(6,mv+k) = ps4(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama4(k)
       enddo
    elseif(np  ==  5)  then      !  gamma + n ==> 4pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w5(n)
          if(rnd  <=  s) go  to  50
       enddo
50     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p5(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n5(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama5(k) = emprot           ! proton
             else
                ama5(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama5(k) = empi0            ! pi0
             else
                ama5(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama5,ps5,w5)
       do k=1,np
          pmemo(4,mv+k) = ps5(1,k)
          pmemo(5,mv+k) = ps5(2,k)
          pmemo(6,mv+k) = ps5(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama5(k)
       enddo
    elseif(np  ==  6)  then      !  gamma + n ==> 5pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w6(n)
          if(rnd  <=  s) go  to  60
       enddo
60     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p6(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n6(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama6(k) = emprot           ! proton
             else
                ama6(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama6(k) = empi0            ! pi0
             else
                ama6(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama6,ps6,w6)
       do k=1,np
          pmemo(4,mv+k) = ps6(1,k)
          pmemo(5,mv+k) = ps6(2,k)
          pmemo(6,mv+k) = ps6(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama6(k)
       enddo
    elseif(np  ==  7)  then      !  gamma + n ==> 6pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w7(n)
          if(rnd  <=  s) go  to  70
       enddo
70     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p7(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n7(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama7(k) = emprot           ! proton
             else
                ama7(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama7(k) = empi0            ! pi0
             else
                ama7(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama7,ps7,w7)
       do k=1,np
          pmemo(4,mv+k) = ps7(1,k)
          pmemo(5,mv+k) = ps7(2,k)
          pmemo(6,mv+k) = ps7(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama7(k)
       enddo
    elseif(np  ==  8)  then      !  gamma + n ==> 7pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w8(n)
          if(rnd  <=  s) go  to  80
       enddo
80     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p8(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n8(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama8(k) = emprot           ! proton
             else
                ama8(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama8(k) = empi0            ! pi0
             else
                ama8(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama8,ps8,w8)
       do k=1,np
          pmemo(4,mv+k) = ps8(1,k)
          pmemo(5,mv+k) = ps8(2,k)
          pmemo(6,mv+k) = ps8(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama8(k)
       enddo
    elseif(np  ==  9)  then      !  gamma + n ==> 8pi + n
       s = zro
       rnd = rndm(-1.0_real64)
       do  n=1,np
          nr = n
          s = s + su_2w9(n)
          if(rnd  <=  s) go  to  90
       enddo
90     continue
       do  k =1,np
          if(ipatne(1) == 1)  then
             imemo(1,mv+k) = isu_2p9(nr,k)  ! target is proton
          else
             imemo(1,mv+k) = isu_2n9(nr,k)  ! target is neutron
          endif
          if(k == 1) then
             imemo(4,mv+k) = 1            ! nucleon
             if(imemo(1,mv+k) == 1)  then
                ama9(k) = emprot           ! proton
             else
                ama9(k) = emneut           ! neutron
             endif
          else
             imemo(4,mv+k) = 0            ! pion
             if(imemo(1,mv+k) == 0)  then
                ama9(k) = empi0            ! pi0
             else
                ama9(k) = empich           ! pi+ or pi-
             endif
          endif
       enddo
       call  genbodl(np,u,ama9,ps9,w9)
       do k=1,np
          pmemo(4,mv+k) = ps9(1,k)
          pmemo(5,mv+k) = ps9(2,k)
          pmemo(6,mv+k) = ps9(3,k)
          pmemo(7,mv+k) = zro
          pmemo(9,mv+k) = ama9(k)
       enddo
    else
       write(*,*) ' gntonh: np>9 ?, np=',np
       stop
    endif
    do  k=1,np
       imemo(2,mv+k) = 0
       imemo(3,mv+k) = 0
       imemo(5,mv+k) = 0
       pmemo(8,mv+k) = sqrt(pmemo(4,mv+k)**2+pmemo(5,mv+k)**2+ &
            & pmemo(6,mv+k)**2+pmemo(9,mv+k)**2)-pmemo(9,mv+k)
    enddo
    if(imemo(1,mv+1) == 1)  then
       idpme(mv+1) = 1120        ! p
    else
       idpme(mv+1) = 1220        ! n
    endif
    do  n=2,np
       if(imemo(1,mv+n) == 0)      then
          idpme(mv+n) =  110        ! pi0
       elseif(imemo(1,mv+n) == 1)  then
          idpme(mv+n) =  120        ! pi+
       else
          idpme(mv+n) = -120        ! pi-
       endif
    enddo
    psx = zro
    psy = zro
    psz = zro
    es  = zro
    ies = 0
    iss = 0
    ibs = 0
    do  k=1,np
       psx = psx + pmemo(4,mv+k)
       psy = psy + pmemo(5,mv+k)
       psz = psz + pmemo(6,mv+k)
       es  = es  + pmemo(8,mv+k) + pmemo(9,mv+k)
       ies = ies + imemo(1,mv+k)
       iss = iss + imemo(3,mv+k)
       ibs = ibs + imemo(4,mv+k)
    enddo
!      check the conservation law
    if(abs(psx) > eps.or.abs(psy) > eps.or.abs(psz) > eps &
         & .or.abs(es-u) > eps.or.ies.ne.ipatne(1).or.iss.ne.0.or.ibs.ne.1) &
         & write(*,*) ' gntonh: np,psx,psy,psz,es,u,ies,iss,ibs=', &
         & np,psx,psy,psz,es,u,ies,iss,ibs
    return
  end
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
  subroutine genbodl(np,tecm,amass,pcms,wt)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
! Revision 1.1.1.1  1996/03/22 16:42:49  mclareni
! Received from Igor Pshenichnov, Oct.,2003
! Mofified by K.K. Gudima for LAQGSM, Nov., 2004
!   SUBROUTINE TO GENERATE N-BODY EVENT
!   ACCORDING TO FERMI LORENTZ-INVARIANT PHASE SPACE
!   ADAPTED FROM FOWL (CERN W505) SEPT. 1974 BY F. JAMES
!   EVENTS ARE GENERATED IN THEIR OWN CENTER-OF-MASS,
!   BUT MAY BE TRANSFORMED TO ANY FRAME USING LOREN4
!
!   INPUT TO SUBROUTINE IS THRU COMMON BLOCK GENIN
!             NP=NUMBER OF OUTGOING PARTICLES (.LT. 19)
!             TECM=TOTAL ENERGY IN CENTER-OF-MASS
!             AMASS(I)=MASS OF ITH OUTGOING PARTICLE
!             KGENEV=1 FOR CONSTANT CROSS SECTION
!                      2 FOR FERMI ENERGY-DEPENDANCE
!
!   OUTPUT FROM SUBROUTINE IS
!             PCMS(1,I)=X-MOMENTUM IF ITH PARTICLE
!             PCMS(2,I)=Y-MOMENTUM IF ITH PARTICLE
!             PCMS(3,I)=Z-MOMENTUM IF ITH PARTICLE
!             PCMS(4,I)=ENERGY OF ITH PARTICLE
!             PCMS(5,I)=MOMENTUM OF ITH PARTICLE
!             WT=WEIGHT OF EVENT
    dimension amass(np),pcms(5,np)
    dimension emm(18)
    dimension rno(50)
! -       PCM1 IS LINEAR EQUIV. OF PCM TO AVOID DOUBLE INDICES
    dimension em(18),pd(18),ems(18),sm(18),ffq(18),pcm1(90),pcm(5,18)
    equivalence (pcm1(1),pcm(1,1))
! FFQ(N)=PI * (TWOPI)**(N-2) / (N-2)FACTORIAL
    data ffq/0.d0,3.141592d0, 19.73921d0, 62.01255d0, 129.8788d0, &
         &   204.0131d0,256.3704d0, 268.4705d0, 240.9780d0, 189.2637d0, &
         &   132.1308d0,  83.0202d0,  47.4210d0,  24.8295d0, &
         &   12.0006d0,   5.3858d0,   2.2560d0,   0.8859d0/
    data knt,twopi,zro,one,two/0,6.2831853073d0,0.0d0,1.0d0,2.0d0/
!#if defined(cernlib_cdc)
    data kgenev/ 2 /
!#endif
!        INITIALIZATION
    nt = np
    do  i=1,np
       em(i) = amass(i)
    enddo
    knt=knt + 1
    if(knt > 1) goto 100
!     WRITE(6,1160)
!     WRITE(6,1200) NP,TECM,(AMASS(JK),JK=1,NP)
100 continue
    if(nt < 2) goto 1001
    if(nt > 18) goto 1002
    ntm1=nt-1
    ntm2=nt-2
    ntp1=nt+1
    ntnm4=3*nt - 4
    emm(1)=em(1)
    tm=zro
!----> do 200 i=1,nt
    do i=1,nt
       ems(i)=em(i)**2
       tm=tm+em(i)
200    sm(i)=tm
    end do
!        CONSTANTS DEPENDING ON TECM
    tecmtm=tecm-tm
    if(tecmtm <= zro) goto 1000
    emm(nt)=tecm
    if(kgenev > 1) goto 400
!        CONSTANT CROSS SECTION AS FUNCTION OF TECM
    emmax=tecmtm+em(1)
    emmin=zro
    wtmax=one
!----> do 350 i=2,nt
    do i=2,nt
       emmin=emmin+em(i-1)
       emmax=emmax+em(i)
350    wtmax=wtmax*pdk(emmax,emmin,em(i))
    end do
    wtmaxq=one/wtmax
    goto 455
! -      FERMI ENERGY DEPENDENCE FOR CROSS SECTION
400 wtmaxq=tecmtm**ntm2*ffq(nt) / tecm
!        CALCULATION OF WT BASED ON EFFECTIVE MASSES EMM
455 continue
! -               FILL RNO WITH 3*NT-4 RANDOM NUMBERS,
! -            OF WHICH THE FIRST NT-2 ARE ORDERED.
!----> do 457 i= 1, ntnm4
    do i= 1, ntnm4
457    rno(i)=rndm(-1.0_real64)
    end do
    if(ntm2) 900,509,460
460 continue
    call flpsor(rno,ntm2)
!----> do 508 j=2,ntm1
    do j=2,ntm1
508    emm(j)=rno(j-1)*(tecmtm)+sm(j)
    end do
509 wt=wtmaxq
    ir=ntm2
!----> do 530 i=1,ntm1
    do i=1,ntm1
       pd(i)=pdk(emm(i+1),emm(i),em(i+1))
530    wt=wt*pd(i)
    end do
! -       COMPLETE SPECIFICATION OF EVENT (RAUBOLD-LYNCH METHOD)
    pcm(1,1)=zro
    pcm(2,1)=pd(1)
    pcm(3,1)=zro
!----> do 570 i=2,nt
    do i=2,nt
       pcm(1,i)=zro
       pcm(2,i)=-pd(i-1)
       pcm(3,i)=zro
       ir=ir+1
       bang=twopi*rno(ir)
       cb=cos(bang)
       sb=sin(bang)
       ir=ir+1
       c=two*rno(ir)-one
       s=sqrt(one-c*c)
       if(i == nt) goto 1567
       esys=sqrt(pd(i)**2+emm(i)**2)
       beta=pd(i)/esys
       gama=esys/emm(i)
!----> do 568 j=1,i
       do j=1,i
          ndx=5*j - 5
          aa= pcm1(ndx+1)**2 + pcm1(ndx+2)**2 + pcm1(ndx+3)**2
          pcm1(ndx+5)=sqrt(aa)
          pcm1(ndx+4)=sqrt(aa+ems(j))
          call rotes2(c,s,cb,sb,pcm,j)
          psave=gama*(pcm(2,j)+beta*pcm(4,j))
568       pcm(2,j)=psave
       end do
       goto 570
!----> do 1568 j=1,i
1567   do j=1,i
          aa=pcm(1,j)**2 + pcm(2,j)**2 + pcm(3,j)**2
          pcm(5,j)=sqrt(aa)
          pcm(4,j)=sqrt(aa+ems(j))
          call rotes2(c,s,cb,sb,pcm,j)
1568      continue
       end do
570    continue
    end do
    do  i=1,nt
       do  j=1,5
          pcms(j,i) = pcm(j,i)
       enddo
    enddo
900 continue
    return
!          ERROR RETURNS
1000 write(*,1100)
    goto 1050
1001 write(*,1101)
    goto 1050
1002 write(*,1102)
1050 write(*,1150) knt
    write(*,1200) np,tecm,(amass(jk),jk=1,np)
    stop
1100 format(' genbod:available energy negative' )
1101 format(' genbod:less than 2 outgoing particles' )
1102 format(' genbod:more than 18 outgoing particles' )
1150 format(' genbod:above error detected in genbod at call number',i7)
1160 format(' genbod: first call to subroutine genbod' )
1200 format(' genbod: input data np=   ',i6/'   tecm=',e16.7, &
         & '  particle masses=',5e15.5/(42x,5e15.5))
  end
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  function pdk(a,b,c)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
! Received from Igor Pshenichnov, Oct.,2003
! Mofified by K.K. Gudima for LAQGSM, Nov., 2004
! -  CALLED FROM -  GENBOD
!     PDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
    a2 = a*a
    b2 = b*b
    c2 = c*c
    sqr_argument=a2 + (b2-c2)**2/a2 - 2.0d0*(b2+c2)
    if(sqr_argument > 0.0d0) then
       pdk = 0.5d0*sqrt(sqr_argument)
    else
       pdk = 0.0d0
    endif
    return
  end
!   ****************************************************************
  subroutine flpsor(rno,ntm2)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
! Received from Igor Pshenichnov, Oct.,2003
! Mofified by K.K. Gudima for LAQGSM, Nov., 2004
! -       ORDER THE FIRST NTM2 RANDOM NUMBERS
! -         TWO IS A SPECIAL CASE (FASTER)
! -       CALLED FROM GENBOD
    dimension rno(50)
    if(ntm2 - 2) 200,160,110
110 km1 = ntm2 - 1
!----> do 150 i= 1, km1
    do i= 1, km1
       iquit = 0
       ni = ntm2 - i
!----> do 140 j= 1, ni
       do j= 1, ni
          if(rno(j) - rno(j+1)) 140,140,120
120       sav = rno(j)
          rno(j) = rno(j+1)
          rno(j+1) = sav
          iquit = 1
140       continue
       end do
       if(iquit) 200,200,150
150    continue
    end do
    goto 200
160 if(rno(1) <= rno(2)) goto 200
    sav = rno(1)
    rno(1) = rno(2)
    rno(2) = sav
200 continue
    return
  end
!   ****************************************************************
  subroutine rotes2(c,s,c2,s2,pr,i)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)
! Received from Igor Pshenichnov, Oct.,2003
! Mofified by K.K. Gudima for LAQGSM, Nov., 2004
! -  CALLED FROM - GENBOD
!         THIS SUBROUTINE NOW DOES TWO ROTATIONS (XY AND XZ)
    dimension pr(50)
    k1 = 5*i - 4
    k2 = k1 + 1
    sa = pr(k1)
    sb = pr(k2)
    a      = sa*c - sb*s
    pr(k2) = sa*s + sb*c
    k2 = k2 + 1
    b = pr(k2)
    pr(k1) = a*c2 - b*s2
    pr(k2) = a*s2 + b*c2
    return
  end
!   ****************************************************************
  subroutine inigamn (egamma)

! ======================================================================
!
!    Main routine to extract ds/do for channel 1-22:
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!    Modified by KKG, Nov., 2004
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMParams, only: zro, two, twpi, degreeToRad
    use modifiedDCMData, only: theta, ctheta, xsectd, ecm, elg

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)

! ======================================================================

    real(real64), parameter :: dtheti = 1.0_real64
    real(real64), dimension(22, 19) :: s, r
    real(real64), dimension(22    ) :: st

! ======================================================================

    real(real64) :: thetai, cthetai, si, ri
    common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)

! ======================================================================

    do j = 1,181
       thetai(j) = dble(j-1)*dtheti
       cthetai(j) = cos(thetai(j)*degreeToRad)
    enddo
    do jch = 1,22
       if (egamma <= elg(jch,2))      then
          ieg1 = 2
          ieg2 = 3
       elseif (egamma >= elg(jch,50))  then
          ieg1 = 49
          ieg2 = 50
       else
          do ie = 3,50
             if (egamma >= elg(jch,ie-1) .and. egamma <= elg(jch,ie)) &
                  & then
                ieg1 = ie-1
                ieg2 = ie
                go to 10
             endif
          enddo
       endif
10     if (ieg1 < 2 .or. ieg1 > 49 .or. ieg2 < 3 .or. ieg2 > 50) &
            & then
          write (*,*) '  stop in inigamn: ieg1, ieg2 = ', ieg1, ieg2
          stop
       endif
       eg1 = elg(jch,ieg1)
       eg2 = elg(jch,ieg2)
       sint = zro
       do j = 1,19
          s1 = xsectd(jch,ieg1,j-1)
          s2 = xsectd(jch,ieg2,j-1)
          s(jch,j) = s1 + ((s2 - s1)/(eg2 - eg1))*(egamma - eg1)
          if (j >= 2) then
             dom = twpi*(ctheta(j-1) - ctheta(j))
             sint = sint + dom*(s(jch,j-1) + s(jch,j))/two
          endif
          r(jch,j) = sint
       enddo
       if(sint < zro)  sint = zro
       st(jch) = sint
       do j = 1,19
          if (sint > zro)  r(jch,j) = r(jch,j)/sint
       enddo
    end do
    do jch = 1,22
       sint = zro
       do j = 1,181
          si(jch,j) = qintxsq(thetai(j), theta, s, jch, 22, 19)
          if (j >= 2) then
             dom = twpi*(cthetai(j-1) - cthetai(j))
             sint = sint + dom*(si(jch,j-1) + si(jch,j))/two
          endif
          ri(jch,j) = sint
       end do
       if(sint < zro)  sint = zro
       si(jch,182) = sint
       do j = 1,181
          if (sint > zro)  ri(jch,j) = ri(jch,j)/sint
       enddo
    end do

    return
  end

! ======================================================================

  function qintxsq (x, th, se, l, m, n)

! ======================================================================
!
!    Interpolation of gamma+N differential cross sections.
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!
! ======================================================================


    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)

! ======================================================================

    dimension th(n), se(m,n)

    data zro, one /0.d0, 1.d0/

! ======================================================================

    do k = 1,n-1
       if (abs(x - th(k)) <= 1.0d-5)  then
          qintxsq = se(l,k)
          qintxsq = max(qintxsq, zro)
          return
       endif
    enddo
    do k = 2,n-1
       if (x < th(k))  then
          k1 = k - 1
          k2 = k
          k3 = k + 1
          go to 10
       endif
    enddo
    k1 = n - 2
    k2 = n - 1
    k3 = n
10  continue
    x1 = th(k1)
    x2 = th(k2)
    x3 = th(k3)
    y1 = se(l,k1)
    y2 = se(l,k2)
    y3 = se(l,k3)
    d  = (x2 - x3)*x1*x1 + (x3 - x1)*x2*x2 + (x1 - x2)*x3*x3
    da = (x2 - x3)*y1    + (x3 - x1)*y2    + (x1 - x2)*y3
    db = (y2 - y3)*x1*x1 + (y3 - y1)*x2*x2 + (y1 - y2)*x3*x3
    dc = (x2*y3 - x3*y2)*x1*x1 + (x3*y1 - x1*y3)*x2*x2 + &
         & (x1*y2 - x2*y1)*x3*x3
    d1 = zro
    if (d.ne.zro) d1 = one/d
    a = da*d1
    b = db*d1
    c = dc*d1
    qintxsq = a*x*x + b*x + c
    qintxsq = max(qintxsq,zro)
    return

! ======================================================================
  end
!  ***************************************************************
  function cosgamm (j0, r1)

! ======================================================================
!
!     Cosine calculation for two body gamma +N ==> hadrons reactions
!     using linear interpolation of function cos(theta)= f(r1),
!     r1=random number, j0=number of channel (see channel1.tab)
!     Energy of gamma is fixed in INIGAM.
!
!   Called by: COSEL COSEX
!
!   Written by K. K. Gudima, Nov. 2003
!   Edited by A. J. Sierk, LANL T-16, March, 2004.
!   Removed call to RNDM to calling subprogram, AJS, March, 2004.
!   Extended to channels 1-22 by KKG, Nov., 2004
! ======================================================================


    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)

! ======================================================================

    common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)

    data zro, one, pi /0.d0, 1.d0, 3.141592d0/

! ======================================================================
    jg = j0
    if (j0 < 1 .or. j0 > 22)  then
       write (*, *)  ' in cosgamm, j0 =', j0
       stop
    endif
    degrad = pi/180.0d0
    do ir = 1,181
       rrri = r1 - ri(jg,ir)
       if (abs(rrri) < 1.0d-5)  then
          cth = cos(thetai(ir)*degrad)
          go to 20
       endif
    enddo
    do ir = 2,181
       rrri = r1 - ri(jg,ir)
       if (rrri < zro) then
          ir1 = ir - 1
          ir2 = ir
          go to 10
       endif
    enddo
    ir1 = 180
    ir2 = 181
10  continue
    x1 = ri(jg,ir1)
    x2 = ri(jg,ir2)
    y1 = thetai(ir1)
    y2 = thetai(ir2)
    th = y1 + (y2 - y1)*((r1 - x1)/(x2 - x1))
    cth = cos(th*degrad)
20  temp1 = abs(cth)
    if (temp1 <= one) then
       cosgamm = cth
    else
       cosgamm = sign(one, cth)
    endif
!   kkg  12/10/04
    if(jg == 17.or.jg == 18)  cosgamm = -cosgamm
!
    return

! ======================================================================
  end
