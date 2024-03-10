
  subroutine pauli1(p1,p2,ip1,ip2,n1,n2,v,np,mv,tint,ip, results)

! ====================================================================
!
!   Check Pauli principle for collision of
!   projectile nucleon (N1,P1,IP1) with target nucleon N2
!   boost NP secondary particles into observer's system, if
!   the collision is allowed
!   Number of cascade particles MV ==> MV + NP
!
!   Edited by CMJ (3/10/17): Added error protection
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMClass, only: mDCMResults

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    type(mDCMResults), intent(inout) :: results

    ! Fragment kinetic energy
    real(real64) :: tfr1 = 0.0_real64, tfr2 = 0.0_real64

    real(real64) ::  masn
    common/ncasca/ncas,ncpri
    common/tprod/tprod(5999)
    common/porig/iori(3,5999)
    common/idn12/id1,id2
    common/nucoll/ nucoll,mvcoll(5999)
    common/cslid/clider(5999)
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2 &
         & /center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    ! This is the hole momentum, position, and excitation/kinetic energy for
    ! proj/targ. Should review physics more, but this is important!
    common /holpt/ ph1(3),ph2(3),rh1(3),rh2(3),ehol1,ehol2,tfp,tft
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common /rint/ rint
    common /idpme/ idpme(5999)
    common /sori/ sori(5999),ssor
    common /hadr1/hadr1(4,5999),hadr2(4,5999),hadi1(4),hadi2(4)
    dimension p1(9),p2(9),ip1(5),ip2(5),ps(3),v(3),pn1(3), &
         & pn2(3),p(3),pl(3),bpr(3),bta(3),rl(3)
    masn=0.940
!----> do  9  l=1,3
    do   l=1,3
       bpr(l)=-vpr(l)
9      bta(l)=-vta(l)
    end do
    tfr1=tfp
    tfr2=tft
    if(np-2)10,10,13
!----> do 11 l=1,9
10  do l=1,9
11     pmemo(l,mv+2)=pmemo(l,mv+3)
    end do
!----> do 12 l=1,5
    do l=1,5
12     imemo(l,mv+2)=imemo(l,mv+3)
    end do
    clider(mv+2)=clider(mv+3)
    idpme(mv+2)=idpme(mv+3)
    go to 16
!----> do 14 l=1,9
13  do l=1,9
       temp=pmemo(l,mv+2)
       pmemo(l,mv+2)=pmemo(l,mv+3)
14     pmemo(l,mv+3)=temp
    end do
!----> do 15 l=1,5
    do l=1,5
       itemp=imemo(l,mv+2)
       imemo(l,mv+2)=imemo(l,mv+3)
15     imemo(l,mv+3)=itemp
    end do
    temc=clider(mv+2)
    clider(mv+2)=clider(mv+3)
    clider(mv+3)=temc
    idte=idpme(mv+2)
    idpme(mv+2)=idpme(mv+3)
    idpme(mv+3)=idte
!----> do 18 l=1,2
16  do l=1,2
       if(imemo(3,mv+l).ne.0.or.imemo(5,mv+l).ne.0.or.imemo(2,mv+l).ne.0) &
            & go  to  18
       ps(1)=pmemo(4,mv+l)
       ps(2)=pmemo(5,mv+l)
       ps(3)=pmemo(6,mv+l)
       call kinemq(ps,v,pl,ctl,stl,cfl,sfl,tl,masn)
       call kinemq(pl,bpr,pn1,ct1,st1,cf1,sf1,t1,masn)
       iforb=1
       if(t1 <= tfr1)  go  to  19
       call kinemq(pl,bta,pn2,ct2,st2,cf2,sf2,t2,masn)
       iforb=2
       if(t2 <= tfr2)  go  to  19 !ton21.07.97
18     continue
    end do
    go  to  20
19  ip=0
    if(ncas >= ncpri) then
       if(iforb == 1) write(16,302) t1,tfr1
       if(iforb == 2) write(16,303) t2,tfr2
302    format(1x,'pauli1 forbided: t1=',f6.3,' < tfr1=',f6.3)
303    format(1x,'pauli1 forbided: t2=',f6.3,' < tfr2=',f6.3)
    endif
    return
20  continue
    nucoll=nucoll+1
    na1=an1+0.1
    na2=an2+0.1
24  e1=ip1(1)
    e2=ip2(1)
    j2=an2+0.1
    xc(2,n2)=xc(2,j2)
    yc(2,n2)=yc(2,j2)
    zc(2,n2)=zc(2,j2)
    iz(2,n2)=iz(2,j2)
    zn2=zn2-e2
    an2=an2-1.
    j1=an1+0.1
    xc(1,n1)=xc(1,j1)
    yc(1,n1)=yc(1,j1)
    zc(1,n1)=zc(1,j1)
    iz(1,n1)=iz(1,j1)
    zn1=zn1-e1
    an1=an1-1.
    mpa(n1)=mpa(j1)
    call  recul(1,ph1(1),ph1(2),ph1(3),rh1(1),rh1(2),rh1(3))
    call  recul(2,ph2(1),ph2(2),ph2(3),rh2(1),rh2(2),rh2(3))
    enext1=enext1+tfr1-ehol1
    enext2=enext2+tfr2-ehol2
    l=1
    np1=np
25  m=mv+l
    el=imemo(1,m)
    ql=imemo(4,m)
    p(1)=pmemo(4,m)
    p(2)=pmemo(5,m)
    p(3)=pmemo(6,m)
    call  kinemq(p,v,pl,ctl,stl,cfl,sfl,tl,pmemo(9,m))
    pmemo(4,m)=pl(1)
    pmemo(5,m)=pl(2)
    pmemo(6,m)=pl(3)
    pmemo(8,m)=tl
    call  kinemr(v,m,rl,taul)
    call  tmatur(np,v,m,p1,p2,mv,p,taul,tmat)
    pmemo(7,m)=tmat
    tprod(m)=tint
    mvcoll(m)=nucoll
    iori(1,m)=id1
    iori(2,m)=id2
    iori(3,m)=0
    do j=1,4
       hadr1(j,m)=hadi1(j)
       hadr2(j,m)=hadi2(j)
    enddo
    sori(m)=ssor
    if(imemo(5,m) == 0)    go  to  125
    taum0= dble(imemo(5,m))/1000.
    temp1 = pmemo(9,m)
!      call err_chk(1,'LaqPauli.f',"144",1,temp1)
    tauml=taum0*(1.+pmemo(8,m)/temp1)
    imemo(5,m)=intg((tauml+pmemo(7,m))*1000.)
    if(imemo(5,m) == 0)  imemo(5,m)=1
125 continue
    elm=pmemo(8,m)+pmemo(9,m)
!      PMEMO(1,M)=(P1(1)+P2(1))/2.+RL(1)-PMEMO(4,M)/ELM*TAUL
!      PMEMO(2,M)=(P1(2)+P2(2))/2.+RL(2)-PMEMO(5,M)/ELM*TAUL
!      PMEMO(3,M)=(P1(3)+P2(3))/2.+RL(3)-PMEMO(6,M)/ELM*TAUL
    dr=rint*rndm(-1.0_real64)**(1./3.) ! rndm assumed > 0; no err_chk.
    ! TODO: Search for 6.28 and replace w/ twpi
    fi=6.283185*rndm(-1.0_real64)
    ct=1.-2.*rndm(-1.0_real64)
    temp1 = 1.-ct**2
!      call err_chk(1,'LaqPauli.f',"157",2,temp1)
    st=sqrt(temp1)
    dx=dr*st*cos(fi)
    dy=dr*st*sin(fi)
    dz=dr*ct
    pmemo(1,m)=(p1(1)+p2(1))/2.+dx
    pmemo(2,m)=(p1(2)+p2(2))/2.+dy
    pmemo(3,m)=(p1(3)+p2(3))/2.+dz
    myp(m)=1
    myt(m)=1
    myy(m)=1
    if(imemo(2,m).ne.0)  myp(m)=0
    if(imemo(2,m).ne.0)  myt(m)=0
    if(imemo(2,m).ne.0)  myy(m)=0
    if(ncas >= ncpri)  write(16,301) &
         & m,(pmemo(k,m),k=1,9),(imemo(k,m),k=1,5),clider(m),idpme(m)
301 format(1x,'newp',i5,9(1x,f 8.3),2x,4i2,i15,f6.3,i5)
    if(l-np1)  36,37,37
36  l=l+1
    go to 25
37  ip=1
    mv=mv+np1
    results%projExc%numTotal=results%projExc%numTotal+1.
    results%projExc%numHoles=results%projExc%numHoles+1.
    results%targExc%numTotal=results%targExc%numTotal+1.
    results%targExc%numHoles=results%targExc%numHoles+1.
!      IF(AN1.LT.ZN1.OR.ZN1.LT.0..OR.AN2.LT.ZN2.OR.ZN2.LT.0.)
!     *write(16,300) AN1,ZN1,AN2,ZN2
300 format(2x,'pauli1',4(2x,f5.0))
    return
  end

! ====================================================
!
! ====================================================

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine pauli2(p1,p2,p3,ip1,ip2,ip3,n1,n2,n3,v,np,mv,tint, &
       & ip,obr1, results)

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMClass, only: mDCMResults

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    type(mDCMResults), intent(inout) :: results

!   Check Pauli principle for collision of
!   cascade particle  (N1,P1,IP1) with projectile nucleon N2
!   (and N3 for absorption by NN pair), and
!   boost NP secondary particles into observer's system, if
!   the collision is allowed
!   Number of cascade particles MV ==> MV + NP - 1
!   Calculate the excitation energy and recoil momentum of
!   projectile nucleus
!

    real(real64) ::  masn
    common/ncasca/ncas,ncpri
    common/tprod/tprod(5999)
    common/porig/iori(3,5999)
    common/idn12/id1,id2
    common/cslid/clider(5999)
    common/nucoll/ nucoll,mvcoll(5999)
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2 &
         & /center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common /rint/ rint
    ! tf* is kinetic energy of proj/target;
    common /holpt/ ph1(3),ph2(3),rh1(3),rh2(3),ehol1,ehol2,tfp,tft
    common /idpme/ idpme(5999)
    common /sori/ sori(5999),ssor
    common /hadr1/hadr1(4,5999),hadr2(4,5999),hadi1(4),hadi2(4)
    common /barpot/ pot
    common/tefabs/ tefabs,ehol3,tfr3
    common /parinc/ pinc(3),einc
    dimension b(3),p1(9),p2(9),ip1(5),ip2(5),ps(3),v(3),pn1(3), &
         & pn2(3),p(3),pl(3),p3(9),ip3(5),rl(3)
    masn=0.940
    an10=an1
    zn10=zn1
    enex10=enext1
    esum=0.
    b(1)=vpr(1)
    b(2)=vpr(2)
    b(3)=vpr(3)
    temp1 = 1.-vpr(1)**2-vpr(2)**2-vpr(3)**2
!      call err_chk(1,'LaqPauli.f',"244",2,temp1)
    temp2 = sqrt(temp1)
!      call err_chk(1,'LaqPauli.f',"246",1,temp2)
    gpr=1./temp2
    temp1 = p2(1)**2+p2(2)**2+p2(3)**2
!      call err_chk(1,'LaqPauli.f',"249",2,temp1)
    r1=sqrt(temp1)
    iabs=0
    tfr1=tfp
    if(np-2)32,10,13
!----> do 11 l=1,9
10  do l=1,9
11     pmemo(l,mv+2)=pmemo(l,mv+3)
    end do
!----> do 12 l=1,5
    do l=1,5
12     imemo(l,mv+2)=imemo(l,mv+3)
    end do
    if(ip1(4) == 0.and.(imemo(4,mv+1)+imemo(4,mv+2)) == 2)  iabs=1
    clider(mv+2)=clider(mv+3)
    idpme(mv+2)=idpme(mv+3)
    go to 16
!----> do 14 l=1,9
13  do l=1,9
       temp=pmemo(l,mv+2)
       pmemo(l,mv+2)=pmemo(l,mv+3)
14     pmemo(l,mv+3)=temp
    end do
!----> do 15 l=1,5
    do l=1,5
       itemp=imemo(l,mv+2)
       imemo(l,mv+2)=imemo(l,mv+3)
15     imemo(l,mv+3)=itemp
    end do
    temc=clider(mv+2)
    clider(mv+2)=clider(mv+3)
    clider(mv+3)=temc
    idte=idpme(mv+2)
    idpme(mv+2)=idpme(mv+3)
    idpme(mv+3)=idte
16  if(imemo(4,mv+1).ne.1)     go  to  17
    if(imemo(3,mv+1).ne.0.or.imemo(5,mv+1).ne.0.or.imemo(2,mv+1).ne.0) &
         & go  to  17
    ps(1)=pmemo(4,mv+1)
    ps(2)=pmemo(5,mv+1)
    ps(3)=pmemo(6,mv+1)
    call cinema(ps,v,pn1,ct1,st1,cf1,sf1,t1,masn)
    tn=t1
    if(t1-tfr1)   31,31,17 !ton21.07.97
17  if(imemo(4,mv+2).ne.1)  go  to  32
18  if(imemo(3,mv+2).ne.0.or.imemo(5,mv+2).ne.0.or.imemo(2,mv+2).ne.0) &
         & go  to  32
    ps(1)=pmemo(4,mv+2)
    ps(2)=pmemo(5,mv+2)
    ps(3)=pmemo(6,mv+2)
    call cinema(ps,v,pn2,ct2,st2,cf2,sf2,t2,masn)
    tn=t2
    if(t2-tfr1) 31,31,32 !ton21.07.97
31  ip=0
    an1=an10
    zn1=zn10
    enext1=enex10
    if(ncas >= ncpri) write(16,302) tn,tfr1
302 format(1x,'pauli2 forbided: tn=',f6.3,' < tfr1=',f6.3)
    return
32  continue
    e2=ip2(1)
    if((an1-1.) < (zn1-e2).or.(zn1-e2) < 0.) go  to  31
    if(iabs == 1.and.np.ne.1)  then
       e3=ip3(1)
       if((an1-1.) < (zn1-e3).or.(zn1-e3) < 0.) go to 31
    endif
    nucoll=nucoll+1
    ibar=0
    if(ip1(4) == -1.and.ip1(3) == 0.and.ip1(5) == 0) ibar=1
    j1=an1+0.1
    if(iabs == 1.and.j1 == n3)   j1=j1-1
    xc(1,n2)=xc(1,j1)
    yc(1,n2)=yc(1,j1)
    zc(1,n2)=zc(1,j1)
    iz(1,n2)=iz(1,j1)
    mpa(n2)=mpa(j1)
    zn1=zn1-e2
    an1=an1-1.
    call  recul(1,ph1(1),ph1(2),ph1(3),rh1(1),rh1(2),rh1(3))
!     ENEXT1=ENEXT1+TFR1-EHOL1
    if(np == 1)  go  to  51
    if(iabs.ne.1)  go  to  51
! *************************************************
! 16.01.97
!     ENEXT1=ENEXT1+TFR3-EHOL3-(TEFABS-EHOL1-EHOL3)-VPI
    results%projExc%numTotal=results%projExc%numTotal+1.
    results%projExc%numHoles=results%projExc%numHoles+1.
    e3=ip3(1)
    call recul(1,p3(4),p3(5),p3(6),p3(1),p3(2),p3(3))
    j1=an1+0.1
    if(j1 < n3)   j1=n3
    xc(1,n3)=xc(1,j1)
    yc(1,n3)=yc(1,j1)
    zc(1,n3)=zc(1,j1)
    iz(1,n3)=iz(1,j1)
    mpa(n3)=mpa(j1)
    zn1=zn1-e3
    an1=an1-1.
51  continue
    l=1
    np1=np
    vpi1=vpi
    nabsn=0
    iflag=0
33  m=mv+l
    el=imemo(1,m)
    ql=imemo(4,m)
    p(1)=pmemo(4,m)
    p(2)=pmemo(5,m)
    p(3)=pmemo(6,m)
    call cinema(p,v,pl,ctl,stl,cfl,sfl,tl,pmemo(9,m))
!
    if(ncas >= ncpri)  write(16,299) &
         & m,pl,tl+pmemo(9,m),(imemo(k,m),k=1,5),idpme(m)
299 format(1x,'newp2 in cms:',i5,4(1x,f8.3),2x,4i2,i15,i5)
!
    tl0=tl
! 16.01.97
    if(imemo(4,m) < 0.or.imemo(3,m).ne.0.or.imemo(2,m).ne.0) &
         & go  to  46
    if(imemo(4,m).ne.0)                go  to  105
    if(imemo(5,m).ne.0)                go  to  46
    if(imemo(1,m))  103,105,104
103 if((zn1*(an1-1.)/2.) < 1.1)       go  to  46
    go  to  105
104 if(((an1-1.)*(an1-zn1)/2.) < 1.1) go  to  46
105 continue
    tl0=tl-ql*(tfr1+eps1)-(1.-ql)*vpi
    if(iabs == 1.and.iflag == 1) tl0=tl-ql*(tfr3+eps1)-(1.-ql)*vpi
    iflag=iflag+1
!   kkg  17.06.02
    if(imemo(5,m).ne.0)  then
       if(tl0 <= 0) then
          if(ncas >= ncpri) write(16,303) tl,tfr1+eps1
303       format(1x,'pauli2 forbided: tdelta=',f6.3,'< tfr1+eps1=',f6.3)
          go  to  31
! NVM          return
       else
          go  to  46
       endif
    endif
!
35  continue
    temp1 = znucl1
!      call err_chk(1,'LaqPauli.f',"386",1,temp1)
    cut=eps1*ql+obr1*el*zn1/temp1           !  04/29/04
!     CUT=0.001*QL+OBR1*EL*ZN1/ZNUCL1         !  06/26/02, 12/13/04
    if(el)  101,102,102
101 cut=0.
102 continue
    if(tl0)  36,100,100
100 continue
    if(tl0-cut) 36,36,46
36  continue
    if((an1+ql) < (zn1+el).or.(zn1+el) < 0.)  then
       tl0=tl
       go  to  46
    endif
!     ENEXT1=ENEXT1+TL0+EPS1*QL+(.14+VPI)*(1.-QL)
    call  recul(1,pl(1),pl(2),pl(3),p1(1),p1(2),p1(3))
    nabsn=nabsn+imemo(4,m)
    if((an1+ql) < (zn1+el).or.(zn1+el) < 0.) go  to  31
    an1=an1+ql
    zn1=zn1+el
    if(imemo(4,m)-1) 141,37,141
37  if(nabsn-1) 38,38,39
38  j1=an1+.1
    xc(1,j1)=p2(1)
    yc(1,j1)=p2(2)
    zc(1,j1)=p2(3)
    go to 40
141 j1=an1+0.1
!----> do  108  k=1,j1
    do   k=1,j1
       if(imemo(1,m))  106,41,107
106    if(iz(1,k).ne.1)   go  to  108
       iz(1,k)=0
       go  to  41
107    if(iz(1,k).ne.0)   go  to  108
       iz(1,k)=1
       go  to  41
108    continue
    end do
    go  to  41
39  j1=an1+.1
    xc(1,j1)=p1(1)
    yc(1,j1)=p1(2)
    zc(1,j1)=p1(3)
40  iz(1,j1)=imemo(1,m)
    mpa(j1)=1
    results%projExc%numTotal=results%projExc%numTotal+1.
    results%projExc%numProtons=results%projExc%numProtons+el
41  if(l-np1) 42,45,45
!----> do 43 k=1,9
42  do k=1,9
43     pmemo(k,m)=pmemo(k,mv+np1)
    end do
!----> do 44 k=1,5
    do k=1,5
44     imemo(k,m)=imemo(k,mv+np1)
    end do
    clider(m)=clider(mv+np1)
    idpme(m)=idpme(mv+np1)
    np1=np1-1
    go to 33
45  np1=np1-1
    go to 48
46  continue
    if(imemo(4,m) == -1.and.imemo(3,m) == 0.and.imemo(5,m) == 0.and. &
         & ibar == 1)       then
       tl0=tl-pot
       if(tl0 <= 0.)  then
!        ENEXT1=ENEXT1+TL+0.940
          an1=an1-1.
          zn1=zn1+el
          call  recul(1,pl(1),pl(2),pl(3),p1(1),p1(2),p1(3))
          go  to  41
       endif
    endif
    temp1 = tl0*(tl0+2.*pmemo(9,m))
!      call err_chk(1,'LaqPauli.f',"456",2,temp1)
    pm=sqrt(temp1)
    pn2(1)=pm*stl*cfl
    pn2(2)=pm*stl*sfl
    pn2(3)=pm*ctl
!
    esum=esum+tl0+pmemo(9,m)
!
    call cinema(pn2,b,ps,ct,st,cf,sf,ttl,pmemo(9,m))
    pmemo(4,m)=ps(1)
    pmemo(5,m)=ps(2)
    pmemo(6,m)=ps(3)
    pmemo(8,m)=ttl
    call  kinemr(b,m,rl,taul)
    call  tmatur(np,v,m,p1,p2,mv,p,taul,tmat)
    pmemo(7,m)=tmat
    tprod(m)=tint
    mvcoll(m)=nucoll
    if(idpme(n1) == idpme(m).and.np == 2) then
! rescattering (elastic or charge exchange)  of the N1 particle
       iori(1,m)=iori(1,n1)
       iori(2,m)=iori(2,n1)
       iori(3,m)=iori(3,n1)+1
       do j=1,4
          hadr1(j,m)=hadr1(j,n1)
          hadr2(j,m)=hadr2(j,n1)
       enddo
       sori(m)=sori(n1)
    else
! production of a new particle M
       iori(1,m)=id1
       iori(2,m)=id2
       iori(3,m)=0
       do j=1,4
          hadr1(j,m)=hadi1(j)
          hadr2(j,m)=hadi2(j)
       enddo
       sori(m)=ssor
    endif
    myp(m)=1
    myt(m)=1
    myy(m)=1
    if(imemo(2,m).ne.0)  myp(m)=0
    if(imemo(2,m).ne.0)  myt(m)=0
    if(imemo(2,m).ne.0)  myy(m)=0
    if(imemo(5,m) == 0)  go  to  146
    taum0= dble(imemo(5,m))/1000.
    temp1 = pmemo(9,m)
!      call err_chk(1,'LaqPauli.f',"504",1,temp1)
    tauml=taum0*(1.+pmemo(8,m)/temp1)
    imemo(5,m)=intg((tauml+pmemo(7,m))*1000.)
    if(imemo(5,m) == 0)  imemo(5,m)=1
146 continue
    x=(p1(1)+p2(1))/2.
    y=(p1(2)+p2(2))/2.
    z=(p1(3)+p2(3))/2.
    vrk=x*vpr(1)+y*vpr(2)+z*vpr(3)
    dr=rint*rndm(-1.0_real64)**(1./3.)
    fi=6.283185*rndm(-1.0_real64)
    ct=1.-2.*rndm(-1.0_real64)
    temp1 = 1.-ct**2
!      call err_chk(1,'LaqPauli.f',"517",2,temp1)
    st=sqrt(temp1)
    dx=dr*st*cos(fi)
    dy=dr*st*sin(fi)
    dz=dr*ct
    temp2 = gpr+1.
!      call err_chk(1,'LaqPauli.f',"523, 524, 525",1,temp2)
    pmemo(1,m)=radp(1)+x-vpr(1)*vrk*gpr/temp2+dx
    pmemo(2,m)=radp(2)+y-vpr(2)*vrk*gpr/temp2+dy
    pmemo(3,m)=radp(3)+z-vpr(3)*vrk*gpr/temp2+dz
    if(ncas >= ncpri)  write(16,301) &
         & m,(pmemo(k,m),k=1,9),(imemo(k,m),k=1,5),clider(m),idpme(m)
301 format(1x,'newp',i5,9(1x,f 8.3),2x,4i2,i15,f6.3,i5)
    if(l-np1) 47,48,48
47  l=l+1
    go to 33
48  ip=1
    mv=mv+np1-1
!----> do 53 k=1,9
52  do k=1,9
53     pmemo(k,n1)=pmemo(k,mv+1)
    end do
!----> do 54 k=1,5
    do k=1,5
54     imemo(k,n1)=imemo(k,mv+1)
    end do
    myp(n1)=myp(mv+1)
    myt(n1)=myt(mv+1)
    myy(n1)=myy(mv+1)
    tprod(n1)=tprod(mv+1)
    mvcoll(n1)=mvcoll(mv+1)
    iori(1,n1)=iori(1,mv+1)
    iori(2,n1)=iori(2,mv+1)
    iori(3,n1)=iori(3,mv+1)
    do j=1,4
       hadr1(j,1)=hadr1(j,mv+1)
       hadr2(j,1)=hadr2(j,mv+1)
    enddo
    sori(n1)=sori(mv+1)
    clider(n1)=clider(mv+1)
    idpme(n1)=idpme(mv+1)
    call eraij(n1)
    if(np1 == 0)              call repij(n1,mv+1)
    results%projExc%numTotal=results%projExc%numTotal+1.
    results%projExc%numHoles=results%projExc%numHoles+1.
!
    enext1=enex10+einc-esum+(an10-an1)*(0.940-eps1)
!
    return
  end

! =========================================
!
! =========================================

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine pauli3(p1,p2,p3,ip1,ip2,ip3,n1,n2,n3,v,np,mv,tint, &
       & ip,obr2, results)

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMClass, only: mDCMResults

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    type(mDCMResults), intent(inout) :: results

!   Check Pauli principle for collision of
!   cascade particle  (N1,P1,IP1) with target nucleon N2
!   (and N3 for absorption by NN pair), and
!   boost NP secondary particles into observer's system, if
!   the collision is allowed
!   Number of cascade particles MV ==> MV + NP - 1
!   Calculate the excitation energy and recoil momentum of
!   target nucleus
!
!
    real(real64) ::  masn
    common/ncasca/ncas,ncpri
    common/tprod/tprod(5999)
    common/porig/iori(3,5999)
    common/idn12/id1,id2
    common/cslid/clider(5999)
    common/nucoll/ nucoll,mvcoll(5999)
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2 &
         & /center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    common/memorylaq/pmemo(9,5999),imemo(5,5999)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common /rint/ rint
    common /holpt/ ph1(3),ph2(3),rh1(3),rh2(3),ehol1,ehol2,tfp,tft
    common /idpme/ idpme(5999)
    common /sori/ sori(5999),ssor
    common /hadr1/hadr1(4,5999),hadr2(4,5999),hadi1(4),hadi2(4)
    common /barpot/ pot
    common/tefabs/ tefabs,ehol3,tfr3
    common /parinc/ pinc(3),einc
    dimension b(3),p1(9),p2(9),ip1(5),ip2(5),ps(3),v(3),pn1(3), &
         & pn2(3),p(3),pt(3),pl(3),p3(9),ip3(5),rl(3)
    masn=0.940
    an20=an2
    zn20=zn2
    enex20=enext2
    esum=0.
    b(1)=vta(1)
    b(2)=vta(2)
    b(3)=vta(3)
    temp1 = 1.-vta(1)**2-vta(2)**2-vta(3)**2
!      call err_chk(1,'LaqPauli.f',"619",2,temp1)
    temp2 = sqrt(temp1)
!      call err_chk(1,'LaqPauli.f',"621",1,temp2)
    gta=1./temp2
    temp3 = p2(1)**2+p2(2)**2+p2(3)**2
!      call err_chk(1,'LaqPauli.f',"624",2,temp3)
    r2=sqrt(temp3)
    iabs=0
    tfr2=tft
    if(np-2)32,10,13
!----> do 11 l=1,9
10  do l=1,9
11     pmemo(l,mv+2)=pmemo(l,mv+3)
    end do
!----> do 12 l=1,5
    do l=1,5
12     imemo(l,mv+2)=imemo(l,mv+3)
    end do
    if(ip1(4) == 0.and.(imemo(4,mv+1)+imemo(4,mv+2)) == 2)  iabs=1
    clider(mv+2)=clider(mv+3)
    idpme(mv+2)=idpme(mv+3)
    go to 16
!----> do 14 l=1,9
13  do l=1,9
       temp=pmemo(l,mv+2)
       pmemo(l,mv+2)=pmemo(l,mv+3)
14     pmemo(l,mv+3)=temp
    end do
!----> do 15 l=1,5
    do l=1,5
       itemp=imemo(l,mv+2)
       imemo(l,mv+2)=imemo(l,mv+3)
15     imemo(l,mv+3)=itemp
    end do
    temc=clider(mv+2)
    clider(mv+2)=clider(mv+3)
    clider(mv+3)=temc
    idte=idpme(mv+2)
    idpme(mv+2)=idpme(mv+3)
    idpme(mv+3)=idte
16  if(imemo(4,mv+1).ne.1)   go  to  17
    if(imemo(3,mv+1).ne.0.or.imemo(5,mv+1).ne.0.or.imemo(2,mv+1).ne.0) &
         & go  to  17
    ps(1)=pmemo(4,mv+1)
    ps(2)=pmemo(5,mv+1)
    ps(3)=pmemo(6,mv+1)
    call cinema(ps,v,pn1,ct1,st1,cf1,sf1,tn,masn)
    if(tn-tfr2) 31,31,17 !ton21.07.97
17  if(imemo(4,mv+2).ne.1)  go  to  32
18  if(imemo(3,mv+2).ne.0.or.imemo(5,mv+2).ne.0.or.imemo(2,mv+2).ne.0) &
         & go  to  32
    ps(1)=pmemo(4,mv+2)
    ps(2)=pmemo(5,mv+2)
    ps(3)=pmemo(6,mv+2)
    call cinema(ps,v,pn2,ct2,st2,cf2,sf2,tn,masn)
    if(tn-tfr2) 31,31,32 !ton21.07.97
31  ip=0
    an2=an20
    zn2=zn20
    enext2=enex20
    if(ncas >= ncpri) write(16,302) tn,tfr2
302 format(1x,'pauli3 forbided: tn=',f6.3,' < tfr2=',f6.3)
    return
32  continue
    ibar=0
    if(ip1(4) == -1.and.ip1(3) == 0.and.ip1(5) == 0) ibar=1
    e2=ip2(1)
    if((an2-1.) < (zn2-e2).or.(zn2-e2) < 0.) go  to  31
    if(iabs == 1.and.np.ne.1)  then
       e3=ip3(1)
       if((an2-1.) < (zn2-e3).or.(zn2-e3) < 0.) go to 31
    endif
    nucoll=nucoll+1
    j2=an2+0.1
    if(iabs == 1.and.j2 == n3)   j2=j2-1
25  xc(2,n2)=xc(2,j2)
    yc(2,n2)=yc(2,j2)
    zc(2,n2)=zc(2,j2)
    iz(2,n2)=iz(2,j2)
    zn2=zn2-e2
    an2=an2-1.
    call  recul(2,ph2(1),ph2(2),ph2(3),rh2(1),rh2(2),rh2(3))
!     ENEXT2=ENEXT2+TFR2-EHOL2
    if(np == 1) go  to  51
    if(iabs.ne.1)  go  to  51
! *************************************************
! 16.01.97
!     ENEXT2=ENEXT2+TFR3-EHOL3-(TEFABS-EHOL2-EHOL3)-VPI
    results%targExc%numTotal=results%targExc%numTotal+1.
    results%targExc%numHoles=results%targExc%numHoles+1.
    call recul(2,p3(4),p3(5),p3(6),p3(1),p3(2),p3(3))
    e3=ip3(1)
    j2=an2+0.1
    if(j2 < n3)   j2=n3
132 xc(2,n3)=xc(2,j2)
    yc(2,n3)=yc(2,j2)
    zc(2,n3)=zc(2,j2)
    iz(2,n3)=iz(2,j2)
    zn2=zn2-e3
    an2=an2-1.
51  continue
    l=1
    np1=np
    vpi1=vpi
    nabsn=0
    iflag=0
33  m=mv+l
    el=imemo(1,m)
    ql=imemo(4,m)
    p(1)=pmemo(4,m)
    p(2)=pmemo(5,m)
    p(3)=pmemo(6,m)
    call cinema(p,v,pl,ctl,stl,cfl,sfl,tl,pmemo(9,m))
!
    if(ncas >= ncpri)  write(16,299) &
         & m,pl,tl+pmemo(9,m),(imemo(k,m),k=1,5),idpme(m)
299 format(1x,'newp3 in cms:',i5,4(1x,f8.3),2x,4i2,i15,i5)
!
    tl0=tl
! 16.01.97
    if(imemo(4,m) < 0.or.imemo(3,m).ne.0.or.imemo(2,m).ne.0) &
         & go  to  46
    if(imemo(4,m).ne.0)                go  to  105
    if(imemo(5,m).ne.0)                go  to  46
    if(imemo(1,m))  103,105,104
103 if((zn2*(an2-1.)/2.) < 1.1)       go  to  46
    go  to  105
104 if(((an2-1.)*(an2-zn2)/2.) < 1.1) go  to  46
105 continue
    tl0=tl-ql*(tfr2+eps2)-(1.-ql)*vpi
    if(iabs == 1.and.iflag == 1) &
         & tl0=tl-ql*(tfr3+eps2)-(1.-ql)*vpi
    iflag=iflag+1
!   kkg  17.06.02
    if(imemo(5,m).ne.0)  then
       if(tl0 <= 0) then
          if(ncas >= ncpri) write(16,303) tl,tfr2+eps2
303       format(1x,'pauli3 forbided: tdelta=',f6.3,'< tfr2+eps2=',f6.3)
          go  to  31
! NVM:          return
       else
          go  to  46
       endif
    endif
!
35  continue
    temp1 = znucl2
!      call err_chk(1,'LaqPauli.f',"758",1,temp1)
    cut=eps2*ql+obr2*el*zn2/temp1         !  04/29/04
!     CUT=0.001*QL+OBR2*EL*ZN2/ZNUCL2      ! 06/26/02, 12/13/04
    if(el)  101,102,102
101 cut=0.
102 continue
    if(tl0)  36,100,100
100 continue
    if(tl0-cut) 36,36,46
36  continue
    if((an2+ql) < (zn2+el).or.(zn2+el) < 0.)  then
       tl0=tl
       go  to  46
    endif
!     ENEXT2=ENEXT2+TL0+EPS2*QL+(.14+VPI)*(1.-QL)
    call recul(2,pl(1),pl(2),pl(3),p1(1),p1(2),p1(3))
    nabsn=nabsn+imemo(4,m)
    if((an2+ql) < (zn2+el).or.(zn2+el) < 0.) go  to  31
    an2=an2+ql
    zn2=zn2+el
    if(imemo(4,m)-1)  141,37,141
37  if(nabsn-1) 38,38,39
38  j2=an2+.1
    xc(2,j2)=p2(1)
    yc(2,j2)=p2(2)
    zc(2,j2)=p2(3)
    go to 40
141 j2=an2+0.1
!----> do  108  k=1,j2
    do   k=1,j2
       if(imemo(1,m))  106,41,107
106    if(iz(2,k).ne.1)   go  to  108
       iz(2,k)=0
       go  to  41
107    if(iz(2,k).ne.0)   go  to  108
       iz(2,k)=1
       go  to  41
108    continue
    end do
    go  to  41
39  j2=an2+.1
    xc(2,j2)=p1(1)
    yc(2,j2)=p1(2)
    zc(2,j2)=p1(3)
40  iz(2,j2)=imemo(1,m)
    results%targExc%numTotal=results%targExc%numTotal+1.
    results%targexc%numProtons=results%targexc%numProtons+el
41  if(l-np1) 42,45,45
!----> do 43 k=1,9
42  do k=1,9
43     pmemo(k,m)=pmemo(k,mv+np1)
    end do
!----> do 44 k=1,5
    do k=1,5
44     imemo(k,m)=imemo(k,mv+np1)
    end do
    clider(m)=clider(mv+np1)
    idpme(m)=idpme(mv+np1)
    np1=np1-1
    go to 33
45  np1=np1-1
    go to 48
46  continue
    if(imemo(4,m) == -1.and.imemo(3,m) == 0.and.imemo(5,m) == 0.and. &
         & ibar == 1)       then
       tl0=tl-pot
       if(tl0 <= 0.)  then
!        ENEXT2=ENEXT2+TL+0.940
          an2=an2-1.
          zn2=zn2+el
          call  recul(2,pl(1),pl(2),pl(3),p1(1),p1(2),p1(3))
          go  to  41
       endif
    endif
    temp1 = tl0*(tl0+2.*pmemo(9,m))
!      call err_chk(1,'LaqPauli.f',"827",2,temp1)
    pm0=sqrt(temp1)
    temp1 = tl*(tl+2.*pmemo(9,m))
!      call err_chk(1,'LaqPauli.f',"830",2,temp1)
    temp2 = sqrt(temp1)
!      call err_chk(1,'LaqPauli.f',"832",1,temp2)
    dpm=pm0/temp2
    pt(1)=pl(1)*dpm
    pt(2)=pl(2)*dpm
    pt(3)=pl(3)*dpm
!
    esum=esum+tl0+pmemo(9,m)
!
    call cinema(pt,b,ps,ctl,stl,cfl,sfl,ttl,pmemo(9,m))
    pmemo(4,m)=ps(1)
    pmemo(5,m)=ps(2)
    pmemo(6,m)=ps(3)
    pmemo(8,m)=ttl
    call  kinemr(b,m,rl,taul)
    call  tmatur(np,v,m,p1,p2,mv,p,taul,tmat)
    pmemo(7,m)=tmat
    tprod(m)=tint
    mvcoll(m)=nucoll
    if(idpme(n1) == idpme(m).and.np == 2) then
! rescattering (elastic or charge exchange)  of the N1 particle
       iori(1,m)=iori(1,n1)
       iori(2,m)=iori(2,n1)
       iori(3,m)=iori(3,n1)+1
       do j=1,4
          hadr1(j,m)=hadr1(j,n1)
          hadr2(j,m)=hadr2(j,n1)
       enddo
       sori(m)=sori(n1)
    else
! production of a new particle M
       iori(1,m)=id1
       iori(2,m)=id2
       iori(3,m)=0
       do j=1,4
          hadr1(j,m)=hadi1(j)
          hadr2(j,m)=hadi2(j)
       enddo
       sori(m)=ssor
    endif
    myp(m)=1
    myt(m)=1
    myy(m)=1
    if(imemo(2,m).ne.0)  myp(m)=0
    if(imemo(2,m).ne.0)  myt(m)=0
    if(imemo(2,m).ne.0)  myy(m)=0
    if(imemo(5,m) == 0)  go  to  146
    taum0= dble(imemo(5,m))/1000.
    temp1 = pmemo(9,m)
!      call err_chk(1,'LaqPauli.f',"880",1,temp1)
    tauml=taum0*(1.+pmemo(8,m)/temp1)
    imemo(5,m)=intg((tauml+pmemo(7,m))*1000.)
    if(imemo(5,m) == 0)  imemo(5,m)=1
146 continue
    x=(p1(1)+p2(1))/2.
    y=(p1(2)+p2(2))/2.
    z=(p1(3)+p2(3))/2.
    vrk=x*vta(1)+y*vta(2)+z*vta(3)
    dr=rint*rndm(-1.0_real64)**(1./3.)
    fi=6.283185*rndm(-1.0_real64)
    ct=1.-2.*rndm(-1.0_real64)
    temp1 = 1.-ct**2
!      call err_chk(1,'LaqPauli.f',"893",2,temp1)
    st=sqrt(temp1)
    dx=dr*st*cos(fi)
    dy=dr*st*sin(fi)
    dz=dr*ct
    temp1 = gta+1.
!      call err_chk(1,'LaqPauli.f',"899, 900, 901",1,temp1)
    pmemo(1,m)=radt(1)+x-vta(1)*vrk*gta/temp1+dx
    pmemo(2,m)=radt(2)+y-vta(2)*vrk*gta/temp1+dy
    pmemo(3,m)=radt(3)+z-vta(3)*vrk*gta/temp1+dz
    if(ncas >= ncpri)  write(16,301) &
         & m,(pmemo(k,m),k=1,9),(imemo(k,m),k=1,5),clider(m),idpme(m)
301 format(1x,'newp',i5,9(1x,f 8.3),2x,4i2,i15,f6.3,i5)
    if(l-np1) 47,48,48
47  l=l+1
    go to 33
48  ip=1
    mv=mv+np1-1
!----> do 53 k=1,9
52  do k=1,9
53     pmemo(k,n1)=pmemo(k,mv+1)
    end do
!----> do 54 k=1,5
    do k=1,5
54     imemo(k,n1)=imemo(k,mv+1)
    end do
    myt(n1)=myt(mv+1)
    myy(n1)=myy(mv+1)
    tprod(n1)=tprod(mv+1)
    mvcoll(n1)=mvcoll(mv+1)
    iori(1,n1)=iori(1,mv+1)
    iori(2,n1)=iori(2,mv+1)
    iori(3,n1)=iori(3,mv+1)
    do j=1,4
       hadr1(j,n1)=hadr1(j,mv+1)
       hadr2(j,n1)=hadr2(j,mv+1)
    enddo
    sori(n1)=sori(mv+1)
    clider(n1)=clider(mv+1)
    idpme(n1)=idpme(mv+1)
    call eraij(n1)
    if(np1 == 0)              call repij(n1,mv+1)
    results%targExc%numTotal=results%targExc%numTotal+1.
    results%targExc%numHoles=results%targExc%numHoles+1.
!
    enext2=enex20+einc-esum+(an20-an2)*(0.940-eps2)
!
    return
  end subroutine pauli3

! =========================================
!
! =========================================


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine  paulid(mv,md,nre,ifd,np,ip)
!
!   Check Pauli principle for decay of
!   cascade particle  MD
!   boost NP secondary particles into observer's system, if
!   the collision is allowed
!   Number of cascade particles MV ==> MV + NP - 1
!

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
    common/ncasca/ncas,ncpri
    common/actim/tint/tprod/tprod(5999)
    common/porig/iori(3,5999)
    common/idn12/id1,id2
    common/cslid/clider(5999)
    common/nucoll/ nucoll,mvcoll(5999)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common /memorylaq/ pme(9,5999),ime(5,5999)
    common /hcasc/ an1,an2,zn1,zn2,t0,eps1,eps2,vpi,a1,a2,c1,c2,d1,d2, &
         & r0n1,r0n2,tf01,tf02,rm1,rm2
    common /nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre, &
         & vep(3),vet(3),gep,get
    common/idpme/ idpme(5999)
    common /sori/ sori(5999),ssor
    common /hadr1/hadr1(4,5999),hadr2(4,5999),hadi1(4),hadi2(4)
    dimension  ps(3),vn(3),rn(3),pl(3),pin(9),iin(5)
    mvn=0
    resm=pme(9,md)
!----> do  108  l=1,np
    do   l=1,np
!----> do   8  k=1,3
       do    k=1,3
          pme(k,mv+l)=pme(k,md)
8         continue
       end do
       if(ime(3,mv+l) == 0.and.ime(4,mv+l) == 1)  mvn=mv+l
108    continue
    end do
    if(ifd == 0)  go  to  113
    if(mvn == 0)  go  to  113
!----> do  9  k=1,9
    do   k=1,9
9      pin(k)=pme(k,mvn)
    end do
!----> do 10  k=1,5
    do  k=1,5
10     iin(k)=ime(k,mvn)
    end do
!----> do  13  nu=1,2
    do   nu=1,2
       if(nu == 1.and.an1 < 2.1)   go  to  13
       if(nu == 2.and.an2 < 2.1)   go  to  13
!----> do  12  k=1,3
       do   k=1,3
          pl(k)=pme(k+3,mvn)
          if(nu == 2)  go  to  11
          vn(k)=-vpr(k)
          rn(k)=radp(k)
          go  to  12
11        vn(k)=-vta(k)
          rn(k)=radt(k)
12        continue
       end do
       temp1 = 1.-vn(1)**2-vn(2)**2-vn(3)**2
!      call err_chk(1,'LaqPauli.f',"997",2,temp1)
       temp2 = sqrt(temp1)
!      call err_chk(1,'LaqPauli.f',"999",1,temp2)
       gn=1./temp2
       temp3 = gn+1.
!      call err_chk(1,'LaqPauli.f',"1002",1,temp3)
       gg=gn*gn/temp3
       vr=(pme(1,mvn)-rn(1))*vn(1)+ &
            & (pme(2,mvn)-rn(2))*vn(2)+ &
            & (pme(3,mvn)-rn(3))*vn(3)
       pin(1)=pme(1,mvn)-rn(1)+vn(1)*vr*gg
       pin(2)=pme(2,mvn)-rn(2)+vn(2)*vr*gg
       pin(3)=pme(3,mvn)-rn(3)+vn(3)*vr*gg
       temp1 =pin(1)**2+pin(2)**2+pin(3)**2
!      call err_chk(1,'LaqPauli.f',"1011",2,temp1)
       rmod=sqrt(temp1)
       if(nu == 1.and.rmod > rm1)  go  to  13
       if(nu == 2.and.rmod > rm2)  go  to  13
       call  cinema(pl,vn,ps,cts,sts,cfs,sfs,ts,pin(9))
       pin(4)=ps(1)
       pin(5)=ps(2)
       pin(6)=ps(3)
       pin(8)=ts
       if(nu == 1) tf=potenq(pin,iin,a1,c1,d1,tf01,vpi,eps1)-eps1
       if(nu == 2) tf=potenq(pin,iin,a2,c2,d2,tf02,vpi,eps2)-eps2
       if(ts > tf)   go  to  13
       ip=0
       tau0= taun(nre)
       temp2 = pme(9,md)
!      call err_chk(1,'LaqPauli.f',"1026",1,temp2)
       taul=tau0*(1.+pme(8,md)/temp2)
       pme(7,md)=0.
       ime(5,md)=intg((taul+pme(7,md))*1000.)
       if(ime(5,md) == 0)  ime(5,md)=1
       if(ncas >= ncpri) write(16,600) ts,tf,nu,ime(5,md)
600    format(1x,'paulid forbided: tn=',f6.3,' < tf=',f6.3,' nu=',i2,1x, &
            & 'ime(5,md)=',i10)
       return
13     continue
    end do
113 ip=1
    nucoll=nucoll+1
!----> do  114  l=1,np
    do   l=1,np
       pme(7,mv+l)=0.
       tprod(mv+l)=tint
! for decay of L particle the IORI(2,...)=0 !
       iori(1,mv+l)=id1
       iori(2,mv+l)=0
       iori(3,mv+l)=0
       do j=1,3
          hadr1(j,mv+l)=pme(j+3,md)
          hadr2(j,mv+l)=0.
       enddo
       hadr1(4,mv+l)=pme(9,md)
       hadr2(4,mv+l)=0.
       sori(mv+l)=resm
       clider(mv+l)=1.
       mvcoll(mv+l)=nucoll
       myp(mv+l)=ifd
       myt(mv+l)=ifd
       myy(mv+l)=ifd
       if(ime(2,mv+l).ne.0) myp(mv+l)=0
       if(ime(2,mv+l).ne.0) myt(mv+l)=0
       if(ime(2,mv+l).ne.0) myy(mv+l)=0
       mvl=mv+l
       if(ncas >= ncpri)  write(16,301) &
            & mvl,(pme(k,mv+l),k=1,9),(ime(k,mv+l),k=1,5),idpme(mv+l)
301    format(1x,'newp',i5,9(1x,f 8.3),2x,4i2,i15,i5)
114    continue
    end do
!----> do  14  k=1,9
    do   k=1,9
14     pme(k,md)=pme(k,mv+np)
    end do
!----> do  15  k=1,5
    do   k=1,5
15     ime(k,md)=ime(k,mv+np)
    end do
    myp(md)=ifd
    myt(md)=ifd
    myy(md)=ifd
    if(ime(2,md).ne.0)  myp(md)=0
    if(ime(2,md).ne.0)  myt(md)=0
    if(ime(2,md).ne.0)  myy(md)=0
    tprod(md)=tint
! for decay of MD particle the IORI(2,...)=0 !
    iori(1,md)=id1
    iori(2,md)=0
    iori(3,md)=0
    do j=1,4
       hadr1(j,md)=hadr1(j,mv+np)
       hadr2(j,md)=hadr2(j,mv+np)
    enddo
    sori(md)=resm
    mvcoll(md)=nucoll
    clider(md)=1.
    idpme(md)=idpme(mv+np)
    call eraij(md)
    mv=mv+np-1
    return
  end

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  subroutine pauli4(p1,p2,ip1,ip2,n1,n2,v,np,mv,tint,ip)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h,o-z), integer(int32) (i-n)
!
!   Check Pauli principle for collision of
!   cascade particle  (N1,P1,IP1) with cascade particle  (N2,P2,IP2)
!   (and N3 for absorption by NN pair), and
!   boost NP secondary particles into observer's system, if
!   the collision is allowed
!   Number of cascade particles MV ==> MV + NP - 2
!
    real(real64) ::  masn
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/ncasca/ncas,ncpri
    common/tprod/tprod(5999)
    common/porig/iori(3,5999)
    common/idn12/id1,id2
    common/cslid/clider(5999)
    common/nucoll/ nucoll,mvcoll(5999)
    common/memorylaq/pme(9,5999),ime(5,5999)
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    common/nucsp/vpr(3),vta(3),radp(3),radt(3),vev(3),vre(3),gev,gre &
         & ,vep(3),vet(3),gep,get
    common /rint/ rint
    common /idpme/ idpme(5999)
    common /sori/ sori(5999),ssor
    common /hadr1/hadr1(4,5999),hadr2(4,5999),hadi1(4),hadi2(4)
    dimension p1(9),p2(9),ip1(5),ip2(5),ps(3),v(3),pn1(3), &
         & pn2(3),p(3),pl(3),bpr(3),bta(3),rl(3)
    masn=0.940
!----> do  9  l=1,3
    do   l=1,3
       bpr(l)=-vpr(l)
9      bta(l)=-vta(l)
    end do
    x=(p1(1)+p2(1))/2.
    y=(p1(2)+p2(2))/2.
    z=(p1(3)+p2(3))/2.
    if((iabs(ip1(4))+iabs(ip2(4))) <= 0)  go  to  19
    if(np-2)19,10,13
!----> do 11 l=1,9
10  do l=1,9
11     pme(l,mv+2)=pme(l,mv+3)
    end do
!----> do 12 l=1,5
    do l=1,5
12     ime(l,mv+2)=ime(l,mv+3)
    end do
    clider(mv+2)=clider(mv+3)
    idpme(mv+2)=idpme(mv+3)
    go to 16
!----> do 14 l=1,9
13  do l=1,9
       temp=pme(l,mv+2)
       pme(l,mv+2)=pme(l,mv+3)
14     pme(l,mv+3)=temp
    end do
!----> do 15 l=1,5
    do l=1,5
       itemp=ime(l,mv+2)
       ime(l,mv+2)=ime(l,mv+3)
15     ime(l,mv+3)=itemp
    end do
    temc=clider(mv+2)
    clider(mv+2)=clider(mv+3)
    clider(mv+3)=temc
    idte=idpme(mv+2)
    idpme(mv+2)=idpme(mv+3)
    idpme(mv+3)=idte
!----> do 17 l=1,2
16  do l=1,2
       if(ime(3,mv+l).ne.0.or.ime(5,mv+l).ne.0.or.ime(4,mv+l).ne.1) &
            & go  to  17
       ps(1)=pme(4,mv+l)
       ps(2)=pme(5,mv+l)
       ps(3)=pme(6,mv+l)
       call kinemq(ps,v,pl,ctl,stl,cfl,sfl,tl,masn)
       if(anucl1 < 2.1)  go  to  117
       call kinemq(pl,bpr,pn1,ct1,st1,cf1,sf1,t1,masn)
       call rpts(x,y,z,xp,yp,zp,1)
       tfr1=tfermiq(xp,yp,zp,1)
       iforb=1
       if(t1 <= tfr1)  go  to  18
117    call kinemq(pl,bta,pn2,ct2,st2,cf2,sf2,t2,masn)
       call rpts(x,y,z,xt,yt,zt,2)
       tfr2=tfermiq(xt,yt,zt,2)
       iforb=2
       if(t2 <= tfr2)  go  to  18
17     continue
    end do
    go  to  19
18  ip=0
    if(ncas >= ncpri) then
       if(iforb == 1) write(16,302) t1,tfr1
       if(iforb == 2) write(16,303) t2,tfr2
302    format(1x,'pauli4 forbided: t1=',f6.3,' < tfr1=',f6.3)
303    format(1x,'pauli4 forbided: t2=',f6.3,' < tfr2=',f6.3)
    endif
    return
19  ip=1
    nucoll=nucoll+1
!----> do  21  l=1,np
    do   l=1,np
       m=mv+l
       p(1)=pme(4,m)
       p(2)=pme(5,m)
       p(3)=pme(6,m)
       call  kinemq(p,v,pl,ctl,stl,cfl,sfl,tl,pme(9,m))
       pme(4,m)=pl(1)
       pme(5,m)=pl(2)
       pme(6,m)=pl(3)
       pme(8,m)=tl
       call  kinemr(v,m,rl,taul)
       call  tmatur(np,v,m,p1,p2,mv,p,taul,tmat)
       pme(7,m)=tmat
       tprod(m)=tint
       iori(1,m)=id1
       iori(2,m)=id2
       iori(3,m)=0
       do j=1,4
          hadr1(j,m)=hadi1(j)
          hadr2(j,m)=hadi2(j)
       enddo
       sori(m)=ssor
       if(idpme(n1) == idpme(m).and.np == 2) then
! for rescattering of the N1 particle
          iori(1,m)=iori(1,n1)
          iori(2,m)=iori(2,n1)
          iori(3,m)=iori(3,n1)+1
          do j=1,4
             hadr1(j,m)=hadr1(j,n1)
             hadr2(j,m)=hadr2(j,n1)
          enddo
          sori(m)=sori(n1)
       endif
       if(idpme(n2) == idpme(m).and.np == 2) then
! for rescattering of the N2 particle
          iori(1,m)=iori(1,n2)
          iori(2,m)=iori(2,n2)
          iori(3,m)=iori(3,n2)+1
          do j=1,4
             hadr1(j,m)=hadr1(j,n2)
             hadr2(j,m)=hadr2(j,n2)
          enddo
          sori(m)=sori(n2)
       endif
       mvcoll(m)=nucoll
       myp(m)=1
       myt(m)=1
       myy(m)=1
       if(ime(2,m).ne.0)    myp(m)=0
       if(ime(2,m).ne.0)    myt(m)=0
       if(ime(2,m).ne.0)    myy(m)=0
       if(ime(5,m) == 0)    go  to  20
       taum0= dble(ime(5,m))/1000.
       temp1 = pme(9,m)
!      call err_chk(1,'LaqPauli.f','1236',1,temp1)
       tauml=taum0*(1.+pme(8,m)/temp1)
       ime(5,m)=intg((tauml+pme(7,m))*1000.)
       if(ime(5,m) == 0)  ime(5,m)=1
20     continue
       dr=rint*rndm(-1.0_real64)**(1./3.)
       fi=6.283185*rndm(-1.0_real64)
       ct=1.-2.*rndm(-1.0_real64)
       temp1 = 1.-ct**2
!      call err_chk(1,'LaqPauli.f','1245',2,temp1)
       st=sqrt(temp1)
       dx=dr*st*cos(fi)
       dy=dr*st*sin(fi)
       dz=dr*ct
       pme(1,m)=(p1(1)+p2(1))/2.+dx
       pme(2,m)=(p1(2)+p2(2))/2.+dy
       pme(3,m)=(p1(3)+p2(3))/2.+dz
       if(ncas >= ncpri)  write(16,301) &
            & m,(pme(k,m),k=1,9),(ime(k,m),k=1,5),clider(m),idpme(m)
301    format(1x,'newp',i5,9(1x,f 8.3),2x,4i2,i15,f6.3,i5)
21     continue
    end do
    m1=mv+np-1
    m2=mv+np
!----> do  22  k=1,9
    do   k=1,9
       pme(k,n1)=pme(k,m1)
22     pme(k,n2)=pme(k,m2)
    end do
    myp(n1)=myp(m1)
    myt(n1)=myt(m1)
    myy(n1)=myy(m1)
    tprod(n1)=tprod(m1)
    iori(1,n1)=iori(1,m1)
    iori(2,n1)=iori(2,m1)
    iori(3,n1)=iori(3,m1)
    sori(n1)=sori(m1)
    mvcoll(n1)=mvcoll(m1)
    clider(n1)=clider(m1)
    idpme(n1)=idpme(m1)
!----> do  23  k=1,5
    do   k=1,5
       ime(k,n1)=ime(k,m1)
23     ime(k,n2)=ime(k,m2)
    end do
    myp(n2)=myp(m2)
    myt(n2)=myt(m2)
    myy(n2)=myy(m2)
    tprod(n2)=tprod(m2)
    iori(1,n2)=iori(1,m2)
    iori(2,n2)=iori(2,m2)
    iori(3,n2)=iori(3,m2)
    sori(n2)=sori(m2)
    mvcoll(n2)=mvcoll(m2)
    clider(n2)=clider(m2)
    idpme(n2)=idpme(m2)
    call eraij(n1)
    call eraij(n2)
    if(np == 1) call repij(n1,m1)
    mv=mv+np-2
    return
  end
