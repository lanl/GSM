
  subroutine  cenum1(na,p,ip4,delta,nc,k1,k2,ik,nu)

! ==============================================================================
!
! Calculates all possible interaction partners in cylinder of
! radius DELTA along the direction of incoming particle (P,IP) in
! projectile (NU=1) or target (NU=2) nucleus;
! all partners are ordered in "time" and its numbers are stored in
! NUKC array
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none 
    integer(int32), intent(in   ) :: na
    real(real64),   intent(in   ) :: p(9)
    integer(int32), intent(in   ) :: ip4
    real(real64),   intent(in   ) :: delta
    integer(int32), intent(  out) :: nc
    integer(int32), intent(in   ) :: k1
    integer(int32), intent(  out) :: k2
    integer(int32), intent(in   ) :: ik
    integer(int32), intent(in   ) :: nu

    integer(int32) :: j, k, ka, kc
    real(real64)   :: del2, delta2, dev, dr1, pm, r, rj, z, z1
    real(real64), dimension(  3) :: rin1, rc1
    real(real64), dimension(100) :: zkc

! ==============================================================================

    real(real64) :: xc, yc, zc
    integer(int32) :: iz
    common/center/xc(2,300),yc(2,300),zc(2,300),iz(2,300)
    integer(int32) :: nukc
    common/cenpar/nukc(100)
    real(real64) :: anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2, &
         & vpi,a1,a2,c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    integer(int32) :: ncas, ncpri
    common /ncasca/ ncas,ncpri

! ==============================================================================

    pm=sqrt(p(4)**2+p(5)**2+p(6)**2)
    if(ik)  9,9,16
9   continue

    nc=0
    k=1
    k2=0
    delta2=delta**2
    r=rm2
    if(nu == 1)   r=rm1
    rin1(3)=(p(1)*p(4)+p(2)*p(5)+p(3)*p(6))/pm
    dev=sqrt(p(1)**2+p(2)**2+p(3)**2-rin1(3)**2)
    if(dev > (r+delta))   go  to  26
104 rc1(3)=(xc(nu,k)*p(4)+yc(nu,k)*p(5)+zc(nu,k)*p(6))/pm
    dr1=rc1(3)-rin1(3)
    if(dr1 < 0.)   go  to  14
    if(abs(dr1) < 0.0001)   go  to  14
    del2=(xc(nu,k)-p(1))**2+(yc(nu,k)-p(2))**2+(zc(nu,k)-p(3))**2- &
         & dr1**2
    if(del2-delta2)   11,11,14
11  z1=rc1(3)
    ka=k
    if(nc == 0)  go  to  113
    j=1
12  if(z1 > zkc(j))  go  to  13
    z=zkc(j)
    zkc(j)=z1
    z1=z
    kc=nukc(j)
    nukc(j)=ka
    ka=kc
13  if(j == nc)   go  to  113
    j=j+1
    go  to  12
113 if(nc == 100)   go  to  114
    nc=nc+1
    zkc(nc)=z1
    nukc(nc)=ka
    go  to  14
114 k2=1

14  if(k-na)  15,26,26
15  k=k+1
    go  to  104
16  if(na-1) 107,107,108
107 k2=0
    go to 26
108 if(ip4 == 1)  go  to  107
17  j=1
18  if(k1-1)  20,19,20
19  j=2
20  r=sqrt((xc(nu,k1)-xc(nu,j))**2+(yc(nu,k1)-yc(nu,j))**2+ &
         & (zc(nu,k1)-zc(nu,j))**2)
    k2=j
21  if(j-k1)  22,24,22
22  rj=sqrt((xc(nu,k1)-xc(nu,j))**2+(yc(nu,k1)-yc(nu,j))**2+ &
         & (zc(nu,k1)-zc(nu,j))**2)
    if(rj-r)  23,24,24
23  r=rj
    k2=j
24  if(j-na)  25,26,26
25  j=j+1

!   kkg 11/05/03
26  continue
    if(ncas >= ncpri)  then
       write(*,*) ' cenum1: nu,nc,k1,k2=',nu,nc,k1,k2
    endif
    return
! ==============================================================================
  end subroutine cenum1
