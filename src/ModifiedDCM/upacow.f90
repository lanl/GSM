
  subroutine  upacow(mv,mvu,iu)

! ====================================================================
!
!
!
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(inout) :: mv
    integer(int32), intent(inout) :: mvu
    integer(int32), intent(in   ) :: iu

    integer(int32) :: i, iper, j, jj, m, mv1

    logical, save :: first = .true.
    logical :: validpart

! ====================================================================

    real(real64) :: an1,an2,zn1,zn2,enext1,enext2,pnucl1, &
         & pnucl2,amnuc1,amnuc2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3), &
         & pnucl2(3),amnuc1(3),amnuc2(3)
    real(real64) :: tint, tprod
    common/actim/tint/tprod/tprod(5999)
    integer(int32) :: iori
    common/porig/iori(3,5999)
    integer(int32) :: nucoll, mvcoll
    common/nucoll/ nucoll,mvcoll(5999)
    real(real64) :: clider
    common/cslid/clider(5999)
    integer(int32) :: kobr
    real(real64) :: blab, glab
    common/bglab/blab,glab,kobr
    real(real64) :: pm
    integer(int32) :: im
    common/memorylaq/pm(9,5999),im(5,5999)
    integer(int32) :: mpa, myp, myt, myy
    common/activ/mpa(300),myp(5999),myt(5999),myy(5999)
    integer(int32) :: npimi
    common/npimi/npimi
    integer(int32) :: ncas, ncpri
    common /ncasca/ ncas,ncpri
    integer(int32) :: intcc
    common/intcc/intcc
    integer(int32) :: idpme
    common /idpme/ idpme(5999)
    real(real64) :: tlimit
    common/tlimit/tlimit
    real(real64) :: sori, ssor
    common /sori/ sori(5999),ssor

! ====================================================================

    mv1=mv
    if(mv == 0)   return
    m=1
9   continue
    if(iu == 1)   go  to  10
    if(an1 < 1.) myp(m)=0
    if(an2 < 1.) myt(m)=0
    if(myp(m).ne.0.or.myt(m).ne.0)                 go to 20
    if(im(5,m).ne.0)                               go to 20
    if(intcc.ne.0.and.im(2,m) == 0)                go to 20
!  kkg  11/05/03
    if(im(2,m).ne.0.and.nucoll == 0)               go to 20
!  kkg  11/05/03
10  continue
    i=mvu+1
    if((11*i) < 66000)       go  to  201
    iper=1
    write(16,'(15x,''array  up1 is exceeded in upacow''/)')
    write( *,'(15x,''array  up1 is exceeded in upacow''/)')
    return
201 continue


! Collect fragment information (transfers INC information to GSM's variables)
    call collectfragment(m)

    if(im(4,m) == 0.and.im(1,m) == -1.and.im(3,m) == 0) &
         & npimi=npimi+1
    mvu=mvu+1
    mv1=mv1-1
    if(ncas >= ncpri) &
         & write(16,202) m,mvu,mv1,(pm(j,m),j=1,9),(im(jj,m),jj=1,5)
    if(m > mv1)   go  to  21
    do j=1,9
       pm(j,m)=pm(j,mv1+1)
    end do
    do j=1,5
       im(j,m)=im(j,mv1+1)
    end do
    myp(m)=myp(mv1+1)
    myt(m)=myt(mv1+1)
    myy(m)=myy(mv1+1)
    tprod(m)=tprod(mv1+1)
    iori(1,m)=iori(1,mv1+1)
    iori(2,m)=iori(2,mv1+1)
    iori(3,m)=iori(3,mv1+1)
    sori(m)=sori(mv1+1)
    mvcoll(m)=mvcoll(mv1+1)
    clider(m)=clider(mv1+1)
    idpme(m)=idpme(mv1+1)
    call repij(m,mv1+1)
    go  to  9
20  if(m == mv1)   go  to  21
    m=m+1
    go  to  9
21  mv=mv1

    return
! ====================================================================
202 format(5x,'error: upacov: m=',i5,2x,'mvu=',i5,2x,'mv1=',i5/ &
         & 1x,'==>',9(1x,f 8.3),2x,4i2,i15,'<==')
1200 format(5x, "error: ", a)
9999 format(3x, "comment: using ", a, " var. in upacow for tallies")
! ====================================================================
  end subroutine upacow
