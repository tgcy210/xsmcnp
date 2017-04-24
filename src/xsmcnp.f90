!*************************************************************************************
!  xsmcnp : 
!  convert to/from PENTRAN format (row-wise P0) xs file to/from MCNP multigroup format
!  Higher order (for odd P_N) supported  for xs-> MCNP
!*************************************************************************************
module mParaset
    character(LEN=40) :: ver_num='xsmcnp v2.45 '
    character(LEN=20) :: logfile='xsmcnp.log'
    integer ::  LOGUNIT=101

    character(LEN=80) :: xsfile='' , grpfile='grp_erg.bnd', chifile='grp_mat.chi', nvfile='grp_mat.nu', curfile
    character(LEN=80) :: mgxs='mgxs' , form_str=''
    character(LEN=30) :: prbname=''
    integer :: xsunit=21, grpunit=22, outunit=23, xsout=30
    integer :: IsGetName=0, IsNubarFileRequired=0, IsChiFileRequired=0

    real*8, parameter :: MySmallestNumber=1e-10
end module

module mLineParse
    integer, parameter :: leg_line=10000
    character (LEN=3) :: cmtchar_end='/!#'
    character (LEN=54) :: cmtchar='/!abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    !character (LEN=16)   :: digchar='0123456789.eE +-'
    character (LEN=13)   :: digchar='0123456789.+-'
    character (LEN=2)    :: scichar='eE'
    integer :: IsCmtLine=0

    character(LEN=leg_line) :: line_buf=''
    !maxium length of keyword
    integer :: len_keyword=15
    !number of numbers in a line
    integer :: itm=0
    !number of words in aline
    integer :: iwd=0
end module


!xs data module
module mXslib


    integer :: num_mat=0, num_grp=0, num_line=0, num_col=0, num_cmt=0
    !for n-gamma coupled mcnp library
    !p_atom is index array for photo-atomic reaction 
    !.e.g. p_atom(1)=2 means H-1 gamma library is located in mat talbe 2 
    integer :: num_mat_p=0, num_mat_n=0, p_atom(150)=0
    integer :: num_grp_p=0, num_grp_n=0, num_grp_i
    integer :: num_mat_r
    integer, dimension(:), pointer :: n_map
    real, dimension(:), pointer :: upper_bnd_n, upper_bnd_p

    integer :: IsUpSca=0 , leg_order=0, leg_mgxs=0, xs_type=0
    !neutron (1) or photon (1) xs file
    integer :: np_flag=1  

    type tPorder
      real,dimension(:),pointer :: sig_a, mv_f, sig_t
      real,dimension(:,:),pointer :: sig_gp2g
      character(len=80) :: mtag=""
    end type tPorder

    type tMat
      integer num
      character(Len=20) name
      character(Len=20) mat   !for MCNP mat file names
      character(Len=10) zaid
      real :: mass=0
      integer :: leg_m=0, z_num=0, a_num=0,leg_p=0, leg_n=0
      integer :: sec_type=0, sec_ngrp=0, IsSkip=0
      integer :: nxs(16)=0, jxs(32)=0
      !integer :: IsFissionable
      real,dimension(:),pointer :: sig_s, chi, nubar, sig_f, upper_bnd
      type(tPorder),dimension(:),pointer :: xs_p
    end type tMat

    type(tMat), dimension(:), allocatable :: mat_info 

    integer, allocatable :: fis_flag(:)
    character(LEN =8) str_date
    character(LEN =12) str_time
    character(LEN=10) str_date_mg
    !MCNP data array
    !chi(num_grp, num_mat)
    real,  allocatable :: grp_erg_bnd(:), grp_mat_chi(:,:), grp_mat_nu(:,:)
    real :: nu_const=2.35, low_erg_bin=1e-9	
    character(len=7) :: zaid="000.22m"
    integer :: block_len=1, block_len_act, block_len_p0, p0_locator
    real, allocatable :: stream_b2(:)
    integer :: line_start, line_seg
    integer :: sca_loc, IsFisFormat=1

    integer :: xs_ihm, xs_iht=3, xs_ihs=4
    !to store scattering moments 
    real, allocatable :: fpl(:)
end module

!Moudle mPolynomials
module mPoly
    integer :: num_poly  !num_ploy(Qpoly order) =(leg_mgxs+1)/2 
    integer :: ord_poly  !reduced orders for Pn downgrade if not positive definite
    type tPoly
        double precision, allocatable :: a_coeff(:), root(:), wgt(:)
        integer :: rank=0
    end type

    type (tPoly), dimension(:), allocatable :: q_poly
    double precision, allocatable :: pnl_rev(:,:)       !inverse Pn coeffs
    double precision, allocatable :: M_mor(:), N_mor(:) !M_mor: adjusted scattering moments, N_mor: q_poly norm 
    double precision, allocatable :: L_mor(:), u_mor(:), sig2_mor(:)  !intermediate arrays for calculating q_poly coeffs
    double precision, allocatable :: Lmat(:,:)  !for Positive definite check
    
contains
    !Pn(x): legendre polynomial
    real*8 function legendre(x,n)
        real*8 x
        integer n

        integer i
        real leg1,leg2

        if(n .eq. 0) then 
           legendre=1
           return
        endif

        if(n .eq. 1) then 
           legendre=x
           return
        endif

        legendre=0
        leg2=1
        leg1=x

        do i=2,n
            legendre=(2.0*i-1)*x/i*leg1-(i-1)*leg2/i
            leg2=leg1
            leg1=legendre
        enddo

        return

    end function
    !P'n(x): derivative
    real*8 function dlegendre(x,n)
        real*8 x
        integer n

        integer i,m
        real leg1,leg2

        m=n-1
        if(m .eq. 0) then 
         dlegendre=1
         return
        endif

        if(m .eq. 1) then 
         dlegendre=2*x
         return
        endif

        dlegendre=0
        leg2=1
        leg1=2*x

        do i=2,m
            dlegendre=(2.0*i+1)*(i+1)*x*leg1/(i*(i+2))-(i+1)*leg2/(i+2)
            leg2=leg1
            leg1=dlegendre
        enddo

        dlegendre=0.5*(n+1)*dlegendre

        return

    end function
    !calculate pnl, called with n=leg_mgxs
    subroutine InitPoly (n)
        integer n

        real*8 areal, breal
        integer j, k

        allocate(pnl_rev(0:n, 0:n))

        pnl_rev=0
        pnl_rev(0,0)=2

        do k=1, n
          do j=0, k
          if(j+1 .gt. k-1) then
            areal=0
          else
            areal=pnl_rev(k-1,j+1)
          endif
          if(j .eq. 0 .or. j-1 .gt. k-1) then
           breal=0
          else
           breal=pnl_rev(k-1,j-1)
          endif
  
          pnl_rev(k,j)=((j+1)*areal+j*breal)/(2*j+1)

          enddo
        enddo
        !allocate
        num_poly=(n+1)/2
        allocate(M_mor(0:n), N_mor(0:n))
        allocate(L_mor(num_poly), u_mor(num_poly), sig2_mor(0:num_poly))
        M_mor=0.0d0; N_mor=0.0d0
        L_mor=0.0d0; u_mor=0.0d0;  sig2_mor=0.0d0
        
        k=num_poly-1
        allocate( Lmat(0:k, 0:k) ) ; Lmat=0.0d0
        allocate(q_poly(0:num_poly))
        do k=0, num_poly
           allocate(q_poly(k)%a_coeff(0:k))
           q_poly(k)%a_coeff=0
        enddo 
    end subroutine


end module

program xsmcnp
use mXslib
use mParaset
use mLineParse

logical :: ex=.FALSE.
integer :: cmt_err=0
integer :: i

write(*, *) ver_num
write(*,*) " Convert xs format for MCNP multigroup "

i = IARGC( )
if(i .gt. 0)  then 
 call GETARG(1, xsfile)
 inquire(file=xsfile, exist=ex)
 if(.not. ex) then 
  write(*, "('file doesnot exist: ',A)") trim(xsfile)
  IsGetName=1
 endif
else
 IsGetName=1
endif
if(IsGetName .eq. 1) then
  do while(.not. ex )
    write(*,*) "Enter PENTRAN or MCNP multigroup xs filename (Cirl+C exit):"
    read(*,*) xsfile
    inquire(file=xsfile, exist=ex)
    if ( .not. ex) then
      write(*, "('file doesnot exist: ',A)") trim(xsfile)
    endif
  enddo
endif

xsfile=adjustl(xsfile)
i=index(xsfile,'.')

if( i .gt. 0) then
prbname=xsfile(1:i-1)
else
prbname=xsfile
endif

open(unit=xsunit, file=trim(xsfile), err=1001)
open(unit=LOGUNIT, file=trim(logfile), err=1001)
write(form_str,"(A, '  Log generated on ' )") trim(ver_num)
call TimeStamp(LOGUNIT) 

read (xsunit, * ,  err=1001, end=1001)
IsCmtLine=0
!check five blank lines for mcnp xs file
do i=1, 5  
read (xsunit, "(A)",  err=1001,end=1001) line_buf
if (len_trim(line_buf) .eq. 0) then
 IsCmtLine=IsCmtLine-1
else
 exit
endif
enddo

rewind(xsunit) 

if(IsCmtLine .le. -1) then
write(LOGUNIT,"('Input cross section file: ',A, ' (mcnp multigroup format) ' )") trim(xsfile)
call mcnp2xs
write(*,*) 'conversion done: from mcnp multigroup to pentran format'
write(LOGUNIT,*) 'conversion done: from mcnp multigroup to pentran format'
else
write(LOGUNIT,"('Input cross section file: ',A, ' (pentran format) ' )") trim(xsfile)
call xs2mcnp
write(LOGUNIT,*) 'conversion done: from pentran to mcnp multigroup format'
endif

goto 2000
1001 write(*,"('read xs file error in file: ', A )") trim(xsfile)
2000 close(LOGUNIT)

end

! convert xs to mcnp
subroutine xs2mcnp
use mXslib
use mParaset
use mLineParse
use mPoly

logical :: ex=.FALSE., ex2=.FALSE.
integer :: cmt_err=0, ierr
integer :: i, j, k,m
integer :: dflag=0, cflag=0
integer :: jxs_16, jxs_17, idx_lpnd, loc_lpnd
integer, allocatable :: line_mat(:), lent_mat(:)


!counting lines
i=0
do while(cmt_err .eq. 0) 
i=i+1
read (xsunit, "(A)",  iostat=cmt_err) line_buf
if(cmt_err .ne. 0) exit

if(len_trim(line_buf) .ge. leg_line) then
write(LOGUNIT, "('number of columns >' , I0)" ) leg_line
write(*, "('number of columns >' , I0)" ) leg_line
stop 'XS conversion failed'
endif 

line_buf= adjustl(line_buf)
IsCmtLine=0

call ParseLine(0)
if(itm .le. 2)  IsCmtLine=1

!counting comment lines and data lines
if(IsCmtLine .eq. 1) then
  cflag=cflag+1
  if(dflag .ne. 0) then 
    
    if(num_grp .eq. 0) then
      num_grp=dflag
    elseif(dflag .ne. num_grp) then
      write(*,"('ERROR when reading line number:', I0)") i
      write(*, "('Expecting ', I0, ' data lines for each mat, but ', I0, ' lines read')") num_grp, dflag
	  write(LOGUNIT,"('ERROR when reading line number:', I0)") i
      write(LOGUNIT, "('Expecting ', I0, ' data lines for each mat, but ', I0, ' lines read')") num_grp, dflag
      stop 'XS conversion failed'
    endif

    num_mat=num_mat+1
    dflag=0
  endif
  cycle
  
else
  dflag=dflag+1
  if(num_cmt .eq. 0) then
    num_cmt=cflag
  elseif(num_cmt .ne. cflag .and. dflag .eq. 0) then
      write(*,"('ERROR when reading line number:', I0)") i
      write(*, "('Expecting ', I0, ' comment lines for each mat, but ', I0, ' lines read')") num_cmt, cflag
	  write(LOGUNIT,"('ERROR when reading line number:', I0)") i
      write(LOGUNIT, "('Expecting ', I0, ' comment lines for each mat, but ', I0, ' lines read')") num_cmt, cflag
      stop 'XS conversion failed'    
  endif
  cflag=0
  call ParseLine(0)
  if(num_col .eq. 0) then
    num_col=itm
  elseif( num_col .ne. itm) then
    write(*,"('ERROR when reading line number:', I0)") i
    write(*, "('Expecting ', I0, ' columns , but ', I0, ' columns read')") num_col, itm
	write(LOGUNIT,"('ERROR when reading line number:', I0)") i
    write(LOGUNIT, "('Expecting ', I0, ' columns , but ', I0, ' columns read')") num_col, itm
    stop 'XS conversion failed'
  endif
endif
enddo  !do while

num_line=i 
if(dflag .ne. 0) num_mat=num_mat+1

ierr=1
  do while(ierr .ne. 0) 
    write(*,*) "Enter xs particle type (2 for photon, 1 for neutron) : "
    read(*,*, iostat=ierr) np_flag
    if(ierr .ne. 0) write (*,*) 'Wrong particle type, enter again (Ctrl+C to exit)'
    
	if( np_flag .eq. 2 ) then
      zaid='000.22g'
     else if(np_flag .eq. 1)  then
      zaid='000.22m'
     else
      write (*,*) 'only 2 or 1 allowed, enter again (Ctrl+C to exit)'
      ierr=1
     endif
  enddo
  
  write(LOGUNIT,"('xs particle type entered (2 for photon, 1 for neutron): ', I0 )") np_flag
  ierr=1
  do while(ierr .ne. 0) 
    write(*,*) "Enter Pn order of xs : "
    read(*,*, iostat=ierr) leg_order
    if(ierr .ne. 0) write (*,*) 'Wrong Pn order, enter again (Ctrl+C to exit)'
    if( leg_order .lt. 0)  then
     write (*,*) 'Pn ordder<0, enter again (Ctrl+C to exit)'
      ierr=1
     endif
  enddo
  write(LOGUNIT,"('Pn order of xs entered: ', I0 )") leg_order
  ierr=1
  if(leg_order .eq. 0) then
    leg_mgxs=0
  else
  do while(ierr .ne. 0) 
    write(*,*) "Enter Pn order of MCNP xs output: "
    read(*,*, iostat=ierr) leg_mgxs
    if(ierr .ne. 0) write (*,*) 'Wrong mcnp xs Pn order, enter again (Ctrl+C to exit)'
    if( leg_mgxs .gt. leg_order) then
      write (*,*) 'mcnp xs Pn ordder > PENTRAN xs Pn order, enter again (Ctrl+C to exit)'
      ierr=1
    endif
    if( leg_mgxs .ne. 0 .and. mod(leg_mgxs,2) .eq. 0) then
      write (*,*) 'mcnp xs Pn ordder has to be an odd number, enter again (Ctrl+C to exit)'
      ierr=1
     endif
  enddo
  endif
  write(LOGUNIT,"('Pn order of MCNP xs output : ', I0 )") leg_mgxs
  if(leg_mgxs .gt. 0 ) then
   ierr=1
   do while(ierr .ne. 0) 
    write(*,*) "Enter PENTRAN xs type (Enter 0 if (2l+1) is not factored in, otherwise Enter 1): "
    read(*,*, iostat=ierr) xs_type
    if(ierr .ne. 0) write (*,*) 'Wrong PENTRAN xs type, enter again (Ctrl+C to exit)'
    if( xs_type .ne. 0 .and. xs_type .ne. 1) then
     write (*,*) 'wrong xs type, only 0 or 1 allowed, enter again (Ctrl+C to exit)'
     ierr=1
    endif
  enddo
  endif
  write(LOGUNIT,"('PENTRAN xs type entered (0: (2l+1) not factored in, 1: factored in): : ', I0 )") xs_type

if(mod(num_mat, leg_order+1) .ne. 0) then
write(LOGUNIT,*) 'number of sections in xs file not equal to num_mat*(PnOrder+1)'
stop 'number of sections in xs file not equal to num_mat*(PnOrder+1)'
endif
num_mat=num_mat/(leg_order+1)

rewind(xsunit)

write(*, "('Number of Columns in xs file: ', I0)") num_col
write(LOGUNIT,*) '---------------checking xs data----------------------'
write(LOGUNIT, "('Number of Columns found in xs file: ', I0)") num_col
if(num_col .le. 3) then
write(LOGUNIT,*) 'number of columns in xs file < 3'
 stop 'number of columns in xs file < 3'
endif
if(num_col .gt. num_grp*2-1+xs_iht) then
  write(LOGUNIT,"('number of columns > ', I0, 'columns beyond ',I0, ' ignored' )") num_grp*2+2, num_grp*2+2
  write(*,"('number of columns > ', I0, 'columns beyond ',I0, ' ignored' )") num_grp*2+2, num_grp*2+2
endif

if(num_grp+xs_iht .eq. num_col) then
  IsUpSca=0
  xs_ihs=xs_iht+1
  xs_ihm=num_grp+xs_iht
elseif(num_grp*2-1+xs_iht .eq. num_col) then
  IsUpSca=1
  xs_ihs=xs_iht+num_grp
  xs_ihm=num_grp*2-1+xs_iht
elseif(num_col .gt. num_grp+xs_iht) then
  IsUpSca=1
  ierr=1
  do while(ierr .ne. 0) 
    write(*,*) "Enter self scattering position (ihs): "
    read(*,*, iostat=ierr) xs_ihs
    if(ierr .ne. 0) write (*,*) 'Wrong ihs number, enter again (Ctrl+C to exit)'
    if( xs_ihs .gt. num_col) write (*,*) 'ihs>num_col, enter again (Ctrl+C to exit)'
    if( xs_ihs .lt. 4) write (*,*) 'ihs>4, enter again (Ctrl+C to exit)'
  enddo
  xs_ihm=num_col
else
  write(LOGUNIT,*) 'wrong format of PENTRAN xs file'
  stop 'wrong format of PENTRAN xs file'
endif

write(LOGUNIT,"('self-scattering column located at: ', I0)") xs_ihs
if(num_mat .eq. 0 .or. num_grp .eq. 0 .or. num_cmt .eq. 0 .or. num_col .eq. 0) then
write(LOGUNIT,*) 'No data found in PENTRAN xs file'
stop 'No data found in PENTRAN xs file'
endif
!allocate xs data array
allocate (mat_info(num_mat))
allocate (fis_flag(num_mat))
fis_flag=0

do i=1, num_mat
! 
allocate(mat_info(i)%sig_s(num_grp) )
allocate(mat_info(i)%xs_p(leg_order+1))
do j=1, leg_order+1
allocate(mat_info(i)%xs_p(j)%sig_a(num_grp),&
         mat_info(i)%xs_p(j)%mv_f(num_grp),&
	 mat_info(i)%xs_p(j)%sig_t(num_grp) )
	 
allocate(mat_info(i)%xs_p(j)%sig_gp2g(num_grp,num_grp))
mat_info(i)%xs_p(j)%sig_a=0
mat_info(i)%xs_p(j)%mv_f=0
mat_info(i)%xs_p(j)%sig_t=0
mat_info(i)%xs_p(j)%sig_gp2g=0
enddo
enddo

!load xs
call LoadXs

do i=1, num_mat
do j=1, num_grp
  if(mat_info(i)%xs_p(1)%mv_f(j) .gt. 0) then
    fis_flag(i)=1
    IsNubarFileRequired=1
    IsChiFileRequired=1
    exit
  endif
enddo
enddo


write(*,"('XS file read sucessful: ', A )") trim(xsfile)
write(*, "('   Num of groups    : ', I0)")  num_grp
write(*, "('   Num of materials : ', I0)")  num_mat
if ( IsUpSca .eq. 0) then
write(*, "('   Upsacttering     : No')")
else
write(*, "('   Upsacttering     : yes')")
endif

write(LOGUNIT,*)
write(LOGUNIT,"('XS file read sucessful: ', A )") trim(xsfile)
write(LOGUNIT, "('   Num of groups    : ', I0)")  num_grp
write(LOGUNIT, "('   Num of materials : ', I0)")  num_mat
if ( IsUpSca .eq. 0) then
write(LOGUNIT, "('   Upsacttering     : No')")
else
write(LOGUNIT, "('   Upsacttering     : yes')")
endif

!write a full upscattering matrix xs file
if(IsUpSca .eq. 1 .and. num_grp*2-1+xs_iht .ne. num_col) then
call OutputXS
endif
!allocate MCNP array
allocate(grp_erg_bnd(num_grp), grp_mat_chi(num_grp, num_mat), grp_mat_nu(num_grp, num_mat))


!check group boundary file (upper boundary for each group)
inquire(file=trim(prbname)//'.grp' , exist=ex)
inquire(file='grp_erg.bnd', exist=ex2) 

if( ex) then !prbname.grp exists
curfile=trim(prbname)//'.grp'
open(unit=grpunit, file=trim(curfile), err=1002)
do i=1, 4
read(grpunit,*,iostat=ierr)
if(ierr .ne. 0) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
  write(LOGUNIT,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif 
enddo
do i=1, num_grp
read(grpunit,*,iostat=ierr) j, grp_erg_bnd(i)
if(ierr .ne. 0) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
  write(LOGUNIT,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
endif
enddo

 write(*,"('group energy boundary file read sucessful: ', A )") trim(curfile)
 write(LOGUNIT,"('group energy boundary file read sucessful: ', A )") trim(curfile)
 
 if(grp_erg_bnd(num_grp) .lt. low_erg_bin) low_erg_bin=grp_erg_bnd(num_grp)/2.0
 do j=2, num_grp
   if(grp_erg_bnd(j) .ge. grp_erg_bnd(j-1)) then
   
     write(*,"(A,':', 'grp upper boundary not in increasing order, from grp=', I0, ' to ', I0)") &
	   trim(curfile), j-1, j
	 write(LOGUNIT,"(A,':', 'grp upper boundary not in increasing order, from grp=', I0, ' to ', I0)") &
	   trim(curfile), j-1, j
	 write(LOGUNIT,"('Conversion Failed ' )") 
     stop 'Conversion Failed' 
   endif
 enddo
if(grp_erg_bnd(num_grp) .lt. low_erg_bin) low_erg_bin=grp_erg_bnd(num_grp)/2.0

elseif(ex2) then !grp_erg.bnd exists

 curfile='grp_erg.bnd'
 cmt_err=0
 i=1
 open(unit=grpunit, file=trim(curfile), err=1002) 
 do while(cmt_err .eq. 0) 
 
 read (grpunit, "(A)",  iostat=cmt_err) line_buf
 if(cmt_err .ne. 0) exit

 call ParseLine(1)
 if(itm .le. 0 .or. iwd .ge. 1) cycle
 k=min(num_grp, i+itm-1)
 read(line_buf, *, iostat=ierr) (grp_erg_bnd(j), j=i, k)
 if(ierr .ne. 0) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
   write(LOGUNIT,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif 
 if(k .eq. num_grp) exit
 i=k+1
 end do 
 close(grpunit)
 if(k .ne. num_grp) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
   write(LOGUNIT,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif
 
 write(*,"('group energy boundary file read sucessful: ', A )") trim(curfile)
 write(LOGUNIT,"('group energy boundary file read sucessful: ', A )") trim(curfile)
 if(grp_erg_bnd(num_grp) .lt. low_erg_bin) low_erg_bin=grp_erg_bnd(num_grp)/2.0
 do j=2, num_grp
   if(grp_erg_bnd(j) .ge. grp_erg_bnd(j-1)) then
      write(*,"('test01', 10(ES12.5,2x))") grp_erg_bnd
     write(*,"(A,':', 'grp upper boundary not in increasing order, from grp=', I0, ' to ', I0)") &
	   trim(curfile), j-1, j
	 write(LOGUNIT,"(A,':', 'grp upper boundary not in increasing order, from grp=', I0, ' to ', I0)") &
	   trim(curfile), j-1, j
	 write(LOGUNIT,"('Conversion Failed ' )") 
     stop 'Conversion Failed' 
   endif
 enddo
if(grp_erg_bnd(num_grp) .lt. low_erg_bin) low_erg_bin=grp_erg_bnd(num_grp)/2.0

!neither grp_erg.bnd nor prbname.grp exists
else
  do i=1, num_grp
    grp_erg_bnd(i)=(num_grp-i+1)*1.0
  enddo
endif

!check group chi (chi for each group and mat)
inquire(file=trim(prbname)//'.chi', exist=ex)
if(ex) then
curfile=trim(prbname)//'.chi'
else
inquire(file=trim(chifile), exist=ex)
curfile=trim(chifile)
endif 


if( ex) then
 cmt_err=0

 open(unit=grpunit, file=trim(curfile), err=1002) 
 matloop: do m=1, num_mat
 i=1
 lineloop: do while(cmt_err .eq. 0) 
 
 read (grpunit, "(A)",  iostat=cmt_err) line_buf
 if(cmt_err .ne. 0) exit matloop
 call ParseLine(1)
 if(itm .le. 0 .or. iwd .ge. 1) cycle lineloop
 k=min(num_grp, i+itm-1)

 read(line_buf,*, iostat=ierr,end=1003)  (grp_mat_chi(j,m), j=i,k)

 if(ierr .ne. 0) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp*num_mat
  write(LOGUNIT,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp*num_mat
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif
 if(k .eq. num_grp) cycle matloop
 i=k+1
 enddo lineloop !do while
 enddo matloop  !m
 close(grpunit)

 if(k .ne. num_grp .or. m .ne. num_mat+1) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp*num_mat
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif
 write(*,"('mat chi file read sucessful: ', A )") trim(curfile)
 write(LOGUNIT,"('mat chi file read sucessful: ', A )") trim(curfile)
elseif(IsChiFileRequired .eq. 0) then 
  do j=1 , num_mat
    grp_mat_chi(:,j)=0.0
    if(fis_flag(j) .ne. 0) grp_mat_chi(1,j)=1.0
  enddo

else
   write(*,"('Chi file missing: ', A )") trim(prbname)//'.chi or '//trim(chifile)
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed'    
endif

!check group nu (nu for each group and mat)
inquire(file=trim(prbname)//'.nub', exist=ex)
if(ex) then
curfile=trim(prbname)//'.nub'
else
inquire(file=nvfile, exist=ex)
curfile=trim(nvfile)
endif

if( ex) then

 cmt_err=0
 open(unit=grpunit, file=trim(curfile), err=1002) 
 nu_matloop: do m=1, num_mat
 i=1
 nu_lineloop: do while(cmt_err .eq. 0)
 read (grpunit, "(A)",  iostat=cmt_err) line_buf
 if(cmt_err .ne. 0) exit nu_matloop
 call ParseLine(1)
 if(itm .le. 0 .or. iwd .ge. 1) cycle nu_lineloop
 k=min(num_grp, i+itm-1)

 read(line_buf, *, iostat=ierr ) (grp_mat_nu(j,m), j=i,k)
 if(ierr .ne. 0) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp*num_mat
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif 
 if(k .eq. num_grp) cycle nu_matloop
 i=k+1
 enddo nu_lineloop
 enddo nu_matloop
 close(grpunit)
 if( k .ne. num_grp .or. m .ne. num_mat+1) then
  write(*,"(A,':', 'read error, ', I0, ' numbers required')") trim(curfile), num_grp*num_mat
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
 endif 

 write(*,"('mat nu file read sucessful: ', A )") trim(curfile)
 write(LOGUNIT,"('mat nu file read sucessful: ', A )") trim(curfile)
elseif(IsNubarFileRequired .eq. 0) then 
  grp_mat_nu=nu_const
else
  write(*,"('nubar file missing: ', A )") trim(prbname)//'.nub or '//trim(nvfile)
  write(LOGUNIT,"('Conversion Failed ' )") 
  stop 'Conversion Failed' 
endif

!output MCNP xs lib file
call date_and_time (str_date, str_time)
write(str_date_mg,"(A2,'/',A2,'/', A4)")  str_date(5:6),str_date(7:8), str_date(1:4)

!accumulate the length of the second block
!7 sections each with 
block_len_p0=num_grp*(num_grp+7)+1
if(leg_mgxs .gt. 0) then
block_len=block_len_p0+num_grp*num_grp*(leg_mgxs+1)+2
else
block_len=block_len_p0
endif
allocate(stream_b2(block_len))
allocate(line_mat(num_mat),lent_mat(num_mat) )

line_mat=0
lent_mat=0

allocate(fpl(0:leg_mgxs))

call InitPoly(leg_mgxs)

open(unit=outunit, file=trim(mgxs))
do i=1, num_mat
  !mat title
  write(outunit, "(I3,A7, 1x, f11.6, 1x, ES12.5, 1x, A12)") i, zaid,i*2.0, 0.0, str_date_mg 
  do j=1, 5
    write(outunit,*) 
  enddo
  if(IsFisFormat .eq. 0 .or. fis_flag(i) .eq. 1) then
    block_len_act=block_len
    p0_locator=block_len_p0
  else
    block_len_act=block_len-num_grp*3
    p0_locator=block_len_p0-num_grp*3
  endif
  !first block
  write(outunit, "(8(I8,1x))") block_len_act, i*1000, leg_mgxs, 0, num_grp, num_grp, num_grp, 0
  if(leg_mgxs .gt. 1) then  
  write(outunit, "(8(I8,1x))")  1, 1, 0, np_flag, 0, 0, 0, 0  !MORSE discrete angle scattering
  else
  write(outunit, "(8(I8,1x))")  0, 1, 0, np_flag, 0, 0, 0, 0  !for P1, cosine bin scattering, and isotropic scattering
  endif

  if(leg_mgxs .gt. 0) then
    jxs_16=block_len_act-1
  else
    jxs_16=0
  endif

  if(fis_flag(i) .ne. 0) then
    write(outunit, "(8(I8,1x))") 1, (num_grp*(k+1)+1, k=1, 5), 0, 0
    write(outunit, "(8(I8,1x))")  0, 0, 0, 0,  num_grp*7+1, 0, 0, jxs_16
  elseif (IsFisFormat .eq. 0 ) then
    write(outunit, "(8(I8,1x))")  1, num_grp*2+1, 0, 0, 0, num_grp*6+1, 0, 0
    write(outunit, "(8(I8,1x))")  0, 0, 0, 0,  num_grp*7+1, 0, 0, jxs_16
  else
    write(outunit, "(8(I8,1x))")  1, num_grp*2+1, 0, 0, 0, num_grp*3+1, 0, 0
    write(outunit, "(8(I8,1x))")  0, 0, 0, 0,  num_grp*4+1, 0, 0, jxs_16
  endif
 
  if(leg_mgxs .gt. 0) then
    jxs_17=block_len_act
  else
    jxs_17=0
  endif

  write(outunit, "(8(I8,1x))")  jxs_17, 0, 0, 0,  0, 0, 0, 0
  write(outunit, "(8(I8,1x))")  0, 0, 0, 0,  0, 0, 0, 0
  !second block
  m=0
  loc_lpnd=0
  !group center
  do j=1, num_grp-1
    m=m+1
    stream_b2(m)=(grp_erg_bnd(j)+grp_erg_bnd(j+1))/2
  enddo
  m=m+1
  stream_b2(m)=(grp_erg_bnd(num_grp)+low_erg_bin)/2

  !group width
  do j=1, num_grp-1
    m=m+1
    stream_b2(m)=grp_erg_bnd(j)-grp_erg_bnd(j+1)
  enddo
  m=m+1
  stream_b2(m)=grp_erg_bnd(num_grp)-low_erg_bin
  !total xs
  do j=1, num_grp
    m=m+1
    stream_b2(m)=mat_info(i)%xs_p(1)%sig_t(j)
  enddo

  !sigma fission
  if(IsFisFormat .eq. 0 .or. fis_flag(i) .eq. 1) then
  do j=1, num_grp
   m=m+1
   if(mat_info(i)%xs_p(1)%mv_f(j) .eq. 0 .or. grp_mat_nu(j,i) .eq. 0) then
     stream_b2(m)=0
   else
     stream_b2(m)=mat_info(i)%xs_p(1)%mv_f(j)/grp_mat_nu(j,i)
   endif
  enddo
  !nu
  do j=1, num_grp
    m=m+1
    stream_b2(m)=grp_mat_nu(j,i)
  enddo
  !chi
  do j=1, num_grp
    m=m+1
    stream_b2(m)=grp_mat_chi(j,i)
  enddo
  endif
  !sigma_a
  do j=1, num_grp
    m=m+1
    if(fis_flag(i) .eq. 1) then
    !mcnp sig_a excludes fission, but seems not used, so commented out
    stream_b2(m)=mat_info(i)%xs_p(1)%sig_a(j)  !& 
    !         - mat_info(i)%xs_p(1)%mv_f(j)/grp_mat_nu(j,i)
    else
    stream_b2(m)=mat_info(i)%xs_p(1)%sig_a(j)
    endif
  enddo
  m=m+1
  sca_loc=m
  stream_b2(m)=m+1
  
  !scattering
  do j=1, num_grp
    do k=1, num_grp
     m=m+1
     stream_b2(m)=mat_info(i)%xs_p(1)%sig_gp2g(j,k)

    enddo
  enddo
  if(leg_mgxs .gt. 0) then
    !LPND block higher Pn locator (see manual)
    idx_lpnd=0
    loc_lpnd=m
    do j=1, num_grp
     do k=1, num_grp
       m=m+1
       idx_lpnd=idx_lpnd+1
       stream_b2(m)=1+(idx_lpnd-1)*leg_mgxs

     enddo
    enddo

    !higher pn order scatering matrix
    do j=1, num_grp
      do k=1, num_grp
      m=m+1 
      fpl(0)=mat_info(i)%xs_p(1)%sig_gp2g(j,k)
      do idx_lpnd=1, leg_mgxs
         !get rid of 2l+1
         if(xs_type .eq. 1) &
             mat_info(i)%xs_p(idx_lpnd+1)%sig_gp2g(j,k)=&
            mat_info(i)%xs_p(idx_lpnd+1)%sig_gp2g(j,k)/(2*idx_lpnd+1)
         fpl(idx_lpnd)=mat_info(i)%xs_p(idx_lpnd+1)%sig_gp2g(j,k)
      enddo

      select case(leg_mgxs)
      case(1)
      if(mat_info(i)%xs_p(1)%sig_gp2g(j,k) .ne. 0) then
        stream_b2(m)=mat_info(i)%xs_p(2)%sig_gp2g(j,k)/mat_info(i)%xs_p(1)%sig_gp2g(j,k)
      else
        stream_b2(m)=0
      endif
      case(2:)
      if(mat_info(i)%xs_p(1)%sig_gp2g(j,k) .ne. 0) then
        call ConvertMnt(i,j,k,dflag)
        !if scattering moments does not satisfy the positive-definite condition, downgrade scattering to isotropic by setting LPND=0  
        if(dflag .gt. 0) then 
            stream_b2(loc_lpnd+j*k)=0
        endif
        do idx_lpnd=1, leg_mgxs-1
          stream_b2(m)=fpl(idx_lpnd)
		  m=m+1
		enddo
        stream_b2(m)=fpl(leg_mgxs)
	  else
        do idx_lpnd=1, leg_mgxs-1
          stream_b2(m)=0
		  m=m+1
		enddo
        stream_b2(m)=0
      endif
      case default
	   write(LOGUNIT,"('legendre order is not supported, conversion failed' )") 
       stop 'legendre order is not supported, conversion failed'
      end select
      

      enddo
    enddo
    !last two are locators indexed by jxs_16 and jxs_17
    m=m+1
    stream_b2(m)=p0_locator + 1
    m=m+1
    stream_b2(m)=stream_b2(m-1)+num_grp*num_grp
  endif  !leg_mgxs

  
  if(block_len_act .ne. m) stop 'm counting error'
!  block_len_act=m
  line_mat(i)=12+int((block_len_act-1)/4)+1
  lent_mat(i)=block_len_act

!output second block
k=int((sca_loc-1)/4)*4
j=mod(sca_loc-1,4)+1
write(outunit, "(4(ES19.12,x))")  stream_b2(1:k) 
select case(j)
  case(1)
   write(outunit, "(I19, x,3(ES19.12,x))")  sca_loc+1, stream_b2(k+2:k+4)
  case(2)
   write(outunit, "(ES19.12,x, I19, x,2(ES19.12, x))")  stream_b2(k+1),sca_loc+1,stream_b2(k+3:k+4)
  case(3)
   write(outunit, "(2(ES19.12,x), I19, x,ES19.12,x)")  stream_b2(k+1:k+2),sca_loc+1,stream_b2(k+4)
  case(4)
   write(outunit, "(3(ES19.12,x), I19, x)")   stream_b2(k+1:k+3),sca_loc+1
end select
write(outunit, "(4(ES19.12,x))")  stream_b2(k+5:block_len_act)
enddo  !num_mat
close(outunit)
write(*,*)
write(*, "('output file (mncp multigroup xs lib):  ', A)") trim(mgxs)
write(LOGUNIT,*) '----------------xs output files--------------------------'
write(LOGUNIT, "('output file (mncp multigroup xs lib):  ', A)") trim(mgxs)
!output xsdir file
open(unit=outunit, file='xsdir')
write(outunit, "('ATOMIC WEIGHT RATIOS')")

do i=1, num_mat
  write(outunit, "(I0,A3, f11.6 )") i, zaid, i*2.0
enddo
!Ref. MCNP5 manual Appendix F
write(outunit,"('DIRECTORY')")
line_seg=0
line_start=1
do i=1, num_mat
!  line_start=line_start+ line_seg
!  block_len_act=num_grp*(num_grp+4+3*(1-IsFisFormat*(1-fis_flag(i)) ) )+1
!  line_seg=12+int((block_len_act-1)/4)+1
  write(outunit,  "(I0,A7,1x,f11.6,1x, A, 1x,  I0, 1x, I0, 1x,  I0, 1x,    I0, 1x,  2(I0,1x), f6.2)") &
                    i, zaid, i*2.0, trim(mgxs),   0,     1,   line_start, lent_mat(i) , 0, 0, 0.00
 line_start=line_start+line_mat(i)
enddo
close(outunit)
write(*, "('output file (xsdir file):  ', A)") 'xsdir'
write(LOGUNIT, "('output file (xsdir file):  ', A)") 'xsdir'

!mcnp sample input deck
open(unit=outunit, file='mcnp.inp')
write (outunit, "('title: sample multigroup input deck with mcard')" ) 
write (outunit, "('c generated by xsmcnp, a code to convert PENTRAN xs (row format P0 only) to MCNP mg format' )" )

write (outunit, "('c group structure, fission chi and nu can be defined by optional input files of xsmcnp')")
write (outunit, "('c group upper bound defined in file : grp_erg.bnd (note: group 1 first)' )")
write (outunit, "('c fission chi defined in file : grp_mat.chi or valid gmix format; note: chi(num_grp, num_mat)' )")
write (outunit, "('c fission  nu defined in file : grp_mat.nu  or valid gmix format; note: nu(num_grp, num_mat)' )")
write (outunit, "('c group structure default is 0.001MeV 1.0 2.0 3.0 ... or valid gmix format ' )")
write (outunit, "('c chi by default is 1.0 for the first group, all other groups are zeros')")
! write (outunit, "('c nu by default is 2.40 if nu_sigmaf (second column) in pentran xs not equal zero')")
write (outunit, "('c ')")
write (outunit, "('c  cell cards')")
write (outunit, "('c  set cell density to 1.0, to make micro = macro xs')")
write (outunit, "('1    1   1.0  -1   imp:n=1')")
write (outunit, "('2      0       1   imp:n=0')")
write (outunit, *)
write (outunit, "('c  surface cards')")
write (outunit, "('1  rpp  -1.0  1.0  -1.0  1.0  -1.0  1.0')")
write (outunit, *)

if(np_flag .eq. 2) then
write (outunit, "('mode p ')") 
else
write (outunit, "('mode n ')") 
endif

write (outunit, "('c  source defination')")
write (outunit, "('c  point source at (0 0 0) in the first group')")
write (outunit, "('sdef erg=', ES12.5)") stream_b2(1)
write (outunit, "('c  material cards')")
do i=1, num_mat
  write (outunit, "('m', I0, 2x, I0, A7, '  1.0')") i, i, zaid
enddo
write (outunit, "('c  ')")
write (outunit, "('c  tally cards')")
if(np_flag .eq. 2) then
write (outunit, "('f4:p 1')")
else
write (outunit, "('f4:n 1')")
endif
i=num_grp
do while (i .gt. 0) 
  if (i .eq. num_grp) then
     write (outunit, "('e0 ', 3x,6(ES10.3,1x)  )") (grp_erg_bnd(k), k=i, max(i-5,1), -1)
  else 
     write (outunit, "( 6x,6(ES10.3,1x)  )") (grp_erg_bnd(k), k=i, max(i-5,1), -1)
  endif
  i=i-6
enddo
write (outunit, "('c  NOTE: MCNP energy bin from lower energy to higher, so the first bin is the last group')")
write(outunit, "('nps  1e3 ')") 
write(outunit, "('mgopt   f  ', I0)") num_grp
close (outunit)
write(*, "('output file (sample mcnp inp with material card):  ', A)") 'mcnp.inp'
write(LOGUNIT, "('output file (sample mcnp inp with material card):  ', A)") 'mcnp.inp'
return
1001 write(LOGUNIT, "('xsmcnp: read xs file error')")
stop "xsmcnp: read xs file error"
1002 write(LOGUNIT,"('open file error:', A)") trim(curfile)
write(*,"('open file error:', A)") trim(curfile)
stop 'Conversion failed'
1003 write(LOGUNIT,"('chi file error:', A)") trim(curfile)
write(*,"('chi file error:', A)") trim(curfile)
stop 'Conversion failed'

end subroutine


!read xs 
subroutine LoadXs
use mXslib
use mParaset

integer i, j, k
integer cmt,sgrp
real temp
integer xs_dist, xs_first,xs_last, xs_empty

xs_dist=xs_ihs - xs_iht-1
do k=1, num_mat
 do j=1, leg_order+1
   do cmt=1,num_cmt
     if(cmt .eq. 1) then   
       read(xsunit,"(A)",err=1002, end=1002)  mat_info(k)%xs_p(j)%mtag
     else
        read(xsunit,*,err=1002, end=1002)
     endif
   enddo
   
   do i=1, num_grp
  xs_empty=xs_dist-(num_grp-i)
  xs_first=min(num_grp, xs_empty+num_grp)
  xs_last=max(1,i-(xs_ihm-xs_ihs))
     if(IsUpSca .eq. 0 ) then
       read(xsunit,*,err=1002,end=1002) mat_info(k)%xs_p(j)%sig_a(i),mat_info(k)%xs_p(j)%mv_f(i),&
	  mat_info(k)%xs_p(j)%sig_t(i),(mat_info(k)%xs_p(j)%sig_gp2g(sgrp,i),sgrp=i,1,-1)         
     else
       read(xsunit,*,err=1002,end=1002) mat_info(k)%xs_p(j)%sig_a(i),mat_info(k)%xs_p(j)%mv_f(i),&
	  mat_info(k)%xs_p(j)%sig_t(i),(temp, sgrp=1,xs_empty),&
	  (mat_info(k)%xs_p(j)%sig_gp2g(sgrp,i),sgrp= xs_first,xs_last,-1) 

     endif 
   enddo
  enddo
 enddo

return
1002 write(LOGUNIT,"('xsmcnp: read pentran xs file error')")
stop "xsmcnp: read pentran xs file error"
end subroutine		


!write xs 
subroutine OutputXs
use mXslib
use mParaset

integer i, j, k
integer cmt,sgrp, fu_front, fu_end
integer fu_ihs, fu_iht, fu_ihm
character(LEN=120) temp, ioform

fu_iht=3
fu_ihm=2*num_grp-1+fu_iht
fu_ihs=fu_iht+num_grp

write(ioform, "('(',I0,'(ES12.5,x) )')" ) fu_ihm
open (unit=xsout, file='full.xs')

rewind(xsunit)
do k=1, num_mat
 do j=1, leg_order+1
   do cmt=1,num_cmt
     read(xsunit,"(A)",err=1002, end=1002)  temp   
     write(xsout,"(A)",err=1003)  temp
   enddo
       
   do i=1, num_grp
     fu_front=i-1
	 fu_end=fu_ihm-fu_iht-fu_front-num_grp
	 write(xsout, ioform,err=1003) &
      mat_info(k)%xs_p(j)%sig_a(i),mat_info(k)%xs_p(j)%mv_f(i),&
	  mat_info(k)%xs_p(j)%sig_t(i),(0.0, sgrp=1,fu_front),&
	  (mat_info(k)%xs_p(j)%sig_gp2g(sgrp,i),sgrp= num_grp,1,-1), &
	  (0.0, sgrp=1,fu_end) 
      read(xsunit,"(A)",err=1002, end=1002) temp
   
   enddo
  enddo
 enddo

close (xsout)
write(*,"('full scattering matrix file:', A,2x, 'iht=', I0, 2x, 'ihs=', I0, 2x, 'ihm=', I0)") &
    'full.xs', fu_iht, fu_ihs, fu_ihm
write(LOGUNIT,"('full scattering matrix file:', A,2x, 'iht=', I0, 2x, 'ihs=', I0, 2x, 'ihm=', I0)") &
    'full.xs', fu_iht, fu_ihs, fu_ihm
return
1002 write(LOGUNIT,"('xsmcnp: read pentran xs file error')")
stop "xsmcnp: read pentran xs file error"

1003 write(LOGUNIT,"('xsmcnp: write full scattering matrix file error: full.xs')")
stop "xsmcnp: write full scattering matrix file error: full.xs"
end subroutine		

!convert mnt to morse angles and weights (see morse code manual) 
subroutine ConvertMnt(mat_num, g_from, g_to,dflag)
use mXslib
use mParaset
use mPoly

integer, intent(in) :: mat_num, g_from, g_to
integer, intent(inout) :: dflag
integer i,j,k, iflag

double precision  areal, breal, term(3)


dflag=leg_mgxs
if(mat_num .eq. 3 .and. g_from .eq. 48 .and. g_to .eq. 48) then
    write (*,*) "test 01"
endif
! calculate M_mor
areal=fpl(0)
do k=0, leg_mgxs
    M_mor(k)=0
    fpl(k)=fpl(k)/areal
    do j=0, k
      M_mor(k)=M_mor(k)+(2*j+1)*pnl_rev(k,j)*fpl(j)/2
    enddo !j
enddo !k

!check if Gram matrix positive definite
!The Cholesky–Banachiewicz and Cholesky–Crout algorithm
!https://en.wikipedia.org/wiki/Cholesky_decomposition

if (mat_num .eq. 3 .and. g_from .eq. 14 .and. g_to .eq. 14) then
    write(LOGUNIT, "('test 03: ', 'mat=',I0, 2x, I0,'->', I0 )") mat_num, g_from, g_to
    write(LOGUNIT, "('M(0:3)= ', 4(f10.4, 2x) )")  M_mor(0:3)
    write(LOGUNIT, "('a(2:0)= ', 3(f10.4, 2x) )")  q_poly(ord_poly)%a_coeff(2:0:-1)
endif
do j=0, num_poly-1
  areal=M_mor(2*j)
  do k=0, j-1
    areal=areal-Lmat(j,k)**2
  enddo
  if(areal .lt. 0) exit
  Lmat(j,j)=dsqrt(areal)
  do i=j+1, num_poly-1
      breal=M_mor(j+i)
      do k=0, j-1
          breal=breal-Lmat(i,k)*Lmat(j,k)
      enddo
      Lmat(i,j)=breal/Lmat(j,j)
  enddo
enddo
ord_poly=num_poly

if(areal .lt. 0) then
    write(LOGUNIT,"('xs2mcnp: Gram Matrix is not positive definite: mat=',I0, 2x, I0,'->', I0 )") mat_num, g_from, g_to
    write(LOGUNIT,"('f(0:', I0,')=', 10(ES12.5, 2x))") leg_mgxs, fpl
    write(LOGUNIT,"('M(0:', I0,')=', 10(ES12.5, 2x))") leg_mgxs, M_mor
    
    if (j .ge. 1) then
       dflag=2*j-1
       ord_poly=j
       write(LOGUNIT,"('scattering downgraded to Pn=', I0 )") 2*j-1 
    else 
       dflag=0 
       write(LOGUNIT,"('scattering downgraded to isotropic' )")
    endif
endif

do i=0, ord_poly
    q_poly(i)%a_coeff(i)=1
    q_poly(i)%rank=i
enddo

if (mat_num .eq. 3 .and. g_from .eq. 14 .and. g_to .eq. 14) then
    write(LOGUNIT, "('test 03: ', 'mat=',I0, 2x, I0,'->', I0 )") mat_num, g_from, g_to
endif

do i=0, ord_poly-1
  N_mor(i)=0
  do k=0, i
    N_mor(i)=N_mor(i)+q_poly(i)%a_coeff(k)*M_mor(k+i)
  enddo
  if(N_mor(i) .le. 0) then
    write(LOGUNIT,"( 'N_mor(', I0,')=', ES12.5, ' for mat=',I0, 2x, I0,'->', I0 )" ) i, N_mor(i),mat_num, g_from, g_to
  endif
  L_mor(i+1)=0
  do k=0,i
   L_mor(i+1)=L_mor(i+1)+q_poly(i)%a_coeff(k)*M_mor(k+i+1)
  enddo
 
!  u_mor(i+1)=L_mor(i+1)/N_mor(i)-L_mor(i)/N_mor(i-1)
  if(i .eq. 0) then
   sig2_mor(i)=0
   if(N_mor(i) .ne. 0) then
    u_mor(i+1)=L_mor(i+1)/N_mor(i)
    ! u_mor(i+1)=fpl(1)
   else 
     u_mor(i+1)=0
   endif
  else
   if(N_mor(i-1) .ne. 0) then
     sig2_mor(i)=N_mor(i)/N_mor(i-1)
   else
    sig2_mor(i)=0
    u_mor(i+1)=0
   endif
   if(N_mor(i) .ne. 0) then
    u_mor(i+1)=L_mor(i+1)/N_mor(i)-L_mor(i)/N_mor(i-1)
   else 
     u_mor(i+1)=0  
   endif
   endif
  
  q_poly(i+1)%a_coeff(i+1)=1
  do j=0, i
      if(j .eq. 0) then
        term(1)=0
       else
        term(1)=q_poly(i)%a_coeff(j-1)
      endif
      term(2)=-u_mor(i+1)*q_poly(i)%a_coeff(j)
      if(i .eq. 0 .or. j .gt. i-1) then
       term(3)=0
      else
       term(3)=-sig2_mor(i)*q_poly(i-1)%a_coeff(j)
      endif

      q_poly(i+1)%a_coeff(j)=term(1)+term(2)+term(3)
  enddo !j

enddo  !i 

call poly_roots(mat_num, g_from, g_to,dflag)
if (mat_num .eq. 3 .and. g_from .eq. 14 .and. g_to .eq. 14) then
    write(LOGUNIT, "('test 02: ', 'mat=',I0, 2x, I0,'->', I0 )") mat_num, g_from, g_to
    write(LOGUNIT, "('M(0:3)= ', 4(f10.4, 2x) )")  M_mor(0:3)
    write(LOGUNIT, "('a(2:0)= ', 3(f10.4, 2x) )")  q_poly(ord_poly)%a_coeff(2:0:-1)
endif
    

term(2)=0
do i=1, ord_poly !num of weight

    term(1)=0
    areal=q_poly(ord_poly)%root(i)

    do k=0, ord_poly-1 !q summation
        breal=0             !q value
        do j=0, k
            breal=breal+q_poly(k)%a_coeff(j)* (areal**j)
        enddo     !j

        if(N_mor(k) .ne. 0) then
            term(1)=term(1)+breal*breal/N_mor(k)
        endif
    enddo  !k

    if( term(1) .ne. 0) then
        q_poly(ord_poly)%wgt(i)=1.0/term(1)
    else
        q_poly(ord_poly)%wgt(i)=0
    endif

    term(2)=term(2)+q_poly(ord_poly)%wgt(i)
enddo  !i

!fpl output
!weight
!if downgraded, put zero weight
j=num_poly-ord_poly
do i=1, j
    fpl(i)=0
enddo
k=0
do i=j+1, num_poly-1
    k=k+1
    if(term(2) .ne. 0) then
      areal=q_poly(ord_poly)%wgt(k)/term(2)
    else
      areal=0
    endif
    if( i .eq. 1) then
       fpl(i)=areal
    else 
       fpl(i)=fpl(i-1)+areal 
    endif
!    if(fpl(i) .lt. 0 .or. fpl(i) .gt. 1) then
!        write (LOGUNIT,"('fpl(', I0, ')=', ES12.5,' for mat=',I0, 2x, I0,'->', I0 )" ) i, fpl(i), mat_num, g_from, g_to
!    endif
enddo

!u
do i=num_poly, num_poly+j-1
    fpl(i)=0
enddo
k=0
do i=num_poly+j, leg_mgxs
 k=k+1   
 fpl(i)=q_poly(num_poly)%root(k)
enddo

if(dflag .gt. leg_mgxs) then
    write(LOGUNIT,"('fpl(1:', I0,')=', 10(ES12.5, 2x))") leg_mgxs, fpl(1:leg_mgxs)
    write(LOGUNIT,"('**********************************')") 
endif

end subroutine

!finding polynomial roots
subroutine poly_roots(mat_num, g_from, g_to,dflag)
use mPoly
use mParaset

integer, intent(in) :: mat_num, g_from, g_to
integer, intent(inout) :: dflag

integer i,k, iter, max_iter
double precision  z,s
double precision :: tol=1E-9, delta
double precision :: q_val, qp_val  !poly value and p' value

if(allocated(q_poly(ord_poly)%root)) deallocate(q_poly(ord_poly)%root)
if(allocated(q_poly(ord_poly)%wgt)) deallocate(q_poly(ord_poly)%wgt)

allocate(q_poly(ord_poly)%root(q_poly(ord_poly)%rank)) 
allocate(q_poly(ord_poly)%wgt(q_poly(ord_poly)%rank)) 

delta=1.0
max_iter=10000
z=-1.0
do k=1, q_poly(ord_poly)%rank
! z=sin(phi+(k-1)*pi/(2.0*sn_2))  !**2
! if(k .ne. 1) z=(z+mu_set(k-1))/2
 delta=1.0
 z=z+tol*1000
 iter=0
 do while(abs(delta) .gt. tol .and. iter .le. max_iter)
  iter=iter+1
  s=0
  do i=1, k-1
   s=s+1.0/(z-q_poly(ord_poly)%root(i) )
  enddo
  !x=sqrt(z)
    
  q_val=0
  do i=0, q_poly(ord_poly)%rank
    q_val=q_val+z**i*q_poly(ord_poly)%a_coeff(i)
    
  enddo
  qp_val=0
  do i=1, q_poly(ord_poly)%rank
  qp_val=qp_val+z**(i-1) * i*q_poly(ord_poly)%a_coeff(i)
  enddo
  
  delta=-q_val/(qp_val-q_val*s)
 
  z=z+delta

 enddo
 q_poly(ord_poly)%root(k)=z

end do 

do k=1, q_poly(ord_poly)%rank
  if(dabs(q_poly(ord_poly)%root(k)) .gt. 1) then
      write(LOGUNIT,"('xs2mcnp: root out of range [-1,1]: mat=',I0, 2x, I0,'->', I0, ' root=', ES12.5 )") mat_num, g_from, g_to, q_poly(ord_poly)%root(k)
  endif
  if(q_poly(ord_poly)%root(k) .gt. 1) then
    q_poly(ord_poly)%root(k)=1.0
  else if(q_poly(ord_poly)%root(k) .lt. -1) then
    q_poly(ord_poly)%root(k)=-1.0
  endif
enddo


return
end subroutine

!Parsing a line
!cmtflag=0: to tell if the line is composed numbers only,
!           yes: itm=non-zero
!           no:  itm=0
!cmtflag=1:  to count how many numbers in the line
!            itm=num of numbers
!            num_word=number of words if present
! if it is a blank line, itm=-1 num_word=-1
! if it is a comment line, itm=0 num_word=0
subroutine ParseLine(cmtflag)
use mLineParse


integer, intent(in) :: cmtflag

integer i,j
integer buf_len
integer IsNum, IsSciFormat, IsWord
integer num_word, num_val

character(len=1) :: mychar=''

line_buf= adjustl(line_buf)
buf_len=len_trim(line_buf)

!blank line
if(buf_len .eq. 0) then
itm=-1
num_word=-1
return
endif

!comment line
if(index (cmtchar_end,line_buf(1:1)) .ne. 0) then 
  itm=0
  num_word=0
  return
endif

!remove comment at the end of line
do j=1, buf_len
 if(index(cmtchar_end,line_buf(j:j)) .ne. 0) then
  buf_len=j-1
  exit
 endif
enddo


!remove prevailing keyword
do j=1, min(buf_len-1,len_keyword)
  if( line_buf(j:j) .eq. '=') then
  do i=1, j
  line_buf(i:i)=' '
  enddo
  line_buf=adjustl(line_buf)
  buf_len=len_trim(line_buf)
  exit
 endif
enddo

itm=0
iwd=0
IsNum=0
IsSciFormat=0
IsWord=0
num_val=0
num_word=0

do j=1, buf_len

mychar=line_buf(j:j)
if (index(' ,', mychar ) .ne. 0) then
  if(IsWord .eq. 0 .and. IsNum .eq. 1) num_val=num_val+1 
  if(IsWord .eq. 1 ) num_word=num_word+1
  IsWord=0
  IsNum=0
  IsSciFormat=0
elseif(IsWord .eq. 1) then
 cycle
elseif(index(digchar, mychar) .ne. 0  .or. &
  ( index(scichar, mychar) .ne. 0 .and. IsSciFormat .eq. 0)   ) then
   IsNum=1
   if( index(scichar, mychar) .ne. 0) IsSciFormat=1
   cycle
else
  IsWord=1
  if(cmtflag .eq. 0 ) then
    itm=0
	iwd=1
    return
  endif	
endif

enddo  !j

!handle if last entry is a number
if(IsWord .eq. 0 .and. IsNum .eq. 1)  num_val=num_val+1
if(IsWord .eq. 1 .and. IsNum .eq. 0)  num_word=num_word+1


itm=num_val
iwd=num_word
return

end subroutine


Subroutine mcnp2xs
use mParaset
use mLineParse
use mXslib
use mPoly
use IEEE_ARITHMETIC

integer :: nxs_block(16)=1
integer :: jxs_block(32)=1
real*8, allocatable :: data_block(:)
logical ex
character(len=1):: mclib_letter(4)
character(len=12) ::  matname
integer i,j,k, ipos, jpos,kpos, m
integer  g_from, g_to, gto(2)
integer  :: cmt_err=0

real*8 :: cos_set(35)=0, pn_set(35)=0, mval,wval   !cosine bin boundaries up to 35 bins
ex=.FALSE.
mclib_letter(1)='m'
mclib_letter(2)='n'
mclib_letter(3)='g'
mclib_letter(4)='.'
num_mat=0
num_mat_p=0
num_mat_n=0

!counting lines
num_line=0
cmt_err=0

do while(cmt_err .eq. 0) 
num_line=num_line+1
read (xsunit, "(A)", iostat=cmt_err) line_buf
if(cmt_err .ne. 0) exit

if(len_trim(line_buf) .ge. leg_line) then
write(*, "('number of columns >' , I0)" ) leg_line
stop 'XS conversion failed'
endif 

line_buf= adjustl(line_buf)
matname=line_buf(1:12)
i=index(matname,mclib_letter(1) )
j=index(matname,mclib_letter(3) )
! i=i+j
k=index(matname,mclib_letter(4) )
if( i .ne. 0 .and. i-k .eq. 3) then
   num_mat_n=num_mat_n+1
else if(j .ne. 0 .and. j-k .eq. 3) then
   num_mat_p=num_mat_p+1
endif
enddo
num_mat=num_mat_n+num_mat_p


rewind(xsunit)

if( num_mat .eq. 0) then
write(*, "('no xs section found in ' , A)" ) trim(xsfile)
write(*, "('each mat section should start with a zzzaaa.xxm format, e.g. 92000.22m ' )" ) 
stop 'XS conversion failed'
else

write(*,  "(I0 ,' materials found in multigroup MCNP xs file: ' , A)" ) num_mat, trim(xsfile)
write(LOGUNIT,  "(I0 ,' materials found in multigroup MCNP xs file: ' , A)" ) num_mat, trim(xsfile)
if(num_mat_n .gt. 0 .and. num_mat_p .gt. 0) then
write(*, "('neutron gamma coupled library detected (n/g): ', I0, '/', I0 )" )  num_mat_n, num_mat_p
write(LOGUNIT, "('neutron gamma coupled library detected (n/g): ', I0, '/', I0 )" )  num_mat_n, num_mat_p
endif

endif

allocate (mat_info(num_mat))
!for now only P0 and P1 can be converted
leg_order=0
!reading lines
m=0
cmt_err=0
do while(cmt_err .eq. 0) 

read (xsunit, "(A)",  err=1002,end=1002) line_buf
num_grp_i=0

line_buf= adjustl(line_buf)
matname=line_buf(1:12)
i=index(matname,mclib_letter(1) )
j=index(matname,mclib_letter(3) )

k=index(matname,mclib_letter(4) )
if( (i+j)*k .ne. 0 .and. i+j-k .eq. 3) then
   m=m+1
   read(line_buf,*,iostat=ipos) mat_info(m)%name, mat_info(m)%mass
   if(ipos .ne. 0)  read(line_buf,*) mat_info(m)%name
   curfile=mat_info(m)%name
   do i=1,5
     read(xsunit, *,  err=1002,end=1002)
   enddo
   read(xsunit, *,  err=1002,end=1002) nxs_block
   mat_info(m)%nxs=nxs_block
   mat_info(m)%z_num=int(nxs_block(2)/1000)
   mat_info(m)%a_num=nxs_block(2)-mat_info(m)%z_num*1000

!check Pn order
!   if( nxs_block(3) .ne. 0 .and. nxs_block(3) .ne. 1) then
!	 write(*,"('Pn order located at nxs(3) is greater than 1 for material section : ', I0)") m 
!	 write(*,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
!   endif   
   
!   if( m .eq. 1) then
!    num_grp=nxs_block(5)
!    if( num_grp .eq. 0) then
!	 write(*,"('number of groups=0 for: ',A)") trim(curfile) 
!	 stop 'mcnp2xs conversion failed'
!	endif 
!   else
!    if(nxs_block(5) .ne. num_grp) then
!	  write(*,"('number of groups doesnot match other mat section: ',A)") trim(curfile) 
!	  stop 'mcnp2xs conversion failed'
!	endif
!  endif 

 !if(nxs_block(6) .ne. num_grp .or. nxs_block(7) .ne. num_grp) then
!	  write(*,"('partial scattering matrix format not supported for mat section: ',A)") trim(curfile) 
!	  stop 'mcnp2xs conversion failed'
!	endif
 read(xsunit, *,  err=1002,end=1002) jxs_block
 allocate ( data_block (nxs_block(1) )  )  
 read(xsunit, *,  err=1002,end=1002) data_block

 if(nxs_block(12) .eq. 2) then 
 if(p_atom(mat_info(m)%z_num) .eq. 0) then
 p_atom(mat_info(m)%z_num)=m
 else
  write(*,"('photoatomic xs for atom number ',I0, ' is already defined before: ', A)") mat_info(m)%z_num, trim(curfile) 
  write(*,  "('  will try to use photoatomic data in section: ', A )" ) trim(curfile)

  write(LOGUNIT,"('photoatomic xs for atom number ',I0, ' is already defined before: ', A)") mat_info(m)%z_num, trim(curfile) 
  write(LOGUNIT,  "('  will try to use photoatomic data in section: ', A )" ) trim(curfile)
 endif
 endif
 
 if(nxs_block(3) .le. 1 .or. nxs_block(9) .eq. 0) then
  mat_info(m)%leg_n=nxs_block(3)
 else 
  mat_info(m)%leg_n=0
  write(*,"('Higher Pn order xs in discrte format: nxs(9)=1 in material section : ', A)") trim(curfile) 
  write(*,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
  write(LOGUNIT,"('Higher Pn order xs in discrte format: nxs(9)=1 in material section : ', A)") trim(curfile) 
  write(LOGUNIT,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
 endif

 mat_info(m)%leg_p=0
 if(nxs_block(8) .ne. 0) then  !for secondary particle
 if(int(data_block(jxs_block(15))) .le. 1 .or. int(data_block(jxs_block(14))) .eq. 0) then
 mat_info(m)%leg_p=int(data_block(jxs_block(15)))
 else
  mat_info(m)%leg_p=0
  write(*,"('Higher Pn order for secondary xs in discrte format: jxs(14)=1 in material section : ', A)") trim(curfile) 
  write(*,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
  write(LOGUNIT,"('Higher Pn order for secondary xs in discrte format: jxs(14)=1 in material section : ', A)") trim(curfile) 
  write(LOGUNIT,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
 endif
 endif
 mat_info(m)%leg_m=max(mat_info(m)%leg_n, mat_info(m)%leg_p)

 if (mat_info(m)%leg_m .gt. 1 .and. mod(mat_info(m)%leg_m,2) .eq. 0) then
    mat_info(m)%leg_m=mat_info(m)%leg_m+1   !leg_order =num_angle=1
 endif
 
 if(mat_info(m)%leg_m .gt. leg_order) leg_order=mat_info(m)%leg_m 
 
 num_grp=nxs_block(5)
 num_grp_i=num_grp
 if(nxs_block(8) .ne. 0) then
   i=jxs_block(11)
   if( i .ne. 0) i=int(data_block(i))
   if (i .ne. 0) mat_info(m)%sec_type=i
   
   i=jxs_block(12)
   if( i .ne. 0) i=int(data_block(i))
   if (i .ne. 0) i=int(data_block(i))
   
   num_grp=num_grp+i
   if( i .ne. 0) mat_info(m)%sec_ngrp=i

   if(nxs_block(8) .gt. 1) then
    write(*,"('More than one secondary particle present in: ',A)") trim(curfile) 
    write(*,"('only the first secondary particle will be processed: ')")  
	write(LOGUNIT,"('More than one secondary particle present in: ',A)") trim(curfile) 
    write(LOGUNIT,"('only the first secondary particle will be processed: ')")  
  endif
 endif
 
 
 allocate(mat_info(m)%xs_p(mat_info(m)%leg_m+1))

 allocate(mat_info(m)%sig_s(num_grp),&
	 mat_info(m)%chi(num_grp),&
	 mat_info(m)%nubar(num_grp),& 
	 mat_info(m)%sig_f(num_grp))
allocate (mat_info(m)%upper_bnd(num_grp))
do j=1, mat_info(m)%leg_m+1
allocate(mat_info(m)%xs_p(j)%sig_a(num_grp),&
         mat_info(m)%xs_p(j)%mv_f(num_grp),&
	 mat_info(m)%xs_p(j)%sig_t(num_grp))

allocate(mat_info(m)%xs_p(j)%sig_gp2g(num_grp,num_grp))

mat_info(m)%xs_p(j)%sig_a=0
mat_info(m)%xs_p(j)%mv_f=0
mat_info(m)%xs_p(j)%sig_t=0
mat_info(m)%xs_p(j)%sig_gp2g=0
enddo  ! leg_order

mat_info(m)%chi=0
mat_info(m)%nubar=0
mat_info(m)%sig_f=0

do i=1, nxs_block(5)
mat_info(m)%upper_bnd(i)=data_block(i)+data_block(i+nxs_block(5))/2
enddo

do i=nxs_block(5)+1, num_grp
j=int(data_block(jxs_block(12)))+i-nxs_block(5)
mat_info(m)%upper_bnd(i)=data_block(j)+data_block(j+num_grp-nxs_block(5))/2
enddo

i=jxs_block(2)
mat_info(m)%xs_p(1)%sig_t(1:num_grp_i)=data_block(i:i+num_grp_i-1)
i=jxs_block(3)
if (i .ne. 0) mat_info(m)%sig_f(1:num_grp_i)=data_block(i:i+num_grp_i-1) 
i=jxs_block(4)
if (i .ne. 0) mat_info(m)%nubar(1:num_grp_i)=data_block(i:i+num_grp_i-1)
i=jxs_block(5)
if (i .ne. 0) mat_info(m)%chi(1:num_grp_i)=data_block(i:i+num_grp_i-1)
i=jxs_block(6)
if (i .ne. 0) then
 if(jxs_block(3) .ne. 0) then
  mat_info(m)%xs_p(1)%sig_a(1:num_grp_i)=data_block(i:i+num_grp_i-1)+&
   mat_info(m)%sig_f(1:num_grp_i)
 else
  mat_info(m)%xs_p(1)%sig_a(1:num_grp_i)=data_block(i:i+num_grp_i-1)
 endif
endif


if(mat_info(m)%leg_m .gt. 0 .and. jxs_block(16) .eq. 0) then
  write(*,"('MCNP xs format conflict: nxs(3) >0 and jxs(16)=0, for mat section ', I0)") m 
  stop 'mcnp2xs conversion failed'
endif

if(mat_info(m)%leg_m .gt. 0 .and. jxs_block(17) .eq. 0) then
  write(*,"('MCNP xs format conflict: nxs(3) >0 and jxs(17)=0, for mat section ', I0)") m 
  stop 'mcnp2xs conversion failed'
endif
!P0 data primary
i=jxs_block(13)
i=int(data_block(i))
do g_from=1, num_grp_i
gto(1)=max(1,g_from-nxs_block(6))
gto(2)=min(nxs_block(5), g_from+nxs_block(7))
do g_to=gto(1), gto(2)
if(data_block(i) .lt. 0) then 
write(LOGUNIT,"('P0 scattering xs <0 in ', A, ' Grp ', I0, ' to ', I0, ' xs=',ES12.5, ' at location ',I0, ' **set to zero**')") &
   trim(mat_info(m)%name) , g_from, g_to,data_block(i), i
data_block(i)=0
endif
mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)=data_block(i)
i=i+1
enddo
enddo

!primary to secondary P0
if(nxs_block(8) .gt. 0) then

i=jxs_block(13)+1
i=int(data_block(i))
do g_from=1, num_grp_i
do g_to=num_grp_i+1, num_grp
if(data_block(i) .lt. 0) then 
write(LOGUNIT,"('P0 scattering xs <0 in ', A, ' Grp ', I0, ' to ', I0, ' xs=',ES12.5, ' at location ',I0, ' **set to zero**')") &
  trim(mat_info(m)%name) , g_from, g_to,data_block(i), i
data_block(i)=0
endif

mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)=data_block(i)
i=i+1
enddo
enddo
endif

!Pn data primary
if( mat_info(m)%leg_n .gt. 0) then
ipos=jxs_block(17)
ipos=int(data_block(ipos))
if(ipos .le. 0) then
write(*,"('Higher Pn order block absent: data(jxs(17))=0 while nxs(3)>0 in material section : ', A)") trim(curfile) 
write(*,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
write(LOGUNIT,"('Higher Pn order block absent: data(jxs(17))=0 while nxs(3)>0 in material section : ', A)") trim(curfile) 
write(LOGUNIT,  "('  Warning : Only P0 scattering data will be converted to PENTRAN format ' )" )
else

i=0
do g_from=1, num_grp_i
gto(1)=max(1,g_from-nxs_block(6))
gto(2)=min(nxs_block(5), g_from+nxs_block(7))
do g_to=gto(1), gto(2)

kpos=int(data_block(int(data_block(jxs_block(16))) +i))
i=i+1
!   if(mat_info(m)%z_num .eq. 63 .and. mat_info(m)%a_num .eq. 0 .and. g_from .eq. 1 .and. g_to .eq. 21) then
!   write(*,*) mat_info(m)%z_num
!   endif
if (kpos .le. 0)  cycle

jpos=ipos+kpos-1
cos_set=0
pn_set=0
do k=1, nxs_block(3)
j=jpos+k-1
cos_set(k)=data_block(j)
enddo

if(mat_info(m)%leg_m .eq. 1) then
mat_info(m)%xs_p(2)%sig_gp2g(g_from,g_to)=mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)*cos_set(1)
cycle
endif

do k=1, nxs_block(3)-1
if(cos_set(k+1)-cos_set(k) .lt. MySmallestNumber) then
write(LOGUNIT,"('zero or negative size cosine bin ', A, ' Grp ', I0, ' to ', I0,  ' at location ',I0, '  **ignored**' )") &
   trim(mat_info(m)%name) , g_from, g_to, jpos+k-1
endif
enddo


do j=2, mat_info(m)%leg_m+1
mval=0

do k=1, nxs_block(3)
pn_set(k)=dlegendre(cos_set(k),j-1)*(1-cos_set(k)**2)/j/(j-1)
enddo
do k=1, nxs_block(3)-1

if(cos_set(k+1)-cos_set(k) .lt. MySmallestNumber) then
wval=0
else
wval=1.0/(nxs_block(3)-1)/(cos_set(k+1)-cos_set(k))
endif
mval=mval+(pn_set(k)-pn_set(k+1))*wval
enddo
mat_info(m)%xs_p(j)%sig_gp2g(g_from,g_to)=mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)*mval

enddo  !j


enddo  !g_to
enddo  !g_from

endif !ipos
endif


!Pn data secondary
if( mat_info(m)%leg_p .gt. 0) then
ipos=jxs_block(17)+1
ipos=int(data_block(ipos))
if(ipos .ne. 0) then
i=0
do g_from=1, num_grp_i
gto(1)=max(1,g_from-nxs_block(6))
gto(2)=min(nxs_block(5), g_from+nxs_block(7))
do g_to=gto(1), gto(2)

kpos=int(data_block(int(data_block(jxs_block(16))) +i))
i=i+1

if (kpos .le. 0)  cycle

jpos=ipos+kpos-1
cos_set=0
pn_set=0
do k=1, nxs_block(3)
j=jpos+k-1
cos_set(k)=data_block(j)
enddo

if(mat_info(m)%leg_m .eq. 1) then
mat_info(m)%xs_p(2)%sig_gp2g(g_from,g_to)=mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)*cos_set(1)
cycle
endif

do k=1, nxs_block(3)-1
if(cos_set(k+1)-cos_set(k) .lt. MySmallestNumber) then
write(LOGUNIT,"('zero or negative size cosine bin ', A, ' Grp ', I0, ' to ', I0,  ' at location ',I0,'  **ignored**')") &
   trim(mat_info(m)%name) , g_from, g_to, jpos+k-1
endif

enddo

do j=2, mat_info(m)%leg_m+1

mval=0
do k=1, nxs_block(3)
pn_set(k)=dlegendre(cos_set(k),j-1)*(1-cos_set(k)**2)/j/(j-1)
enddo
do k=1, nxs_block(3)-1
if(cos_set(k+1)-cos_set(k) .lt. MySmallestNumber) then
wval=0
else
wval=1.0/(nxs_block(3)-1)/(cos_set(k+1)-cos_set(k))
endif

mval=mval+(pn_set(k)-pn_set(k+1))*wval
enddo
mat_info(m)%xs_p(j)%sig_gp2g(g_from,g_to)=mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)*mval

enddo  !j


enddo  !g_to
enddo  !g_from

endif !ipos
endif

!if(mat_info(m)%leg_m .eq. 1) then
!i=jxs_block(17)
!i=int(data_block(i))
!!for P1 only 
!!mu_average=p1/p0
!do g_from=1, num_grp
!do g_to=1, num_grp

!mat_info(m)%xs_p(2)%sig_gp2g(g_from,g_to)=mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)*data_block(i)
!i=i+1
!enddo
!enddo

!endif

deallocate (data_block)
endif !end mat block

if(m .eq. num_mat) cmt_err=1
enddo  !while
close(xsunit)

Call CheckData
call OutputXs_m2x

return
1001  write(LOGUNIT,"('open mcnp xs file error')")
stop "open mcnp xs file error"
write(LOGUNIT,"('read file error in material section :', A)") trim(curfile)
1002 write(*,"('read file error in material section :', A)") trim(curfile) 
end subroutine

subroutine CheckData
use mLineParse
use mXslib
use mParaset


integer :: nxs_block(16)=1
integer :: jxs_block(32)=1

integer i,j,k, ipos, jpos,kpos, m
integer  g_from, g_to, gto(2), m_p, gb(2)
real :: mval=0

num_grp_n=0
num_grp_p=0
write(LOGUNIT,*) '--------------checking data warnings---------------------'

do m=1, num_mat
mat_info(m)%IsSkip=0
!primary group number check
if( mat_info(m)%nxs(12) .eq. 1) then

if(num_grp_n .eq. 0) then 
  num_grp_n=mat_info(m)%nxs(5)
elseif(num_grp_n .ne. mat_info(m)%nxs(5)) then
 write(*,"('number of neutron group inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
 write(LOGUNIT,"('number of neutron group inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name

 mat_info(m)%IsSkip=mat_info(m)%IsSkip+1
endif

else if( mat_info(m)%nxs(12) .eq. 2) then

if(num_grp_p .eq. 0) then 
  num_grp_p=mat_info(m)%nxs(5)
elseif(num_grp_p .ne. mat_info(m)%nxs(5)) then
 write(*,"('number of photon group inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
 write(LOGUNIT,"('number of photon group inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
 mat_info(m)%IsSkip=mat_info(m)%IsSkip+1
endif
endif  !nxs(12)


!secondary group number check
if ( mat_info(m)%nxs(8) .ge. 1) then
if( mat_info(m)%sec_type .eq. 1) then
 write(*,"('secondary particle is set to neutron in : ', A)") mat_info(m)%name
 write(LOGUNIT,"('secondary particle is set to neutron in : ', A)") mat_info(m)%name
 if(num_grp_n .eq. 0) then 
  num_grp_n=mat_info(m)%sec_ngrp
 elseif(num_grp_n .ne. mat_info(m)%sec_ngrp) then
  write(*,"('number of neutron group inconsistant in : ', A)") mat_info(m)%name
  write(*,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('number of neutron group inconsistant in : ', A)") mat_info(m)%name
  write(LOGUNIT,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
  mat_info(m)%IsSkip=mat_info(m)%IsSkip+1
 endif 
else if  ( mat_info(m)%sec_type .eq. 2) then
 
if(num_grp_p .eq. 0) then 
  num_grp_p=mat_info(m)%sec_ngrp
elseif(num_grp_p .ne. mat_info(m)%sec_ngrp) then
 write(*,"('number of photon group inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('number of photon group inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : will be skpped : ', A )" ) mat_info(m)%name
 mat_info(m)%IsSkip=mat_info(m)%IsSkip+1
endif

endif !sec_type
endif !nxs(8) number of sec. particle
enddo !m
write(LOGUNIT,  "('  primary particle number of group : checked ' )" ) 
write(LOGUNIT,  "('  secondary particle number of group : checked ' )" ) 
 
if (mat_info(m)%nxs(12) .ne. 0) num_grp=num_grp_n+num_grp_p

num_mat_r=num_mat
allocate(n_map(num_mat_r))
do m=1,num_mat_r
n_map(m)=m
enddo

!group boundary check
if(num_grp_n .gt. 0)  then
 allocate(upper_bnd_n(num_grp_n))
 upper_bnd_n=0
endif

if(num_grp_p .gt. 0)  then
 allocate(upper_bnd_p(num_grp_p))
 upper_bnd_p=0
endif

do m=1, num_mat

if(mat_info(m)%IsSkip .ne. 0) cycle
!primary (incident particle)
if( mat_info(m)%nxs(12) .eq. 1) then

do i=1, num_grp_n
if(upper_bnd_n(i) .le. MySmallestNumber) then 
  upper_bnd_n(i)=mat_info(m)%upper_bnd(i)
elseif(abs(upper_bnd_n(i) - mat_info(m)%upper_bnd(i)) .gt. MySmallestNumber) then
 write(*,"('neutron energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('neutron energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
endif
enddo !i

elseif( mat_info(m)%nxs(12) .eq. 2) then !primary particle is photon

do i=1, num_grp_p
if(upper_bnd_p(i) .le. MySmallestNumber) then 
  upper_bnd_p(i)=mat_info(m)%upper_bnd(i)
elseif(abs(upper_bnd_p(i) - mat_info(m)%upper_bnd(i)) .gt. MySmallestNumber) then
 write(*,"('photon energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('photon energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
endif
enddo !i

endif !nxs(12)



if ( mat_info(m)%nxs(8) .ge. 1) then  !secondary particle present
if( mat_info(m)%sec_type .eq. 1) then !secondary particle is neutron

do i=1, num_grp_n
j=i+mat_info(m)%nxs(5)
if(upper_bnd_n(i) .le. MySmallestNumber) then 
  upper_bnd_n(i)=mat_info(m)%upper_bnd(j)
elseif(abs(upper_bnd_n(i) - mat_info(m)%upper_bnd(j)) .gt. MySmallestNumber) then
 write(*,"('neutron energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('neutron energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
endif
enddo !i

else if( mat_info(m)%sec_type .eq. 2) then
do i=1, num_grp_p
j=i+mat_info(m)%nxs(5)
if(upper_bnd_p(i) .le. MySmallestNumber) then 
  upper_bnd_p(i)=mat_info(m)%upper_bnd(j)
elseif(abs(upper_bnd_p(i) - mat_info(m)%upper_bnd(j)) .gt. MySmallestNumber) then
 write(*,"('photon energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(*,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('photon energy boundaries inconsistant in : ', A)") mat_info(m)%name
 write(LOGUNIT,  "('  Warning : boundary values ignored : ', A )" ) mat_info(m)%name
endif
enddo !i
endif !sec_type
endif !nxs(8) number of sec. particle

enddo !m
write(LOGUNIT,  "('  primary pariticle energy group boundary: checked ' )" ) 
write(LOGUNIT,  "('  secondary pariticle energy group boundary: checked ' )" ) 

if(num_grp_p .eq. 0 .or. num_grp_n .eq. 0) return

num_mat_r=0
n_map=0
do m=1, num_mat
if(mat_info(m)%IsSkip .ne. 0) cycle
if(mat_info(m)%nxs(12) .eq. 2) cycle
if(mat_info(m)%nxs(12) .eq. 1 .and. mat_info(m)%nxs(8) .eq. 0) cycle
num_mat_r=num_mat_r+1
n_map(num_mat_r)=m
m_p=p_atom(mat_info(m)%z_num)
if (m_p .eq. 0) then
  write(*,"('photoatompic xs data not found for : ', A)") mat_info(m)%name
  write(*,  "('  Warning : photon scattering xs will be set to zero : ', A )" ) mat_info(m)%name
    write(LOGUNIT,"('photoatompic xs data not found for : ', A)") mat_info(m)%name
  write(LOGUNIT,  "('  Warning : photon scattering xs will be set to zero : ', A )" ) mat_info(m)%name
else if(mat_info(m_p)%IsSkip .ne. 0) then
m_p=0
do i=1, num_mat
if(mat_info(i)%nxs(12) .eq. 2 .and. mat_info(i)%z_num .eq. mat_info(m)%z_num .and. mat_info(i)%IsSkip .eq. 0) then
m_p=i
exit
endif
enddo !i
if (m_p .eq. 0) then
  write(*,"('photoatompic xs data found for : ', A, ' but number of group not match')") mat_info(m)%name
  write(*,  "('  Warning : photon scattering xs will be set to zero : ', A )" ) mat_info(m)%name
  write(LOGUNIT,"('photoatompic xs data found for : ', A, ' but number of group not match')") mat_info(m)%name
  write(LOGUNIT,  "('  Warning : photon scattering xs will be set to zero : ', A )" ) mat_info(m)%name
endif !m_p still .eq. 0 

endif !m_p .eq. 0

!add photoatomic data to neutron data where primary is photon
do k=1, mat_info(m_p)%leg_n+1
do g_from=1, num_grp_p
gto(1)=max(1,g_from-mat_info(m_p)%nxs(6))
gto(2)=min(mat_info(m_p)%nxs(5), g_from+mat_info(m_p)%nxs(7))
gb(1)=g_from+num_grp_n
if(k .eq.1) then
mat_info(m)%xs_p(k)%sig_t(gb(1))=mat_info(m_p)%xs_p(k)%sig_t(g_from)
mat_info(m)%xs_p(k)%sig_a(gb(1))=mat_info(m_p)%xs_p(k)%sig_a(g_from)
endif
do g_to=gto(1), gto(2)
gb(2)=g_to+num_grp_n
mat_info(m)%xs_p(k)%sig_gp2g(gb(1),gb(2))=mat_info(m_p)%xs_p(k)%sig_gp2g(g_from,g_to)
enddo  !go_to

enddo  !go_from

enddo !k leg_n

enddo  !m
write(LOGUNIT,  "('  photoatomic data match with neutron library: checked ' )" ) 
return
end subroutine

!output xs
subroutine OutputXs_m2x
use mLineParse
use mXslib
use mParaset

integer i, j, k
integer g_from, g_to, m, order_uni,mp,mi
integer cmt,sgrp, fu_front, fu_end
integer fu_ihs, fu_iht, fu_ihm
character(LEN=120) temp, ioform
real*8 :: tiny_num=10E-16, sig



call date_and_time (str_date, str_time)
write(str_date_mg,"(A2,'/',A2,'/', A4)")  str_date(5:6),str_date(7:8), str_date(1:4)


IsUpSca=0
m_loop : do mp=1, num_mat_r
m=n_map(mp)
do g_from=1, num_grp
do g_to=1, g_from-1
sig=mat_info(m)%xs_p(1)%sig_gp2g(g_from,g_to)
if(sig .gt. 0) then
 IsUpSca=1
 exit m_loop
endif
enddo
enddo
enddo  m_loop !m

do mp=1, num_mat_r
m=n_map(mp)
do i=1, num_grp
mat_info(m)%xs_p(1)%mv_f(i)=mat_info(m)%sig_f(i)*mat_info(m)%nubar(i)
enddo
enddo !m

fu_iht=3
if(IsUpSca .eq. 1) then
fu_ihm=2*num_grp-1+fu_iht
fu_ihs=fu_iht+num_grp
else
fu_ihm=num_grp+fu_iht
fu_ihs=fu_iht+1
endif

write(ioform, "('(',I0,'(ES12.5,x) )')" ) fu_ihm
curfile=trim(prbname)//'.xs'
open (unit=xsunit, file=trim(curfile) )
 
if(leg_order .gt. 1) then
 order_uni=0
else
 order_uni=leg_order
endif

order_uni=leg_order
 
do mp=1, num_mat_r
 k=n_map(mp)
! order_uni=mat_info(k)%leg_m
 do j=1, order_uni+1
  write(xsunit,"(A, ' G=', I0, ' iht=', I0, ' ihs=',I0, ' ihm=',I0, '  p', I0, ' data  from ', A)",err=1003) &
   trim(mat_info(k)%name), num_grp,fu_iht, fu_ihs, fu_ihm, j-1, trim(prbname)
  write(xsunit,"('#',13X, A10,15x,'ampx', 8x, I2.2,I3.3)",err=1003) str_date_mg, mat_info(k)%z_num, mat_info(k)%a_num 
  if(j .gt. mat_info(k)%leg_m+1) then
     do i=1, num_grp
     write(xsunit, ioform, err=1003) (0.0, sgrp=1, fu_ihm)
     enddo 
    cycle
  endif 
  if(IsUpsca .eq. 1) then   
   do i=1, num_grp
     fu_front=i-1
     fu_end=fu_ihm-fu_iht-fu_front-num_grp
     write(xsunit, ioform, err=1003) &
      mat_info(k)%xs_p(j)%sig_a(i),mat_info(k)%xs_p(j)%mv_f(i),&
          mat_info(k)%xs_p(j)%sig_t(i),(0.0, sgrp=1,fu_front),&
          (mat_info(k)%xs_p(j)%sig_gp2g(sgrp,i)*(2*j-1),sgrp= num_grp,1,-1), &
          (0.0, sgrp=1,fu_end) 

   enddo
   
   else
     do i=1, num_grp

	  write(xsunit, ioform,err=1003) &
      mat_info(k)%xs_p(j)%sig_a(i),mat_info(k)%xs_p(j)%mv_f(i),&
          mat_info(k)%xs_p(j)%sig_t(i),&
          (mat_info(k)%xs_p(j)%sig_gp2g(sgrp,i)*(2*j-1),sgrp= i,1,-1), (0.0, sgrp=1, num_grp-i)
	 enddo
   endif
 enddo
 enddo

close (xsunit)
write(LOGUNIT,*) '---------------------xs output files-------------------'
write(*,"('full scattering matrix file: ', A,2x, 'iht=', I0, 2x, 'ihs=', I0, 2x, 'ihm=', I0, ' PnOrder=', I0)") &
    trim(curfile), fu_iht, fu_ihs, fu_ihm, order_uni
write(LOGUNIT,"('full scattering matrix file: ', A,2x, 'iht=', I0, 2x, 'ihs=', I0, 2x, 'ihm=', I0, ' PnOrder=', I0)") &
    trim(curfile), fu_iht, fu_ihs, fu_ihm, order_uni
if(order_uni .gt. 0) then
write(*,"('2L+1 factor is pre-multiplied, PENTRAN xs type=1')")
write(LOGUNIT,"('2L+1 factor is pre-multiplied, PENTRAN xs type=1')") 
endif
!chi file
curfile=trim(prbname)//'.chi'
open (unit=xsunit, file=trim(curfile) )
 
do mp=1, num_mat_r
k=n_map(mp)
write (xsunit, "('/# fission chi ', A)") mat_info(k)%name
write(xsunit, "(10(ES12.5,2x) )") mat_info(k)%chi
enddo !k
close (xsunit)
write(*,"('chi file:', A)") trim(curfile)
write(LOGUNIT,"('chi file:', A)") trim(curfile)
!nubar file
curfile=trim(prbname)//'.nub'
open (unit=xsunit, file=trim(curfile) )
 
do mp=1, num_mat_r
k=n_map(mp)
write (xsunit, "('/# fission nubar ', A)") mat_info(k)%name
write(xsunit, "(10(ES12.5,2x) )") mat_info(k)%nubar
enddo !k
close (xsunit)
write(*,"('nubar file:', A)") trim(curfile)
write(LOGUNIT,"('nubar file:', A)") trim(curfile)
!grp file
curfile=trim(prbname)//'.grp'
open (unit=xsunit, file=trim(curfile) )
 

write (xsunit, "( '# Upper Bins in MeV, problem name: ', A )") trim(prbname)
write (xsunit, "( ' num_grps   FissNeutTemp(t,e,f)' )") 
if(num_grp_p .eq. 0 .or. num_grp_n .eq. 0) then 
  write (xsunit, "( ' ' , I8, '      t  ' )") num_grp
else
   write (xsunit, "( ' ' , I8, '      t  ', ' with ', I0, ' photon groups')") num_grp_n, num_grp_p
endif
if(num_grp_n*num_grp_p .eq. 0) then
 write (xsunit, "( ' Grp#     Upper_MeV' )")
 do k=1, num_grp_n
  write(xsunit, "(I6, 2x, ES12.5 )") k, upper_bnd_n(k)
 enddo !k

 do k=1, num_grp_p
  write(xsunit, "(I6, 2x, ES12.5 )") k, upper_bnd_p(k) 
 enddo !k
else 
 write (xsunit, "( '/ Grp#     Upper_MeV  Particle_type' )")
 do k=1, num_grp_n
  write(xsunit, "(I6, 2x, ES12.5,'      neutron' )") k, upper_bnd_n(k)
 enddo !k

 do k=1, num_grp_p
  write(xsunit, "(I6, 2x, ES12.5, '      photon  ',  I4 )") k+num_grp_n, upper_bnd_p(k),k 
 enddo !k

endif


close (xsunit)
write(*,"('group structure file:', A)") trim(curfile)
write(LOGUNIT,"('group structure file:', A)") trim(curfile)
!xrf file
curfile=trim(prbname)//'.xrf'
open (unit=xsunit, file=trim(curfile) )
 

write (xsunit, "( '/# File created by xsmcnp for GMIX: ', A )") trim(prbname)
write (xsunit, "( 'entry     nuclide     P0-id        mass          MI   Comment' )") 

do mp=1, num_mat_r
k=n_map(mp)
m=mat_info(k)%nxs(2)  !zzaaa number
mi=0
do i=1,mp-1
j=n_map(i)
if(mat_info(j)%nxs(2) .eq. m ) then
mi=mi+1
endif
enddo
write (xsunit, "(I4,6x, I2.2,I3.3,7x,I4, 8x, f12.7,3x, I2, 3x)") mp, &
      mat_info(k)%z_num,mat_info(k)%a_num, (mp-1)*(leg_order+1)+1, mat_info(k)%mass, mi
enddo !mp

close (xsunit)
write(*,"('material list file .xrf:', A)") trim(curfile)
write(LOGUNIT,"('material list file .xrf:', A)") trim(curfile)
return

1003 write(*,"('Output xs file error')") 
stop "output xs file error"
end subroutine                

subroutine TimeStamp(unit_num)
  use mParaset

  integer unit_num

  character(LEN =8) str_date
  character(LEN =12) str_time

  call date_and_time (str_date, str_time)

  if(unit_num .eq. 0) then

  write(*, "(A, 4x,A2,':',A2,':',A2,2x,'on ',A2,'/',A2,'/', A4)" )  trim(form_str), &
     str_time(1:2),str_time(3:4),str_time(5:6),str_date(5:6),str_date(7:8), str_date(1:4)
  write(*,*)

  else
  write(unit_num, "(A, 4x,A2,':',A2,':',A2,2x,'on ',A2,'/',A2,'/', A4)" )  trim(form_str), &
     str_time(1:2),str_time(3:4),str_time(5:6),str_date(5:6),str_date(7:8), str_date(1:4)
  write(unit_num,*)
  endif
end subroutine
