! ------------------------------------------------------------------------------
module mod_leconf
! ------------------------------------------------------------------------------
  use mod_precision
  use mod_UnitUtil
  use mod_readcir
  !.
  type t_msk
    private
    integer(ip), allocatable :: dat_i(:)
    real(rp),    allocatable :: dat_r(:)
  end type t_msk
  !.
  type t_stm
    private
    integer(ip)          :: ntrain
    real(rp), allocatable    :: data(:)
  end type t_stm
  !.
  type, public:: t_cnf
    integer(ip)           :: ict
    integer(ip)           :: post(2)
    real(rp)              :: step(3)
    type(t_msk)           :: prm
    type(t_stm)           :: stm
    integer(ip), allocatable :: file_parameters_i(:)
    integer(ip), allocatable :: file_currents_i(:)
  end type t_cnf
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine readconfile (cnf)
! ------------------------------------------------------------------------------
  implicit none
  type(t_cnf)           :: cnf
  !.
  
  integer(ip), parameter   :: MAX=100
  integer(ip)              :: lu_cfg, n, m, io_err
  integer(ip)              :: dat_i(MAX)
  real(rp)                 :: dat_r(MAX)
  character (len = 256)    :: fname,strout
  logical                  :: flag
  integer(ip)              :: file_parameters_i(MAX)
  integer(ip)              :: file_currents_i(MAX)

  call a_unit(lu_cfg)
  call fun_getarg (1,'-i',fname)  !.Nombre del archivo de configuracion
  m=len_trim(fname)
  open(unit=lu_cfg, file=fname(1:m), status='old')
  !.
  do 
    call srch_lec (lu_cfg,'#',strout, flag, io_err)
    !.
    if (flag) then
      if     (index(strout,'#MODEL')  == 1) then
        !.
        read (lu_cfg, * ,iostat=io_err) cnf%ict
        if (io_err > 0) stop '>>>.Error reading the cell type'
        !.
      elseif (index(strout,'#PARAMETERS') == 1) then
        !.
        read (lu_cfg, * ,iostat=io_err) n,(dat_i(m),m=1,n),(dat_r(m),m=1,n)
        if (io_err > 0) stop '>>>.Error reading parameter mask and values'
        allocate (cnf%prm%dat_i(n))
        allocate (cnf%prm%dat_r(n))
        cnf%prm%dat_i(1:n)=dat_i(1:n)
        cnf%prm%dat_r(1:n)=dat_r(1:n)
        !.
      elseif (index(strout,'#STEP') == 1) then
        !.
        read (lu_cfg, * ,iostat=io_err) (cnf%step(m),m=1,3)
        if (io_err > 0) stop '>>>.Error reading step information'
        !.
      elseif (index(strout,'#STIMULUS')    == 1) then
        !.
        read (lu_cfg,*,iostat=io_err) n
        cnf%stm%ntrain = n
        if (io_err > 0) stop '>>>.Error reading the number of stimulus'
        allocate (cnf%stm%data(4*n))
        read (lu_cfg, * ,iostat=io_err) (cnf%stm%data(m),m=1,4*n)
        if (io_err > 0) stop '>>>.Error reading stimulus'
      elseif (index(strout,'#POST')     == 1) then 
        !.
        read (lu_cfg, * ,iostat=io_err) (cnf%post(m),m=1,2)
        if (io_err > 0) stop '>>>.Error reading postprocess frequency'
        !.
      elseif (index(strout,'#FILE_PARAMETERS') == 1) then
        !.
        read (lu_cfg, * ,iostat=io_err) n,(file_parameters_i(m),m=1,n)
        if (io_err > 0) stop '>>>.Error reading file parameters mask and values'
        allocate (cnf%file_parameters_i(n))
        cnf%file_parameters_i(1:n)=file_parameters_i(1:n)
        !cnf%file_parameters_i(1:n)=dat_i
        !.
       elseif (index(strout,'#FILE_CURRENTS') == 1) then
        !.
        read (lu_cfg, * ,iostat=io_err) n,(file_currents_i(m),m=1,n)
        if (io_err > 0) stop '>>>.Error reading file currents mask and values'
        allocate (cnf%file_currents_i(n))
        cnf%file_currents_i(1:n)=file_currents_i(1:n)
        !.
      else
          if (index(strout,'#END')     == 1) exit
      endif
    else
      if (io_err == -1) exit
    endif
  end do
  !.
  call close_unit (lu_cfg)
  !.
  return
end subroutine readconfile
! ------------------------------------------------------------------------------
function get_celltype (cnf) result(ict)
! ------------------------------------------------------------------------------
 type(t_cnf)            :: cnf     
 integer(ip)            :: ict

 ict = cnf%ict
 return
end function get_celltype
! ------------------------------------------------------------------------------
function get_number_of_prm(prm) result(n)
! ------------------------------------------------------------------------------
  implicit none
  type(t_msk), intent(in) :: prm
  integer(ip)             :: n

  n = size(prm%dat_i)
  if (n<0) then
    n=0
  endif
end function get_number_of_prm
! ------------------------------------------------------------------------------
subroutine get_mask(prm, dat_i, dat_r)
! ------------------------------------------------------------------------------
  implicit none
  type(t_msk), intent(in) :: prm
  integer(ip), intent(out):: dat_i(:)
  real(rp),    intent(out):: dat_r(:)
  integer(ip)             :: n

  n = size(prm%dat_i)
  if (n>0) then
      dat_i(1:n) = prm%dat_i
      dat_r(1:n) = prm%dat_r
  endif
end subroutine get_mask
! ------------------------------------------------------------------------------
subroutine get_step(cnf, vec)
! ------------------------------------------------------------------------------
  implicit none
  type(t_cnf),intent(in):: cnf
  real(rp),intent(out)  :: vec(:)

  vec(:) = cnf%step
  
  return
end subroutine get_step
! ------------------------------------------------------------------------------
function get_init_save(cnf) result(init_save)
! ------------------------------------------------------------------------------
  implicit none
  type(t_cnf),intent(in) :: cnf
  integer(ip)            :: init_save

  init_save = cnf%post(1)
  
  return
end function get_init_save
! ------------------------------------------------------------------------------
function get_post_freq(cnf) result(outfreq)
! ------------------------------------------------------------------------------
  implicit none
  type(t_cnf),intent(in) :: cnf
  integer(ip)            :: outfreq

  outfreq = cnf%post(2)
  
  return
end function get_post_freq
! ------------------------------------------------------------------------------
function get_stimulus(stm, t, tend) result  (Istm)
! ------------------------------------------------------------------------------
  implicit none
  type(t_stm),intent(in) :: stm
  real(rp),intent(in)    :: t, tend 
  real(rp)               :: Istm
  integer(ip)            :: j,ntrain
  real(rp)               :: tist, tfst

  Istm = 0.0_rp
  ntrain = stm%ntrain
  j=0
  tist = stm%data(1)
  do while (tist .le. tend)
    tfst = tist+stm%data(4*j+3)
    if((t>=tist).and.(t<=tfst)) then
         Istm = -stm%data(4*j+4)
         exit
    endif
    tist=tist+stm%data(4*j+2)
    if(j<(ntrain-1)) then
        if (tist >= stm%data(4*(j+1)+1)) then
            j=j+1
            tist=stm%data(4*j+1)
        endif
    endif
  enddo
  return
end function get_stimulus
! ------------------------------------------------------------------------------
end module mod_leconf
