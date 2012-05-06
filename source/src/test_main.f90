! -------------------------------------------------------------------------------
program test_main
! -------------------------------------------------------------------------------
  use mod_precision
  use mod_leconf
  use mod_UnitUtil
  use mod_Iion
  !.
  implicit none

  type(t_cnf)               :: cnf
  character(len=256)        :: outname
  integer(ip), parameter    :: mxprm=100
  integer(ip)               :: lu_curr,lu_stat,iteration
  integer(ip)               :: aux, nprm, ict, outfreq, init_save
  real(rp)                  :: step(3), Istm, U, time, Qion
  real(rp)                  :: t1,t2,v_prm(mxprm)
  real(rp), allocatable     :: v_ve(:)
  real(rp), allocatable     :: v_cr(:)
  real(rp), allocatable     :: dat_r(:)
  integer(ip), allocatable  :: dat_i(:)
  integer                   :: numeroDeParametro
  integer                   :: numeroDeCurrent
  integer(ip)               :: valorDeParametro

  !.

  write(*,*) 'Begining reading ...'
  call readconfile(cnf)
  write(*,*) 'Done reading ...'
  ict = get_celltype(cnf)           ! Retrieve cell type
  call PrintNumCellEM ('Solving ',ict,ict)
  nprm = get_number_of_prm(cnf%prm) ! Retrieve mask
  write(*,*) 'Number of param to modified : ',nprm
  if (nprm>0) then
      allocate(dat_r(nprm))
      allocate(dat_i(nprm))
      call get_mask(cnf%prm, dat_i, dat_r)
      write(*,*) 'Mask:'
      write(*,10) 'dat_i: ', (dat_i(aux),aux=1,nprm)
      write(*,20) 'dat_r: ', (dat_r(aux),aux=1,nprm)
  endif
  call get_step(cnf,step)   ! Retrieve simulation time and increment
  write(*,*) 'Step'
  write(*,30) 'Tini: ', step(1)
  write(*,30) 'DT  : ', step(2)
  write(*,30) 'Tfin: ', step(3)
  outfreq =  get_post_freq(cnf) 
  init_save= get_init_save(cnf) 
  write(*,40) 'Saving starts at increment: ',init_save
  write(*,40) 'OutputFreq: ',outfreq

  call fun_getarg(3,'-o',outname) !.Retrieve root name for outputfile
  aux=len_trim(outname)
  call a_unit(lu_curr) 
  call a_unit(lu_stat) 
  open(lu_curr,FILE=outname(1:aux)//'_curr.dat');
  open(lu_stat,FILE=outname(1:aux)//'_stat.dat');

  iteration =0;
  time = step(1);
  ! Initialicing state variables and parameters
  aux = get_numstatevar(ict)
write(*,*) "get_numstatevar=>",aux
  allocate(v_ve(aux))
  v_ve(1:aux)=get_ic_em(ict,aux)
write(*,*) "get_ic_em=>",aux
  aux = get_numcur(ict)
write(*,*) "get_numcur=>",aux
  allocate(v_cr(aux))
  call get_parameter_nd(ict,v_prm,aux);
write(*,*) "get_parameter_nd=>",aux
  if(nprm>0) then
      v_prm(dat_i(1:nprm))=dat_r(1:nprm)
  endif
  U = initcond_elecmodel (ict) 
  if( allocated(cnf%file_parameters_i)) then
          write(*,*) "Numero de Parametros a Disco=>",size(cnf%file_parameters_i)
          do numeroDeParametro=1 , size(cnf%file_parameters_i)
          write(*,*) "Parametro=>", cnf%file_parameters_i(numeroDeParametro:numeroDeParametro)
          enddo
  endif
  if( allocated(cnf%file_currents_i)) then
          write(*,*) "Numero de Corrientes=>",size(cnf%file_currents_i)
          do numeroDeCurrent=1 , size(cnf%file_currents_i)
          write(*,*) "Parametro=>", cnf%file_currents_i(numeroDeCurrent:numeroDeCurrent)
          enddo
  endif


  call get_time(t1)
  do while (time < step(3))
    Istm = get_stimulus(cnf%stm, time, step(3))
!
    if(iteration.ge.init_save) then
            if(mod(iteration,outfreq) == 0) then
                    write(*,101) iteration
                    write(*,*) "time=>" , time
                    write(*,*) "U=> " , U
                    write(*,*) "v_ve(1:1)=> ", v_ve(1:1)
                    ! write to file only the selected params or all if any.
                    ! PARAMETERS
                    if(allocated(cnf%file_parameters_i))then
                            if(size(cnf%file_parameters_i) == 0) then
                                    !no guardo nada
                                    !write(lu_stat,50) time, U, v_ve
                            else
                                    do numeroDeParametro=1 , size(cnf%file_parameters_i)
                                    valorDeParametro = cnf%file_parameters_i(numeroDeParametro)
                                    if(numeroDeParametro==1) then
                                            if(size(cnf%file_parameters_i) == 1) then
                                                    write(lu_stat,50) time, U , v_ve(valorDeParametro:valorDeParametro)
                                            else
                                                    write(lu_stat,90) time, U , v_ve(valorDeParametro:valorDeParametro)
                                            end if
                                    else
                                            write(lu_stat,95) v_ve(valorDeParametro:valorDeParametro)
                                    end if
                                    enddo ! loop for every param.
                                    write(lu_stat,'()')   ! SALTO LINEA
                            end if ! params or not.
                    endif
                    ! CURRENTS
                    if(allocated(cnf%file_currents_i))then
                            if(size(cnf%file_currents_i) == 0) then
                                    !no guardo nada
                                    !write(lu_curr,50) time, v_cr
                            else
                                    if(size(cnf%file_currents_i) == 1 .AND. cnf%file_currents_i(1)==0) then
                                            !guardo todo
                                            write(lu_curr,50) time, v_cr
                                    else
                                            !guardo solo lo indicado.
                                            do numeroDeParametro=1 , size(cnf%file_currents_i)
                                            valorDeParametro = cnf%file_currents_i(numeroDeParametro)
                                            if(numeroDeParametro==1) then
                                                    if(size(cnf%file_currents_i) == 1) then
                                                            write(lu_curr,50) time,     v_cr(valorDeParametro:valorDeParametro)
                                                    else
                                                            write(lu_curr,90) time,     v_cr(valorDeParametro:valorDeParametro)
                                                    end if
                                            else
                                                    write(lu_curr,95) v_cr(valorDeParametro:valorDeParametro)
                                            end if
                                            enddo ! loop for every param.
                                            write(lu_curr,'()')   ! SALTO LINEA
                                    end if ! params or not.
                            endif
                    endif
            end if
    endif
!
    call get_Iion (ict,step(2),Istm,0.0_rp,U,0.0_rp,Qion, &
                   v_prm(1:mxprm),v_ve, v_cr)
    U = U - step(2) * (Qion+Istm);
    time = time + step(2);
    iteration = iteration + 1;
!
  enddo
  call get_time(t2)
  write(*,60) 'Elapsed time:',t2-t1,'s'
  close(lu_curr)
  close(lu_stat)

10 format (A,40(2X,I3))
20 format (A,40(2X,F12.6))
30 format (A,2X,F16.2)
40 format (A,2X,I15)
50 format (F15.2,50(2X,F15.10))
60 format (A,2X,F16.3,2X,A)
90 format (F15.2,50(2X,F15.10,$))
!95 format (50(2X,F15.10))
95 format (50(2X,F15.10,$))
!60 format (F14.6,F12.6)
101 format (10(I10,$))
end program test_main


