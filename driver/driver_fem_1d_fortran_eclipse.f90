! ============================================================================
! name        : driver_fem_1d_fortran_eclipse.f90
! author      : amir
! version     :
! copyright   : your copyright notice
! description : hello world in fortran
! ============================================================================
program driver_fem_1d_fortran_eclipse
use linearsolvers
use fem_geometry
use material_behavior
use fem_libs
    implicit none
!     ================================================================
!                           the variables
!     ================================================================
      integer :: icalib !calibration parameter
!     ================================================================
!                      solution variables
!     ================================================================
      real*8,allocatable::glk(:,:)
      real*8,allocatable::glf(:) !the global source vector
      real*8,allocatable::glu(:) !the global solution
      real*8,allocatable::glp(:) !the previus degrees of freedom vector
      real*8,allocatable::glb(:) !the before previus degrees of freedom

      real*8,allocatable::glu_ref(:) !the global solution from refrence configuration
      real*8,allocatable::glp_ref(:) !the previus degrees of freedom vector
      real*8,allocatable::glb_ref(:) !the before previus degrees of freedom
!     ================================================================
!                      solver variables
!     ================================================================
      real*8,allocatable::gls(:) !the sulution vector out of linear solver
      real*8,allocatable::glq(:) !internal force vector
      real*8,allocatable::glt(:) !external fource vector
      real*8,allocatable::glr(:) !total residual vector
!     ============================iteration variables
      real*8::  eps,error,normal !the error tolerance and error
      integer:: itmax !maxumim number of iteration
      integer:: iter ! iteration counter
      logical:: successful,converge
!     ================================================================
!                      file unit
!     ================================================================
      integer::curved_node !curved node for
!     ================================================================
!                      time variables
!     ================================================================
!      integer::ntime
!      real*8::deltime
      real*8::loadfactor    ! time
      real*8::freq
      real*4::timearray, telapse
      integer ( kind = 4 )start_stamp_v(8),end_stamp_v(8),values(8)
      character*8::ctime(2)
      integer :: clck_counts_beg, clck_counts_end, clck_rate
      real ::  beg_cpu_time, end_cpu_time
!     ================================================================
      real*8::k00,k01,labda01
      real*8::labda_e_01,k_e_00,k_e_01
      real*8::labda_k_01,k_k_00,k_k_01
!     ================================================================
      integer::inode,jnode,ibond,jbond ,i,j,id,in
!     ================================================================
      integer::igauss,ngauss
!     =======================trivial meshing arrays the meshing seed parameters
      integer ::neldirectional(dimen)
      real*8  ::length(dimen)
!     ================================================================
      common/iteration/iter
      common/visco_elastic/k00,k01,labda01,            &
                          k_e_00,k_e_01,labda_e_01,     &
                          k_k_00,k_k_01,labda_k_01
!     ===============================================================
!                       p r e p r o c e s s o r   u n i t
!     ===============================================================
      call system_clock ( clck_counts_beg, clck_rate)
      call cpu_time (beg_cpu_time)
      call timestamp()
!     ===============================================================
!     reading iteration and convergence critria's
!     ===============================================================
      itmax=30;eps=1e-2
!     ======================== ve formulation
      k00=1.0 ; k01=0.0; labda01=0.0
      k_e_00=0.4 ; k_e_01=0.25; labda_e_01=0.5 ;
      k_k_00=1.0 ; k_k_01=0.0; labda_k_01=0.0 ;
!     ===============================================================
!     reading time increment varibales
!     ===============================================================
      dtime=1.0;freq=0.1
      ntime=100

    length=10.0
    neldirectional(1)=100
      call geometry_1D(neldirectional,length)
!      call printMatrix(coords,'c')
!      call printMatrix_int(nod,'nod')
!     ===============================================================
!     define the solution parameters
!     ===============================================================
      nnm=size(coords,dim=1)
      nem=size(nod,dim=1)
      neq  = nnm*ndf;     nn=npe*ndf;
      ipdf = iel+1;    ngauss=(ipdf**dimen)*nem
!     ===============================================================
!     reading the boundary conditions
!     ===============================================================
      call bounday_1D()
      allocate(glk(neq,neq),glu(neq),glq(neq),glt(neq),glr(neq))
      glu=0.0d0;glt=0.0d0;glr=0.0d0;
      call  form_history(ngauss,dimen,ipdf,nem)
!     ===============================================================
!                        time increment starts here
!     ===============================================================
      vspvt=vspv
      vssvt=vssv

      do itime=0,ntime
      time(1)=dtime;time(2)=itime*dtime

     loadfactor=1.0d0*time(2)

!      loadfactor=250*sawtooth(time(2)*freq)*1.0

!     loadfactor=250*sin(2.0*pi*time(2)*freq)*1.0e3

      vspv=loadfactor*vspvt
      vssv=loadfactor*vssvt

      glt=0.0d0;
      glu(bnd_no_pr_vec)=0.0d0;
!     ===============================================================
!                 nonlinear solution iteration starts here
!     ===============================================================
      iter=0;CONVERGE=.false.;   ;
130   iter=iter+1;
!      write(1,*)'iter',iter
      if(iter.gt.itmax)then; write(*,330);goto 500;endif
!     ===============================================================
!                        forming the global matrices
!     ===============================================================
      glk=0.0d0;glq=0.0d0;!
!      write(1,*)'glu=',(glu(i),i=1,neq);
      call glbmatrcs(glu,glk,glq)
!      write(1,*)'glq=',(glq(i),i=1,neq);
!     ===============================================================
!                             newton raphson sprocedure
!     ===============================================================
      glt(bnd_no_se_vec) = vssv ;
      glr=glt-glq ;

!      write(1,*)'glt=',(glt(i),i=1,neq);
!      write(1,*)'glr=',(glr(i),i=1,neq);

!     ===============================================================
!                        solving the linear system of equations
!     ===============================================================
      call symmetric_primary_bounday(glk,glr)
      vspv=0.0d0
!!      write(1,*)'glr=',(glr(i),i=1,neq);

      call gaussian_elimination_solver(glk,glr);  gls=glr
!!      write(1,*)'gls=',(gls(i),i=1,neq);
!     call lapack_gesv( glk, glr ); gls=glr
!     ===============================================================
!     call result_printer_1D(iter,glr)
      glu=glu+gls
!
      normal=norm_vect(glu);if(normal.eq.0)then;normal=1;endif
      error=norm_vect(gls)
      error=error/normal

      write(out,*)'time',time
      write(out,*)'error',error
      write(out,*)'iter',iter
!

      if(error.gt.eps)then;converge=.false. ;goto 130;endif;
      converge=.true.
!     ===============================================================
!                        time increment ends here
!     ===============================================================
!      glu_ref=glu+glu_ref
!      glb_ref=glp_ref
!      glp_ref=glu_ref
!
!      glb=glp
!      glp=glu

!      call update_history(ngauss,dimen)
!     ======================updatng coordinates  updated lagrangian
!      forall(i=1:nnm) coords(1:dimen,i)=
!     1 coords(1:dimen,i)+glu((i-1)*ndf+[1:dimen])
!     ===============================================================
!     plotting the converged result
!     ===============================================================
    call result_printer_1D(iter,glu)


    enddo !itime=0,ntime


500   call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()

999 FORMAT (11(:,1X,F8.5) )
330   FORMAT(/,5X,'***** CONVERGENCE CRITERION IS NOT SATISFIED *****')
end program

