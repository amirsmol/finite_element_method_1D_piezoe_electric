      module material_behavior

      use linearsolvers

!     =========================element center point
      real*8 :: el_center_coord(dimen)
      real*8 :: qtransform(dimen,dimen) ! this is the transformtion matrix
!     ================================================================
!                      history variables
!     ================================================================
!      real*8,allocatable::gauss_coords(:,:) !gauss_coords(igauss,dimen)
      real*8,allocatable::ele_pol_vector(:,:)  !gauss_coords(ele,dimen)
!      real*8,allocatable::shape_function_weight(:,:) ! shape_function_weight(nem,igauss_element)
!     ================================================================
!                      history variables
!     ================================================================
      real*8,allocatable::q_t(:),q_t_p(:),q_t_b(:)
      real*8,allocatable::                      &
        q_sigma_c_t(:,:,:),q_sigma_e_t(:,:,:),  &
        q_sigma_c_p(:,:,:),q_sigma_e_p(:,:,:),   &
        q_sigma_c_b(:,:,:),q_sigma_e_b(:,:,:),    &
        q_d_e_t(:,:),q_d_k_t(:,:),                 &
        q_d_e_p(:,:),q_d_k_p(:,:),                  &
        q_d_e_b(:,:),q_d_k_b(:,:)
!     ================================================================
!                      polarization history variables
!     ================================================================
      real*8,allocatable::      &
        q_rempol_t(:,:),         &
        q_rempol_p(:,:),          &
        q_rempol_b(:,:)

!     =======================mechanical properties variables
      real*8::cmate_voight(9,9)
      real*8::ctens(3,3,3,3)
      real*8::ctenst(3,3,3,3)
      real*8::cmat_voight(6,6)
      real*8::stense(3,3,3,3)
      real*8::e1,e2,e3,nu12,nu13,nu23,g12,g13,g23
      real*8::lambda,nu,mu,kbulk
!     =======================electric variables
      real*8,parameter::eps_vac_permit=8.58e-12
      real*8::epz(3,3,3)
      real*8,dimension(3,3)::ktense,ktense_trans
      real*8::ktense_inv(3,3) !the material matrix
      real*8::xeta(3,3,3)! the non linear electric permibility tensor

!     ===========================materials tensors
      real*8,dimension(3,3,3):: epz_t,epz_trans,epz_nl !the electro coupling tensor
      real*8,dimension(3,3,3,3)      :: b_tilt,b_tilt_trans ! the nonlinear electro coupling tensor

!     ===================================remanent polarization
      real*8,dimension(dimen):: rem_pol_t,drem_pol_t
      real*8,dimension(dimen):: rem_pol_p,drem_pol_p
      real*8,dimension(dimen):: rem_pol_b
!     ================================================================
!                      polarization history variables local
!     ================================================================
      real*8,dimension(dimen)::  &
        q_rempol_t_local,        &
        q_rempol_p_local,        &
        q_rempol_b_local
!     =======================polarization switching variables
      real*8::e_coercive ! coercive field
      integer::ipol
      real*8::p_per_sat
      real*8::pol_rate
      real*8::der_total(3)
!     ================================================================
!                      time variables
!     ================================================================
      real*8::time(2),dtime
      integer::ntime,itime  ! time
!     ================================================================
      integer::noelem
!     ================================================================
!     ================================================================
!     ================================================================
!     ================================================================
!     ================================================================
      contains

      subroutine material_properties()
      implicit none
!     ================================================================
!                      material variables
!     ================================================================

!     ================================================================
      integer::eye(dimen,dimen)
      integer::i,j,k,l,m,n

!     =======================polarization switching variables
      real*8 ::ecoup_norma
      real*8 ::ecoup_perpe
      real*8 ::ecoup_shear
!     =======================electric polarization (displacement in refrence confiquaration)
      real*8 :: polarization_vec(dimen)
      real*8 :: der_total(dimen),pol(3)
      real*8 :: der_unit(3)
      real*8,dimension(dimen)::  elr_t,elr_p,elr_b
      real*8 :: alpha_de(3,3)
      real*8 :: pol_intesity ,sat_pol
      real*8 :: kappa_0
!     ================================================================
      real*8::ey
      real*8::k00,k01,labda01
      real*8::labda_e_01,k_e_00,k_e_01
      real*8::labda_k_01,k_k_00,k_k_01
      common/visco_elastic/k00,k01,labda01,            &
                          k_e_00,k_e_01,labda_e_01,     &
                          k_k_00,k_k_01,labda_k_01
!     ================================================================
      eye = 0.0; do i = 1,dimen; eye(i,i) = 1;enddo

!      e_coercive=0.67e6  ![v/m] coercive field
!      pol_rate=5.4   !remnant polarization rate
!      p_per_sat= 0.26 ! [c/m^2]
!      kappa_0=0.04e-6
!!
!      e1= 6.06e10; ; ![pa]
!      nu23=0.3;
!     ==============================general polarization
      pol_intesity=norm_vect(polarization_vec);

      e1= 1 ![n/m^2] pzt 5a;
!      e1= 62e9 ![n/m^2] pzt 5a;
      nu23=0.3;

      ktense=0.0d0
      ktense(1,1)=(4.0d0)*1e-9;
      ktense(2,2)=ktense(1,1)
      ktense(3,3)=ktense(1,1)/2.0;

      ecoup_norma =-44.69;   ![n/(v.m)]
      ecoup_perpe =18.93;    ![n/(v.m)]
      ecoup_shear =-23.27;   ![n/(v.m)]

!      k00=1.0 ; k01=0.0; labda01=0.0
!      k_e_00=1.0 ; k_e_01=1.0; labda_e_01=1.0 ;
!      k_k_00=1.0 ; k_k_01=0.0; labda_k_01=0.0 ;

      e2=e1;e3=e2
      nu13=nu23;nu12=nu13;
      g23 = e1/(1+nu12)/2.0;
      g13 = g23;g12=g23;
      ey=e1
      nu=nu12

      mu=ey/(1+nu)/2.0
      lambda=ey*nu/(1+nu)/(1-2.0*nu)
      kbulk=lambda+2.0*mu/3.0

      der_unit=unit_vect(polarization_vec)
!     ==============================test
       der_unit=0.0d0
       der_unit(1)=1.0d0
!     ==============================
      do i=1,dimen
      do j=1,dimen
      alpha_de(i,j)=eye(i,j)-der_unit(i)*der_unit(j)
      enddo;enddo

      epz=0.0d0
      do i=1,dimen
      do j=1,dimen
      do k=1,dimen
      epz(k,i,j)=ecoup_norma*der_unit(k)*der_unit(i)*der_unit(j)  &
               +ecoup_perpe*der_unit(k)*alpha_de(i,j)             &
               +ecoup_shear*der_unit(i)*alpha_de(j,k)*0.5         &
               +ecoup_shear*der_unit(j)*alpha_de(i,k)*0.5
      enddo;enddo;enddo


     ktense_inv=0.0d0
     forall(k=1:3) ktense_inv(k,k)=1/ktense(k,k)

      xeta=0.0d0;
      b_tilt=0.0d0

      do i=1,dimen;do j=1,dimen;do k=1,dimen;do l=1,dimen
      ctens(i,j,k,l)=lambda*eye(i,j)*eye(k,l)+mu*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))
      enddo;enddo;enddo;enddo

      do i=1,dimen;do j=1,dimen;do k=1,dimen;do l=1,dimen
      stense(i,j,k,l)=(1.0/kbulk/9.0)*eye(i,j)*eye(k,l)+(1/mu/4.0)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)-(2.0/3.0)*eye(i,j)*eye(k,l))
      enddo;enddo;enddo;enddo
!     ==============================x3 polarized
!     epz=0.0d0
!     epz(3,3,3)=45.744;
!     epz(3,1,1)=-19.375;
!     epz(1,1,3)=23.81
!
!     epz(3,2,2)=epz(3,1,1)
!
!     epz(2,2,3)=epz(1,1,3)
!     epz(1,3,1)=epz(1,1,3)
!     epz(2,3,2)=epz(1,1,3)
!     b_tilt(3,3,3,3)=-3.35e-2;
!     b_tilt(3,3,1,1)=3.35e-2;
!     b_tilt(3,3,2,2)=3.35e-2;
!     ==============================x3 polarized
      ktense=0.0d0
      ktense(1,1)=(4.0d0)*1e-9;
      ktense(2,2)=ktense(1,1)
      ktense(3,3)=ktense(1,1)/2.0;


!      call material_transofrmation_tensor_x2()
!      call tensor_transform_rank_3(qtransform,epz,epz_trans)
!      call tensor_transform_rank_4(qtransform,b_tilt,b_tilt_trans)
!      call tensor_transform_rank_2(qtransform,ktense,ktense_trans)
!
!      epz=epz_trans
!      b_tilt=b_tilt_trans
!      ktense=ktense_trans

      epz_t=0.0;
      epz_nl=0.0;

      do 50 i=1,dimen; do 50 j=1,dimen;do 50 k=1,dimen;do 50 l=1,dimen
      epz_nl(k,i,j)=epz(k,i,j)+0.5*b_tilt(k,l,i,j)*abs(elr_t(l))
50    epz_t(k,i,j)=epz(k,i,j)+b_tilt(k,l,i,j)*abs(elr_t(l))


!    1D formulation
ctens=0; ctens(1,1,1,1)=1;
epz=0; epz(1,1,1)=-1;
ctenst=ctens;
epz_t=epz;
ktense=0.0d0  ;ktense(1,1)=1;

      endsubroutine material_properties


      subroutine  update_history(ngauss,dimen)
      integer::ngauss,dimen

      q_sigma_c_b=q_sigma_c_p;
      q_sigma_c_p=q_sigma_c_t;

      q_sigma_e_b=q_sigma_e_p;
      q_sigma_e_p=q_sigma_e_t;

      q_d_e_b=q_d_e_p
      q_d_e_p=q_d_e_t

      q_d_k_b=q_d_k_p
      q_d_k_p=q_d_k_t

      q_rempol_b=q_rempol_p
      q_rempol_p=q_rempol_t

      endsubroutine  update_history

!     ===============================================================
!           forming history variables
!     ===============================================================
      subroutine  form_history(ngauss,dimen,ipdf,nem)
      integer::ngauss,dimen,ipdf,nem
!     ===================allocate history variables
      allocate(                                                                  &
      q_sigma_c_t(dimen,dimen,ngauss),q_sigma_e_t(dimen,dimen,ngauss),            &
      q_sigma_c_p(dimen,dimen,ngauss),q_sigma_e_p(dimen,dimen,ngauss),            &
      q_sigma_c_b(dimen,dimen,ngauss),q_sigma_e_b(dimen,dimen,ngauss),             &
        q_d_e_t(dimen,ngauss),q_d_k_t(dimen,ngauss),                                &
        q_d_e_p(dimen,ngauss),q_d_k_p(dimen,ngauss),                                 &
        q_d_e_b(dimen,ngauss),q_d_k_b(dimen,ngauss),                                  &
      q_t(ngauss),q_t_p(ngauss),q_t_b(ngauss),                                         &
      q_rempol_t(dimen,ngauss),                                                         &
      q_rempol_p(dimen,ngauss),                                                          &
      q_rempol_b(dimen,ngauss),          &
      ele_pol_vector(nem,dimen))



      q_sigma_c_t=0.0d0;q_sigma_e_t=0.0d0;
      q_sigma_c_p=0.0d0;q_sigma_e_p=0.0d0;
      q_sigma_c_b=0.0d0;q_sigma_e_b=0.0d0;

       q_d_e_t=0.0d0;q_d_k_t=0.0d0;
       q_d_e_p=0.0d0;q_d_k_p=0.0d0;
       q_d_e_b=0.0d0;q_d_k_b=0.0d0;

     q_rempol_t=0.0d0;
     q_rempol_p=0.0d0;
     q_rempol_b=0.0d0;
     ele_pol_vector=0.0d0

      q_t=0.0d0 ; q_t_p=0.0d0; q_t_b=0.0d0

      end subroutine  form_history

subroutine stress_elect_displacement(ur,der,sigma)
!     ================================================================
!                          input variables
!     ================================================================
real*8,intent(in) :: ur(:,:)! global coefficient
!     ================================================================
!                         output variables
!     ================================================================
integer :: i !integer counters
real*8  :: der(:),sigma(:,:)
!     ================================================================
real*8::strain(3,3);
real*8::electric_feild(3);
!     ================================================================
strain=0.0d0;electric_feild=0.0d0;
strain(1:dimen,1:dimen)=0.5d0*(   ur(1:dimen,:)+transpose(ur(1:dimen,:))+      &
                              matmul(  transpose(   ur(1:dimen,:)  ),ur(1:dimen,:))    );

electric_feild(1:dimen)=-ur(dimen+1,:)

!write(1,*)'strain=',strain(1,1)
!write(1,*)'electric_feild=',electric_feild(1)

sigma =  0.0d0 ;
der   =  0.0d0 ;

sigma(1,1)=ctens(1,1,1,1)*strain(1,1)-epz(1,1,1)*electric_feild(1)
der(1)=epz(1,1,1)*strain(1,1)+ktense(1,1)*electric_feild(1)


!write(1,*)'ktense(1,1)=',ktense(1,1)
!write(1,*)'sigma=',sigma(1,1)
!write(1,*)'der=',der(1)

endsubroutine stress_elect_displacement

      end module material_behavior
