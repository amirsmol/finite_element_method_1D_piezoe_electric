      module fem_libs
      use fem_geometry
      use material_behavior
      use linearsolvers
      implicit none
      real*8,allocatable:: elcrds(:,:) !the geometries of nodes of the element
      contains
!     ================================================================
!                 forming the global martixes and vetors
!     ================================================================
      subroutine glbmatrcs(glu,glk,glq)
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      real*8,intent(in)::glu(:) ! the dof vector, current time
!     ================================================================
!                          output variables
!     ================================================================
      real*8,intent(out)::glq(:) ! residual vector
      real*8,intent(out)::glk(:,:) ! the global stiffness matrix
!     ================================================================
!                          trivial variables
!     ================================================================
      real*8,allocatable:: elq(:),elk(:,:) !the force and coeffi
      real*8,allocatable:: elu(:) ! degrees of freedom vector in current time

!     ================================================================
!                      time variables
!     ================================================================
!      real*8::time(2),dtime
!      integer::ntime,itime  !time
!     ================================================================
!                      history variables
!     ================================================================
      integer::igauss,elneq
!     ================================================================
!      integer::i,j,ni,l,li,nn,n,neq,
      integer::inpe,iglobal,idof,jdof
      integer::inode,jnode
      integer::noelem
!     =========================element global node and dof map
      integer,allocatable::el_nod_vec(:) !el_nod_vec(npe)
      integer,allocatable::el_dof_vec(:)  !el_dof_vec(npe*ndf)
!     ================================================================

      elneq=npe*ndf;
      allocate(elk(elneq,elneq),elq(elneq),elu(elneq))
!     ================================================================
      elu=0.0d0;elq=0.0d0;elk=0.0d0;
      igauss=0


!     starting the element

      allocate(el_nod_vec(npe),el_dof_vec(elneq))

      do noelem=1,nem
      el_nod_vec=nod(noelem,:)
      elcrds=transpose(coords(el_nod_vec,:))
!      call printMatrix(elcrds,'elcrds');

      do inpe=1,npe
      iglobal=nod(noelem,inpe)
      idof=inpe*ndf
      do jdof=0,ndf-1
      el_dof_vec(idof-ndf+1+jdof)=(iglobal*ndf-ndf+1+jdof)

      enddo !j=1,ndf
      enddo !inpe=1,npe

!     call printVector_int(el_dof_vec,'el_dof_vec')
      elu=glu(el_dof_vec)
!      write(1,*)noelem;
!      write(1,999)(elu(idof),idof=1,elneq);
!     ===============================================================
!     finding the center point of element
!     ===============================================================

!     ===============================================================
!     forming the coefficitne matrix and source vector for each element
!     ===============================================================
      call elmatrcs_1d(elu,elk,elq)
!     ===============================================================
!     assembling the global source vector and coefficient matrix
!     ===============================================================
!      write(3,*)noelem;
!      write(3,*)(elq(iglobal),iglobal=1,elneq);
      glq(el_dof_vec)=glq(el_dof_vec)+elq
      glk(el_dof_vec,el_dof_vec)=glk(el_dof_vec,el_dof_vec)+elk


!     call printMatrix(elk,'elk');

     enddo ! noelem=1,nem
!    write(3,*);write(3,*)(glq(iglobal),iglobal=1,neq);
999 FORMAT (11(:,1X,F8.5) )
     endsubroutine glbmatrcs


!     ===============================================================
!                       the element coefficient matrix
!     ===============================================================
      subroutine elmatrcs_1d(elu,elk,elq)
!     ================================================================
!     element  calculations based on  linear and quadratic rectangular
!     relements with isoparametric  formulation.
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      real*8,intent(in):: elu(:) !dof vector
!     ================================================================
!                          otput variables
!     ================================================================
      real*8,intent(out):: elq(:) !element source vector
      real*8,intent(out):: elk(:,:) !the element coefficient matrix
!     ================================================================
      integer::ni,nj,nk,i,j,ii,jj,k,li,lj
      integer::inode,jnode
!     ================================================================
      real*8:: elu_tense(npe,ndf) !dof current time
      real*8:: elp_tense(npe,ndf) !dof current time
      real*8::coord(dimen) !the global coordinate of each node
      real*8::xi(dimen)! transformed coordinates
      real*8::sf(npe),gdsf(dimen,npe),ur(ndf,dimen)
      real*8::k_coef(ndf,ndf) ! the auxilary matrix for cuefficitne matr
      real*8::res_vect(ndf)
      real*8::cnst,det
      real*8::der(dimen) !electric displacement
      real*8::sigma(dimen,dimen) !stress
!     ================================================================
      real*8:: gauspt(5,5),gauswt(5,5) ! gauss points and gauss weights
!     ================================================================
       data gauspt/5*0.0d0, -0.57735027d0, 0.57735027d0, 3*0.0d0,            &
       -0.77459667d0, 0.0d0, 0.77459667d0, 2*0.0d0, -0.86113631d0,           &
       -0.33998104d0, 0.33998104d0, 0.86113631d0, 0.0d0, -0.90617984d0,      &
       -0.53846931d0,0.0d0,0.53846931d0,0.90617984d0/
!
      data gauswt/2.0d0, 4*0.0d0, 2*1.0d0, 3*0.0d0, 0.55555555d0,        &
        0.88888888d0, 0.55555555d0, 2*0.0d0, 0.34785485d0,               &
      2*0.65214515d0, 0.34785485d0, 0.0d0, 0.23692688d0,                 &
        0.47862867d0, 0.56888888d0, 0.47862867d0, 0.23692688d0/
!     ===============================================================
      elq=0.0d0
      elk=0.0d0

!     ===============================================================

    do i=1,npe
    elu_tense(i,:)=elu(1+(i-1)*ndf:i*ndf)
    enddo !i=1,npe

    do 200 ni = 1,ipdf
    xi(1) = gauspt(ni,ipdf)

      call shapefuncs_1d(xi,det,sf,gdsf)
      cnst = det*gauswt(ni,ipdf)
      coord=0
    do i=1,npe
    coord(:)=coord(:)+elcrds(i,:)*sf(i)
    enddo !i=1,npe
!    write(1,*)'coord=',coord
!     ===============================================================
!                       post process
!     ===============================================================

    ur=transpose(matmul(gdsf,elu_tense))
    call material_properties()
    call stress_elect_displacement(ur,der,sigma)

!     ===============================================================
!                       filling the elemt matrixes
!     ===============================================================
    ii=1
      do 180 inode=1,npe
!     inode,jnode : are the node iterator
      call residual(inode,gdsf,ur,der,sigma,res_vect)
      do 162 li=1,ndf
!     li, lj : are the degree of freedom iterator
 162  elq(ii+li-1)=elq(ii+li-1)+res_vect(li)*cnst
      jj=1
      do 160 jnode=1,npe
      call  k_gen(inode,jnode,gdsf,sf,ur,k_coef)
      do 161 li=1,ndf
      do 161 lj=1,ndf
161   elk(ii+li-1,jj+lj-1)=elk(ii+li-1,jj+lj-1)+k_coef(li,lj)*cnst

!     ii,jj: define the location of information in the element matrix

  160 jj = ndf*jnode+1
  180 ii = ndf*inode+1

  200 continue
!write(1,*)'elq=',elq

      endsubroutine elmatrcs_1d


!     ===============================================================
!                       THE SHAPE FUNCTION SUBROUTINE
!     ===============================================================
      subroutine shapefuncs_1d(xi,det,sf,gdsf)
      real*8::det ! the jacobian determinant
      real*8,intent(out) ::sf(:),gdsf(:,:) ! shape function sf(npe),gdsf(dimen,npe)

      real*8::dsf(dimen,npe) ! the differentiatian of shape function
      real*8::gj(dimen,dimen)! global jacobian matrix
      real*8::gjinv(dimen,dimen) ! inverse of global jacobian matrix
      real*8::xnode(20,3) ! coordinate of nodes in refrence
      real*8::xp(dimen),xi(dimen) !locad coordinate in element
      integer::i
!     ================================================================
      data xnode/                                         &
     -1,1,1,-1,-1,1,1,-1,0,1,0,-1,-1,1,1,-1,0,1,0,-1,   &
     -1,-1,1,1,-1,-1,1,1,-1,0,1,0,-1,-1,1,1,-1,0,1,0,   &
     -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,0,0,0,0,1,1,1,1/
!     ================================================================
      do 10 i=1,npe
      xp(1:dimen)=xnode(i,1:dimen)
      sf(i)=0.5*(1+xi(1)*xnode(i,1))
      dsf(1,i)=0.5*(xnode(i,1))
10    continue

!     ===================================================================
!         compute the jacobian matrix [gj] and its inverse [gjinv]
!     ===================================================================
      gj=matmul(dsf,transpose(elcrds))
      det=dabs(gj(1,1))
      gjinv=1/gj(1,1)
      gdsf=matmul(gjinv,dsf)
      end subroutine shapefuncs_1d

!     ===============================================================
!                       residula vector subroutine
!     ===============================================================
      subroutine residual(inode,gdsf,ur,der,sigma,res_vect)
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      integer,intent(in):: inode !node number
      real*8,intent(in)::gdsf(dimen,npe) !  the shape functions
      real*8,intent(in),dimension(ndf,dimen):: ur
      real*8,intent(in),dimension(dimen)::der ! electric displacement tensor
      real*8,intent(in),dimension(dimen,dimen)::sigma !stress tensor
!     ================================================================
!                          output variables
!     ================================================================
      real*8,intent(out)::res_vect(ndf) !the coefficient matri
!     ================================================================
      integer::i,j,k,l,m,n,p,s,t
      integer::del(dimen,dimen),eye(ndf,ndf)
!     ================================================================
!     deformmation and strain tensors
!     ================================================================
      real*8,dimension(dimen,dimen)::udsp_t,udsp_p,udsp_b
!     ===========================shape function tensors
      real*8::bepsilon(dimen,dimen,ndf),belec(dimen,ndf)
!     ================================================================
      del = 0.0; do i = 1,dimen; del(i,i) = 1;enddo
      eye = 0.0; do i = 1,ndf; eye(i,i) = 1;enddo
!!     ================================================================
!!            c     shape function tensor
!!     ================================================================
      belec=0.0d0;
      bepsilon=0.0d0;
!

      belec(:,dimen+1) = -gdsf(:,inode)
      do 10 k=1,ndf
      do 10 i=1,dimen
      do 10 j=1,dimen
      bepsilon(i,j,k)= 0.5d0*gdsf(j,inode)*eye(k,i)  &
                     + 0.5d0*gdsf(i,inode)*eye(k,j)
      do 10 m=1,dimen
      bepsilon(i,j,k)=bepsilon(i,j,k)                       &
                   + 0.5d0*gdsf(i,inode)*eye(k,m)*ur(m,j)   &
                   + 0.5d0*gdsf(j,inode)*eye(k,m)*ur(m,i)
 10   continue
!write(1,*)'gdsf=',gdsf

!     ================================================================
      res_vect=0.0d0

      do 40 k=1,ndf
      do 40 i=1,dimen
      res_vect(k)=res_vect(k)-der(i)*belec(i,k)
      do 40 j=1,dimen
40    res_vect(k)=res_vect(k)+sigma(i,j)*bepsilon(i,j,k)

!write(1,*)'belec=',belec
!write(1,*)'bepsilon=',bepsilon


!write(1,*)'sigma=',sigma(1,1)
!write(1,*)'der=',der(1)
      endsubroutine residual

!     ===============================================================
!                       coefficient matrix subroutine
!     ===============================================================
      subroutine  k_gen(inode,jnode,gdsf,sf,ur,k_coef)
!     ================================================================
      implicit none
!     ================================================================
!                          input variables
!     ================================================================
      integer,intent(in)::inode,jnode !node number iterators
      real*8,intent(in)::sf(:),gdsf(:,:) !  the shape functions
      real*8,intent(in)::ur(:,:)! displacement
!     ================================================================
!                          output variables
!     ================================================================
      real*8,intent(out)::k_coef(:,:) !the coefficient matrix
!     ================================================================
      integer::i,j,k,l,m,n,p,s,t
      integer,save::counter
      integer::del(dimen,dimen),eye(ndf,ndf)
      real*8::dsgr(dimen,dimen)

!     electric field variables
!     ================================================================
      real*8,dimension(dimen)::der,elr_t,elr_p,elr_b
      real*8,dimension(dimen)::d_elr_t,d_elr_p,d_elr_b
!     ===========================shape function tensors
      real*8::bepsiloni(dimen,dimen,ndf),beleci(dimen,ndf)
      real*8::bepsilonj(dimen,dimen,ndf),belecj(dimen,ndf)
!!     ================================================================
      counter=counter+1
      del = 0.0; do i = 1,dimen; del(i,i) = 1;enddo
      eye = 0.0; do i = 1,ndf;   eye(i,i) = 1;enddo
!!     ================================================================
!!            c     material properties
!!     ================================================================
      ctenst=ctens;

      beleci=0.0d0;
      belecj=0.0d0;
      bepsiloni=0.0d0;
      bepsilonj=0.0d0;

      beleci(:,dimen+1) = -gdsf(:, inode)
      belecj(:,dimen+1) = -gdsf(:, jnode)

      do 10 k=1,ndf
      do 10 i=1,dimen
      do 10 j=1,dimen

      bepsilonj(i,j,k)= 0.5*gdsf(j,jnode)*eye(k,i)+0.5*gdsf(i,jnode)*eye(k,j)

      bepsiloni(i,j,k)= 0.5*gdsf(j,inode)*eye(k,i)+0.5*gdsf(i,inode)*eye(k,j)

      do 10 m=1,dimen

      bepsilonj(i,j,k)=bepsilonj(i,j,k)                   &
                   + 0.5*gdsf(i,jnode)*eye(k,m)*ur(m,j)   &
                   + 0.5*gdsf(j,jnode)*eye(k,m)*ur(m,i)
      bepsiloni(i,j,k)=bepsiloni(i,j,k)                   &
                   + 0.5*gdsf(i,inode)*eye(k,m)*ur(m,j)   &
                   + 0.5*gdsf(j,inode)*eye(k,m)*ur(m,i)

 10   continue
!       write(*,*)ktense
!!     ===========================componenst of tangent matrixe
      k_coef=0.0d0;
!!     ===========================the tangent matrix
      do 40 k=1,ndf;
      do 40 l=1,ndf;

      do 40 i=1,dimen;
      do 40 j=1,dimen;

      k_coef(k,l)=k_coef(k,l)-beleci(i,l)*ktense(i,j)*belecj(j,k)

      do 40 m=1,dimen;

      k_coef(k,l)=k_coef(k,l)-epz_t(m,i,j)*beleci(m,l)*bepsilonj(i,j,k)
      k_coef(k,l)=k_coef(k,l)-epz_t(m,i,j)*belecj(m,k)*bepsiloni(i,j,l)

      do 40 n=1,dimen;

40    k_coef(k,l)=k_coef(k,l)+ctenst(i,j,m,n)*bepsilonj(i,j,k)*bepsiloni(m,n,l)

      endsubroutine k_gen

   !     ===============================================================
!                  the boundary condition subroutine
!     ===============================================================
     subroutine symmetric_primary_bounday(glk,glr)
     implicit none
!     ================================================================
!                          input variables
!     ================================================================
      real*8,intent(inout) :: glk(:,:),glr(:)! global coefficient
      integer::i !integer counters
      real*8,allocatable::gls(:)
!     ================================================================
!     ================================================================
!                          primary variables
!     for u(i)=alpha, it puts k(i,j)=1 and f(i)=alpha
!     ================================================================
      allocate(gls(size(glr)));
        gls=0.0d0
        gls(bnd_no_pr_vec)=vspv
        glr=glr-matmul(glk,gls)

        glk(bnd_no_pr_vec,:)=0.0d0;
        glk(:,bnd_no_pr_vec)=0.0d0;

        do i=1,nspv;
        glk(bnd_no_pr_vec(i),bnd_no_pr_vec(i))=1.0d0; enddo
        glr(bnd_no_pr_vec)=vspv


    endsubroutine symmetric_primary_bounday


subroutine result_printer_1D(iter,glu)
      real*8,INTENT(IN)::glu(:);
      integer::iter;
!     ================================================================
!                          trivial variables
!     ================================================================
      integer::i,j,pdf,k;
      integer::curved_node;
      integer,parameter::gnuplot=5
      integer,save::counter
      character(len=5)::str,x1
      character(len=20) :: fmt,filename ! format descriptor
      fmt = '(i4.4)' ! an integer of width 5 with zeros at the left

!       write(*,*)dimen

      curved_node=nnm
      counter=counter+1


      i=1
      write(x1,fmt)i
      filename='1d_afc'//trim(x1)//'.csv'
!
      open (out,file='out_1d_data_outpu.txt')
      open (gnuplot,file='gnu_out_put.gnu')
!
      write(out,910);write(out,*)'time',time
      write(out,910);write(out,910);write(out,670);write(out,910)
      do i=1,nnm
      pdf=(i-1)*ndf
      write(out,950)i,(coords(i,j),j=1,dimen),(glu(pdf+k),k=1,ndf)
      enddo
      if(counter.eq.1)then
      open (csv,file=filename)
      write(csv,672)
      endif
      pdf=(curved_node-1)*ndf
      write(csv,951)iter,time(2),glu(pdf+1)
      write(gnuplot,*)glu(pdf+2),time(2),glu(pdf+1)

!     ===============================================================
!                              formats
!     ===============================================================
910   format (2x,80('_'),/)
670   format(5x,'node     x-coord. u-solut.')
950   format(5x,i5,10e14.5)
951   format(5x,i5,' , ',10(e14.5,' , '))
672   format(5x,'itr, time[sec],u[m],  elec_pot[v]')

end subroutine result_printer_1D

      end module fem_libs


