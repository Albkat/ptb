module purify 
   use timing_utilities
   use iso_fortran_env, only : wp => real64, out => output_unit
   use checks
   use cuda_
   use cli, only : pur
   use accel_lib 
   implicit none
   public :: purification 
   private 
contains

!> non-orthogonal purification
subroutine purification(H, ndim, S, P, P_purified,method,&
      & verbose,chempot,metr,fix,upper_limit,step,iter,cycle,sparse_cli,cuda)
   
   use timing_utilities
   implicit none
   
   !> number of basis functions
   integer, intent(in)      :: ndim                  
   
   !> overlap matrix packed
   real(wp), intent(in)     :: S(ndim*(ndim+1)/2) 
   
   !> density matrix packed
   real(wp), intent(in)     :: P(ndim*(ndim+1)/2)
   
   !> Hamiltonian matrix packed
   real(wp), intent(in)     :: H(ndim*(ndim+1)/2) 
   
   !> purified density matrix
   real(wp), intent(out)    :: P_purified(ndim*(ndim+1)/2)

   !> settigns
   character(len=*), intent(in) :: method
   
   !> intial guess
   real(wp), intent(inout) :: chempot

   !> if verbose
   logical, intent(in) :: verbose

   !> metric
   character(len=*), intent(in) :: metr   

   !> fixed chemical potential
   logical, intent(in) :: fix

   !> step size
   real(wp),intent(in) :: step

   !> upper limit for chempot scan
   real(wp),intent(in) :: upper_limit

   !> fixed chemical potential
   logical, intent(in) :: iter

   !> number of purification cycles
   integer, intent(in) :: cycle

   !> sparse 
   logical, intent(in) :: sparse_cli

   !> CUDA support
   logical, intent(in) :: cuda


   ! local vars !

   !> density matrix 
   real(wp) :: Psym(ndim,ndim)
   
   !> overlap matrix
   real(wp) :: Ssym(ndim,ndim)
   
   !> hamiltonian matrix 
   real(wp) :: Hsym(ndim,ndim)
   
   !> metric used for initial guess
   real(wp) :: metric (ndim,ndim)
   
   !> cash density matrix 
   real(wp) :: tmp(ndim,ndim)
   
   !> S^(1/2)
   real(wp) :: root(ndim,ndim)

   !> S^(-1/2)
   real(wp), dimension(ndim,ndim) :: root_inv
   
   !> upper and lower spectrum bounds of H 
   real(wp) :: hmax, hmin 

   !> debug mode
   logical, parameter :: time(2) = [.false.,.true.]
   
   
   !-------!
   ! SETUP !
   !-------!
   hmax  = 0.0_wp
   hmin  = 0.0_wp
   Hsym(:,:)=0.0_wp
   Ssym(:,:)=0.0_wp
   Psym(:,:)=0.0_wp

    
   ! transform array to matrix ! 
   call blowsym(ndim,P,Psym)
   call blowsym(ndim,H,Hsym)  
   call blowsym(ndim,S,Ssym)
   
   ! max and min eigenvalues !
   if (method .eq. "mcweeny".or. method .eq. "cp") &
      call eigs(ndim, Hsym, hmax, hmin)
   
   ! obtain S^(-1) !
   if (metr .eq. "-1") then
      
      call inverse_(ndim,Ssym,metric)
   
   ! obtain the S^(1/2) !
   else if (metr == "0.5") then

      call sq_root(ndim,Ssym,root)
      metric = root

   ! obtain the S^(-1/2) !
   else if (metr == "-0.5") then
      
      call sq_root(ndim,Ssym,root)
      call inverse_(ndim,root,root_inv)
      metric = root_inv
   
   endif
   
   
   ! print the purification settings !
   write(out,'(/,a,1x,a,1x,a,/)') repeat('*',24), "PURIFICATION", repeat('*',25)
   write(out,'(2x,a)') "__SETTINGS__" 
   write(out,'(2x,a,5x,a)') "Purification type:            ", method
   if (method.eq."sign") then
      if (iter) then
         write(out,'(2x,a,5x,a)') "Sign mode:                    ", "iterative"
      else
         write(out,'(2x,a,5x,a)') "Sign mode:                    ", "diagonalisation"
      endif
   endif
   write(out,'(2x,a,5x,a)') "Overlap metric (S^metric):    ", metr 
   write(out,'(2x,a,5x,I0)') "Purification cycle number:    ", cycle

   if (method .ne. 'cp') &
      &  write(out,'(2x,a,2x,f14.8)') "Initial chemical potential:   ", chempot
   if (.not.fix .and. method.ne."cp") then
      write(out,'(2x,a,2x,f14.8)') "Upper limit               :   ", upper_limit
      write(out,'(2x,a,2x,f14.8)') "Step (initial -> limit)   :   ", step
   else if(method.ne."cp") then
      write(out,'(2x,a,6x,L1)') "Fixed mode:                   ",fix
   endif
   if (verbose) & 
      & write(out,'(2x,a,6x,L1)') "Verbose:                      ", verbose
   write(out,'(2x,a,5x,l)') "Cuda support:                 ", pur%cuda

   write(out,'()')

!------------------------------------------!
!-------------- Purification --------------!
!------------------------------------------! 
      
   
   if (method .eq. "cp")then
      call cp_purification(ndim, metric, Ssym, Psym, Hsym, hmax, hmin, tmp, cycle, verbose)
   
   else
      if (method .eq. "mcweeny") then
         call gcp_purification(ndim,metric,Ssym,Hsym,hmax,hmin,chempot, &
            & upper_limit,step,tmp,cycle,fix,verbose)
      else
         call sign_purification(ndim,metric,Ssym,Hsym,chempot, &
            & upper_limit,step,tmp,verbose,fix,iter,cycle,sparse_cli)
      endif
   endif

   call packsym(ndim,tmp,P_purified)


end subroutine purification

!> Sign Method 
subroutine sign_purification(ndim, metric, S, H, chempot, upl, &
      & step, P_purified, verbose, fix, iterative, cycle_,sparse_cli)
   
   use ieee_arithmetic, only : ieee_is_NaN
   use gtb_lapack_eig, only : la_sygvd
   use gtb_la, only : la_gemm, la_symm

   !> number of basis functions
   integer, intent(in) :: ndim

   !> upper limit for chempot scan
   real(wp), intent(in) :: upl

   !> metric used for Hamiltonian diagonalization
   real(wp), intent(in) :: metric(ndim,ndim)

   !> hamiltonian matrix 
   real(wp), intent(in) :: H(ndim,ndim)
   
   !> overlap matrix 
   real(wp), intent(in) :: S(ndim,ndim)
   
   !> chemical potential
   real(wp), intent(inout) :: chempot
   
   !> step size for chempot scan
   real(wp), intent(in) :: step

   !> purified density matrix 
   real(wp), intent(out) :: P_purified(ndim,ndim)
   
   !> verbosity
   logical, intent(in) :: verbose

   !> if print the iterations inside purification cycle
   logical, intent(in) :: fix

   !> iterative or diagonalization
   logical,intent(in) :: iterative 
   
   !> cycles
   integer,intent(in) :: cycle_

   !> farced sparse mode
   logical, intent(in) :: sparse_cli

   !> iterators
   integer :: i
   
   !> normalization
   real(wp) :: frob, gersh, Ascale
   
   !> terms for initial guess 
   real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, tmp
   
   !> initial guess
   real(wp) :: sign_A(ndim,ndim)
   
   !> final guess
   real(wp) :: final_(ndim,ndim)
   
   !> I
   real(wp) :: identity(ndim,ndim)
   
   !> if converged
   logical :: conv

   !> if sparse
   logical :: sparse

   !> if NaN
   logical :: NaN_

   !> if debug mode
   logical :: debug

   !> err
   integer :: err

!---------------------------------------------------!
!-------------- initial configuration --------------!
!---------------------------------------------------! 

   debug = verbose .and. ndim < 20

   ! identity matrix !
   identity=0.0_wp
   do i=1,ndim
      identity(i,i) = 1.0_wp
   enddo

!---------------------------------------------------!
!------------------ chempot scan -------------------!
!---------------------------------------------------!

   ! find optimal chemical potential !
   chemp: do while (chempot <= upl) 
         
      ! sign(A) ! 
      term1 = chempot * identity
      
      if (pur%cuda) then
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, H, metric, 0.0_wp,term2,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, metric, term2, 0.0_wp,term3,err)
      else
         !call la_symm('L','U',ndim,ndim,1.0_wp,H,ndim,metric,ndim,0.0_wp,term2,ndim) ! H * metric
         call la_gemm('N','N',ndim,ndim,ndim,1.0_wp,H,ndim,metric,ndim,0.0_wp,term2,ndim) ! metric * H * metric 
         call la_gemm('N','N',ndim,ndim,ndim,1.0_wp,metric,ndim,term2,ndim,0.0_wp,term3,ndim) ! metric * H * metric 
      endif 
      
      term4 = term3 - term1
      sign_A = term4 ! symmetric !
         
      ! frobenius norm !
      frob=sqrt(sum(sign_A**2))
      
      ! gershgorin norm !
      gersh=0.0_wp
      do i=1,ndim
        gersh=max(gersh,sum(abs(Sign_A(:,i))))
      enddo

      ! print the normalization !
      if (debug) then
         write(out,'(a, 1x, F12.6)') "Gershgorin norm", gersh
         write(out,'(a, 1x, F12.6)') "Frobenuis norm",  frob
      endif

      Ascale=min(frob,gersh)
      sign_A=sign_A/Ascale

      ! for debugging !
      if (debug) &
         & call print_blowed_matrix(ndim,sign_A,'sign_A')

      ! check sparisty!
      if (.not. allocated(pur%sparse)) then
         call chk_sparsity(ndim,sign_A,verbose,sparse) 
      else 
         sparse = pur%sparse
      endif

      ! Submatrix method !
      if(sparse .and. .not.iterative) then
         
         tmp = sign_A
         call apply_submatrix_method(ndim,tmp,sign_A,verbose)

         ! for debugging !
         if (debug) &
            & call print_blowed_matrix(ndim,sign_A,'after sub')
      
      ! Dense matrix case !
      else

         if (iterative) then
            call sign_iterative(ndim,chempot,sign_A,H,metric,identity,cycle_,fix,final_,conv)                     
         else
            call sign_diagonalization(ndim,sign_A)
         endif
      
      endif
         
      ! construct density matrix !
      call sign_density(ndim,final_, sign_A, identity, metric, NaN_)
      
      ! print iterarion summary !
      call print_summary(ndim,chempot,sign_A,final_,H,identity,iterative,&
                              &  NaN_, conv)
      
      ! exit if chempot is user-fix !
      if (fix) exit chemp
      
      ! next iteration parameters !
      chempot=chempot+step
   
   enddo chemp

   P_purified = final_ * 2.0

end subroutine sign_purification

!> iterative purification of sign matrix
subroutine sign_iterative(ndim, chempot, sign_A, H, metric, identity, cycle_, fix, final_, conv)

   use gtb_la, only : la_gemm
   use iso_fortran_env, only : wp => real64, out => output_unit
   implicit none

   !> number of basis functions
   integer, intent(in) :: ndim
   
   !> chemical potential
   real(wp), intent(in) :: chempot

   !> sign matrix
   real(wp), intent(inout) :: sign_A(ndim,ndim)
   
   !> hamiltonian
   real(wp), intent(in) :: H(ndim,ndim)

   !> metric
   real(wp), intent(in) :: metric(ndim,ndim)

   !> identity matrix
   real(wp), intent(in) :: identity(ndim,ndim)

   !> number of cycles
   integer, intent(in) :: cycle_

   !> if fixed mode
   logical, intent(in) :: fix
   
   !> density matrix
   real(wp), intent(out) :: final_(ndim,ndim)
   
   !> convergence
   logical, intent(out) :: conv
   
   ! local vars !
   logical :: cuda

   !> iterator
   integer :: purificator
   
   !> temporary matrices
   real(wp), dimension(ndim,ndim) :: x2k, x4k, tmp, tmp2
   real(wp) :: norm

   !> if NaN
   logical :: NaN_

   !> error handle
   integer :: err

   conv=.false.
   cuda = pur%cuda 

   ! Pade approximation !
   pur: do purificator=1, cycle_
      
      if (cuda) then
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, sign_A, sign_A, 0.0_wp, x2k, err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, x2k, x2k, 0.0_wp, x4k, err)
      else
         call la_gemm(sign_A,sign_A,x2k) 
         call la_gemm(x2k,x2k,x4k)
      endif

      tmp = (15.0_wp * identity) - (10.0_wp * x2k) + (3.0_wp * x4k)
      
      if (cuda) then 
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 0.125_wp, sign_A, tmp, 0.0_wp, tmp2, err)
      else
         call la_gemm(sign_A,tmp,tmp2, alpha=0.125_wp)
      endif
      
      sign_A = tmp2
      
      call sign_density(ndim,final_, sign_A, identity, metric, NaN_) ! create density matrix from sign matrix
      if (NaN_) exit pur ! exit the current loop if NaN is found

      ! convergence check ! 
      if (purificator>1) then
         if (abs(norm2(sign_A)-norm) < 1.0e-5) then
            conv=.true.
            write(out,'(2x,a,I0)') "cycle number: ", purificator
            exit pur
         else
            norm=norm2(sign_A)
         endif
      else
         norm=norm2(sign_A)
      endif

      ! printout for fix modus !
      if (fix) then

         if (purificator==1) &
            write(out,'(2x,a)') "_ITERATIVE PURIFICATION_" 
         
         call print_summary(ndim,chempot,sign_A,final_,H,identity,.true.,NaN_,.true.)
         
         if (purificator==cycle_) &
            write(out,'(2x,a,I0)') "cycle number: ", purificator
      endif

   enddo pur 

end subroutine sign_iterative

!> construct density matrix from sign matrix
subroutine sign_density(ndim, P, sign_A, identity, metric, NaN_)

   use gtb_la, only : la_gemm
   use ieee_arithmetic, only : ieee_is_NaN
   use iso_fortran_env, only : wp => real64
   implicit none

   !> number of basis functions
   integer, intent(in) :: ndim

   !> density matrix
   real(wp), intent(out) :: P(ndim,ndim)

   !> sign matrix
   real(wp), intent(in) :: sign_A(ndim,ndim)

   !> identity matrix
   real(wp), intent(in) :: identity(ndim,ndim)

   !> metric
   real(wp), intent(in) :: metric(ndim,ndim)

   !> NaN
   logical, intent(out) :: NaN_

   !> temporary matrices
   real(wp) :: tmp(ndim,ndim)

   integer :: err
   

   ! (16) density matrix !
   if (pur%cuda) then 
      call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, identity-sign_A, metric, 0.0_wp,tmp,err)
      call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 0.5_wp, metric,tmp, 0.0_wp,P,err)
   else
      call la_gemm(identity-sign_A,metric,tmp)
      call la_gemm(metric,tmp,P,alpha=0.5_wp)
   endif
   NaN_ = any((ieee_is_NaN(P))) ! check for NaN
      
end subroutine sign_density

!> purification of the sign matrix via diagonalization
subroutine sign_diagonalization(ndim, sign_A)

   use gtb_la, only : la_gemm
   use iso_fortran_env, only : wp => real64, out => output_unit, int => int64
   implicit none

   !> number of basis functions
   integer, intent(in) :: ndim

   !> sign matrix
   real(wp), intent(inout) :: sign_A(ndim,ndim)

   !> eigenvectors
   real(wp) :: eigvec(ndim,ndim)

   !> eigenvalues
   real(wp) :: eigval(ndim), eigval2(ndim,ndim)

   !> temporary vars
   integer :: i
   real(wp) :: tmp(ndim,ndim)

   integer :: err

   call eigendecompostion(ndim, sign_A, eigvec, eigval)
            
   eigval2 = 0.0_wp
   tmp = 0.0_wp
   
   ! signum() !
   do i = 1, ndim
      if (abs(eigval(i)) < 0.1E-10_wp) then
         eigval2(i,i) = 0.0_wp
      else 
         if (eigval(i) > 0.0) then 
            eigval2(i,i) = 1.0_wp 
         else 
            eigval2(i,i) = -1.0_wp 
         endif
      endif
   enddo

   if(pur%cuda) then
      call cuda_dgemm(ctx, 'N', 'T', ndim, ndim, ndim, 1.0_wp, eigval2, eigvec, 0.0_wp,tmp,err)
      call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, eigvec, tmp, 0.0_wp,sign_A,err)
   else
      call la_gemm(eigval2,eigvec,tmp,transb='T')
      call la_gemm(eigvec,tmp,sign_A)
   endif 

end subroutine sign_diagonalization

!> grand-canonical purification
subroutine gcp_purification(ndim, metric, S, H, hmax, hmin, chempot, &
      & upl, step, P_purified, cycle_, fix, verbose)
   
   use gtb_la, only : la_gemm
   use ieee_arithmetic, only : ieee_is_NaN
   use iso_fortran_env, only : wp => real64, out => output_unit 

   !> number of basis functions
   integer, intent(in) :: ndim
   
   !> upper limit for chempot scan
   real(wp), intent(in) :: upl

   !> metric (S version)
   real(wp), intent(in) :: metric(ndim,ndim)
   
   !> hamiltonian matrix
   real(wp), intent(in) :: H(ndim,ndim)
   
   !> overlap matrix
   real(wp), intent(in) :: S(ndim,ndim)
   
   !> upper and of the eigenvalue spectrum of H
   real(wp), intent(in) :: hmax, hmin 

   !> chemical potential
   real(wp), intent(inout) :: chempot
   
   !> step size
   real(wp), intent(in) :: step   

   !> purified density 
   real(wp), intent(out) :: P_purified(ndim,ndim)

   !> max purification number
   integer, intent(in) :: cycle_

   !> if chempot user-defined
   logical, intent(in) :: fix

   !> printlevel
   logical, intent(in) :: verbose

!----------------------------------------------------------
   
   !> iterators
   integer :: num, purificator, i, j, ic, jc, iter
   
   !> scaling parameters
   real(wp) :: alpha_max, alpha_min, alpha
   
   !> l2 matrix norm
   real(wp) :: norm 
   
   !> band structure energy
   real(wp) :: band 
   
   !> number of electrons
   real(wp) :: nel, number_of_electrons 
   
   !> normalization
   real(wp) :: frob, gersh, Ascale
   
   !> terms for (11b)
   real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, term5, term6
   
   !> terms for (10) McWeeny purification 
   real(wp), dimension(ndim,ndim) :: SP, PSP, PSPSP, res
   
   !> initial guess
   real(wp) :: P0(ndim,ndim)
   
   !> identity matrix
   real(wp) :: identity(ndim,ndim)
   
   !> P*H
   real(wp) :: PH(ndim,ndim)
   
   !> NaN catch
   logical :: error
   
   !> if converged
   logical :: conv
   
   integer :: err
   
   logical :: cuda
!-----------------------! 
! initial configuration !
!-----------------------! 
   
   error = .false.
   conv  = .false. 

   ! construct identity !
   identity =0.0_wp
   do i=1,ndim
      identity(i,i) = 1.0_wp
   enddo
   cuda = pur%cuda
!--------------!
! chempot scan !
!--------------!

   chemp: do while (chempot <= upl)
   
      ! construct initial guess !   
      ! (11a) !
      term1 = chempot*metric

      if (pur%cuda) then

         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, H, metric, 0.0_wp,term2,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, metric, term2, 0.0_wp,term3,err)
      else 
         call la_gemm(H,metric,term2)
         call la_gemm(metric,term2,term3)

      endif

      term4 = term1 - term3            
      
      ! frobenius norm !
      frob=sqrt(sum(term4**2))
   
      ! gershgorin norm !
      gersh=0.0_wp
      do i=1,ndim
         gersh=max(gersh,sum(abs(term4(:,i))))
      enddo
      
      ! print norms !
      if (verbose) then
         write(out,'(2x,a,f14.9)'), "Gershgorin Norm  = ", gersh
         write(out,'(2x,a,f14.9)'), "Frobenius Norm   = ", frob
      endif

      alpha=1.0_wp/min(gersh,frob)

      term5 = alpha * 0.5_wp * term4   ! scaling !
      P0 = term5 + 0.5_wp * metric     ! shift !
      
      purify: do purificator=1,pur%cycle
         if(cuda) then 
            call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, S, P0, 0.0_wp,SP,err)
            call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, P0, SP, 0.0_wp,PSP,err)
            call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, PSP, SP, 0.0_wp,PSPSP,err)
         else
         ! McWeeny purification !
            call la_gemm(S,P0,SP)
            call la_gemm(P0,SP,PSP)
            call la_gemm(PSP,SP,PSPSP)
         endif

         res=3*PSP-2*PSPSP

         ! band-str energy !
         if(cuda) then 
            call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, res, H, 0.0_wp,PH,err)
         else
            call la_gemm(res,H,PH)
         endif

         band=0.0_wp
         
         ! NaN catch !
         do ic=1,ndim
            do jc=1, ndim
               
               ! if NaN !
               if (ieee_is_NaN(res(ic,jc))) then
                  error=.true.
                  exit purify
               endif
               
               ! band-structure energy !
               if (ic == jc) then
                  band=band+PH(ic,jc)
               endif

            enddo
         enddo
         
         if (.not.fix) then
            
            if (purificator>1) then 

               if (abs(norm2(res)-norm)<1.0E-6_wp) then
                  conv=.true.
                  iter=purificator
                  exit purify
               else
                  norm=norm2(res)
               endif

            else
               norm=norm2(res)               
            endif
         
         else
            
            nel=check_electrons(ndim,S,res)

            if (purificator==1) write(*,'(a)') "Fixed chempot purification"
            write(out,111) "chempot = ", chempot, ", iterations = ", purificator, ", electrons = ", nel*2 ,", Band-str E = ", band*2,", norm2=", norm2(res)
            111 format(a,F14.9,a,I0,a,F14.9,a,F14.9,a,F14.9)
         
         endif
         
         P0=res
      
      enddo purify


      if (.not.fix) then
         if (error) then
            write(out,104) "chempot = ", chempot, ', NaN '
            104 format(a,F14.9,a)
         else 
            if (conv) then
               nel=check_electrons(ndim,S,res)
               write(out,112) "chempot = ", chempot, ", iterations = ", iter, ", electrons = ", nel*2.0_wp ,", Band-str E = ", band
               112 format(a,F14.9,a,I0,a,F14.9,a,F14.9)
            else
               write(out,113) "chempot = ", chempot, ', not converged '
               113 format(a,F14.5,a)
            endif
         endif
      endif

      ! parameters for next iteration !
      chempot = chempot + step
      error=.false.
      conv=.false.
      
      ! exit after 1 iter if density is given !
      if (fix) exit chemp
   
   enddo chemp

   P_purified = res * 2.0_wp
   
end subroutine gcp_purification


subroutine cp_purification(ndim, metric, S, P_ptb, H, hmax, hmin, P_purified, cycle_, ver)
   
   use gtb_la, only : la_gemm
   use ieee_arithmetic, only : ieee_is_NaN
   use iso_fortran_env, only : wp => real64, out => output_unit
   implicit none 
   
   !> dimensionality (nbf)
   integer, intent(in) :: ndim
   
   !> metric
   real(wp), intent(in) :: metric(ndim,ndim)

   !> overlap matrix 
   real(wp), intent(in) :: S(ndim,ndim)
   
   !> PTB density matrix 
   real(wp), intent(in) :: P_ptb(ndim,ndim)
   
   !> hamiltonian 
   real(wp), intent(in) :: H(ndim,ndim)
   
   !> spectral bounds  of H
   real(wp), intent(in) :: hmax, hmin 
   
   !> purified density
   real(wp), intent(out) :: P_purified(ndim,ndim)

   !> number of cycles
   integer, intent(in) :: cycle_
   
   !> verbosity
   logical, intent(in) :: ver

!---------------------------------------------------------
   
   !> iterators
   integer :: num, purificator, i, j, iter
   
   !> scaling parameter
   real(wp) :: alpha_max, alpha_min, alpha
   
   !> Gershgorin's circle theorem 
   real(wp) :: h_max, h_min, buff_min, buff_max
   
   !> l2 norm 
   real(wp) :: norm 
   
   !> trace
   real(wp) :: tr, tr1, tr2
   
   !> TODO
   real(wp) :: hmin1, hmax1
   
   !> purification factor 
   real(wp) :: cn

   !> chemical potential 
   real(wp) :: chempot
   
   !> norms 
   real(wp) :: gersh, frob
   
   !> desired number of electrons
   real(wp) :: nel, N_e, N
   
   !> equation terms
   real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, term5
   
   !> McWeeny purification terms
   real(wp), dimension(ndim,ndim) :: SP, PSP, PSPSP, res 
   
   !> temporary matrices
   real(wp) :: matr1(ndim,ndim),matr2(ndim,ndim)
   
   !> initial guess
   real(wp) :: p0(ndim,ndim)
   
   !> density in iteration 
   real(wp) :: pn(ndim,ndim)
   
   !> identity matrix
   real(wp) :: identity(ndim,ndim)
   
   !> if NAN present
   logical :: error
   
   !> if purification converged
   logical :: conv

   !> if frobenius/gershgorin norm should be used
   logical :: frob_norm
   
   integer :: err
   
   !-----------------------!
   ! INITIAL CONFIGURATION !
   !-----------------------!
   
   frob_norm   = .true.
   error       = .false.
   
   N_e = 14.0_wp
   N   = ndim

   ! chempot !
   tr  = 0.0_wp
   h_min  = 0.0_wp
   h_max  = 0.0_wp
   do i=1,ndim
      tr = tr + H(i,i)
      buff_min = H(i,i)
      buff_max = H(i,i)
      do j=1,ndim
         if (i.ne.j) then
            buff_min = buff_min - abs(H(i,j))
            buff_max = buff_max + abs(H(i,j))
         endif
      enddo
      if (i.eq.1) then
         h_min = buff_min
         h_max = buff_max
      else 
         if (buff_max > h_max) h_max = buff_max
         if (buff_min < h_min) h_min = buff_min
      endif

   enddo


   chempot = tr / N

   if (ver) then 
      write(out,"(a, 2x, F14.6)") "Chemical potential", chempot   
      write(out,"(a, 2x, F14.6)") "hmax ", hmax
      write(out,"(a, 2x, F14.6)") "hmin ", hmin
      write(out,"(a, 2x, F14.6)") "h_max ", h_max
      write(out,"(a, 2x, F14.6)") "h_min ", h_min
   endif

   identity = 0.0_wp
   do i=1, ndim
      identity(i,i) = 1.0_wp
   enddo
   
   !-------------!  
   ! CALCULATION !
   !-------------!
   ! initial guess !            
   if (frob_norm) then
      
      term1 = chempot * metric
      
      if (pur%cuda) then  
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, metric, H, 0.0_wp,term2,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, term2, metric, 0.0_wp,term2,err)
      else
         call la_gemm(metric,H,term2)
         call la_gemm(term2,metric,term2)
      endif

      term3 = term1 - term2
      
      ! frobenius norm !
      frob=sqrt(sum(term3**2))

      ! gershgorin norm !
      gersh=0.0_wp
      do i=1,ndim
         gersh=max(gersh,sum(abs(term3(:,i))))
      enddo
      
      ! print norms !
      if (ver) then
         write(out, '(2x,a,f14.9)'), "Gershgorin Norm: ", gersh
         write(out, '(2x,a,f14.9)'), "Frobenius Norm:  ", frob
      endif
      
      alpha=N_e/min(gersh,frob)

      term4 =  metric * N_e    
      p0 = ( (term3 * alpha)  + term4 ) / N_e

   else
   
      hmin1 = ( N - N_e) / (chempot - hmin)
      hmax1 =  N_e / (hmax - chempot)
      
      alpha = min(hmin1,hmax1)
      
      term1 = alpha * chempot * metric
      if(pur%cuda) then  
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, metric,H, 0.0_wp,term2,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, alpha, term2,metric, 0.0_wp,term2,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, N_e, identity,metric, 0.0_wp,term3,err)
      else
         call la_gemm(metric,H,term2)
         call la_gemm(term2,metric,term2,alpha=alpha)
         call la_gemm(identity,metric,term3,alpha=N_e)
      endif
      term4 = (term1 - term2 + term3) / N_e ! 1/N 
      p0    = term4 
      
   endif
   
   if (ver) then
      write(out,'(2x,a,f14.9)') , "number of electrons: ", check_electrons(ndim,S,p0)
      call print_blowed_matrix(ndim,p0, "initial guess")
   endif
   pn = p0*2 

   ! purification loop !
   do num=1, cycle_ 

      if(pur%cuda) then  
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, S, pn, 0.0_wp,SP,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, pn,SP, 0.0_wp,PSP,err)
         call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, PSP,SP, 0.0_wp,PSPSP,err)
      else
         call la_gemm(S,pn,SP)
         call la_gemm(pn,SP,PSP)
         call la_gemm(PSP,SP,PSPSP)
      endif
      matr1=PSP-PSPSP
      matr2=pn-PSP

      tr1=0.0_wp
      tr2=0.0_wp

      do i=1, ndim
         tr1=tr1+matr1(i,i)
         tr2=tr2+matr2(i,i)
      enddo
   
      cn=tr1/tr2

      
      if(cn >= 0.5_wp) then
         res=( (1.0_wp - 2.0_wp * cn ) * pn + (1.0_wp + cn) * PSP - PSPSP ) / (1.0_wp-cn)
      endif
      
      nel = check_electrons(ndim, S, res)
      if (num==1) write(*,'(a)') "Purification loop"
      write(out,114) "iterations = ", num, ", electrons = ", nel, ", band-str E = ", bandStrEnergy(ndim,H,res),", norm2=", norm2(res)
      114 format(a, I0, a, F14.6, a, F14.6, a, F14.6)
      
      if (ver) then   
         call print_blowed_matrix(ndim,res,"P")
         write(out,'(a,2x,f14.6)') "The cn factor", cn
      endif
      
      pn=res
   
   enddo
    
   P_purified = pn

end subroutine cp_purification

subroutine eigs(ndim,matrix,maxim,minim)
   
   use gtb_la, only : la_syev
   use iso_fortran_env, only : wp => real64
   
   integer, intent(in)     :: ndim
   real(wp), intent(in)    :: matrix(ndim,ndim)
   real(wp), intent(out)   :: maxim, minim
   
   real(wp), allocatable :: work(:)
      !! workspace
   integer  :: lwork
      !! dimension of workspace
   real(wp) :: w(ndim)
      !! eigenvalues

   real(wp) :: matrix_syev(ndim,ndim), maxim1,minim1, hsumm 
   real(wp) :: hmaxim1,hminim1,hmaxim,hminim
   integer  :: ipiv(ndim), info
   integer :: i,j
   logical :: lapack 
   
   lapack=.true.
   
   if (lapack) then
      matrix_syev=matrix
      if(pur%cuda) then
         call cuda_dsyevd(ctx,ndim,matrix_syev,w,info)
      else
         allocate(work(1))
         lwork=-1
         call la_syev('N','U',ndim,matrix_syev,ndim,w,work,lwork,info)
         lwork=idint(work(1))
         deallocate(work)
         allocate(work(lwork))
         call la_syev('N','U',ndim,matrix_syev,ndim,w,work,lwork,info)
      endif
      maxim = maxval(w)
      minim = minval(w)

   else
    
      !-------------------------------------------------
      !                       (13a/b)
      !-------------------------------------------------

      !> Gershgorin formulas to obtain lower and upper bounds of the spectrum of H 
      do i=1,ndim      
         do j=1,ndim
            if (i.ne.j) then 
               hsumm=hsumm+abs(matrix(i,j))
            else
               hmaxim1=matrix(i,j)
               hminim1=matrix(i,j)
            endif
         enddo

         hmaxim1 = hmaxim1 + hsumm
         hminim1 = hminim1 - hsumm
         hsumm=0.0_wp

         if (i==1) then
            hmaxim=hmaxim1
            hminim=hminim1
         else
            if (hmaxim1>hmaxim) maxim=maxim1
            if (hminim1>hminim) minim=minim1
        endif
      enddo
  
  endif

end subroutine eigs

!> matrix^(1/2)
subroutine sq_root(ndim,matrix,root)
   
   use gtb_la, only : la_syev, la_gemm
   use iso_fortran_env, only : wp => real64
   use accel_lib
   implicit none 
   integer, intent(in) :: ndim
   real(wp), intent(in) :: matrix(ndim,ndim)
   real(wp), intent(inout) :: root(ndim,ndim)
   
   real(wp), allocatable :: work(:)
      !! workspace
   integer  :: lwork
      !! dimension of workspace
   real(wp) :: w(ndim), sqrtW(ndim)
      !! eigenvalues

   real(wp), dimension(ndim,ndim) :: matrix_syev, matrix_syev_T, root2
   real(wp), dimension(ndim,ndim) :: matrix_syev_inv, D, tmp
   real(wp),dimension(ndim) :: left,right
   integer  :: ipiv(ndim), info
   integer :: i,j
   

   matrix_syev=matrix
   if(pur%cuda) then
      call cuda_dsyevd(ctx,ndim,matrix_syev,w,info)
   else
      allocate(work(1))
      w = 0.0_wp
      lwork = -1
      call la_syev('V','U',ndim,matrix_syev,ndim,w,work,lwork,info)
      lwork=idint(work(1))
      deallocate(work)
      allocate(work(lwork))

      ! V, D !
      call la_syev('V','U',ndim,matrix_syev,ndim,w,work,lwork,info)
   endif 
   ! D^(1/2) !
   sqrtW=sqrt(w)
   D = 0.0_wp
   do i=1,ndim
      D(i,i) = sqrtW(i) 
   enddo

   ! V^(T) !
   do i=1,ndim
      do j=1,ndim
         matrix_syev_T(j,i)=matrix_syev(i,j)
      enddo
   enddo

   ! V^(-1) !
   matrix_syev_inv = 0.0_wp
   call inverse_(ndim,matrix_syev,matrix_syev_inv)
   
   ! r = V * D^(1/2) * V^(-1)
   if(pur%cuda) then
      call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, D, matrix_syev_inv, 0.0_wp,tmp,info)
      call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, matrix_syev, tmp, 0.0_wp,root,info)
      call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, root, root, 0.0_wp,root2,info)
   else
      call la_gemm(D,matrix_syev_inv,tmp)
      call la_gemm(matrix_syev,tmp,root)
      call la_gemm(root,root,root2)
   endif

end subroutine sq_root

!> get inverse of the matrix 
subroutine inverse_(ndim,matrix,inverse) 
   
   use iso_fortran_env, only : wp => real64, out => output_unit
   use multicharge_lapack, only : sytri, sytrf, getrf, getri
   use gtb_la, only : la_gemm

   implicit none
   
   !> number of AOs
   integer,  intent(in)    :: ndim                                  
   
   !> A matrix
   real(wp), intent(in)    :: matrix (ndim,ndim)             
   
   !> inverse
   real(wp), intent(out)   :: inverse (ndim,ndim)             
   
   
   !> A*A^(-1) check
   real(wp) :: check (ndim,ndim), tmp (ndim,ndim) 

   !> identity matrix
   real(wp) :: identity (ndim,ndim)
   
   !> lower or upper triangular part
   character(len=1) :: uplo
   
   !> exit status of routine
   integer :: info
   
   integer :: ipiv(ndim)
   integer :: i, j, ic, jc

   !> use lapack to calculate square root
   logical :: lapack
   
   !> symmetric or general matrix
   logical :: sy

   !> spectral norm vs eigenvalue
   logical :: spect_norm
   
   !> max absolute column sum of matrix
   real(wp) :: a_1(ndim)
   
   !> max absolute row sum of matrix
   real(wp) :: a_inf(ndim)
   
   real(wp),allocatable :: eigvals(:), eigvecs(:,:) 

   lapack   = .true.
   check    = 0.0_wp 

   ! check if matrix sym !
   sy = matrix(1,2) == matrix(2,1)

   ! get inverse of matrix using LAPACK !
   if (lapack) then

      if (sy) then
         if(pur%cuda) then
            allocate(eigvecs(ndim,ndim), eigvals(ndim))
            eigvecs=matrix
            call cuda_dsyevd(ctx, ndim, eigvecs, eigvals, info)
            do j = 1, ndim
               do i = 1, ndim
                  inverse(i, j) = 1.0_wp/eigvals(j)*eigvecs(i, j)
               end do
            end do
            call cuda_dgemm(ctx, 'N', 'T', ndim, ndim, ndim, 1.0_wp, inverse, eigvecs, &
                     0.0_wp, inverse, info)
            call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, matrix, inverse, &
                     0.0_wp, check, info)
         else
            inverse=matrix

            call sytrf(inverse, ipiv, info=info, uplo='l')
            if (info == 0) then
                  call sytri(inverse, ipiv, info=info, uplo='l')
                  if (info == 0) then
                     do ic = 1, ndim
                        do jc = ic+1, ndim
                           inverse(ic, jc) = inverse(jc, ic)
                        end do
                     end do
                  end if
            end if
            call la_gemm(matrix,inverse,check)
         endif

      else

         inverse=matrix
         call getrf(ndim,ndim,inverse,ndim,ipiv,info)
         if (info == 0) then
            call getri(inverse, ipiv, info)  
            if (info /= 0) then
               print *, "The inverse matrix cannot be computed"
               stop
            endif
         endif
         call la_gemm(matrix,inverse,check)
      
      endif

      if (any(abs(check)>1.5_wp)) & 
         stop "Not identity, error during inversion"

   else
       
      ! Identity matrix !
      identity = 0
      do i=1,ndim
         identity(i,i) = 1
      enddo
     
      ! initial guess (27) !
      do j=1,ndim
         a_1(j)=sum(abs(matrix(:,j)))
         a_inf(j)=sum(abs(matrix(j,:)))
      enddo

      write (out,'(a,f12.6)') "a1: ", maxval(a_1)  
      write (out,'(a,f12.6)') "ainf: ", maxval(a_inf)
      
      inverse = (1.0_wp/(maxval(a_1) * maxval(a_inf))) * matrix

      
      !> check  
      call la_gemm(inverse,matrix,tmp)
      spect_norm = 1 > maxval(identity - tmp)
      write(out,'(a,l)'),"Is the largest spectral norm of initial guess < 1? ", spect_norm
      tmp=0.0_wp

      !-------------------------------------------------
      !                       (25)
      !-------------------------------------------------
      
      !> The Newton-Schulz iteration
      do i=1,2
         call la_gemm(inverse,inverse,tmp)
         call la_gemm(tmp,matrix,tmp)
         inverse= (2.0_wp * inverse - tmp) 
      enddo

      call la_gemm(matrix,inverse,check)
      
      if (any(abs(check)>1.5_wp)) then 
         write(out,'(a,l)') "Not identity, error during inversion"
         stop   
      endif

   endif 
   
end subroutine inverse_

subroutine eigendecompostion(ndim,matrix,eigvec,eigval)
  
   use gtb_lapack_eig, only : la_syevd
   use gtb_la, only : la_gemm
   use iso_fortran_env, only : wp => real64, out => output_unit
   use accel_lib
   
   !> dimension
   integer, intent(in)  :: ndim
   
   !> original matrix
   real(wp), intent(in) :: matrix(ndim,ndim)
   
   !> eigenvectors
   real(wp), intent(inout):: eigvec(ndim,ndim)
   
   !> eigenvalues
   real(wp), intent(inout):: eigval(ndim)
   
   !> workspace
   real(wp), allocatable :: work(:)
   integer, allocatable  :: iwork(:)
   !> workspace dim
   integer  :: lwork, liwork

   logical, parameter :: debug = .false.

   real(wp), dimension(ndim,ndim) :: matrix_syev, check, tmp
   integer  :: ipiv(ndim)
   integer  :: info
   integer  :: i,j
   integer  :: err

   matrix_syev=matrix
   
   ! CUDA support !
   if (pur%cuda) then
      
      call cuda_dsyevd(ctx,ndim,matrix_syev,eigval,err)
      
      if (err == 0) then
         eigvec = matrix_syev
      else
         stop "An error occured in the CUDA accelerated code."
      end if

   else
      
      allocate(work(1))
      allocate(iwork(1))
      lwork=-1
      
      call la_syevd('V', 'U', ndim, matrix_syev, ndim, eigval, &
               & work, lwork, iwork, liwork, info)
      
      lwork=idint(work(1))
      liwork=iwork(1)
      deallocate(work)
      deallocate(iwork)
      allocate(work(lwork))
      allocate(iwork(liwork))

      call la_syevd('V', 'U', ndim, matrix_syev, ndim, eigval, &
               & work, lwork, iwork, liwork, info)
      
      !-------!
      ! CHECK !         
      !-------!
      if (info==0) then
         
         eigvec = matrix_syev
         
         if (debug) then
            check = 0.0_wp 
         
            do i=1, ndim 
               check(i,i) = eigval(i)
            enddo

            call la_gemm(check,eigvec,tmp,transb='T')
            call la_gemm(eigvec,tmp,check)

            if (any(abs(matrix-check) > 1E-8_wp)) then
               stop "Diagonalization failure!"
            endif
         endif

      else 
         stop "Matrix can not be diagonalized."
      endif

   endif

end subroutine eigendecompostion

!> submatrix method
subroutine apply_submatrix_method(ndim, input_matrix, result, verbose)
   
   use iso_fortran_env, only : wp => real64, int => int64, out => output_unit
   implicit none
   integer, intent(in) :: ndim
   real(wp), dimension(ndim,ndim), intent(in) :: input_matrix
   real(wp), dimension(ndim,ndim), intent(out) :: result
   logical, intent(in) :: verbose

   integer(int) :: N(2), dim, i, j, k,  result_col
   integer :: nonzero
   real(wp), dimension(:,:), allocatable :: submatrix, submatrix_result
   integer(int), dimension(:), allocatable :: indices
   character(len=20), allocatable :: msg
   logical :: debug
   real(wp) :: thr = 1E-6
   
   debug = ndim < 20 .and. verbose
   N = shape(input_matrix)
   result=0.0_wp
   if (verbose) &
      write(out,'(a,2x,i0)') "Applying submatrix method to matrix of dimension", ndim   

   ! We generate and process a submatrix for each column i of the input matrix
   do i = 1, ndim

      ! The size of the submatrix is determined by the number of nonzero elements
      !  in the inducing column i.
      nonzero = count(abs(input_matrix(:,i)) > thr) ! Criteruim ???
      allocate(submatrix(nonzero, nonzero), submatrix_result(nonzero, nonzero), &
               &indices(nonzero))!, eigvec(nonzero,nonzero), eigval(nonzero))


      ! We need to determine which elements in that column are nonzero because
      ! they determine the elements copied into the submatrix.
      k = 1
      do j = 1, ndim
         if (abs(input_matrix(j,i)) > thr) then

            indices(k) = j

            ! We should take note which column in the submatrix contains elements from the i'th column
            ! of the original matrix because that will be our result column.
            ! Note: Elements on the diagonal of the input matrix always have to be nonzero!
            if (j == i) then
               result_col = k
            end if
            
            k = k+1

         end if

      end do
      
      ! Building the submatrix now just means selecting the right rows and columns of the original matrix.
      submatrix = input_matrix(indices,indices)
      
      ! print submatrices if it is small enough !
      if (debug) then
         write(out,'(a,x,i0)') "Submatrix", i
         call print_blowed_matrix(nonzero,submatrix)
         write(out, *)
      end if


      ! Apply the function of interest to the submatrix.
      call sign_diagonalization(nonzero,submatrix)
      
      ! print submatrices if it is small enough !
      if (debug) then
         write(out,'(a,x,i0)') "diagonalized submatrix", i
         call print_blowed_matrix(nonzero,submatrix)
         write(out, *)
      end if

      ! Now we copy the matching result column into the approximate result matrix.
      result(indices,i) = submatrix(:,result_col)

      deallocate(submatrix, submatrix_result, indices)
   end do

end subroutine apply_submatrix_method

!> check if matrix should be decomposed 
subroutine chk_sparsity(ndim,A,verbose,sparse)

   use ISO_FORTRAN_ENV, only : wp => real64, out => output_unit, int => int64
   implicit none
  
   !> thresholds 
   real(wp), parameter :: thr = 1.0e-5_wp
   real(wp), parameter :: sparse_thr = 65.0_wp
   
   !> dimension
   integer, intent(in) :: ndim
   
   !> matrix
   real(wp), dimension(ndim,ndim), intent(in) :: A

   !> print level
   logical, intent(in) :: verbose

   !> if sparse
   logical, intent(out) :: sparse

   !> sparsity ratio
   real(wp) :: ratio
   
   !> temp
   integer(int) :: i, nonzero, full, dense

   !> max number of nonzero elements in a column
   real(wp) :: max
   
   dense = 0.0_wp
   max   = 0
    
   ! check overall sparsity !
   ratio = count(abs(A(:,:)) < thr) 
   ratio = ratio / (ndim**2) * 100.0_wp
   sparse = ratio > sparse_thr

   if (verbose) &
      write(out,'(a, 1x, f9.2, a)') "Sparsity ratio:", ratio , "%"

   ! check sparsity of columns !
   do i=1,ndim
      nonzero = count(abs(A(:,i)) > thr)
      dense = dense + (nonzero * nonzero * nonzero) ! comp time for submatrices
      if (nonzero > max) max = nonzero ! swap if nonzero is larger
   enddo

   ratio = real(max)/real(ndim) * 100.0_wp
   sparse = ratio < sparse_thr
   if(verbose) &
      write(out,'(a, 1x, f9.2, a)') "(max nonzero column / dimension) ratio ", ratio, "%"
   
   !  check comp time !
   full = ndim * ndim * ndim ! comp time for full matrix
   ratio = real(dense)/real(full) 
   sparse = ratio < 1.0_wp
   
   if (verbose) then
      if (sparse) then
         write(out,'(a,f12.1)') "The approximate decrease of comp time: ", 1.0_wp/ratio
      else
         write(out,'(a,f12.1)') "The approximate increase of comp time: ", ratio
      endif
   endif


end subroutine chk_sparsity

subroutine print_summary(ndim, chempot, sign_A, P, H, identity, iterative, NaN_, conv)
   
   use gtb_la, only : la_symm
   use iso_fortran_env, only : wp => real64, out => output_unit
   implicit none

   !> dimension
   integer, intent(in) :: ndim

   !> chemical potential
   real(wp), intent(in) :: chempot
   
   !> sign function
   real(wp), intent(in) :: sign_A(ndim,ndim)

   !> density matrix
   real(wp), intent(in) :: P(ndim,ndim)

   !> hamiltonian
   real(wp), intent(in) :: H(ndim,ndim)

   !> identity matrix
   real(wp), intent(in) :: identity(ndim,ndim)

   !> band-structure energy
   logical, intent(in) :: NaN_

   !> if converged
   logical, intent(in) :: conv

   !> if iterative
   logical, intent(in) :: iterative

   !> band-structure energy
   real(wp) :: PH(ndim,ndim), band   
   
   !> number of electrons
   real(wp) :: nel

   !> temp
   integer :: i
   real(wp) :: tmp(ndim,ndim)

   !> format holder
   character(len=*), parameter :: fmt_conv = '(2x,a,f14.9,a,f14.9,a,f14.9)'
   character(len=*), parameter :: fmt_not_conv = '(2x,a,F14.9,a)'

   if (NaN_) then
      write(out,'(a,f9.4,a)') "NaN in the density matrix (", chempot, ")"
      return
   endif

   call la_symm(P,H,PH) 
   
   band=0.0_wp
   do i=1,ndim
      band=band+PH(i,i)
   enddo

   nel = 2.0_wp * check_el2(ndim,sign_A,identity)
   
   if(.not.iterative) then
      write(out,'(a,f9.4,a,f9.4,a,f9.4)') &
         "chempot = ", chempot, ", electrons = ",nel, ", band-str E = ", band*2
   else
      if (conv) then
         write(out,fmt_conv) "chempot = ", chempot, ", electrons = ", nel ,", Band-str E = ", band*2
      else
         write(out,fmt_not_conv) "chempot = ", chempot, ', not converged'
      endif
   endif

   write(out, *)

end subroutine print_summary

!> the chempot prediction algorithm
subroutine chempot_scan(ndim,eigvec,eigval2,nel,chempot)

   use iso_fortran_env, only : wp => real64
   implicit none

   !> dimensionality
   integer, intent(in) :: ndim
   
   !> eigenvalues of sign function
   real(wp), intent(in) :: eigval2(ndim,ndim)

   !> eigenvectors of sign function
   real(wp), intent(in) :: eigvec(ndim,ndim)

   !> exact  number of electrons
   real(wp), intent(in) :: nel

   !> chemical potential
   real(wp), intent(inout) :: chempot
   

   !> iterators
   integer :: i

   !> number of electrons
   real(wp) :: nel_guess

   !> coorection to chempot
   real(wp) :: chempot_corr
   
   do while( nel-nel_guess < 1.0E-6_wp) 
      
     ! nel_guess=0.0_wp
     ! do i=1,dim
     !    k 



   enddo
end subroutine chempot_scan

end module purify

