module purify 
   use timing_utilities
   use iso_fortran_env, only : wp => real64, out => output_unit
   use checks
   use cuda_
   use cli
   use accel_lib 
   use gtb_la, only : la_gemm
   use ieee_arithmetic, only : ieee_is_NaN
   implicit none
   public :: purification 
   private 
   real(wp), allocatable :: identity(:,:)
   real(wp), parameter :: conv_thrs = 1.0e-7_wp
contains

!> non-orthogonal purification
subroutine purification(H, ndim, S, P, P_purified)
   
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
   
   !> debug mode
   logical, parameter :: time(2) = [.false.,.true.]

   !> counter 
   integer :: i 
   
   !-------!
   ! SETUP !
   !-------!
   
   ! print the purification settings !
   call timer%click(2,'settings')
   call print_set(pur,out)
   call timer%click(2)
   
   ! preparation for purification !
   call timer%click(3,'transformation matrix')
   
   Hsym(:,:) = 0.0_wp
   Ssym(:,:) = 0.0_wp
   Psym(:,:) = 0.0_wp
    
   call blowsym(ndim,H,Hsym)  
   call blowsym(ndim,S,Ssym)
   allocate(identity(ndim,ndim))
   do i=1,ndim
      identity(i,i) = 1.0_wp
   enddo

   call timer_gpu%new(1,.false.)

   call get_transformation_matrix(ndim,Ssym,metric,pur%metric)

   call timer%click(3)
!------------------------------------------!
!-------------- Purification --------------!
!------------------------------------------! 
   call timer%click(4,'only purification')
   call iterator(ndim, metric, Ssym, Hsym, tmp)
   call packsym(ndim, tmp, P_purified)
   call timer%click(4)
   call timer%write(out)
   

end subroutine purification

!> global iterator  
subroutine iterator(ndim, X, S, H, P)

   real(wp), parameter :: thrs = 1.0e-8_wp

   !> dimensionality (number of BFs)
   integer, intent(in) :: ndim

   !> transformation matrix
   real(wp), intent(in) :: X(ndim,ndim)

   !> overlap matrix
   real(wp), intent(in) :: S(ndim,ndim) 

   !> hamiltonian
   real(wp), intent(in) :: H(ndim,ndim)

   !> purified density 
   real(wp), intent(out) :: P(ndim,ndim)
   
   
   !> chemical potential
   real(wp) :: chempot

   !> use bisection
   logical :: bisection

   !> lower/upper search bounds
   real(wp) :: lower_bound, upper_bound
   logical :: is_lower_bound, is_upper_bound

   !> increment
   real(wp) :: incr

   !> actual number of electrons
   real(wp) :: nel

   !> calculated number of electrons
   real(wp) :: nelCalc 

   !> counters
   integer :: i
   
   !> upper and lower spectrum bounds of H 
   real(wp) :: hmax, hmin 


   ! DEFS !

   P = 0.0_wp
   incr = pur%step
   nel = pur%nel
   chempot = pur%chempot
   bisection = .false.
   
   ! max and min eigenvalues !
   if (pur%method .eq. mcweeny .or. pur%method .eq. mcweeny_cp) then
      hmax  = 0.0_wp
      hmin  = 0.0_wp
      call eigs(ndim, H, hmax, hmin)
   endif

   ! Adjustable chempot !
   do i = 1, 100

      
      if (bisection) then
         chempot = (lower_bound + upper_bound) / 2.0_wp
         if (abs(upper_bound-lower_bound) < thrs) exit
      endif
      select case(pur%method)
      case(sign, sign_mm)
         
         call sign_purification(ndim, X ,S, H, chempot, &
            & P, pur%verbose, pur%sparse, nelCalc)
      
      case(mcweeny) 
         
         call gcp_purification(ndim, X, S, H, hmax, hmin, chempot, &
            & P, pur%cycle, pur%verbose, nelCalc) 
         
      endselect

      
      ! exit condition !
      if (abs(nelCalc-nel) < 0.001_wp .or. pur%fixed) exit
      
      ! adjust search boundaries !
      if (nelCalc < nel) then
         is_lower_bound = .true.
         lower_bound = chempot
         chempot = chempot + incr
      else                                 
         is_upper_bound = .true.
         upper_bound = chempot
         chempot = chempot - incr
      endif
      
      bisection = is_lower_bound .and. is_upper_bound
      incr = incr * 2.0_wp

   enddo

   P = P * 2.0_wp

endsubroutine iterator

!> Sign Method 
subroutine sign_purification(ndim, X, S, H, chempot, &
      & P_purified, verbose, sparse_cli, nelCalc)
   
   use ieee_arithmetic, only : ieee_is_NaN
   use gtb_lapack_eig, only : la_sygvd
   use gtb_la, only : la_gemm, la_symm

   !> number of basis functions
   integer, intent(in) :: ndim

   !> metric used for Hamiltonian diagonalization
   real(wp), intent(in) :: X(ndim,ndim)

   !> overlap matrix 
   real(wp), intent(in) :: S(ndim,ndim)
   
   !> hamiltonian matrix 
   real(wp), intent(in) :: H(ndim,ndim)
   
   !> chemical potential
   real(wp), intent(in) :: chempot
   
   !> purified density matrix 
   real(wp), intent(out) :: P_purified(ndim,ndim)
   
   !> verbosity
   logical, intent(in) :: verbose

   !> forced sparse mode
   logical, intent(in) :: sparse_cli

   !> calculated number of electrons
   real(wp), intent(out) :: nelCalc


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
   
   !> if sparse
   logical :: sparse

   !> if debug mode
   logical :: debug = .false.

   !> err
   integer :: err

   !> track diagonalization
   logical, save :: first = .true.
   logical :: make_a_guess
   logical :: conv

   ! calculate sign only if it is first iteration !
   ! or MM purification ! 
   make_a_guess = first .or. pur%method .eq. sign_mm

   if (make_a_guess) then

      ! sign(A) ! 
      term1 = chempot * identity
      call la_gemm(H, X, term2) !  H * metric 
      call la_gemm(X, term2, term3) ! metric * H * metric 
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

      Ascale = min(frob,gersh)
      sign_A = sign_A/Ascale

      ! for debugging !
      if (debug) &
         & call print_blowed_matrix(ndim, sign_A, 'sign_A')

      ! check sparisty!
      if (.not. allocated(pur%sparse)) then
         allocate(pur%sparse)
         call chk_sparsity(ndim, sign_A, verbose, pur%sparse) 
      endif
   
   endif


   ! Submatrix method !
   if (pur%sparse) then
      
      tmp = sign_A
      call apply_submatrix_method(ndim, tmp, sign_A, verbose)
 
      ! for debugging !
      if (debug) &
         & call print_blowed_matrix(ndim, sign_A, 'after submatrix method')
   
   ! Dense matrix case !
   else

      if (pur%method == sign_mm) then

         call sign_iterative(ndim, chempot, sign_A, H, X, pur%cycle, final_)  

      else
         call sign_diagonalization(ndim, sign_A, chempot, nelCalc, first, conv)
         first = .false.               
         if (.not.conv) return
      endif
   
   endif
      
   ! construct density matrix !
   call sign_density(ndim, P_purified, sign_A, X)
   
   ! print iterarion summary !
   call print_summary(ndim, chempot, P_purified, H, nelCalc, sign_A = sign_a)


end subroutine sign_purification

!> iterative purification of sign matrix
subroutine sign_iterative(ndim, chempot, sign_A, H, metric, cycles, final_)


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

   !> number of cycles
   integer, intent(in) :: cycles

   !> density matrix
   real(wp), intent(out) :: final_(ndim,ndim)
   
   !> iterator
   integer :: purificator
   
   !> temporary matrices
   real(wp), dimension(ndim,ndim) :: x2k, x4k, tmp, tmp2
   real(wp) :: norm


   !> error handle
   integer :: err

   !> number of electrons
   real(wp) :: nelCalc

   !> debug
   logical :: debug = .false.

   ! Pade approximation !
   purify: do purificator=1, cycles
      
      norm=norm2(sign_A)
      
      call la_gemm(sign_A, sign_A, x2k) 
      call la_gemm(x2k, x2k, x4k)
      tmp = (15.0_wp * identity) - (10.0_wp * x2k) + (3.0_wp * x4k)
      call la_gemm(sign_A, tmp, tmp2, alpha=0.125_wp)
      
      sign_A = tmp2
      
      call sign_density(ndim, final_, sign_A, metric) ! create density matrix from sign matrix

      ! convergence check ! 
      if (abs(norm2(sign_A)-norm) < conv_thrs) then
         write(out,'(1x,I2,a,1x)', advance = 'no') purificator, " cycles || "
         exit purify
      endif

      ! printout for fix modus !
      if (pur%fixed .or. debug) then

         if (purificator == 1) &
            write(out,'(//,2x,a)') "_ITERATIVE PURIFICATION_" 
         
         write(out,'(a,I0,2x)', advance='no') 'cycle ', purificator
         call print_summary(ndim, chempot, final_, H, nelCalc, sign_A=sign_A)
         
      endif

   enddo purify 

end subroutine sign_iterative

!> construct density matrix from sign matrix
subroutine sign_density(ndim, P, sign_A, X)

   use ieee_arithmetic, only : ieee_is_NaN
   
   implicit none

   !> number of basis functions
   integer, intent(in) :: ndim

   !> density matrix
   real(wp), intent(out) :: P(ndim,ndim)

   !> sign matrix
   real(wp), intent(in) :: sign_A(ndim,ndim)

   !> metric
   real(wp), intent(in) :: X(ndim,ndim)

   !> temporary matrices
   real(wp) :: tmp(ndim,ndim)


   ! (16) density matrix !
   call la_gemm(identity-sign_A,X,tmp)
   call la_gemm(X,tmp,P,alpha=0.5_wp)
   
   if (any((ieee_is_NaN(P)))) & ! check for NaN
      error stop "Not a Number encountered"
      
end subroutine sign_density

!> purification of the sign matrix via diagonalization

subroutine sign_diagonalization(ndim, sign_A, chempot, nelCalc, first, conv)

   use iso_fortran_env, only : int => int64
   use gtb_lapack_eig, only : la_syevd
   implicit none

   !> number of basis functions
   integer, intent(in) :: ndim

   !> sign matrix
   real(wp), intent(inout) :: sign_A(ndim,ndim)
   
   !> chempot
   real(wp), intent(in) :: chempot

   !> number of electrons
   real(wp), intent(out) :: nelCalc

   !> if first diagonalization 
   logical, intent(in) :: first

   !> if converged
   logical, intent(out) :: conv

   !> eigenvectors
   real(wp), allocatable, save :: eigvec(:,:)

   !> eigenvalues
   real(wp), allocatable, save :: eigval(:)
   
   real(wp) :: eigval_guess(ndim), eigval2(ndim,ndim)

   !> temporary vars
   integer :: i, j, info
   real(wp) :: tmp(ndim,ndim), tmp2

   logical :: debug = .false.

   if(first .or. debug) then 
      if (.not. allocated(eigval)) allocate(eigvec(ndim,ndim), eigval(ndim))
      eigvec = sign_A
      call la_syevd(eigvec, eigval, info)
   endif
   
   conv = .false.
   eigval2 = 0.0_wp
   nelCalc = 0.0_wp
   
   ! signum() !
   do i = 1, ndim

      tmp2 = 0.0_wp
      if (.not. debug) then
         eigval_guess(i) = eigval(i) - chempot
      else
         eigval_guess(i) = eigval(i)
      endif
      
      if (abs(eigval_guess(i)) < conv_thrs) then
         eigval2(i,i) = 0.0_wp
      else 
         if (eigval_guess(i) > 0.0_wp) then 
            eigval2(i,i) = 1.0_wp 
         else 
            eigval2(i,i) = -1.0_wp 
         endif
      endif

      if (.not. debug) then
         ! calculate nelCalc
         do j=1, ndim
            tmp2 = tmp2 + eigvec(j,i) * eigvec(j,i)
         enddo
         nelCalc = nelCalc + (1.0_wp - (tmp2 * eigval2(i,i)))
      endif

   enddo

   if (abs(pur%nel-nelCalc) < conv_thrs) conv = .true.

   if (conv .or.debug) then 
      tmp = 0.0_wp
      call la_gemm(eigval2, eigvec, tmp, transb='T')
      call la_gemm(eigvec, tmp, sign_A)
   endif

end subroutine sign_diagonalization

!> grand-canonical purification
subroutine gcp_purification(ndim, X, S, H, hmax, hmin, chempot, &
      & P_purified, cycles, verbose, nelCalc)
   
   !> number of basis functions
   integer, intent(in) :: ndim

   !> metric (S version)
   real(wp), intent(in) :: X(ndim,ndim)
   
   !> hamiltonian matrix
   real(wp), intent(in) :: H(ndim,ndim)
   
   !> overlap matrix
   real(wp), intent(in) :: S(ndim,ndim)
   
   !> upper and of the eigenvalue spectrum of H
   real(wp), intent(in) :: hmax, hmin 

   !> chemical potential
   real(wp), intent(inout) :: chempot
   
   !> purified density 
   real(wp), intent(out) :: P_purified(ndim,ndim)

   !> max purification number
   integer, intent(in) :: cycles

   !> printlevel
   logical, intent(in) :: verbose

   real(wp), intent(out) :: nelCalc
   
   !> iterators
   integer ::  purificator, i 
   
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
   
   !> debug mode
   logical, parameter :: debug = .false.
   
   
   ! construct initial guess !   
   ! (11a) !
   term1 = chempot * X
   call la_gemm(H, X, term2)
   call la_gemm(X, term2, term3)
   term4 = term1 - term3            
      
   ! frobenius norm !
   frob = sqrt(sum(term4**2))

   ! gershgorin norm !
   gersh = 0.0_wp
   do i = 1, ndim
      gersh = max(gersh,sum(abs(term4(:,i))))
   enddo
   
   ! print norms !
   if (verbose .and. debug) then
      write(out,'(2x,a,f14.9)'), "Gershgorin Norm  = ", gersh
      write(out,'(2x,a,f14.9)'), "Frobenius Norm   = ", frob
   endif

   alpha = 1.0_wp / min(gersh,frob)
   term5 = alpha * 0.5_wp * term4   ! scaling !
   P0 = term5 + 0.5_wp * X     ! shift !
   norm = norm2(P0)
   
   purify: do purificator = 1, cycles

      call la_gemm(S, P0, SP)
      call la_gemm(P0, SP, PSP)
      call la_gemm(PSP, SP, PSPSP)

      res = 3 * PSP - 2 * PSPSP

      ! if NaN !
      if (any(ieee_is_NaN(res))) &
         error stop "Not a Number encountered during Mcweeny Purification"  
         
      ! convergence check !
      if (abs(norm2(res)-norm) < conv_thrs) &
         exit purify
      
      if (debug .or. pur%fixed) then


         if (purificator == 1) &
            write(out,'(//,2x,a)') "_ITERATIVE PURIFICATION_" 
         
         write(out,'(1x,I0,a)', advance='no') purificator, ' cycle || '
         call print_summary(ndim, chempot, res, H, nelCalc, S=S)
      
      endif

      P0 = res
      norm = norm2(res)
   
   enddo purify
   
   call print_summary(ndim, chempot, res, H, nelCalc, S=S)
   P_purified = res
   
end subroutine gcp_purification


! subroutine cp_purification(ndim, metric, S, P_ptb, H, hmax, hmin, P_purified, cycle_, ver)
   
!    use ieee_arithmetic, only : ieee_is_NaN
!    use iso_fortran_env, only : wp => real64, out => output_unit
!    implicit none 
   
!    !> dimensionality (nbf)
!    integer, intent(in) :: ndim
   
!    !> metric
!    real(wp), intent(in) :: metric(ndim,ndim)

!    !> overlap matrix 
!    real(wp), intent(in) :: S(ndim,ndim)
   
!    !> PTB density matrix 
!    real(wp), intent(in) :: P_ptb(ndim,ndim)
   
!    !> hamiltonian 
!    real(wp), intent(in) :: H(ndim,ndim)
   
!    !> spectral bounds  of H
!    real(wp), intent(in) :: hmax, hmin 
   
!    !> purified density
!    real(wp), intent(out) :: P_purified(ndim,ndim)

!    !> number of cycles
!    integer, intent(in) :: cycle_
   
!    !> verbosity
!    logical, intent(in) :: ver

! !---------------------------------------------------------
   
!    !> iterators
!    integer :: num, purificator, i, j, iter
   
!    !> scaling parameter
!    real(wp) :: alpha_max, alpha_min, alpha
   
!    !> Gershgorin's circle theorem 
!    real(wp) :: h_max, h_min, buff_min, buff_max
   
!    !> l2 norm 
!    real(wp) :: norm 
   
!    !> trace
!    real(wp) :: tr, tr1, tr2
   
!    !> TODO
!    real(wp) :: hmin1, hmax1
   
!    !> purification factor 
!    real(wp) :: cn

!    !> chemical potential 
!    real(wp) :: chempot
   
!    !> norms 
!    real(wp) :: gersh, frob
   
!    !> desired number of electrons
!    real(wp) :: nel, N_e, N
   
!    !> equation terms
!    real(wp), dimension(ndim,ndim) :: term1, term2, term3, term4, term5
   
!    !> McWeeny purification terms
!    real(wp), dimension(ndim,ndim) :: SP, PSP, PSPSP, res 
   
!    !> temporary matrices
!    real(wp) :: matr1(ndim,ndim),matr2(ndim,ndim)
   
!    !> initial guess
!    real(wp) :: p0(ndim,ndim)
   
!    !> density in iteration 
!    real(wp) :: pn(ndim,ndim)
   
   
!    !> if NAN present
!    logical :: error
   
!    !> if purification converged
!    logical :: conv

!    !> if frobenius/gershgorin norm should be used
!    logical :: frob_norm
   
!    integer :: err
   
!    !-----------------------!
!    ! INITIAL CONFIGURATION !
!    !-----------------------!
   
!    frob_norm   = .true.
!    error       = .false.
   
!    N_e = 14.0_wp
!    N   = ndim

!    ! chempot !
!    tr  = 0.0_wp
!    h_min  = 0.0_wp
!    h_max  = 0.0_wp
!    do i=1,ndim
!       tr = tr + H(i,i)
!       buff_min = H(i,i)
!       buff_max = H(i,i)
!       do j=1,ndim
!          if (i.ne.j) then
!             buff_min = buff_min - abs(H(i,j))
!             buff_max = buff_max + abs(H(i,j))
!          endif
!       enddo
!       if (i.eq.1) then
!          h_min = buff_min
!          h_max = buff_max
!       else 
!          if (buff_max > h_max) h_max = buff_max
!          if (buff_min < h_min) h_min = buff_min
!       endif

!    enddo


!    chempot = tr / N

!    if (ver) then 
!       write(out,"(a, 2x, F14.6)") "Chemical potential", chempot   
!       write(out,"(a, 2x, F14.6)") "hmax ", hmax
!       write(out,"(a, 2x, F14.6)") "hmin ", hmin
!       write(out,"(a, 2x, F14.6)") "h_max ", h_max
!       write(out,"(a, 2x, F14.6)") "h_min ", h_min
!    endif

   
!    !-------------!  
!    ! CALCULATION !
!    !-------------!
!    ! initial guess !            
!    if (frob_norm) then
      
!       term1 = chempot * metric
      
!       if (pur%cuda) then  
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, metric, H, 0.0_wp,term2,err)
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, term2, metric, 0.0_wp,term2,err)
!       else
!          call la_gemm(metric,H,term2)
!          call la_gemm(term2,metric,term2)
!       endif

!       term3 = term1 - term2
      
!       ! frobenius norm !
!       frob=sqrt(sum(term3**2))

!       ! gershgorin norm !
!       gersh=0.0_wp
!       do i=1,ndim
!          gersh=max(gersh,sum(abs(term3(:,i))))
!       enddo
      
!       ! print norms !
!       if (ver) then
!          write(out, '(2x,a,f14.9)'), "Gershgorin Norm: ", gersh
!          write(out, '(2x,a,f14.9)'), "Frobenius Norm:  ", frob
!       endif
      
!       alpha=N_e/min(gersh,frob)

!       term4 =  metric * N_e    
!       p0 = ( (term3 * alpha)  + term4 ) / N_e

!    else
   
!       hmin1 = ( N - N_e) / (chempot - hmin)
!       hmax1 =  N_e / (hmax - chempot)
      
!       alpha = min(hmin1,hmax1)
      
!       term1 = alpha * chempot * metric
!       if(pur%cuda) then  
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, metric,H, 0.0_wp,term2,err)
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, alpha, term2,metric, 0.0_wp,term2,err)
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, N_e, identity,metric, 0.0_wp,term3,err)
!       else
!          call la_gemm(metric,H,term2)
!          call la_gemm(term2,metric,term2,alpha=alpha)
!          call la_gemm(identity,metric,term3,alpha=N_e)
!       endif
!       term4 = (term1 - term2 + term3) / N_e ! 1/N 
!       p0    = term4 
      
!    endif
   
!    if (ver) then
!       write(out,'(2x,a,f14.9)') , "number of electrons: ", check_electrons(ndim,p0,S)
!       call print_blowed_matrix(ndim,p0, "initial guess")
!    endif
!    pn = p0*2 

!    ! purification loop !
!    do num=1, cycle_ 

!       if(pur%cuda) then  
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, S, pn, 0.0_wp,SP,err)
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, pn,SP, 0.0_wp,PSP,err)
!          call cuda_dgemm(ctx, 'N', 'N', ndim, ndim, ndim, 1.0_wp, PSP,SP, 0.0_wp,PSPSP,err)
!       else
!          call la_gemm(S,pn,SP)
!          call la_gemm(pn,SP,PSP)
!          call la_gemm(PSP,SP,PSPSP)
!       endif
!       matr1=PSP-PSPSP
!       matr2=pn-PSP

!       tr1=0.0_wp
!       tr2=0.0_wp

!       do i=1, ndim
!          tr1=tr1+matr1(i,i)
!          tr2=tr2+matr2(i,i)
!       enddo
   
!       cn=tr1/tr2

      
!       if(cn >= 0.5_wp) then
!          res=( (1.0_wp - 2.0_wp * cn ) * pn + (1.0_wp + cn) * PSP - PSPSP ) / (1.0_wp-cn)
!       endif
      
!       nel = check_electrons(ndim, res, S)
!       if (num==1) write(*,'(a)') "Purification loop"
!       write(out,114) "iterations = ", num, ", electrons = ", nel, ", band-str E = ", bandStrEnergy(ndim,H,res),", norm2=", norm2(res)
!       114 format(a, I0, a, F14.6, a, F14.6, a, F14.6)
      
!       if (ver) then   
!          call print_blowed_matrix(ndim,res,"P")
!          write(out,'(a,2x,f14.6)') "The cn factor", cn
!       endif
      
!       pn=res
   
!    enddo
    
!    P_purified = pn

! end subroutine cp_purification

!> get eigenvalues and 
subroutine eigs(ndim,matrix,maxim,minim)
   
   use gtb_lapack_eig, only : la_syevd
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
      call la_syevd(matrix_syev,w,info)
      
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
subroutine get_transformation_matrix(ndim,matrix,root,metric)
   
   use gtb_la, only : la_gemm
   use gtb_lapack_eig, only : la_syevd
   use iso_fortran_env, only : wp => real64
   use accel_lib
   use cli, only : inv, inv_sqrt, sqrt_
   implicit none 
   integer, intent(in) :: ndim
   real(wp), intent(in) :: matrix(ndim,ndim)
   real(wp), intent(inout) :: root(ndim,ndim)

   !> S power
   integer, intent(in) :: metric
   
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

   !> Debug mode
   logical, parameter :: debug = .false.

   
   matrix_syev = matrix

   
   ! Eigendecomposition !
   call la_syevd(matrix_syev,w,info)

   D = 0.0_wp
   sqrtW=sqrt(w)

   selectcase(metric)
   case(sqrt_)
      do i = 1, ndim
         D(i,i) = sqrtW(i)
      enddo
   case(inv_sqrt)
      do i = 1, ndim
         D(i,i) = 1.0_wp/sqrtW(i)
      enddo
   case(inv)
      do i = 1, ndim
         D(i,i) = 1.0_wp/w(i)
      enddo
   endselect
   
   ! r = V * D^(1/2) * V^(-1)
   call la_gemm(matrix_syev,D,tmp)
   call la_gemm(tmp, matrix_syev,root,transb='T')
   
   ! check root
   if(debug) then
      selectcase(metric)
      case(inv)
         call print_blowed_matrix(ndim,root,'inverse')
         call la_gemm(root, matrix, tmp)
         call print_blowed_matrix(ndim,tmp,'identity')
      case(sqrt_)
         call la_gemm(root,root,tmp)
         call print_blowed_matrix(ndim,matrix,'original')
         call print_blowed_matrix(ndim,tmp,'check')
      case(inv_sqrt) 
         call la_gemm(root,root,tmp)
         call la_gemm(tmp,matrix,root2)
         call print_blowed_matrix(ndim,root2,'identity')
      endselect
   endif

   block
      logical :: apprx = .false.
      real(wp), dimension(ndim) :: a_1, a_inf
      real(wp), dimension(ndim,ndim) :: inverse, check
      real(wp) :: spect_norm


      ! using numerical techniques !
      if (apprx) then

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

   end block
end subroutine get_transformation_matrix

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
   logical, save :: first = .false.
   logical :: conv
   real(wp) :: nelCalc
   real(wp):: chempot
   
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
      call sign_diagonalization(nonzero,submatrix,chempot,nelCalc,first, conv)
      
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

subroutine print_summary(ndim, chempot, P, H, nelCalc, S, sign_A)
   implicit none

   !> dimension
   integer, intent(in) :: ndim

   !> chemical potential
   real(wp), intent(in) :: chempot
   
   !> density matrix
   real(wp), intent(in) :: P(ndim,ndim)

   !> hamiltonian
   real(wp), intent(in) :: H(ndim,ndim)

   !> number of electrons
   real(wp),  intent(out) :: nelCalc
   !> overlap
   real(wp), optional, intent(in) :: S(ndim,ndim)
   
   !> sign function
   real(wp), optional, intent(in) :: sign_A(ndim,ndim)
   
   !> temporary variables 
   real(wp) :: band, nel
   integer :: i

   band = bandStrEnergy(ndim, P, H)
  
   ! GCP case !
   if (present(S)) then
      nel = 2.0_wp * check_electrons(ndim, P, S)
   else if (present(sign_A)) then
      nel = 2.0_wp * check_el2(ndim, sign_A, identity)
   endif

   write(out,'(a,f9.4,a,f9.4,a,f13.6)') &
      "chempot = ", chempot, ", electrons = ",nel, ", band-str E = ", band * 2.0_wp

   nelCalc = nel

end subroutine print_summary

end module purify

