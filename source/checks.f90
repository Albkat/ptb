module checks
   use iso_fortran_env, only : wp => real64, out => output_unit
   use accel_lib
   use cli
   use cuda_
   use gtb_la, only : la_gemm 
   implicit none
   public :: check_density, check_electrons, bandStrEnergy, check_el2, calculate_rmsd
   public :: trace
   private
   real(wp),save :: bndEnergy(2)
   real(wp),save :: electrons(2)
   logical, save :: run1 = .true.
contains

!> check density for number of electrons and band structure energy
subroutine check_density(ndim, P, S, H, str)

   use ISO_FORTRAN_ENV, only : wp => real64
   implicit none

   !> number of basis functions
   integer, intent(in) :: ndim
   
   !> density 
   real(wp), intent(in) :: P(ndim*(ndim+1)/2)

   !> overlap 
   real(wp), intent(in) :: S(ndim*(ndim+1)/2)

   !> hamiltonian
   real(wp), intent(in) :: H(ndim*(ndim+1)/2)

   !> message
   character(len=*), intent(in) :: str

   !> tmp matrices
   real(wp), dimension(ndim,ndim) :: Psym, Ssym, Hsym
   
   !> format for I/O
   character(len=*), parameter :: form = '( /, a, 3x, *, /, a, 1x, *)'

   real(wp) :: bnd, el

   ! print all input matrices !
   call print_packed_matrix(ndim,P,str)
   
   ! transform to matrices !
   call blowsym(ndim,P,Psym)
   call blowsym(ndim,S,Ssym)
   call blowsym(ndim,H,Hsym)

   el=check_electrons(ndim,Ssym,Psym)
   if(run1)then
      electrons(1) = el
   else
      electrons(2) = el
   endif
   bnd=bandStrEnergy(ndim,Hsym,Psym) 
   if(run1)then 
      bndEnergy(1) = bnd
   else
      bndEnergy(2) = bnd
   endif
   write(6,*) "Number of electrons: ",el 
   write(6,*)"Band structure energy: ",bnd 
   write(6, "(a)") repeat('-',72)
   run1 =.false.

end subroutine check_density

real(wp) function bandStrEnergy(ndim,P,H) result(band)
   
   implicit none
   
   integer, intent(in) :: ndim
   real(wp), intent(in) :: H(ndim,ndim)
   real(wp), intent(in) :: P(ndim,ndim)
   
   real(wp) :: PH(ndim,ndim)
   integer :: i, err 
   

   call la_gemm(P, H, PH)
   band = trace(PH)

end function bandStrEnergy

!> number of electrns from (18)
real(wp) function check_el2(ndim, sign_, I) result(nel)
   
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in) :: ndim
   real(wp), intent(in) :: sign_(ndim,ndim) 
   real(wp), intent(in) :: I(ndim,ndim)
  
   real(wp) :: tmp(ndim,ndim)
   integer :: ii, jj 
   
   
   tmp = 0.5_wp * (I - sign_)
   nel = trace(tmp)

end function check_el2

real(wp) function check_electrons(ndim, P, S) result(check_nel)
   
   implicit none
   integer, intent(in) :: ndim
   real(wp), intent(in) :: S(ndim,ndim) 
   real(wp), intent(in) :: P(ndim,ndim)
   
   real(wp) :: PS(ndim,ndim) 
   integer :: i 

   call la_gemm(P, S, PS) 
   check_nel = trace(PS)

end function check_electrons

!> idempotency check
subroutine idempotency(ndim, mat, S, msg)
   
   use gtb_la, only : la_symm
   use iso_fortran_env, only : wp => real64
   implicit none
   
   real(wp), intent(in) :: mat(ndim,ndim)
   real(wp), intent(in) :: S(ndim,ndim)
   integer, intent(in) :: ndim
   
   character(len=*), intent(in) :: msg
   real(wp) ::  matSmat(ndim,ndim) 
   
   matSmat=0.0_wp 
   
   ! P*S*P !
   call la_symm(mat,S,matSmat)
   call la_symm(matSmat,S,matSmat,alpha=0.5_wp)
   
   call print_blowed_matrix(ndim,matSmat,msg)

end subroutine idempotency


subroutine calculate_rmsd(ndim, A, B)
   use ISO_FORTRAN_ENV, only : wp => real64
   implicit none
   integer, intent(in) :: ndim
   real(wp), intent(in) :: A(ndim*(ndim+1)/2), B(ndim*(ndim+1)/2)
   real(wp) :: rmsd
   real(wp), dimension(ndim,ndim) :: P1, P2
   integer :: i, j
   real(wp) :: sum_sq_diff
   
   call blowsym(ndim,A,P1)
   call blowsym(ndim,B,P2)
   sum_sq_diff = 0.0_wp
   do i = 1, ndim
       do j = 1, ndim
           sum_sq_diff = sum_sq_diff + (P1(i, j) - P2(i, j))**2
       end do
   end do
   
   rmsd = sqrt(sum_sq_diff / (ndim**2))
   write(6,*)"RMSD: ", rmsd 
   write(6,*)"d_el: ", electrons(1)-electrons(2)
   write(6,*)"d_bnd:", bndEnergy(1)-bndEnergy(2)

end subroutine calculate_rmsd

pure function trace(matrix) result(tr)
   real(wp), intent(in) :: matrix(:,:)
   real(wp) :: tr
   integer :: i

   tr = 0.0_wp
   do i = 1, size(matrix,2)
      tr = tr + matrix(i,i)
   enddo
end function
end module checks