module checks
   use iso_fortran_env, only : wp => real64, out => output_unit
   implicit none
   public :: check_density, check_electrons, bandStrEnergy, check_el2
   private
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
   character(len=*), parameter :: form = '( /, a, 3x, f12.6, /, a, 1x, f12.6)'
   
   ! print all input matrices !
   call print_packed_matrix(ndim,P,str)
   
   ! transform to matrices !
   call blowsym(ndim,P,Psym)
   call blowsym(ndim,S,Ssym)
   call blowsym(ndim,H,Hsym)

   write(6,form) "Number of electrons:", check_electrons(ndim,Ssym,Psym), "Band structure energy:", bandStrEnergy(ndim,Hsym,Psym) 
   write(6, "(a)") repeat('-',72)

end subroutine check_density

real(wp) function bandStrEnergy(ndim,H,P) result(band)
   
   use gtb_la, only : la_symm 
   use iso_fortran_env, only : wp => real64
   implicit none
   
   integer, intent(in) :: ndim
   real(wp), intent(in) :: H(ndim,ndim)
   real(wp), intent(in) :: P(ndim,ndim)
   
   real(wp) :: PH(ndim,ndim)
   integer :: i 
   
   band=0.0_wp

   call la_symm(P,H,PH)

   do i=1, ndim
      band = band + PH(i,i)
   enddo

end function bandStrEnergy

!> number of electorns from (18)
real(wp) function check_el2(ndim, sign_, I) result(nel)
   
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in) :: ndim
   real(wp), intent(in) :: sign_(ndim,ndim) 
   real(wp), intent(in) :: I(ndim,ndim)
  
   real(wp) :: tmp(ndim,ndim)
   integer :: ii, jj 
   
   nel=0.0_wp
   
   tmp = 0.5_wp * (I - sign_)
   do ii=1, ndim
     nel = nel + tmp(ii,ii)
   enddo

end function check_el2

real(wp) function check_electrons(ndim,S,P) result(check_nel)
   
   use gtb_la, only : la_symm 
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in) :: ndim
   real(wp), intent(in) :: S(ndim,ndim) 
   real(wp), intent(in) :: P(ndim,ndim)
   
   real(wp) ::  PS(ndim,ndim) 
   integer :: i 
   
   check_nel=0.0_wp
   call la_symm(P,S,PS) 
   
   do i=1, ndim
      check_nel=check_nel+PS(i,i)
   enddo

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
end module checks