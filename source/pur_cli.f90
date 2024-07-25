
module cli
   use cuda_
   use iso_fortran_env, only : wp => real64
   implicit none
   
   !> modes 
   integer, parameter :: sign = 0 ! def
   integer, parameter :: sign_mm = 1 
   integer, parameter :: mcweeny = 2
   integer, parameter :: mcweeny_cp = 3

   !>  metrics
   integer, parameter :: inv = 4
   integer, parameter :: inv_sqrt = 5 ! def 
   integer, parameter :: sqrt_ = 6


   
   
   type :: pur_set
   ! member variables !

      !> S metric
      integer :: metric = 5

      !> purification type
      integer :: method  = 0
      
      !> for sign: iterative or diagonalisation
      logical :: iter = .false.

      !> chemical potetial
      real(wp) :: chempot = -1.5_wp

      !> upper limit for chempot search
      real(wp) :: limit = 30.0_wp

      !> step
      real(wp) :: step = 0.1_wp

      !> printout level
      logical :: verbose = .true.

      !> if chemical potential fixed
      logical :: fixed = .false.
      
      !> if sparsity should be forced
      logical, allocatable :: sparse 

      !> purification iteration number
      integer :: cycle = 40

      !> cuda support
      logical :: cuda = .false.

      !> nel
      integer :: nel 
   contains
   ! type-bound procedures !
      procedure :: settings 
   end type pur_set
      
   
   type(pur_set), allocatable  :: pur

contains

subroutine settings(this, io)

   character(len=*), parameter :: eq = "="

   !> purifictaion holder
   class(pur_set), intent(inout) :: this

   !> I/O unit
   integer, intent(in) :: io

   !> number of CLI arguments
   integer :: numArgs

   !> raw value of flag
   character*80 :: arg

   !> raw value of flag argument
   character(len=80) :: line 
   character(len=:), allocatable :: arg1, arg3

   !> length of flag argument
   integer :: leng

   !> variables for readline
   real(wp) :: floats(10)
   character(len=80) :: str(10)
   integer :: ns, nf 
   
   integer :: ios

   !> debug mode
   logical, parameter :: debug = .true.

   this%chempot=-1.5_wp
   readsettings: do
      read(io,'(a)',iostat=ios) line
      if (ios .ne. 0) exit
      call readline(line, floats,str,ns,nf)
      arg1 = trim(str(1))
      if (arg1(1:1) .eq. '!') cycle

      ! do only if equal sign is present as a second argument !
      if (trim(str(2)) == eq) then
         select case(arg1)         

         case('mode', 'method')
            arg3 = trim(str(3))
            
            select case(arg3)  
            case('mcw') 
               this%method = mcweeny
            
            case('mcwcp')
               this%method = mcweeny_cp
               this%fixed = .true.
            
            case('sign')
               this%method = sign

            case('signmm')
               this%method = sign_mm
            endselect

         case('metric')
            arg3 = trim(str(3))
            
            select case(arg3)  
            case('inverse') 
               this%metric = inv
            
            case('inverse_sqrt')
               this%metric = inv_sqrt
            
            case('sqrt')
               this%metric = sqrt_
            endselect
                     
         case('cycle')
            if (nf > 0) this%cycle = int(floats(1))

         
         case('chempot')
            if (nf > 0) this%chempot = floats(1)

         case('step')
            if (nf > 0) this%step = floats(1)

         case('sparse')
            arg3 = trim(str(3))
            if (arg3 .eq. 'true') then
               this%sparse = .true.
            else if (arg3 .eq. 'false') then
               this%sparse = .false.
            endif
         endselect

      else

         select case(arg1)
         case('fix')
            this%fixed = .true.
         case('verbose')
            this%verbose = .true.
         case('cuda')
            this%cuda = .true.
            call initialize_ctx()
         endselect         
      endif

   enddo readsettings


end subroutine settings

subroutine print_set(pur, out)

   !> purification holder
   type(pur_set), intent(in) :: pur

   !> IO unit
   integer,  intent(in) :: out

   ! HEADER !
   write(out,'(/,a)') repeat('*',72)
   write(out,'(a,1x,a,1x,a)') repeat('*',29), "PURIFICATION", repeat('*',29)
   write(out,'(a,/)') repeat('*',72)

   write(out,'(2x,a)') "__settings__" 
   write(out,'(2x,a,5x)',advance='no') "Purification type:            "
   selectcase(pur%method)
   case(sign)
      write(out,'(a)') 'Sign Diagonalization'
   case(sign_mm)
      write(out,'(a)') 'Sign Matrix Multiplication'
   case(mcweeny)
      write(out,'(a)') 'McWeeny Grand Canonical Purification'
   case(mcweeny_cp)
      write(out,'(a)') 'McWeeny Canonical Purification'
   endselect
   
   write(out,'(2x,a,5x)',advance='no') "S power:                      "
   selectcase(pur%metric)
   case(inv)
      write(out,'(a)') '-1'
   case(sqrt_)
      write(out,'(a)') '0.5'
   case(inv_sqrt)
      write(out,'(a)') '-0.5'
   endselect
   if (pur%method.ne.mcweeny_cp) &
   write(out,'(2x,a,3x,f13.8)')        "Initial Chemical Potential:   ", pur%chempot
   if (pur%method.ne.sign) &
   write(out,'(2x,a,5x,I0)')           "Matrix M*M Cycle Number   :   ", pur%cycle
   if (pur%fixed) then
      write(out,'(2x,a,5x,L1)') "Fixed Mode:                   ", pur%fixed
   else
      write(out,'(2x,a,2x,f14.8)') "Initial Increment:           ", pur%step
   endif
   if (pur%verbose) &
      write(out,'(2x,a,5x,L1)') "Verbose Printlevel:           ", pur%verbose
   if (allocated(pur%sparse)) &
      write(out,'(2x,a,5x,L1)') "Sparse mode:                  ", pur%sparse 
   if (pur%cuda) &
      write(out,'(2x,a,5x,L1)') "CUDA Mode:                    ", pur%cuda

   write(out,'()')

end subroutine print_set

end module cli