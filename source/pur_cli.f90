
module cli
   use cuda_
   use iso_fortran_env, only : wp => real64
   public :: pur, get_purification_cli 
   private
   
   type :: pur_set
   ! member variables !

      !> S metric
      character(len=:), allocatable :: metric

      !> purification type
      character(len=:), allocatable :: type
      
      !> for sign: iterative or diagonalisation
      logical :: iter = .false.

      !> chemical potetial
      real(wp) :: chempot = 0.0_wp

      !> upper limit for chempot search
      real(wp) :: limit = 30.0_wp

      !> step
      real(wp) :: step = 0.1_wp

      !> printout level
      logical :: verbose = .false.

      !> if chemical potential fixed
      logical :: fixed = .false.
      
      !> if sparsity should be forced
      logical, allocatable :: sparse 

      !> purification iteration number
      integer :: cycle = 20

      !> cuda support
      logical :: cuda = .false.

   contains
   ! type-bound procedures !
      
   end type pur_set
      
   
   type(pur_set) :: pur

contains

subroutine get_purification_cli()

 !> number of CLI arguments
   integer :: numArgs

   !> raw value of flag
   character*80 :: arg

   !> raw value of flag argument
   character(len=:), allocatable :: line 

   !> length of flag argument
   integer :: leng

   !> variables for readline
   real(wp) :: floats(10)
   character(len=80) :: str(10)
   integer :: ns, nf 

   pur%metric = "-0.5"

   ! Parser !
   numArgs = command_argument_count()

   do i=2, numArgs
      call get_command_argument(i,arg)
      
      ! verbose printout !
      if(index(arg,'-V').ne.0) pur%verbose = .true. 
         
      ! chempot !
      if(index(arg,'--chempot').ne.0) then

         ! initial guess !
         call get_command_argument(i+1,length=leng) 
         allocate( character(len=leng) :: line)
         call get_command_argument(i+1,value=line)
         call readline(line,floats,str,ns,nf)
         deallocate(line)

         if (nf.ne.0) then
            
            pur%chempot=floats(1)
            if (pur%chempot>=pur%limit) pur%limit = pur%chempot + 30.0_wp 
            
            ! limit !
            call get_command_argument(i+2,length=leng) 
            allocate( character(len=leng) :: line)
            call get_command_argument(i+2,value=line)
            call readline(line,floats,str,ns,nf) 
            deallocate(line)

            if (nf.ne.0) then
               
               if(floats(1)>pur%chempot) pur%limit = floats(1)
               
               ! step !
               call get_command_argument(i+3,length=leng) 
               allocate( character(len=leng) :: line)
               call get_command_argument(i+3,value=line)
               call readline(line,floats,str,ns,nf)
               deallocate(line)

               if(nf.ne.0 .and. floats(1)<(pur%limit-pur%chempot).and.floats(1)>0.0_wp) then
                  pur%step=floats(1)
               endif

            endif 
         endif
      endif
      
      ! mode sign(default) | mcweeny | cp !
      if(index(arg,'--pur').ne.0) then
         call get_command_argument(i+1,length=leng) 
         allocate( character(len=leng) :: line)
         call get_command_argument(i+1,value=line)
         if (line.eq."sign" .or. line.eq."mcweeny" .or. line.eq."cp") then
            pur%type = line

            if (line.eq."sign") then
               deallocate(line)
               call get_command_argument(i+2,length=leng) 
               allocate( character(len=leng) :: line)
               call get_command_argument(i+2,value=line)
               if (line.eq."iter") then
                  pur%iter = .true.
               endif
            endif
         
         endif

         deallocate(line)
      endif
      
      if(index(arg,"--metric").ne.0) then 
         call get_command_argument(i+1,length=leng) 
         allocate( character(len=leng) :: line)
         call get_command_argument(i+1,value=line)
         if (line.eq."-1" .or. line.eq."-0.5" .or. line .eq. "+1") then
            pur%metric = line
         endif
         deallocate(line)
      endif

      if(index(arg,"--fix").ne.0) pur%fixed = .true.

      if(index(arg,"--cycle").ne.0) then
         call get_command_argument(i+1,length=leng)
         allocate( character(len=leng) :: line)
         call get_command_argument(i+1,value=line)
         call readline(line,floats,str,ns,nf)
         if (nf.ne.0 .and. floats(1)>0) then
            pur%cycle = int(floats(1))
         endif
         deallocate(line)
      endif
      
      if(index(arg,"--sparse").ne.0) then
         allocate(pur%sparse)
         call get_command_argument(i+1,length=leng)
         allocate( character(len=leng) :: line)
         call get_command_argument(i+1,value=line)
         
         if (line.eq."true") then
            pur%sparse = .true.
         elseif(line .eq. "false") then
            pur%sparse = .false.
         else 
            pur%sparse = .true.
         endif

         deallocate(line)
      endif

      if(index(arg,"--cuda").ne.0) then
         pur%cuda = .true.
         call initialize_ctx()
      endif

   enddo

end subroutine get_purification_cli
end module cli