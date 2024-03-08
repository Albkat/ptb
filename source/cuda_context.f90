module cuda_
   use accel_lib
   implicit none
   public :: ctx, initialize_ctx
   private
   type(cuda_context) :: ctx

contains
subroutine initialize_ctx()
   integer :: err

   err = load_cuda()
   if (err /= 0) then
      stop "CUDA acceleration library could not be loaded."
   end if
   
   ctx = cuda_init(err)
   if (err /= 0) then
      stop "Could not initialize GPU acceleration library."
   end if

end subroutine initialize_ctx
end module cuda_
