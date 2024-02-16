module timing_utilities

   use iso_fortran_env, only: wp => real64, int64
   implicit none
   private
   public :: start_timer, stop_timer, print_elapsed_time
   integer(int64) :: start_count, stop_count, count_rate
   real(wp) :: start_cpu, stop_cpu

contains
   
   !> start timer button
   subroutine start_timer()

      call system_clock(start_count, count_rate)
      call cpu_time(start_cpu)
   
   end subroutine start_timer

   !> stop timer button
   subroutine stop_timer(str)

      character(len=*), intent(in) :: str
      call system_clock(stop_count)
      call cpu_time(stop_cpu)
      call print_elapsed_time(str)
      
   end subroutine stop_timer

   !> print elapsed time 
   subroutine print_elapsed_time(str)
      
      use iso_fortran_env, only: out => output_unit

      character(len=*), intent(in) :: str
      real(wp) :: walltime, cputime   
      integer(int64) ::  cpudays, cpuhours, cpumins
      integer(int64) :: walldays,wallhours,wallmins
      
      walltime = real(stop_count - start_count, wp) / real(count_rate, wp)
      cputime = stop_cpu - start_cpu

      ! convert overall elapsed CPU time into days, hours, minutes !
      cpudays = int(cputime/86400._wp)
      cputime = cputime - cpudays*86400._wp
      cpuhours = int(cputime/3600._wp)
      cputime = cputime - cpuhours*3600._wp
      cpumins = int(cputime/60._wp)
      cputime = cputime - cpumins*60._wp

      ! convert overall elapsed wall time into days, hours, minutes !
      walldays = int(walltime/86400._wp)
      walltime = walltime - walldays*86400._wp
      wallhours = int(walltime/3600._wp)
      walltime = walltime - wallhours*3600._wp
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp
      
      write(out, '(a)') repeat('=', 80)
      write(out, '(a, a)') 'Elapse time for ', trim(str)
      write(out,'(1x,a26,i5," d, ",i2," h, ",i2," min, ",f15.8," sec")') &
         'Wall time: ',walldays,wallhours,wallmins,walltime
      write(out,'(1x,a26,i5," d, ",i2," h, ",i2," min, ",f15.8," sec")') &
         'CPU time: ', cpudays, cpuhours, cpumins, cputime
      write(out, '(a)') repeat('=', 80)
  
   end subroutine print_elapsed_time

end module timing_utilities
