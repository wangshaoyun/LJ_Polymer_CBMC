program main
use global_variables
use input_output
use initialize_update
use compute_energy
implicit none

!##########Data Dictionary############!
  integer :: i, j, k
  real*8  :: EE, EE1=0, st, fn
  real*8  :: DeltaE, time(3)
!#####################################!

!#############Initialize##############!
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
  !
  !
  if (restart_or_continue /= 1 ) then
    !
    !initialize position
    call Initialize_position
    call write_pos
    call write_pos1(1)
    call write_hist
    !
    !initialize energy and parameters of potential
    call initialize_energy_parameters
    !
    !Error analysis of Ewald sum
    if ( qq /= 0 ) then
      call error_analysis
    end if
    !
    !Compute total energy
    call total_energy(EE)
    i=1
  else
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !initialize energy and parameters of potential
    call initialize_energy_parameters
    !
    !Error analysis of Ewald sum
    if ( qq /= 0 ) then
      call error_analysis
    end if
    !
    !Compute total energy
    call total_energy(EE)
  end if
!#####################################!

call Monte_Carlo_Move_and_Time(EE, DeltaE, time)
write(*,*) 'time in lj            :', time(1)
write(*,*) 'time in real space    :', time(2)
write(*,*) 'time in fourier space :', time(3)

!##############Preheation#############!
 if ( i <= StepNum0 ) then
  do step = i, StepNum0
    if ( mod(step,DeltaStep1) == 0 ) then
      call Monte_Carlo_Move_and_Time(EE, DeltaE, time)
      call compute_physical_quantities
      call total_energy(EE1)
      call write_physical_quantities( step, EE, EE1, DeltaE, time )
!       EE = EE1
      write(*,*) 'time in lj            :', time(1)
      write(*,*) 'time in real space    :', time(2)
      write(*,*) 'time in fourier space :', time(3)
    else
      call Monte_Carlo_Move(EE, DeltaE)
    end if    
    call update_verlet_list
    if ( mod(step,DeltaStep3) == 0 ) then
      call write_pos1(step)
      if ( qq /= 0 ) then
        call error_analysis
      end if
    end if
  end do
  i = step
end if
!#####################################!

!###############Running###############!
  do step=i, StepNum+StepNum0
    if ( mod(step,DeltaStep1) == 0 ) then 
      call Monte_Carlo_Move_and_Time(EE, DeltaE, time)
      call compute_physical_quantities
      call total_energy(EE1)
      call write_physical_quantities( step, EE, EE1, DeltaE, time )
    else
      call Monte_Carlo_Move(EE, DeltaE)
    end if
    call update_verlet_list
    if ( mod(step, DeltaStep2) == 0 ) then
      call histogram
    end if
    if ( mod(step, DeltaStep3) == 0 ) then
      call write_pos1(step)
      call write_hist
    end if
  end do
!#####################################!

  call cpu_time(finished)
  total_time=finished-started+total_time
  call write_time(total_time)
  write(*,*) 'finished!'

end program main








