module initialize_update
  !------------------------------------!
  !Use system parameters to initialize
  !the positions and update the positions.
  !------------------------------------!
  implicit none

  contains

subroutine Initialize_position
  !------------------------------------!
  !Initialize position
  !   This program is used to initialize the position of
  !   Polyelectrolytes and ions, and parameters, energy of
  !   the potential.
  !Input
  !   pos, random_or_uniform
  !Output
  !   pos
  !External Variables
  !   pos, random_or_uniform
  !Routine Referenced:
  !1.subroutine random_grafted
  !   initialize chains by randomly grafting on the plate
  !2.subroutine uniform_grafted
  !   initialize chains by uniformly grafting on the plate
  !3.subroutine initialize_ions
  !   initialize ions in the system
  !------------------------------------!
  use global_variables
  implicit none

  pos=0
  
  if ( random_or_uniform == 0 ) then
    call random_grafted     ! Don't forget periodic condition
  else 
    call uniform_grafted
  end if

  if ( qq /= 0 ) then
    call initialize_ions
  end if

end subroutine Initialize_position


subroutine Monte_Carlo_Move( EE, DeltaE )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out)   :: DeltaE
  integer :: j
  real*8 :: EE1, EE2

  do j = 1, NN-Ngl
!     call total_energy(EE1)

    call Choose_Particle
    call New_Position
    call Delta_Energy(DeltaE)
    call Move_or_not(EE, DeltaE)

    !
    !test EE2-EE1 = DeltaE
!     call total_energy(EE2)
!     write(*,*) EE2 - EE1, DeltaE, EE2, EE1  
  end do

end subroutine Monte_Carlo_Move


subroutine Monte_Carlo_Move_and_Time( EE, DeltaE, time )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(out) :: DeltaE
  real*8, intent(inout)  :: EE
  real*8, dimension(3), intent(out) :: time
  integer :: j
  
  time = 0
  do j = 1, NN-Ngl
    call Choose_Particle
    call New_Position
    call Delta_Energy_time(DeltaE,time)
    call Move_or_not(EE, DeltaE)
  end do

end subroutine Monte_Carlo_Move_and_Time


subroutine choose_particle
  !------------------------------------!
  !This subroutine is used to choose a particle ip to move.
  !   
  !Input
  !   
  !Output
  !   ip
  !External Variables
  !   NN, Nm, Npe, ip
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd                 

  call random_number(rnd)
  ip = int(rnd*NN) + 1
  !
  !The monomer anchored on the plate can't move, so we need to choose again.
  do while( mod(ip,Nml) == 1 .and. ip <= Npe )
    call random_number(rnd)
    ip = int(rnd*NN) + 1
  end do

end subroutine choose_particle


subroutine New_Position
  !--------------------------------------!
  !This program is used to generate new position.
  !   
  !Input
  !   ip
  !Output
  !   pos1
  !External Variables
  !   pos, pos1, dr
  !Routine Referenced:
  !1. Periodic_condition( rr(2) )
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd(3)

  pos_ip0 = pos(ip,:)
  call random_number(rnd)
  pos_ip1(1:3) = pos_ip0(1:3) + (rnd - 0.5D0) * dr
  pos_ip1(4)   = pos_ip0(4)
  call periodic_condition( pos_ip1(1:2) )

end subroutine New_Position


subroutine Move_or_not(EE, DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   pos, pos_ip0, pos_ip1, ip, Beta
  !Routine Referenced:
  !1.
  !Reference:
  !
  !--------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8,  intent(in)   :: DeltaE
  real*8,  intent(inout) :: EE
  real*8  :: rnd
  !
  !Judge whether move or not
  if ( DeltaE < 0 ) then
    pos(ip,1:3) = pos(ip,1:3) + pos_ip1(1:3) - pos_ip0(1:3)
    EE = EE + DeltaE
    if ( pos_ip0(4) /= 0 ) then
      call update_rhok
    end if
  else 
    call random_number(rnd)
    if ( rnd < Exp(-Beta*DeltaE) ) then
      pos(ip,1:3) = pos(ip,1:3) + pos_ip1(1:3) - pos_ip0(1:3)
      EE = EE + DeltaE
      if ( pos_ip0(4) /= 0 ) then
        call update_rhok
      end if
    endif
  endif
end subroutine Move_or_not

end module initialize_update



