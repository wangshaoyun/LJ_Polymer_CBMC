module input_output
  implicit none
  
  save

  real*8, allocatable, dimension(:,:)  :: hist1   !brushes distribution
  real*8, allocatable, dimension(:,:)  :: hist2   !height distribution
  real*8, allocatable, dimension(:,:)  :: hist3   !ions distribution
  real*8, allocatable, dimension(:,:)  :: hist4   !cos distribution
  real*8, allocatable, dimension(:,:)  :: hist5   !max height distribution
  real*8, allocatable, dimension(:,:)  :: hist6   !polar angle distribution
  integer, allocatable, dimension(:,:) :: hist7   !top view of monomer
  integer, allocatable, dimension(:,:) :: hist8   !side view of monomer
  integer, allocatable, dimension(:,:) :: hist9   !side view of counterion

  real*8 :: hh
  real*8 :: hh_max
  real*8 :: Rg
  real*8 :: Rgz
  real*8 :: RR2
  real*8 :: RR2z

  contains

subroutine initialize_parameters
  !------------------------------------!
  !Input and initialize system, timing and histogram parameters.
  !Allocate arrays.
  !------------------------------------!
  use global_variables
  implicit none
  logical alive

  !
  !Judge whether restart or continue
  Inquire(file='start_time.txt',exist=alive)
  if (alive) then
    open(11,file='./start_time.txt')
      read(11,*) restart_or_continue
    close(11)
  else
    restart_or_continue=0
  end if

  !
  !Input parameters
  call read_data

  ! data operation
  Nx  = nint( Ngl ** (1./2) )
  Ny  = Nx
  Npe = Nml * Ngl
  if ( qq == 0 ) then
    Nq = 0
    NN = Npe
  else
    if ( man == 1 ) then !the anchored points are not charged
      Nq = Ngl * (Nml - 1) * ( 1 + nint(abs(qq)) ) 
      NN = Npe + Ngl * (Nml-1) * nint(abs(qq))
    else   ! if man/=1, mod(Nml,man) must be zero
      Nq = Ngl * Nml / man * ( 1 + nint(abs(qq)) )
      NN = Npe + Ngl * Nml / man * nint(abs(qq))
    end if  
  end if

  Lx = sqrt(Ngl / sigmag * ratio_xy)
  Ly = Lx / ratio_xy
  Z_empty  = 1 + Z_empty*Lx/Lz

  !Write data
  call write_data_to_screen

  !
  !Allocate arrays and initialize them
  call allocatte_arrays_and_initialize

end subroutine initialize_parameters


subroutine read_data
  implicit none
  use global_variables

  open(unit=100, file='system_data.txt')
    read(100,*) random_or_uniform
    read(100,*) Lz
    read(100,*) ratio_xy
    read(100,*) sigmag
    read(100,*) Beta
    read(100,*) qq
    read(100,*) Z_empty
    read(100,*) Nml
    read(100,*) Ngl
    read(100,*) man
    read(100,*) StepNum0
    read(100,*) StepNum
    read(100,*) DeltaStep1
    read(100,*) DeltaStep2
    read(100,*) DeltaStep3
    read(100,*) multistep
    read(100,*) dr
    read(100,*) SizeHist
  close(100)
end subroutine read_data


subroutine write_data_to_screen
  use global_variables
  implicit none

  write(*,*) '******************system_data***********************'
  write(*,*) 'Total anchored particles, Ngl:', Ngl
  write(*,*) 'Particles of each chain,  Nml:', Nml
  write(*,*) 'Total particles,          NN :', NN
  write(*,*) 'Total charged particles,  Nq :', Nq
  write(*,*) 'Total brushes particles,  Npe:', Npe
  write(*,*) 'Manning effect parameter, man:', man
  write(*,*) 'Length of the box,        Lx :', Lx
  write(*,*) 'Width of the box,         Ly :', Ly
  write(*,*) 'Height of the box,        Lz :', Lz
  write(*,*) '****************************************************'

  write(*,*) '******************running_steps*********************'
  write(*,*) 'Preheating steps             :', StepNum0
  write(*,*) 'Running steps                :', StepNum
  write(*,*) 'Total steps                  :', (StepNum0+StepNum)
  write(*,*) 'DeltaStep1                   :', DeltaStep1
  write(*,*) 'DeltaStep2                   :', DeltaStep2
  write(*,*) 'DeltaStep3                   :', DeltaStep3
  write(*,*) 'Multisteps of coulomb energy :', MultiStep
  write(*,*) 'Distance of each move        :', dr
  write(*,*) '****************************************************'
end subroutine write_data_to_screen


subroutine allocatte_arrays_and_initialize
  use global_variables
  implicit none

  allocate( pos(NN, 4)                )
  allocate( hist1(SizeHist, 2)        )
  allocate( hist2(SizeHist, 2)        )
  allocate( hist3(SizeHist, 2)        )
  allocate( hist4(Nml, 2)             )
  allocate( hist5(SizeHist, 2)        )
  allocate( hist6(SizeHist, 2)        )
  allocate( hist7(SizeHist, SizeHist) )
  allocate( hist8(SizeHist, SizeHist) )
  allocate( hist9(SizeHist, SizeHist) )

  hist1=0
  hist2=0
  hist3=0
  hist4=0
  hist5=0
  hist6=0
  hist7=0
  hist8=0
  hist9=0
end subroutine allocatte_arrays_and_initialize


subroutine continue_read_data(l)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none
  integer, intent(out) :: l
  integer :: i, j 

  open(20,file='./data/pos1.txt')
    read(20,*) ((pos(i,j),j=1,4),i=1,NN)
  close(20)
  open(19,file='./start_time.txt')
    read(19,*)
    read(19,*) l
    read(19,*) total_time
  close(19)
  
  open(22,file='./data/hist1.txt')
  open(23,file='./data/hist2.txt')
  open(24,file='./data/hist3.txt')
  open(25,file='./data/hist4.txt')
  open(26,file='./data/hist5.txt')
  open(27,file='./data/hist6.txt')
  open(28,file='./data/hist7.txt')
  open(29,file='./data/hist8.txt')
  open(30,file='./data/hist9.txt')
    read(22,*) ((hist1(i,j),j=1,2),i=1,SizeHist)
    read(23,*) ((hist2(i,j),j=1,2),i=1,SizeHist)
    read(24,*) ((hist3(i,j),j=1,2),i=1,SizeHist)
    read(25,*) ((hist4(i,j),j=1,2),i=1,Nml)
    read(26,*) ((hist5(i,j),j=1,2),i=1,SizeHist)
    read(27,*) ((hist6(i,j),j=1,2),i=1,SizeHist)
    read(28,*) ((hist7(i,j),j=1,SizeHist),i=1,SizeHist)
    read(29,*) ((hist8(i,j),j=1,SizeHist),i=1,SizeHist)
    read(30,*) ((hist9(i,j),j=1,SizeHist),i=1,SizeHist)
  close(30)
  close(29)
  close(28)
  close(27)
  close(26)
  close(25)
  close(24)
  close(23)
  close(22)
end subroutine continue_read_data


subroutine compute_physical_quantities
  !----------------------------------------!
  !
  !input:
  !  pos
  !output:
  !  Rg, Rgz, RR2, RR2z, hh, hh_max
  !External Variables:
  !  Ngl, Nml, Npe, NN,
  !----------------------------------------!
  use global_variables
  implicit none
  integer i,j,k
  real*8 :: rr,h_max
  real*8, dimension(3) :: rij
  
  Rg     = 0
  Rgz    = 0
  RR2    = 0
  RR2z   = 0
  hh     = 0
  hh_max = 0
  do i = 1, Ngl
    h_max = 0
    do j = 1, Nml
      do k = 1, Nml
        if ( j == k ) cycle
        call rij_and_rr(rij, rr, j, k)
        Rg  = Rg + rr
        Rgz = Rgz + rij(3)*rij(3)
      end do
      hh = hh + pos(Nml*(i-1)+j,3)
      if ( h_max < pos((i-1)*Nml+j,3) ) then
        h_max = pos((i-1)*Nml+j,3)
      end if
    end do
    hh_max = hh_max + h_max
    call rij_and_rr(rij, rr, Nml*(i-1)+1, Nml*i)
    RR2  = RR2 + rr
    RR2z = RR2z + rij(3)*rij(3)
  end do
  hh_max = hh_max / Ngl
  hh     = hh / Npe
  RR2    = RR2 / Ngl
  RR2z   = RR2z / Ngl
  Rg     = Rg / Nml / (Nml-1) / Ngl / 2
  Rgz    = Rgz / Nml / (Nml-1) / Ngl / 2
  
end subroutine compute_physical_quantities


subroutine histogram
  !----------------------------------------!
  !input:
  !  pos
  !output:
  !  hist1(distribution hisotgram from PE to rod)
  !External Variants:
  !  Npe, Ngl, NN, Nml, Sizehist, Lz 
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k
  real*8, dimension(3) :: rij
  real*8 rsqr, max_h, theta, rr
  
  !
  !hist1
  do i=1, Npe
    if ( mod(i,Nml)==1 .and. i<=Npe ) cycle
    k = ceiling( pos(i,3) / (Lz/SizeHist) )
    if ( k<=0 .or. k>SizeHist ) then
      write(*,*) 'Wrong in histogram1'
      cycle
    end if
    hist1(k,2) = hist1(k,2) + 1
  end do
  !
  !hist2
  do i = 1, Ngl
    k = ceiling( pos(i*Nml,3) / (Lz/SizeHist) )
    if ( k<=0 .or. k>SizeHist ) then
      write(*,*) 'Wrong in histogram2'
      cycle
    end if
    hist2(k,2) = hist2(k,2) + 1
  end do
  !
  !hist3
  do i = 1, Ngl*nint(abs(qq))
    k = ceiling( pos(Ngl*Nml+i,3) / (Lz/SizeHist) )
    if ( k<=0 .or. k>SizeHist ) then
      write(*,*) 'Wrong in histogram3'
      cycle
    end if
    hist3(k,2) = hist3(k,2) + 1
  end do
  !
  !hist4
  do i=1, Ngl
    do j=2, Nml
      call rij_and_rr(rij, rsqr, (i-1)*Nml+j, (i-1)*Nml+j-1)
      hist4(j,2) = hist4(j,2) + rij(3)/sqrt(rsqr)/Ngl
    end do
  end do
  !
  !hist5
  do i = 1, Ngl
    max_h = 0
    do j = 1, Nml
      if ( max_h < pos((i-1)*Nml+j,3) ) then
        max_h = pos((i-1)*Nml+j,3)
      end if
    end do
    k = ceiling( max_h / (Lz/SizeHist) )
    if ( k<=0 .or. k>SizeHist ) then
      write(*,*) 'Wrong in histogram5'
      cycle
    end if
    hist5(k,2) = hist5(k,2) + 1
  end do
  !
  !hist6
  do i = 1, Ngl
    call rij_and_rr(rij,rr,i*Nml,(i-1)*Nml+1)
    theta = acos( pos(i*Nml,3) / sqrt(rr) )
    k = ceiling( theta / (pi/2/SizeHist) )
    if ( k<=0 .or. k>SizeHist ) then
      write(*,*) 'Wrong in histogram6'
      cycle
    end if
    hist6(k,2) = hist6(k,2) + 1
  end do
  !
  !hist7
  do k = 1, Npe
    if ( mod(k,Nml)==1 .and. k<=Npe ) cycle
    i = ceiling( (pos(k,1)+Lx/2) / (Lx/SizeHist) )
    j = ceiling( (pos(k,2)+Ly/2) / (Ly/SizeHist) )
    if (i<=0 .or. i>SizeHist .or. j<=0 .or. j>SizeHist) then
      write(*,*) 'Wrong in histogram7'
      cycle
    end if
    hist7(i,j)=hist7(i,j)+1
  end do
  !
  !hist8
  do k = 1, Npe
    if ( mod(k,Nml)==1 .and. k<=Npe ) cycle
    i = ceiling( (pos(k,2)+Ly/2) / (Ly/SizeHist) )
    j = ceiling( pos(k,3) / (Lz/SizeHist) )
    if ( i<=0 .or. i>SizeHist .or. j<=0 .or. j>SizeHist ) then
      write(*,*) 'Wrong in histogram8'
      write(*,*) i, j, k, pos(k,1:3)
      cycle
    end if
    hist8(i,j) = hist8(i,j) + 1
  end do
  !
  !hist9
  do k = Npe+1, NN
    i = ceiling( (pos(k,2)+Ly/2) / (Ly/SizeHist) )
    j = ceiling( pos(k,3) / (Lz/SizeHist) )
    if ( i<=0 .or. i>SizeHist .or. j<=0 .or. j>SizeHist ) then
      write(*,*) 'Wrong in histogram9'
      cycle
    end if
    hist9(i,j) = hist9(i,j) + 1
  end do

end subroutine histogram


subroutine write_pos
  !----------------------------------------!
  !write position to pos.txt
  !input:
  !  pos
  !External Variants:
  !  NN
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: i

  open(30,file='./data/pos.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3), pos(i,4)
      300 format(4F15.6)
    end do
  close(30)

end subroutine write_pos


subroutine write_pos1(j)
  !----------------------------------------!
  !write position to pos1.txt
  !input:
  !  pos
  !External Variants:
  !  NN
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  integer :: i

  open(30,file='./data/pos1.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3), pos(i,4)
      300 format(4F15.6)
    end do
  close(30)

  open(32,file='./start_time.txt')
    write(32,*) 1
    write(32,*) j
    call cpu_time(finished)
    total_time=total_time+finished-started
    call cpu_time(started)
    write(32,*) total_time
  close(32)

end subroutine write_pos1


subroutine write_hist
  !----------------------------------------!
  !Write distribution histogram to the file hist1.txt ... hist4.txt
  !input:
  !  hist1, hist2, hist3, hist4,...hist9
  !External Variants: 
  !  SizeHist, Lz, 
  !----------------------------------------!
  use global_variables
  implicit none
  integer i,j
  
  open(34,file='./data/hist1.txt')
  open(35,file='./data/hist2.txt')
  open(36,file='./data/hist3.txt')
  open(37,file='./data/hist5.txt')
  open(38,file='./data/hist6.txt')
    do i=1,SizeHist
      hist1(i,1)=i*Lz/SizeHist
      write(34,340) hist1(i,1), hist1(i,2)
      hist2(i,1)=i*Lz/SizeHist
      write(35,340) hist2(i,1), hist2(i,2)
      hist3(i,1)=i*Lz/SizeHist
      write(36,340) hist3(i,1), hist3(i,2)
      hist5(i,1)=i*Lz/SizeHist
      write(37,340) hist5(i,1), hist5(i,2)
      hist6(i,1)=i*pi/2/SizeHist
      write(38,340) hist6(i,1), hist6(i,2)
    end do
    340 format(2F15.6)
  close(34)
  close(35)
  close(36)
  close(37)
  close(38)

  open(39,file='./data/hist4.txt')
    do i=1,Nml
      hist4(i,1)=i*1.
      write(39,340) hist4(i,1), hist4(i,2)
    end do
  close(39)

  open(40,file='./data/hist7.txt')
  open(41,file='./data/hist8.txt')
  open(42,file='./data/hist9.txt')
    do i=1,SizeHist
      write(40,'(500I10)') (hist7(i,j),j=1,SizeHist) 
      write(41,'(500I10)') (hist8(i,j),j=1,SizeHist)  
      write(42,'(500I10)') (hist9(i,j),j=1,SizeHist)  
    end do
  close(40)
  close(41)
  close(42)

end subroutine write_hist


subroutine write_physical_quantities(j, EE, EE1, DeltaE, time)
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  real*8, intent(in) :: EE
  real*8, intent(in) :: EE1
  real*8, intent(in) :: DeltaE
  real*8, dimension(3), intent(in) :: time
  real*8 :: prob

  if ( DeltaE < 0 ) then
    prob = 1
  else
    prob = exp(-Beta*DeltaE)
  end if 
  
  open(36,position='append', file='./data/height.txt')
    write(36,360) 1.*j, hh, hh_max, Rg, Rgz, RR2, RR2z
    360 format(7F15.6)
  close(36)

  open(37,position='append', file='./data/energy_and_time.txt')
    write(37,370) 1.*j, EE, EE1, DeltaE, prob, &
                  time(1), time(2), time(3)
    370 format(8F15.6)
  close(37)

end subroutine write_physical_quantities


subroutine write_time(time)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none

  real*8, intent(in) :: time
  open(10,file='./data/time.txt')
    write(10,*) 'time:(seconds)', real(total_time)
    write(10,*) 'time:(hours)  ', real(total_time/3600)
    write(10,*) 'time:(days)   ', real(total_time/86400)
    write(10,*) 'Lx:           ', real(Lx)
    write(10,*) 'Ly:           ', real(Ly)
    write(10,*) 'Lz:           ', real(Lz)
    write(10,*) 'Ngl:          ', Ngl
    write(10,*) 'Nml:          ', Nml
    write(10,*) 'Nq:           ', Nq
    write(10,*) 'NN:           ', NN
    write(10,*) 'sigmag:       ', real(sigmag)
    write(10,*) 'qq:           ', nint(qq)
  close(10)

end subroutine write_time


end module input_output
