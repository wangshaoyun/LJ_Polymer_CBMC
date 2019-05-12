module compute_energy
  !--------------------------------------!
  !Input:
  ! pos, pos_ip0, pos_ip1, ip
  ! and the system parameters
  !Output: 
  ! EE, DeltaE 
  !--------------------------------------!
  implicit none

  save

!############coefficient in potential function#############!
!
!lj potential
  real*8,  private :: epsilon     !Energy unit epsilon in lj potential
  real*8,  private :: sigma       !Distance sigma in lj potential
  real*8,  private :: rc_lj       !Cut off radius of LJ potential
  real*8,  private :: rv_lj       !Verlet list radius of LJ potential
  real*8,  private :: rsk_lj      !Skin between cut off sphere and verlet list 
                                  !sphere
  integer, private :: npair1      !number of pairs in the lj verlet sphere
!
!fene potential
  real*8,  private :: R0_2        !Max bond length R0 square in FENE potential
  real*8,  private :: kFENE       !FENE spring constant k
  integer, private :: N_bond      !Number of Chemical bond of polymers
!##########end coefficient in potential function###########!


!##########################arrays##########################!
  integer, allocatable, dimension( : ), private :: lj_pair_list  
                                  !LJ potential verlet list
  integer, allocatable, dimension( : ), private :: lj_point
                                  !the particles near i are from
                                  !lj_pair_list(lj_point(i-1)) to 
                                  !lj_pair_list(lj_point(i))
  integer, allocatable, dimension( : ), private :: fene_list
                                  !list of chemical bonds
  integer, allocatable, dimension( : ), private :: fene_point
                                  !same as lj_point and real_point
!########################end arrays########################!


contains


subroutine initialize_energy_parameters
  !--------------------------------------!
  !Initial parameters are not inputted from file and compute
  !the total energy of the system.
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !Reference:
  !The computation of alpha, rc_real et al are refered to
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.304-306.
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read energy parameters from file
  call read_energy_parameters
  !
  !Initialize lj parameters and array allocate.
  !
  call initialize_lj_parameters
  !
  !build lj_pair_list and lj_point
  call build_lj_verlet_list
  !
  !Initialize fene parameters and array allocate.
  call build_fene_list
  !
  !write energy parameters
  call write_energy_parameters

end subroutine initialize_energy_parameters


subroutine total_energy (EE)
  !--------------------------------------!
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
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: EE

  EE=0

  call LJ_Energy(EE)

  call FENE_Energy(EE)

end subroutine total_energy


subroutine LJ_energy (EE)
  !--------------------------------------!
  !Compute total LJ potential energy,
  !including LJ energy of wall.
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   lj_point, lj_pair_list, pos, 
  !   epsilon, sigma, rc_lj, Lz
  !Routine Referenced:
  !1. rij_and_rr( rij, rr, i, j )
  !Reference:
  !1.In fact, the cut-off radius in good solvent is 2^(1/6), 
  !  which was first obatained by JOHN D. WEEKS, DAVID CHANDLER
  !  and HANS C. ANDERSEN. So it is called WCA potential.
  !  JOHN D. WEEKS, DAVID CHANDLER and HANS C. ANDERSEN, 'Role of 
  !  Repulsive Forces in Determining the Equilibrium Structure of
  !  Simple Liquids', THE JOURNAL OF CHEMICAL PHYSICS, Vol. 54, 
  !  pp.5237-5247, (1971).
  !2.The potential of particle and wall are 9-3 LJ potential which is
  !  cut off at 0.4^(1/6) = 0.86.
  !  Yu-Fan Ho, et al, 'Structure of Polyelectrolyte Brushes Subject
  !  to Normal Electric Fields', Langmuir, 19, pp.2359-2370, (2013).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j, k, l, m
  real*8  :: rr, rij(3), inv_rr2, inv_rr6

  do i = 1, NN
    if ( i == 1) then
      k = 1
      l = lj_point(1)
    else
      k = lj_point(i-1)+1
      l = lj_point(i)
    end if
    do m = k, l
      j = lj_pair_list(m)
      call rij_and_rr( rij, rr, i, j )
      if ( rr < rc_lj * rc_lj ) then
        inv_rr2  = sigma*sigma/rr
        inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
        EE = EE + 4 * epsilon * ( inv_rr6 * inv_rr6 - inv_rr6 + 0.25D0) / 2
      end if
    end do
  end do

  do i = 1, NN
    if ( pos(i,3) > 0 .and. pos(i,3) < 0.86 ) then
      rr = pos(i,3)
      EE = EE + epsilon * ( 2.D0/15 * ( sigma / rr )**9 &
                            - ( sigma / rr )**3 )
    elseif ( pos(i,3) <Lz .and. pos(i,3) > (Lz - 0.86) ) then
      rr = Lz-pos(i,3)
      EE = EE + epsilon * ( 2.D0/15 * ( sigma / rr )**9 &
                            - ( sigma / rr )**3 )
    endif
  end do
end subroutine LJ_energy


subroutine FENE_energy(EE)
  !--------------------------------------!
  !Compute FENE potential energy.
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   fene_list, fene_point
  !   kFENE, Nbond, R0_2
  !Routine Referenced:
  !1. rij_and_rr
  !Reference:
  !1.Michael Murat, Gary S. Grest, 'Structure of a Grafted
  !Polymer Brush: A Molecular Dynamics Simulation', Macromolecules,
  !vol. 22, pp.4054-4059, (1989).
  !The FENE potential was first derived on 1972.
  !2.Harold R. Warner, 'Kinetic Theory and Rheology of Dilute 
  !Suspensions of Finitely Extendible Dumbbells', Ind. Eng. Chem.
  !Fundam., Vol. 11, No. 3, pp.379-387, (1972).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j, k, l, m
  real*8  :: rij(3), rr

  do i = 1, Npe
    if ( i == 1 ) then
      k = 1
      l = fene_point(1)
    else
      k = fene_point(i-1) + 1
      l = fene_point(i)
    end if
    do m = k, l
      j = fene_list(m)
      call rij_and_rr( rij, rr, i, j)
      EE = EE - 1.D0/2 * kFENE * R0_2 * log( 1 - rr/R0_2 ) / 2
    end do
  end do

end subroutine FENE_energy


subroutine update_verlet_list
  !--------------------------------------!
  !Judge whether renew verlet list or not
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   Nq
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none

  if ( mod(step, nint(rsk_lj/dr)*2) == 0 ) then
    call build_lj_verlet_list
  end if

  if ( mod(step, nint(rsk_real/dr)*2) == 0 .and. real_verlet == 1) then
    call build_real_verlet_list
  end if

end subroutine update_verlet_list


subroutine Delta_Energy(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
	real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute energy of LJ potential
  call Delta_LJ_Energy(DeltaE)
  !
  !Compute Delta energy of FENE potential
  if ( ip <= Npe ) then
    call Delta_FENE_Energy(DeltaE)
  end if

end subroutine Delta_Energy


subroutine Delta_Energy_time( DeltaE, time )
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  real*8, dimension(3), intent(inout) :: time
  real*8 :: st_lj, fn_lj, st_real, fn_real
  real*8 :: st_Fourier, fn_Fourier

  DeltaE = 0
  !
  !Compute energy of LJ potential
  call cpu_time(st_lj)
  call Delta_LJ_Energy(DeltaE)
  call cpu_time(fn_lj)
  time(1) = time(1) + fn_lj - st_lj
  !
  !Compute Delta energy of FENE potential
  if ( ip <= Npe ) then
    call Delta_FENE_Energy(DeltaE)
  end if

end subroutine Delta_Energy_time


subroutine Delta_lj_Energy(DeltaE)
  !--------------------------------------!
  !Compute change of LJ potential Energy.
  !   
  !Input
  !   DeltaE
  !Output
  !   DeltaE
  !External Variables
  !   pos, lj_pair_list, lj_point
  !   pos_ip0, pos_ip1, ip
  !   Lx, Ly, Lz, sigma, epsilon, rc_lj
  !Routine Referenced:
  !
  !Reference:
  !In fact, the cut-off radius in good solvent is 2^(1/6), 
  !which was first obatained by JOHN D. WEEKS, DAVID CHANDLER
  !and HANS C. ANDERSEN. So it is called WCA potential.
  !JOHN D. WEEKS, DAVID CHANDLER and HANS C. ANDERSEN, 'Role of 
  !Repulsive Forces in Determining the Equilibrium Structure of
  !Simple Liquids', THE JOURNAL OF CHEMICAL PHYSICS, Vol. 54, 
  !pp.5237-5247, (1971).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: EE, h_lx, h_ly, nh_lx, nh_ly, sigma2, rc_lj2
  real*8  :: rij(3), rr, inv_rr2, inv_rr6, inv_rr12
  integer :: i, j, k, l

  EE     = 0
  h_lx   = Lx / 2
  h_ly   = Ly / 2
  nh_lx  = - h_lx
  nh_ly  = - h_ly
  sigma2 = sigma * sigma
  rc_lj2 = rc_lj *rc_lj
  !
  !ip can't be 1 because it can't 
  !be the anchored particles.
  k = lj_point( ip-1 ) + 1
  l = lj_point( ip )
  do j= k, l
    i = lj_pair_list(j)
    !
    !Energy of old configuration
    !
    rij = pos(i, 1:3) - pos_ip0(1:3)
    !
    !periodic condition
    if ( rij(1) > h_lx ) then
      rij(1) = rij(1) - Lx
    elseif ( rij(1) < nh_lx ) then
      rij(1) = rij(1) + Lx
    end if
    if ( rij(2) > h_ly ) then
      rij(2) = rij(2) - Ly
    elseif ( rij(2) < nh_ly ) then
      rij(2) = rij(2) + Ly
    end if
    !
    !lj energy
    rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
    if ( rr < rc_lj2 ) then
      inv_rr2  = sigma2 / rr
      inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
      inv_rr12 = inv_rr6 * inv_rr6
      EE       = EE + inv_rr6 - inv_rr12 - 0.25D0
    end if
    !
    !Energy of new configuration
    !
    rij = pos(i, 1:3) - pos_ip1(1:3)
    !
    !periodic condition
    if ( rij(1) > h_lx ) then
      rij(1) = rij(1) - Lx
    elseif ( rij(1) < nh_lx ) then
      rij(1) = rij(1) + Lx
    end if
    if ( rij(2) > h_ly ) then
      rij(2) = rij(2) - Ly
    elseif ( rij(2) < nh_ly ) then
      rij(2) = rij(2) + Ly
    end if
    !
    !lj energy
    rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
    if ( rr < rc_lj2 ) then
      inv_rr2  = sigma2 / rr
      inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
      inv_rr12 = inv_rr6 * inv_rr6
      EE       = EE + inv_rr12 - inv_rr6 + 0.25D0
    end if
  end do

  DeltaE = DeltaE + 4 * epsilon * EE

  !
  !Energy of particle and wall of old configuration
  !0.4^(1/6) = 0.86
  if ( pos_ip0(3) < 0.86 ) then
    rr = pos_ip0(3)
    DeltaE = DeltaE - epsilon * ( 2.D0/15 * ( sigma / rr )**9 &
                            - ( sigma / rr )**3 )
  elseif ( pos_ip0(3) > (Lz - 0.86) ) then
    rr = Lz-pos_ip0(3)
    DeltaE = DeltaE - epsilon * ( 2.D0/15 * ( sigma / rr )**9 &
                            - ( sigma / rr )**3 )
  end if
  !
  !Energy of particle and wall of new configuration
  if ( pos_ip1(3) < 0.86 ) then
    rr = pos_ip1(3)
    DeltaE = DeltaE + epsilon * ( 2.D0/15 * ( sigma / rr )**9 &
                            - ( sigma / rr )**3 )
  elseif ( pos_ip1(3) > (Lz - 0.86) ) then
    rr = Lz-pos_ip1(3)
    DeltaE = DeltaE + epsilon * ( 2.D0/15 * ( sigma / rr )**9 &
                            - ( sigma / rr )**3 )
  end if

end subroutine Delta_lj_Energy


subroutine Delta_FENE_Energy(DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   DeltaE
  !Output
  !   DeltaE
  !External Variables
  !   pos, fene_list, fene_point
  !   pos_ip0, pos_ip1, ip
  !   Lx, Ly, Lz, kFENE, R0_2
  !Routine Referenced:
  !
  !Reference:
  !1.Michael Murat, Gary S. Grest, 'Structure of a Grafted
  !Polymer Brush: A Molecular Dynamics Simulation', Macromolecules,
  !vol. 22, pp.4054-4059, (1989).
  !The FENE potential was first derived on 1972.
  !2.Harold R. Warner, 'Kinetic Theory and Rheology of Dilute 
  !Suspensions of Finitely Extendible Dumbbells', Ind. Eng. Chem.
  !Fundam., Vol. 11, No. 3, pp.379-387, (1972).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  integer :: i, j, k, l, m
  real*8  :: rr, rij(3), EE

  EE=0
  !
  !ip can't be 1 because it can't 
  !be the anchored particles.
  l = fene_point(ip-1)+1
  m = fene_point(ip)
  do k= l, m
    i = fene_list(k)
    if (i==ip) then
      write(*,*) i
      stop 
    end if
    !
    !Energy of FENE potential of odd configuration
    !
    rij = pos(i,1:3) - pos_ip0(1:3)
    !
    !Peridoic condition
    if ( rij(1) > Lx/2 ) then
      rij(1) = rij(1) - Lx
    elseif ( rij(1) < -Lx/2 ) then
      rij(1) = rij(1) + Lx
    end if
    if ( rij(2) > Ly/2 ) then
      rij(2) = rij(2) - Ly
    elseif ( rij(2) < -Ly/2 ) then
      rij(2) = rij(2) + Ly
    end if
    rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
    EE = EE + log( 1 - rr / R0_2 )
    !
    !Energy of FENE potential of new configuration
    !
    rij = pos(i,1:3) - pos_ip1(1:3)
    !
    !Periodic condition
    if ( rij(1) > Lx/2 ) then
      rij(1) = rij(1) - Lx
    elseif ( rij(1) < -Lx/2 ) then
      rij(1) = rij(1) + Lx
    end if
    if ( rij(2) > Ly/2 ) then
      rij(2) = rij(2) - Ly
    elseif ( rij(2) < -Ly/2 ) then
      rij(2) = rij(2) + Ly
    end if
    rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
    if ( rr > R0_2 ) then
      write(*,*) 'Chemical bonds are Broken off!'
      stop
    end if
    EE = EE - log( 1 - rr / R0_2 )
  enddo
  DeltaE = DeltaE + 0.5D0 * KFENE * R0_2 * EE

end subroutine Delta_FENE_Energy


subroutine read_energy_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  open(unit=100, file='energy_data.txt')
    read(100,*) epsilon
    read(100,*) sigma
    read(100,*) rc_lj
    read(100,*) rv_lj
    read(100,*) rsk_lj
    read(100,*) R0_2
    read(100,*) kFENE
  close(100)

end subroutine read_energy_parameters


subroutine write_energy_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  write(*,*) '******************  Potential  *********************'
  write(*,*) 'rc_lj      :', rc_lj
  write(*,*) 'rv_lj      :', rv_lj
  write(*,*) '****************************************************'

end subroutine write_energy_parameters


subroutine initialize_lj_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rho, v_verlet

  !
  !allocate verlet list of LJ potential
  if ( allocated(lj_point) ) deallocate(lj_point)
  allocate(  lj_point(NN)  )
  lj_point   = 0
  rho = NN / (Lx * Ly * Lz)
  v_verlet = 8.D0/3 * pi * rv_lj**3
  if ( allocated(lj_pair_list) ) deallocate(lj_pair_list)
  allocate(  lj_pair_list(25*NN*ceiling(rho*v_verlet))  )
  lj_pair_list = 0

end subroutine initialize_lj_parameters


subroutine build_lj_verlet_list
  !--------------------------------------!
  !Construct lj_pair_list and lj_point by the link list
  !method with the complexity of O(N)
  !   
  !Input
  !   pos
  !Output
  !   lj_pair_list, lj_point
  !External Variables
  !   NN, Lx, Ly, Lz, rv_lj, lj_pair_list, lj_point, pos
  !Routine Referenced:
  !1. rij_and_rr(rij, rr, i, j)
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.550-552. 
  !--------------------------------------!
  use global_variables
  implicit none
  integer i,j,k,l,m,n,p,q,r,maxnab
  integer icel,jcel,kcel,ncel1,ncel2,ncel3
  real*8, dimension(3) :: rij
  real*8 :: rsqr,rcel1,rcel2,rcel3
  integer, dimension(NN) :: cell_list
  integer,allocatable,dimension(:,:,:)::hoc

  ncel1=int(Lx/rv_lj)   !number of cell in x direction
  ncel2=int(Ly/rv_lj)   !number of cell in y direction
  ncel3=int(Lz/rv_lj)   !number of cell in z direction
  allocate(hoc(0:ncel1-1,0:ncel2-1,0:ncel3-1))

  maxnab=size(lj_pair_list)
  hoc=0
  rcel1=Lx/ncel1      !Size of each cell in x direction
  rcel2=Ly/ncel2      !Size of each cell in y direction
  rcel3=Lz/ncel3      !Size of each cell in z direction
  do i=1,NN
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    cell_list(i)=hoc(icel,jcel,kcel)
    hoc(icel,jcel,kcel)=i
  end do

  k=0
  do i=1,NN
    icel=int((pos(i,1)+Lx/2)/rcel1)  
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    do l=-1,1
      if (icel+l .ge. ncel1) then
        p=icel+l-ncel1
      elseif(icel+l<0) then
        p=icel+l+ncel1
      else
        p=icel+l
      end if
      do m=-1,1
        if (jcel+m .ge. ncel2) then
          q=jcel+m-ncel2
        elseif(jcel+m<0) then
          q=jcel+m+ncel2
        else
          q=jcel+m
        end if
        do n=-1,1
          if (kcel+n .ge. ncel3) then
            cycle
          elseif(kcel+n<0) then
            cycle
          else
            r=kcel+n
          end if
          j=hoc(p,q,r)
          do while (j /= 0)
            call rij_and_rr(rij,rsqr,i,j)
            if ( i/=j .and. rsqr<(rv_lj*rv_lj) ) then
              k = k + 1
              if ( k > maxnab ) then
                write(*,*) 'maxnab', maxnab
                write(*,*) 'k',  k
                write(*,*) 'lj verlet list is too small!'
                stop
              end if
              lj_pair_list(k)=j
            end if
            j=cell_list(j)
          end do
        end do
      end do
    end do
    lj_point(i)=k
  end do
  npair1=k
  deallocate(hoc)
end subroutine build_lj_verlet_list


subroutine build_real_verlet_list
  !--------------------------------------!
  !Construct real_pair_list and real_point by the link list
  !method with the complexity of O(N)
  !   
  !Input
  !   pos
  !Output
  !   real_pair_list, real_point
  !External Variables
  !   Nq, Lx, Ly, Lz, rv_real,
  !   real_pair_list, real_point, pos, charge
  !Routine Referenced:
  !1. rij_and_rr(rij, rr, i, j)
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.550-552. 
  !--------------------------------------!
  use global_variables
  implicit none

  integer i,j,k,m,n,p,q,r,u,v,w,maxnab
  integer icel,jcel,kcel,ncel1,ncel2,ncel3
  real*8, dimension(3) :: rij
  real*8 :: rsqr,rcel1,rcel2,rcel3
  integer, dimension(Nq) :: cell_list
  integer,allocatable,dimension(:,:,:)::hoc
  ncel1=int(Lx/rv_real)
  ncel2=int(Ly/rv_real)
  ncel3=int(Lz/rv_real)
  allocate(hoc(0:ncel1-1,0:ncel2-1,0:ncel3-1))

  hoc=0
  maxnab=size(real_pair_list)
  rcel1=Lx/ncel1
  rcel2=Ly/ncel2
  rcel3=Lz/ncel3
  do m=1,Nq
    i=charge(m)
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    cell_list(m)=hoc(icel,jcel,kcel)
    hoc(icel,jcel,kcel)=m
  end do

  k=0
  do m=1,Nq
    i=charge(m)
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    do u=-1,1
      if (icel+u>=ncel1) then
        p=icel+u-ncel1
      elseif(icel+u<0) then
        p=icel+u+ncel1
      else
        p=icel+u
      end if
      do v=-1,1
        if (jcel+v>=ncel2) then
          q=jcel+v-ncel2
        elseif(jcel+v<0) then
          q=jcel+v+ncel2
        else
          q=jcel+v
        end if
        do w=-1,1
          if (kcel+w>=ncel3) then
            cycle
          elseif(kcel+w<0) then
            cycle
          else
            r=kcel+w
          end if
          n=hoc(p,q,r)
          do while (n /= 0)
            j=charge(n)
            call rij_and_rr(rij,rsqr,i,j)
            if ( i/=j .and. rsqr<(rv_real*rv_real) ) then
              k=k+1
              if ( k > maxnab ) then
                write(*,*) 'maxnab', maxnab
                write(*,*) 'k',  k
                write(*,*) 'real verlet list is too small!'
                stop
              end if
!               real_pair_list(k,1)=i
!               real_pair_list(k,2)=j
              real_pair_list(k) = j
            end if
            n = cell_list(n)
          end do
        end do
      end do
    end do
    real_point(m) = k
  end do
  npair2 = k
  deallocate(hoc)

end subroutine build_real_verlet_list


subroutine build_fene_list
  !--------------------------------------!
  !Construct the fene_list array.
  !   
  !Input
  !   
  !Output
  !   fene_list, fene_point
  !External Variables
  !   N_bond, Npe, Nml, Ngl
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l

  N_bond = Ngl * (Nml-1) * 2
  allocate( fene_list(N_bond) )
  allocate( fene_point(Npe) )

  l = 0
  do i = 1, Npe
    if ( mod( i, Nml ) == 1 ) then
      l = l + 1
      fene_list(l)  = i + 1 
      fene_point(i) = l
    elseif( mod( i, Nml ) == 0 ) then
      l = l + 1
      fene_list(l)  = i - 1
      fene_point(i) = l
    else
      l = l + 1
      fene_list(l)  = i - 1
      l = l + 1
      fene_list(l)  = i + 1
      fene_point(i) = l
    end if
  end do

end subroutine build_fene_list


subroutine rij_and_rr(rij, rsqr, i, j)
  !-----------------------------------------!
  !compute displacement vector and displacement of two particles
  !input:
  !  post(pos or pos1), i, j(particle number) 
  !output:
  !  rij(displacement vecter), rr(square of displacement)
  !External Variant:
  !  Lz(used in period condition)
  !note:
  !  including period condition
  !-----------------------------------------!
  use global_variables
  implicit none
  real*8, dimension(3), intent(out) :: rij
  real*8, intent(out) :: rsqr
  integer, intent(in) :: i
  integer, intent(in) :: j

  rij = pos(i,1:3) - pos(j,1:3)

  if ( rij(1) > Lx/2 ) then
    rij(1) = rij(1) - Lx
  elseif( rij(1) <= -Lx/2 ) then
    rij(1) = rij(1) + Lx
  end if
  if ( rij(2) > Ly/2 ) then
    rij(2) = rij(2) - Ly
  elseif( rij(2) <= -Ly/2 ) then
    rij(2) = rij(2) + Ly
  end if

  rsqr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

end subroutine rij_and_rr


end module compute_energy


















