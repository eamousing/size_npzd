module mod_phytoplankton

   use kind_parameters
   use parameters
   use utilities

   implicit none

   private
   public :: phytoplankton, init_phyto, tstep_phyto

   type phytoplankton
      !! Phytoplankton type definition
      integer :: nb
      real(dp), allocatable :: length(:) !! Length (um)
      real(dp), allocatable :: width(:) !! Width (um)
      real(dp), allocatable :: volume(:) !! Volume (um**3)
      
      real(dp), allocatable :: biom_c(:) !! Carbon biomass (mmol C (m)-3)
      real(dp), allocatable :: biom_no3(:) !! NO3 biomass (mmol N (m)-3)
      real(dp), allocatable :: chla(:) !! Chlorophyll (mg (m-3))

      real(dp), allocatable :: qno3(:) !! Cellular NO3 reserve (mmol N (mmol C)-1) 

      ! Nutrient uptake
      real(dp), allocatable :: qminno3(:)
      real(dp), allocatable :: qmaxno3(:)
      real(dp), allocatable :: vmaxno3(:) !! Maximum NO3 uptake rate (mmol N (mmol C)-1 d-1)
      real(dp), allocatable :: kno3(:) !! NO3 Half-saturation constant (mmol N m**-3)

      ! Photosynthetic traits
      real(dp), allocatable :: pmax(:) !! Max photosynthetic rate (d-1)
   end type phytoplankton

contains

   subroutine init_phyto(self, nb_groups)
      type(phytoplankton), intent(inout) :: self
      integer, intent(in) :: nb_groups

      self%nb = nb_groups

      ! Allocate arrays
      allocate(self%length(nb_groups))
      allocate(self%width(nb_groups))
      allocate(self%volume(nb_groups))
      allocate(self%biom_c(nb_groups))
      allocate(self%biom_no3(nb_groups))
      allocate(self%chla(nb_groups))
      allocate(self%qminno3(nb_groups))
      allocate(self%qmaxno3(nb_groups))
      allocate(self%vmaxno3(nb_groups))
      allocate(self%kno3(nb_groups))
      allocate(self%pmax(nb_groups))

      ! Set lengths, width and volume
      call set_lengths(nb_groups, 2.0_dp, 200.0_dp, self%length)
      call set_widths(self%length, self%width)
      call set_volumes(self%length, self%width, self%volume)

      ! Initialize carbon and nitrate biomass in Redfield proportions
      self%biom_c = phy_c0
      self%biom_no3 = phy_n0

      self%chla = phy_chla

      ! Set Nutrient uptake traits
      call set_allometric(self%volume, a_vmaxno3, b_vmaxno3, self%vmaxno3)
      call set_allometric(self%volume, a_qminno3, b_qminno3, self%qminno3)
      call set_allometric(self%volume, a_qmaxno3, b_qmaxno3, self%qmaxno3)
      call set_allometric(self%volume, a_kno3, b_kno3, self%kno3)
      call set_allometric(self%volume, a_pmax, b_pmax, self%pmax)

   end subroutine init_phyto

   subroutine tstep_phyto(self, te, no3, irr, dt)
      type(phytoplankton), intent(inout) :: self
      real(dp), intent(in) :: te
      real(dp), intent(inout) :: no3
      real(dp), intent(in) :: irr
      real(dp), intent(in) :: dt !! time step length

      real(dp) :: gamma_t
      real(dp) :: qno3(self%nb)
      real(dp) :: vno3(self%nb)
      real(dp) :: vc(self%nb)
      real(dp) :: vchla(self%nb)

      ! Calculate temperature dependent
      gamma_t = exp(tsens * (te - tref))

      ! Nutrient uptake
      call nutrient_uptake(self, no3, gamma_t, qno3, vno3)

      ! Photosynthesis
      call photosynthesis(self, qno3, irr, vno3, gamma_t, vc, vchla)

      ! Update biomasses
      call update_biom(self, vc, vno3, vchla, no3, dt)

   end subroutine tstep_phyto

   subroutine nutrient_uptake(self, no3, gamma_t, qno3, vno3)
      type(phytoplankton), intent(in) :: self
      real(dp), intent(in) :: no3
      real(dp), intent(in) :: gamma_t
      real(dp), intent(out) :: qno3(:)
      real(dp), intent(out) :: vno3(:)
      
      real(dp) :: qstat_no3(self%nb)
      real(dp) :: no3_lim(self%nb)
      integer :: i

      ! Calculate internal quotas
      qno3 = self%biom_no3 / self%biom_c

      ! Calculate uptake capacity
      qstat_no3 = (self%qmaxno3 - qno3) / (self%qmaxno3 - self%qminno3)

      do i = 1, self%nb
         if (qstat_no3(i) > 1.0_dp) qstat_no3(i) = 1.0_dp
         if (qstat_no3(i) < 0.0_dp) qstat_no3(i) = 0.0_dp
      end do

      ! Calculate nutrient uptake
      no3_lim = (no3 / (no3 + self%kno3))
      do i = 1, self%nb
         if (no3_lim(i) .lt. 0.0_dp) no3_lim(i) = 0.0_dp
      end do
      vno3 = self%vmaxno3 * no3_lim * qstat_no3 * gamma_t
   end subroutine nutrient_uptake

   subroutine photosynthesis(self, qno3, irr, vno3, gamma_t, vc, vchla)
      type(phytoplankton), intent(inout) :: self
      real(dp), intent(in) :: qno3(:)
      real(dp), intent(in) :: irr
      real(dp), intent(in) :: vno3(:)
      real(dp), intent(in) :: gamma_t
      real(dp), intent(out) :: vc(:)
      real(dp), intent(out) :: vchla(:)

      real(dp) :: pchla(self%nb)
      real(dp) :: gamma_no3(self%nb)
      real(dp) :: pc(self%nb)
      real(dp) :: psat(self%nb)
      real(dp) :: qchla(self%nb)

      integer :: i

      ! Calculate nitrogen quota limitation term
      gamma_no3 = (qno3 - self%qminno3) / (self%qmaxno3 - self%qminno3)

      do i = 1, self%nb
         if (gamma_no3(i) .gt. 1.0_dp) gamma_no3(i) = 1.0_dp
         if (gamma_no3(i) .lt. 0.0_dp) gamma_no3(i) = 0.0_dp
      end do

      ! Calculate carbon-specific light-saturated photosynthetic rate
      psat = self%pmax * gamma_no3 * gamma_t

      ! Calculate internal chl a quota
      qchla = self%chla / self%biom_c

      ! Calculate light-limited photosynthesis
      pc = psat * ( 1.0_dp - exp((-alpha * qchla * irr) / psat))

      ! Calculate net carbon uptake
      vc = pc - zeta * vno3 * self%biom_c

      ! Calculate nitrogen-specific Chl a synthesis
      pchla = omaxno3 * (pc / (alpha * qchla * irr))

      ! Calculate carbon-specific chla a synthesis
      vchla = pchla * vno3
   end subroutine photosynthesis

   subroutine update_biom(self, vc, vno3, vchla, no3, dt)
      type(phytoplankton), intent(inout) :: self
      real(dp), intent(in) :: vc(:)
      real(dp), intent(in) :: vno3(:)
      real(dp), intent(in) :: vchla(:)
      real(dp), intent(inout) :: no3
      real(dp), intent(in) :: dt

      self%biom_c = self%biom_c + (self%biom_c * vc * dt)
      self%biom_no3 = self%biom_no3 + (self%biom_c * vno3 * dt)
      self%chla = self%chla + (self%biom_c * vchla * dt)
      no3 = no3 - sum(self%biom_c * vno3 * dt)
   end subroutine update_biom

end module mod_phytoplankton