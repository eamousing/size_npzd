module mod_zooplankton

   use kind_parameters
   use parameters
   use utilities
   use mod_phytoplankton, only: phytoplankton

   implicit none

   private
   public :: zooplankton, init_zoo, tstep_zoo

   type zooplankton
      !! Zooplankton type definition
      integer :: nb
      real(dp), allocatable :: length(:)
      real(dp), allocatable :: width(:)
      real(dp), allocatable :: volume(:)

      real(dp), allocatable :: biom_c(:)
      real(dp), allocatable :: biom_no3(:)

      ! Predation
      real(dp), allocatable :: gmax(:)

   end type zooplankton

contains

   subroutine init_zoo(self, nb_groups)
      type(zooplankton), intent(inout) :: self
      integer, intent(in) :: nb_groups

      self%nb = nb_groups

      ! Allocate arrays
      allocate(self%length(nb_groups))
      allocate(self%width(nb_groups))
      allocate(self%volume(nb_groups))
      allocate(self%biom_c(nb_groups))
      allocate(self%biom_no3(nb_groups))
      allocate(self%gmax(nb_groups))

      ! Set lengths, widths and volumes
      call set_lengths(nb_groups, 20.0_dp, 500.0_dp, self%length)
      call set_widths(self%length, self%width)
      call set_volumes(self%length, self%width, self%volume)

      ! Initialize carbon and nitrate biomass
      self%biom_c = zoo_c0
      self%biom_no3 = zoo_n0

      ! Set allometric parameters
      call set_allometric(self%volume, a_gmax, b_gmax, self%gmax)



   end subroutine init_zoo

   subroutine tstep_zoo(self, phyto, te, no3, dt)
      type(zooplankton), intent(inout) :: self
      type(phytoplankton), intent(inout) :: phyto
      real(dp), intent(in) :: te
      real(dp), intent(inout) :: no3
      real(dp), intent(in) :: dt

      real(dp) :: gc_p(phyto%nb, self%nb)
      real(dp) :: gc_z(self%nb, self%nb)
      real(dp) :: assim_no3(self%nb)
      real(dp) :: assim_c(self%nb)
      real(dp) :: ege_no3(self%nb)

      real(dp) :: gamma_t
      integer :: i, j

      ! Calculate temperature dependent
      gamma_t = exp(tsens * (te - tref))

      call grazing(self, phyto, gamma_t, gc_p, gc_z)
      
      call assimilation(self, gc_p, gc_z, assim_no3, assim_c, ege_no3, dt)

      ! Update
      do i = 1, self%nb
         self%biom_no3(i) = self%biom_no3(i) + assim_no3(i)
         self%biom_c(i) = self%biom_c(i) + assim_c(i)
         no3 = no3 + sum(ege_no3)
         do j = 1, phyto%nb
            phyto%biom_c(j) = phyto%biom_c(j) - (gc_p(i,j) * self%biom_c(i) * dt)
            phyto%biom_no3(j) = phyto%biom_no3(j) - (gc_p(i,j) * self%biom_c(i) * dt) * (phyto%biom_no3(j) / phyto%biom_c(j))
            if (phyto%biom_c(j) .lt. phy_c0) phyto%biom_c = phy_c0
            if (phyto%biom_no3(j) .lt. phy_n0) phyto%biom_no3 = phy_n0
         end do
         do j = 1, self%nb
            self%biom_c(j) = self%biom_c(j) - (gc_p(i,j) * self%biom_c(i) * dt)
            self%biom_no3(j) = self%biom_no3(j) - (gc_p(i,j) * self%biom_c(i) * dt) * (self%biom_no3(j) / self%biom_c(j))
            if (self%biom_c(j) .lt. zoo_c0) self%biom_c = zoo_c0
            if (self%biom_no3(j) .lt. zoo_c0) self%biom_no3 = zoo_n0
         end do
      end do

      ! Remove grazing/predation

   end subroutine tstep_zoo

   subroutine grazing(self, phyto, gamma_t, gc_p, gc_z)
      type(zooplankton), intent(in) :: self
      type(phytoplankton), intent(in) :: phyto
      real(dp), intent(in) :: gamma_t
      real(dp), intent(out) :: gc_p(:,:)
      real(dp), intent(out) :: gc_z(:,:)

      integer :: i, j
      real(dp) :: p_ratio(phyto%nb)
      real(dp) :: z_ratio(self%nb)
      real(dp) :: phi_p(phyto%nb)
      real(dp) :: phi_z(self%nb)
      real(dp) :: pref_p
      real(dp) :: pref_z
      real(dp) :: tot_c

      real(dp) :: phyto_biom
      real(dp) :: gc, gc_tot

      

      ! For each zooplankton predator do
      do i = 1, self%nb
         ! Calculate predator-prey volume ratio
         phyto_biom = sum(phyto%biom_no3(:))
         p_ratio = self%length(i) / phyto%length
         z_ratio = self%length(i) / self%length
         
         ! Calculate prey palatability
         phi_p = exp(-(log((p_ratio)/(opt)))**2.0_dp * (2.0_dp * sigma**2.0_dp)**(-1.0_dp))
         phi_z = exp(-(log((z_ratio)/(opt)))**2.0_dp * (2.0_dp * sigma**2.0_dp)**(-1.0_dp))

         ! Calculate prey preference
         pref_p = sum(phi_p*phyto%biom_c**2.0_dp) / (sum(phi_p*phyto%biom_c**2.0_dp) + sum(phi_z*self%biom_c**2.0_dp))
         pref_z = sum(phi_z*self%biom_c**2.0_dp) /  (sum(phi_p*phyto%biom_c**2.0_dp) + sum(phi_z*self%biom_c**2.0_dp))
         
         ! Calculate total amount of carbon available to the predator
         tot_c = sum(phi_p * phyto%biom_c + phi_z * self%biom_c)

         ! Calculate predation biomass-specific grazing rate for each prey type
         gc_p(:,i) = self%gmax(i) * gamma_t * (phi_p * phyto%biom_c / (tot_c + kc)) * (1 - exp(aref * tot_c)) * pref_p
         gc_z(:,i) = self%gmax(i) * gamma_t * (phi_z * self%biom_c / (tot_c + kc)) * (1 - exp(aref * tot_c)) * pref_z

         gc_tot = 0.0_dp
         do j = 1, phyto%nb
            gc = self%gmax(i) * gamma_t * (phi_p(j) * phyto%biom_c(j) / (tot_c + kc)) * (1 - exp(aref * tot_c)) * pref_p
            gc_tot = gc_tot + gc
            ! if (phyto_biom .gt. 1.0_dp) print *, gc, gc_tot
         end do

         
      end do
   end subroutine grazing

   subroutine assimilation(self, gc_p, gc_z, assim_no3, assim_c, ege_no3, dt)
      type(zooplankton), intent(in) :: self
      real(dp), intent(in) :: gc_p(:,:)
      real(dp), intent(in) :: gc_z(:,:)
      real(dp), intent(out) :: assim_no3(:)
      real(dp), intent(out) :: assim_c(:)
      real(dp), intent(out) :: ege_no3(:)
      real(dp), intent(in) :: dt

      real(dp) :: lambda_no3
      real(dp) :: lambda_c
      real(dp) :: qno3
      integer :: i
      real(dp) :: tot_pred

      
      ! Loop over all predators
      do i = 1, self%nb
         ! Calculate total biomass-specific grazing
         tot_pred = sum(gc_p(i,:)*self%biom_c) + sum(gc_z(i,:)*self%biom_c) * dt
         
         ! Calculate nitrogen assimilation efficiency
         qno3 = self%biom_no3(i) / self%biom_c(i)
         lambda_no3 = lambda_max * ((qmaxno3 - qno3) / (qmaxno3 - qminno3))
         if (lambda_no3 .gt. 1.0_dp) lambda_no3 = 1.0_dp
         if (lambda_no3 .lt. 0.0_dp) lambda_no3 = 0.0_dp

         ! Nitrogen assimilation
         assim_no3(i) = tot_pred * lambda_no3
         ege_no3(i) = tot_pred * (1.0_dp - lambda_no3)

         ! Carbon assimilation efficiency
         lambda_c = lambda_max * ((qno3 - qminno3)/(qmaxno3 - qminno3))
         if (lambda_c .gt. 1.0_dp) lambda_c = 1.0_dp
         if (lambda_c .lt. 0.0_dp) lambda_c = 0.0_dp

         ! Carbon assimilation
         assim_c(i) = tot_pred * lambda_c
      end do
   end subroutine assimilation

end module mod_zooplankton