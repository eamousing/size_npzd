program size_npzd

   use kind_parameters
   use mod_phytoplankton
   use mod_zooplankton

   implicit none

   character(len=2) :: nb_phy_in
   character(len=2) :: nb_zoo_in
   integer :: nb_phy
   integer :: nb_zoo
   type(phytoplankton) :: phyto
   type(zooplankton) :: zoo

   real(dp) :: no3 = 10.0_dp
   real(dp) :: te = 283.15_dp
   real(dp) :: irr = 1.3e6_dp
   real(dp) :: dt = 0.0417

   integer :: i

   ! Check if number of trophic group have been specified
   if (command_argument_count() .ne. 2) then
      print *, 'Error: Number of trophic groups of both phytoplankton and zooplankton required'
      stop
   end if

   ! Read and store arguments
   call get_command_argument(1, nb_phy_in)
   call get_command_argument(2, nb_zoo_in)
   read(nb_phy_in, *) nb_phy
   read(nb_zoo_in, *) nb_zoo

   ! Initialize plankton populations
   call init_phyto(phyto, nb_phy)
   call init_zoo(zoo, nb_zoo)

   do i = 1, 20000
      call tstep_phyto(phyto, te, no3, irr, dt)
      call tstep_zoo(zoo, phyto, te, no3, dt)
      print *, i, no3, sum(phyto%biom_no3(:)), sum(zoo%biom_no3(:)), no3 + sum(phyto%biom_no3(:)) + sum(zoo%biom_no3(:))
   end do

end program size_npzd