module parameters

   use kind_parameters

   implicit none

   real(dp), parameter :: pi = 3.14159265359_dp

   ! Initialization
   real(dp), parameter :: phy_c0 = 1.0e-9_dp !! mmol C m-3
   real(dp), parameter :: phy_n0 = 1.51e-10_dp !! mmol N m-3
   real(dp), parameter :: phy_chla = 6.28e-10_dp !! mg chla m-3
   real(dp), parameter :: zoo_c0 = 1.0e-9_dp !! mmol C m-3
   real(dp), parameter :: zoo_n0 = 1.51e-10_dp !! mmol N m-3

   ! Physiology
   real(dp), parameter :: tref = 293.15_dp !! degC
   real(dp), parameter :: tsens = 0.05_dp !! no unit

   ! Nutrient uptake
   real(dp), parameter :: a_vmaxno3 = 0.51_dp !! Intercept maximum NO3 uptake rate (mmol N (mmol C)-1 d-1)
   real(dp), parameter :: b_vmaxno3 = -0.27_dp !! Coefficient maximum NO3 uptake rate (mmol N (mmol C)-1 d-1)
   real(dp), parameter :: a_kno3 = 0.17_dp !! Intercept half-sat constant NO3 uptake (mmol N m-3)
   real(dp), parameter :: b_kno3 = 0.27_dp !! Coefficient half-sat constant NO3 uptake (mmol N m-3)
   real(dp), parameter :: a_qmaxno3 = 0.25_dp !! Max N:C quota (mmol N (mmol C)-1)
   real(dp), parameter :: b_qmaxno3 = -0.13_dp !! Max N:C quota (mmol N (mmol C)-1)
   real(dp), parameter :: a_qminno3 = 0.07_dp !! Min N:C quota (mmol N (mmol C)-1)
   real(dp), parameter :: b_qminno3 = -0.17_dp !! Min N:C quota (mmol N (mmol C)-1)

   ! Photosynthesis
   real(dp), parameter :: omaxno3 = 3.0_dp !! Max chla:N ratio (mg chla (mmol N)-1)
   real(dp), parameter :: alpha = 3.83e-7_dp !! Initial slope of photosynthesis-irr curve (mmol C (mg chla)-1 (uE m-2)-1)
   real(dp), parameter :: zeta = 2.33 !! Cost of biosynthesis (mmol C (mmol N)-1)
   real(dp), parameter :: a_pmax = 2.1_dp !! Intercept max photosynthetic rate (d-1)
   real(dp), parameter :: b_pmax = -0.15_dp !! Coefficient max photosynthetic rate (d-1)

   ! Predation
   real(dp), parameter :: a_gmax = 21.9_dp !! d-1
   real(dp), parameter :: b_gmax = -0.16 !! d-1
   real(dp), parameter :: qminno3 = 0.075_dp !! mmol N (mmol C)-1
   real(dp), parameter :: qmaxno3 = 0.151_dp !! mmol N (mmol C)-1
   real(dp), parameter :: sigma = 0.5_dp !! no unit?
   real(dp), parameter :: opt = 10.0_dp !! no unit?
   real(dp), parameter :: kc = 1.0_dp
   real(dp), parameter :: aref = -1.0_dp
   real(dp), parameter :: lambda_max = 0.7_dp
end module parameters