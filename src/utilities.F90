module utilities

   use kind_parameters
   use parameters, only: pi

   implicit none

   public

contains

   subroutine set_lengths(nb, smin, smax, lengths)
      integer, intent(in) :: nb
      real(dp), intent(in) :: smin
      real(dp), intent(in) :: smax
      real(dp), intent(out) :: lengths(:)

      real(dp) :: inc
      integer :: i

      lengths(1) = smin
      inc = (smax - smin) / real((nb - 1), dp)
      do i = 2, nb
         lengths(i) = lengths(i-1) + inc
      end do
   end subroutine set_lengths

   subroutine set_widths(lengths, widths)
      real(dp), intent(in) :: lengths(:)
      real(dp), intent(out) :: widths(:)

      widths = 0.49829_dp * lengths + 4.42329_dp
   end subroutine set_widths

   subroutine set_volumes(lengths, widths, volumes)
      real(dp), intent(in) :: lengths(:)
      real(dp), intent(in) :: widths(:)
      real(dp), intent(out) :: volumes(:)

      volumes = (4.0_dp / 3.0_dp) * pi * (lengths / 2.0_dp) * (widths / 2.0_dp)**2.0_dp
   end subroutine set_volumes

   subroutine set_allometric(volumes, a, b, trait)
      real(dp), intent(in) :: volumes(:)
      real(dp), intent(in) :: a
      real(dp), intent(in) :: b
      real(dp), intent(out) :: trait(:)

      trait = a * volumes**b
   end subroutine set_allometric

end module utilities