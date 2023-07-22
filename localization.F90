module localization
    use precision
    implicit none

    private

    ! public subroutines
    public    ::     get_dist, letkf_loc_gc

    contains

        function get_dist(lat,r) result(dist)

            real(r4), intent(in)    :: lat
            real(r4), intent(in)    :: r(2)
            real(r4)                :: dist
            dist = ( (abs(lat)-0.0)*r(2)+(90.0-abs(lat))*r(1) ) / 90.0

        end function get_dist

        pure function letkf_loc_gc(z, l) result(res)
            real(r4), intent(in) :: z  !< value to localize
            real(r4), intent(in) :: l  !< the equivalent to the gaussian standard deviation
            real(r4) :: res
            real(r4) :: c
            real(r4) :: abs_z, z_c

            c = l / sqrt(0.3_r4)
            abs_z = abs(z)
            z_c = abs_z / c

            if (abs_z >= 2*c) then
               res = 0.0
            elseif (abs_z > c) then
               res = &
                      0.08333_r4 * z_c**5 &
                    - 0.50000_r4 * z_c**4 &
                    + 0.62500_r4 * z_c**3 &
                    + 1.66667_r4 * z_c**2 &
                    - 5.00000_r4 * z_c &
                    + 4_r4 &
                    - 0.66667_r4 * c/abs_z
            else
               res = &
                     -0.25000_r4 * z_c**5 &
                    + 0.50000_r4 * z_c**4 &
                    + 0.62500_r4 * z_c**3 &
                    - 1.66667_r4 * z_c**2 &
                    + 1_r4
            end if
          end function letkf_loc_gc

endmodule localization
