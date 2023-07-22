module precision
    implicit none
    private
    public :: r8, r4
    ! r8 for colminput
    ! r4 for grapesinput and omb
    integer, parameter:: r8 = selected_real_kind(8)
    integer, parameter:: r4 = selected_real_kind(4)
endmodule precision

