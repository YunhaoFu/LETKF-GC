    !--------------------------------------------------!
    !                   parameter.h                    !
    integer, parameter            ::   nftune         = 14
    integer, parameter            ::   nfcon          = 119
    integer, parameter            ::   nfvar          = 126
    integer, parameter            ::   lon_points_lnd = 1440
    integer, parameter            ::   lat_points_lnd = 720
    integer, parameter            ::   lon_id_lnd     = 2      ! in fcon
    integer, parameter            ::   lat_id_lnd     = 1      ! in fcon
    integer, parameter            ::   rd_invar_cnt   = 1
    ! todo_colm_invar #1
    integer, parameter            ::   invar_idx_s(rd_invar_cnt)&
                                         = (/ 3 /)              ! id for land water type            [itypwat]
    integer, parameter            ::   invar_len  (rd_invar_cnt)&
                                         = (/ 1 /)              ! length for land water type        [itypwat]
    ! todo_colm_var #1
    integer, parameter            ::   scv_idx            = 81  ! id for scv in case wliq, wice is changed
    integer, parameter            ::   nvar_lnd_raw       = 6   ! 1(tss)+1(wliq)+1(wice)+1(tg)+1(tlsun)+1(tlsha)
    integer, parameter            ::   nvar_lnd           = 9   ! 2(tss)+2(wliq)+2(wice)+1(tg)+1(tlsun)+1(tlsha)
    integer, parameter            ::   var_idx_s(nvar_lnd_raw)&     
                                         = (/ 31, &             ! id for Soil temperature           [tss]  
                                              46, &             ! id for Liquid water in layers     [wliq]
                                              61, &             ! id for Ice lens in layers         [wice]
                                              76, &             ! id for Ground surface temperature [tg]
                                              77, &             ! id for Sunlit leaf temperature    [tlsun]
                                              78  &             ! id for Shaded leaf temperature    [tlsha]
                                              /)             
    integer, parameter            ::   var_len  (nvar_lnd_raw)&     
                                         = (/ 15, &             ! length for Soil temperature           [tss]
                                              15, &             ! length for Liquid water in layers     [wliq]
                                              15, &             ! length for Ice lens in layers         [wice]
                                              1 , &             ! length for Ground surface temperature [tg]
                                              1 , &             ! length for Sunlit leaf temperature    [tlsun]
                                              1   &             ! length for Shaded leaf temperature    [tlsha]
                                              /)            
    integer, parameter            ::   itypwat_max = 3          ! todo_colm_var #2
                                                                ! soil => 0 ; 
                                                                ! urban and built-up => 1; 
                                                                ! wetland => 2; 
                                                                ! land ice => 3
    integer, parameter            ::   lon_points_atm = 1440
    integer, parameter            ::   lat_points_atm = 721
    integer, parameter            ::   level_atm      = 89
    real   , parameter            ::   edges_atm      = -90.00
    real   , parameter            ::   edgen_atm      =  90.00
    real   , parameter            ::   edgew_atm      =  00.00
    real   , parameter            ::   edgee_atm      =  359.75
    integer, parameter            ::   head_len       = 101
    ! todo_grapes #1
    integer, parameter            ::   nvar_atm       = 178 ! 89(th) + 89(qv)
