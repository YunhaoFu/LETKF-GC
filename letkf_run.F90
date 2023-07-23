program letkf_run

    use precision
    use ini_fnl
    use prepost
    use letkf

    implicit none

include './parameter.h'

    character(len=1024)      ::    bkg_dir
    character(len=1024)      ::    ana_dir
    character(len=1024)      ::    obs_dir
    integer                  ::    ens_size

    integer                  ::    nobs
    real    , allocatable    ::    olat(:)
    real    , allocatable    ::    olon(:)
    real    , allocatable    ::    hdxb(:,:)
    real    , allocatable    ::    error(:)
    real    , allocatable    ::    omb(:)
    real                     ::    qc_ens

    real(r4), allocatable    ::    x_ens_atm(:,:,:,:)
    real(r4), allocatable    ::    lonxy_atm(:,:)
    real(r4), allocatable    ::    latxy_atm(:,:)

    real(r8), allocatable    ::    x_ens_lnd(:,:,:)
    integer, allocatable     ::    grid2patch_start(:,:)
    integer, allocatable     ::    grid2patch_count(:,:)
    real(r4), allocatable    ::    lonxy_lnd(:,:)
    real(r4), allocatable    ::    latxy_lnd(:,:)
    integer , allocatable    ::    itypwat(:,:)                  ! todo_colm_invar #1
    integer , allocatable    ::    lb_patch(:,:)                 ! todo_colm_var #

    real                     ::    omb_max
    real                     ::    infl
    real                     ::    radius(2)

    call da_init(bkg_dir,ana_dir,obs_dir,ens_size,    &
                 qc_ens,omb_max,infl,radius           )          

    call readin(bkg_dir,obs_dir,ens_size,                     &
                nobs,olat,olon,hdxb,error,omb,qc_ens,omb_max, &
                x_ens_atm, lonxy_atm, latxy_atm,              &
                x_ens_lnd, grid2patch_start, grid2patch_count,&
                lonxy_lnd, latxy_lnd,                         &
                itypwat,lb_patch                              ) ! todo_colm_invar #2
    
    call letkf_ini(ens_size,nobs,olat,olon,hdxb,error,omb,    &
                   x_ens_atm,lonxy_atm,latxy_atm,             &
                   x_ens_lnd,lonxy_lnd,latxy_lnd,itypwat      ) ! todo_colm_invar #3

    call letkf_drv(ens_size,nobs,olat,olon,hdxb,error,omb,    &
                   x_ens_atm,lonxy_atm,latxy_atm,             &
                   x_ens_lnd,lonxy_lnd,latxy_lnd,itypwat,     & ! todo_colm_invar #4
                   grid2patch_start,grid2patch_count,         &
                   infl,radius)

    call letkf_fnl(nobs,olat,olon,hdxb,error,omb)

    call writeout(nobs,ana_dir,ens_size,x_ens_atm,x_ens_lnd,  &
                  itypwat,lb_patch                            )
    
    call da_final(nobs,lonxy_atm,latxy_atm,           &
                  grid2patch_start, grid2patch_count, &
                  lonxy_lnd,latxy_lnd,                &
                  itypwat,lb_patch                    )     ! todo_colm_invar #5

end program letkf_run
