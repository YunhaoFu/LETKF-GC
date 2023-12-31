module letkf

    use precision
    use netcdf
    use kdtree
    use localization

    implicit none

    private

include './parameter.h'

    ! public    ::     letkf_ini, 
    public    ::     letkf_drv, letkf_fnl

    real(r4), parameter                 ::    stdev2max = sqrt(40.0/3.0)
    integer , parameter                 ::    max_obs = 100
    real(r4)                            ::    max_search_dist
    real(r4)                            ::    dist
    real(r4)                            ::    rd

    contains
        
        ! subroutine letkf_ini(ens_size,numpatch,nobs,olat,olon,hdxb,error,omb,&
        !                      x_ens_atm,lonxy_atm,latxy_atm,             &
        !                      x_ens_lnd,lonxy_lnd,latxy_lnd,itypwat      ) ! todo_colm_invar #1
        !     integer              , intent(in)   ::    ens_size
        !     integer              , intent(in)   ::    numpatch
        !     integer              , intent(in)   ::    nobs
        !     real    , allocatable, intent(in)   ::    olat(:)
        !     real    , allocatable, intent(in)   ::    olon(:)
        !     real    , allocatable, intent(in)   ::    hdxb(:,:)
        !     real    , allocatable, intent(in)   ::    error(:)
        !     real    , allocatable, intent(in)   ::    omb(:)
        !     real(r4), allocatable, intent(in)   ::    x_ens_atm(:,:,:,:)
        !     real(r4), allocatable, intent(in)   ::    lonxy_atm(:,:)
        !     real(r4), allocatable, intent(in)   ::    latxy_atm(:,:)
        !     real(r8), allocatable, intent(in)   ::    x_ens_lnd(:,:,:)
        !     real(r4), allocatable, intent(in)   ::    lonxy_lnd(:,:)
        !     real(r4), allocatable, intent(in)   ::    latxy_lnd(:,:)
        !     integer , allocatable, intent(in)   ::    itypwat(:,:)       ! todo_colm_invar #2

        !     if(nobs == 0) return

        ! endsubroutine letkf_ini

        subroutine letkf_drv(ens_size,numpatch,nobs,olat,olon,hdxb,error,omb, &
                             x_ens_atm,lonxy_atm,latxy_atm,             &
                             x_ens_lnd,lonxy_lnd,latxy_lnd,itypwat,     & ! todo_colm_invar #3
                             grid2patch_start,grid2patch_count,         &
                             infl,radius,lb_patch                       )

            integer , intent(in)    ::    ens_size
            integer , intent(in)    ::    numpatch
            integer , intent(in)    ::    nobs
            real    , intent(in)    ::    olat(nobs)
            real    , intent(in)    ::    olon(nobs)
            real    , intent(in)    ::    hdxb(nobs,ens_size)
            real    , intent(in)    ::    error(nobs)
            real    , intent(in)    ::    omb(nobs)
            real(r4),  intent(inout)::    x_ens_atm(ens_size,lon_points_atm,lat_points_atm,nvar_atm)
            real(r4),  intent(in)   ::    lonxy_atm(lon_points_atm,lat_points_atm)
            real(r4),  intent(in)   ::    latxy_atm(lon_points_atm,lat_points_atm)
            real(r8),  intent(inout)::    x_ens_lnd(ens_size,numpatch,nvar_lnd)
            real(r4),  intent(in)   ::    lonxy_lnd(lon_points_lnd,lat_points_lnd)
            real(r4),  intent(in)   ::    latxy_lnd(lon_points_lnd,lat_points_lnd)
            integer ,  intent(in)   ::    itypwat(numpatch,1)       ! todo_colm_invar #4
            integer,   intent(in)   ::    grid2patch_start(lon_points_lnd,lat_points_lnd)
            integer,   intent(in)   ::    grid2patch_count(lon_points_lnd,lat_points_lnd)
            real   ,   intent(in)   ::    infl
            real   ,   intent(in)   ::    radius(2)
            integer,   intent(in)   ::    lb_patch(ens_size,numpatch)
            logical                 ::    skip(numpatch)
            integer                 ::    np

            if(nobs == 0) return

            ! skip when lb in ensemble are different, tested !
            ! skip(1:numpatch) = .false.
            skip(1:numpatch) = maxval(lb_patch(1:ens_size,1:numpatch),dim=1) &
                            /= minval(lb_patch(1:ens_size,1:numpatch),dim=1)
            ! do np=1,numpatch
            !     if(maxval(lb_patch(1:ens_size,np)) /= minval(lb_patch(1:ens_size,np))) skip(np) = .true.
            ! enddo

            print *, '******************* LETKF RUN **********************'

            call atm_da(ens_size,nobs,olat,olon,hdxb,error,omb,    &
                        x_ens_atm,lonxy_atm,latxy_atm,             &
                        infl,radius                                )
            call lnd_da(ens_size,numpatch,nobs,olat,olon,hdxb,error,omb,&
                        x_ens_lnd,lonxy_lnd,latxy_lnd,itypwat,     & ! todo_colm_invar #5
                        grid2patch_start,grid2patch_count,         &
                        infl,radius,skip                           )

            print *, '****************** LETKF OVER **********************'

        endsubroutine letkf_drv

        subroutine atm_da(ens_size,nobs,olat,olon,hdxb,error,omb,    &
                          x_ens_atm,lonxy,latxy,             &
                          infl,radius                                )

            integer , intent(in)    ::    ens_size
            integer , intent(in)    ::    nobs
            real    , intent(in)    ::    olat(nobs)
            real    , intent(in)    ::    olon(nobs)
            real    , intent(in)    ::    hdxb(nobs,ens_size)
            real    , intent(in)    ::    error(nobs)
            real    , intent(in)    ::    omb(nobs)
            real(r4),  intent(inout)::    x_ens_atm(ens_size,lon_points_atm,lat_points_atm,nvar_atm)
            real(r4),  intent(in)   ::    lonxy(lon_points_atm,lat_points_atm)
            real(r4),  intent(in)   ::    latxy(lon_points_atm,lat_points_atm)
            real    , intent(in)    ::    infl
            real    , intent(in)    ::    radius(2)

            integer                 ::    i, j, nvar, ens, ns, idx
            type(kd_root)           ::    obs_tree
            real(r4)                ::    bkg_mean(lon_points_atm,lat_points_atm,nvar_atm)
            real(r4)                ::    bkg(ens_size,lon_points_atm,lat_points_atm,nvar_atm)

            integer                 ::    obs_ij_cnt
            integer                 ::    obs_ij_idx(max_obs)
            real(r4)                ::    obs_ij_dist(max_obs)
            real(r4)                ::    obs_ij_hdxb(ens_size,max_obs)
            real(r4)                ::    obs_ij_rdiag(max_obs)
            real(r4)                ::    obs_ij_dep(max_obs)
            integer                 ::    obs_lg_cnt
            real(r4)                ::    obs_lg_hdxb(ens_size,max_obs)
            real(r4)                ::    obs_lg_rdiag(max_obs)
            real(r4)                ::    obs_lg_dep(max_obs)
            real(r4)                ::    obs_lg_rloc(max_obs)
            real(r4)                ::    trans(ens_size,ens_size)

            call kd_init(obs_tree,olon,olat)
            
            bkg_mean = sum(x_ens_atm,dim=1)/real(ens_size,kind=r4)
            do ens=1,ens_size
                bkg(ens,:,:,:) = x_ens_atm(ens,:,:,:) - bkg_mean(:,:,:)
            enddo
            x_ens_atm(:,:,:,:) = bkg(:,:,:,:)

            ! do nvar=1,nvar_atm
                ! do i=1,lon_points_atm
                    ! do j=1,lat_points_atm
                        ! bkg_mean(i,j,nvar) = sum(x_ens_atm(:,i,j,nvar))/real(ens_size)
                        ! bkg(:,   i,j,nvar) = x_ens_atm(:,i,j,nvar) - bkg_mean(i,j,nvar)
                    ! enddo
                ! enddo
            ! enddo
            lon_loop: do i=1,lon_points_atm
                lat_loop: do j=1,lat_points_atm

                    dist=get_dist(latxy(i,j),radius)
                    max_search_dist = dist * stdev2max
                    ! cycle ?

                    ! get search radius
                    call kd_search_radius(obs_tree, lonxy(i,j),latxy(i,j),&
                         max_search_dist, obs_ij_idx, obs_ij_dist, obs_ij_cnt, .false.)
                        
                    if(obs_ij_cnt > max_obs) stop 26
                    if(obs_ij_cnt .le. 0) cycle

                    ! prepare obs when there is one "nearing" the ij grid point
                    do ns=1,obs_ij_cnt
                        idx                = obs_ij_idx(ns)
                        obs_ij_hdxb (:,ns) = hdxb(:,idx)
                        obs_ij_rdiag  (ns) = error(idx)
                        obs_ij_dep    (ns) = omb(idx)
                    end do

                    obs_lg_cnt = 0
                    obs_loop: do ns=1,obs_ij_cnt
                        rd = letkf_loc_gc(obs_ij_dist(ns),dist)
                        if (rd > 1.0_r4) then
                            obs_lg_cnt                = obs_lg_cnt + 1
                            obs_lg_rloc  (obs_lg_cnt) = rd
                            obs_lg_hdxb(:,obs_lg_cnt) = obs_ij_hdxb (:,ns)
                            obs_lg_rdiag (obs_lg_cnt) = obs_ij_rdiag  (ns)
                            obs_lg_dep   (obs_lg_cnt) = obs_ij_dep    (ns)
                        end if
                    end do obs_loop

                    if (obs_lg_cnt > 0) then
                        !calculate trans matrix in ij grid point format
                        call letkf_core_solve_e0(ens_size, obs_lg_cnt, obs_lg_hdxb(:,:obs_lg_cnt),&
                             obs_lg_rdiag(:obs_lg_cnt), obs_lg_rloc(:obs_lg_cnt), &
                             obs_lg_dep(:obs_lg_cnt), infl, trans)

                        !apply the trans calculated above into different var within the same lat lon, i.e., novertical localization
                        do nvar=1,nvar_atm
                            call sgemm('n', 'n', 1, ens_size, ens_size, 1.0_r4, bkg(:,i,j,nvar), &
                                1, trans, ens_size, 0.0_r4, x_ens_atm(:,i,j,nvar), 1)
                        enddo
                    end if

                enddo lat_loop
            enddo lon_loop

            do ens=1,ens_size
                  x_ens_atm(ens,:,:,:) &
                = x_ens_atm(ens,:,:,:) + bkg_mean(:,:,:)
            enddo

            call kd_free(obs_tree)

        endsubroutine atm_da

        subroutine lnd_da(ens_size,numpatch,nobs,olat,olon,hdxb,error,omb,&
                          x_ens_lnd,lonxy,latxy,itypwat,             & ! todo_colm_invar #6
                          grid2patch_start,grid2patch_count,         &
                          infl,radius,skip                           )
            integer   , intent(in)    ::    ens_size
            integer   , intent(in)    ::    numpatch
            integer   , intent(in)    ::    nobs
            real      , intent(in)    ::    olat(nobs)
            real      , intent(in)    ::    olon(nobs)
            real      , intent(in)    ::    hdxb(nobs,ens_size)
            real      , intent(in)    ::    error(nobs)
            real      , intent(in)    ::    omb(nobs)
            real(r8)   , intent(inout)::    x_ens_lnd(ens_size,numpatch,nvar_lnd)
            real(r4)   , intent(in)   ::    lonxy(lon_points_lnd,lat_points_lnd)
            real(r4)   , intent(in)   ::    latxy(lon_points_lnd,lat_points_lnd)
            integer    , intent(in)   ::    itypwat(numpatch,1)       ! todo_colm_invar #7
            integer    , intent(in)   ::    grid2patch_start(lon_points_lnd,lat_points_lnd)
            integer    , intent(in)   ::    grid2patch_count(lon_points_lnd,lat_points_lnd)
            real       , intent(in)   ::    infl
            real       , intent(in)   ::    radius(2)
            logical    , intent(in)   ::    skip(numpatch)

            integer                   ::    i, j, nvar, ens, ns, idx
            type(kd_root)             ::    obs_tree
            real(r8)                  ::    bkg_mean(numpatch,nvar_lnd)
            real(r8)                  ::    bkg(ens_size,numpatch,nvar_lnd)
            integer                   ::    np

            integer                   ::    obs_ij_cnt
            integer                   ::    obs_ij_idx(max_obs)
            real(r4)                  ::    obs_ij_dist(max_obs)
            real(r4)                  ::    obs_ij_hdxb(ens_size,max_obs)
            real(r4)                  ::    obs_ij_rdiag(max_obs)
            real(r4)                  ::    obs_ij_dep(max_obs)
            integer                   ::    obs_lg_cnt
            real(r4)                  ::    obs_lg_hdxb(ens_size,max_obs)
            real(r4)                  ::    obs_lg_rdiag(max_obs)
            real(r4)                  ::    obs_lg_dep(max_obs)
            real(r4)                  ::    obs_lg_rloc(max_obs)
            real(r4)                  ::    trans(ens_size,ens_size)

            call kd_init(obs_tree,olon,olat)

            bkg_mean = sum(x_ens_lnd,dim=1)/real(ens_size,kind=r4)
            do ens=1,ens_size
                bkg(ens,:,:) = x_ens_lnd(ens,:,:) - bkg_mean(:,:)
            enddo
            x_ens_lnd(:,:,:) = bkg(:,:,:)

            ! do nvar=1,nvar_lnd
            !     do np=1,numpatch
            !         bkg_mean(np,nvar) = sum(x_ens_lnd(:,np,nvar))/real(ens_size)
            !         bkg(:,   np,nvar) = x_ens_lnd(:,np,nvar) - bkg_mean(np,nvar)
            !     enddo
            ! enddo
            lon_loop: do i=1,lon_points_lnd
                lat_loop: do j=1,lat_points_lnd

                    dist=get_dist(latxy(i,j),radius)
                    max_search_dist = dist * stdev2max
                    ! cycle ?

                    ! get search radius
                    call kd_search_radius(obs_tree, lonxy(i,j),latxy(i,j),&
                         max_search_dist, obs_ij_idx, obs_ij_dist, obs_ij_cnt, .false.)
                        
                    if(obs_ij_cnt > max_obs) stop 26
                    if(obs_ij_cnt .le. 0) cycle

                    ! prepare obs when there is one "nearing" the ij grid point
                    do ns=1,obs_ij_cnt
                        idx                = obs_ij_idx(ns)
                        obs_ij_hdxb (:,ns) = hdxb(:,idx)
                        obs_ij_rdiag  (ns) = error(idx)
                        obs_ij_dep    (ns) = omb(idx)
                    end do

                    obs_lg_cnt = 0
                    obs_loop: do ns=1,obs_ij_cnt
                        rd = letkf_loc_gc(obs_ij_dist(ns),dist)
                        if (rd > 1.0_r4) then
                            obs_lg_cnt                = obs_lg_cnt + 1
                            obs_lg_rloc  (obs_lg_cnt) = rd
                            obs_lg_hdxb(:,obs_lg_cnt) = obs_ij_hdxb (:,ns)
                            obs_lg_rdiag (obs_lg_cnt) = obs_ij_rdiag  (ns)
                            obs_lg_dep   (obs_lg_cnt) = obs_ij_dep    (ns)
                        end if
                    end do obs_loop

                    if (obs_lg_cnt > 0) then
                        !calculate trans matrix in ij grid point format
                        call letkf_core_solve_e0(ens_size, obs_lg_cnt, obs_lg_hdxb(:,:obs_lg_cnt),&
                             obs_lg_rdiag(:obs_lg_cnt), obs_lg_rloc(:obs_lg_cnt), &
                             obs_lg_dep(:obs_lg_cnt), infl, trans)

                        !apply the trans calculated above into different var within the same lat lon, i.e., novertical localization
                        do nvar=1,nvar_lnd
                            do np=grid2patch_start(i,j),grid2patch_start(i,j)+grid2patch_count(i,j)-1
                                if(itypwat(np,1) > itypwat_max .or. skip(np)) cycle
                                call dgemm('n', 'n', 1, ens_size, ens_size, 1.0_r8, bkg(:,np,nvar), &
                                    1, real(trans,kind=r8), ens_size, 0.0_r8, x_ens_lnd(:,np,nvar), 1)
                            enddo
                        enddo
                    end if

                enddo lat_loop
            enddo lon_loop

            do ens=1,ens_size
                  x_ens_lnd(ens,:,:) &
                = x_ens_lnd(ens,:,:) + bkg_mean(:,:)
            enddo

            call kd_free(obs_tree)

        endsubroutine lnd_da

        subroutine letkf_fnl(nobs,olat,olon,hdxb,error,omb)
            integer              , intent(in)   ::    nobs
            real    , allocatable, intent(inout)::    olat(:)
            real    , allocatable, intent(inout)::    olon(:)
            real    , allocatable, intent(inout)::    hdxb(:,:)
            real    , allocatable, intent(inout)::    error(:)
            real    , allocatable, intent(inout)::    omb(:)

            if(nobs == 0) then
                print *, 'no observation assimilated'
            else
                print *, 'valid observation count = ', nobs
                print *, 'lon max = ', maxval(olon), ' min = ', minval(olon)
                print *, 'lat max = ', maxval(olat), ' min = ', minval(olat)
                print *, "Hx' max = ", maxval(hdxb), ' min = ', minval(hdxb)
                print *, 'err max = ', maxval(error),' min = ', minval(error)
                print *, 'omb max = ', maxval(omb) , ' min = ', minval(omb) 
            endif

            deallocate(olat)
            deallocate(olon)
            deallocate(hdxb)
            deallocate(error)
            deallocate(omb)

        endsubroutine letkf_fnl

        subroutine letkf_core_solve_e0(nbv, nobs, hdxb, rdiag, rloc, dep, infl, trans)
            integer, intent(in)   :: nbv             !< ens_size
            integer, intent(in)   :: nobs            !< number of observations
            real(r4), intent(in)  :: hdxb(nbv,nobs)  !< ensemble perturbations in obs space
            real(r4), intent(in)  :: rdiag(nobs)     !< observation error variance
            real(r4), intent(in)  :: rloc(nobs)      !< observation localization weights
            real(r4), intent(in)  :: dep(nobs)       !< observation departures
            real(r4), intent(in)  :: infl            !< covariance inflation
            real(r4), intent(out) :: trans(nbv, nbv) !< ensemble transformation matrix
            ! temporary intermediate values
            real(r4)              :: hdxb_rinv(nbv,nobs)
            real(r4)              :: work1(nbv,nbv)
            ! real(r4)              :: work2(nbv,nobs)
            ! todo_letkf #1
            real(r4)              :: work2(nbv)
            real(r4)              :: work3(nbv)
            real(r4)              :: pa(nbv,nbv)
            real(r4)              :: eival(nbv)
            real(r4)              :: eivec(nbv,nbv)
            real(r4), allocatable :: evwork(:)
            integer               :: i, j, err
            real(r4)              :: r
            allocate(evwork((64+2)*nbv))
            ! hdxb rinv
            do j = 1, nobs
               hdxb_rinv(:,j) = hdxb(:,j) / rdiag(j) * rloc(j)
            end do
            ! hdxb^t rinv hdxb
            call sgemm('n','t', nbv, nbv, nobs, &
                 1.0_r4, hdxb_rinv, nbv, hdxb, nbv, 0.0_r4, work1, nbv)
            ! hdxb^t rinv hdxb + (k-1) i / rho (covariance inflation)
            r = real(nbv-1,r8)*1.0_r4/infl
            do i=1,nbv
               work1(i,i) = work1(i,i) + r
            end do
            ! eigenvalues and eigenvectors of above
            call ssyev('v','u', nbv, work1, nbv, eival, evwork, size(evwork), err)
            eivec = work1
            ! pa = [hdxb^t rinv hdxb + (m-1) i] inv
            do i=1,nbv
               work1(:,i) = eivec(:,i) / eival(i)
            end do
            call sgemm('n','t',nbv,nbv,nbv,1.0_r4, work1, nbv, eivec,&
                 nbv, 0.0_r4, pa, nbv)
!===========================================================================
!=================== complexity: nbv*nbv*nobs+nbv*nobs =====================
            ! ! pa hdxb_rinv^t
            ! call sgemm('n', 'n', nbv, nobs, nbv, 1.0_r4, pa, nbv, hdxb_rinv,&
            !      nbv, 0.0_r4, work2, nbv)
            ! ! pa hdxb_rinv^t dep
            ! work3 = 0
            ! do j=1,nobs
            !    work3 = work3 + work2(:,j) * dep(j)
            ! end do
!===========================================================================
!===================== complexity: nbv*nbv+nbv*nobs ======================== 
            ! todo_letkf #2
            ! hdxb_rinv^t dep
            work2 = 0
            do j=1,nobs
              work2 = work2 + hdxb_rinv(:,j) * dep(j)
            enddo
            ! pa hdxb_rinv^t dep
            work3 = 0
            do j=1,nbv
              work3 = work3 + pa(:,j) * work2(j)
            enddo           
!===========================================================================
            ! t = sqrt[(m-1)pa]
            do j = 1, nbv
               r = sqrt(real(nbv-1,r8) / eival(j))
               work1(:,j) = eivec(:,j) * r
            end do
            call sgemm('n', 't', nbv, nbv, nbv, 1.0_r4, work1, nbv, eivec, &
                 nbv, 0.0_r4, trans, nbv)
            ! relaxation (rtpp or rtps?)
            ! todo
            ! t + pa hdxb_rinv^t dep
            do j=1,nbv
               trans(:,j) = trans(:,j) + work3
            end do
            ! adaptive inflation
            ! todo
            deallocate(evwork)
        end subroutine letkf_core_solve_e0

        subroutine letkf_core_solve_d0(nbv, nobs, hdxb, rdiag, rloc, dep, infl, trans)
            integer, intent(in)   :: nbv             !< ens_size
            integer, intent(in)   :: nobs            !< number of observations
            real(r8), intent(in)  :: hdxb(nbv,nobs)  !< ensemble perturbations in obs space
            real(r8), intent(in)  :: rdiag(nobs)     !< observation error variance
            real(r8), intent(in)  :: rloc(nobs)      !< observation localization weights
            real(r8), intent(in)  :: dep(nobs)       !< observation departures
            real(r8), intent(in)  :: infl            !< covariance inflation
            real(r8), intent(out) :: trans(nbv, nbv) !< ensemble transformation matrix
            ! temporary intermediate values
            real(r8)              :: hdxb_rinv(nbv,nobs)
            real(r8)              :: work1(nbv,nbv)
            ! real(r8)              :: work2(nbv,nobs)
            ! todo_letkf #3
            real(r8)              :: work2(nbv)
            real(r8)              :: work3(nbv)
            real(r8)              :: pa(nbv,nbv)
            real(r8)              :: eival(nbv)
            real(r8)              :: eivec(nbv,nbv)
            real(r8), allocatable :: evwork(:)
            integer               :: i, j, err
            real(r8)              :: r
            allocate(evwork((64+2)*nbv))
            ! hdxb rinv
            do j = 1, nobs
               hdxb_rinv(:,j) = hdxb(:,j) / rdiag(j) * rloc(j)
            end do
            ! hdxb^t rinv hdxb
            call dgemm('n','t', nbv, nbv, nobs, &
                 1.0_r8, hdxb_rinv, nbv, hdxb, nbv, 0.0_r8, work1, nbv)
            ! hdxb^t rinv hdxb + (k-1) i / rho (covariance inflation)
            r = real(nbv-1,r8)*1.0_r8/infl
            do i=1,nbv
               work1(i,i) = work1(i,i) + r
            end do
            ! eigenvalues and eigenvectors of above
            call dsyev('v','u', nbv, work1, nbv, eival, evwork, size(evwork), err)
            eivec = work1
            ! pa = [hdxb^t rinv hdxb + (m-1) i] inv
            do i=1,nbv
               work1(:,i) = eivec(:,i) / eival(i)
            end do
            call dgemm('n','t',nbv,nbv,nbv,1.0_r8, work1, nbv, eivec,&
                 nbv, 0.0_r8, pa, nbv)
!===========================================================================
!=================== complexity: nbv*nbv*nobs+nbv*nobs =====================
            ! ! pa hdxb_rinv^t
            ! call dgemm('n', 'n', nbv, nobs, nbv, 1.0_r8, pa, nbv, hdxb_rinv,&
            !      nbv, 0.0_r8, work2, nbv)
            ! ! pa hdxb_rinv^t dep
            ! work3 = 0
            ! do j=1,nobs
            !    work3 = work3 + work2(:,j) * dep(j)
            ! end do
!===========================================================================
!===================== complexity: nbv*nbv+nbv*nobs ========================
            ! todo_letkf #4
            ! hdxb_rinv^t dep
            work2 = 0
            do j=1,nobs
              work2 = work2 + hdxb_rinv(:,j) * dep(j)
            enddo
            ! pa hdxb_rinv^t dep
            work3 = 0
            do j=1,nbv
              work3 = work3 + pa(:,j) * work2(j)
            enddo
!===========================================================================            
            ! t = sqrt[(m-1)pa]
            do j = 1, nbv
               r = sqrt(real(nbv-1,r8) / eival(j))
               work1(:,j) = eivec(:,j) * r
            end do
            call dgemm('n', 't', nbv, nbv, nbv, 1.0_r8, work1, nbv, eivec, &
                 nbv, 0.0_r8, trans, nbv)
            ! relaxation (rtpp or rtps?)
            ! todo
            ! t + pa hdxb_rinv^t dep
            do j=1,nbv
               trans(:,j) = trans(:,j) + work3
            end do
            ! adaptive inflation
            ! todo
            deallocate(evwork)
        end subroutine letkf_core_solve_d0

endmodule letkf
