module prepost
    use precision
    use io_control

    implicit none

    private

include './parameter.h'

    public      ::    readin  , writeout

    contains
        subroutine readin(bkg_dir,obs_dir,ens_size,numpatch,            &
                          nobs,olat,olon,hdxb,error,omb,qc_ens,omb_max, &
                          x_ens_atm, lonxy_atm, latxy_atm,              &
                          x_ens_lnd, grid2patch_start, grid2patch_count,&
                                     lonxy_lnd, latxy_lnd,              &
                                     itypwat,lb_patch                   & ! todo_colm_invar #1
                                     ) 
            character(*)         , intent(in)       ::      bkg_dir
            character(*)         , intent(in)       ::      obs_dir
            integer              , intent(in)       ::      ens_size
            integer              , intent(out)      ::      numpatch
            integer              , intent(inout)    ::      nobs
            real    , allocatable, intent(inout)    ::      olat         (:)
            real    , allocatable, intent(inout)    ::      olon         (:)
            real    , allocatable, intent(inout)    ::      hdxb         (:,:)
            real    , allocatable, intent(inout)    ::      error        (:)
            real    , allocatable, intent(inout)    ::      omb          (:)
            real    ,              intent(in)       ::      qc_ens 
            real    ,              intent(in)       ::      omb_max
            real(r4), allocatable, intent(inout)    ::      x_ens_atm(:,:,:,:)
            real(r4), allocatable, intent(inout)    ::      lonxy_atm(:,:)
            real(r4), allocatable, intent(inout)    ::      latxy_atm(:,:)
            real(r8), allocatable, intent(inout)    ::      x_ens_lnd(:,:,:)
            integer , allocatable, intent(inout)    ::      grid2patch_start(:,:)
            integer , allocatable, intent(inout)    ::      grid2patch_count(:,:)
            real(r4), allocatable, intent(inout)    ::      lonxy_lnd(:,:)
            real(r4), allocatable, intent(inout)    ::      latxy_lnd(:,:)

            integer , allocatable, intent(inout)    ::      itypwat(:,:) ! todo_colm_invar #2
            integer , allocatable, intent(inout)    ::      lb_patch(:,:)
            !===================================================================
            !                     MODEL STATE
            !===================================================================
            !----------------------- grapes ---------------------------
            real(r4), allocatable        ::      head(:)
            real(r4), allocatable        ::      pi (:, :, :)
            real(r4), allocatable        ::      u  (:, :, :)
            real(r4), allocatable        ::      v  (:, :, :)
            real(r4), allocatable        ::      w  (:, :, :)
            real(r4), allocatable        ::      th (:, :, :)
            real(r4), allocatable        ::      qv (:, :, :)
            real(r4), allocatable        ::      ts (:, :, :)
            real(r4), allocatable        ::      t2 (:, :, :)
            real(r4), allocatable        ::      q2 (:, :, :)
            real(r4), allocatable        ::      u10(:, :, :)
            real(r4), allocatable        ::      v10(:, :, :)
            logical                      ::      first_read_grapes = .true.
            !--------------------------- colm ---------------------------
            integer                      ::      idate(3)
            integer , allocatable        ::      ixy_patch(:)
            integer , allocatable        ::      jxy_patch(:)
            integer , allocatable        ::      mxy_patch(:)
            real(r8), allocatable        ::      wtxy_patch(:)
            real(r8), allocatable        ::      ftune(:)
            real(r8), allocatable        ::      fcon(:,:)
            real(r8), allocatable        ::      fvar(:,:)
            real(r8), allocatable        ::      oro(:)
            real(r8), allocatable        ::      topo(:,:)
            real(r8), allocatable        ::      mask(:,:)
            logical                      ::      first_read_colm = .true.
            !===================================================================
            !                     OBS STATE
            !===================================================================
            integer                      ::      nobs_raw
            real    , allocatable        ::      rlat         (:)
            real    , allocatable        ::      rlon         (:)
            real    , allocatable        ::      tbb          (:)
            integer , allocatable        ::      qc_flag      (:,:)
            real    , allocatable        ::      oberr        (:) 
            real    , allocatable        ::      tbb_bmo      (:,:)
            logical                      ::      first_read_obs  = .true.
            !---------------------- temp use ------------------------
            integer                      ::    ens, num

            print *, '-----------------'
            print *, 'read in start'
            print *, ''
            num = int(log10(real(ens_size)))+1
            do ens=1,ens_size
                call obs_read(obs_dir, ens_size, ens, num,        &
                              nobs_raw, rlat, rlon, tbb, qc_flag, &
                              oberr, tbb_bmo, first_read_obs      )
                first_read_obs    = .false.
            enddo
            call obs_pre(ens_size,nobs_raw,rlat,rlon,tbb,qc_flag,oberr,tbb_bmo, &
                         nobs,olat,olon,hdxb,error,omb,qc_ens,omb_max)
            call obs_deallocate(rlat, rlon, tbb, qc_flag, oberr, tbb_bmo)
            if(nobs == 0) return
            do ens=1,ens_size
                call grapes_input(bkg_dir,ens,num,head,    &
                                  pi, u, v, w, th, qv,     &
                                  ts, t2, q2, u10, v10, first_read_grapes)
                if(first_read_grapes) then
                    call grapes_grid_init(lonxy_atm,latxy_atm)
                endif
                call grapes2x_ens(x_ens_atm, ens, ens_size, first_read_grapes,  &
                                  th,qv)      ! todo_grapes #1 th + qv
                first_read_grapes = .false.
                call colm_input(bkg_dir,ens,num,numpatch,idate,          &
                                ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                                ftune,fcon,fvar,oro,topo,mask,first_read_colm )
                if(first_read_colm) then
                    call colm_grid2patch_init(numpatch,ixy_patch,jxy_patch,    &
                                             grid2patch_start, grid2patch_count)
                    call colm_grid_init(numpatch,fcon,grid2patch_start,lonxy_lnd,latxy_lnd)
                    call colm_invar_read(numpatch,fcon, &
                                        itypwat) ! todo_colm_invar #3
                endif 
                call colm2x_ens(numpatch,x_ens_lnd,fvar,first_read_colm, &
                                ens,ens_size,itypwat,lb_patch  ) ! todo_colm_var #1
                first_read_colm   = .false.
            enddo
            print *, ''
            print *, 'read in over'
            print *, '-----------------'
            call grapes_deallocate(head, pi, u, v, w, th, qv, ts, t2, q2, u10, v10)
            call colm_deallocate(ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                                ftune,fcon,fvar,oro,topo,mask)
        endsubroutine readin

        subroutine writeout(nobs,ana_dir,ens_size,x_ens_atm,x_ens_lnd,    &
                            itypwat,lb_patch)
            integer              , intent(in)       ::      nobs
            character(*)         , intent(in)       ::      ana_dir
            integer              , intent(in)       ::      ens_size
            real(r4), allocatable, intent(inout)    ::      x_ens_atm(:,:,:,:)
            real(r8), allocatable, intent(inout)    ::      x_ens_lnd(:,:,:)
            integer , allocatable, intent(inout)    ::      itypwat(:,:) ! todo_colm_invar #4
            integer , allocatable, intent(inout)    ::      lb_patch(:,:)
            !===================================================================
            !                     MODEL STATE
            !===================================================================
            !----------------------- grapes ---------------------------
            real(r4), allocatable        ::      head(:)
            real(r4), allocatable        ::      pi (:, :, :)
            real(r4), allocatable        ::      u  (:, :, :)
            real(r4), allocatable        ::      v  (:, :, :)
            real(r4), allocatable        ::      w  (:, :, :)
            real(r4), allocatable        ::      th (:, :, :)
            real(r4), allocatable        ::      qv (:, :, :)
            real(r4), allocatable        ::      ts (:, :, :)
            real(r4), allocatable        ::      t2 (:, :, :)
            real(r4), allocatable        ::      q2 (:, :, :)
            real(r4), allocatable        ::      u10(:, :, :)
            real(r4), allocatable        ::      v10(:, :, :)
            logical                      ::      first_read_grapes = .true.
            !----------------------- colm ---------------------------
            integer                      ::      numpatch
            integer                      ::      idate(3)
            integer , allocatable        ::      ixy_patch(:)
            integer , allocatable        ::      jxy_patch(:)
            integer , allocatable        ::      mxy_patch(:)
            real(r8), allocatable        ::      wtxy_patch(:)
            real(r8), allocatable        ::      ftune(:)
            real(r8), allocatable        ::      fcon(:,:)
            real(r8), allocatable        ::      fvar(:,:)
            real(r8), allocatable        ::      oro(:)
            real(r8), allocatable        ::      topo(:,:)
            real(r8), allocatable        ::      mask(:,:)
            logical                      ::      first_read_colm = .true.
            !---------------------- temp use ------------------------
            integer                      ::    ens, num

            if(nobs == 0) return

            print *, '-----------------'
            print *, 'write out start'
            print *, ''
            num = int(log10(real(ens_size)))+1
            do ens=1,ens_size
                call grapes_input (ana_dir,ens,num,head,     &
                                   pi, u, v, w, th, qv,      &
                                   ts, t2, q2, u10, v10, first_read_grapes)
                call x_ens2grapes(x_ens_atm, ens, ens_size,  &
                                th,qv)         ! todo_grapes #2 th + qv
                first_read_grapes = .false.
                call grapes_output(ana_dir,ens,num,head,     &
                                   pi, u, v, w, th, qv,      &
                                   ts, t2, q2, u10, v10                   )
                call colm_input(ana_dir,ens,num,numpatch,idate,          &
                                ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                                ftune,fcon,fvar,oro,topo,mask,first_read_colm )
                call x_ens2colm(numpatch,x_ens_lnd,fvar,first_read_colm, &
                                ens,ens_size,itypwat,lb_patch  ) ! todo_colm_var #2
                first_read_colm   = .false.
                call colm_output(ana_dir,ens,num,numpatch,idate,         &
                                ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                                ftune,fcon,fvar,oro,topo,mask                 )
            enddo
            print *, ''
            print *, 'write out over'
            print *, '-----------------'
            call grapes_deallocate(head, pi, u, v, w, th, qv, ts, t2, q2, u10, v10)
            call colm_deallocate(ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                                ftune,fcon,fvar,oro,topo,mask)
        endsubroutine writeout 

        subroutine grapes_grid_init(lonxy,latxy)
            real(r4), allocatable, intent(inout)  ::   lonxy(:,:)
            real(r4), allocatable, intent(inout)  ::   latxy(:,:)

            real(r4)                              ::   dx, dy
            integer                               ::   i, j

            allocate(latxy(lon_points_atm,lat_points_atm))
            allocate(lonxy(lon_points_atm,lat_points_atm))

            dx = (edgee_atm-edgew_atm)/(lon_points_atm-1)
            dy = (edgen_atm-edges_atm)/(lat_points_atm-1)

            do i=1,lon_points_atm
                lonxy(i,:) = edgew_atm + (i-1) * dx
            enddo
            
            do j=1,lat_points_atm
                latxy(:,j) = edges_atm + (j-1) * dy
            enddo

            ! do i=1,lon_points_atm
                ! do j=1,lat_points_atm
                    ! latxy(i,j) = edges_atm + (j-1) * dy
                    ! lonxy(i,j) = edgew_atm + (i-1) * dx
                ! enddo
            ! enddo
            
        endsubroutine grapes_grid_init

        subroutine grapes2x_ens(x_ens_atm, ens, ens_size, first_read_grapes,  &
                                v1,v2)    ! todo_grapes #3
            real(r4), allocatable, intent(inout)      ::      x_ens_atm(:,:,:,:)
            integer              , intent(in)         ::      ens_size
            integer              , intent(in)         ::      ens
            logical              , intent(in)         ::      first_read_grapes
            ! todo_grapes #4
            real(r4)             , intent(in)         ::      v1(lon_points_atm,lat_points_atm,level_atm)
            real(r4)             , intent(in)         ::      v2(lon_points_atm,lat_points_atm,level_atm)

            integer                                   ::      l, level

            if(first_read_grapes) then
                allocate(x_ens_atm(ens_size,lon_points_atm,lat_points_atm,nvar_atm))
            endif

            l = 0

            level = size(v1,dim=3)
            l = l + level
            x_ens_atm(ens,:,:,l-level+1:l) = v1(:,:,:)

            level = size(v2,dim=3)
            l = l + level
            x_ens_atm(ens,:,:,l-level+1:l) = v2(:,:,:)

            ! todo_grapes #5 (can be added more, v3 v4 v5...)

            if(l /= nvar_atm) then
                print *, 'wrong level setting for atm: '
                print *, 'nvar = ', nvar_atm, ', due nvar = ', l
                stop 20
            endif

        endsubroutine grapes2x_ens

        subroutine x_ens2grapes(x_ens_atm, ens, ens_size ,&
                                v1,v2)    ! todo_grapes #6
            real(r4), allocatable, intent(inout)      ::      x_ens_atm(:,:,:,:)
            integer              , intent(in)         ::      ens
            integer              , intent(in)         ::      ens_size
            ! todo_grapes #7
            real(r4)             , intent(inout)      ::      v1(lon_points_atm,lat_points_atm,level_atm)
            real(r4)             , intent(inout)      ::      v2(lon_points_atm,lat_points_atm,level_atm)

            integer                                   ::      l, level

            l = 0

            level = size(v1,dim=3)
            l = l + level
            v1(:,:,:) = x_ens_atm(ens,:,:,l-level+1:l)

            level = size(v2,dim=3)
            l = l + level
            v2(:,:,:) = x_ens_atm(ens,:,:,l-level+1:l)

            ! todo_grapes #8 (can be added more, v3 v4 v5...)

            if(l /= nvar_atm) then
                print *, 'wrong level setting for atm: '
                print *, 'nvar = ', nvar_atm, ', due nvar = ', l
                stop 20
            endif

            if(ens == ens_size) then
                deallocate(x_ens_atm)
            endif

        endsubroutine x_ens2grapes

        subroutine colm_grid2patch_init(numpatch,ixy_patch,jxy_patch,    &
                                        grid2patch_start, grid2patch_count)
            integer              , intent(in)       ::      numpatch
            integer              , intent(in)       ::      ixy_patch(numpatch)
            integer              , intent(in)       ::      jxy_patch(numpatch)
            integer , allocatable, intent(inout)    ::      grid2patch_start(:,:)
            integer , allocatable, intent(inout)    ::      grid2patch_count(:,:)
            integer                                 ::      i, j, k, i_save, j_save
            logical                                 ::      first_save

            allocate(grid2patch_start(lon_points_lnd,lat_points_lnd))
            allocate(grid2patch_count(lon_points_lnd,lat_points_lnd))

            grid2patch_count = 0
            first_save       = .true.
            do k=1, numpatch
              i = ixy_patch(k)
              j = jxy_patch(k)
              grid2patch_count(i,j) = grid2patch_count(i,j) + 1
              if(first_save .or. i/=i_save .or. j/=j_save) then
                    first_save = .false.
                    i_save = i
                    j_save = j
                    grid2patch_start(i,j) = k
                endif
            enddo
            
        endsubroutine colm_grid2patch_init

        subroutine colm_grid_init(numpatch,fcon,grid2patch_start,lonxy_lnd,latxy_lnd)
            integer              , intent(in)       ::      numpatch
            real(r8)             , intent(in)       ::      fcon(numpatch,nfcon)
            integer              , intent(in)       ::      grid2patch_start(lon_points_lnd,lat_points_lnd)
            real(r4), allocatable, intent(inout)    ::      lonxy_lnd(:,:)
            real(r4), allocatable, intent(inout)    ::      latxy_lnd(:,:)
            integer                                 ::      i, j, k, i_save, j_save
            real(r8), parameter                     ::      pi180 = 180.0_r8/(4.0_r8*atan(1.0_r8))

            allocate(lonxy_lnd(lon_points_lnd,lat_points_lnd))
            allocate(latxy_lnd(lon_points_lnd,lat_points_lnd))

            do j=1,lat_points_lnd
              do i=1,lon_points_lnd
                k=grid2patch_start(i,j)
                lonxy_lnd(i,j) = real(fcon(k,lon_id_lnd)*pi180,kind=r4)
                latxy_lnd(i,j) = real(fcon(k,lat_id_lnd)*pi180,kind=r4)
              enddo
            enddo

        endsubroutine colm_grid_init

        subroutine colm_invar_read(numpatch, fcon, &
                                  v1) ! todo_colm_invar #5
            integer              , intent(in)       ::      numpatch
            real(r8)             , intent(in)       ::      fcon(numpatch,nfcon)
            ! todo_colm_invar #6
            integer , allocatable, intent(inout)    ::      v1  (:,:)
            integer                                 ::      cnt = 0

            cnt = cnt + 1 
            allocate(v1(numpatch,invar_len(cnt)))
            v1(:,:) = int(fcon(:,invar_idx_s(cnt):invar_idx_s(cnt)+invar_len(cnt)-1))

            ! todo_colm_invar #7 (can be added more, v2 v3 v4...)

            if(cnt /= rd_invar_cnt) then
                print *, 'wrong input setting for colm_invar: '
                print *, 'num = ', rd_invar_cnt, ', due num = ', cnt
                stop 22
            endif

        endsubroutine colm_invar_read

        subroutine colm2x_ens(numpatch,x_ens_lnd,fvar,first_read_colm, &
                              ens,ens_size,itypwat,lb_patch  ) ! todo_colm_var #3
            integer              , intent(in)      ::      numpatch
            real(r8), allocatable, intent(inout)   ::      x_ens_lnd(:,:,:)
            real(r8)             , intent(in)      ::      fvar(numpatch,nfvar)
            logical              , intent(in)      ::      first_read_colm
            integer              , intent(in)      ::      ens
            integer              , intent(in)      ::      ens_size
            integer                                ::      cnt, l, np, nl
            integer              , intent(in)      ::      itypwat(numpatch,1) ! todo_colm_var #4
            integer, allocatable , intent(inout)   ::      lb_patch(:,:)! todo_colm_var #
            real(r8)                               ::      lnd_buf(numpatch,sum(var_len))

            if(first_read_colm) then
                allocate(x_ens_lnd(ens_size,numpatch,nvar_lnd))
                allocate(lb_patch(ens_size,numpatch))
            endif

            l = 0
            do cnt=1,nvar_lnd_raw
                l = l + var_len(cnt)
                lnd_buf(:,l-var_len(cnt)+1:l) = fvar(:,var_idx_s(cnt):var_idx_s(cnt)+var_len(cnt)-1)
            enddo
            
            lb_patch(ens,:) = 5
            do np=1,numpatch
                if(itypwat(np,1) > itypwat_max) cycle
                do nl=1,5
                    if(lnd_buf(np,15+nl)+lnd_buf(np,30+nl) > 0.) lb_patch(ens,np) = lb_patch(ens,np) - 1
                enddo
                ! todo_speedy #1
                !if(lnd_buf(np,0+lb_patch(ens,np)+1) /= lnd_buf(np,46)) then
                !    print *, 'patch: ', np, 'itypwat: ', itypwat(np,1)
                !    print *, "t_soil  at  : ", 0+lb_patch(ens,np)+1," :", lnd_buf(np,0+lb_patch(ens,np)+1)
                !    print *, "t_grnd      : ", lnd_buf(np,46)
                !    print *, "t_soil at whole: ", lnd_buf(np,1:15)
                !    stop 23
                !endif
                cnt = 0
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,0+lb_patch(ens,np)+1)    ! todo_colm_var #5 tss[1]
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,0+lb_patch(ens,np)+1+1)  ! todo_colm_var #6 tss[2]
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,15+lb_patch(ens,np)+1)   ! todo_colm_var #7 wliq[1]
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,15+lb_patch(ens,np)+1+1) ! todo_colm_var #8 wliq[2]
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,30+lb_patch(ens,np)+1)   ! todo_colm_var #9 wice[1]
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,30+lb_patch(ens,np)+1+1) ! todo_colm_var #10 wice[2]
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,46)                  ! todo_colm_var #11 tg
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,47)                  ! todo_colm_var #12 tlsun
                cnt = cnt + 1
                x_ens_lnd(ens,np,cnt) = lnd_buf(np,48)                  ! todo_colm_var #13 tlsha
                ! todo_colm_var #14 (can be added more)
            enddo

            if(cnt /= nvar_lnd) then
                print *, 'wrong input setting for colm_var_da: '
                print *, 'num = ', nvar_lnd, ', due num = ', cnt
                stop 24
            endif

        endsubroutine colm2x_ens

        subroutine x_ens2colm(numpatch,x_ens_lnd,fvar,first_read_colm, &
                              ens,ens_size,itypwat,lb_patch            ) ! todo_colm_var #15
            integer              , intent(in)      ::      numpatch
            real(r8), allocatable, intent(inout)   ::      x_ens_lnd(:,:,:)
            real(r8)             , intent(inout)   ::      fvar(numpatch,nfvar)
            logical              , intent(in)      ::      first_read_colm
            integer              , intent(in)      ::      ens
            integer              , intent(in)      ::      ens_size
            integer              , intent(in)      ::      itypwat(numpatch,1) ! todo_colm_var #16
            integer              , intent(in)      ::      lb_patch(ens_size,numpatch)
            integer                                ::      cnt, l, np, nl
            ! todo_colm_var #17
            integer , parameter                    ::      lenth = sum(var_len)+1  ! 1 for scv
            real(r8)                               ::      lnd_buf(numpatch,lenth)
            real(r8)                               ::      scv_check

            l = 0
            do cnt=1,nvar_lnd_raw
                l = l + var_len(cnt)
                lnd_buf(:,l-var_len(cnt)+1:l) = fvar(:,var_idx_s(cnt):var_idx_s(cnt)+var_len(cnt)-1)
            enddo
            lnd_buf(:,lenth) = fvar(:,scv_idx)
            
            do np=1,numpatch
                if(itypwat(np,1) > itypwat_max) cycle
                !do nl=1,5
                !    if(lnd_buf(np,15+nl)+lnd_buf(np,30+nl) > 0.) lb_patch(ens,np) = lb_patch(ens,np) - 1
                !enddo
                !! todo_speedy #2
                !if(lnd_buf(np,0+lb_patch(ens,np)+1) /= lnd_buf(np,46)) then
                !    print *, 'patch: ', np, 'itypwat: ', itypwat(np,1)
                !    print *, "t_soil  at  : ", 0+lb_patch(ens,np)+1," :", lnd_buf(np,0+lb_patch(ens,np)+1)
                !    print *, "t_grnd      : ", lnd_buf(np,46)
                !    print *, "t_soil at whole: ", lnd_buf(np,1:15)
                !    stop 23
                !endif
                cnt = 0
                cnt = cnt + 1
                lnd_buf(np,0+lb_patch(ens,np)+1)    = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #17 tss[1]
                cnt = cnt + 1
                lnd_buf(np,0+lb_patch(ens,np)+1+1)  = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #18 tss[2]
                cnt = cnt + 1
                lnd_buf(np,15+lb_patch(ens,np)+1)   = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #19 wliq[1]
                cnt = cnt + 1
                lnd_buf(np,15+lb_patch(ens,np)+1+1) = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #20 wliq[2]
                cnt = cnt + 1
                lnd_buf(np,30+lb_patch(ens,np)+1)   = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #21 wice[1]
                cnt = cnt + 1
                lnd_buf(np,30+lb_patch(ens,np)+1+1) = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #22 wice[2]
                cnt = cnt + 1
                lnd_buf(np,46)                  = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #23 tg
                cnt = cnt + 1
                lnd_buf(np,47)                  = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #24 tlsun
                cnt = cnt + 1
                lnd_buf(np,48)                  = x_ens_lnd(ens,np,cnt)   ! todo_colm_var #25 tlsha
                ! todo_colm_var #26 (can be added more)
            enddo 

            if(cnt /= nvar_lnd) then
                print *, 'wrong input setting for colm_var_da: '
                print *, 'num = ', nvar_lnd, ', due num = ', cnt
                stop 24
            endif

            do np=1,numpatch
                if(itypwat(np,1) > itypwat_max) cycle
                scv_check = 0.0_r8
                do nl=1,5
                    scv_check = scv_check + lnd_buf(np,15+nl) + lnd_buf(np,30+nl)
                enddo
                lnd_buf(np,lenth) = scv_check
            enddo

            l = 0
            do cnt=1,nvar_lnd_raw
                l = l + var_len(cnt)
                fvar(:,var_idx_s(cnt):var_idx_s(cnt)+var_len(cnt)-1) = lnd_buf(:,l-var_len(cnt)+1:l)
            enddo

            fvar(:,scv_idx) = lnd_buf(:,lenth)

        endsubroutine x_ens2colm

        subroutine obs_pre(ens_size,nobs_raw,rlat,rlon,tbb,qc_flag,oberr,tbb_bmo, &
                           nobs,olat,olon,hdxb,error,omb,qc_ens,omb_max           )
            integer               , intent(in)           ::    ens_size
            integer               , intent(in)           ::    nobs_raw
            real                  , intent(in)           ::    rlat         (nobs_raw)
            real                  , intent(in)           ::    rlon         (nobs_raw)
            real                  , intent(in)           ::    tbb          (nobs_raw)
            integer               , intent(in)           ::    qc_flag      (nobs_raw,ens_size)
            real                  , intent(in)           ::    oberr        (nobs_raw) 
            real                  , intent(in)           ::    tbb_bmo      (nobs_raw,ens_size)
            integer               , intent(inout)        ::    nobs
            real    , allocatable , intent(inout)        ::    olat         (:)
            real    , allocatable , intent(inout)        ::    olon         (:)
            real    , allocatable , intent(inout)        ::    hdxb         (:,:)
            real    , allocatable , intent(inout)        ::    error        (:)
            real    , allocatable , intent(inout)        ::    omb          (:)
            real    ,               intent(in)           ::    qc_ens
            real    ,               intent(in)           ::    omb_max

            real                                         ::    hxb_raw      (nobs_raw,ens_size)
            logical                                      ::    qc_raw       (nobs_raw,ens_size)
            real                                         ::    hxb_mean     (nobs_raw)
            logical                                      ::    qc_list      (nobs_raw)
            logical                                      ::    omb_list     (nobs_raw)
            real                                         ::    hxb1         (ens_size)
            real                                         ::    hxb2         (ens_size)
            integer                                      ::    s, cnt, ens_min

            ens_min = int(ens_size*qc_ens)
            print *, 'ens_size =', ens_size
            print *, 'at least ',ens_min, ' must have qc=0'
            if(ens_min <= 0) then
                stop 25
            endif

            do s=1,nobs_raw
                hxb_raw(s,1:ens_size) = tbb(s)+tbb_bmo(s,1:ens_size)
            enddo
            qc_raw  = qc_flag == 0.
            hxb_mean(1:nobs_raw) = sum(hxb_raw,dim=2,mask=qc_raw)/real(count(mask=qc_raw,dim=2),kind=r4)
            omb_list(1:nobs_raw) = abs(hxb_mean(1:nobs_raw)-tbb(1:nobs_raw)) <= omb_max
            qc_list (1:nobs_raw) = count(mask=qc_raw,dim=2) >= ens_min
            nobs = count(omb_list .and. qc_list)
            allocate(olat(nobs))
            allocate(olon(nobs))
            allocate(hdxb(nobs,ens_size))
            allocate(error(nobs))
            allocate(omb (nobs))
            cnt = 0
            hxb2= 0.0_r4
            do s=1,nobs_raw
                if(.not. omb_list(s) .or. .not. qc_list(s)) cycle
                cnt = cnt + 1
                olat(cnt) = rlat(s)
                olon(cnt) = rlon(s)
                error(cnt)= oberr(s)
                omb (cnt) = tbb(s)-hxb_mean(s)
                hxb1(1:ens_size)=hxb_raw(s,1:ens_size) - hxb_mean(s)
                hdxb(cnt,1:ens_size) = merge(hxb1(1:ens_size),hxb2(1:ens_size),qc_raw(s,:))
            enddo
            if(cnt/=nobs) then
                print *, 'cnt = ', cnt
                print *, 'nobs= ', nobs
                stop 300
            endif

            ! allocate(olat(nobs_raw))
            ! allocate(olon(nobs_raw))
            ! allocate(hdxb(nobs_raw,ens_size))
            ! allocate(error(nobs_raw))
            ! allocate(omb (nobs_raw))
            ! do s=1,nobs_raw
            !     hxb1(:) = tbb (s) + tbb_bmo(s,:)                  ! cal raw ens_size*H(x) 
            !     cnt = 0                                           ! valid ens_size with qc = 0
            !     do ens=1,ens_size                                      
            !         if(qc_flag(s,ens) /= 0) cycle                     
            !         cnt      = cnt + 1                                
            !         hxb2(cnt) = tbb (s) + tbb_bmo(s,ens)          ! cal cnt*H(x) satisfying qc == 0
            !     enddo
            !     hxb_mean = sum(hxb2(1:cnt))/real(ens_size,kind=r4)
            !     ! todo_obs # replace bad hxb with mean(hxb with good quality)
            !     do ens=1,ens_size
            !         if(qc_flag(s,ens) /= 0) hxb1(ens) = hxb_mean
            !     enddo
            !     if(cnt < ens_min .or. &                           ! valid ens < ens_min
            !        abs(tbb(s)-hxb_mean) > omb_max) cycle          ! abs(omb) > omb_max
            !     nobs = nobs + 1
            !     olat (nobs)   =   rlat(s)                         ! put raw lat into 
            !     olon (nobs)   =   rlon(s)                         ! put raw lon into
            !     hdxb (nobs,:) =   hxb1(:) - hxb_mean              ! raw ens_size*H(x) - mean(cnt*H(x))
            !     error(nobs)   =   oberr(s)                        ! put raw err into
            !     omb  (nobs)   =   tbb(s)  - hxb_mean              ! o - mean(cnt*H(x))
            ! enddo

        endsubroutine

endmodule prepost
