module io_control

    use precision
    use netcdf

    implicit none

    private 
include './parameter.h'
    public     ::     grapes_input, grapes_output, grapes_deallocate, &
                      colm_input  , colm_output, colm_deallocate    , &
                      obs_read, obs_deallocate

    contains

        subroutine grapes_input(bkg_dir,ens,num,head,&
                            pi, u, v, w, th, qv,     &
                            ts, t2, q2, u10, v10, first_read_grapes)
            character(*)         ,  intent(in)    ::      bkg_dir
            integer              ,  intent(in)    ::      ens
            integer              ,  intent(in)    ::      num
            real(r4), allocatable, intent(inout)  ::      head(:)
            real(r4), allocatable, intent(inout)  ::      pi (:, :, :)
            real(r4), allocatable, intent(inout)  ::      u  (:, :, :)
            real(r4), allocatable, intent(inout)  ::      v  (:, :, :)
            real(r4), allocatable, intent(inout)  ::      w  (:, :, :)
            real(r4), allocatable, intent(inout)  ::      th (:, :, :)
            real(r4), allocatable, intent(inout)  ::      qv (:, :, :)
            real(r4), allocatable, intent(inout)  ::      ts (:, :, :)
            real(r4), allocatable, intent(inout)  ::      t2 (:, :, :)
            real(r4), allocatable, intent(inout)  ::      q2 (:, :, :)
            real(r4), allocatable, intent(inout)  ::      u10(:, :, :)
            real(r4), allocatable, intent(inout)  ::      v10(:, :, :)
            logical ,              intent(in)     ::      first_read_grapes

            integer     ,  parameter     ::    fu = 200
            character(len=1024)          ::    fn
            character(len=256)           ::    fmt
            character(len=20), parameter ::    fi = 'grapesinput'
            
            write(fmt, '(A,I0,A)') '(A,A,I0.',num,',A,A)' 
            write(fn,fmt) trim(bkg_dir),'/', ens,'/', fi
            open(fu, file=trim(fn), form="unformatted", access="stream", status="old")
            print*, 'read in '//trim(fn)
            if(first_read_grapes) then
                allocate(head(head_len))
                allocate(pi (lon_points_atm,lat_points_atm,level_atm))
                allocate( u (lon_points_atm,lat_points_atm,level_atm))
                allocate( v (lon_points_atm,lat_points_atm,level_atm))
                allocate( w (lon_points_atm,lat_points_atm,level_atm))
                allocate(th (lon_points_atm,lat_points_atm,level_atm))
                allocate(qv (lon_points_atm,lat_points_atm,level_atm))
                allocate( ts(lon_points_atm,lat_points_atm,1))
                allocate( t2(lon_points_atm,lat_points_atm,1))
                allocate( q2(lon_points_atm,lat_points_atm,1))
                allocate(u10(lon_points_atm,lat_points_atm,1))
                allocate(v10(lon_points_atm,lat_points_atm,1))
            endif
            read(fu) head
            read(fu) pi
            read(fu) u
            read(fu) v
            read(fu) w
            read(fu) th
            read(fu) qv
            read(fu) ts
            read(fu) t2
            read(fu) q2
            read(fu) u10
            read(fu) v10
            close(fu)

        endsubroutine grapes_input

        subroutine grapes_output(ana_dir,ens,num,head,&
                            pi, u, v, w, th, qv,      &
                            ts, t2, q2, u10, v10      )
            character(*)         ,  intent(in)    ::      ana_dir
            integer              ,  intent(in)    ::      ens
            integer              ,  intent(in)    ::      num
            real(r4), allocatable, intent(inout)  ::      head(:)
            real(r4), allocatable, intent(inout)  ::      pi (:, :, :)
            real(r4), allocatable, intent(inout)  ::      u  (:, :, :)
            real(r4), allocatable, intent(inout)  ::      v  (:, :, :)
            real(r4), allocatable, intent(inout)  ::      w  (:, :, :)
            real(r4), allocatable, intent(inout)  ::      th (:, :, :)
            real(r4), allocatable, intent(inout)  ::      qv (:, :, :)
            real(r4), allocatable, intent(inout)  ::      ts (:, :, :)
            real(r4), allocatable, intent(inout)  ::      t2 (:, :, :)
            real(r4), allocatable, intent(inout)  ::      q2 (:, :, :)
            real(r4), allocatable, intent(inout)  ::      u10(:, :, :)
            real(r4), allocatable, intent(inout)  ::      v10(:, :, :)

            integer     ,  parameter     ::    fu = 201
            character(len=1024)          ::    fn
            character(len=256)           ::    fmt
            character(len=20), parameter ::    fi = 'grapesinput'
            
            write(fmt, '(a,i0,a)') '(a,a,i0.',num,',a,a)' 
            write(fn,fmt) trim(ana_dir),'/', ens,'/', fi

            open(fu, file=trim(fn), form="unformatted", access="stream", status="replace", action='write')
            print*, 'write out '//trim(fn)
            write(fu) head
            write(fu) pi
            write(fu) u
            write(fu) v
            write(fu) w
            write(fu) th
            write(fu) qv
            write(fu) ts
            write(fu) t2
            write(fu) q2
            write(fu) u10
            write(fu) v10
            close(fu)
            
        endsubroutine grapes_output

        subroutine grapes_deallocate(head, pi, u, v, w, th, qv, ts, t2, q2, u10, v10)

            real(r4), allocatable, intent(out)    ::      head(:)
            real(r4), allocatable, intent(out)    ::      pi (:, :, :)
            real(r4), allocatable, intent(out)    ::      u  (:, :, :)
            real(r4), allocatable, intent(out)    ::      v  (:, :, :)
            real(r4), allocatable, intent(out)    ::      w  (:, :, :)
            real(r4), allocatable, intent(out)    ::      th (:, :, :)
            real(r4), allocatable, intent(out)    ::      qv (:, :, :)
            real(r4), allocatable, intent(out)    ::      ts (:, :, :)
            real(r4), allocatable, intent(out)    ::      t2 (:, :, :)
            real(r4), allocatable, intent(out)    ::      q2 (:, :, :)
            real(r4), allocatable, intent(out)    ::      u10(:, :, :)
            real(r4), allocatable, intent(out)    ::      v10(:, :, :)

            deallocate(head)
            deallocate( pi)
            deallocate(  u)
            deallocate(  v)
            deallocate(  w)
            deallocate( th)
            deallocate( qv)
            deallocate( ts)
            deallocate( t2)
            deallocate( q2)
            deallocate(u10)
            deallocate(v10)
            ! print *, '* deallocate var memory over (grapes)'

        endsubroutine grapes_deallocate

        subroutine colm_input(bkg_dir,ens,num,numpatch,idate,          &
                              ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                              ftune,fcon,fvar,oro,topo,mask,first_read_colm )
            character(*),  intent(in)    ::    bkg_dir
            integer     ,  intent(in)    ::    ens
            integer     ,  intent(in)    ::    num

            integer     ,  intent(inout)   ::    numpatch
            integer     ,  intent(inout)   ::    idate(3)
            integer,  allocatable, intent(inout)  ::   ixy_patch(:)
            integer,  allocatable, intent(inout)  ::   jxy_patch(:)
            integer,  allocatable, intent(inout)  ::   mxy_patch(:)
            real(r8), allocatable, intent(inout)  ::   wtxy_patch(:)
            real(r8), allocatable, intent(inout)  ::   ftune(:)
            real(r8), allocatable, intent(inout)  ::   fcon(:,:)
            real(r8), allocatable, intent(inout)  ::   fvar(:,:)
            real(r8), allocatable, intent(inout)  ::   oro(:)
            real(r8), allocatable, intent(inout)  ::   topo(:,:)
            real(r8), allocatable, intent(inout)  ::   mask(:,:)
            logical              , intent(in)     ::   first_read_colm           

            character(len=1024)          ::    fn
            character(len=256)           ::    fmt
            character(len=20), parameter ::    fi = 'colminput'
            integer     ,  parameter     ::    fu = 100
            integer                      ::    n

            integer                      ::    storage_test
            real(r8)    ,  parameter     ::    pi  = 4.*atan(1.0_r8)
            
            write(fmt, '(A,I0,A)') '(A,A,I0.',num,',A,A)' 
            write(fn,fmt) trim(bkg_dir),'/', ens,'/', fi
            open(fu,file=trim(fn),form='unformatted',convert="little_endian",access='sequential',status='old',action='read')
            print*, 'read in '//trim(fn)
            read(fu) numpatch
            ! print *, 'numpatch=',numpatch
            if(first_read_colm) then
                allocate(ixy_patch        (numpatch))
                allocate(jxy_patch        (numpatch))
                allocate(mxy_patch        (numpatch))
                allocate(wtxy_patch       (numpatch))
                allocate(ftune            (nftune))
                allocate(fcon             (numpatch,nfcon))
                allocate(fvar             (numpatch,nfvar))
                allocate(oro              (numpatch))
                allocate(topo             (lon_points_lnd,lat_points_lnd))
                allocate(mask             (lon_points_lnd,lat_points_lnd))
            endif
            read(fu) idate    
            read(fu) ixy_patch  
            read(fu) jxy_patch  
            read(fu) mxy_patch  
            read(fu) wtxy_patch 
            read(fu) ftune
            do n=1,nfcon
                read(fu) fcon(:,n)
            enddo
            do n=1,nfvar
                read(fu) fvar(:,n)
            enddo
            read(fu) oro             
            read(fu) topo            
            read(fu) mask

            close(fu)

        endsubroutine colm_input

        subroutine colm_output(ana_dir,ens,num,numpatch,idate,          &
                              ixy_patch,jxy_patch,mxy_patch,wtxy_patch ,&
                              ftune,fcon,fvar,oro,topo,mask             )
            character(*),  intent(in)    ::    ana_dir
            integer     ,  intent(in)    ::    ens
            integer     ,  intent(in)    ::    num

            integer     ,  intent(in)   ::    numpatch
            integer     ,  intent(in)   ::    idate(3)
            integer,  allocatable, intent(inout)  ::   ixy_patch(:)
            integer,  allocatable, intent(inout)  ::   jxy_patch(:)
            integer,  allocatable, intent(inout)  ::   mxy_patch(:)
            real(r8), allocatable, intent(inout)  ::   wtxy_patch(:)
            real(r8), allocatable, intent(inout)  ::   ftune(:)
            real(r8), allocatable, intent(inout)  ::   fcon(:,:)
            real(r8), allocatable, intent(inout)  ::   fvar(:,:)
            real(r8), allocatable, intent(inout)  ::   oro(:)
            real(r8), allocatable, intent(inout)  ::   topo(:,:)
            real(r8), allocatable, intent(inout)  ::   mask(:,:)

            character(len=1024)          ::    fn
            character(len=256)           ::    fmt
            character(len=20), parameter ::    fi = 'colminput'
            integer     ,  parameter     ::    fu = 101
            integer                      ::    n
            
            write(fmt, '(A,I0,A)') '(A,A,I0.',num,',A,A)' 
            write(fn,fmt) trim(ana_dir),'/', ens,'/', fi
            print*, 'write out '//trim(fn)
            open(fu,file=trim(fn),form='unformatted',convert="little_endian",access='sequential',status='replace',action='write')
            write(fu) numpatch
            write(fu) idate    
            write(fu) ixy_patch  
            write(fu) jxy_patch  
            write(fu) mxy_patch  
            write(fu) wtxy_patch 
            write(fu) ftune
            do n=1,nfcon
                write(fu) fcon(:,n)
            enddo
            do n=1,nfvar
                write(fu) fvar(:,n)
            enddo
            write(fu) oro             
            write(fu) topo            
            write(fu) mask
            close(fu)

        endsubroutine colm_output

        subroutine colm_deallocate(ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                                ftune,fcon,fvar,oro,topo,mask)

            integer,  allocatable, intent(out)    ::   ixy_patch(:)
            integer,  allocatable, intent(out)    ::   jxy_patch(:)
            integer,  allocatable, intent(out)    ::   mxy_patch(:)
            real(r8), allocatable, intent(out)    ::   wtxy_patch(:)
            real(r8), allocatable, intent(out)    ::   ftune(:)
            real(r8), allocatable, intent(out)    ::   fcon(:,:)
            real(r8), allocatable, intent(out)    ::   fvar(:,:)
            real(r8), allocatable, intent(out)    ::   oro(:)
            real(r8), allocatable, intent(out)    ::   topo(:,:)
            real(r8), allocatable, intent(out)    ::   mask(:,:)

            deallocate(ixy_patch       ) 
            deallocate(jxy_patch       ) 
            deallocate(mxy_patch       ) 
            deallocate(wtxy_patch      ) 
            deallocate(ftune           ) 
            deallocate(fcon            ) 
            deallocate(fvar            ) 
            deallocate(oro             ) 
            deallocate(topo            ) 
            deallocate(mask            ) 
            ! print *, '* deallocate var memory over (colm)'

        endsubroutine colm_deallocate

        subroutine obs_read(obs_dir, ens_size, ens, num,        &
                            nobs_raw, rlat, rlon, tbb, qc_flag, &
                            oberr, tbb_bmo, first_read_obs      )
            character(*),  intent(in)    ::    obs_dir
            integer     ,  intent(in)    ::    ens_size
            integer     ,  intent(in)    ::    ens
            integer     ,  intent(in)    ::    num

            integer               , intent(inout)        ::    nobs_raw
            real    , allocatable , intent(inout)        ::    rlat         (:)
            real    , allocatable , intent(inout)        ::    rlon         (:)
            real    , allocatable , intent(inout)        ::    tbb          (:)
            integer , allocatable , intent(inout)        ::    qc_flag      (:,:)
            real    , allocatable , intent(inout)        ::    oberr        (:) 
            real    , allocatable , intent(inout)        ::    tbb_bmo      (:,:)
            logical               , intent(in)           ::    first_read_obs

            character(len=1024)          ::    fn
            character(len=256)           ::    fmt
            character(len=20), parameter ::    fi = 'omb'
            integer     ,  parameter     ::    fu = 300
            integer                      ::    n
            integer                      ::    ncid, dimid, varid
            
            write(fmt, '(A,I0,".",I0,A)') '(A,A,A,A,I',num,num,',A)' !'(A,A,I0.',num,',A)' 
            write(fn,fmt) trim(obs_dir),'/', trim(fi),'.', ens,'.nc'
            print*, 'read in '//trim(fn)

            call check(nf90_open(trim(fn),nf90_nowrite,ncid))
            call check(nf90_inq_dimid(ncid,'index',dimid))
            call check(nf90_inquire_dimension(ncid,dimid,len=nobs_raw))

            ! print *, 'nobs_raw = ',nobs_raw

            if(first_read_obs) then
                allocate(rlat    (nobs_raw))
                allocate(rlon    (nobs_raw))
                allocate(tbb     (nobs_raw))
                allocate(qc_flag (nobs_raw,ens_size))
                allocate(oberr   (nobs_raw))
                allocate(tbb_bmo (nobs_raw,ens_size))
                ! print *, '* allocate obs space over'
            endif

            if(first_read_obs) then
                call check(nf90_inq_varid(ncid,'rlat',varid))
                call check(nf90_get_var(ncid,varid,rlat(:)))
                call check(nf90_inq_varid(ncid,'rlon',varid))
                call check(nf90_get_var(ncid,varid,rlon(:)))
                call check(nf90_inq_varid(ncid,'tbb',varid))
                call check(nf90_get_var(ncid,varid,tbb(:)))
                call check(nf90_inq_varid(ncid,'oberr',varid))
                call check(nf90_get_var(ncid,varid,oberr(:)))
            endif

            call check(nf90_inq_varid(ncid,'qc_flag',varid))
            call check(nf90_get_var(ncid,varid,qc_flag(:,ens)))
            call check(nf90_inq_varid(ncid,'tbb_BMO',varid))
            call check(nf90_get_var(ncid,varid,tbb_bmo(:,ens)))

            call check(nf90_close(ncid))

        endsubroutine obs_read

        subroutine obs_deallocate(rlat, rlon, tbb, qc_flag, oberr, tbb_bmo)
            real    , allocatable , intent(out)         ::    rlat         (:)
            real    , allocatable , intent(out)         ::    rlon         (:)
            real    , allocatable , intent(out)         ::    tbb          (:)
            integer , allocatable , intent(out)         ::    qc_flag      (:,:)
            real    , allocatable , intent(out)         ::    oberr        (:) 
            real    , allocatable , intent(out)         ::    tbb_bmo      (:,:)
            
            deallocate(rlat    )
            deallocate(rlon    )
            deallocate(tbb     )
            deallocate(qc_flag )
            deallocate(oberr   )
            deallocate(tbb_bmo )
        end subroutine obs_deallocate

        subroutine check(status,str)
            integer  , intent(in)               ::    status
            character(*), optional, intent(in)  ::    str    
            character(:), allocatable           ::    str2
            if(status /= nf90_noerr) then

                str2 = trim(nf90_strerror(status))
                if(present(str)) then
                    str2 = str2//':'//trim(str)
                endif
                write(*,"('Error: ', A)") trim(str2)
                stop 21
            endif
        endsubroutine check

endmodule io_control
