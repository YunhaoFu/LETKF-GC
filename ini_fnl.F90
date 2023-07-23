module ini_fnl

    use precision

    implicit none

    private

include './parameter.h'

    public        ::        da_init, da_final

    contains

        subroutine da_init(bkg_dir,ana_dir,obs_dir,ens_size,    &
                           qc_ens,omb_max,infl,radius)
            character(len=256)                ::    cmd
            character(len=1024),intent(out)   ::    bkg_dir        ! bkg directory
            character(len=1024),intent(out)   ::    ana_dir        ! ana directory
            character(len=1024),intent(out)   ::    obs_dir        ! obs directory
            integer            ,intent(out)   ::    ens_size       ! ensemble size
            real               ,intent(out)   ::    qc_ens         ! (0,1) for ensemble qc
            real               ,intent(out)   ::    omb_max        ! maximum of abs(omb) 
            real               ,intent(out)   ::    infl           ! inflation rate for LETKF
            real               ,intent(out)   ::    radius(2)      ! localization radius for LETKF
            character(len=1024)               ::    nmlfile        ! absolute dir+name for nml file
                                                                   ! which contains the info above
            integer, parameter                ::    argv_cnt  = 1
            integer, parameter                ::    fu = 101

            namelist /letkf/ bkg_dir , &
                             ana_dir , &
                             obs_dir , &
                             ens_size, &
                             qc_ens  , &
                             omb_max , &
                             infl    , &
                             radius

            write(*, "('============================== ** Initialize ** =================================')")
            print *
            if (command_argument_count() < argv_cnt) then
              call get_command_argument(0, cmd)
              write(*, "('Usage: ', A, 'dir+name for nml')") trim(cmd)
              stop 18
            else
              call get_command_argument(1, nmlfile)
              ! write(*, "('nml: ', A)") trim(nmlfile)
            endif
            
            open(unit=fu,file=trim(nmlfile))
            read(unit=fu,nml=letkf)
            close(unit=fu)

            write(*, "('bkg dir  : ', A)") trim(bkg_dir)
            write(*, "('ana dir  : ', A)") trim(ana_dir)
            write(*, "('obs dir  : ', A)") trim(obs_dir)
            write(*, "('ens_size : ', I3.3)") ens_size
            write(*, "('qc_ens   : ', F10.5)") qc_ens
            write(*, "('omb_max  : ', F10.5)") omb_max
            write(*, "('infl     : ', F10.5)") infl
            write(*, "('radius   : ', 2F14.5)") radius

        endsubroutine da_init

        subroutine da_final(nobs,lonxy_atm,latxy_atm,              &
                               grid2patch_start, grid2patch_count, &
                               lonxy_lnd,latxy_lnd,                &
                               v1,lb_patch)                              ! todo_colm_invar #1 
            integer              , intent(in)       ::      nobs
            real(r4), allocatable, intent(out)      ::      lonxy_atm(:,:)
            real(r4), allocatable, intent(out)      ::      latxy_atm(:,:)
            integer , allocatable, intent(out)      ::      grid2patch_start(:,:)
            integer , allocatable, intent(out)      ::      grid2patch_count(:,:)
            real(r4), allocatable, intent(out)      ::      lonxy_lnd(:,:)
            real(r4), allocatable, intent(out)      ::      latxy_lnd(:,:)
            integer , allocatable, intent(out)      ::      v1(:,:) ! todo_colm_invar #2
            integer , allocatable, intent(out)      ::      lb_patch(:,:)
            integer                                 ::      cnt1=0

            if(nobs > 0) then
                deallocate(lonxy_atm)
                deallocate(latxy_atm)
                deallocate(grid2patch_start)
                deallocate(grid2patch_count)
                deallocate(lonxy_lnd)
                deallocate(latxy_lnd)
                cnt1 = cnt1 + 1
                deallocate(v1)
                ! todo_colm_invar #3 (can be added more, v2 v3 v4...)
                if(cnt1 /= rd_invar_cnt) then
                    print *, 'wrong deallocate setting for colm_invar: '
                    print *, 'num = ', rd_invar_cnt, ', due num = ', cnt1
                    stop 22
                endif
                deallocate(lb_patch)
            endif
            print *, ''
            write(*, "('============================== ** Finalized ** =================================')")

        endsubroutine da_final

endmodule ini_fnl
