! These diagnostics are implemented for MPAS7 vs MPAS7-acc verification only
! and should not be used for production.
!
module twc_diagnostics

    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines

    public :: twc_diagnostics_setup, &
              twc_diagnostics_update, &
              twc_diagnostics_compute, &
              twc_diagnostics_reset

    private

    type (MPAS_pool_type), pointer :: mesh, state, diag, phys
    real (kind=RKIND), pointer :: dt
    integer, pointer :: nCells, nVertLevels
    logical :: need_rain, need_conv, need_snow, need_snow_total, need_snow_ratio

    contains


!
! Setup
!
    subroutine twc_diagnostics_setup(configs, all_pools)

        use mpas_atm_diagnostics_utils, only : mpas_stream_inclusion_count

        implicit none

        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: all_pools

        call mpas_pool_get_config(configs, 'config_dt', dt)
        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', phys)

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

        need_rain = (mpas_stream_inclusion_count('rain_bucket', direction=MPAS_STREAM_OUTPUT) > 0)
        need_conv = (mpas_stream_inclusion_count('conv_bucket', direction=MPAS_STREAM_OUTPUT) > 0)
        need_snow = (mpas_stream_inclusion_count('snow_bucket', direction=MPAS_STREAM_OUTPUT) > 0)
        need_snow_total = (mpas_stream_inclusion_count('snow_total', direction=MPAS_STREAM_OUTPUT) > 0)
        need_snow_ratio = need_snow .or. need_snow_total .or. &
                          (mpas_stream_inclusion_count('snow_ratio', direction=MPAS_STREAM_OUTPUT) > 0)
   
    end subroutine twc_diagnostics_setup

!
! Update
!
    subroutine twc_diagnostics_update()
        implicit none

        real (kind=RKIND), dimension(:), pointer :: rainncv, raincv, snowncv, graupelncv
        real (kind=RKIND), dimension(:), pointer :: rain_bucket, conv_bucket, snow_bucket
        real (kind=RKIND), dimension(:), pointer :: snow_total, snow_ratio

        integer :: iCell
        real :: rain_dt, conv_dt, snow_dt

        if (need_snow_ratio) then
            call update_kuchera_snow_ratio()
        end if

        if (need_rain .or. need_conv .or. need_snow .or. need_snow_total) then
            ! Inputs
            call mpas_pool_get_array(phys, 'rainncv', rainncv)
            call mpas_pool_get_array(phys, 'raincv', raincv)
            call mpas_pool_get_array(phys, 'snowncv', snowncv)
            call mpas_pool_get_array(phys, 'graupelncv', graupelncv)
            call mpas_pool_get_array(diag, 'snow_ratio', snow_ratio)
            ! Outputs
            call mpas_pool_get_array(diag, 'rain_bucket', rain_bucket)
            call mpas_pool_get_array(diag, 'conv_bucket', conv_bucket)
            call mpas_pool_get_array(diag, 'snow_bucket', snow_bucket)
            call mpas_pool_get_array(diag, 'snow_total', snow_total)

            do iCell = 1, nCells
                ! Timestep totals
                rain_dt = max(0._RKIND, rainncv(iCell) - snowncv(iCell) - graupelncv(iCell))
                conv_dt = max(0._RKIND, raincv(iCell))
                snow_dt = max(0._RKIND, snowncv(iCell) + graupelncv(iCell))

                ! Accumulate
                rain_bucket(iCell) = rain_bucket(iCell) + rain_dt
                conv_bucket(iCell) = conv_bucket(iCell) + conv_dt
                snow_bucket(iCell) = snow_bucket(iCell) + snow_dt
                snow_total(iCell) = snow_total(iCell) + snow_dt * snow_ratio(iCell) * 0.001_RKIND
            end do
        end if

    end subroutine twc_diagnostics_update


    subroutine update_kuchera_snow_ratio()
        use mpas_constants, only : rvord

        implicit none

        integer, pointer :: index_qv
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: exner
        real (kind=RKIND), dimension(:,:), pointer :: pressure_b
        real (kind=RKIND), dimension(:,:), pointer :: pressure_p
        real (kind=RKIND), dimension(:), pointer :: snowncv
        real (kind=RKIND), dimension(:), pointer :: graupelncv
        real (kind=RKIND), dimension(:), pointer :: snow_ratio

        real (kind=RKIND) :: snow_ratio_tmp, temperature_max, pressure, temperature
        integer :: iCell, k

        real (kind=RKIND), parameter :: snow_ratio_default = 1._RKIND
        real (kind=RKIND), parameter :: snow_ratio_min = 8._RKIND
        real (kind=RKIND), parameter :: snow_ratio_max = 18._RKIND

        ! Inputs
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)
        call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
        call mpas_pool_get_array(state, 'scalars', scalars, 1)
        call mpas_pool_get_array(diag, 'exner', exner)
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
        call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
        call mpas_pool_get_array(phys, 'snowncv', snowncv)
        call mpas_pool_get_array(phys, 'graupelncv', graupelncv)
        ! Outputs
        call mpas_pool_get_array(diag, 'snow_ratio', snow_ratio)

        do iCell = 1, nCells
            snow_ratio_tmp = snow_ratio_default
            if (snowncv(iCell) + graupelncv(iCell) > 0) then
                temperature_max = 0._RKIND
                kloop: do k = 1, nVertLevels-1
                   pressure = (pressure_p(k,iCell)+pressure_b(k,iCell))/100._RKIND
                   if (pressure >= 500) then
                      temperature = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell)
                      temperature_max = max(temperature_max,temperature)
                   else
                      exit kloop
                   end if
                end do kloop
                if (temperature_max > 271.16) then
                    snow_ratio_tmp = 12.0+2.0*(271.16-temperature_max)
                else
                    snow_ratio_tmp = 12.0+(271.16-temperature_max)
                end if

            end if
            snow_ratio(iCell) = max(snow_ratio_min, min(snow_ratio_tmp, snow_ratio_max))
        end do

    end subroutine update_kuchera_snow_ratio


!
! Compute
!
    subroutine twc_diagnostics_compute()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        if (MPAS_field_will_be_written('dewpoint_2m')) then
            call compute_dewpoint_2m()
        end if

        if (MPAS_field_will_be_written('refl10cm_prate')) then
            call compute_refl10cm_prate()
        end if

    end subroutine twc_diagnostics_compute



    ! Dewpoint temperature (K) following Bolton (1980)
    subroutine compute_dewpoint_2m()
        implicit none

        real (kind=RKIND), dimension(:), pointer :: dewpoint_2m, q2, t2m, psfc
        real (kind=RKIND) :: evp
        integer           :: iCell

        ! Inputs
        call mpas_pool_get_array(phys, 'q2', q2)
        call mpas_pool_get_array(phys, 't2m', t2m)
        call mpas_pool_get_array(diag, 'surface_pressure', psfc)
        ! Outputs
        call mpas_pool_get_array(diag, 'dewpoint_2m', dewpoint_2m)

        do iCell = 1,nCells
            evp = 0.01_RKIND * psfc(iCell) * q2(iCell) / (q2(iCell) + 0.622_RKIND)
            evp = max(evp, 1.0e-8_RKIND)
            dewpoint_2m(iCell) = (243.5_RKIND * log(evp/6.112_RKIND)) / (17.67_RKIND - log(evp/6.112_RKIND))
            dewpoint_2m(iCell) = dewpoint_2m(iCell) + 273.15_RKIND
            dewpoint_2m(iCell) = min(dewpoint_2m(iCell), t2m(iCell))
        end do
    end subroutine compute_dewpoint_2m
        
    ! From cpox_diagnostics
    subroutine compute_refl10cm_prate()
        implicit none

        real (kind=RKIND), dimension(:), pointer :: rainncv, raincv, refl
        real (kind=RKIND) :: prate
        integer :: iCell

        ! Inputs
        call mpas_pool_get_array(phys, 'rainncv', rainncv)
        call mpas_pool_get_array(phys, 'raincv', raincv)
        ! Outputs
        call mpas_pool_get_array(diag, 'refl10cm_prate', refl)

        do iCell = 1,nCells
            refl(iCell) = 0.0_RKIND
            prate = max((rainncv(iCell) + raincv(iCell)) / dt, 0.0_RKIND)
            if (prate .ge. 1.e-6_RKIND) then
                refl(iCell) = 10.*log10(300.*((prate*3600.)**1.4))
            end if
        end do

    end subroutine compute_refl10cm_prate


!
! Reset
!
    subroutine twc_diagnostics_reset()
        use mpas_pool_routines, only : mpas_pool_get_array
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        real (kind=RKIND), dimension(:), pointer :: p

        if (MPAS_field_will_be_written('rain_bucket')) then
            call mpas_pool_get_array(diag, 'rain_bucket', p)
            p(1:nCells) = 0.0_RKIND
        end if
        if (MPAS_field_will_be_written('conv_bucket')) then
            call mpas_pool_get_array(diag, 'conv_bucket', p)
            p(1:nCells) = 0.0_RKIND
        end if
        if (MPAS_field_will_be_written('snow_bucket')) then
            call mpas_pool_get_array(diag, 'snow_bucket', p)
            p(1:nCells) = 0.0_RKIND
        end if
    end subroutine twc_diagnostics_reset



end module twc_diagnostics
