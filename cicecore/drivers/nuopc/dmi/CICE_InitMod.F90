!=======================================================================
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
! 2008: E. Hunke moved ESMF code to its own driver

      module CICE_InitMod

      use ice_kinds_mod
      use ice_exit, only: abort_ice
      use ice_fileunits, only: init_fileunits, nu_diag
      use ice_memusage, only: ice_memusage_init, ice_memusage_print
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc, only: icepack_init_itd, icepack_init_itd_hist
      use icepack_intfc, only: icepack_init_fsd_bounds, icepack_init_wave
      use icepack_intfc, only: icepack_init_snow, icepack_init_radiation
      use icepack_intfc, only: icepack_configure
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_flags, &
          icepack_query_tracer_indices, icepack_query_tracer_sizes

      implicit none
      private
      public :: CICE_Initialize, cice_init

!=======================================================================

      contains

!=======================================================================

!  Initialize the basic state, grid and all necessary parameters for
!  running the CICE model.  Return the initial state in routine
!  export state.
!  Note: This initialization driver is designed for standalone and
!        CESM-coupled applications.  For other
!        applications (e.g., standalone CAM), this driver would be
!        replaced by a different driver that calls subroutine cice_init,
!        where most of the work is done.

      subroutine CICE_Initialize(mpi_comm)

      integer (kind=int_kind), optional, intent(in) :: mpi_comm ! communicator from nuopc
      character(len=*), parameter :: subname='(CICE_Initialize)'
   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      if (present(mpi_comm)) then
          call cice_init(mpi_comm)
      else
          call cice_init()
      endif

      end subroutine CICE_Initialize

!=======================================================================
!
!  Initialize CICE model.

      subroutine cice_init(mpi_comm)

      use ice_arrays_column, only: hin_max, c_hi_range, alloc_arrays_column
      use ice_arrays_column, only: floe_rad_l, floe_rad_c, &
          floe_binwidth, c_fsd_range
      use ice_state, only: alloc_state
      use ice_flux_bgc, only: alloc_flux_bgc
      use ice_calendar, only: dt, dt_dyn, write_ic, &
          init_calendar, advance_timestep, calc_timesteps
      use ice_communicate, only: init_communicate, my_task, master_task
      use ice_diagnostics, only: init_diags
      use ice_domain, only: init_domain_blocks
      use ice_domain_size, only: ncat, nfsd
      use ice_dyn_eap, only: init_eap
      use ice_dyn_evp, only: init_evp
      use ice_dyn_vp, only: init_vp
      use ice_dyn_shared, only: kdyn
      use ice_flux, only: init_coupler_flux, init_history_therm, &
          init_history_dyn, init_flux_atm, init_flux_ocn, alloc_flux
      use ice_forcing, only: init_forcing_ocn, init_forcing_atmo, &
          get_forcing_atmo, get_forcing_ocn, get_wave_spec, init_snowtable
      use ice_forcing_bgc, only: get_forcing_bgc, get_atm_bgc, &
          faero_default, alloc_forcing_bgc, fiso_default
      use ice_grid, only: init_grid1, init_grid2, alloc_grid, dealloc_grid
      use ice_history, only: init_hist, accum_hist
      use ice_restart_shared, only: restart, runtype
      use ice_init, only: input_data, init_state
      use ice_init_column, only: init_thermo_vertical, init_shortwave, init_zbgc, input_zbgc, count_tracers
      use ice_kinds_mod
      use ice_restoring, only: ice_HaloRestore_init
      use ice_timers, only: timer_total, init_ice_timers, ice_timer_start
      use ice_transport_driver, only: init_transport

      integer (kind=int_kind), optional, intent(in) :: &
         mpi_comm ! communicator for sequential ccsm

      logical(kind=log_kind) :: tr_aero, tr_zaero, skl_bgc, z_tracers, &
         tr_iso, tr_fsd, wave_spec, tr_snow
      character(len=char_len) :: snw_aging_table
      character(len=*), parameter :: subname = '(cice_init)'

      if (present(mpi_comm)) then
          call init_communicate(mpi_comm)     ! initial setup for message passing
      else
          call init_communicate     ! initial setup for message passing
      endif
      call init_fileunits       ! unit numbers

      ! tcx debug, this will create a different logfile for each pe
      ! if (my_task /= master_task) nu_diag = 100+my_task

      call icepack_configure()  ! initialize icepack
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(trim(subname), &
          file=__FILE__,line= __LINE__)

      call input_data           ! namelist variables
      call input_zbgc           ! vertical biogeochemistry namelist
      call count_tracers        ! count tracers

      ! Call this as early as possible, must be after memory_stats is read
      if (my_task == master_task) then
         call ice_memusage_init(nu_diag)
         call ice_memusage_print(nu_diag,subname//':start')
      endif

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
      call alloc_grid           ! allocate grid arrays
      call alloc_arrays_column  ! allocate column arrays
      call alloc_state          ! allocate state arrays
      call alloc_flux_bgc       ! allocate flux_bgc arrays
      call alloc_flux           ! allocate flux arrays
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables
      call init_zbgc            ! vertical biogeochemistry initialization
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file

      if (kdyn == 1) then
         call init_evp
      else if (kdyn == 2) then
         call init_eap          ! define eap dynamics parameters, variables
      else if (kdyn == 3) then
         call init_vp           ! define vp dynamics parameters, variables
      endif

      call init_coupler_flux    ! initialize fluxes exchanged with coupler

      call init_thermo_vertical ! initialize vertical thermodynamics

      call icepack_init_itd(hin_max=hin_max)  ! ice thickness distribution
      if (my_task == master_task) then
         call icepack_init_itd_hist(hin_max=hin_max, c_hi_range=c_hi_range) ! output
      endif

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(trim(subname), &
          file=__FILE__,line= __LINE__)

      if (tr_fsd) call icepack_init_fsd_bounds ( &
         floe_rad_l_out = floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c_out = floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth_out = floe_binwidth, &  ! fsd size bin width in m (radius)
         c_fsd_range_out = c_fsd_range,   &  ! string for history output
         write_diags=(my_task == master_task))  ! write diag on master only
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

#ifndef CICE_IN_NEMO
      call init_forcing_ocn(dt) ! initialize sss and sst from data
#endif
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
      call ice_HaloRestore_init ! restored boundary conditions

      call icepack_query_parameters(skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, &
          wave_spec_out=wave_spec, snw_aging_table_out=snw_aging_table)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(trim(subname), &
          file=__FILE__,line= __LINE__)

      if (skl_bgc .or. z_tracers) call alloc_forcing_bgc ! allocate biogeochemistry arrays

      call init_restart         ! initialize restart variables
      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables
      call calc_timesteps       ! update timestep counter if not using npt_unit="1"
      call icepack_init_radiation ! initialize icepack shortwave tables

      call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_zaero_out=tr_zaero)
      call icepack_query_tracer_flags(tr_iso_out=tr_iso, tr_snow_out=tr_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(trim(subname), &
          file=__FILE__,line= __LINE__)

      ! Initialize shortwave components using swdn from previous timestep
      ! if restarting. These components will be scaled to current forcing
      ! in prep_radiation.
      if (trim(runtype) == 'continue' .or. restart) &
         call init_shortwave    ! initialize radiative transfer

      call advance_timestep()

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

#ifndef CICE_IN_NEMO
      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)
#endif

#ifndef coupled
#ifndef CESMCOUPLED
      if (tr_fsd .and. wave_spec) call get_wave_spec ! wave spectrum in ice
      call get_forcing_atmo     ! atmospheric forcing from data
#ifndef CICE_DMI
      call get_forcing_ocn(dt)  ! ocean forcing from data
#endif

      ! snow aging lookup table initialization
      if (tr_snow) then         ! advanced snow physics
         call icepack_init_snow()
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)
         if (snw_aging_table(1:4) /= 'test') then
            call init_snowtable()
         endif
      endif

      ! isotopes
      if (tr_iso)     call fiso_default                 ! default values

      ! aerosols
      ! if (tr_aero)  call faero_data                   ! data file
      ! if (tr_zaero) call fzaero_data                  ! data file (gx1)
      if (tr_aero .or. tr_zaero)  call faero_default    ! default values
      if (skl_bgc .or. z_tracers) call get_forcing_bgc  ! biogeochemistry
#endif
#endif
      if (z_tracers) call get_atm_bgc                   ! biogeochemistry

      if (runtype == 'initial' .and. .not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler

      if (write_ic) call accum_hist(dt) ! write initial conditions

      call dealloc_grid         ! deallocate temporary grid arrays
      if (my_task == master_task) then
         call ice_memusage_print(nu_diag,subname//':end')
      endif

      end subroutine cice_init

!=======================================================================

      subroutine init_restart

      use ice_arrays_column, only: dhsn
      use ice_blocks, only: nx_block, ny_block
      use ice_calendar, only: calendar
      use ice_constants, only: c0
      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat, n_iso, n_aero, nfsd, nslyr
      use ice_dyn_eap, only: read_restart_eap
      use ice_dyn_shared, only: kdyn
      use ice_grid, only: tmask
      use ice_init, only: ice_ic
      use ice_init_column, only: init_age, init_FY, init_lvl, init_snowtracers, &
          init_meltponds_lvl, init_meltponds_topo, &
          init_isotope, init_aerosol, init_hbrine, init_bgc, init_fsd
      use ice_restart_column, only: restart_age, read_restart_age, &
          restart_FY, read_restart_FY, restart_lvl, read_restart_lvl, &
          restart_pond_lvl, read_restart_pond_lvl, &
          restart_pond_topo, read_restart_pond_topo, &
          restart_snow, read_restart_snow, &
          restart_fsd, read_restart_fsd, &
          restart_iso, read_restart_iso, &
          restart_aero, read_restart_aero, &
          restart_hbrine, read_restart_hbrine, &
          restart_bgc
      use ice_flux, only: Tf
      use ice_restart_driver, only: restartfile
      use ice_restart_shared, only: runtype, restart
      use ice_state ! almost everything

      integer(kind=int_kind) :: &
         i, j        , & ! horizontal indices
         iblk            ! block index
      logical(kind=log_kind) :: &
          tr_iage, tr_FY, tr_lvl, tr_pond_lvl, &
          tr_pond_topo, tr_snow, tr_fsd, tr_iso, tr_aero, tr_brine, &
          skl_bgc, z_tracers
      integer(kind=int_kind) :: &
          ntrcr
      integer(kind=int_kind) :: &
          nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, &
          nt_smice, nt_smliq, nt_rhos, nt_rsnw, &
          nt_iage, nt_FY, nt_aero, nt_fsd, nt_isosno, nt_isoice

      character(len=*), parameter :: subname = '(init_restart)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_parameters(skl_bgc_out=skl_bgc, &
           z_tracers_out=z_tracers)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
           tr_lvl_out=tr_lvl, tr_pond_lvl_out=tr_pond_lvl, &
           tr_pond_topo_out=tr_pond_topo, tr_aero_out=tr_aero, tr_brine_out=tr_brine, &
           tr_snow_out=tr_snow, tr_fsd_out=tr_fsd, tr_iso_out=tr_iso)
      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, &
           nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, &
           nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_aero_out=nt_aero, nt_fsd_out=nt_fsd, &
           nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, &
           nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw, &
           nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (trim(runtype) == 'continue') then
         ! start from core restart file
         call restartfile()           ! given by pointer in ice_in
         call calendar()              ! update time parameters
         if (kdyn == 2) call read_restart_eap ! EAP
      else if (restart) then          ! ice_ic = core restart file
         call restartfile (ice_ic)    !  or 'default' or 'none'
         !!! uncomment to create netcdf
         ! call restartfile_v4 (ice_ic)  ! CICE v4.1 binary restart file
         !!! uncomment if EAP restart data exists
         ! if (kdyn == 2) call read_restart_eap
      endif

      ! tracers
      ! ice age tracer
      if (tr_iage) then
         if (trim(runtype) == 'continue') restart_age = .true.
         if (restart_age) then
            call read_restart_age
         else
            do iblk = 1, nblocks
               call init_age(trcrn(:,:,nt_iage,:,iblk))
            enddo ! iblk
         endif
      endif
      ! first-year area tracer
      if (tr_FY) then
         if (trim(runtype) == 'continue') restart_FY = .true.
         if (restart_FY) then
            call read_restart_FY
         else
            do iblk = 1, nblocks
               call init_FY(trcrn(:,:,nt_FY,:,iblk))
            enddo ! iblk
         endif
      endif
      ! level ice tracer
      if (tr_lvl) then
         if (trim(runtype) == 'continue') restart_lvl = .true.
         if (restart_lvl) then
            call read_restart_lvl
         else
            do iblk = 1, nblocks
               call init_lvl(iblk,trcrn(:,:,nt_alvl,:,iblk), &
                             trcrn(:,:,nt_vlvl,:,iblk))
            enddo ! iblk
         endif
      endif
      ! level-ice melt ponds
      if (tr_pond_lvl) then
         if (trim(runtype) == 'continue') restart_pond_lvl = .true.
         if (restart_pond_lvl) then
            call read_restart_pond_lvl
         else
            do iblk = 1, nblocks
               call init_meltponds_lvl(trcrn(:,:,nt_apnd,:,iblk), &
                                       trcrn(:,:,nt_hpnd,:,iblk), &
                                       trcrn(:,:,nt_ipnd,:,iblk), &
                                       dhsn(:,:,:,iblk))
            enddo ! iblk
         endif
      endif
      ! topographic melt ponds
      if (tr_pond_topo) then
         if (trim(runtype) == 'continue') restart_pond_topo = .true.
         if (restart_pond_topo) then
            call read_restart_pond_topo
         else
            do iblk = 1, nblocks
               call init_meltponds_topo(trcrn(:,:,nt_apnd,:,iblk), &
                                        trcrn(:,:,nt_hpnd,:,iblk), &
                                        trcrn(:,:,nt_ipnd,:,iblk))
            enddo ! iblk
         endif ! .not. restart_pond
      endif

      ! snow redistribution/metamorphism
      if (tr_snow) then
         if (trim(runtype) == 'continue') restart_snow = .true.
         if (restart_snow) then
            call read_restart_snow
         else
            do iblk = 1, nblocks
               call init_snowtracers(trcrn(:,:,nt_smice:nt_smice+nslyr-1,:,iblk), &
                                     trcrn(:,:,nt_smliq:nt_smliq+nslyr-1,:,iblk), &
                                     trcrn(:,:,nt_rhos :nt_rhos +nslyr-1,:,iblk), &
                                     trcrn(:,:,nt_rsnw :nt_rsnw +nslyr-1,:,iblk))
            enddo ! iblk
         endif
      endif

      ! floe size distribution
      if (tr_fsd) then
         if (trim(runtype) == 'continue') restart_fsd = .true.
         if (restart_fsd) then
            call read_restart_fsd
         else
            call init_fsd(trcrn(:,:,nt_fsd:nt_fsd+nfsd-1,:,:))
         endif
      endif

      ! isotopes
      if (tr_iso) then
         if (trim(runtype) == 'continue') restart_iso = .true.
         if (restart_iso) then
            call read_restart_iso
         else
            do iblk = 1, nblocks
               call init_isotope(trcrn(:,:,nt_isosno:nt_isosno+n_iso-1,:,iblk), &
                                 trcrn(:,:,nt_isoice:nt_isoice+n_iso-1,:,iblk))
            enddo ! iblk
         endif
      endif

      if (tr_aero) then ! ice aerosol
         if (trim(runtype) == 'continue') restart_aero = .true.
         if (restart_aero) then
            call read_restart_aero
         else
            do iblk = 1, nblocks
               call init_aerosol(trcrn(:,:,nt_aero:nt_aero+4*n_aero-1,:,iblk))
            enddo ! iblk
         endif ! .not. restart_aero
      endif

      if (trim(runtype) == 'continue') then
         if (tr_brine) restart_hbrine = .true.
         if (skl_bgc .or. z_tracers) restart_bgc = .true.
      endif

      if (tr_brine .or. skl_bgc) then ! brine height tracer
         call init_hbrine
         if (tr_brine .and. restart_hbrine) call read_restart_hbrine
      endif

      if (skl_bgc .or. z_tracers) then ! biogeochemistry
         if (tr_fsd) then
            write (nu_diag,*) 'FSD implementation incomplete for use with BGC'
            call icepack_warnings_flush(nu_diag)
            if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
               file=__FILE__, line=__LINE__)
         endif
         call init_bgc
      endif

      !-----------------------------------------------------------------
      ! aggregate tracers
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j,iblk)) then
            call icepack_aggregate(aicen = aicen(i,j,:,iblk),     &
                                   trcrn = trcrn(i,j,:,:,iblk),   &
                                   vicen = vicen(i,j,:,iblk),     &
                                   vsnon = vsnon(i,j,:,iblk),     &
                                   aice  = aice (i,j,  iblk),     &
                                   trcr  = trcr (i,j,:,iblk),     &
                                   vice  = vice (i,j,  iblk),     &
                                   vsno  = vsno (i,j,  iblk),     &
                                   aice0 = aice0(i,j,  iblk),     &
                                   trcr_depend   = trcr_depend,   &
                                   trcr_base     = trcr_base,     &
                                   n_trcr_strata = n_trcr_strata, &
                                   nt_strata     = nt_strata,     &
                                   Tf            = Tf(i,j,iblk))
         else
            ! tcraig, reset all tracer values on land to zero
            trcrn(i,j,:,:,iblk) = c0
         endif
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine init_restart

!=======================================================================

      end module CICE_InitMod

!=======================================================================
