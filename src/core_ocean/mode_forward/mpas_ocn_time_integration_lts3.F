! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts3
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This module contains the LTS3 time integration and LTS init routines
!
!-----------------------------------------------------------------------

module ocn_time_integration_lts3 

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled
   use ocn_wetting_drying

   use ocn_effective_density_in_land_ice
   use ocn_time_varying_forcing 

   use ocn_time_integration_ssprk3

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_time_integrator_lts3, ocn_time_integration_lts3_init, ocn_lts_fast_tends

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_lts3
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This routine integrates one timestep (dt) using an LTS3 time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_lts3(domain,dt)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    !   Local Time Stepping of order 3 scheme
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    type (domain_type), intent(inout) :: domain
    real (kind=RKIND), intent(in) :: dt

    integer :: iCell, iEdge, iRegion, k, ic, ie, im, M, nRegions, err

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: tendPool
    type (mpas_pool_type), pointer :: tracersTendPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: scratchPool
    type (mpas_pool_type), pointer :: swForcingPool

    ! LTS Pools
    type (mpas_pool_type), pointer :: LTSPool
    type (mpas_pool_type), pointer :: tendSlowPool
    type (mpas_pool_type), pointer :: tendSum1stPool
    type (mpas_pool_type), pointer :: tracersTendSum1stPool
    type (mpas_pool_type), pointer :: tendSum2ndPool
    type (mpas_pool_type), pointer :: tracersTendSum2ndPool
    type (mpas_pool_type), pointer :: tendSum3rdPool
    type (mpas_pool_type), pointer :: tracersTendSum3rdPool
    type (mpas_pool_type), pointer :: prevTendSum1stPool, nextTendSum1stPool
    type (mpas_pool_type), pointer :: prevTendSum2ndPool, nextTendSum2ndPool
    type (mpas_pool_type), pointer :: prevTendSum3rdPool, nextTendSum3rdPool
    type (mpas_pool_type), pointer :: prevTendSlowPool, nextTendSlowPool

    ! Dimensions
    integer, pointer :: nCells, nEdges, nEdgesSolve, nVertLevels, haveFine, haveCoarse, haveInterface

    ! Config options
    integer, pointer :: config_dt_scaling_LTS
    logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
    logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
    logical, pointer :: config_use_GM
    logical, pointer :: config_use_cvmix_kpp
    logical, pointer :: config_disable_vel_vmix
    logical, pointer :: config_disable_tr_vmix
    logical, pointer :: config_use_tracerGroup
    logical, pointer :: config_disable_thick_all_tend
    logical, pointer :: config_disable_vel_all_tend
    logical, pointer :: config_disable_tr_all_tend
    logical, pointer :: config_LTS_split
    logical, pointer :: config_LTS_inline_fast_tends
    real (kind=RKIND), pointer :: config_mom_del4
    real (kind=RKIND), pointer :: config_drying_min_cell_height
    logical, pointer :: config_use_wetting_drying
    logical, pointer :: config_verify_not_dry
    character (len=StrKIND), pointer :: config_land_ice_flux_mode
    logical :: activeTracersOnlyIn

    ! State indices
    integer, pointer :: indexTemperature
    integer, pointer :: indexSalinity

    ! Diagnostics Indices
    integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
    integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

    ! Mesh array pointers
    integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop
    integer, dimension(:,:), pointer :: cellsOnEdge
    real (kind=RKIND), dimension(:), pointer :: bottomDepth

    ! Tend Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, layerThicknessTend, normalVelocityTendSlow
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTendSum1st, lowFreqDivergenceTendSum1st, normalVelocityTendSum1st, layerThicknessTendSum1st
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTendSum2nd, lowFreqDivergenceTendSum2nd, normalVelocityTendSum2nd, layerThicknessTendSum2nd
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTendSum3rd, lowFreqDivergenceTendSum3rd, normalVelocityTendSum3rd, layerThicknessTendSum3rd
    real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend, tracersGroupTendSum1st, tracersGroupTendSum2nd, tracersGroupTendSum3rd

    ! Diagnostics Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
    real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
    real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity, normalGMBolusVelocity
    real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
    real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
    real (kind=RKIND), dimension(:), pointer :: gradSSH
    real (kind=RKIND), dimension(:), pointer :: gradSSHX, gradSSHY, gradSSHZ
    real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
    real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, sshGradient

    ! State Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew, normalVelocityFirstStage, normalVelocitySecondStage
    real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew, layerThicknessFirstStage, layerThicknessSecondStage
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew, highFreqThicknessFirstStage, highFreqThicknessSecondStage
    real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew, lowFreqDivergenceFirstStage, lowFreqDivergenceSecondStage 
    real (kind=RKIND), dimension(:,:), pointer :: wettingVelocity
    real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew, sshFirstStage, sshSecondStage
    real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupCur, tracersGroupNew, tracersGroupFirstStage, tracersGroupSecondStage

    ! Diagnostics Field Pointers
    type (field1DReal), pointer :: boundaryLayerDepthField, effectiveDensityField
    type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

    ! Tracer Group Iteartion
    type (mpas_pool_iterator_type) :: groupItr
    character (len=StrKIND) :: modifiedGroupName
    character (len=StrKIND) :: configName

    ! Tidal boundary condition
    logical, pointer :: config_use_tidal_forcing
    character (len=StrKIND), pointer :: config_tidal_forcing_type
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue
    real (kind=RKIND), dimension(:,:), pointer :: restingThickness
    real (kind=RKIND) :: totalDepth

    ! LTS objects
    real (kind=RKIND) :: dtFine, alpha, alphaHat, beta, betaHat, gam, gamHat
    integer, dimension(:,:), pointer :: nCellsInLTSHalo, nEdgesInLTSHalo
    integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
    real (kind=RKIND) :: weightOld, weightNew, weightTend
    real (kind=RKIND) :: weightTendSum1st, weightTendSum2nd, weightTendSum3rd
    
    real (kind=RKIND), pointer :: forcingTimeIncrement

    ! Get config options
    call mpas_pool_get_config(domain % configs, 'config_dt_scaling_LTS', config_dt_scaling_LTS)
    call mpas_pool_get_config(domain % configs, 'config_LTS_split', config_LTS_split)
    call mpas_pool_get_config(domain % configs, 'config_LTS_inline_fast_tends', config_LTS_inline_fast_tends)
    call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
    call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
    call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
    call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
    call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
    call mpas_pool_get_config(domain % configs, 'config_use_GM', config_use_GM)
    call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
    call mpas_pool_get_config(domain % configs, 'config_disable_vel_vmix', config_disable_vel_vmix)
    call mpas_pool_get_config(domain % configs, 'config_disable_tr_vmix', config_disable_tr_vmix)
    call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
    call mpas_pool_get_config(domain % configs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_use_wetting_drying', config_use_wetting_drying)
    call mpas_pool_get_config(domain % configs, 'config_verify_not_dry', config_verify_not_dry)
    call mpas_pool_get_config(domain % configs, 'config_drying_min_cell_height', config_drying_min_cell_height)
    call mpas_pool_get_config(domain % configs, 'config_use_tidal_forcing', config_use_tidal_forcing)
    call mpas_pool_get_config(domain % configs, 'config_tidal_forcing_type', config_tidal_forcing_type)

    activeTracersOnlyIn = .false.
  
    M = config_dt_scaling_LTS
    nRegions = 2

    err = 0

    dtFine = dt / M

    weightOld = 0.75
    weightNew = 0.25
    weightTend = 0.25

    weightTendSum1st = 1.0 / 6.0
    weightTendSum2nd = 1.0 / 6.0
    weightTendSum3rd = 2.0 / 3.0

    call mpas_timer_start("LTS3-init variables and create extra pools for LTS")
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'state', statePool)
       call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_create_pool(tendSum1stPool)
       call mpas_pool_clone_pool(tendPool, tendSum1stPool, 1)
       call mpas_pool_create_pool(tendSum2ndPool)
       call mpas_pool_clone_pool(tendPool, tendSum2ndPool, 1)
       call mpas_pool_create_pool(tendSum3rdPool)
       call mpas_pool_clone_pool(tendPool, tendSum3rdPool, 1)
       call mpas_pool_create_pool(tendSlowPool)
       call mpas_pool_clone_pool(tendPool, tendSlowPool, 1)

       call mpas_pool_add_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
       call mpas_pool_add_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
       call mpas_pool_add_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
       call mpas_pool_add_subpool(block % structs, 'tend_slow', tendSlowPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

       call mpas_pool_get_array(tendSum1stPool, 'normalVelocity', normalVelocityTendSum1st)
       call mpas_pool_get_array(tendSum1stPool, 'layerThickness', layerThicknessTendSum1st)
       call mpas_pool_get_subpool(tendSum1stPool, 'tracersTend', tracersTendSum1stPool)

       call mpas_pool_get_array(tendSum2ndPool, 'normalVelocity', normalVelocityTendSum2nd)
       call mpas_pool_get_array(tendSum2ndPool, 'layerThickness', layerThicknessTendSum2nd)
       call mpas_pool_get_subpool(tendSum2ndPool, 'tracersTend', tracersTendSum2ndPool)

       call mpas_pool_get_array(tendSum3rdPool, 'normalVelocity', normalVelocityTendSum3rd)
       call mpas_pool_get_array(tendSum3rdPool, 'layerThickness', layerThicknessTendSum3rd)
       call mpas_pool_get_subpool(tendSum3rdPool, 'tracersTend', tracersTendSum3rdPool)

       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)

       call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

       call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
       forcingTimeIncrement = 0.0_RKIND

       !$omp do schedule(runtime) private(k)
       do iEdge = 1, nEdges
          do k = 1, maxLevelEdgeTop(iEdge)
             normalVelocityFirstStage(k, iEdge) = normalVelocityCur(k, iEdge)
             normalVelocitySecondStage(k, iEdge) = normalVelocityCur(k, iEdge)
             normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
          end do
       end do
       !$omp end do

       !$omp do schedule(runtime) private(k)
       do iCell = 1, nCells
          do k = 1, maxLevelCell(iCell)
             layerThicknessFirstStage(k, iCell) = layerThicknessCur(k, iCell)
             layerThicknessSecondStage(k, iCell) = layerThicknessCur(k, iCell)
             layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
          end do
       end do
       !$omp end do

       normalVelocityTendSum1st(:,:) = 0.0_RKIND
       layerThicknessTendSum1st(:,:) = 0.0_RKIND

       normalVelocityTendSum2nd(:,:) = 0.0_RKIND
       layerThicknessTendSum2nd(:,:) = 0.0_RKIND

       normalVelocityTendSum3rd(:,:) = 0.0_RKIND
       layerThicknessTendSum3rd(:,:) = 0.0_RKIND

       call mpas_pool_begin_iteration(tracersPool)
       do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
          if ( groupItr % memberType == MPAS_POOL_FIELD ) then
             configName = 'config_use_' // trim(groupItr % memberName)
             call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

             if ( config_use_tracerGroup ) then
                modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                call mpas_pool_get_array(tracersTendSum1stPool, modifiedGroupName, tracersGroupTendSum1st)
                call mpas_pool_get_array(tracersTendSum2ndPool, modifiedGroupName, tracersGroupTendSum2nd)
                call mpas_pool_get_array(tracersTendSum3rdPool, modifiedGroupName, tracersGroupTendSum3rd)
                tracersGroupTendSum1st(:,:,:) = 0.0_RKIND
                tracersGroupTendSum2nd(:,:,:) = 0.0_RKIND
                tracersGroupTendSum3rd(:,:,:) = 0.0_RKIND

                call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupCur, 1)
                call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupNew, 2)
                call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupFirstStage, 3)
                call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupSecondStage, 4)
                !$omp do schedule(runtime) private(k)
                do iCell = 1, nCells
                   do k = 1, maxLevelCell(iCell)
                      tracersGroupFirstStage(:, k, iCell) = tracersGroupCur(:, k, iCell)
                      tracersGroupSecondStage(:, k, iCell) = tracersGroupCur(:, k, iCell)
                      tracersGroupNew(:, k, iCell) = tracersGroupCur(:, k, iCell)
                   end do
                end do
                !$omp end do
             end if
          end if
       end do

       block => block % next
    end do

    block => domain % blocklist
    do while(associated(block))
       if (associated(block % prev)) then
          call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_1st', tendSum1stPool)
          call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_2nd', tendSum2ndPool)         
          call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_3rd', tendSum3rdPool)         
          call mpas_pool_get_subpool(block % prev % structs, 'tend_slow', tendSlowPool)         
       else
          nullify(prevTendSum1stPool)
          nullify(prevTendSum2ndPool)
          nullify(prevTendSum3rdPool)
          nullify(prevTendSlowPool)
       end if

       if (associated(block % next)) then
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_1st', nextTendSum1stPool)
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_2nd', nextTendSum2ndPool)
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_3rd', nextTendSum3rdPool)
            call mpas_pool_get_subpool(block % next % structs, 'tend_slow', nextTendSlowPool)
       else
          nullify(nextTendSum1stPool)
          nullify(nextTendSum2ndPool)
          nullify(nextTendSum3rdPool)
          nullify(nextTendSlowPool)
       end if

       call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
       call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
       call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
       call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)

       if (associated(prevTendSum1stPool) .and. associated(nextTendSum1stPool)) then
          call mpas_pool_link_pools(tendSum1stPool, prevTendSum1stPool, nextTendSum1stPool)
       else if (associated(prevTendSum1stPool)) then
          call mpas_pool_link_pools(tendSum1stPool, prevTendSum1stPool)
       else if (associated(nextTendSum1stPool)) then
          call mpas_pool_link_pools(tendSum1stPool,nextPool=nextTendSum1stPool)
       else
          call mpas_pool_link_pools(tendSum1stPool)
       end if

       if (associated(prevTendSum2ndPool) .and. associated(nextTendSum2ndPool)) then
          call mpas_pool_link_pools(tendSum2ndPool, prevTendSum2ndPool, nextTendSum2ndPool)
       else if (associated(prevTendSum2ndPool)) then
          call mpas_pool_link_pools(tendSum2ndPool, prevTendSum2ndPool)
       else if (associated(nextTendSum2ndPool)) then
          call mpas_pool_link_pools(tendSum2ndPool,nextPool=nextTendSum2ndPool)
       else
          call mpas_pool_link_pools(tendSum2ndPool)
       end if

       if (associated(prevTendSum3rdPool) .and. associated(nextTendSum3rdPool)) then
          call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool, nextTendSum3rdPool)
       else if (associated(prevTendSum3rdPool)) then
          call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool)
       else if (associated(nextTendSum3rdPool)) then
          call mpas_pool_link_pools(tendSum3rdPool,nextPool=nextTendSum3rdPool)
       else
          call mpas_pool_link_pools(tendSum3rdPool)
       end if

       if (associated(prevTendSlowPool) .and. associated(nextTendSlowPool)) then
          call mpas_pool_link_pools(tendSlowPool, prevTendSlowPool, nextTendSlowPool)
       else if (associated(prevTendSlowPool)) then
          call mpas_pool_link_pools(tendSlowPool, prevTendSlowPool)
       else if (associated(nextTendSlowPool)) then
          call mpas_pool_link_pools(tendSlowPool,nextPool=nextTendSlowPool)
       else
          call mpas_pool_link_pools(tendSlowPool)
       end if

       call mpas_pool_link_parinfo(block, tendSum1stPool)
       call mpas_pool_link_parinfo(block, tendSum2ndPool)
       call mpas_pool_link_parinfo(block, tendSum3rdPool)
       call mpas_pool_link_parinfo(block, tendSlowPool)

       block => block % next
    end do

    call mpas_timer_stop("LTS3-init variables and create extra pools for LTS")

    call mpas_timer_start("LTS3-main loop")

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! BEGIN LTS SCHEME
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! -------------------------------------- IF SPLITTING IS ON  --------------------------------------
      ! -------------------------------- COMPUTE SLOW TEND FOR EVERYBODY  -------------------------------
      if (config_LTS_split) then

         ! set configs to compute only the slow terms on the momentum equation
         call setConfigsForSlowTends(domain)

         call mpas_timer_start("LTS3-slow tend computation")

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
            call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
            forcingTimeIncrement = 0.0_RKIND

            ! TENDENCIES COMPUTATION --- (this is only for the momentum equation)
            call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
            call ocn_tend_vel(tendSlowPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1, dt)

            block => block % next
         end do

         ! set configs to compute only the fast terms on the momentum equation
         call setConfigsForFastTends(domain)

         call mpas_timer_stop("LTS3-slow tend computation")
      end if
     ! --------------------------------------  STEP  ------------------------------------------
     ! --- compute the first stage of SSPRK for interface layer 1, interface layer 2 and coarse

      ! --- update halos for diagnostic variables
      if (config_use_cvmix_kpp) then
         call mpas_timer_start("LTS3-boundary layer depth halo update")
         call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
         call mpas_timer_stop("LTS3-boundary layer depth halo update")
      end if

      call mpas_timer_start("LTS3-diagnostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
      if (config_mom_del4 > 0.0_RKIND) then
         call mpas_dmpar_field_halo_exch(domain, 'divergence')
         call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
      end if
      call mpas_timer_stop("LTS3-diagnostic halo update")

      ! --- compute tendencies for velocity, thickness, and tracers
      call mpas_timer_start("LTS3-compute tendencies vel/thick/tracer")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
         call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
         forcingTimeIncrement = 0.0_RKIND
         
         ! VELOCITY: VERTICAL TRANSPORT COMPUTATION --- 
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, layerThicknessCur, &
                                              layerThicknessEdge, normalVelocityCur, sshCur, dt, vertAleTransportTop, err)

         ! VELOCITY AND THICKNESS: TENDENCY COMPUTATION ---
         if (.not. config_LTS_split) then
            call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
         end if
         if ((.not. config_LTS_inline_fast_tends) .or. (.not. config_LTS_split)) then
            call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1, dt)
            call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)
         else if (config_LTS_split .and. config_LTS_inline_fast_tends) then
            call ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, 1, 0, 1, 1, 1)
         end if

         ! TRACERS: TENDENCY COMPUTATION ---
         call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dt, activeTracersOnlyIn, 1)

         block => block % next
      end do
      call mpas_timer_stop("LTS3-compute tendencies vel/thick/tracer")

      call mpas_timer_start("LTS3-advance soln")
      ! --- advance solution for interface layer 1, interface layer 2 and coarse with first stage of SSPRK
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3) !time level 3 is first stage
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3) !time level 3 is first stage

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
         call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

         if (config_LTS_split) then

            normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

         end if

         ! --- NORMAL VELOCITY

         ! --- interface layers
         do iRegion =1,nRegions
            do ie = 1, nEdgesInLTSHalo(iRegion,2)
               iEdge = edgesInLTSHalo(iRegion,2,ie)
               normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge)
            end do
         end do
         ! --- coarse
         do ie = 1, nEdgesInLTSHalo(2,1)
            iEdge = edgesInLTSHalo(2,1,ie)
            normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge)
         end do
         ! --- fine layers close to interface layers
         do ie = 1, nEdgesInLTSHalo(1,3)
            iEdge = edgesInLTSHalo(1,3,ie)
            normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge)
         end do
    
 
         ! --- LAYER THICKNESS

         ! --- interface layers
         do iRegion =1,nRegions
            do ic = 1, nCellsInLTSHalo(iRegion,2)
               iCell = cellsInLTSHalo(iRegion,2,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
               end do
            end do
         end do
         ! --- coarse
         do ic = 1, nCellsInLTSHalo(2,1)
            iCell = cellsInLTSHalo(2,1,ic)
            do k = 1, maxLevelCell(iCell)
               layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
            end do
         end do
         ! --- fine layers close to interface layers
         do ic = 1, nCellsInLTSHalo(1,3)
            iCell = cellsInLTSHalo(1,3,ic)
            do k = 1, maxLevelCell(iCell)
               layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
            end do
         end do


         ! --- TRACERS

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               configName = 'config_use_' // trim(groupItr % memberName)
               call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

               if ( config_use_tracerGroup ) then
                  call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupFirstStage, 3) !time level 3 is first stage
                  call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

                  modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                  call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                  !$omp do schedule(runtime) private(k)
                  ! --- interface layers
                  do iRegion =1,nRegions
                     do ic = 1, nCellsInLTSHalo(iRegion,2)
                        iCell = cellsInLTSHalo(iRegion,2,ic)
                        do k = 1, maxLevelCell(iCell)
                            tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) & 
                                                                  +   dt * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell)
                        end do
                     end do
                  end do
                  ! --- coarse
                  do ic = 1, nCellsInLTSHalo(2,1)
                     iCell = cellsInLTSHalo(2,1,ic)
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                +   dt * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell)
                     end do
                  end do
                  ! --- fine layers close to interface layers
                  do ic = 1, nCellsInLTSHalo(1,3)
                     iCell = cellsInLTSHalo(1,3,ic)
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                +   dt * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell)
                     end do
                  end do
                  !$omp end do
               end if
            end if
         end do

         block => block % next
      end do
      call mpas_timer_stop("LTS3-advance soln")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      
      call mpas_timer_start("LTS3-vel/thick prognostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=3)
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=3)
      call mpas_timer_stop("LTS3-vel/thick prognostic halo update")
      
      
      call mpas_timer_start("LTS3-tracer prognostic halo update")
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=3)
         end if
      end do
      call mpas_timer_stop("LTS3-tracer prognostic halo update")

      ! --------------------------------------  STEP  ------------------------------------------
      ! --- compute the second stage of SSPRK for interface layer 1, interface layer 2 and coarse

      ! --- update halos for diagnostic variables.
      if (config_use_cvmix_kpp) then
         call mpas_timer_start("LTS3-boundary layer depth halo update")
         call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
         call mpas_timer_stop("LTS3-boundary layer depth halo update")
      end if

      call mpas_timer_start("LTS3-diagnostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
      if (config_mom_del4 > 0.0_RKIND) then
         call mpas_dmpar_field_halo_exch(domain, 'divergence')
         call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
      end if
      call mpas_timer_stop("LTS3-diagnostic halo update")

      ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
      call mpas_timer_start("LTS3-update diagnostic variables")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
         call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityFirstStage(:, iEdge)
         end do
         !$omp end do

         if(haveCoarse == 1 .or. haveInterface == 1) then

            ! DIAGNOSTICS UPDATE ---

            if ( (.not. config_LTS_split) .or. (.not. config_LTS_inline_fast_tends) ) then
               call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 3)
            end if

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-update diagnostic variables")

      
      call mpas_timer_start("LTS3-compute tendencies vel/thick/tracer")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
         call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)
         call mpas_pool_get_array(statePool, 'ssh', sshFirstStage, 3)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
        
         call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
         forcingTimeIncrement = dt
 
         if(haveCoarse == 1 .or. haveInterface == 1) then


            ! VERTICAL TRANSPORT COMPUTATION ---
            call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, layerThicknessFirstStage, &
                                              layerThicknessEdge, normalTransportVelocity, sshFirstStage, dt, vertAleTransportTop, err)

            ! TENDENCIES COMPUTATION ---
            if (.not. config_LTS_split) then
               call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
            end if
            if ((.not. config_LTS_inline_fast_tends) .or. (.not. config_LTS_split)) then
               call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 3, dt)
               call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)
            else if (config_LTS_split .and. config_LTS_inline_fast_tends) then
               call ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, 3, 0, 0, 1, 1)
            end if
            call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dt, activeTracersOnlyIn, 3)
           
         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-compute tendencies vel/thick/tracer")

      call mpas_timer_start("LTS3-advance soln")
      ! --- advance solution for interface layer 1, interface layer 2 and coarse with second stage of SSPRK
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4) !time level 4 is second stage
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4) !time level 4 is second stage

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
         call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

         if(haveCoarse == 1 .or. haveInterface == 1) then

            if (config_LTS_split) then

               normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

            end if
            ! --- NORMAL VELOCITY

            ! --- coarse
            do ie = 1, nEdgesInLTSHalo(2,1)
               iEdge = edgesInLTSHalo(2,1,ie)
               normalVelocitySecondStage(:,iEdge) = weightOld * normalVelocityCur(:,iEdge) &
                                                  + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                  + weightTend * dt * normalVelocityTend(:,iEdge)
            end do

            ! --- interface layers
            do iRegion = 1, nRegions
               do ie = 1, nEdgesInLTSHalo(iRegion,2)
                  iEdge = edgesInLTSHalo(iRegion,2,ie)
                  normalVelocitySecondStage(:,iEdge) = weightOld * normalVelocityCur(:,iEdge) &
                                                  + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                  + weightTend * dt * normalVelocityTend(:,iEdge)
               end do
            end do

            ! --- LAYER THICKNESS

            ! --- coarse
            do ic = 1, nCellsInLTSHalo(2,1)
               iCell = cellsInLTSHalo(2,1,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) & 
                                                     + weightNew * layerThicknessFirstStage(k,iCell) &
                                                     + weightTend * dt * layerThicknessTend(k,iCell)
               end do
            end do

            ! --- interface layers
            do iRegion = 1, nRegions
               do ic = 1, nCellsInLTSHalo(iRegion,2)
                  iCell = cellsInLTSHalo(iRegion,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) &
                                                        + weightNew * layerThicknessFirstStage(k,iCell) &
                                                        + weightTend * dt * layerThicknessTend(k,iCell)
                  end do
               end do
            end do

            ! --- TRACERS

            call mpas_pool_begin_iteration(tracersPool)
            do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
               if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                  configName = 'config_use_' // trim(groupItr % memberName)
                  call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)
                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupFirstStage, 3)
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupSecondStage, 4)
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

                     modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                     !$omp do schedule(runtime) private(k)
                     ! --- coarse
                     do ic = 1, nCellsInLTSHalo(2,1)
                        iCell = cellsInLTSHalo(2,1,ic)
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupSecondStage(:,k,iCell) = ( weightOld * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                  + weightNew * tracersGroupFirstStage(:,k,iCell) * layerThicknessFirstStage(k,iCell)  &
                                                                  + weightTend * dt * tracersGroupTend(:,k,iCell) ) / layerThicknessSecondStage(k,iCell)
                        end do
                     end do
                     ! --- interface layers
                     do iRegion =1,nRegions
                        do ic = 1, nCellsInLTSHalo(iRegion,2)
                           iCell = cellsInLTSHalo(iRegion,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupSecondStage(:,k,iCell) = ( weightOld * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                        + weightNew * tracersGroupFirstStage(:,k,iCell) * layerThicknessFirstStage(k,iCell)  &
                                                                        + weightTend * dt * tracersGroupTend(:,k,iCell) ) / layerThicknessSecondStage(k,iCell)
                           end do
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end do

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-advance soln")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      
      call mpas_timer_start("LTS3-vel/thick prognostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=4)
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=4)
      call mpas_timer_stop("LTS3-vel/thick prognostic halo update")
      
      call mpas_timer_start("LTS3-tracer prognostic halo update")
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=4)
         end if
      end do
      call mpas_timer_stop("LTS3-tracer prognostic halo update")

      ! --------------------------------------  STEP  ------------------------------------------
      ! --- substepping for fine soln advancement begins here. we also compute
      ! the correction terms for the interface to be used after the substepping

      do im =1, M

         call mpas_timer_start("LTS3-interface prediction")
         ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln cur in soln new
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

            if (haveInterface == 1) then
               alpha = REAL((im - 1)) / M
               alphaHat = REAL((im - 1) * (im - 1)) / (M * M)
               do ie = 1, nEdgesInLTSHalo(1,2) !interface layer 1
                  iEdge = edgesInLTSHalo(1,2,ie)
                  normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
                  normalVelocityCur(:,iEdge) = (1.0 - alpha - alphaHat) * normalVelocityCur(:,iEdge) + (alpha - alphaHat) * normalVelocityFirstStage(:,iEdge) &
                                             + 2.0 * alphaHat *normalVelocitySecondStage(:,iEdge)
               end do
               do ic = 1, nCellsInLTSHalo(1,2) !interface layer 1
                  iCell = cellsInLTSHalo(1,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
                     layerThicknessCur(k,iCell) = (1.0 - alpha - alphaHat) * layerThicknessCur(k,iCell) + (alpha - alphaHat) * layerThicknessFirstStage(k,iCell) &
                                               + 2.0 * alphaHat * layerThicknessSecondStage(k,iCell) 
                  end do
               end do

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupNew, 2)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupFirstStage, 3)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupSecondStage, 4)
                        !$omp do schedule(runtime) private(k)
                        do ic = 1, nCellsInLTSHalo(1,2)
                           iCell = cellsInLTSHalo(1,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupNew(:, k, iCell) = tracersGroupCur(:, k, iCell)
                              tracersGroupCur(:, k, iCell) = (1.0 - alpha - alphaHat) * tracersGroupCur(:,k,iCell) + &
                                                            (alpha - alphaHat) * tracersGroupFirstStage(:,k,iCell) + 2.0 * alphaHat * tracersGroupSecondStage(:,k,iCell)
                           end do
                        end do
                        !$omp end do
                     end if
                  end if
               end do

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-interface prediction")

         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         
         call mpas_timer_start("LTS3-vel/thick prognostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=1)
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=1)
         call mpas_timer_stop("LTS3-vel/thick prognostic halo update")

         call mpas_timer_start("LTS3-tracer prognostic halo update")
         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=1)
            end if
         end do
         call mpas_timer_stop("LTS3-tracer prognostic halo update")

         ! --- update halos for diagnostic variables
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("LTS3-boundary layer depth halo update")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("LTS3-boundary layer depth halo update")
         end if

         call mpas_timer_start("LTS3-diagnostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence')
            call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("LTS3-diagnostic halo update")


         ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
         call mpas_timer_start("LTS3-update diagnostic variables")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalVelocityCur(:, iEdge)
            end do
            !$omp end do

            if(haveFine == 1 .or. haveInterface == 1) then

               ! DIAGNOSTICS UPDATE ---

               if ( (.not. config_LTS_split) .or. (.not. config_LTS_inline_fast_tends) ) then
                  call ocn_diagnostic_solve(dtFine, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)
               end if

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-update diagnostic variables")

         call mpas_timer_start("LTS3-compute tendencies vel/thick/tracer")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
            call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
           
            call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
            forcingTimeIncrement = REAL((im-1))*dtFine + 0.0_RKIND

            if(haveFine == 1 .or. haveInterface == 1) then

               ! VERTICAL TRANSPORT COMPUTATION ---
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, layerThicknessCur, &
                                                    layerThicknessEdge, normalTransportVelocity, sshCur, dtFine, vertAleTransportTop, err)

               ! TENDENCIES COMPUTATION ---
               if (.not. config_LTS_split) then
                  call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
               end if
               if ((.not. config_LTS_inline_fast_tends) .or. (.not. config_LTS_split)) then
                  call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 1, dtFine)
                  call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)
               else if (config_LTS_split .and. config_LTS_inline_fast_tends) then
                  call ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, 1, 1, 1, 0, 1)
               end if
               call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dtFine, activeTracersOnlyIn, 1)

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-compute tendencies vel/thick/tracer")


         call mpas_timer_start("LTS3-advance soln")
         ! --- advance fine solution (first stage of SSPRK) AND sum up the tendecies above into tendSum2nd for interface layer 1 and interface
         ! layer 2 also restore soln cur by copying it back from soln new
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
            call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

            call mpas_pool_get_array(tendSum2ndPool, 'normalVelocity', normalVelocityTendSum2nd)
            call mpas_pool_get_array(tendSum2ndPool, 'layerThickness', layerThicknessTendSum2nd)

            call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

            if(haveFine == 1 .or. haveInterface == 1) then

               if (config_LTS_split) then

                  normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

               end if

               ! --- NORMAL VELOCITY

               ! --- interface layer 1 (correction)
               do ie = 1, nEdgesInLTSHalo(1,2)
                  iEdge = edgesInLTSHalo(1,2,ie)
                  normalVelocityCur(:,iEdge) = normalVelocityNew(:,iEdge)
                  normalVelocityTendSum2nd(:,iEdge) = normalVelocityTendSum2nd(:,iEdge) + normalVelocityTend(:,iEdge)
               end do
               ! --- interface layer 2 (correction)
               do ie = 1, nEdgesInLTSHalo(2,2)
                  iEdge = edgesInLTSHalo(2,2,ie)
                  normalVelocityTendSum2nd(:,iEdge) = normalVelocityTendSum2nd(:,iEdge) + normalVelocityTend(:,iEdge)
               end do
      
               ! --- fine (soln advancement)  
               do ie = 1, nEdgesInLTSHalo(1,1)
                  iEdge = edgesInLTSHalo(1,1,ie)
                  normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dtFine * normalVelocityTend(:,iEdge) !soln update for the fine
               end do
               do ie = 1, nEdgesInLTSHalo(1,3)
                  iEdge = edgesInLTSHalo(1,3,ie)
                  normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dtFine * normalVelocityTend(:,iEdge) !soln update for the fine
               end do
 
               ! --- LAYER THICKNESS

               ! --- interface layer 1 (correction)
               do ic = 1, nCellsInLTSHalo(1,2)
                  iCell = cellsInLTSHalo(1,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessCur(k,iCell) = layerThicknessNew(k,iCell)
                     layerThicknessTendSum2nd(k,iCell) = layerThicknessTendSum2nd(k,iCell) +  layerThicknessTend(k,iCell)
                  end do
               end do

               ! --- interface layer 2 (correction)
               do ic = 1, nCellsInLTSHalo(2,2)
                  iCell = cellsInLTSHalo(2,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessTendSum2nd(k,iCell) = layerThicknessTendSum2nd(k,iCell) +  layerThicknessTend(k,iCell)
                  end do
               end do

               ! --- fine (soln advancement)
               do ic = 1, nCellsInLTSHalo(1,1)
                  iCell = cellsInLTSHalo(1,1,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dtFine * layerThicknessTend(k,iCell) !soln update for the fine
                  end do
               end do
               do ic = 1, nCellsInLTSHalo(1,3)
                  iCell = cellsInLTSHalo(1,3,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dtFine * layerThicknessTend(k,iCell) !soln update for the fine
                  end do
               end do

               ! --- TRACERS

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)
                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupFirstStage, 3)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

                        modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                        call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                        call mpas_pool_get_array(tracersTendSum2ndPool, modifiedGroupName, tracersGroupTendSum2nd)
                        !$omp do schedule(runtime) private(k)

                        ! --- interface layer 1 (correction)
                        do ic = 1, nCellsInLTSHalo(1,2)
                           iCell = cellsInLTSHalo(1,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupCur(:,k,iCell) = tracersGroupNew(:,k,iCell)
                              tracersGroupTendSum2nd(:,k,iCell) = tracersGroupTendSum2nd(:,k,iCell) +  tracersGroupTend(:,k,iCell)
                           end do
                        end do
            
                        ! --- interface layer 2 (correction)
                        do ic = 1, nCellsInLTSHalo(2,2)
                           iCell = cellsInLTSHalo(2,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupTendSum2nd(:,k,iCell) = tracersGroupTendSum2nd(:,k,iCell) +  tracersGroupTend(:,k,iCell)
                           end do
                        end do

                        ! --- fine (soln advancement)
                        do ic = 1, nCellsInLTSHalo(1,1)
                           iCell = cellsInLTSHalo(1,1,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                    + dtFine * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell) !soln update for the fine
                           end do
                        end do
                        do ic = 1, nCellsInLTSHalo(1,3)
                           iCell = cellsInLTSHalo(1,3,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                    + dtFine * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell) !soln update for the fine
                           end do
                        end do

                        !$omp end do
                     end if
                  end if
               end do

            end if

            block => block % next
         end do

         call mpas_timer_stop("LTS3-advance soln")

         call mpas_timer_start("LTS3-interface prediction")
         ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln cur in soln new
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

            if (haveInterface == 1) then
               beta = REAL(im) / M
               betaHat = REAL((im - 1) * (im + 1)) / (M * M)
               do ie = 1, nEdgesInLTSHalo(1,2) !interface layer 1
                  iEdge = edgesInLTSHalo(1,2,ie)
                  normalVelocityNew(:,iEdge) = normalVelocityFirstStage(:,iEdge)
                  normalVelocityFirstStage(:,iEdge) = (1.0 - beta - betaHat) * normalVelocityCur(:,iEdge) + (beta - betaHat) * normalVelocityFirstStage(:,iEdge) &
                                             + 2.0 * betaHat *normalVelocitySecondStage(:,iEdge)
               end do
               do ic = 1, nCellsInLTSHalo(1,2) !interface layer 1
                  iCell = cellsInLTSHalo(1,2,ic)
                  do k = 1, maxLevelCell(iCell)
                    layerThicknessNew(k,iCell) = layerThicknessFirstStage(k,iCell)
                    layerThicknessFirstStage(k,iCell) = (1.0 - beta - betaHat) * layerThicknessCur(k,iCell) + (beta - betaHat) * layerThicknessFirstStage(k,iCell) &
                                                         + 2.0 * betaHat * layerThicknessSecondStage(k,iCell) 
                  end do
               end do

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupNew, 2)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupFirstStage, 3)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupSecondStage, 4)
                        !$omp do schedule(runtime) private(k)
                        do ic = 1, nCellsInLTSHalo(1,2)
                           iCell = cellsInLTSHalo(1,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupNew(:, k, iCell) = tracersGroupFirstStage(:, k, iCell)
                              tracersGroupFirstStage(:, k, iCell) = (1.0 - beta - betaHat) * tracersGroupCur(:,k,iCell) + &
                                                                    (beta - betaHat) * tracersGroupFirstStage(:,k,iCell) + 2.0 * betaHat * tracersGroupSecondStage(:,k,iCell)
                           end do
                        end do
                        !$omp end do
                     end if
                  end if
               end do

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-interface prediction")

         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool) 
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         
         call mpas_timer_start("LTS3-vel/thick prognostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=3)
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=3)
         call mpas_timer_stop("LTS3-vel/thick prognostic halo update")

         call mpas_timer_start("LTS3-tracer prognostic halo update")
         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=3)
            end if
         end do
         call mpas_timer_stop("LTS3-tracer prognostic halo update")

         ! --- update halos for diagnostic variables
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("LTS3-boundary layer depth halo update")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("LTS3-boundary layer depth halo update")
         end if

         call mpas_timer_start("LTS3-diagnostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence')
            call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("LTS3-diagnostic halo update")


         ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
         call mpas_timer_start("LTS3-update diagnostic variables")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalVelocityFirstStage(:, iEdge)
            end do
            !$omp end do

            if(haveFine == 1 .or. haveInterface == 1) then

               ! DIAGNOSTICS UPDATE ---

               if ( (.not. config_LTS_split) .or. (.not. config_LTS_inline_fast_tends) ) then
                  call ocn_diagnostic_solve(dtFine, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 3)
               end if

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-update diagnostic variables")

         call mpas_timer_start("LTS3-compute tendencies vel/thick/tracer")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)
            call mpas_pool_get_array(statePool, 'ssh', sshFirstStage, 3)
            call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
            call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
            
            call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
            forcingTimeIncrement = REAL((im-1))*dtFine + dtFine

            if(haveFine == 1 .or. haveInterface == 1) then

               ! VERTICAL TRANSPORT COMPUTATION ---
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, layerThicknessFirstStage, &
                                                    layerThicknessEdge, normalTransportVelocity, sshFirstStage, dtFine, vertAleTransportTop, err)
 
               ! TENDENCIES COMPUTATION ---
               if (.not. config_LTS_split) then
                  call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
               end if
               if ((.not. config_LTS_inline_fast_tends) .or. (.not. config_LTS_split)) then
                  call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 3, dtFine)
                  call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)
               else if (config_LTS_split .and. config_LTS_inline_fast_tends) then 
                  call ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, 3, 1, 1, 0, 1)
               end if
               call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dtFine, activeTracersOnlyIn, 3)

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-compute tendencies vel/thick/tracer")


         call mpas_timer_start("LTS3-advance soln")
         ! --- advance fine solution (second stage of SSPRK) AND sum up the tendecies above into tendSum1st for interface layer 1 and interface
         ! layer 2 also restore soln first stage by copying it back from soln new
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
            call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

            call mpas_pool_get_array(tendSum1stPool, 'normalVelocity', normalVelocityTendSum1st)
            call mpas_pool_get_array(tendSum1stPool, 'layerThickness', layerThicknessTendSum1st)

            call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

            if(haveFine == 1 .or. haveInterface == 1) then

               if (config_LTS_split) then

                  normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

               end if
    
               ! --- NORMAL VELOCITY

               ! --- interface layer 1 (correction)
               do ie = 1, nEdgesInLTSHalo(1,2)
                  iEdge = edgesInLTSHalo(1,2,ie)
                  normalVelocityFirstStage(:,iEdge) = normalVelocityNew(:,iEdge)
                  normalVelocityTendSum1st(:,iEdge) = normalVelocityTendSum1st(:,iEdge) + normalVelocityTend(:,iEdge)
               end do
               ! --- interface layer 2 (correction)
               do ie = 1, nEdgesInLTSHalo(2,2)
                  iEdge = edgesInLTSHalo(2,2,ie)
                  normalVelocityTendSum1st(:,iEdge) = normalVelocityTendSum1st(:,iEdge) + normalVelocityTend(:,iEdge)
               end do
      
               ! --- fine (soln advancement)  
               do ie = 1, nEdgesInLTSHalo(1,1)
                  iEdge = edgesInLTSHalo(1,1,ie)
                  normalVelocitySecondStage(:,iEdge) =  weightOld * normalVelocityCur(:,iEdge) + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                        + weightTend *  dtFine * normalVelocityTend(:,iEdge) !soln update for the fine
               end do
               do ie = 1, nEdgesInLTSHalo(1,3)
                  iEdge = edgesInLTSHalo(1,3,ie)
                  normalVelocitySecondStage(:,iEdge) =  weightOld * normalVelocityCur(:,iEdge) + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                        + weightTend *  dtFine * normalVelocityTend(:,iEdge) !soln update for the fine
               end do
 
               ! --- LAYER THICKNESS

               ! --- interface layer 1 (correction)
               do ic = 1, nCellsInLTSHalo(1,2)
                  iCell = cellsInLTSHalo(1,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessFirstStage(k,iCell) = layerThicknessNew(k,iCell)
                     layerThicknessTendSum1st(k,iCell) = layerThicknessTendSum1st(k,iCell) +  layerThicknessTend(k,iCell)
                  end do
               end do

               ! --- interface layer 2 (correction)
               do ic = 1, nCellsInLTSHalo(2,2)
                  iCell = cellsInLTSHalo(2,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessTendSum1st(k,iCell) = layerThicknessTendSum1st(k,iCell) +  layerThicknessTend(k,iCell)
                  end do
               end do

               ! --- fine (soln advancement)
               do ic = 1, nCellsInLTSHalo(1,1)
                  iCell = cellsInLTSHalo(1,1,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) + weightNew * layerThicknessFirstStage(k,iCell) &
                                                          + weightTend * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
                  end do
               end do
               do ic = 1, nCellsInLTSHalo(1,3)
                  iCell = cellsInLTSHalo(1,3,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) + weightNew * layerThicknessFirstStage(k,iCell) &
                                                          + weightTend * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
                  end do
               end do

               ! --- TRACERS

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)
                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupFirstStage, 3)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupSecondStage, 4)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
   
                        modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                        call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                        call mpas_pool_get_array(tracersTendSum1stPool, modifiedGroupName, tracersGroupTendSum1st)
                        !$omp do schedule(runtime) private(k)
   
                        ! --- interface layer 1 (correction)
                        do ic = 1, nCellsInLTSHalo(1,2)
                           iCell = cellsInLTSHalo(1,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupFirstStage(:,k,iCell) = tracersGroupNew(:,k,iCell)
                              tracersGroupTendSum1st(:,k,iCell) = tracersGroupTendSum1st(:,k,iCell) +  tracersGroupTend(:,k,iCell)
                           end do
                        end do
               
                        ! --- interface layer 2 (correction)
                        do ic = 1, nCellsInLTSHalo(2,2)
                           iCell = cellsInLTSHalo(2,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupTendSum1st(:,k,iCell) = tracersGroupTendSum1st(:,k,iCell) +  tracersGroupTend(:,k,iCell)
                           end do
                        end do

                        ! --- fine (soln advancement)
                        do ic = 1, nCellsInLTSHalo(1,1)
                           iCell = cellsInLTSHalo(1,1,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupSecondStage(:,k,iCell) = ( weightOld * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + &
                                                                     weightNew * tracersGroupFirstStage(:,k,iCell) * layerThicknessFirstStage(k,iCell) + &
                                                                     weightTend * dtFine * tracersGroupTend(:,k,iCell) ) / layerThicknessSecondStage(k,iCell)
                           end do
                        end do
                        do ic = 1, nCellsInLTSHalo(1,3)
                           iCell = cellsInLTSHalo(1,3,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupSecondStage(:,k,iCell) = ( weightOld * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + &
                                                                     weightNew * tracersGroupFirstStage(:,k,iCell) * layerThicknessFirstStage(k,iCell) + &
                                                                     weightTend * dtFine * tracersGroupTend(:,k,iCell) ) / layerThicknessSecondStage(k,iCell)
                           end do
                        end do

                        !$omp end do
                     end if
                  end if
               end do

            end if

            block => block % next
         end do

         call mpas_timer_stop("LTS3-advance soln")


         call mpas_timer_start("LTS3-interface prediction")
         ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln cur in soln new
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

            if (haveInterface == 1) then
               gam = REAL(2 * (im - 1) + 1) / ( 2 * M )
               gamHat = REAL(2 * (im - 1) * (im - 1) + 2 * (im - 1) + 1) / (2 * M * M)
               do ie = 1, nEdgesInLTSHalo(1,2) !interface layer 1
                  iEdge = edgesInLTSHalo(1,2,ie)
                  normalVelocityNew(:,iEdge) = normalVelocitySecondStage(:,iEdge)
                  normalVelocitySecondStage(:,iEdge) = (1.0 - gam - gamHat) * normalVelocityCur(:,iEdge) + (gam - gamHat) * normalVelocityFirstStage(:,iEdge) &
                                             + 2.0 * gamHat *normalVelocitySecondStage(:,iEdge)
               end do
               do ic = 1, nCellsInLTSHalo(1,2) !interface layer 1
                  iCell = cellsInLTSHalo(1,2,ic)
                  do k = 1, maxLevelCell(iCell)
                    layerThicknessNew(k,iCell) = layerThicknessSecondStage(k,iCell)
                    layerThicknessSecondStage(k,iCell) = (1.0 - gam - gamHat) * layerThicknessCur(k,iCell) + (gam - gamHat) * layerThicknessFirstStage(k,iCell) &
                                                         + 2.0 * gamHat * layerThicknessSecondStage(k,iCell) 
                  end do
               end do

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupNew, 2)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupFirstStage, 3)
                        call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupSecondStage, 4)
                        !$omp do schedule(runtime) private(k)
                        do ic = 1, nCellsInLTSHalo(1,2)
                           iCell = cellsInLTSHalo(1,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupNew(:, k, iCell) = tracersGroupSecondStage(:, k, iCell)
                              tracersGroupSecondStage(:, k, iCell) = (1.0 - gam - gamHat) * tracersGroupCur(:,k,iCell) + &
                                                                     (gam - gamHat) * tracersGroupFirstStage(:,k,iCell) + 2.0 * gamHat * tracersGroupSecondStage(:,k,iCell)
                           end do
                        end do
                        !$omp end do
                     end if
                  end if
               end do

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-interface prediction")

         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_timer_start("LTS3-vel/thick prognostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=4)
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=4)
         call mpas_timer_stop("LTS3-vel/thick prognostic halo update")

         call mpas_timer_start("LTS3-tracer prognostic halo update")
         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
              call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=4)
           end if
         end do
         call mpas_timer_stop("LTS3-tracer prognostic halo update")


         ! --- update halos for diagnostic variables
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("LTS3-boundary layer depth halo update")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("LTS3-boundary layer depth halo update")
         end if

         call mpas_timer_start("LTS3-diagnostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence')
            call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("LTS3-diagnostic halo update")


         ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
         call mpas_timer_start("LTS3-update diagnostic variables")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalVelocitySecondStage(:,iEdge)
            end do
            !$omp end do

            if(haveFine == 1 .or. haveInterface == 1) then

               ! DIAGNOSTICS UPDATE ---

               if ( (.not. config_LTS_split) .or. (.not. config_LTS_inline_fast_tends) ) then
                  call ocn_diagnostic_solve(dtFine, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 4)
               end if

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-update diagnostic variables")

         call mpas_timer_start("LTS3-compute tendencies vel/thick/tracer")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)
            call mpas_pool_get_array(statePool, 'ssh', sshSecondStage, 4)
            call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
            call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
           
            call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
            forcingTimeIncrement = REAL((im-1))*dtFine + dtFine/2.0_RKIND

           if(haveFine == 1 .or. haveInterface == 1) then

               ! VERTICAL TRANSPORT COMPUTATION ---
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, layerThicknessSecondStage, &
                                                    layerThicknessEdge, normalTransportVelocity, sshSecondStage, dtFine, vertAleTransportTop, err)
 
               ! TENDENCIES COMPUTATION ---
               if (.not. config_LTS_split) then
                  call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
               end if
               if ((.not. config_LTS_inline_fast_tends) .or. (.not. config_LTS_split)) then
                  call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 4, dtFine)
                  call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)
               else if (config_LTS_split .and. config_LTS_inline_fast_tends) then
                  call ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, 4, 1, 1, 0, 1)
               end if
               call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dtFine, activeTracersOnlyIn, 4)

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-compute tendencies vel/thick/tracer")


         call mpas_timer_start("LTS3-advance soln")
         ! --- advance fine solution (third stage of SSPRK) AND sum up the tendecies above into tendSum3rd for interface layer 1 and interface
         ! layer 2 also restore soln second stage by copying it back from soln new
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
            call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

            call mpas_pool_get_array(tendSum3rdPool, 'normalVelocity', normalVelocityTendSum3rd)
            call mpas_pool_get_array(tendSum3rdPool, 'layerThickness', layerThicknessTendSum3rd)

            call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityFirstStage, 3)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessFirstStage, 3)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

            if(haveFine == 1 .or. haveInterface == 1) then

               if (config_LTS_split) then

                  normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

               end if
    
               ! --- NORMAL VELOCITY

               ! --- interface layer 1 (correction)
               do ie = 1, nEdgesInLTSHalo(1,2)
                  iEdge = edgesInLTSHalo(1,2,ie)
                  normalVelocitySecondStage(:,iEdge) = normalVelocityNew(:,iEdge)
                  normalVelocityTendSum3rd(:,iEdge) = normalVelocityTendSum3rd(:,iEdge) + normalVelocityTend(:,iEdge)
               end do
               ! --- interface layer 2 (correction)
               do ie = 1, nEdgesInLTSHalo(2,2)
                  iEdge = edgesInLTSHalo(2,2,ie)
                  normalVelocityTendSum3rd(:,iEdge) = normalVelocityTendSum3rd(:,iEdge) + normalVelocityTend(:,iEdge)
               end do
      
               ! --- fine (soln advancement)  
               do ie = 1, nEdgesInLTSHalo(1,1)
                  iEdge = edgesInLTSHalo(1,1,ie)
                  normalVelocityCur(:,iEdge) =  (1.0 / 3.0) * normalVelocityCur(:,iEdge) + (2.0 / 3.0) * normalVelocitySecondStage(:,iEdge) &
                                                        + (2.0 / 3.0) *  dtFine * normalVelocityTend(:,iEdge) !soln update for the fine
               end do
               do ie = 1, nEdgesInLTSHalo(1,3)
                  iEdge = edgesInLTSHalo(1,3,ie)
                  normalVelocityCur(:,iEdge) =  (1.0 / 3.0) * normalVelocityCur(:,iEdge) + (2.0 / 3.0) * normalVelocitySecondStage(:,iEdge) &
                                                        + (2.0 / 3.0) *  dtFine * normalVelocityTend(:,iEdge) !soln update for the fine
               end do
 
               ! --- LAYER THICKNESS

               ! --- interface layer 1 (correction)
               do ic = 1, nCellsInLTSHalo(1,2)
                  iCell = cellsInLTSHalo(1,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessSecondStage(k,iCell) = layerThicknessNew(k,iCell)
                     layerThicknessTendSum3rd(k,iCell) = layerThicknessTendSum3rd(k,iCell) +  layerThicknessTend(k,iCell)
                  end do
               end do

               ! --- interface layer 2 (correction)
               do ic = 1, nCellsInLTSHalo(2,2)
                  iCell = cellsInLTSHalo(2,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessTendSum3rd(k,iCell) = layerThicknessTendSum3rd(k,iCell) +  layerThicknessTend(k,iCell)
                  end do
               end do

               ! --- fine (soln advancement)
               do ic = 1, nCellsInLTSHalo(1,1)
                  iCell = cellsInLTSHalo(1,1,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessNew(k,iCell) = (1.0 / 3.0)  * layerThicknessCur(k,iCell) + (2.0 / 3.0) * layerThicknessSecondStage(k,iCell) &
                                                          + (2.0 / 3.0) * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
                  end do
               end do
               do ic = 1, nCellsInLTSHalo(1,3)
                  iCell = cellsInLTSHalo(1,3,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessNew(k,iCell) = (1.0 / 3.0)  * layerThicknessCur(k,iCell) + (2.0 / 3.0) * layerThicknessSecondStage(k,iCell) &
                                                          + (2.0 / 3.0) * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
                  end do
               end do

               ! --- TRACERS

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)
                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupFirstStage, 3)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupSecondStage, 4)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
   
                        modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                        call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                        call mpas_pool_get_array(tracersTendSum3rdPool, modifiedGroupName, tracersGroupTendSum3rd)
                        !$omp do schedule(runtime) private(k)
   
                        ! --- interface layer 1 (correction)
                        do ic = 1, nCellsInLTSHalo(1,2)
                           iCell = cellsInLTSHalo(1,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupSecondStage(:,k,iCell) = tracersGroupNew(:,k,iCell)
                              tracersGroupTendSum3rd(:,k,iCell) = tracersGroupTendSum3rd(:,k,iCell) +  tracersGroupTend(:,k,iCell)
                           end do
                        end do
               
                        ! --- interface layer 2 (correction)
                        do ic = 1, nCellsInLTSHalo(2,2)
                           iCell = cellsInLTSHalo(2,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupTendSum3rd(:,k,iCell) = tracersGroupTendSum3rd(:,k,iCell) +  tracersGroupTend(:,k,iCell)
                           end do
                        end do

                        ! --- fine (soln advancement)
                        do ic = 1, nCellsInLTSHalo(1,1)
                           iCell = cellsInLTSHalo(1,1,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupCur(:,k,iCell) = ( (1.0 / 3.0) * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + &
                                                             (2.0 / 3.0) * tracersGroupSecondStage(:,k,iCell) * layerThicknessSecondStage(k,iCell) + &
                                                             (2.0 / 3.0) * dtFine * tracersGroupTend(:,k,iCell) ) / layerThicknessNew(k,iCell)
                           end do
                        end do
                        do ic = 1, nCellsInLTSHalo(1,3)
                           iCell = cellsInLTSHalo(1,3,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupCur(:,k,iCell) = ( (1.0 / 3.0) * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + &
                                                             (2.0 / 3.0) * tracersGroupSecondStage(:,k,iCell) * layerThicknessSecondStage(k,iCell) + &
                                                             (2.0 / 3.0) * dtFine * tracersGroupTend(:,k,iCell) ) / layerThicknessNew(k,iCell)
                           end do
                        end do

                        !$omp end do
                     end if
                  end if
               end do

               do ic = 1, nCellsInLTSHalo(1,1)
                  iCell = cellsInLTSHalo(1,1,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessCur(k,iCell) = layerThicknessNew(k,iCell)
                  end do
               end do
               do ic = 1, nCellsInLTSHalo(1,3)
                  iCell = cellsInLTSHalo(1,3,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessCur(k,iCell) = layerThicknessNew(k,iCell) 
                  end do
               end do

            end if

            block => block % next
         end do

         call mpas_timer_stop("LTS3-advance soln")


         if ( (.not. (config_disable_vel_vmix .or. config_disable_vel_all_tend) ) .or. & 
              (.not. (config_disable_tr_vmix .or. config_disable_tr_all_tend) ) ) then

            ! --- TODO: check if this halo update is needed or not
            call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

            call mpas_timer_start("LTS3-vel/thick prognostic halo update")
            call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=1)
            call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=1)
            call mpas_timer_stop("LTS3-vel/thick prognostic halo update")

            call mpas_timer_start("LTS3-tracer prognostic halo update")
            call mpas_pool_begin_iteration(tracersPool)
            do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
               if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                 call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=1)
              end if
            end do
            call mpas_timer_stop("LTS3-tracer prognostic halo update")

            ! --- update halos for diagnostic variables
            if (config_use_cvmix_kpp) then
               call mpas_timer_start("LTS3-boundary layer depth halo update")
               call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
               call mpas_timer_stop("LTS3-boundary layer depth halo update")
            end if

            call mpas_timer_start("LTS3-diagnostic halo update")
            call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
            if (config_mom_del4 > 0.0_RKIND) then
               call mpas_dmpar_field_halo_exch(domain, 'divergence')
              call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
            end if
            call mpas_timer_stop("LTS3-diagnostic halo update")

            ! --- compute diagnostics for implicit vert diff for velocity, thickness, and tracers
            call mpas_timer_start("LTS3-update diagnostic variables for vert diff")
            block => domain % blocklist
            do while (associated(block))

               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
               call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
               call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

               call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
               call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
               call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)

               !$omp do schedule(runtime)
               do iEdge = 1, nEdges
                  normalTransportVelocity(:, iEdge) = normalVelocityCur(:, iEdge)
               end do
               !$omp end do

               if (haveFine == 1) then

                  ! DIAGNOSTICS UPDATE ---

                  call ocn_diagnostic_solve(dtFine, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)

               end if

               block => block % next
            end do
            call mpas_timer_stop("LTS3-update diagnostic variables for vert diff")

         end if

         call mpas_timer_start("LTS3-compute vmix for fine")
         !TODO: this has not been tested yet
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)

            if (haveFine == 1) then
               call ocn_vmix_implicit(dtFine, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 1)
            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-compute vmix for fine")

      end do ! --- end substepping
 

      call mpas_timer_start("LTS3-copy soln into new soln")
      ! --- copy old soln into new soln for the fine cells and edges
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)

         if(haveFine == 1) then


            ! --- normal velocity
            do ie = 1, nEdgesInLTSHalo(1,1) 
               iEdge = edgesInLTSHalo(1,1,ie)
               normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
            end do
            do ie = 1, nEdgesInLTSHalo(1,3)
               iEdge = edgesInLTSHalo(1,3,ie)
               normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
            end do

            ! --- layer thickness
            do ic = 1, nCellsInLTSHalo(1,1) 
               iCell = cellsInLTSHalo(1,1,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
               end do
            end do
            do ic = 1, nCellsInLTSHalo(1,3) 
               iCell = cellsInLTSHalo(1,3,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
               end do
            end do

            ! --- tracers
            call mpas_pool_begin_iteration(tracersPool)
            do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
               if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                  configName = 'config_use_' // trim(groupItr % memberName)
                  call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then

                     call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupCur, 1)
                     call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupNew, 2)
                     !$omp do schedule(runtime) private(k)
                     do ic = 1, nCellsInLTSHalo(1,1)
                        iCell = cellsInLTSHalo(1,1,ic)
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupNew(:, k, iCell) = tracersGroupCur(:, k, iCell)
                        end do
                     end do
                     do ic = 1, nCellsInLTSHalo(1,3)
                        iCell = cellsInLTSHalo(1,3,ic)
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupNew(:, k, iCell) = tracersGroupCur(:, k, iCell)
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end do

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-copy soln into new soln")

      ! --------------------------------------  STEP  ------------------------------------------
      ! --- compute the third stage of SSPRK for coarse
      !

      ! --- update halos for diagnostic variables
      if (config_use_cvmix_kpp) then
         call mpas_timer_start("LTS3-boundary layer depth halo update")
         call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
         call mpas_timer_stop("LTS3-boundary layer depth halo update")
      end if

      call mpas_timer_start("LTS3-diagnostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
      if (config_mom_del4 > 0.0_RKIND) then
         call mpas_dmpar_field_halo_exch(domain, 'divergence')
         call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
      end if
      call mpas_timer_stop("LTS3-diagnostic halo update")

      ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
      call mpas_timer_start("LTS3-update diagnostic variables")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool) 
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

         ! TODO: I THINK THIS CAN BE REMOVED
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocitySecondStage(:,iEdge)
         end do
         !$omp end do

         if(haveCoarse == 1) then

            ! DIAGNOSTICS UPDATE ---
            if ( (.not. config_LTS_split) .or. (.not. config_LTS_inline_fast_tends) ) then
               call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 4)
            end if

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-update diagnostic variables")

      call mpas_timer_start("LTS3-compute tendencies vel/thick/tracer")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)
         call mpas_pool_get_array(statePool, 'ssh', sshSecondStage, 4)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
           
         call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
         forcingTimeIncrement = dt/2.0_RKIND

         if(haveCoarse == 1) then

            ! VERTICAL TRANSPORT COMPUTATION ---
            call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, layerThicknessSecondStage, &
                                              layerThicknessEdge, normalTransportVelocity, sshSecondStage, dt, vertAleTransportTop, err)

            ! TENDENCIES COMPUTATION ---
            if (.not. config_LTS_split) then
               call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
            end if
            if ((.not. config_LTS_inline_fast_tends) .or. (.not. config_LTS_split)) then
                  call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, 4, dt)
                  call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)
            else if (config_LTS_split .and. config_LTS_inline_fast_tends) then
                  call ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, 4, 0, 0, 1, 0)
            end if
            call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dt, activeTracersOnlyIn, 4)

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-compute tendencies vel/thick/tracer")

      call mpas_timer_start("LTS3-advance soln")
      ! --- advance solution for coarse with third stage of SSPRK
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocitySecondStage, 4)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessSecondStage, 4)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
         call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

         call mpas_pool_get_array(diagnosticsPool, 'wettingVelocity', wettingVelocity)

         if (haveCoarse == 1) then
            if (config_LTS_split) then

               normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

            end if 
            ! --- NORMAL VELOCITY

            ! --- coarse
            do ie = 1, nEdgesInLTSHalo(2,1)
               iEdge = edgesInLTSHalo(2,1,ie)
               normalVelocityNew(:,iEdge) = (1.0 / 3.0)  * normalVelocityCur(:,iEdge) + (2.0 / 3.0) * normalVelocitySecondStage(:,iEdge) &
                                            + (2.0 / 3.0) * dt * normalVelocityTend(:,iEdge)
            end do

            ! --- LAYER THICKNESS

            ! --- coarse
            do ic = 1, nCellsInLTSHalo(2,1)
               iCell = cellsInLTSHalo(2,1,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessNew(k,iCell) = (1.0 / 3.0)  * layerThicknessCur(k,iCell) + (2.0 / 3.0) * layerThicknessSecondStage(k,iCell) &
                                               + (2.0 / 3.0) * dt * layerThicknessTend(k,iCell)
               end do
            end do

            ! --- TRACERS

            call mpas_pool_begin_iteration(tracersPool)
            do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
               if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                  configName = 'config_use_' // trim(groupItr % memberName)
                  call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupSecondStage, 4)
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

                     modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                     !$omp do schedule(runtime) private(k)
                     ! --- coarse
                     do ic = 1, nCellsInLTSHalo(2,1)
                        iCell = cellsInLTSHalo(2,1,ic)
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupNew(:,k,iCell) = ( (1.0 / 3.0) * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + &
                                                          (2.0 / 3.0) * tracersGroupSecondStage(:,k,iCell) * layerThicknessSecondStage(k,iCell) + &
                                                          (2.0 / 3.0) * dt * tracersGroupTend(:,k,iCell) ) / layerThicknessNew(k,iCell)
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end do

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-advance soln")

      ! --------------------------------------  STEP  ------------------------------------------
      ! --- interface correction
      !

      call mpas_timer_start("LTS3-interface correction")
      ! --- correct solution for interface layer 1 and interface layer 2
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
         call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
         call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

         call mpas_pool_get_array(tendSum1stPool, 'normalVelocity', normalVelocityTendSum1st)
         call mpas_pool_get_array(tendSum1stPool, 'layerThickness', layerThicknessTendSum1st)

         call mpas_pool_get_array(tendSum2ndPool, 'normalVelocity', normalVelocityTendSum2nd)
         call mpas_pool_get_array(tendSum2ndPool, 'layerThickness', layerThicknessTendSum2nd)

         call mpas_pool_get_array(tendSum3rdPool, 'normalVelocity', normalVelocityTendSum3rd)
         call mpas_pool_get_array(tendSum3rdPool, 'layerThickness', layerThicknessTendSum3rd)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         if (haveInterface == 1) then

            do iRegion = 1,nRegions

               ! --- NORMAL VELOCITY

               ! --- interface layers

               do ie = 1, nEdgesInLTSHalo(iRegion,2)
                  iEdge = edgesInLTSHalo(iRegion,2,ie)
                  normalVelocityNew(:,iEdge) =  normalVelocityCur(:,iEdge) + weightTendSum2nd * dtFine * normalVelocityTendSum2nd(:,iEdge) &
                                                + weightTendSum1st * dtFine * normalVelocityTendSum1st(:,iEdge) + weightTendSum3rd * dtFine * normalVelocityTendSum3rd(:,iEdge)
               end do

               ! --- LAYER THICKNESS

               ! --- interface layers

               do ic = 1, nCellsInLTSHalo(iRegion,2)
                  iCell = cellsInLTSHalo(iRegion,2,ic)
                  do k= 1, maxLevelCell(iCell)
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) + weightTendSum2nd * dtFine * layerThicknessTendSum2nd(k,iCell) &
                                                  + weightTendSum1st * dtFine * layerThicknessTendSum1st(k,iCell) &
                                                  + weightTendsum3rd * dtFine * layerThicknessTendSum3rd(k,iCell)
                  end do
               end do
            end do
   
            call mpas_pool_begin_iteration(tracersPool)
            do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
               if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                  configName = 'config_use_' // trim(groupItr % memberName)
                  call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

                     modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array(tracersTendSum1stPool, modifiedGroupName, tracersGroupTendSum1st)
                     call mpas_pool_get_array(tracersTendSum2ndPool, modifiedGroupName, tracersGroupTendSum2nd)
                     call mpas_pool_get_array(tracersTendSum3rdPool, modifiedGroupName, tracersGroupTendSum3rd)

                     !$omp do schedule(runtime) private(k)

                     ! --- TRACERS

                     ! --- interface layers

                     do iRegion = 1, nRegions
                        do ic = 1, nCellsInLTSHalo(iRegion,2)
                           iCell = cellsInLTSHalo(iRegion,2,ic)
                           do k= 1, maxLevelCell(iCell)
                              tracersGroupNew(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                             + weightTendSum2nd * dtFine * tracersGroupTendSum2nd(:,k,iCell) &
                                                             + weightTendSum1st * dtFine * tracersGroupTendSum1st(:,k,iCell) &
                                                             + weightTendsum3rd * dtFine * tracersgroupTendSum3rd(:,k,iCell) ) / layerThicknessNew(k,iCell)
                           end do
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end do

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-interface correction")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_timer_start("LTS3-vel/thick prognostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=2)
      call mpas_timer_stop("LTS3-vel/thick prognostic halo update")

      call mpas_timer_start("LTS3-tracer prognostic halo update")
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do
      call mpas_timer_stop("LTS3-tracer prognostic halo update")

      ! --------------------------------------  STEP  ------------------------------------------
      ! --- implicit vert diff for coarse and interface
      !

      if ( (.not. (config_disable_vel_vmix .or. config_disable_vel_all_tend) ) .or. & 
           (.not. (config_disable_tr_vmix .or. config_disable_tr_all_tend) ) ) then

         ! --- update halos for diagnostic variables
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("LTS3-boundary layer depth halo update")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("LTS3-boundary layer depth halo update")
         end if

         call mpas_timer_start("LTS3-diagnostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("LTS3-diagnostic halo update")

         ! --- compute diagnostics for implicit vert diff for velocity, thickness, and tracers
         call mpas_timer_start("LTS3-update diagnostic variables for vert diff")
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
            call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
            end do
            !$omp end do

            if ( haveCoarse == 1 .or. haveInterface == 1 ) then

               ! DIAGNOSTICS UPDATE ---

               call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

            end if

            block => block % next
         end do
         call mpas_timer_stop("LTS3-update diagnostic variables for vert diff")

      end if


      call mpas_timer_start("LTS3-compute vmix for coarse and interface")
      !TODO: this has not been tested yet
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
         call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

         if( haveCoarse == 1 .or. haveInterface == 1 ) then
            call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)
         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-compute vmix for coarse and interface")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END LTS SCHEME
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_stop("LTS3-main loop")
      
      !
      !  A little clean up at the end
      !

      ! reset configs
      if (config_LTS_split) then

         call resetConfigs(domain)

      end if

      call mpas_timer_start("LTS3-cleanup phase")
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1) 
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityZonal', indexSurfaceVelocityZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityMeridional', indexSurfaceVelocityMeridional)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientZonal', indexSSHGradientZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientMeridional', indexSSHGradientMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX) 
         call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH) 
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)
         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
         call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
         call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)

         if (config_prescribe_velocity) then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
            end do
            !$omp end do
         end if

         if (config_prescribe_thickness) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
            end do
            !$omp end do
         end if

         ! ------------------------------------------------------------------
         ! Accumulating various parametrizations of the transport velocity
         ! ------------------------------------------------------------------
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
         end do
         !$omp end do
         ! ------------------------------------------------------------------
         ! End: Accumulating various parametrizations of the transport velocity
         ! ------------------------------------------------------------------

         ! direct application of tidal boundary condition
         if (config_use_tidal_forcing .and. trim(config_tidal_forcing_type) == 'direct') then
           do iCell=1, nCells
             ! artificially assumes boolean mask for now, could generalize to
             ! tappered sponge layer
             if (tidalInputMask(iCell) == 1.0_RKIND) then
               ! compute total depth for relative thickness contribution
               totalDepth = 0.0_RKIND
               do k = 1, maxLevelCell(iCell)
                 totalDepth = totalDepth + restingThickness(k,iCell)
               end do

               ! only modify layer thicknesses on tidal boundary
               do k = 1, maxLevelCell(iCell)
                 layerThicknessNew(k, iCell) = tidalInputMask(iCell)*(tidalBCValue(iCell) + bottomDepth(iCell))*(restingThickness(k,iCell)/totalDepth)
                 !(1.0_RKIND - tidalInputMask(iCell))*layerThicknessNew(k,
                 !iCell)  ! generalized tappered assumption code
               end do
             end if
           end do
         end if


         ! Update the effective density in land ice if we're coupling to land ice
         call ocn_effective_density_in_land_ice_update(meshPool, forcingPool, statePool, scratchPool, err)

         call mpas_reconstruct(meshPool,  normalVelocityNew, &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

         call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do

         call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 2)

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("LTS3-effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("LTS3-effective density halo")
      end if

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         ! DIAGNOSTICS UPDATE ---

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

         block => block % next
      end do

      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
         call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
         call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
         call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)

         call mpas_pool_destroy_pool(tendSum1stPool)
         call mpas_pool_destroy_pool(tendSum2ndPool)
         call mpas_pool_destroy_pool(tendSum3rdPool)
         call mpas_pool_destroy_pool(tendSlowPool)

         call mpas_pool_remove_subpool(block % structs, 'tend_sum_1st')
         call mpas_pool_remove_subpool(block % structs, 'tend_sum_2nd')
         call mpas_pool_remove_subpool(block % structs, 'tend_sum_3rd')
         call mpas_pool_remove_subpool(block % structs, 'tend_slow')

         block => block % next
      end do
      call mpas_timer_stop("LTS3-cleanup phase")


   end subroutine ocn_time_integrator_lts3

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts3_init
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This routine integrates one timestep (dt) using an LTS3 time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integration_lts3_init(domain)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set local time stepping (LTS) regions and lists
    !
    ! Output: LTS instances are written
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    type (domain_type), intent(inout) :: domain

    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: statePool, meshPool, tendPool, LTSPool
    integer, dimension(:), allocatable :: isLTSRegionEdgeAssigned
    integer :: i, iCell, iEdge, iRegion, iHalo, coarseRegions, fineRegions, fineRegionsM1 
    integer, pointer :: nCells, nCellsSolve, nEdges, nEdgesSolve, haveFine, haveCoarse, haveInterface
    integer, dimension(:), pointer :: LTSRegion, LTSRegionLocal, nEdgesOnCell
    integer, dimension(:,:), pointer :: edgesOnCell, nCellsInLTSHalo, nEdgesInLTSHalo
    integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
    integer, dimension(:), pointer :: mpiRankCell, blockPartitionCell
    integer, dimension(2) :: minMaxLTSRegion

    minMaxLTSRegion(1) = 1
    minMaxLTSRegion(2) = 2

    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

       call mpas_pool_get_array(meshPool, 'mpiRankCell', mpiRankCell)
       call mpas_pool_get_array(meshPool, 'blockPartitionCell', blockPartitionCell)

       call mpas_pool_get_array(LTSPool, 'LTSRegion', LTSRegion)
       call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
       call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
       call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
       call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)
       call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

       ! LTS Regions code:
       ! 1 = fine
       ! 2 = coarse
       ! 3 = interface layer 1
       ! 4 = interface layer 2
       ! 5 = fine (to advance when doing 1st stage on interface)
       
       nCellsInLTSHalo(:,:) = 0
       nEdgesInLTSHalo(:,:) = 0

       haveFine = 1
       haveCoarse = 1
       haveInterface = 1

       ! this is to visualize MPI partitions on Paraview
       mpiRankCell(1:nCells) = domain % dminfo % my_proc_id

       ! this is to visualize block partitions on Paraview
       do iCell = 1, nCells
          blockPartitionCell(iCell) = block % blockID
       end do

      block => block % next
    end do

    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
       call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)

       call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
       call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
       call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
       call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)
       call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

       ! this is a loop to build the lists of elements in the fine, coarse, and interface regions
       do iCell = 1, nCellsSolve !this loops up to nCellsSolve because in the time stepping we only want to advance the cells owned by the processor (or block)
          do iRegion = 1,2
             if (iRegion == minMaxLTSRegion(iRegion)) then
                if(LTSRegionLocal(iCell) == minMaxLTSRegion(iRegion)) then
                   nCellsInLTSHalo(iRegion,1) = nCellsInLTSHalo(iRegion,1) + 1
                   cellsInLTSHalo(iRegion,1,nCellsInLTSHalo(iRegion,1)) = iCell
                end if
                if(LTSRegionLocal(iCell) == (minMaxLTSRegion(iRegion) + 2) ) then
                   nCellsInLTSHalo(iRegion,2) = nCellsInLTSHalo(iRegion,2) + 1
                   cellsInLTSHalo(iRegion,2,nCellsInLTSHalo(iRegion,2)) = iCell
                end if
             end if
          end do
          if (LTSRegionLocal(iCell) == 5) then
             nCellsInLTSHalo(1,3) = nCellsInLTSHalo(1,3) + 1
             cellsInLTSHalo(1,3,nCellsInLTSHalo(1,3)) = iCell
          end if
       end do

       if ((nCellsInLTSHalo(1,1) == 0) .and. (nCellsInLTSHalo(1,3) == 0)) then
             haveFine = 0
       end if

       if (nCellsInLTSHalo(2,1) == 0) then
             haveCoarse = 0
       end if

       if ((nCellsInLTSHalo(1,2) == 0) .and. (nCellsInLTSHalo(2,2) == 0)) then
          haveInterface = 0
       end if

       ! below we fill out the lists for the edges, according to the
       ! LTSRegionLocal that have been assigned to the cells.
       ! we move from the fine to the coarse (i.e. from the fine to the nearest
       ! LTS region in the direction of the coarse).
       ! Note that edges shared between cells of different LTS regions are owned
       ! by the cell in the LTS region closest to the fine region,
       ! see Figure 3 in "Conservative explicit local time-stepping schemes for
       ! the shallow water equations" by Hoang et al. (halo edges however
       ! are owned by whatever processor they are initially assigned to)

       allocate(isLTSRegionEdgeAssigned(nEdgesSolve))
       isLTSRegionEdgeAssigned(:) = 0

       do iCell = 1, nCellsInLTSHalo(1,1)
          do i = 1, nEdgesOnCell(cellsInLTSHalo(1,1,iCell))
             iEdge = edgesOnCell(i,cellsInLTSHalo(1,1,iCell))
             if (iEdge .le. nEdgesSolve) then
                if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                   nEdgesInLTSHalo(1,1) = nEdgesInLTSHalo(1,1) + 1
                   edgesInLTSHalo(1,1, nEdgesInLTSHalo(1,1)) = iEdge
                   isLTSRegionEdgeAssigned(iEdge) = 1
                end if
             end if
          end do
       end do

       fineRegions = 3
       fineRegionsM1 = 2
       do iHalo = 1, fineRegionsM1
          do iCell = 1, nCellsInLTSHalo(1, fineRegions - iHalo + 1)
             do i = 1, nEdgesOnCell(cellsInLTSHalo(1, fineRegions - iHalo + 1, iCell))
                iEdge = edgesOnCell(i,cellsInLTSHalo(1, fineRegions - iHalo + 1, iCell))
                if (iEdge .le. nEdgesSolve) then
                   if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                      nEdgesInLTSHalo(1, fineRegions - iHalo + 1) = nEdgesInLTSHalo(1, fineRegions - iHalo + 1) + 1
                      edgesInLTSHalo(1, fineRegions - iHalo + 1, nEdgesInLTSHalo(1, fineRegions - iHalo + 1)) = iEdge
                      isLTSRegionEdgeAssigned(iEdge) = 1 
                   end if
                end if
             end do
          end do
       end do

       coarseRegions = 2
       do iHalo = 1, coarseRegions
          do iCell = 1, nCellsInLTSHalo(2, coarseRegions - iHalo + 1)
             do i = 1, nEdgesOnCell(cellsInLTSHalo(2,coarseRegions - iHalo + 1,iCell))
                iEdge = edgesOnCell(i,cellsInLTSHalo(2,coarseRegions - iHalo + 1,iCell))
                if (iEdge .le. nEdgesSolve) then
                   if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                      nEdgesInLTSHalo(2, coarseRegions - iHalo + 1) = nEdgesInLTSHalo(2, coarseRegions - iHalo + 1) + 1
                      edgesInLTSHalo(2, coarseRegions - iHalo + 1, nEdgesInLTSHalo(2, coarseRegions - iHalo + 1)) = iEdge
                      isLTSRegionEdgeAssigned(iEdge) = 1
                   end if
                end if
             end do
          end do
       end do

       deallocate(isLTSRegionEdgeAssigned)

       block => block % next
    end do

    ! BEGIN: CHECK TO SEE IF THE INTERFACE BLOCKS HAVE ENOUGH CELLS
    block => domain % blocklist
    do while (associated(block))

       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

       call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
       call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
       call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

       call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

       if(haveInterface == 1) then
          if (nCellsInLTSHalo(1,2) + nCellsInLTSHalo(2,2) < 100) then
             print*, 'WARNING: interface block', block % blockID, 'has', nCellsInLTSHalo(1,2) + nCellsInLTSHalo(2,2), 'cells'
          end if
       end if

       block => block % next
    end do

    ! END: CHECK TO SEE IF THE INTERFACE BLOCKS HAVE ENOUGH CELLS

    end subroutine ocn_time_integration_lts3_init

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ocn_lts_fast_tends(statePool, diagnosticsPool, LTSPool, meshPool, tendPool, timeLevelIn, computeOnFineBig, computeOnFineSmall, computeOnCoarse, computeOnInterface)

   integer, intent(in) :: timeLevelIn, computeOnFineBig, computeOnFineSmall, computeOnCoarse, computeOnInterface

   type (mpas_pool_type), intent(in) :: &
      statePool          !< Input: state variables

   type (mpas_pool_type), intent(in) :: &
      diagnosticsPool          !< Input: diagnostics variables

   type (mpas_pool_type), intent(in) :: &
      LTSPool          !< Input: LTS data

   type (mpas_pool_type), intent(in) :: &
      meshPool          !< Input: mesh information

   type (mpas_pool_type), intent(inout) :: &
      tendPool          !< Input: tendency variables

   integer, pointer :: nEdges, nVertLevels
   integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask, edgesOnCell, edgeSignOnCell
   real (kind=RKIND), dimension(:), pointer :: dcEdge, ssh, dvEdge, areaCell, bottomDepth
   real (kind=RKIND), pointer :: config_self_attraction_and_loading_beta
   integer, dimension(:,:), pointer :: nCellsInLTSHalo, nEdgesInLTSHalo
   integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
   real (kind=RKIND), dimension(:,:), pointer :: normalVelocityTend, layerThicknessTend, layerThicknessEdge, normalVelocity, layerThickness, zMid, zTop
   integer, dimension(:), pointer :: maxLevelEdgeBot, maxLevelEdgeTop, nEdgesOnCell, maxLevelCell, nCellsArray

   integer :: iEdge, cell1, cell2, k, ie, iRegion, nRegions, ic, iCell, i, nCells
   real (kind=RKIND) :: invdcEdge, betaSelfAttrLoad, flux, invAreaCell

   call mpas_pool_get_config(ocnConfigs, 'config_self_attraction_and_loading_beta', config_self_attraction_and_loading_beta) 

   betaSelfAttrLoad = config_self_attraction_and_loading_beta
   nRegions = 2

   call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
   call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
   call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
   call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
   call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
   call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
   call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
   call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
   call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
   call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
   call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
   call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
   call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
   call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
   call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

   call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
   call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
   call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
   call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)

   call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
   call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

   call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
   call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
   call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)

   call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)
   call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevelIn)
   call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevelIn)

   ! inline computation of the diagnostics
   ! TODO: check if this can be done using LTS regions

   ! layerThicknessEdge
   do iEdge = 1, nEdges
      layerThicknessEdge(:, iEdge) = -1.0e34_RKIND
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)
      do k = 1, maxLevelEdgeTop(iEdge)
         layerThicknessEdge(k,iEdge) = 0.5_RKIND * (layerThickness(k,cell1) + layerThickness(k,cell2))
      end do
   end do

   ! ssh
   nCells = nCellsArray( size(nCellsArray) )
   do iCell = 1, nCells
      k = maxLevelCell(iCell)
      zMid(k:nVertLevels,iCell) = -bottomDepth(iCell) + 0.5_RKIND*layerThickness(k,iCell)
      zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) + layerThickness(k,iCell)

      do k = maxLevelCell(iCell)-1, 1, -1
         zMid(k,iCell) = zMid(k+1,iCell) + 0.5_RKIND*(  layerThickness(k+1,iCell) + layerThickness(k  ,iCell))
         zTop(k,iCell) = zTop(k+1,iCell) + layerThickness(k  ,iCell)
      end do

      ssh(iCell) = zTop(1,iCell)

   end do

   ! inline computation of the tendencies

   normalVelocityTend(:,:) = 0.0_RKIND
   layerThicknessTend(:, :) = 0.0_RKIND

   ! interface
   if (computeOnInterface == 1) then
      do iRegion = 1, nRegions
         ! velocity tendency
         do ie = 1, nEdgesInLTSHalo(iRegion,2)
            iEdge = edgesInLTSHalo(iRegion,2,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            do k=1,maxLevelEdgeTop(iEdge)
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  - betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
         ! thickness tendency
         do ic = 1, nCellsInLTSHalo(iRegion,2)
            iCell = cellsInLTSHalo(iRegion,2,ic)
            invAreaCell = 1.0_RKIND / areaCell(iCell)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = 1, maxLevelEdgeBot(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
               end do
            end do
         end do
      end do
   end if

   if (computeOnFineBig == 1) then
      ! velocity tendency
      do ie = 1, nEdgesInLTSHalo(1,1)
         iEdge = edgesInLTSHalo(1,1,ie)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) & 
                                          - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  - betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
         end do
      end do
      ! thickness tendency
      do ic = 1, nCellsInLTSHalo(1,1)
         iCell = cellsInLTSHalo(1,1,ic)
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeBot(iEdge)
               flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
               layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
            end do
         end do
      end do
   end if

   if (computeOnFineSmall == 1) then
      ! velocity tendency
      do ie = 1, nEdgesInLTSHalo(1,3)
         iEdge = edgesInLTSHalo(1,3,ie)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) & 
                                          - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  - betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
         end do
      end do
      ! thickness tendency
      do ic = 1, nCellsInLTSHalo(1,3)
         iCell = cellsInLTSHalo(1,3,ic)
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeBot(iEdge)
               flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
               layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
            end do
         end do
      end do
   end if

   if (computeOnCoarse == 1) then
      ! velocity tendency
      do ie = 1, nEdgesInLTSHalo(2,1)
         iEdge = edgesInLTSHalo(2,1,ie)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         do k=1,maxLevelEdgeTop(iEdge)
            normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                          - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  - betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
         end do
      end do
      ! thickness tendency
      do ic = 1, nCellsInLTSHalo(2,1)
         iCell = cellsInLTSHalo(2,1,ic)
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelEdgeBot(iEdge)
               flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
               layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
            end do
         end do
      end do
   end if

   end subroutine ocn_lts_fast_tends

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on
   ! a sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance

end module ocn_time_integration_lts3
