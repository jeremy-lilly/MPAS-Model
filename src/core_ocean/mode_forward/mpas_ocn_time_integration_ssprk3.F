! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_ssprk3
!
!> \brief MPAS ocean SSPRK3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This module contains the SSPRK3 time integration
!
!-----------------------------------------------------------------------

module ocn_time_integration_ssprk3

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

   public :: ocn_time_integrator_ssprk3, ocn_time_integrator_compute_vel_tends
   public :: ocn_time_integrator_compute_thick_tends, ocn_time_integrator_compute_tracer_tends

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_ssprk3
!
!> \brief MPAS ocean SSPRK3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This routine integrates one timestep (dt) using an SSPRK3 time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_ssprk3(domain,dt)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    !   Strong Stability Preserving Runge Kutta of order 3 scheme
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    type (domain_type), intent(inout) :: domain
    real (kind=RKIND), intent(in) :: dt

    integer :: iCell, iEdge, k, err

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

    integer :: ssprk_step

    real (kind=RKIND), dimension(3) :: weights_old, weights_new, weights_tend, forcingTimeIncrementSSPRK3

    ! Dimensions
    integer, pointer :: nCells, nEdges, nVertLevels, num_tracers

    ! Config options
    logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
    logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
    logical, pointer :: config_use_GM
    logical, pointer :: config_use_cvmix_kpp
    logical, pointer :: config_use_tracerGroup
    logical, pointer :: config_disable_vel_vmix
    logical, pointer :: config_disable_thick_all_tend
    logical, pointer :: config_disable_vel_all_tend
    logical, pointer :: config_disable_tr_all_tend
    real (kind=RKIND), pointer :: config_mom_del4
    real (kind=RKIND), pointer :: config_drying_min_cell_height
    logical, pointer :: config_use_wetting_drying
    logical, pointer :: config_verify_not_dry
    character (len=StrKIND), pointer :: config_land_ice_flux_mode

    ! State indices
    integer, pointer :: indexTemperature
    integer, pointer :: indexSalinity

    ! Diagnostics Indices
    integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
    integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

    ! Mesh array pointers
    integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop
    real (kind=RKIND), dimension(:), pointer :: bottomDepth

    ! Tend Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, layerThicknessTend
    real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend

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
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew
    real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew
    real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew, wettingVelocity
    real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew
    real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupCur, tracersGroupNew

    ! Diagnostics Field Pointers
    type (field1DReal), pointer :: boundaryLayerDepthField, effectiveDensityField
    type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

    ! Tracer Group Iteration
    type (mpas_pool_iterator_type) :: groupItr
    character (len=StrKIND) :: modifiedGroupName
    character (len=StrKIND) :: configName

    ! Tidal boundary condition
    logical, pointer :: config_use_tidal_forcing
    character (len=StrKIND), pointer :: config_tidal_forcing_type
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue
    real (kind=RKIND), dimension(:,:), pointer :: restingThickness
    real (kind=RKIND) :: totalDepth
    
    real (kind=RKIND), pointer :: forcingTimeIncrement

    ! Get config options
    call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
    call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
    call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
    call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
    call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
    call mpas_pool_get_config(domain % configs, 'config_use_GM', config_use_GM)
    call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
    call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
    call mpas_pool_get_config(domain % configs, 'config_disable_vel_vmix', config_disable_vel_vmix)
    call mpas_pool_get_config(domain % configs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_use_wetting_drying', config_use_wetting_drying)
    call mpas_pool_get_config(domain % configs, 'config_verify_not_dry', config_verify_not_dry)
    call mpas_pool_get_config(domain % configs, 'config_drying_min_cell_height', config_drying_min_cell_height)
    call mpas_pool_get_config(domain % configs, 'config_use_tidal_forcing', config_use_tidal_forcing)
    call mpas_pool_get_config(domain % configs, 'config_tidal_forcing_type', config_tidal_forcing_type)  

    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'state', statePool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
       call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
       call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

       call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
       call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
       call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
       call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

       call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)

       forcingTimeIncrement = 0.0_RKIND

       !$omp do schedule(runtime) private(k)
       do iEdge = 1, nEdges
          do k = 1, maxLevelEdgeTop(iEdge)
             normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
          end do
       end do
       !$omp end do

       !$omp do schedule(runtime) private(k)
       do iCell = 1, nCells
          sshNew(iCell) = sshCur(iCell)
          do k = 1, maxLevelCell(iCell)
             layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
          end do
       end do
       !$omp end do

       call mpas_pool_begin_iteration(tracersPool)
       do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

          if ( groupItr % memberType == MPAS_POOL_FIELD ) then

             call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupCur, 1)
             call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracersGroupNew, 2)

             if ( associated(tracersGroupCur) .and. associated(tracersGroupNew) ) then
                !$omp do schedule(runtime) private(k)
                do iCell = 1, nCells  
                   do k = 1, maxLevelCell(iCell)
                      tracersGroupNew(:, k, iCell) = tracersGroupCur(:, k, iCell) 
                   end do
                end do
                !$omp end do
             end if
          end if
       end do

       if (associated(highFreqThicknessCur)) then
          !$omp do schedule(runtime)
          do iCell = 1, nCells
             highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
           end do
           !$omp end do
       end if

       if (associated(lowFreqDivergenceCur)) then
          !$omp do schedule(runtime)
          do iCell = 1, nCells
             lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
          end do
          !$omp end do
       end if

       block => block % next
    end do

    call mpas_threading_barrier() 

    weights_old(1) = 1.0
    weights_old(2) = 0.75 !use this for SSPRK3
    !weights_old(2) = 0.5   !use this for SSPRK2
    weights_old(3) = 1./3.

    weights_new(1) = 0.0
    weights_new(2) = 0.25 !use this for SSPRK3
    !weights_new(2) = 0.5   !use this for SSPRK2
    weights_new(3) = 2./3.

    weights_tend(1) = dt * 1.0
    weights_tend(2) = dt * 0.25 !use this for SSPRK3
    !weights_tend(2) = dt * 0.5   !use this for SSPRK2
    weights_tend(3) = dt * 2./3.

    forcingTimeIncrementSSPRK3(1) = 0 * dt 
    forcingTimeIncrementSSPRK3(2) = dt 
    forcingTimeIncrementSSPRK3(3) = dt/2 

    err = 0 ! for MPI_Barrier tests

    call mpas_timer_start("SSPRK3-main loop")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BEGIN SSPRK loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ssprk_step = 1, 3

      if (config_disable_thick_all_tend .and. config_disable_vel_all_tend .and. config_disable_tr_all_tend) then
          exit ! don't compute in loop meant to update velocity, thickness, and tracers
      end if


      ! Update halos for diagnostic variables.
      if (config_use_cvmix_kpp) then
         call mpas_timer_start("SSPRK3-boundary layer depth halo update")
         call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
         call mpas_timer_stop("SSPRK3-boundary layer depth halo update")
      end if


      call mpas_timer_start("SSPRK3-diagnostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
      if (config_mom_del4 > 0.0_RKIND) then
         call mpas_dmpar_field_halo_exch(domain, 'divergence')
         call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
      end if
      call mpas_timer_stop("SSPRK3-diagnostic halo update")
      call mpas_threading_barrier()


      ! Compute tendencies for high frequency thickness
      if (config_use_freq_filtered_thickness) then
         call mpas_timer_start("SSPRK3-tendency computations")
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)

            call ocn_tend_freq_filtered_thickness(tendPool, statePool, diagnosticsPool, meshPool, 2)
            call mpas_threading_barrier()
            block => block % next
         end do
         call mpas_timer_stop("SSPRK3-tendency computations")

         call mpas_timer_start("SSPRK3-prognostic halo update")
         call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
         call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')
         call mpas_timer_stop("SSPRK3-prognostic halo update")
         call mpas_threading_barrier()

         ! Compute next substep state for high frequency thickness.
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
            call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

            !$omp do schedule(runtime)
            do iCell = 1, nCells
               highFreqThicknessNew(:, iCell) = weights_old(ssprk_step) * highFreqThicknessCur(:,iCell) &
                                                + weights_new(ssprk_step) * highFreqThicknessNew(:,iCell) &
                                                + weights_tend(ssprk_step) * highFreqThicknessTend(:, iCell)
            end do
            !$omp end do
            block => block % next
         end do

      end if


      ! require that cells don't dry out
      if (config_use_wetting_drying) then
         !TODO: MISSING FOR NOW
      end if 

      block => domain % blocklist
      do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
           call mpas_pool_get_array(forcingPool, 'forcingTimeIncrement', forcingTimeIncrement)
           forcingTimeIncrement = forcingTimeIncrementSSPRK3(ssprk_step)
           block => block % next
      end do

      ! Compute tendencies for velocity, thickness, and tracers.
      call mpas_timer_start("SSPRK3 vel/thick/tracers tendency computations")
      !print*, "ssprk_step---------->", ssprk_step
      call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)
      block => domain % blocklist
      do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call ocn_time_integrator_compute_vel_tends(block, dt, err, statePool, 2)
           call ocn_time_integrator_compute_thick_tends(block, dt, err, statePool, 2)
           call ocn_time_integrator_compute_tracer_tends( block, dt, err, statePool, 2)
         block => block % next
      end do
      call mpas_timer_stop("SSPRK3 vel/thick/tracers tendency computations")


      ! Update halos for prognostic variables.
      call mpas_timer_start("SSPRK3 vel/thick/tracers prognostic halo update")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'tend',tendPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

      call mpas_dmpar_field_halo_exch(domain, 'tendNormalVelocity')
      call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

      call mpas_pool_begin_iteration(tracersTendPool)
      do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, trim(groupItr % memberName))
         end if
      end do
      call mpas_timer_stop("SSPRK3 vel/thick/tracers prognostic halo update")
      call mpas_threading_barrier() 


      call mpas_timer_start("SSPRK3 advance soln")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend) 
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

         call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
         call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)
         call mpas_pool_get_array(diagnosticsPool, 'wettingVelocity', wettingVelocity)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalVelocityNew(:,iEdge) = weights_old(ssprk_step) * normalVelocityCur(:,iEdge) & 
                                         + weights_new(ssprk_step) * normalVelocityNew(:,iEdge) & 
                                         + weights_tend(ssprk_step) * normalVelocityTend(:,iEdge)
            normalVelocityNew(:, iEdge) = normalVelocityNew(:, iEdge) * (1.0_RKIND - wettingVelocity(:, iEdge))
         end do
         !$omp end do

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               configName = 'config_use_' // trim(groupItr % memberName)
               call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

               if ( config_use_tracerGroup ) then 
                  call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

                  if ( associated(tracersGroupNew) ) then
                     !$omp do schedule(runtime) private(k)
                     do iCell = 1, nCells
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupNew(:, k, iCell) = tracersGroupNew(:, k,iCell) * layerThicknessNew(k,iCell)
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end if
         end do


         !$omp do schedule(runtime)
         do iCell = 1, nCells
            do k = 1, maxLevelCell(iCell)
               layerThicknessNew(k,iCell) = weights_old(ssprk_step) * layerThicknessCur(k,iCell) &
                                         + weights_new(ssprk_step) * layerThicknessNew(k,iCell) &
                                         + weights_tend(ssprk_step) * layerThicknessTend(k,iCell)
            end do
         end do
         !$omp end do

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               configName = 'config_use_' // trim(groupItr % memberName)
               call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

               if ( config_use_tracerGroup ) then
                  call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)
                  call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

                  modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                  call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                  if ( associated(tracersGroupNew) .and. associated(tracersGroupTend) ) then
                     !$omp do schedule(runtime) private(k)
                     do iCell = 1, nCells
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupNew(:,k,iCell) = ( weights_old(ssprk_step) * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell)  &
                                                          + weights_new(ssprk_step) * tracersGroupNew(:,k,iCell) &
                                                          + weights_tend(ssprk_step) * tracersGroupTend(:,k,iCell) ) / layerThicknessNew(k,iCell)
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end if
         end do

         !if (associated(highFreqThicknessNew)) then    !TODO: WHY IS THIS ADVANCED HERE BUT ALSO ABOVE FOR THE SAME TENDENCY?
         !   !$omp do schedule(runtime)
         !   do iCell = 1, nCells
         !      highFreqThicknessNew(:, iCell) = weights_old(ssprk_step) * highFreqThicknessCur(:,iCell) &
         !                                       + weights_new(ssprk_step) * highFreqThicknessNew(:,iCell) &
         !                                       + weights_tend(ssprk_step) * highFreqThicknessTend(:, iCell)
         !   end do
         !   !$omp end do
         !end if

         if (associated(lowFreqDivergenceNew)) then ! TODO: WHERE IS lowFreqDivergenceTend COMPUTED?
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               lowFreqDivergenceNew(:, iCell) = weights_old(ssprk_step) * lowFreqDivergenceCur(:,iCell) &
                                                + weights_new(ssprk_step) * lowFreqDivergenceNew(:,iCell) &
                                                + weights_tend(ssprk_step) * lowFreqDivergenceTend(:, iCell)
            end do
            !$omp end do
         end if

         block => block % next
      end do
      call mpas_timer_stop("SSPRK3 advance soln")
      call mpas_threading_barrier()


      ! Compute diagnostics for velocity, thickness, and tracers.
      call mpas_timer_start("SSPRK3-update diagnostic variables")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)
         call mpas_threading_barrier()


         ! ------------------------------------------------------------------
         ! Accumulating various parametrizations of the transport velocity
         ! ------------------------------------------------------------------
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
         end do
         !$omp end do

         ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if
         ! respective flags are turned on
         if (config_use_GM) then
            call ocn_gm_compute_Bolus_velocity(statePool, diagnosticsPool, meshPool, scratchPool, timeLevelIn=2)
         end if
         call mpas_threading_barrier()

         if (config_use_GM) then
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalTransportVelocity(:,iEdge) + normalGMBolusVelocity(:,iEdge)
         end do
         !$omp end do
         end if
         ! ------------------------------------------------------------------
         ! End: Accumulating various parametrizations of the transport velocity
         ! ------------------------------------------------------------------

         block => block % next
      end do
      call mpas_timer_stop("SSPRK3-update diagnostic variables")
      call mpas_threading_barrier()

   end do
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! END SSPRK loop
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! verify that cells are not dry at conclusion of time step
      if (config_use_wetting_drying) then
         !TODO: MISSING FOR NOW
      end if

      call mpas_timer_stop("SSPRK3-main loop")
      call mpas_threading_barrier()

      !
      !  A little clean up at the end
      !

      call mpas_timer_start("SSPRK3-compute vmix")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)


         call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)
         call mpas_threading_barrier()

         block => block % next
      end do
      call mpas_timer_stop("SSPRK3-compute vmix")
      call mpas_threading_barrier()

      call mpas_timer_start("SSPRK3-cleanup phase") 

      call mpas_timer_start("SSPRK3-implicit vert mix")
      ! Update halo on u and tracers, which were just updated for implicit
      ! vertical mixing.  If not done, this leads to lack of volume conservation.  It is required because halo
      ! updates in SSPRK3 are only  conducted on tendencies, not on the velocity and tracer fields.  So this
      ! update is required to communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("SSPRK3-implicit vert mix halos")

      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do

      call mpas_timer_stop("SSPRK3-implicit vert mix halos")
      call mpas_timer_stop("SSPRK3-implicit vert mix")
      call mpas_threading_barrier()


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
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
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

         ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if
         ! respective flags are turned on
         if (config_use_GM) then
            call ocn_gm_compute_Bolus_velocity(statePool, diagnosticsPool, meshPool, scratchPool, timeLevelIn=2)
         end if
         call mpas_threading_barrier()

         if (config_use_GM) then
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalTransportVelocity(:,iEdge) + normalGMBolusVelocity(:,iEdge)
         end do
         !$omp end do
         end if
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

         if (.not. config_disable_vel_vmix) then
            call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)
         end if
         call mpas_threading_barrier()

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
         call mpas_threading_barrier()

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do

         call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 2)

         if (config_use_GM) then
            call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         end if
         call mpas_threading_barrier()

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("SSPRK3-effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("SSPRK3-effective density halo")
      end if

      call mpas_timer_stop("SSPRK3-cleanup phase")

      call mpas_threading_barrier()

    end subroutine ocn_time_integrator_ssprk3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


    subroutine ocn_time_integrator_compute_vel_tends(block, dt, err, statePool, timeLevelIn)   

      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err
      type (mpas_pool_type), intent(in) :: statePool 
      integer, intent(in), optional :: timeLevelIn            

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool

      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, normalVelocity
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThickness

      logical, pointer :: config_filter_btr_mode

      integer :: timeLevel

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, timeLevel)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity. 
      if (associated(highFreqThickness)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThickness,layerThicknessEdge, normalVelocity, &
            ssh, dt, vertAleTransportTop, err, highFreqThickness)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThickness,layerThicknessEdge, normalVelocity, &
            ssh, dt, vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, timeLevel, dt)
      call mpas_threading_barrier()

    end subroutine ocn_time_integrator_compute_vel_tends


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


    subroutine ocn_time_integrator_compute_thick_tends(block, dt, err, statePool, timeLevelIn) 

      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err
      type (mpas_pool_type), intent(in) :: statePool
      integer, intent(in), optional :: timeLevelIn

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool

      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity
      real (kind=RKIND), dimension(:, :), pointer :: highFreqThickness

      logical, pointer :: config_filter_btr_mode

      integer :: timeLevel

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity) 

      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, timeLevel)


      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity. 
      if (associated(highFreqThickness)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThickness, layerThicknessEdge, normalTransportVelocity, &
            ssh, dt, vertAleTransportTop, err, highFreqThickness)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThickness, layerThicknessEdge, normalTransportVelocity, &
            ssh, dt, vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_compute_thick_tends


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


   subroutine ocn_time_integrator_compute_tracer_tends(block, dt, err, statePool, timeLevelIn)

      type (block_type), intent(in) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err
      type (mpas_pool_type), intent(in) :: statePool
      integer, intent(in), optional :: timeLevelIn

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool
      type (mpas_pool_type), pointer :: scratchPool, tendPool
      type (mpas_pool_type), pointer :: swForcingPool

      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, normalVelocity
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessEdge, vertAleTransportTop
      real (kind=RKIND), dimension(:, :), pointer :: normalTransportVelocity, highFreqThickness

      logical, pointer :: config_filter_btr_mode

      integer :: timeLevel

      logical :: activeTracersOnlyIn

      activeTracersOnlyIn = .false.

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, timeLevel)

      ! advection of u uses u, while advection of layerThickness and tracers use normalTransportVelocity.
      if (associated(highFreqThickness)) then
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThickness, layerThicknessEdge, normalTransportVelocity, &
            ssh, dt, vertAleTransportTop, err, highFreqThickness)
      else
         call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
            layerThickness, layerThicknessEdge, normalTransportVelocity, &
            ssh, dt, vertAleTransportTop, err)
      endif
      call mpas_threading_barrier()

      if (config_filter_btr_mode) then
          call ocn_filter_btr_mode_tend_vel(tendPool, statePool, diagnosticsPool, meshPool, timeLevel)
      endif
      call mpas_threading_barrier()

      call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dt, activeTracersOnlyIn, timeLevel)
      call mpas_threading_barrier()

   end subroutine ocn_time_integrator_compute_tracer_tends


end module ocn_time_integration_ssprk3
