! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module sw_time_integration

   use mpas_vector_reconstruction
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_log
   use mpas_timer

   use sw_constants

   contains


   subroutine sw_timestep(domain, dt, timeStamp, totalTimeAdvance, totalTimeTendencies, totalTimeBarrier, totalTimeInterface)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp
      real (kind=RKIND), intent(inout) :: totalTimeAdvance, totalTimeTendencies, totalTimeBarrier, totalTimeInterface

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool

      character (len=StrKIND), pointer :: xtime
      character (len=StrKIND), pointer :: config_time_integration

      call mpas_pool_get_config(domain % configs, 'config_time_integration', config_time_integration)

      if (trim(config_time_integration) == 'RK4') then
         call sw_rk4(domain, dt) 
      else if (trim(config_time_integration) == 'SSPRK2') then
         call sw_ssprk(domain, dt, 2)
      else if (trim(config_time_integration) == 'SSPRK3') then
         call sw_ssprk(domain, dt, 3)
      else if (trim(config_time_integration) == 'FWDEULR') then
         call sw_ssprk(domain, dt, 1)
      else if (trim(config_time_integration) == 'LTS2') then
         call sw_lts(domain, dt, 2, totalTimeAdvance, totalTimeTendencies, totalTimeBarrier, totalTimeInterface)
      else if (trim(config_time_integration) == 'LTS3') then
         call sw_lts(domain, dt, 3, totalTimeAdvance, totalTimeTendencies, totalTimeBarrier, totalTimeInterface)
      else
         call mpas_log_write('Unknown time integration option '//trim(config_time_integration), MPAS_LOG_ERR)
         call mpas_log_write('Currently, only ''RK4'', ''FWDEULR'', ''SSPRK2'', ''SSPRK3'', ''LTS2'' and ''LTS3'' are supported.', MPAS_LOG_CRIT)
      end if

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_array(statePool, 'xtime', xtime, 2)
         xtime = timeStamp 
         block => block % next
      end do

   end subroutine sw_timestep


   subroutine sw_rk4(domain, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   4th order Runge-Kutta
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      integer :: iCell, k
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: provisStatePool
      type (mpas_pool_type), pointer :: prevProvisPool, nextProvisPool

      integer :: rk_step

      real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels

      real (kind=RKIND), dimension(:,:), pointer :: uOld, uNew, uProvis, uTend
      real (kind=RKIND), dimension(:,:), pointer :: hOld, hNew, hProvis, hTend
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersOld, tracersNew, tracersProvis, tracersTend

      integer, pointer :: config_test_case
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4


      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)
      call mpas_pool_get_config(domain % configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)

     !
     ! Initialize time_levs(2) with state at current time
     ! Initialize first RK state
     ! Couple tracers time_levs(2) with h in time-levels
     ! Initialize RK weights
     !
     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)

        allocate(provisStatePool)
        call mpas_pool_create_pool(provisStatePool)
        call mpas_pool_clone_pool(statePool, provisStatePool, 1)

        call mpas_pool_add_subpool(block % structs, 'provis_state', provisStatePool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(statePool, 'u', uOld, 1)
        call mpas_pool_get_array(statePool, 'u', uNew, 2)
        call mpas_pool_get_array(statePool, 'h', hOld, 1)
        call mpas_pool_get_array(statePool, 'h', hNew, 2)
        call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)
        call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

        uNew(:,:) = uOld(:,:)
        hNew(:,:) = hOld(:,:)

        call mpas_pool_initialize_time_levels(statePool)

        do iCell = 1, nCells  ! couple tracers to h
          do k = 1, nVertLevels
            tracersNew(:,k,iCell) = tracersOld(:,k,iCell) * hOld(k,iCell)
           end do
        end do

        block => block % next
     end do

     block => domain % blocklist
     do while(associated(block))
        if (associated(block % prev)) then
           call mpas_pool_get_subpool(block % prev % structs, 'provis_state', prevProvisPool)
        else
           nullify(prevProvisPool)
        end if

        if (associated(block % next)) then
           call mpas_pool_get_subpool(block % next % structs, 'provis_state', nextProvisPool)
        else
           nullify(nextProvisPool)
        end if

        call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

        if (associated(prevProvisPool) .and. associated(nextProvisPool)) then
           call mpas_pool_link_pools(provisStatePool, prevProvisPool, nextProvisPool)
        else if (associated(prevProvisPool)) then
           call mpas_pool_link_pools(provisStatePool, prevProvisPool)
        else if (associated(nextProvisPool)) then
           call mpas_pool_link_pools(provisStatePool, nextPool=nextProvisPool)
        else
           call mpas_pool_link_pools(provisStatePool)
        end if

        call mpas_pool_link_parinfo(block, provisStatePool)

        block => block % next
     end do


     rk_weights(1) = dt/6.
     rk_weights(2) = dt/3.
     rk_weights(3) = dt/3.
     rk_weights(4) = dt/6.

     rk_substep_weights(1) = dt/2.
     rk_substep_weights(2) = dt/2.
     rk_substep_weights(3) = dt
     rk_substep_weights(4) = 0.


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     ! BEGIN RK loop 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     do rk_step = 1, 4

! --- update halos for diagnostic variables
        !call mpas_dmpar_field_halo_exch(domain, 'pv_edge', timeLevel=1)

        if (config_h_mom_eddy_visc4 > 0.0) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
            call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
        end if

! --- compute tendencies

       block => domain % blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
          call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

          call sw_compute_tend(tendPool, provisStatePool, meshPool, 1)
          call sw_compute_scalar_tend(tendPool, provisStatePool, meshPool, 1)
          call sw_enforce_boundary_edge(tendPool, meshPool)
          block => block % next
       end do

! --- update halos for prognostic variables

       call mpas_dmpar_field_halo_exch(domain, 'tend_u')
       call mpas_dmpar_field_halo_exch(domain, 'tend_h')
       call mpas_dmpar_field_halo_exch(domain, 'tend_tracers')

! --- compute next substep state

       if (rk_step < 4) then
          block => domain % blocklist
          do while (associated(block))
             call mpas_pool_get_subpool(block % structs, 'state', statePool)
             call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
             call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

             call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
             call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

             call mpas_pool_get_array(provisStatePool, 'u', uProvis)
             call mpas_pool_get_array(provisStatePool, 'h', hProvis)
             call mpas_pool_get_array(provisStatePool, 'tracers', tracersProvis)

             call mpas_pool_get_array(statePool, 'u', uOld, 1)
             call mpas_pool_get_array(statePool, 'h', hOld, 1)
             call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

             call mpas_pool_get_array(tendPool, 'u', uTend)
             call mpas_pool_get_array(tendPool, 'h', hTend)
             call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

             uProvis(:,:) = uOld(:,:) + rk_substep_weights(rk_step) * uTend(:,:)
             hProvis(:,:) = hOld(:,:) + rk_substep_weights(rk_step) * hTend(:,:)
             do iCell = 1, nCells
                do k = 1, nVertLevels
                   tracersProvis(:,k,iCell) = ( hOld(k,iCell) * tracersOld(:,k,iCell) +  &
                                                rk_substep_weights(rk_step) * tracersTend(:,k,iCell) ) / hProvis(k,iCell)
                end do
             end do
             if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uProvis(:,:) = uOld(:,:)
             end if
             call sw_compute_solve_diagnostics(dt, provisStatePool, meshPool)
             block => block % next
          end do
       end if

!--- accumulate update (for RK4)

       block => domain % blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(statePool, 'u', uOld, 1)
          call mpas_pool_get_array(statePool, 'h', hOld, 1)
          call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

          call mpas_pool_get_array(statePool, 'u', uNew, 2)
          call mpas_pool_get_array(statePool, 'h', hNew, 2)
          call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

          call mpas_pool_get_array(tendPool, 'u', uTend)
          call mpas_pool_get_array(tendPool, 'h', hTend)
          call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

          uNew(:,:) = uNew(:,:) + rk_weights(rk_step) * uTend(:,:) 
          hNew(:,:) = hNew(:,:) + rk_weights(rk_step) * hTend(:,:) 
          do iCell = 1, nCells
             do k = 1, nVertLevels
                tracersNew(:,k,iCell) = tracersNew(:,k,iCell) + rk_weights(rk_step) * tracersTend(:,k,iCell)
             end do
          end do
          block => block % next
       end do

      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! END RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      !
      !  A little clean up at the end: decouple new scalar fields and compute diagnostics for new state
      !
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(statePool, 'u', uOld, 1)
         call mpas_pool_get_array(statePool, 'u', uNew, 2)
         call mpas_pool_get_array(statePool, 'h', hNew, 2)
         call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

         call mpas_pool_get_array(statePool, 'uReconstructX', uReconstructX, 2)
         call mpas_pool_get_array(statePool, 'uReconstructY', uReconstructY, 2)
         call mpas_pool_get_array(statePool, 'uReconstructZ', uReconstructZ, 2)
         call mpas_pool_get_array(statePool, 'uReconstructZonal', uReconstructZonal, 2)
         call mpas_pool_get_array(statePool, 'uReconstructMeridional', uReconstructMeridional, 2)

         do iCell = 1, nCells
            do k = 1, nVertLevels
               tracersNew(:,k,iCell) = tracersNew(:,k,iCell) / hNew(k,iCell)
            end do
         end do

         if (config_test_case == 1) then    ! For case 1, wind field should be fixed
            uNew(:,:) = uOld(:,:)
         end if

         call sw_compute_solve_diagnostics(dt, statePool, meshPool, 2)

         call mpas_dmpar_field_halo_exch(domain, 'pv_edge', timeLevel=1)

         call mpas_reconstruct(meshPool, uNew,          &
                          uReconstructX, uReconstructY, uReconstructZ, &
                          uReconstructZonal, uReconstructMeridional )

         block => block % next
      end do

      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

         call mpas_pool_destroy_pool(provisStatePool)

         call mpas_pool_remove_subpool(block % structs, 'provis_state')
         block => block % next
      end do



   end subroutine sw_rk4

   subroutine sw_ssprk(domain, dt, ssprk_order)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   a Strong Stability Preserving Runge Kutta method of order ssprk_order
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: ssprk_order

      integer :: iCell, k, err
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tendPool

      ! BEGIN: TO REMOVE
      type (mpas_pool_type), pointer :: LTSPool
      integer, dimension(:,:), pointer :: nCellsInLTSHalo, nEdgesInLTSHalo
      integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
      ! END: TO REMOVE

      integer :: ssprk_step

      real (kind=RKIND), dimension(ssprk_order) :: weights_old, weights_new, weights_tend

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels

      real (kind=RKIND), dimension(:,:), pointer :: uOld, uNew, uTend
      real (kind=RKIND), dimension(:,:), pointer :: hOld, hNew, hTend
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersOld, tracersNew, tracersTend

      integer, pointer :: config_test_case
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4

      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)
      call mpas_pool_get_config(domain % configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)

     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(statePool, 'h', hOld, 1)
        call mpas_pool_get_array(statePool, 'h', hNew, 2)
        call mpas_pool_get_array(statePool, 'u', uOld, 1)
        call mpas_pool_get_array(statePool, 'u', uNew, 2)
        call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)
        call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

        uNew(:,:) = uOld(:,:)
        hNew(:,:) = hOld(:,:)
        do iCell = 1, nCells  
          do k = 1, nVertLevels
            tracersNew(:,k,iCell) = tracersOld(:,k,iCell) 
           end do
        end do

        call mpas_pool_initialize_time_levels(statePool)

        block => block % next
     end do

     weights_old(1) = 1.
     if (ssprk_order == 2) then 
        weights_old(2) = 0.5
     else if (ssprk_order == 3) then
        weights_old(2) = 0.75
        weights_old(3) = 1./3.
     end if 


     weights_new(1) = 0.
     if (ssprk_order == 2) then
        weights_new(2) = 0.5
     else if (ssprk_order == 3) then
        weights_new(2) = 0.25
        weights_new(3) = 2./3.
     end if 


     weights_tend(1) = dt * 1.
     if (ssprk_order == 2) then 
        weights_tend(2) = dt * 0.5
     else if (ssprk_order == 3) then 
     weights_tend(2) = dt * 0.25
     weights_tend(3) = dt * 2./3.
     end if

     err = 0 ! for MPI_Barrier tests

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     ! BEGIN SSPRK loop 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     do ssprk_step = 1, ssprk_order
   
call mpas_timer_start("halo updates SSPRK")
! --- update halos for diagnostic variables
        !call mpas_dmpar_field_halo_exch(domain, 'pv_edge', timeLevel=1)

        if (config_h_mom_eddy_visc4 > 0.0) then
            call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
            call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
        end if
call mpas_timer_stop("halo updates SSPRK")

! --- compute tendencies

call mpas_timer_start("tendencies SSPRK")
       block => domain % blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block % structs, 'state', statePool)

          call sw_compute_tend(tendPool, statePool, meshPool, 2)
          call sw_compute_scalar_tend(tendPool, statePool, meshPool, 2)
          call sw_enforce_boundary_edge(tendPool, meshPool)
          block => block % next
       end do
call mpas_timer_stop("tendencies SSPRK")

! --- update halos for prognostic variables

call mpas_timer_start("barrier time")
call MPI_Barrier(domain % dminfo % comm, err)
call mpas_timer_stop("barrier time")
call mpas_timer_start("halo updates SSPRK")
       call mpas_dmpar_field_halo_exch(domain, 'tend_u')
       call mpas_dmpar_field_halo_exch(domain, 'tend_h')
       call mpas_dmpar_field_halo_exch(domain, 'tend_tracers')
call mpas_timer_stop("halo updates SSPRK")

! --- advance solution

call mpas_timer_start("solution advancement SSPRK")
       block => domain % blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

          ! BEGIN: TO REMOVE
          call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
          call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
          call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
          call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
          call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
          ! END: TO REMOVE

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(statePool, 'u', uNew, 2)
          call mpas_pool_get_array(statePool, 'h', hNew, 2)
          call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

          call mpas_pool_get_array(statePool, 'u', uOld, 1)
          call mpas_pool_get_array(statePool, 'h', hOld, 1)
          call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

          call mpas_pool_get_array(tendPool, 'u', uTend)
          call mpas_pool_get_array(tendPool, 'h', hTend)
          call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

          if (config_test_case /= 1) then ! For case 1, wind should be fixed
          uNew(:,:) = weights_old(ssprk_step) * uOld(:,:) + weights_new(ssprk_step) * uNew(:,:) + weights_tend(ssprk_step) * uTend(:,:) 
          end if

          do iCell = 1, nCells
             do k = 1, nVertLevels
                tracersNew(:,k,iCell) = tracersNew(:,k,iCell) * hNew(k,iCell)
                hNew(k,iCell) = weights_old(ssprk_step) * hOld(k,iCell) + weights_new(ssprk_step) * hNew(k,iCell) + weights_tend(ssprk_step) * hTend(k,iCell)
                tracersNew(:,k,iCell) = ( weights_old(ssprk_step) * tracersOld(:,k,iCell) * hOld(k,iCell)  + weights_new(ssprk_step) * tracersNew(:,k,iCell) + & 
                                          weights_tend(ssprk_step) * tracersTend(:,k,iCell) ) / hNew(k,iCell)
             end do
          end do


          ! BEGIN: TO REMOVE
          !if (ssprk_step == 3 .and. block % blockID == 2) then
          !   do k =1, nVertLevels
          !      do iCell = 1, nCellsInLTSHalo(1,2)
          !         print*, hNew(k,cellsInLTSHalo(1,2,iCell))
          !      end do
          !print*, '-----------------------------------'
          !      do iCell = 1, nCellsInLTSHalo(2,2)
          !         print*, hNew(k,cellsInLTSHalo(2,2,iCell))
          !      end do
          !   end do
          !end if
          ! END: TO REMOVE


call mpas_timer_start("tendencies SSPRK")
          call sw_compute_solve_diagnostics(dt, statePool, meshPool, 2)
call mpas_timer_stop("tendencies SSPRK")
          block => block % next
       end do

call mpas_timer_stop("solution advancement SSPRK")
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! END SSPRK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


call mpas_timer_start("clean-up time SSPRK")
      !
      !  A little clean up at the end: compute diagnostics for new state
      !
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(statePool, 'h', hNew, 2)
         call mpas_pool_get_array(statePool, 'u', uNew, 2)
         call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

         call mpas_pool_get_array(statePool, 'uReconstructX', uReconstructX, 2)
         call mpas_pool_get_array(statePool, 'uReconstructY', uReconstructY, 2)
         call mpas_pool_get_array(statePool, 'uReconstructZ', uReconstructZ, 2)
         call mpas_pool_get_array(statePool, 'uReconstructZonal', uReconstructZonal, 2)
         call mpas_pool_get_array(statePool, 'uReconstructMeridional', uReconstructMeridional, 2)

         !call sw_compute_solve_diagnostics(dt, statePool, meshPool, 2)


call mpas_timer_start("barrier time")
call MPI_Barrier(domain % dminfo % comm, err)
call mpas_timer_stop("barrier time")
call mpas_timer_start("last halo update SSPRK")         
         call mpas_dmpar_field_halo_exch(domain, 'pv_edge', timeLevel=1)
call mpas_timer_stop("last halo update SSPRK")

         call mpas_reconstruct(meshPool, uNew,          &
                          uReconstructX, uReconstructY, uReconstructZ, &
                          uReconstructZonal, uReconstructMeridional )

         block => block % next
      end do
call mpas_timer_stop("clean-up time SSPRK")

   end subroutine sw_ssprk

   subroutine sw_lts(domain, dt, lts_order, totalTimeAdvance, totalTimeTendencies, totalTimeBarrier, totalTimeInterface)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using
   !   the local time stepping (LTS) from "Conservative explicit local
   !   time-stepping schemes for the shallow water equations", by Hoang at al.
   !   Journal of Computational Physics 382, 152-176, (2019)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: lts_order
      real (kind=RKIND), intent(inout) :: totalTimeAdvance, totalTimeTendencies, totalTimeBarrier, totalTimeInterface

      integer :: iCell, iEdge, iRegion, k, ic, ie, im, M, nRegions, err 
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool, meshPool, tendPool, LTSPool
      type (mpas_pool_type), pointer :: tendSum1stPool, tendSum2ndPool, tendSum3rdPool
      type (mpas_pool_type), pointer :: prevTendSum1stPool, nextTendSum1stPool, prevTendSum2ndPool, nextTendSum2ndPool, prevTendSum3rdPool, nextTendSum3rdPool
 
      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, haveFine, haveCoarse, haveInterface

      real (kind=RKIND) :: startTendencies, finishTendencies, startAdvance, finishAdvance, startBarrier, finishBarrier, startInterface, finishInterface 
      real (kind=RKIND), dimension(:,:), pointer :: uOld, uNew, uTend, uTendSum1st, uTendSum2nd, uTendSum3rd, uFirstStage, uSecondStage
      real (kind=RKIND), dimension(:,:), pointer :: hOld, hNew, hTend, hTendSum1st, hTendSum2nd, hTendSum3rd, hFirstStage, hSecondStage
      real (kind=RKIND), dimension(:,:), pointer :: uMPIHaloUpdate, hMPIHaloUpdate
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersMPIHaloUpdate
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersOld, tracersNew, tracersTend, tracersTendSum1st, tracersTendSum2nd, tracersTendSum3rd
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersFirstStage, tracersSecondStage
      real (kind=RKIND) :: dtFine, alpha, alphaHat, beta, betaHat, gam, gamHat

      integer, dimension(:,:), pointer :: nCellsInLTSHalo, nEdgesInLTSHalo
      integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo

      real (kind=RKIND) :: weightOld, weightNew, weightTend
      real (kind=RKIND) :: weightTendSum1st, weightTendSum2nd, weightTendSum3rd

      integer, pointer :: config_test_case, config_dt_scaling_LTS
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4

      call mpas_pool_get_config(domain % configs, 'config_dt_scaling_LTS', config_dt_scaling_LTS)
      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)
      call mpas_pool_get_config(domain % configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)

      M = config_dt_scaling_LTS
      nRegions = 2

      err = 0

      dtFine = dt / M

      weightOld = 0.5
      weightNew = 0.5
      weightTend = 0.5

      weightTendSum1st = 0.5
      weightTendSum2nd = 0.5
      weightTendSum3rd = 0.0

      if (lts_order == 3) then
         weightOld = 0.75
         weightNew = 0.25
         weightTend = 0.25

         weightTendSum1st = 1.0 / 6.0
         weightTendSum2nd = 1.0 / 6.0
         weightTendSum3rd = 2.0 / 3.0
      end if

call mpas_timer_start("create extra pools LTS")
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

         allocate(tendSum1stPool)
         allocate(tendSum2ndPool)
         allocate(tendSum3rdPool)

         call mpas_pool_create_pool(tendSum1stPool)
         call mpas_pool_clone_pool(tendPool, tendSum1stPool, 1)
         call mpas_pool_create_pool(tendSum2ndPool)
         call mpas_pool_clone_pool(tendPool, tendSum2ndPool, 1)
         call mpas_pool_create_pool(tendSum3rdPool)
         call mpas_pool_clone_pool(tendPool, tendSum3rdPool, 1)

         call mpas_pool_add_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
         call mpas_pool_add_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
         call mpas_pool_add_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)

         call mpas_pool_get_array(statePool, 'h', hOld, 1)
         call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
         call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
         call mpas_pool_get_array(statePool, 'u', uOld, 1)
         call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
         call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
         call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)
         call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)
         call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)

         call mpas_pool_get_array(tendSum1stPool, 'u', uTendSum1st)
         call mpas_pool_get_array(tendSum1stPool, 'h', hTendSum1st)
         call mpas_pool_get_array(tendSum1stPool, 'tracers', tracersTendSum1st)

         call mpas_pool_get_array(tendSum2ndPool, 'u', uTendSum2nd)
         call mpas_pool_get_array(tendSum2ndPool, 'h', hTendSum2nd)
         call mpas_pool_get_array(tendSum2ndPool, 'tracers', tracersTendSum2nd)
         
         call mpas_pool_get_array(tendSum3rdPool, 'u', uTendSum3rd)
         call mpas_pool_get_array(tendSum3rdPool, 'h', hTendSum3rd)
         call mpas_pool_get_array(tendSum3rdPool, 'tracers', tracersTendSum3rd)

         hFirstStage(:,:) = hOld(:,:)
         hSecondStage(:,:) = hOld(:,:)

         uFirstStage(:,:) = uOld(:,:)
         uSecondStage(:,:) = uOld(:,:)

         tracersFirstStage(:,:,:) = tracersOld(:,:,:)
         tracersSecondStage(:,:,:) = tracersOld(:,:,:)

         uTendSum1st(:,:) = 0
         hTendSum1st(:,:) = 0
         tracersTendSum1st(:,:,:) = 0
         
         uTendSum2nd(:,:) = 0
         hTendSum2nd(:,:) = 0
         tracersTendSum2nd(:,:,:) = 0
         
         uTendSum3rd(:,:) = 0
         hTendSum3rd(:,:) = 0
         tracersTendSum3rd(:,:,:) = 0

         call mpas_pool_initialize_time_levels(statePool)

         block => block % next
      end do

      block => domain % blocklist
      do while(associated(block))
         if (associated(block % prev)) then
            call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_1st', tendSum1stPool)
            call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_2nd', tendSum2ndPool)
         else
            nullify(prevTendSum1stPool)
            nullify(prevTendSum2ndPool)
         end if

         if (associated(block % next)) then
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_1st', nextTendSum1stPool)
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_2nd', nextTendSum2ndPool)
         else
            nullify(nextTendSum1stPool)
            nullify(nextTendSum2ndPool)
         end if

         call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
         call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)

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

         call mpas_pool_link_parinfo(block, tendSum1stPool)
         call mpas_pool_link_parinfo(block, tendSum2ndPool)

         if(lts_order == 3) then

         if (associated(block % prev)) then
            call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_3rd', prevTendSum3rdPool)
         else
            nullify(prevTendSum3rdPool)
         end if

         if (associated(block % next)) then
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_3rd', nextTendSum3rdPool)
         else
            nullify(nextTendSum3rdPool)
         end if

         call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)

         if (associated(prevTendSum3rdPool) .and. associated(nextTendSum3rdPool)) then
            call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool, nextTendSum3rdPool)
         else if (associated(prevTendSum3rdPool)) then
            call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool)
         else if (associated(nextTendSum3rdPool)) then
            call mpas_pool_link_pools(tendSum3rdPool,nextPool=nextTendSum3rdPool)
         else
            call mpas_pool_link_pools(tendSum3rdPool)
         end if

         call mpas_pool_link_parinfo(block, tendSum3rdPool)

         end if

         block => block % next
     end do

call mpas_timer_stop("create extra pools LTS")

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! BEGIN LTS SCHEME
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! --- compute the first stage of SSPRK for interface layer 1, interface layer 2 and coarse 

call mpas_timer_start("halo updates  LTS")
     ! --- update halos for diagnostic variables

     if (config_h_mom_eddy_visc4 > 0.0) then  ! ALWAYS USED 0.0 HERE FOR THIS PARAMETER (so it has not been tested)
        call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
        call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
     end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
     ! --- compute tendencies for interface layer 1, interface layer 2 and coarse
     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

        call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
        call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
        call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

        if(haveCoarse == 1 .or. haveInterface == 1) then
        call sw_compute_solve_diagnostics(dt, statePool, meshPool, 1) 
        call sw_compute_tend(tendPool, statePool, meshPool, 1) 
        call sw_compute_scalar_tend(tendPool, statePool, meshPool, 1)  
        call sw_enforce_boundary_edge(tendPool, meshPool) 
        end if
        if (lts_order == 3 .and. haveFine == 1) then
        call sw_compute_solve_diagnostics(dt, statePool, meshPool, 1) 
        call sw_compute_tend(tendPool, statePool, meshPool, 1) 
        call sw_compute_scalar_tend(tendPool, statePool, meshPool, 1)  
        call sw_enforce_boundary_edge(tendPool, meshPool) 
        end if

        block => block % next
     end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
     ! --- advance solution for interface layer 1, interface layer 2 and coarse with first stage of SSPRK
     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
        call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
        call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

        call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
        call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
        call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

        call mpas_pool_get_array(statePool, 'u', uOld, 1)
        call mpas_pool_get_array(statePool, 'h', hOld, 1)
        call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

        call mpas_pool_get_array(tendPool, 'u', uTend)
        call mpas_pool_get_array(tendPool, 'h', hTend)
        call mpas_pool_get_array(tendPool, 'tracers', tracersTend)
       
        if(haveFine == 1 .or. haveCoarse == 1 .or. haveInterface == 1) then
        do k =1, nVertLevels
           do iRegion =1,nRegions
              do ie = 1, nEdgesInLTSHalo(iRegion,2)
                 iEdge = edgesInLTSHalo(iRegion,2,ie)
                 uFirstStage(k,iEdge) =  uOld(k,iEdge) + dt * uTend(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(iRegion,2)
                 iCell = cellsInLTSHalo(iRegion,2,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dt * hTend(k,iCell)
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dt * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
              end do
           end do
           do ie = 1, nEdgesInLTSHalo(2,1)
              iEdge = edgesInLTSHalo(2,1,ie)
              uFirstStage(k,iEdge) =  uOld(k,iEdge) + dt * uTend(k,iEdge)
           end do
           do ic = 1, nCellsInLTSHalo(2,1)
                 iCell = cellsInLTSHalo(2,1,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dt * hTend(k,iCell)
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dt * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
           end do

           if (lts_order == 3) then
              do ie = 1, nEdgesInLTSHalo(2,3) 
                 iEdge = edgesInLTSHalo(2,3,ie)
                 uFirstStage(k,iEdge) =  uOld(k,iEdge) + dt * uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(2,4) 
                 iEdge = edgesInLTSHalo(2,4,ie)
                 uFirstStage(k,iEdge) =  uOld(k,iEdge) + dt * uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(1,3) 
                 iEdge = edgesInLTSHalo(1,3,ie)
                 uFirstStage(k,iEdge) =  uOld(k,iEdge) + dt * uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(1,4) 
                 iEdge = edgesInLTSHalo(1,4,ie)
                 uFirstStage(k,iEdge) =  uOld(k,iEdge) + dt * uTend(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(2,3) !these are coarse cells that have a region of their own if nLTSHalos=3
                 iCell = cellsInLTSHalo(2,3,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dt * hTend(k,iCell)
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dt * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(2,4) !these are coarse cells that have a region of their own if nLTSHalos=3
                 iCell = cellsInLTSHalo(2,4,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dt * hTend(k,iCell)
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dt * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(1,3) !these are fine cells that have a region of their own if nLTSHalos=3
                 iCell = cellsInLTSHalo(1,3,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dt * hTend(k,iCell)
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dt * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(1,4) !these are fine cells that have a region of their own if nLTSHalos=3
                 iCell = cellsInLTSHalo(1,4,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dt * hTend(k,iCell)
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dt * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
              end do
           end if

        end do

        if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uFirstStage(:,:) = uOld(:,:)
        end if

        end if

        !BEGIN: TO REMOVE
        !if (block % blockID == 2) then
        !print*, 'LTS - First Stage - Coarse'
        !do k =1, nVertLevels
        !   do ie = 1, nEdgesInLTSHalo(2,1)
        !      iEdge = edgesInLTSHalo(2,1,ie)
        !      !print*, 'u=', uFirstStage(k,iEdge), 'edge=', iEdge
        !   end do
        !   do ic = 1, nCellsInLTSHalo(2,2)
        !       iCell = cellsInLTSHalo(2,2,ic)
        !      print*,  hFirstStage(k,iCell)
        !   end do 
        !   do ic = 1, nCellsInLTSHalo(1,2)
        !       iCell = cellsInLTSHalo(1,2,ic)
        !      print*,  hFirstStage(k,iCell)
        !   end do
        !end do
        !end if
        !END: TO REMOVE

        block => block % next
     end do
call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")

call mpas_timer_start("halo updates  LTS")
     call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=3)
     call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=3)
     call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=3)
call mpas_timer_stop("halo updates  LTS")

     ! ---------------------------------------------------------------------------------------------------------------------

     ! --- compute the second stage of SSPRK for coarse

call mpas_timer_start("halo updates  LTS")
     ! --- update halos for diagnostic variables

     if (config_h_mom_eddy_visc4 > 0.0) then
        call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
        call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
     end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
     ! --- compute tendencies for coarse
     block => domain % blocklist
     do while (associated(block))
 
        call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

        call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

        if(haveCoarse == 1) then
        call sw_compute_solve_diagnostics(dt, statePool, meshPool, 3) 
        call sw_compute_tend(tendPool, statePool, meshPool, 3) 
        call sw_compute_scalar_tend(tendPool, statePool, meshPool, 3) 
        call sw_enforce_boundary_edge(tendPool, meshPool)
        end if

        block => block % next
     end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
     ! --- advance solution (2nd stage of SSPRK for the coarse using dt)
     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

        call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
        call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
        call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

        call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
        call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
        call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)
        
        call mpas_pool_get_array(statePool, 'u', uOld, 1)
        call mpas_pool_get_array(statePool, 'h', hOld, 1)
        call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

        call mpas_pool_get_array(tendPool, 'u', uTend)
        call mpas_pool_get_array(tendPool, 'h', hTend)
        call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

        if(haveCoarse == 1) then
        do k =1, nVertLevels
           do ie = 1, nEdgesInLTSHalo(2,1)
              iEdge = edgesInLTSHalo(2,1,ie)
              uSecondStage(k,iEdge) = weightOld * uOld(k,iEdge) + weightNew * uFirstStage(k,iEdge) + weightTend * dt * uTend(k,iEdge)
           end do
           do ic = 1, nCellsInLTSHalo(2,1)
              iCell = cellsInLTSHalo(2,1,ic)
              hSecondStage(k,iCell) = weightOld * hOld(k,iCell) + weightNew * hFirstStage(k,iCell) + weightTend * dt * hTend(k,iCell)
              tracersSecondStage(:,k,iCell) = ( weightOld * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                                weightNew * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + &
                                                weightTend * dt * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell) 
           end do

           if (lts_order == 3) then !these are coarse cells that have a region of their own if nLTSHalos=2
              do ie = 1, nEdgesInLTSHalo(2,3)
                 iEdge = edgesInLTSHalo(2,3,ie)
                 uSecondStage(k,iEdge) = weightOld * uOld(k,iEdge) + weightNew * uFirstStage(k,iEdge) + weightTend * dt * uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(2,4)
                 iEdge = edgesInLTSHalo(2,4,ie)
                 uSecondStage(k,iEdge) = weightOld * uOld(k,iEdge) + weightNew * uFirstStage(k,iEdge) + weightTend * dt * uTend(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(2,3)
                 iCell = cellsInLTSHalo(2,3,ic)
                 hSecondStage(k,iCell) = weightOld * hOld(k,iCell) + weightNew * hFirstStage(k,iCell) + weightTend * dt * hTend(k,iCell)
                 tracersSecondStage(:,k,iCell) = ( weightOld * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                                   weightNew * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + &
                                                   weightTend * dt * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(2,4)
                 iCell = cellsInLTSHalo(2,4,ic)
                 hSecondStage(k,iCell) = weightOld * hOld(k,iCell) + weightNew * hFirstStage(k,iCell) + weightTend * dt * hTend(k,iCell)
                 tracersSecondStage(:,k,iCell) = ( weightOld * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                                   weightNew * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + &
                                                   weightTend * dt * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell)
              end do

           end if

        end do

        if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uSecondStage(:,:) = uOld(:,:)
        end if

        end if

        !BEGIN: TO REMOVE
        !if (block % blockID == 1) then
        !print*, 'LTS - Second Stage - Coarse'
        !do k =1, nVertLevels
        !   do ie = 1, nEdgesInLTSHalo(2,1)
        !      iEdge = edgesInLTSHalo(2,1,ie)
        !      print*, 'u=', uSecondStage(k,iEdge), 'edge=', iEdge
        !   end do
        !   do ic = 1, nCellsInLTSHalo(2,1)
        !      iCell = cellsInLTSHalo(2,1,ic)
        !      print*, hSecondStage(k,iCell)
        !   end do
        !   do ic = 1, nCellsInLTSHalo(2,3)
        !      iCell = cellsInLTSHalo(2,3,ic)
        !      print*, hSecondStage(k,iCell)
        !   end do
        !end do
        !end if
        !END: TO REMOVE

        block => block % next
     end do

call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")

call mpas_timer_start("halo updates  LTS")
     call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=4)
     call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=4)
     call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=4)
call mpas_timer_stop("halo updates  LTS")

     ! ---------------------------------------------------------------------------------------------------------------------

     ! --- compute the second stage of SSPRK only for interface 1 and interface 2 (this is needed here only if lts_order = 3)
     if (lts_order == 3) then

call mpas_timer_start("halo updates  LTS")
     ! --- update halos for diagnostic variables

     if (config_h_mom_eddy_visc4 > 0.0) then
        call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
        call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
     end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
     ! --- compute tendencies for interface
     block => domain % blocklist
     do while (associated(block))

        call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

        call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

       if(haveInterface == 1) then
        call sw_compute_solve_diagnostics(dt, statePool, meshPool, 3)
        call sw_compute_tend(tendPool, statePool, meshPool, 3)
        call sw_compute_scalar_tend(tendPool, statePool, meshPool, 3)
        call sw_enforce_boundary_edge(tendPool, meshPool)
        end if

        block => block % next
     end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
        ! --- advance solution for interface layer 1, interface layer 2 with second stage of SSPRK
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
           call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
           call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

           call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
           call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
           call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)       

           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           call mpas_pool_get_array(tendPool, 'u', uTend)
           call mpas_pool_get_array(tendPool, 'h', hTend)
           call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

           if (haveInterface == 1) then
           do k =1, nVertLevels
              do iRegion =1,nRegions
                 do ie = 1, nEdgesInLTSHalo(iRegion,2)
                    iEdge = edgesInLTSHalo(iRegion,2,ie)
                    uSecondStage(k,iEdge) =  0.75 * uOld(k,iEdge) + 0.25 * uFirstStage(k,iEdge) + 0.25 * dt * uTend(k,iEdge)
                 end do
                 do ic = 1, nCellsInLTSHalo(iRegion,2)
                    iCell = cellsInLTSHalo(iRegion,2,ic)
                    hSecondStage(k,iCell) = 0.75 * hOld(k,iCell) + 0.25 * hFirstStage(k,iCell) + 0.25 * dt * hTend(k,iCell)
                    tracersSecondStage(:,k,iCell) = ( 0.75 * tracersOld(:,k,iCell) * hOld(k,iCell) + &
                                                    0.25 * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + &
                                                    0.25 * dt * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell)
                 end do
              end do
           end do

           if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uSecondStage(:,:) = uOld(:,:)
           end if

           end if

        !BEGIN: TO REMOVE
        !if (block % blockID == 2) then
        !print*, 'LTS - Second Stage - Interface'
        !do k =1, nVertLevels
        !   do ie = 1, nEdgesInLTSHalo(2,1)
        !      iEdge = edgesInLTSHalo(2,1,ie)
        !      print*, 'u=', uSecondStage(k,iEdge), 'edge=', iEdge
        !   end do
        !   do ic = 1, nCellsInLTSHalo(2,2)
        !      iCell = cellsInLTSHalo(2,2,ic)
        !      print*, hSecondStage(k,iCell)
        !   end do
        !   do ic = 1, nCellsInLTSHalo(1,2)
        !      iCell = cellsInLTSHalo(1,2,ic)
        !      print*, hSecondStage(k,iCell)
        !   end do
        !end do
        !end if
        !END: TO REMOVE
          
           block => block % next
        end do

call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")

call mpas_timer_start("halo updates  LTS")
     call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=4)
     call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=4)
     call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=4)
call mpas_timer_stop("halo updates  LTS")

     end if

     ! ---------------------------------------------------------------------------------------------------------------------

     ! --- FINE AND COARSE ADVANCEMENT STEP: all stages for fine, second (and third if needed) for coarse
  
     ! --- begin substepping for fine solution advancement
     do im =1, M


call mpas_timer_start("interface prediction LTS")
call cpu_time(startInterface)
        ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln old in soln new
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
           call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
           call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

           call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
           call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
           call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)
 
           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           if (haveInterface == 1) then
           alpha = REAL((im - 1)) / M
           alphaHat = 0.0
           if (lts_order == 3) then
              alphaHat = REAL((im - 1) * (im - 1)) / (M * M)
           end if
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(1,2) !interface layer 1
                 iEdge = edgesInLTSHalo(1,2,ie)
                 uNew(k,iEdge) = uOld(k,iEdge) 
                 uOld(k,iEdge) = (1.0 - alpha - alphaHat) * uOld(k,iEdge) + (alpha - alphaHat) * uFirstStage(k,iEdge) + 2.0 * alphaHat * uSecondStage(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(1,2) !interface layer 1
                 iCell = cellsInLTSHalo(1,2,ic)
                 hNew(k,iCell) = hOld(k,iCell)
                 tracersNew(:,k,iCell) = tracersOld(:,k,iCell)
                 hOld(k,iCell) = (1.0 - alpha - alphaHat) * hOld(k,iCell) + (alpha - alphaHat) * hFirstStage(k,iCell) + 2.0 * alphaHat * hSecondStage(k,iCell)
                 tracersOld(:,k,iCell) = (1.0 - alpha - alphaHat) * tracersOld(:,k,iCell) + & 
                                         (alpha - alphaHat) * tracersFirstStage(:,k,iCell) + & 
                                          2.0 * alphaHat * tracersSecondStage(:,k,iCell) 
              end do
           end do
           if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uOld(:,:) = uNew(:,:)
           end if
           end if

           block => block % next
        end do
call cpu_time(finishInterface)
totalTimeInterface = totalTimeInterface + (finishInterface - startInterface)
call mpas_timer_stop("interface prediction LTS")

call mpas_timer_start("halo updates  LTS")
        ! --- update halos for diagnostic variables

        if (config_h_mom_eddy_visc4 > 0.0) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
           call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
        end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
        ! --- compute tendencies for interface layer 1, interface layer 2 and fine
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

           call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           if(haveFine == 1 .or. haveInterface == 1) then
           call sw_compute_solve_diagnostics(dt, statePool, meshPool, 1) 
           call sw_compute_tend(tendPool, statePool, meshPool, 1) 
           call sw_compute_scalar_tend(tendPool, statePool,meshPool, 1) 
           call sw_enforce_boundary_edge(tendPool, meshPool)
           end if

           block => block % next
        end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
        ! --- advance fine solution (first stage of SSPRK) AND sum up the tendecies above into tendSum2nd for interface layer 1 and interface layer 2
        ! --- also restore soln old by copying it back from soln new
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
           call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           call mpas_pool_get_array(tendPool, 'u', uTend)
           call mpas_pool_get_array(tendPool, 'h', hTend)
           call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

           call mpas_pool_get_array(tendSum2ndPool, 'u', uTendSum2nd)
           call mpas_pool_get_array(tendSum2ndPool, 'h', hTendSum2nd)
           call mpas_pool_get_array(tendSum2ndPool, 'tracers', tracersTendSum2nd)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
           call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
           call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

           if(haveFine == 1 .or. haveInterface == 1) then
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(1,2)
                 iEdge = edgesInLTSHalo(1,2,ie)
                 uOld(k,iEdge) = uNew(k,iEdge)
                 uTendSum2nd(k,iEdge) = uTendSum2nd(k,iEdge) + uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(2,2)
                 iEdge = edgesInLTSHalo(2,2,ie)
                 uTendSum2nd(k,iEdge) = uTendSum2nd(k,iEdge) + uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(1,1)
                 iEdge = edgesInLTSHalo(1,1,ie)
                 uFirstStage(k,iEdge) =  uOld(k,iEdge) + dtFine * uTend(k,iEdge) !soln update for the fine
              end do
              do ic = 1, nCellsInLTSHalo(1,2)
                 iCell = cellsInLTSHalo(1,2,ic)
                 hOld(k,iCell) = hNew(k,iCell)
                 tracersOld(:,k,iCell) = tracersNew(:,k,iCell)
                 hTendSum2nd(k,iCell) = hTendSum2nd(k,iCell) +  hTend(k,iCell)
                 tracersTendSum2nd(:,k,iCell) = tracersTendSum2nd(:,k,iCell) + tracersTend(:,k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(2,2)
                 iCell = cellsInLTSHalo(2,2,ic)
                 hTendSum2nd(k,iCell) = hTendSum2nd(k,iCell) +  hTend(k,iCell)
                 tracersTendSum2nd(:,k,iCell) = tracersTendSum2nd(:,k,iCell) + tracersTend(:,k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(1,1)
                 iCell = cellsInLTSHalo(1,1,ic)
                 hFirstStage(k,iCell) = hOld(k,iCell) + dtFine * hTend(k,iCell) !soln update for the fine
                 tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dtFine * tracersTend(:,k,iCell) ) / hFirstStage(k,iCell)
              end do

              if (lts_order == 3) then
                 do ie = 1, nEdgesInLTSHalo(1,3)
                    iEdge = edgesInLTSHalo(1,3,ie)
                    uFirstStage(k,iEdge) =  uOld(k,iEdge) + dtFine * uTend(k,iEdge) !soln update for the fine
                 end do
                 do ie = 1, nEdgesInLTSHalo(1,4)
                    iEdge = edgesInLTSHalo(1,4,ie)
                    uFirstStage(k,iEdge) =  uOld(k,iEdge) + dtFine * uTend(k,iEdge) !soln update for the fine
                 end do
                 do ic = 1, nCellsInLTSHalo(1,3)
                    iCell = cellsInLTSHalo(1,3,ic)
                    hFirstStage(k,iCell) = hOld(k,iCell) + dtFine * hTend(k,iCell) !soln update for the fine
                    tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dtFine * tracersTend(:,k,iCell)) / hFirstStage(k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,4)
                    iCell = cellsInLTSHalo(1,4,ic)
                    hFirstStage(k,iCell) = hOld(k,iCell) + dtFine * hTend(k,iCell) !soln update for the fine
                    tracersFirstStage(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + dtFine * tracersTend(:,k,iCell)) / hFirstStage(k,iCell)
                 end do
              end if

           end do

           if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uFirstStage(:,:) = uOld(:,:)
           end if

           end if   

           !BEGIN: TO REMOVE
           !if (block % blockID == 0) then
           !print*, 'LTS - First Stage - Fine'
           !do k =1, nVertLevels
           !   do ie = 1, nEdgesInLTSHalo(1,1)
           !      iEdge = edgesInLTSHalo(1,1,ie)
           !      print*, 'u=', uNew(k,iEdge), 'edge=', iEdge
           !   end do
           !   do ic = 1, nCellsInLTSHalo(1,1)
           !      iCell = cellsInLTSHalo(1,1,ic)
           !      print*, hFirstStage(k,iCell)
           !   end do
           !   do ic = 1, nCellsInLTSHalo(1,3)
           !      iCell = cellsInLTSHalo(1,3,ic)
           !      print*, hFirstStage(k,iCell)
           !   end do
           !end do
           !end if
           !END: TO REMOVE

           block => block % next
        end do

call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")

call mpas_timer_start("halo updates  LTS")
     call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=3)
     call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=3)
     call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=3)
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("interface prediction LTS")
call cpu_time(startInterface)
        ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln first stage in soln new
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
           call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
           call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

           call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
           call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
           call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)
 
           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           if (haveInterface == 1) then
           beta = REAL(im) / M
           betaHat = 0.0
           if (lts_order == 3) then
              betaHat = REAL((im - 1) * (im + 1)) / (M * M)
           end if
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(1,2) !interface layer 1
                 iEdge = edgesInLTSHalo(1,2,ie)
                 uNew(k,iEdge) = uFirstStage(k,iEdge)
                 uFirstStage(k,iEdge) = (1.0 - beta - betaHat) * uOld(k,iEdge) + (beta - betaHat) * uFirstStage(k,iEdge) + 2.0 * betaHat * uSecondStage(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(1,2) !interface layer 1
                 iCell = cellsInLTSHalo(1,2,ic)
                 hNew(k,iCell) = hFirstStage(k,iCell)
                 tracersNew(:,k,iCell) = tracersFirstStage(:,k,iCell)
                 hFirstStage(k,iCell) = (1.0 - beta - betaHat) * hOld(k,iCell) + (beta - betaHat) * hFirstStage(k,iCell) + 2.0 * betaHat * hSecondStage(k,iCell)
                 tracersFirstStage(:,k,iCell) = (1.0 - beta - betaHat) * tracersOld(:,k,iCell) + &
                                         (beta - betaHat) * tracersFirstStage(:,k,iCell) + & 
                                          2.0 * betaHat * tracersSecondStage(:,k,iCell)
              end do
           end do
           if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uFirstStage(:,:) = uOld(:,:)
           end if
           end if

           block => block % next
        end do
call cpu_time(finishInterface)
totalTimeInterface = totalTimeInterface + (finishInterface - startInterface)
call mpas_timer_stop("interface prediction LTS")

call mpas_timer_start("halo updates  LTS")
        ! --- update halos for diagnostic variables

        if (config_h_mom_eddy_visc4 > 0.0) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
           call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
        end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
        ! --- compute tendencies for interface layer 1, interface layer 2 and fine
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)

           call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           if(haveFine == 1 .or. haveInterface == 1) then
           call sw_compute_solve_diagnostics(dt, statePool, meshPool, 3) 
           call sw_compute_tend(tendPool, statePool, meshPool, 3) 
           call sw_compute_scalar_tend(tendPool, statePool, meshPool, 3) 
           call sw_enforce_boundary_edge(tendPool, meshPool)
           end if
           block => block % next
        end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
        ! --- advance fine solution (second stage of SSPRK) AND sum up the tendecies above into tendSum1st for interface layer 1 and interface layer 2
        ! --- also restore soln first stage by copying it back from soln new
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
           call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
           call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

           call mpas_pool_get_array(tendPool, 'u', uTend)
           call mpas_pool_get_array(tendPool, 'h', hTend)
           call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

           call mpas_pool_get_array(tendSum1stPool, 'u', uTendSum1st)
           call mpas_pool_get_array(tendSum1stPool, 'h', hTendSum1st)
           call mpas_pool_get_array(tendSum1stPool, 'tracers', tracersTendSum1st)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
           call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
           call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

           call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
           call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
           call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)
           
           if(haveFine == 1 .or. haveInterface == 1) then
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(1,2)
                 iEdge = edgesInLTSHalo(1,2,ie)
                 uFirstStage(k,iEdge) = uNew(k,iEdge)
                 uTendSum1st(k,iEdge) = uTendSum1st(k,iEdge) + uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(2,2)
                 iEdge = edgesInLTSHalo(2,2,ie)
                 uTendSum1st(k,iEdge) = uTendSum1st(k,iEdge) + uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(1,1)
                 iEdge = edgesInLTSHalo(1,1,ie)
                 uSecondStage(k,iEdge) =  weightOld * uOld(k,iEdge) + weightNew * uFirstStage(k,iEdge) + weightTend *  dtFine * uTend(k,iEdge) !soln update for the fine
              end do
              do ic = 1, nCellsInLTSHalo(1,2)
                 iCell = cellsInLTSHalo(1,2,ic)
                 hFirstStage(k,iCell) = hNew(k,iCell)
                 tracersFirstStage(:,k,iCell) = tracersNew(:,k,iCell)
                 hTendSum1st(k,iCell) = hTendSum1st(k,iCell) +  hTend(k,iCell)
                 tracersTendSum1st(:,k,iCell) = tracersTendSum1st(:,k,iCell) + tracersTend(:,k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(2,2)
                 iCell = cellsInLTSHalo(2,2,ic)
                 hTendSum1st(k,iCell) = hTendSum1st(k,iCell) + hTend(k,iCell)
                 tracersTendSum1st(:,k,iCell) = tracersTendSum1st(:,k,iCell) + tracersTend(:,k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(1,1)
                 iCell = cellsInLTSHalo(1,1,ic)
                 hSecondStage(k,iCell) = weightOld * hOld(k,iCell) + weightNew * hFirstStage(k,iCell) + weightTend * dtFine * hTend(k,iCell) !soln update for the fine
                 tracersSecondStage(:,k,iCell) = ( weightOld * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                                   weightNew * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + & 
                                                   weightTend * dtFine * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell)
              end do
           end do

           if (lts_order == 3) then
              do k =1, nVertLevels
                 do ie = 1, nEdgesInLTSHalo(1,3)
                   iEdge = edgesInLTSHalo(1,3,ie)
                   uSecondStage(k,iEdge) =  weightOld * uOld(k,iEdge) + weightNew * uFirstStage(k,iEdge) + weightTend * dtFine * uTend(k,iEdge)!soln update for the fine
                 end do
                 do ie = 1, nEdgesInLTSHalo(1,4)
                   iEdge = edgesInLTSHalo(1,4,ie)
                   uSecondStage(k,iEdge) =  weightOld * uOld(k,iEdge) + weightNew * uFirstStage(k,iEdge) + weightTend * dtFine * uTend(k,iEdge)!soln update for the fine
                 end do
                 do ic = 1, nCellsInLTSHalo(1,3)
                    iCell = cellsInLTSHalo(1,3,ic)
                    hSecondStage(k,iCell) = weightOld * hOld(k,iCell) + weightNew * hFirstStage(k,iCell) + weightTend * dtFine * hTend(k,iCell) !soln update for the fine
                    tracersSecondStage(:,k,iCell) = ( weightOld * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                                      weightNew * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + & 
                                                      weightTend * dtFine * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,4)
                    iCell = cellsInLTSHalo(1,4,ic)
                    hSecondStage(k,iCell) = weightOld * hOld(k,iCell) + weightNew * hFirstStage(k,iCell) + weightTend * dtFine * hTend(k,iCell) !soln update for the fine
                    tracersSecondStage(:,k,iCell) = ( weightOld * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                                      weightNew * tracersFirstStage(:,k,iCell) * hFirstStage(k,iCell) + & 
                                                      weightTend * dtFine * tracersTend(:,k,iCell) ) / hSecondStage(k,iCell)
                 end do
              end do
           end if

           if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uSecondStage(:,:) = uOld(:,:)
           end if

           end if

           !BEGIN: TO REMOVE
           !if (block % blockID == 0) then
           !print*, 'LTS - Second Stage - Fine'
           !do k =1, nVertLevels
           !   do ie = 1, nEdgesInLTSHalo(1,1)
           !      iEdge = edgesInLTSHalo(1,1,ie)
           !      print*, 'u=', uNew(k,iEdge), 'edge=', iEdge
           !   end do
           !   do ic = 1, nCellsInLTSHalo(1,1)
           !      iCell = cellsInLTSHalo(1,1,ic)
           !      print*, hSecondStage(k,iCell)
           !   end do
           !   do ic = 1, nCellsInLTSHalo(1,3)
           !      iCell = cellsInLTSHalo(1,3,ic)
           !      print*, hSecondStage(k,iCell)
           !   end do
           !end do
           !end if
           !END: TO REMOVE

           block => block % next
        end do

call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")

call mpas_timer_start("halo updates  LTS")
     call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=4)
     call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=4)
     call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=4)
call mpas_timer_stop("halo updates  LTS")

! --- copy second stage soln into old soln for the fine
call mpas_timer_start("copy second stage soln into old soln  LTS")
        if(lts_order == 2) then

           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

              call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

              call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)

              call mpas_pool_get_array(statePool, 'u', uOld, 1)
              call mpas_pool_get_array(statePool, 'h', hOld, 1)
              call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

              call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
              call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
              call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)

              if(haveFine == 1) then
                 do k =1, nVertLevels
                    do ie = 1, nEdgesInLTSHalo(1,1)
                       iEdge = edgesInLTSHalo(1,1,ie)
                       uOld(k,iEdge) =  uSecondStage(k,iEdge)
                    end do
                    do ic = 1, nCellsInLTSHalo(1,1)
                       iCell = cellsInLTSHalo(1,1,ic)
                       hOld(k,iCell) = hSecondStage(k,iCell)
                       tracersOld(:,k,iCell) = tracersSecondStage(:,k,iCell)
                    end do
                 end do
              end if

              block => block % next
           end do

call mpas_timer_start("halo updates  LTS")
           call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=1)
           call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=1)
           call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=1)
call mpas_timer_stop("halo updates  LTS")


         end if
call mpas_timer_stop("copy second stage soln into old soln  LTS")

        if (lts_order == 3) then

call mpas_timer_start("interface prediction LTS")
call cpu_time(startInterface)
           ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln second stage in soln new
           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

              call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

              call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

              call mpas_pool_get_array(statePool, 'u', uFirstStage, 3)
              call mpas_pool_get_array(statePool, 'h', hFirstStage, 3)
              call mpas_pool_get_array(statePool, 'tracers', tracersFirstStage, 3)

              call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
              call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
              call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)

              call mpas_pool_get_array(statePool, 'u', uOld, 1)
              call mpas_pool_get_array(statePool, 'h', hOld, 1)
              call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

              call mpas_pool_get_array(statePool, 'u', uNew, 2)
              call mpas_pool_get_array(statePool, 'h', hNew, 2)
              call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)
 
              if (haveInterface == 1) then
              gam = REAL(2 * (im - 1) + 1) / ( 2 * M )
              gamHat = REAL(2 * (im - 1) * (im - 1) + 2 * (im - 1) + 1) / (2 * M * M)
              do k =1, nVertLevels
                 do ie = 1, nEdgesInLTSHalo(1,2) !interface layer 1
                    iEdge = edgesInLTSHalo(1,2,ie)
                    uNew(k,iEdge) = uSecondStage(k,iEdge)
                    uSecondStage(k,iEdge) = (1.0 - gam - gamHat) * uOld(k,iEdge) + (gam - gamHat) * uFirstStage(k,iEdge) + 2.0 * gamHat * uSecondStage(k,iEdge)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,2) !interface layer 1
                    iCell = cellsInLTSHalo(1,2,ic)
                    hNew(k,iCell) = hSecondStage(k,iCell)
                    tracersNew(:,k,iCell) = tracersSecondStage(:,k,iCell)
                    hSecondStage(k,iCell) = (1.0 - gam - gamHat) * hOld(k,iCell) + (gam - gamHat) * hFirstStage(k,iCell) + 2.0 * gamHat * hSecondStage(k,iCell)
                    tracersSecondStage(:,k,iCell) = (1.0 - gam - gamHat) * tracersOld(:,k,iCell)  + &
                                            (gam - gamHat) * tracersFirstStage(:,k,iCell) + & 
                                             2.0 * gamHat * tracersSecondStage(:,k,iCell)
                 end do
              end do
              if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uSecondStage(:,:) = uNew(:,:)
              end if
              end if

              block => block % next
           end do
call cpu_time(finishInterface)
totalTimeInterface = totalTimeInterface + (finishInterface - startInterface)
call mpas_timer_stop("interface prediction LTS")

call mpas_timer_start("halo updates  LTS")
           ! --- update halos for diagnostic variables

           if (config_h_mom_eddy_visc4 > 0.0) then
              call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
              call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
           end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
           ! --- compute tendencies for interface layer 1, interface layer 2 and fine
           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
              call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)

              call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
              call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

              if(haveFine == 1 .or. haveInterface == 1) then
              call sw_compute_solve_diagnostics(dt, statePool, meshPool, 4) 
              call sw_compute_tend(tendPool, statePool, meshPool, 4) 
              call sw_compute_scalar_tend(tendPool, statePool, meshPool, 4) 
              call sw_enforce_boundary_edge(tendPool, meshPool)
              end if
              block => block % next
           end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
           ! --- advance fine solution (third stage of SSPRK) AND sum up the tendecies above into tendSum3rd for interface layer 1 and interface layer 2
           ! -- also restore soln second stage by copying it back from soln new
           block => domain % blocklist
           do while (associated(block))
              call mpas_pool_get_subpool(block % structs, 'tend', tendPool) 
              call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
              call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
              call mpas_pool_get_subpool(block % structs, 'state', statePool)
              call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

              call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

              call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
              call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
              call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

              call mpas_pool_get_array(tendPool, 'u', uTend)
              call mpas_pool_get_array(tendPool, 'h', hTend)
              call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

              call mpas_pool_get_array(tendSum3rdPool, 'u', uTendSum3rd)
              call mpas_pool_get_array(tendSum3rdPool, 'h', hTendSum3rd)
              call mpas_pool_get_array(tendSum3rdPool, 'tracers', tracersTendSum3rd)

              call mpas_pool_get_array(statePool, 'u', uSecondStage, 4) 
              call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
              call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)

              call mpas_pool_get_array(statePool, 'u', uOld, 1)
              call mpas_pool_get_array(statePool, 'h', hOld, 1)
              call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

              call mpas_pool_get_array(statePool, 'u', uNew, 2)
              call mpas_pool_get_array(statePool, 'h', hNew, 2)
              call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

              if(haveFine == 1 .or. haveInterface == 1) then
              if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uNew(:,:) = uOld(:,:)
              end if
              do k =1, nVertLevels
                 do ie = 1, nEdgesInLTSHalo(1,2)
                    iEdge = edgesInLTSHalo(1,2,ie)
                    uSecondStage(k,iEdge) = uNew(k,iEdge)
                    uTendSum3rd(k,iEdge) = uTendSum3rd(k,iEdge) + uTend(k,iEdge)
                 end do
                 do ie = 1, nEdgesInLTSHalo(2,2)
                    iEdge = edgesInLTSHalo(2,2,ie)
                    uTendSum3rd(k,iEdge) = uTendSum3rd(k,iEdge) + uTend(k,iEdge)
                 end do
                 do ie = 1, nEdgesInLTSHalo(1,1)
                    iEdge = edgesInLTSHalo(1,1,ie)
                    uOld(k,iEdge) =  (1.0 / 3.0) * uOld(k,iEdge) + (2.0 / 3.0) * uSecondStage(k,iEdge) + (2.0 / 3.0) *  dtFine * uTend(k,iEdge)
                 end do
                 do ie = 1, nEdgesInLTSHalo(1,3)
                    iEdge = edgesInLTSHalo(1,3,ie)
                    uOld(k,iEdge) =  (1.0 / 3.0) * uOld(k,iEdge) + (2.0 / 3.0) * uSecondStage(k,iEdge) + (2.0 / 3.0) *  dtFine * uTend(k,iEdge)
                 end do
                 do ie = 1, nEdgesInLTSHalo(1,4)
                    iEdge = edgesInLTSHalo(1,4,ie)
                    uOld(k,iEdge) =  (1.0 / 3.0) * uOld(k,iEdge) + (2.0 / 3.0) * uSecondStage(k,iEdge) + (2.0 / 3.0) *  dtFine * uTend(k,iEdge)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,2)
                    iCell = cellsInLTSHalo(1,2,ic)
                    hSecondStage(k,iCell) = hNew(k,iCell)
                    tracersSecondStage(:,k,iCell) = tracersNew(:,k,iCell)
                    hTendSum3rd(k,iCell) = hTendSum3rd(k,iCell) + hTend(k,iCell)
                    tracersTendSum3rd(:,k,iCell) = tracersTendSum3rd(:,k,iCell) + tracersTend(:,k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(2,2)
                    iCell = cellsInLTSHalo(2,2,ic)
                    hTendSum3rd(k,iCell) = hTendSum3rd(k,iCell) + hTend(k,iCell)
                    tracersTendSum3rd(:,k,iCell) = tracersTendSum3rd(:,k,iCell) + tracersTend(:,k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,1)
                    iCell = cellsInLTSHalo(1,1,ic)
                    hNew(k,iCell) = (1.0 / 3.0) * hOld(k,iCell) + (2.0 / 3.0) * hSecondStage(k,iCell) + (2.0 / 3.0) * dtFine * hTend(k,iCell) !soln update for the fine
                    tracersOld(:,k,iCell) = ( (1.0 / 3.0) * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                              (2.0 / 3.0) * tracersSecondStage(:,k,iCell) * hSecondStage(k,iCell)+ & 
                                              (2.0 / 3.0) * dtFine * tracersTend(:,k,iCell) ) / hNew(k,iCell)
                    hOld(k,iCell) = hNew(k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,3)
                    iCell = cellsInLTSHalo(1,3,ic)
                    hNew(k,iCell) = (1.0 / 3.0) * hOld(k,iCell) + (2.0 / 3.0) * hSecondStage(k,iCell) + (2.0 / 3.0) * dtFine * hTend(k,iCell) !soln update for the fine
                    tracersOld(:,k,iCell) = ( (1.0 / 3.0) * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                              (2.0 / 3.0) * tracersSecondStage(:,k,iCell) * hSecondStage(k,iCell) + &
                                              (2.0 / 3.0) * dtFine * tracersTend(:,k,iCell) ) / hNew(k,iCell)
                    hOld(k,iCell) = hNew(k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,4)
                    iCell = cellsInLTSHalo(1,4,ic)
                    hNew(k,iCell) = (1.0 / 3.0) * hOld(k,iCell) + (2.0 / 3.0) * hSecondStage(k,iCell) + (2.0 / 3.0) * dtFine * hTend(k,iCell) !soln update for the fine
                    tracersOld(:,k,iCell) = ( (1.0 / 3.0) * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                              (2.0 / 3.0) * tracersSecondStage(:,k,iCell) * hSecondStage(k,iCell) + &
                                              (2.0 / 3.0) * dtFine * tracersTend(:,k,iCell) ) / hNew(k,iCell)
                    hOld(k,iCell) = hNew(k,iCell)
                 end do
              end do
              if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uOld(:,:) = uNew(:,:)
              end if
              end if

        !BEGIN: TO REMOVE
        !if (block % blockID == 0) then
        !print*, 'LTS - New Soln - Fine'
        !do k =1, nVertLevels
        !   do ie = 1, nEdgesInLTSHalo(2,1)
        !      iEdge = edgesInLTSHalo(2,1,ie)
        !      print*, 'u=', uSecondStage(k,iEdge), 'edge=', iEdge
        !   end do
        !   do ic = 1, nCellsInLTSHalo(1,1)
        !      iCell = cellsInLTSHalo(1,1,ic)
        !      print*, hOld(k,iCell)
        !   end do
        !print*, '--------------------------------'
        !   do ic = 1, nCellsInLTSHalo(1,3)
        !      iCell = cellsInLTSHalo(1,3,ic)
        !      print*, hOld(k,iCell)
        !  end do
        !end do
        !end if
        !END: TO REMOVE
 
              block => block % next
           end do

call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")

call mpas_timer_start("halo updates  LTS")
     call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=1)
     call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=1)
     call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=1)
call mpas_timer_stop("halo updates  LTS")

        end if

     end do  ! end of substepping

     ! ---------------------------------------------------------------------------------------------------------------------

call mpas_timer_start("copy soln into new soln LTS")
        ! --- copy old soln into new soln for the fine cells and edges
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           if(haveFine == 1) then
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(1,1)
                 iEdge = edgesInLTSHalo(1,1,ie)
                 uNew(k,iEdge) =  uOld(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(1,1)
                 iCell = cellsInLTSHalo(1,1,ic)
                 hNew(k,iCell) = hOld(k,iCell)
                 tracersNew(:,k,iCell) = tracersOld(:,k,iCell)
              end do

              if (lts_order == 3) then
                 do ie = 1, nEdgesInLTSHalo(1,3)
                    iEdge = edgesInLTSHalo(1,3,ie)
                    uNew(k,iEdge) =  uOld(k,iEdge)
                 end do
                 do ie = 1, nEdgesInLTSHalo(1,4)
                    iEdge = edgesInLTSHalo(1,4,ie)
                    uNew(k,iEdge) =  uOld(k,iEdge)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,3)
                    iCell = cellsInLTSHalo(1,3,ic)
                    hNew(k,iCell) = hOld(k,iCell)
                    tracersNew(:,k,iCell) = tracersOld(:,k,iCell)
                 end do
                 do ic = 1, nCellsInLTSHalo(1,4)
                    iCell = cellsInLTSHalo(1,4,ic)
                    hNew(k,iCell) = hOld(k,iCell)
                    tracersNew(:,k,iCell) = tracersOld(:,k,iCell)
                 end do
              end if

           end do
           end if

           block => block % next
        end do
call mpas_timer_stop("copy soln into new soln LTS")

     ! --- copy secondStage soln into newSoln for coarse

     if (lts_order == 2) then

call mpas_timer_start("copy soln into new soln LTS")
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
           call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
           call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)

           if(haveCoarse == 1) then
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(2,1)
                 iEdge = edgesInLTSHalo(2,1,ie)
                 uNew(k,iEdge) =  uSecondStage(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(2,1)
                 iCell = cellsInLTSHalo(2,1,ic)
                 hNew(k,iCell) = hSecondStage(k,iCell)
                 tracersNew(:,k,iCell) = tracersSecondStage(:,k,iCell) 
              end do
           end do
           end if
           block => block % next
        end do

call mpas_timer_stop("copy soln into new soln LTS")
     end if

     ! --- compute third stage of SSPRK for coarse if needed

     if (lts_order == 3) then

call mpas_timer_start("halo updates  LTS")
        ! --- update halos for diagnostic variables

        if (config_h_mom_eddy_visc4 > 0.0) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence', timeLevel=2)
           call mpas_dmpar_field_halo_exch(domain, 'vorticity', timeLevel=2)
        end if
call mpas_timer_stop("halo updates  LTS")

call mpas_timer_start("tendencies LTS")
call cpu_time(startTendencies)
        ! --- compute tendencies for coarse
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

           call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

           if(haveCoarse == 1) then
           call sw_compute_solve_diagnostics(dt, statePool, meshPool, 4) 
           call sw_compute_tend(tendPool, statePool, meshPool, 4) 
           call sw_compute_scalar_tend(tendPool, statePool, meshPool, 4) 
           call sw_enforce_boundary_edge(tendPool, meshPool)
           end if
           block => block % next
        end do
call cpu_time(finishTendencies)
totalTimeTendencies = totalTimeTendencies + (finishTendencies - startTendencies)
call mpas_timer_stop("tendencies LTS")

call mpas_timer_start("advance soln LTS")
call cpu_time(startAdvance)
        ! --- advance solution (3rd stage of SSPRK for the coarse using dt)
        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

           call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
  
           call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
           call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

           call mpas_pool_get_array(statePool, 'u', uSecondStage, 4)
           call mpas_pool_get_array(statePool, 'h', hSecondStage, 4)
           call mpas_pool_get_array(statePool, 'tracers', tracersSecondStage, 4)

           call mpas_pool_get_array(statePool, 'u', uOld, 1)
           call mpas_pool_get_array(statePool, 'h', hOld, 1)
           call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

           call mpas_pool_get_array(statePool, 'u', uNew, 2)
           call mpas_pool_get_array(statePool, 'h', hNew, 2)
           call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

           call mpas_pool_get_array(tendPool, 'u', uTend)
           call mpas_pool_get_array(tendPool, 'h', hTend)
           call mpas_pool_get_array(tendPool, 'tracers', tracersTend)

           if(haveCoarse == 1) then
           do k =1, nVertLevels
              do ie = 1, nEdgesInLTSHalo(2,1)
                 iEdge = edgesInLTSHalo(2,1,ie)
                 uNew(k,iEdge) = (1.0 / 3.0) * uOld(k,iEdge) + (2.0 / 3.0) * uSecondStage(k,iEdge) + (2.0 / 3.0) * dt * uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(2,3)
                 iEdge = edgesInLTSHalo(2,3,ie)
                 uNew(k,iEdge) = (1.0 / 3.0) * uOld(k,iEdge) + (2.0 / 3.0) * uSecondStage(k,iEdge) + (2.0 / 3.0) * dt * uTend(k,iEdge)
              end do
              do ie = 1, nEdgesInLTSHalo(2,4)
                 iEdge = edgesInLTSHalo(2,4,ie)
                 uNew(k,iEdge) = (1.0 / 3.0) * uOld(k,iEdge) + (2.0 / 3.0) * uSecondStage(k,iEdge) + (2.0 / 3.0) * dt * uTend(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(2,1)
                 iCell = cellsInLTSHalo(2,1,ic)
                 hNew(k,iCell) = (1.0 / 3.0) * hOld(k,iCell) + (2.0 / 3.0) * hSecondStage(k,iCell) + (2.0 / 3.0) * dt * hTend(k,iCell)
                 tracersNew(:,k,iCell) = ( (1.0 / 3.0) * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                           (2.0 / 3.0) * tracersSecondStage(:,k,iCell) * hSecondStage(k,iCell) + &
                                           (2.0 / 3.0) * dt * tracersTend(:,k,iCell) ) / hNew(k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(2,3)
                 iCell = cellsInLTSHalo(2,3,ic)
                 hNew(k,iCell) = (1.0 / 3.0) * hOld(k,iCell) + (2.0 / 3.0) * hSecondStage(k,iCell) + (2.0 / 3.0) * dt * hTend(k,iCell)
                 tracersNew(:,k,iCell) = ( (1.0 / 3.0) * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                           (2.0 / 3.0) * tracersSecondStage(:,k,iCell) * hSecondStage(k,iCell) + &
                                           (2.0 / 3.0) * dt * tracersTend(:,k,iCell) ) / hNew(k,iCell)
              end do
              do ic = 1, nCellsInLTSHalo(2,4)
                 iCell = cellsInLTSHalo(2,4,ic)
                 hNew(k,iCell) = (1.0 / 3.0) * hOld(k,iCell) + (2.0 / 3.0) * hSecondStage(k,iCell) + (2.0 / 3.0) * dt * hTend(k,iCell)
                 tracersNew(:,k,iCell) = ( (1.0 / 3.0) * tracersOld(:,k,iCell) * hOld(k,iCell) + & 
                                           (2.0 / 3.0) * tracersSecondStage(:,k,iCell) * hSecondStage(k,iCell) + &
                                           (2.0 / 3.0) * dt * tracersTend(:,k,iCell) ) / hNew(k,iCell)
              end do
           end do
           if (config_test_case == 1) then    ! For case 1, wind field should be fixed
                uNew(:,:) = uOld(:,:)
           end if
           end if

        !BEGIN: TO REMOVE
        !if (block % blockID == 1) then
        !print*, 'LTS - New Soln - Coarse'
        !do k =1, nVertLevels
        !   do ie = 1, nEdgesInLTSHalo(2,1)
        !      iEdge = edgesInLTSHalo(2,1,ie)
        !      print*, 'u=', uSecondStage(k,iEdge), 'edge=', iEdge
        !   end do
        !   do ic = 1, nCellsInLTSHalo(2,1)
        !      iCell = cellsInLTSHalo(2,1,ic)
        !      print*, hNew(k,iCell)
        !   end do
        !   do ic = 1, nCellsInLTSHalo(2,3)
        !      iCell = cellsInLTSHalo(2,3,ic)
        !      print*, hNew(k,iCell)
        !  end do
        !end do
        !end if
        !END: TO REMOVE

           block => block % next
        end do
call cpu_time(finishAdvance)
totalTimeAdvance = totalTimeAdvance + (finishAdvance - startAdvance)
call mpas_timer_stop("advance soln LTS")
     end if

     ! --- INTERFACE CORRECTION STEP: correct interface layer 1 and interface layer 2 using tendSum2nd and tendSum1st (and tendSum3rd)  

call mpas_timer_start("interface correction LTS")
call cpu_time(startInterface)
     ! --- correct solution for interface layer 1 and interface layer 2
     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
        call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
        call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
        call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
        call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

        call mpas_pool_get_array(tendSum3rdPool, 'u', uTendSum3rd)
        call mpas_pool_get_array(tendSum3rdPool, 'h', hTendSum3rd)
        call mpas_pool_get_array(tendSum3rdPool, 'tracers', tracersTendSum3rd)

        call mpas_pool_get_array(tendSum2ndPool, 'u', uTendSum2nd)
        call mpas_pool_get_array(tendSum2ndPool, 'h', hTendSum2nd)
        call mpas_pool_get_array(tendSum2ndPool, 'tracers', tracersTendSum2nd)
 
        call mpas_pool_get_array(tendSum1stPool, 'u', uTendSum1st)
        call mpas_pool_get_array(tendSum1stPool, 'h', hTendSum1st)
        call mpas_pool_get_array(tendSum1stPool, 'tracers', tracersTendSum1st)
          
        call mpas_pool_get_array(statePool, 'u', uOld, 1)
        call mpas_pool_get_array(statePool, 'h', hOld, 1)
        call mpas_pool_get_array(statePool, 'tracers', tracersOld, 1)

        call mpas_pool_get_array(statePool, 'u', uNew, 2)
        call mpas_pool_get_array(statePool, 'h', hNew, 2)
        call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

        if (haveInterface == 1) then
        do k =1, nVertLevels
           do iRegion = 1,nRegions
              do ie = 1, nEdgesInLTSHalo(iRegion,2)
                 iEdge = edgesInLTSHalo(iRegion,2,ie)
                 uNew(k,iEdge) =  uOld(k,iEdge) + weightTendSum2nd * dtFine * uTendSum2nd(k,iEdge) + weightTendSum1st * dtFine * uTendSum1st(k,iEdge) + &
                                  weightTendSum3rd * dtFine * uTendSum3rd(k,iEdge)
              end do
              do ic = 1, nCellsInLTSHalo(iRegion,2)
                 iCell = cellsInLTSHalo(iRegion,2,ic)
                 hNew(k,iCell) = hOld(k,iCell) + weightTendSum2nd * dtFine * hTendSum2nd(k,iCell) + weightTendSum1st * dtFine * hTendSum1st(k,iCell) + &
                                 weightTendsum3rd * dtFine * hTendSum3rd(k,iCell)
                 tracersNew(:,k,iCell) = ( tracersOld(:,k,iCell) * hOld(k,iCell) + &
                                           weightTendSum2nd * dtFine * tracersTendSum2nd(:,k,iCell) + &
                                           weightTendSum1st * dtFine * tracersTendSum1st(:,k,iCell) + &
                                           weightTendSum3rd * dtFine * tracersTendSum3rd(:,k,iCell) ) / hNew(k,iCell)
              end do
           end do
        end do
        if (config_test_case == 1) then    ! For case 1, wind field should be fixed
           uNew(:,:) = uOld(:,:)
        end if
        end if

        !BEGIN: TO REMOVE
        !if (block % blockID == 2) then
        !print*, 'LTS - New Soln - Interface Layers'
        !do k =1, nVertLevels
        !  do ie = 1, nEdgesInLTSHalo(2,2)
        !     iEdge = edgesInLTSHalo(2,2,ie)
        !     print*, 'u=',  uNew(k,iEdge), 'edge=', iEdge
        !  end do
        !  do ic = 1, nCellsInLTSHalo(1,2)
        !     iCell = cellsInLTSHalo(1,2,ic)
        !     print*, hNew(k,iCell)
        !  end do
        !  print*,'----------------------------------'
        !  do ic = 1, nCellsInLTSHalo(2,2)
        !     iCell = cellsInLTSHalo(2,2,ic)
        !     print*, hNew(k,iCell)
        !  end do
        !end do
        !end if
        !END: TO REMOVE

        block => block % next
     end do
call cpu_time(finishInterface)
totalTimeInterface = totalTimeInterface + (finishInterface - startInterface)
call mpas_timer_stop("interface correction LTS")

call mpas_timer_start("barrier time")
call cpu_time(startBarrier)
call MPI_Barrier(domain % dminfo % comm, err)
call cpu_time(finishBarrier)
totalTimeBarrier = totalTimeBarrier + (finishBarrier - startBarrier)
call mpas_timer_stop("barrier time")
call mpas_timer_start("halo updates  LTS")
        ! --- update halos for prognostic variables
        call mpas_dmpar_field_halo_exch(domain, 'u', timeLevel=2)
        call mpas_dmpar_field_halo_exch(domain, 'h', timeLevel=2)
        call mpas_dmpar_field_halo_exch(domain, 'tracers', timeLevel=2)
call mpas_timer_stop("halo updates  LTS")

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! END LTS SCHEME
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mpas_timer_start("clean-up time LTS")
     !
     !  A little clean up at the end: compute diagnostics for new state
     !

     block => domain % blocklist
     do while (associated(block))
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(statePool, 'u', uOld, 1)
        call mpas_pool_get_array(statePool, 'u', uNew, 2)
        call mpas_pool_get_array(statePool, 'h', hNew, 2)
        call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

        call mpas_pool_get_array(statePool, 'uReconstructX', uReconstructX, 2)
        call mpas_pool_get_array(statePool, 'uReconstructY', uReconstructY, 2)
        call mpas_pool_get_array(statePool, 'uReconstructZ', uReconstructZ, 2)
        call mpas_pool_get_array(statePool, 'uReconstructZonal', uReconstructZonal, 2)
        call mpas_pool_get_array(statePool, 'uReconstructMeridional', uReconstructMeridional, 2)
        
        if (config_test_case == 1) then    ! For case 1, wind field should be fixed
           uNew(:,:) = uOld(:,:)
        end if

        !call sw_compute_solve_diagnostics(dt, statePool, meshPool, 2)

call mpas_timer_start("barrier time")
call cpu_time(startBarrier)
call MPI_Barrier(domain % dminfo % comm, err)
call cpu_time(finishBarrier)
totalTimeBarrier = totalTimeBarrier + (finishBarrier - startBarrier)
call mpas_timer_stop("barrier time")
call mpas_timer_start("last halo update  LTS")
        call mpas_dmpar_field_halo_exch(domain, 'pv_edge', timeLevel=1)
call mpas_timer_stop("last halo update  LTS")

        call mpas_reconstruct(meshPool, uNew,          &
                          uReconstructX, uReconstructY, uReconstructZ, &
                          uReconstructZonal, uReconstructMeridional )

        block => block % next
     end do

     block => domain % blocklist
     do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
        call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
        call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)

        call mpas_pool_destroy_pool(tendSum1stPool)
        call mpas_pool_destroy_pool(tendSum2ndPool)
        call mpas_pool_destroy_pool(tendSum3rdPool)

        call mpas_pool_remove_subpool(block % structs, 'tend_sum_1st')
        call mpas_pool_remove_subpool(block % structs, 'tend_sum_2nd')
        call mpas_pool_remove_subpool(block % structs, 'tend_sum_3rd')

        block => block % next
     end do
call mpas_timer_stop("clean-up time LTS")

   end subroutine sw_lts

   subroutine sw_compute_tend(tendPool, statePool, meshPool, timeLevelIn)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed tendencies for prognostic variables
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in), optional :: timeLevelIn

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, timeLevel
      real (kind=RKIND) :: flux, vorticity_abs, workpv, q, upstream_bias

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle, &
                                                  meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, divergence, h_vertex
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell
      real (kind=RKIND) :: r, u_diffusion

      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity

      real (kind=RKIND), dimension(:,:), pointer :: u_src
      real (kind=RKIND), parameter :: rho_ref = 1000.0
      real (kind=RKIND) :: ke_edge

      logical, pointer :: config_wind_stress, config_bottom_drag
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc2, config_h_mom_eddy_visc4

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(swConfigs, 'config_bottom_drag', config_bottom_drag)
      call mpas_pool_get_config(swConfigs, 'config_wind_stress', config_wind_stress)
      call mpas_pool_get_config(swConfigs, 'config_h_mom_eddy_visc2', config_h_mom_eddy_visc2)
      call mpas_pool_get_config(swConfigs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)

      call mpas_pool_get_array(statePool, 'h', h, timeLevel)
      call mpas_pool_get_array(statePool, 'u', u, timeLevel)
      call mpas_pool_get_array(statePool, 'v', v, timeLevel)
      call mpas_pool_get_array(statePool, 'h_edge', h_edge, timeLevel)
      call mpas_pool_get_array(statePool, 'circulation', circulation, timeLevel)
      call mpas_pool_get_array(statePool, 'vorticity', vorticity, timeLevel)
      call mpas_pool_get_array(statePool, 'divergence', divergence, timeLevel)
      call mpas_pool_get_array(statePool, 'ke', ke, timeLevel)
      call mpas_pool_get_array(statePool, 'pv_edge', pv_edge, timeLevel)
      call mpas_pool_get_array(statePool, 'vh', vh, timeLevel)

      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'u_src', u_src)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      call mpas_pool_get_array(tendPool, 'h', tend_h)
      call mpas_pool_get_array(tendPool, 'u', tend_u)
                  
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      !
      ! Compute height tendency for each cell
      !

      ! THIS IS CORRECT BUT DOES NOT PROVIDE BFB MATCHING AMONG DIFFERENT MPI PARTITIONS
      ! tend_h(:,:) = 0.0
      ! do iEdge = 1, nEdges
      !    cell1 = cellsOnEdge(1,iEdge)
      !    cell2 = cellsOnEdge(2,iEdge)
      !    do k = 1, nVertLevels
      !       flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
      !       tend_h(k,cell1) = tend_h(k,cell1) - flux
      !       tend_h(k,cell2) = tend_h(k,cell2) + flux
      !    end do
      ! end do 
      ! do iCell = 1, nCellsSolve
      !    do k = 1, nVertLevels
      !       tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
      !    end do
      ! end do

      tend_h(:,:) = 0.0
      do iCell = 1, nCells
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, nVertLevels
               flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
               tend_h(k, iCell) = tend_h(k, iCell) + edgeSignOnCell(i, iCell) * flux
            end do
         end do
         do k = 1, nVertLevels
            tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
         end do
      end do

      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !
      tend_u(:,:) = 0.0
      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)
         
         do k = 1, nVertLevels
            q = 0.0
            do j = 1, nEdgesOnEdge(iEdge)
               eoe = edgesOnEdge(j,iEdge)
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
               !workpv = pv_edge(k,iEdge) 
               q = q + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * h_edge(k,eoe) 
            end do

            tend_u(k,iEdge) =       &
                              q     &
                              - (   ke(k,cell2) - ke(k,cell1) + &
                                    gravity * (h(k,cell2) + h_s(cell2) - h(k,cell1) - h_s(cell1)) &
                                  ) / dcEdge(iEdge)
         end do
      end do


     ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
     !                    only valid for visc == constant
     if (config_h_mom_eddy_visc2 > 0.0) then
        do iEdge = 1, nEdgesSolve
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)

           do k = 1, nVertLevels
              u_diffusion =   ( divergence(k,cell2)  -  divergence(k,cell1) ) / dcEdge(iEdge) &
                   -(vorticity(k,vertex2)  - vorticity(k,vertex1) ) / dvEdge(iEdge)
              u_diffusion = meshScalingDel2(iEdge) * config_h_mom_eddy_visc2 * u_diffusion
              tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion
           end do
        end do
     end if

     !
     ! velocity tendency: del4 dissipation, -\nu_4 \nabla^4 u
     !   computed as \nabla^2 u = \nabla divergence + k \times \nabla vorticity
     !   applied recursively.
     !   strictly only valid for h_mom_eddy_visc4 == constant
     !
     if (config_h_mom_eddy_visc4 > 0.0) then
        allocate(delsq_divergence(nVertLevels, nCells+1))
        allocate(delsq_u(nVertLevels, nEdges+1))
        allocate(delsq_circulation(nVertLevels, nVertices+1))
        allocate(delsq_vorticity(nVertLevels, nVertices+1))

        delsq_u(:,:) = 0.0

        ! Compute \nabla^2 u = \nabla divergence + k \times \nabla vorticity
        do iEdge = 1, nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)

           do k = 1, nVertLevels

              delsq_u(k,iEdge) = ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                   -( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)

           end do
        end do

        ! vorticity using \nabla^2 u
        delsq_circulation(:,:) = 0.0
        do iEdge = 1, nEdges
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)
           do k=1,nVertLevels
              delsq_circulation(k,vertex1) = delsq_circulation(k,vertex1) &
                   - dcEdge(iEdge) * delsq_u(k,iEdge)
              delsq_circulation(k,vertex2) = delsq_circulation(k,vertex2) &
                   + dcEdge(iEdge) * delsq_u(k,iEdge)
           end do
        end do
        do iVertex = 1, nVertices
           r = 1.0 / areaTriangle(iVertex)
           do k = 1, nVertLevels
              delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
           end do
        end do

        ! Divergence using \nabla^2 u
        delsq_divergence(:,:) = 0.0
        do iEdge = 1, nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           do k=1,nVertLevels
              delsq_divergence(k,cell1) = delsq_divergence(k,cell1) &
                   + delsq_u(k,iEdge)*dvEdge(iEdge)
              delsq_divergence(k,cell2) = delsq_divergence(k,cell2) &
                   - delsq_u(k,iEdge)*dvEdge(iEdge)
           end do
        end do
        do iCell = 1,nCells
           r = 1.0 / areaCell(iCell)
           do k = 1, nVertLevels
              delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
           end do
        end do

        ! Compute - \kappa \nabla^4 u 
        ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
        do iEdge = 1, nEdgesSolve
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           vertex1 = verticesOnEdge(1,iEdge)
           vertex2 = verticesOnEdge(2,iEdge)

           do k = 1, nVertLevels

              u_diffusion = (  delsq_divergence(k,cell2) &
                   - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                   -(  delsq_vorticity(k,vertex2) &
                   - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)

              u_diffusion = meshScalingDel4(iEdge) * config_h_mom_eddy_visc4 * u_diffusion
              tend_u(k,iEdge) = tend_u(k,iEdge) - u_diffusion

           end do
        end do

        deallocate(delsq_divergence)
        deallocate(delsq_u)
        deallocate(delsq_circulation)
        deallocate(delsq_vorticity)

     end if

     ! Compute u (velocity) tendency from wind stress (u_src)
     if(config_wind_stress) then
         do iEdge = 1, nEdges
            tend_u(1,iEdge) =  tend_u(1,iEdge) &
                  + u_src(1,iEdge) / rho_ref / h_edge(1,iEdge)
         end do
     endif

     if (config_bottom_drag) then
         do iEdge = 1,  nEdges
             ! bottom drag is the same as POP:
             ! -c |u| u  where c is unitless and 1.0e-3.
             ! see POP Reference guide, section 3.4.4.
             ke_edge = 0.5 * ( ke(1,cellsOnEdge(1,iEdge)) &
                   + ke(1,cellsOnEdge(2,iEdge)))

             tend_u(1,iEdge) = tend_u(1,iEdge)  &
                  - 1.0e-3*u(1,iEdge) &
                  *sqrt(2.0*ke_edge)/h_edge(1,iEdge)
         end do
     endif
 
   end subroutine sw_compute_tend


   subroutine sw_compute_scalar_tend(tendPool, statePool, meshPool, timeLevelIn)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in), optional :: timeLevelIn

      integer :: iCell, iEdge, k, iTracer, cell1, cell2, i
      real (kind=RKIND) :: flux, tracer_edge, r
      real (kind=RKIND) :: invAreaCell1, invAreaCell2, tracer_turb_flux
      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), allocatable :: boundaryMask
      real (kind=RKIND), dimension(:,:,:), allocatable:: delsq_tracer
      
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers, tracer_tend
      integer, dimension(:,:), pointer :: cellsOnEdge, boundaryCell
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order
      real (kind=RKIND), dimension(:,:), pointer :: u, h_edge

      integer, pointer :: config_tracer_adv_order
      logical, pointer :: config_monotonic

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nTracers
      integer :: timeLevel

      real (kind=RKIND), pointer :: config_h_tracer_eddy_diff2, config_h_tracer_eddy_diff4

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      endif

      call mpas_pool_get_config(swConfigs, 'config_tracer_adv_order', config_tracer_adv_order)
      call mpas_pool_get_config(swConfigs, 'config_monotonic', config_monotonic)
      call mpas_pool_get_config(swConfigs, 'config_h_tracer_eddy_diff2', config_h_tracer_eddy_diff2)
      call mpas_pool_get_config(swConfigs, 'config_h_tracer_eddy_diff4', config_h_tracer_eddy_diff4)

      call mpas_pool_get_array(statePool, 'u', u, timeLevel)
      call mpas_pool_get_array(statePool, 'h_edge', h_edge, timeLevel)
      call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'deriv_two', deriv_two)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

      call mpas_pool_get_array(tendPool, 'tracers', tracer_tend)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
      call mpas_pool_get_dimension(meshPool, 'nTracers', nTracers)

      coef_3rd_order = 0.
      if (config_tracer_adv_order == 3) coef_3rd_order = 1.0
      if (config_tracer_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      tracer_tend(:,:,:) = 0.0

      if (config_tracer_adv_order == 2) then

      !TODO: change this as for the tendency of h in order to have bfb matching
      !between MPI partitons. Right now it is correct but it does not provide
      !bfb matching.
      do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCells .and. cell2 <= nCells) then
               do k = 1, nVertLevels
                  do iTracer = 1, nTracers
                     tracer_edge = 0.5 * (tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))
                     flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) * tracer_edge
                     tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux / areaCell(cell1)
                     tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux / areaCell(cell2)
                  end do 
               end do 
            end if
      end do 

      else if (config_tracer_adv_order == 3) then

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= nCells .and. cell2 <= nCells) then

               do k = 1, nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  do iTracer = 1, nTracers
 
                     !-- if not a boundary cell
                     if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then

                        d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1)
                        d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2)

                        !-- all edges of cell 1
                        do i = 1, nEdgesOnCell(cell1)
                                d2fdx2_cell1 = d2fdx2_cell1 + &
                                deriv_two(i+1,1,iEdge) * tracers(iTracer,k, CellsOnCell(i,cell1))
                        end do

                        !-- all edges of cell 2
                        do i = 1, nEdgesOnCell(cell2)
                                d2fdx2_cell2 = d2fdx2_cell2 + &
                                deriv_two(i+1,2,iEdge) * tracers(iTracer,k,CellsOnCell(i,cell2))
                        end do

                     endif

                     !-- if u > 0:
                     if (u(k,iEdge) > 0) then
                        flux = dvEdge(iEdge) * u(k,iEdge) * h_edge(k,iEdge) * (          &
                             0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                             -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                     !-- else u <= 0:
                     else
                        flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                             0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                             +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
                     end if

                     !-- update tendency
                     tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux / areaCell(cell1)
                     tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux / areaCell(cell2)
                  enddo
               end do
            end if
         end do

      else  if (config_tracer_adv_order == 4) then

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if an edge is not on the outer-most ring of the halo
            if (cell1 <= nCells .and. cell2 <= nCells) then

               do k = 1, nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  do iTracer = 1, nTracers

                     !-- if not a boundary cell
                     if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then

                        d2fdx2_cell1 = deriv_two(1,1,iEdge) * tracers(iTracer,k,cell1)
                        d2fdx2_cell2 = deriv_two(1,2,iEdge) * tracers(iTracer,k,cell2)

                        !-- all edges of cell 1
                        do i = 1, nEdgesOnCell(cell1)
                                d2fdx2_cell1 = d2fdx2_cell1 + &
                                deriv_two(i+1,1,iEdge) * tracers(iTracer,k, cellsOnCell(i,cell1))
                        end do

                        !-- all edges of cell 2
                        do i = 1, nEdgesOnCell(cell2)
                                d2fdx2_cell2 = d2fdx2_cell2 + &
                                deriv_two(i+1,2,iEdge) * tracers(iTracer,k, cellsOnCell(i,cell2))
                        end do

                     endif

                     flux = dvEdge(iEdge) *  u(k,iEdge) * h_edge(k,iEdge) * (          &
                          0.5*(tracers(iTracer,k,cell1) + tracers(iTracer,k,cell2))      &
                             -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                     !-- update tendency
                     tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux / areaCell(cell1)
                     tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux / areaCell(cell2)
                  enddo
               end do
            end if
         end do

      endif   ! if (config_tracer_adv_order == 2 )

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      if ( config_h_tracer_eddy_diff2 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(nVertLevels, nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invAreaCell1 = 1.0/areaCell(cell1)
            invAreaCell2 = 1.0/areaCell(cell2)

            do k = 1, nVertLevels
              do iTracer = 1,  nTracers
                 ! \kappa_2 \nabla \phi on edge
                 tracer_turb_flux = config_h_tracer_eddy_diff2 &
                    *( tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge)

                 ! div(h \kappa_2 \nabla \phi) at cell center
                 flux = dvEdge(iEdge) * h_edge(k,iEdge) * tracer_turb_flux * boundaryMask(k, iEdge)
                 tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) + flux * invAreaCell1
                 tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) - flux * invAreaCell2
              end do
            end do

         end do

        deallocate(boundaryMask)

      end if

      !
      ! tracer tendency: del4 horizontal tracer diffusion, &
      !    div(h \kappa_4 \nabla [div(h \nabla \phi)])
      !
      if ( config_h_tracer_eddy_diff4 > 0.0 ) then

         !
         ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
         !
         allocate(boundaryMask(nVertLevels, nEdges+1))
         boundaryMask = 1.0
         where(boundaryEdge.eq.1) boundaryMask=0.0

         allocate(delsq_tracer(nTracers, nVertLevels, nCells+1))

         delsq_tracer(:,:,:) = 0.

         ! first del2: div(h \nabla \phi) at cell center
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k = 1, nVertLevels
              do iTracer = 1, nTracers
                 delsq_tracer(iTracer,k,cell1) = delsq_tracer(iTracer,k,cell1) &
                    + dvEdge(iEdge) * h_edge(k,iEdge) * (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge) * boundaryMask(k,iEdge)
                 delsq_tracer(iTracer,k,cell2) = delsq_tracer(iTracer,k,cell2) &
                    - dvEdge(iEdge) * h_edge(k,iEdge) * (tracers(iTracer,k,cell2) - tracers(iTracer,k,cell1)) / dcEdge(iEdge) * boundaryMask(k,iEdge)
              end do
            end do

         end do

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k = 1, nVertLevels
            do iTracer = 1, nTracers
               delsq_tracer(iTracer,k,iCell) = delsq_tracer(iTracer,k,iCell) * r
            end do
            end do
         end do

         ! second del2: div(h \nabla [delsq_tracer]) at cell center
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invAreaCell1 = 1.0 / areaCell(cell1)
            invAreaCell2 = 1.0 / areaCell(cell2)

            do k = 1, nVertLevels
            do iTracer = 1, nTracers
               tracer_turb_flux = config_h_tracer_eddy_diff4 * (delsq_tracer(iTracer,k,cell2) - delsq_tracer(iTracer,k,cell1)) / dcEdge(iEdge)
               flux = dvEdge(iEdge) * tracer_turb_flux
               tracer_tend(iTracer,k,cell1) = tracer_tend(iTracer,k,cell1) - flux * invAreaCell1 * boundaryMask(k,iEdge)
               tracer_tend(iTracer,k,cell2) = tracer_tend(iTracer,k,cell2) + flux * invAreaCell2 * boundaryMask(k,iEdge)
            end do
            enddo

         end do

         deallocate(delsq_tracer)
         deallocate(boundaryMask)

      end if

   end subroutine sw_compute_scalar_tend


   subroutine sw_compute_solve_diagnostics(dt, statePool, meshPool, timeLevelIn)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (mpas_pool_type), intent(inout) :: statePool
      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in), optional :: timeLevelIn

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
      real (kind=RKIND) :: flux, vorticity_abs, workpv, invAreaTri, r_tmp

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree
      real (kind=RKIND), dimension(:), pointer :: h_s, fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, tend_h, tend_u, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, divergence, &
                                                    h_vertex, vorticity_cell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex, boundaryEdge, boundaryCell, & 
                                          cellsOnCell, edgeSignOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND) :: r, h1, h2
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND) :: coef_3rd_order

      integer :: timeLevel

      logical, pointer :: config_monotonic
      integer, pointer :: config_thickness_adv_order
      real (kind=RKIND), pointer :: config_apvm_upwinding

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(swConfigs, 'config_monotonic', config_monotonic)
      call mpas_pool_get_config(swConfigs, 'config_thickness_adv_order', config_thickness_adv_order)
      call mpas_pool_get_config(swConfigs, 'config_apvm_upwinding', config_apvm_upwinding)

      call mpas_pool_get_array(statePool, 'h', h, timeLevel)
      call mpas_pool_get_array(statePool, 'u', u, timeLevel)
      call mpas_pool_get_array(statePool, 'v', v, timeLevel)
      call mpas_pool_get_array(statePool, 'vh', vh, timeLevel)
      call mpas_pool_get_array(statePool, 'h_edge', h_edge, timeLevel)
      call mpas_pool_get_array(statePool, 'h_vertex', h_vertex, timeLevel)
      call mpas_pool_get_array(statePool, 'h', tend_h, timeLevel)
      call mpas_pool_get_array(statePool, 'u', tend_u, timeLevel)
      call mpas_pool_get_array(statePool, 'circulation', circulation, timeLevel)
      call mpas_pool_get_array(statePool, 'vorticity', vorticity, timeLevel)
      call mpas_pool_get_array(statePool, 'divergence', divergence, timeLevel)
      call mpas_pool_get_array(statePool, 'ke', ke, timeLevel)
      call mpas_pool_get_array(statePool, 'pv_edge', pv_edge, timeLevel)
      call mpas_pool_get_array(statePool, 'pv_vertex', pv_vertex, timeLevel)
      call mpas_pool_get_array(statePool, 'pv_cell', pv_cell, timeLevel)
      call mpas_pool_get_array(statePool, 'vorticity_cell', vorticity_cell, timeLevel)
      call mpas_pool_get_array(statePool, 'gradPVn', gradPVn, timeLevel)
      call mpas_pool_get_array(statePool, 'gradPVt', gradPVt, timeLevel)

      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'deriv_two', deriv_two)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
                  
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)

      !
      ! Find those cells that have an edge on the boundary
      !
      boundaryCell(:,:) = 0
      do iEdge = 1, nEdges
       do k = 1, nVertLevels
         if(boundaryEdge(k,iEdge).eq.1) then
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           boundaryCell(k,cell1) = 1
           boundaryCell(k,cell2) = 1
         endif
       enddo
      enddo

      !
      ! Compute height on cell edges at velocity locations
      !   Namelist options control the order of accuracy of the reconstructed h_edge value
      !

      coef_3rd_order = 0.
      if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
      if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

      if (config_thickness_adv_order == 2) then

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCells .and. cell2 <= nCells) then
               do k = 1, nVertLevels
                  h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
               end do 
            end if
         end do 

      else if (config_thickness_adv_order == 3) then

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= nCells .and. cell2 <= nCells) then

               do k = 1, nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  !-- if not a boundary cell
                  if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then

                     d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                     d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                     !-- all edges of cell 1
                     do i = 1, nEdgesOnCell(cell1)
                             d2fdx2_cell1 = d2fdx2_cell1 + &
                             deriv_two(i+1,1,iEdge) * h(k, cellsOnCell(i,cell1))
                     end do

                     !-- all edges of cell 2
                     do i = 1, nEdgesOnCell(cell2)
                             d2fdx2_cell2 = d2fdx2_cell2 + &
                             deriv_two(i+1,2,iEdge) * h(k, cellsOnCell(i,cell2))
                     end do

                  endif

                  !-- if u > 0:
                  if (u(k,iEdge) > 0) then
                     h_edge(k,iEdge) =     &
                          0.5*(h(k,cell1) + h(k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                          -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  !-- else u <= 0:
                  else
                     h_edge(k,iEdge) =     &
                          0.5*(h(k,cell1) + h(k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                          +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
                  end if

               end do   ! do k
            end if      ! if (cell1 <=
         end do         ! do iEdge

      else  if (config_thickness_adv_order == 4) then

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            !-- if a cell not on the most outside ring of the halo
            if (cell1 <= nCells .and. cell2 <= nCells) then

               do k = 1, nVertLevels

                  d2fdx2_cell1 = 0.0
                  d2fdx2_cell2 = 0.0

                  !-- if not a boundary cell
                  if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then

                     d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                     d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                     !-- all edges of cell 1
                     do i = 1, nEdgesOnCell(cell1)
                             d2fdx2_cell1 = d2fdx2_cell1 + &
                             deriv_two(i+1,1,iEdge) * h(k, cellsOnCell(i,cell1))
                     end do

                     !-- all edges of cell 2
                     do i = 1, nEdgesOnCell(cell2)
                             d2fdx2_cell2 = d2fdx2_cell2 + &
                             deriv_two(i+1,2,iEdge) * h(k, cellsOnCell(i,cell2))
                     end do

                  endif

                  h_edge(k,iEdge) =   &
                       0.5*(h(k,cell1) + h(k,cell2))      &
                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               end do   ! do k
            end if      ! if (cell1 <=
         end do         ! do iEdge

      endif   ! if(config_thickness_adv_order == 2)

      !
      ! set the velocity in the nEdges+1 slot to zero, this is a dummy address
      !    used to when reading for edges that do not exist
      !
      u(:,nEdges+1) = 0.0

      !
      ! Compute circulation and relative vorticity at each vertex
      !

      do iVertex = 1, nVertices
        ! circulation(:, iVertex) = 0.0_RKIND
         vorticity(:, iVertex) = 0.0_RKIND
         invAreaTri = 1.0_RKIND / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, nVertLevels
              r_tmp = dcEdge(iEdge) * u(k, iEdge)
              !circulation(k, iVertex) = circulation(k, iVertex) +
              !edgeSignOnVertex(i, iVertex) * r_tmp
             vorticity(k, iVertex) = vorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri
            end do
         end do
      end do

      ! THIS IS CORRECT BUT DOES NOT PROVIDE BFB MATCHING AMONG DIFFERENT MPI PARTITIONS
      !circulation(:,:) = 0.0
      !do iEdge = 1, nEdges
      !   do k = 1, nVertLevels
      !      circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
      !      circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
      !   end do
      !end do
      !do iVertex = 1, nVertices
      !   do k = 1, nVertLevels
      !      vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
      !   end do
      !end do


      !
      ! Compute the divergence at each cell center
      !
      divergence(:,:) = 0.0
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         if (cell1 <= nCells) then
            do k = 1, nVertLevels
              divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
            enddo
         endif
         if(cell2 <= nCells) then
            do k = 1, nVertLevels
              divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
            enddo
         end if
      end do
      do iCell = 1, nCells
        r = 1.0 / areaCell(iCell)
        do k = 1, nVertLevels
           divergence(k,iCell) = divergence(k,iCell) * r
        enddo
      enddo

      !
      ! Compute kinetic energy in each cell
      !
      ke(:,:) = 0.0
      do iCell = 1, nCells
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k = 1, nVertLevels
               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
            end do
         end do
         do k = 1, nVertLevels
            ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
         end do
      end do

      !
      ! Compute v (tangential) velocities
      !
      v(:,:) = 0.0
      do iEdge = 1,nEdges
         do i = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1,nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do

#ifdef NCAR_FORMULATION
      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      vh(:,:) = 0.0
      do iEdge = 1, nEdgesSolve
         do j = 1, nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j,iEdge)
            do k = 1, nVertLevels
               vh(k,iEdge) = vh(k,iEdge) + weightsOnEdge(j,iEdge) * u(k,eoe) * h_edge(k,eoe)
            end do
         end do
      end do
#endif


      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
      !
      do iVertex = 1,nVertices
         do k = 1, nVertLevels
            h_vertex(k,iVertex) = 0.0
            do i = 1, vertexDegree
               h_vertex(k,iVertex) = h_vertex(k,iVertex) + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex(k,iVertex) = h_vertex(k,iVertex) / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex(k,iVertex)
         end do
      end do


      !
      ! Compute gradient of PV in the tangent direction
      !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
      !
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
           gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / &
                              dvEdge(iEdge)
         enddo
      enddo

      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
        do i = 1, vertexDegree
           iEdge = edgesOnVertex(i,iVertex)
           do k = 1, nVertLevels
              pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
           end do
        end do
      end do


      !
      ! Modify PV edge with upstream bias. 
      !
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
         enddo
      enddo


      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
      !
      pv_cell(:,:) = 0.0
      vorticity_cell(:,:) = 0.0
      do iVertex = 1, nVertices
       do i = 1, vertexDegree
         iCell = cellsOnVertex(i,iVertex)
         if (iCell <= nCells) then
           do k = 1, nVertLevels
             pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
             vorticity_cell(k,iCell) = vorticity_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * vorticity(k, iVertex) / areaCell(iCell)
           enddo
         endif
       enddo
      enddo


      !
      ! Compute gradient of PV in normal direction
      !   ( this computes gradPVn for all edges bounding real cells )
      !
      gradPVn(:,:) = 0.0
      do iEdge = 1, nEdges
        if( cellsOnEdge(1,iEdge) <= nCells .and. cellsOnEdge(2,iEdge) <= nCells) then
          do k = 1, nVertLevels
            gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / &
                                 dcEdge(iEdge)
          enddo
        endif
      enddo

      ! Modify PV edge with upstream bias.
      !
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
           pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * u(k,iEdge) * dt * gradPVn(k,iEdge)
         enddo
      enddo

      !
      ! set pv_edge = fEdge / h_edge at boundary points
      !
   !  if (maxval(boundaryEdge).ge.0) then
   !  do iEdge = 1,nEdges
   !     cell1 = cellsOnEdge(1,iEdge)
   !     cell2 = cellsOnEdge(2,iEdge)
   !     do k = 1,nVertLevels
   !       if(boundaryEdge(k,iEdge).eq.1) then
   !         v(k,iEdge) = 0.0
   !         if(cell1.gt.0) then
   !            h1 = h(k,cell1)
   !            pv_edge(k,iEdge) = fEdge(iEdge) / h1
   !            h_edge(k,iEdge) = h1
   !         else
   !            h2 = h(k,cell2)
   !            pv_edge(k,iEdge) = fEdge(iEdge) / h2
   !            h_edge(k,iEdge) = h2
   !         endif
   !       endif
   !     enddo
   !  enddo
   !  endif


   end subroutine sw_compute_solve_diagnostics


   subroutine sw_enforce_boundary_edge(tendPool, meshPool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Enforce any boundary conditions on the normal velocity at each edge
   !
   ! Input: grid - grid metadata
   !
   ! Output: tend_u set to zero at boundaryEdge == 1 locations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool
      type (mpas_pool_type), intent(in) :: meshPool

      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), pointer :: tend_u
      integer, pointer :: nCells, nEdges, nVertices, nVertLevels
      integer :: iEdge, k

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)
      call mpas_pool_get_array(tendPool, 'u', tend_u)

      if(maxval(boundaryEdge).le.0) return

      do iEdge = 1, nEdges
        do k = 1, nVertLevels

          if(boundaryEdge(k,iEdge).eq.1) then
             tend_u(k,iEdge) = 0.0
          endif

        enddo
       enddo

   end subroutine sw_enforce_boundary_edge

   subroutine sw_enforce_boundary_edge_lts(tendPool, meshPool, LTSPool, LTSRegion, LTSHalo)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Enforce any boundary conditions on the normal velocity at each edge
   !
   ! Input: grid - grid metadata
   !
   ! Output: tend_u set to zero at boundaryEdge == 1 locations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool
      type (mpas_pool_type), intent(in) :: meshPool, LTSPool
      type (integer), intent(in) :: LTSRegion, LTSHalo

      integer, dimension(:,:), pointer :: boundaryEdge
      real (kind=RKIND), dimension(:,:), pointer :: tend_u
      integer, pointer :: nVertices, nVertLevels
      integer :: iEdge, k, i

            integer, dimension(:,:), pointer :: nEdgesInLTSHalo
      integer, dimension(:,:,:), pointer :: edgesInLTSHalo
  
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)
      call mpas_pool_get_array(tendPool, 'u', tend_u)

      call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
      call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)

      if(maxval(boundaryEdge).le.0) return

      do i = 1, nEdgesInLTSHalo(LTSRegion,LTSHalo)
         iEdge = edgesInLTSHalo(LTSRegion,LTSHalo,i)
         do k = 1, nVertLevels

            if(boundaryEdge(k,iEdge).eq.1) then
               tend_u(k,iEdge) = 0.0
            end if

         end do
      end do

   end subroutine sw_enforce_boundary_edge_lts

end module sw_time_integration
