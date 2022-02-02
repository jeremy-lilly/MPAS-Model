! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_explicit_rayleigh_drag
!
!> \brief MPAS ocean explicit Rayleigh drag
!> \author Mark Petersen
!> \date   August 2017
!> \details
!>  This module contains the routine for computing
!>  tendencies from explicit Rayleigh drag.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_explicit_rayleigh_drag

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer

   use ocn_constants
   use ocn_forcing

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

   public :: ocn_vel_forcing_explicit_rayleigh_drag_tend, &
             ocn_vel_forcing_explicit_rayleigh_drag_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: explicitRayleighDragOn
   real (kind=RKIND) :: explicitRayleighDragCoef

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_rayleigh_drag_tend
!
!> \brief   Computes tendency term from explicit Rayleigh drag
!> \author  Mark Petersen
!> \date    15 August 2017
!> \details
!>  This routine computes the explicit Rayleigh drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_explicit_rayleigh_drag_tend(meshPool, normalVelocity, & !{{{
                                     kineticEnergyCell, layerThicknessEdge, tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,       &!< Input: velocity
         kineticEnergyCell,    &!< Input: kinetic energy at cell
         layerThicknessEdge     !< Input: thickness at edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, cell1, cell2, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      err = 0

      if ( .not. explicitRayleighDragOn ) return

      call mpas_timer_start('vel explicit rayleigh drag')

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      !$omp do schedule(runtime) private(k, cell1, cell2)
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        k =  maxLevelEdgeTop(iEdge)

        ! Explicit Rayleigh drag term:
        ! du/dt = ... - R u
        ! appied to bottom layer only.
        ! This term comes from the Rayleigh boundary condition in the vertical
        ! momentum mixing 

        tend(k,iEdge) = tend(k,iEdge) - explicitRayleighDragCoef * &
           normalVelocity(k,iEdge) 
      enddo
      !$omp end do

      call mpas_timer_stop('vel explicit rayleigh drag')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_rayleigh_drag_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_rayleigh_drag_init
!
!> \brief   Initializes ocean explicit Rayleigh drag forcing
!> \author  Mark Petersen
!> \date    August 2017
!> \details
!>  This routine initializes quantities related to explicit Rayleigh drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_explicit_rayleigh_drag_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_disable_vel_explicit_Rayleigh_drag
      logical, pointer :: config_use_explicit_Rayleigh_drag
      real (kind=RKIND), pointer :: config_explicit_Rayleigh_drag_coeff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_explicit_Rayleigh_drag', config_use_explicit_Rayleigh_drag)
      call mpas_pool_get_config(ocnConfigs, 'config_explicit_Rayleigh_drag_coeff', config_explicit_Rayleigh_drag_coeff)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_explicit_Rayleigh_drag', config_disable_vel_explicit_Rayleigh_drag)

      explicitRayleighDragCoef = 0.0_RKIND

      if (config_use_explicit_Rayleigh_drag) then
          explicitRayleighDragOn = .true.
          explicitRayleighDragCoef = config_explicit_Rayleigh_drag_coeff
      endif

      if (config_disable_vel_explicit_Rayleigh_drag) explicitRayleighDragOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_rayleigh_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_explicit_rayleigh_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
