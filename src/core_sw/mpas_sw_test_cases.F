! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module sw_test_cases

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_log
   use mpas_timer


   contains


   subroutine setup_sw_test_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the shallow water test case 
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: statePool, meshPool, LTSPool

      integer, pointer :: config_test_case
      logical, pointer:: config_use_local_time_stepping      

      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)
      call mpas_pool_get_config(domain % configs, 'config_use_local_time_stepping', config_use_local_time_stepping)

      if (config_test_case == 0) then
         call mpas_log_write('Using initial conditions supplied in input file')

      else if (config_test_case == 1) then
         call mpas_log_write('Setting up shallow water test case 1')
         call mpas_log_write(' -- Advection of Cosine Bell over the Pole')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_1(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 2) then
         call mpas_log_write('Setting up shallow water test case 2')
         call mpas_log_write(' -- Setup shallow water test case 2: Global Steady State Nonlinear Zonal Geostrophic Flow')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            !call compute_weightsOnEdgeFortran(meshPool) decomment to use in-house weights

            call sw_test_case_2(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 5) then
         call mpas_log_write('Setting up shallow water test case 5')
         call mpas_log_write(' -- Setup shallow water test case 5: Zonal Flow over an Isolated Mountain')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_5(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 6) then
         call mpas_log_write('Setting up shallow water test case 6')
         call mpas_log_write(' -- Rossby-Haurwitz Wave')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_6(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else
         call mpas_log_write('Only test case 1, 2, 5, and 6 are currently supported.')
         stop
      end if

      if (config_use_local_time_stepping) then
         call mpas_log_write('Setting up local time step regions and halos')

call mpas_timer_start("init LTS")
         call mpas_init_LTS(domain)
call mpas_timer_stop("init LTS")

      end if

   end subroutine setup_sw_test_case


   subroutine sw_test_case_1(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 1: Advection of Cosine Bell over the Pole
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: h0 = 1000.0
      real (kind=RKIND), parameter :: theta_c = 0.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: alpha = pii/4.0

      integer :: iCell, iEdge, iVtx
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: r, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)

      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(latVertex(iVtx)) * cos(alpha) - &
                                       cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize cosine bell at (theta_c, lambda_c)
      !
      do iCell = 1, nCells
         r = sphere_distance(theta_c, lambda_c, latCell(iCell), lonCell(iCell), a) 
         if (r < a/3.0) then
            h(1,iCell) = (h0 / 2.0) * (1.0 + cos(pii*r*3.0/a))
         else
            h(1,iCell) = h0 / 2.0
         end if
      end do


   end subroutine sw_test_case_1


   subroutine sw_test_case_2(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 2: Global Steady State Nonlinear Zonal 
   !                                  Geostrophic Flow
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: gh0 = 29400.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h


      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
    
      !xCell = xCell * a
      !yCell = yCell * a
      !zCell = zCell * a
      !xVertex = xVertex * a
      !yVertex = yVertex * a
      !zVertex = zVertex * a
      !xEdge = xEdge * a
      !yEdge = yEdge * a
      !zEdge = zEdge * a
      !dvEdge = dvEdge * a
      !dcEdge = dcEdge * a
      !areaCell = areaCell * a**2.0
      !areaTriangle = areaTriangle * a**2.0
      !kiteAreasOnVertex = kiteAreasOnVertex * a**2.0
      

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                 sin(latVertex(iVtx)) * cos(alpha) - &
                                 cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                               )
      end do
      do iEdge = 1,nEdges
         u(1,iEdge) = -1.0 * ( &
                                 psiVertex(verticesOnEdge(2,iEdge)) - &
                                 psiVertex(verticesOnEdge(1,iEdge)) &
                               ) / dvEdge(iEdge)

      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      do iEdge = 1, nEdges
         fEdge(iEdge) = 2.0 * omega * &
                            ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha) + &
                              sin(latEdge(iEdge)) * cos(alpha) &
                            )
      end do
      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0 * omega * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) + &
                                          sin(latVertex(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                             (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                                              sin(latCell(iCell)) * cos(alpha) &
                                             )**2.0 &
                                      ) / &
                                      gravity
      end do

   end subroutine sw_test_case_2


   subroutine sw_test_case_5(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 5: Zonal Flow over an Isolated Mountain
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 20.
      real (kind=RKIND), parameter :: gh0 = 5960.0*gravity
      real (kind=RKIND), parameter :: hs0 = 2000.
      real (kind=RKIND), parameter :: theta_c = pii/6.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rr = pii/9.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx, nTracers, k
      integer, pointer :: nVertices, nEdges, nCells, nVertLevels
      real (kind=RKIND) :: r, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, h_s, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
      nTracers = size(tracers, dim=1)

      !xCell = xCell * a
      !yCell = yCell * a
      !zCell = zCell * a
      !xVertex = xVertex * a
      !yVertex = yVertex * a
      !zVertex = zVertex * a
      !xEdge = xEdge * a
      !yEdge = yEdge * a
      !zEdge = zEdge * a
      !dvEdge = dvEdge * a
      !dcEdge = dcEdge * a
      !areaCell = areaCell * a**2.0
      !areaTriangle = areaTriangle * a**2.0
      !kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(latVertex(iVtx)) * cos(alpha) - &
                                       cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                     )
      end do
      
      do k = 1, nVertLevels 
         do iEdge = 1, nEdges
            u(k,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
         end do
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      do iEdge = 1, nEdges
         fEdge(iEdge) = 2.0 * omega * &
                                        (-cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha) + &
                                          sin(latEdge(iEdge)) * cos(alpha) &
                                        )
      end do
      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0 * omega * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) + &
                                           sin(latVertex(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize mountain
      !
      do iCell = 1, nCells
         if (lonCell(iCell) < 0.0) lonCell(iCell) = lonCell(iCell) + 2.0 * pii
         r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
         h_s(iCell) = hs0 * (1.0 - r/rr)
      end do

      !
      ! Initialize tracer fields
      !
      do k = 1, nVertLevels
         do iCell = 1, nCells
               r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
               tracers(1,k,iCell) = 1.0 - r/rr
         end do
         if (nTracers > 1) then
            do iCell = 1, nCells
               r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + &
                         (latCell(iCell) - theta_c - pii/6.0)**2.0 &
                        ) &
                    )
               tracers(2,k,iCell) = 1.0 - r/rr
            end do
         end if
      end do

      !
      ! Initialize height field (actually, fluid thickness field)
      !

      do k = 1, nVertLevels
         do iCell = 1, nCells
            h(k,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                         (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                                           sin(latCell(iCell)) * cos(alpha) &
                                         )**2.0 &
                                      ) / &
                                      gravity
            h(k,iCell) = h(k,iCell) - h_s(iCell)
         end do
      end do

   end subroutine sw_test_case_5


   subroutine sw_test_case_6(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 6: Rossby-Haurwitz Wave
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: h0 = 8000.0
      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      integer :: iCell, iEdge, iVtx, nTracers
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, h_s
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers



      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
      nTracers = size(tracers, dim=1)


      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * a * w * sin(latVertex(iVtx)) + &
                            a *a * K * (cos(latVertex(iVtx))**R) * &
                            sin(latVertex(iVtx)) * cos(R * lonVertex(iVtx))
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gravity * h0 + a*a*aa(latCell(iCell)) + &
                                                      a*a*bb(latCell(iCell)) * cos(R*lonCell(iCell)) + &
                                                      a*a*cc(latCell(iCell)) * cos(2.0*R*lonCell(iCell)) &
                                      ) / gravity
      end do

   end subroutine sw_test_case_6

   subroutine mpas_init_LTS(domain)
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
   integer :: i, iCell, iEdge, iVertex, iRegion, iHalo, iHaloP1, iCellHalo, cell1, cell2, nLTSHalosCopy, haveFineCounter, haveCoarseCounter
   integer, pointer :: nLTSHalos, nCells, nCellsSolve, nEdges, nEdgesSolve, haveFine, haveCoarse, haveInterface, doExtraLayers, nLTSHalosPExtra, moreCellsOnInterface
   real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, areaCell
   integer, dimension(:), pointer :: LTSRegion, LTSRegionLocal, nEdgesOnCell, minMaxLTSRegion, minMaxAfterHaloUpdate
   integer, dimension(:,:), pointer :: cellsOnEdge,  edgesOnCell, verticesOnCell, nCellsInLTSHalo, nEdgesInLTSHalo
   integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
   integer, dimension(:), pointer :: mpiRankCell, blockPartitionCell
   real (kind=RKIND) :: latPoint, lonPoint, distLTS
   
   call mpas_pool_get_subpool(domain % blocklist % structs, 'LTS', LTSPool)
   call mpas_pool_get_dimension(LTSPool, 'nLTSHalos', nLTSHalos)
   call mpas_pool_get_dimension(LTSPool, 'moreCellsOnInterface', moreCellsOnInterface)
   call mpas_pool_get_dimension(LTSPool, 'nLTSHalosPExtra', nLTSHalosPExtra)

   ! set to 1 if need more cells in interface layers (for load balancing)
   !moreCellsOnInterface = 0

   nLTSHalosCopy = nLTSHalos
   if (moreCellsOnInterface == 1) then
      !nLTSHalosCopy = 78 !for CPU test with multiple interface layers for test 5 (it's 78 for LTS2, 80 for LTS3)
      nLTSHalosCopy = nLTSHalosPExtra !for speed up test with multiple interface layers for test 5 (it's 18 (or 38) for LTS2, 20 (or 40) for LTS3)
   end if

   !NOTE:
   ! the word "halo" in LTS instances is not related to the MPI halo.

   !NOTE:
   ! to run LTS2, need to set nLTSHalos=1 in the Registry
   ! to run LTS3, need to set nLTSHalos=3 in the Registry

   !NOTE:
   ! for performance, LTS should to be run with either one of the following setups:
   ! 1) N MPI partitions and 3N blocks, each blocks owns an LTS region (fine, coarse or interface)
   ! 2) N MPI partitions and N blocks,  each blocks owns an LTS region (fine, coarse or interface) 
   ! the first case is better theoretically because it optimizes the parallel efficiency by minimizing the wait time among processors.
   ! There are python scripts to modify the graph.info and the graph.info.part.X files to achieve 1) or 2).
   ! Nevertheless, this code is written so that LTS can run with virtually any graph.info file.

   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'mpiRankCell', mpiRankCell)
      call mpas_pool_get_array(meshPool, 'blockPartitionCell', blockPartitionCell)

      call mpas_pool_get_array(LTSPool, 'LTSRegion', LTSRegion)
      call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
      call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
      call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
      call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
      call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)
      call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
      call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
      call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
      call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

      ! the loop below initializes the LTS regions for the current processor (or block)
      ! LTS Regions code:
      ! 1 = fine
      ! 2 = coarse
      ! 3 = interface layer 1
      ! 4 = interface layer 2
      ! then if nLTSHalos=3 we also have
      ! 5 = fine (to advance when doing 1st stage on interface)
      ! 6 = coarse
      ! 7 = fine (to advance when doing 1st stage on interface)
      ! 8 = coarse
      
      LTSRegion(:) = 1
      do iCell = 1, nCells
         !if (2.0 * sqrt(areaCell(iCell) / pii) .ge. 300000) then
         !   LTSRegion(iCell) = 2
         !end if

         latPoint = pii / 6.0
         lonPoint = 3.0 * pii * 0.5
         distLTS = sphere_distance(latCell(iCell), lonCell(iCell), latPoint, lonPoint, 1.0)
         ! NOTE: the number on the right hand side of the inequality dictates where is the interface between coarse and fine on the mesh
         !if (distLTS > 1.2) then  !this was for the convergence test for test 5
         !if (distLTS > 0.35) then   !this is for the CPU ratio test for test 5
         !if (distLTS > 0.7) then   !this is for the CPU ratio test with multiple interface layers for test 5
         if (distLTS > 0.55) then   !this is for the speed up test for test 5
         !if (distLTS > 0.55) then   !this is for the speed up test with multiple interface layers for test 5 (this goes with 38 or 40 layers)
         !if (distLTS > 0.3) then   !this is for the speed up test with multiple interface layers for test 5 (this goes with 18 or 20 layers)
            LTSRegion(iCell) = 2
         end if
         !BEGIN: test initializaton
         !if (latCell(iCell) < 0.) then
         !   LTSRegion(iCell) = 1
         !else
         !   LTSRegion(iCell) = 2
         !end if
         !END: test initialization
      end do

      nCellsInLTSHalo(:,:) = 0
      nEdgesInLTSHalo(:,:) = 0
      LTSRegionLocal(:) = LTSRegion(:) 

      haveFine = 1
      haveCoarse = 1
      haveInterface = 1
      doExtraLayers = 0

      ! this is to visualize MPI partitions on Paraview
      mpiRankCell(1:nCells) = domain % dminfo % my_proc_id

      ! this is to visualize block partitions on Paraview
      do iCell = 1, nCells
         blockPartitionCell(iCell) = block % blockID
      end do
      
      minMaxLTSRegion(1) = minval(LTSRegionLocal(1:nCells))  
      minMaxLTSRegion(2) = maxval(LTSRegionLocal(1:nCells))  

      if (minMaxLTSRegion(2) - minMaxLTSRegion(1) > 1) then
         call mpas_log_write('minMaxLTSRegion(2) - minMaxLTSRegion(1) > 1 is not allowed. Change LTSRegion initialization in init_LTS subroutine.')
         stop
      end if

      if (minMaxLTSRegion(2) - minMaxLTSRegion(1) == 1) then  !this means there are two regions, a fine and a coarse

         doExtraLayers = 1

         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            ! here we write the first layers of the LTS halo (interface layer 1 and interface layer 2)
            if (LTSRegionLocal(cell1) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell2) == minMaxLTSRegion(2)) then
               if (cell1 .le. nCells) then
                  LTSRegionLocal(cell1) = minMaxLTSRegion(1) + 2
               end if
               if (cell2 .le. nCells) then
                  LTSRegionLocal(cell2) = minMaxLTSRegion(2) + 2
               end if
            else if (LTSRegionLocal(cell2) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell1) == minMaxLTSRegion(2)) then
               if (cell2 .le. nCells) then
                  LTSRegionLocal(cell2) = minMaxLTSRegion(1) + 2
               end if
               if (cell1 .le. nCells) then
                  LTSRegionLocal(cell1) = minMaxLTSRegion(2) + 2
               end if
            ! the next if statements are there because above we are possibly changing two LTSRegionLocal at the
            ! time, so we have to take into consideration the fact that an edge
            ! might be shared by two elements, one of which had its LTSRegionLocal already changed
            else if (LTSRegionLocal(cell1) == minMaxLTSRegion(1) + 2 .and. LTSRegionLocal(cell2) == minMaxLTSRegion(2)) then
               if (cell2 .le. nCells) then
                  LTSRegionLocal(cell2) = minMaxLTSRegion(2) + 2
               end if
            else if (LTSRegionLocal(cell2) == minMaxLTSRegion(1) + 2 .and. LTSRegionLocal(cell1) == minMaxLTSRegion(2)) then
               if (cell1 .le. nCells) then
                  LTSRegionLocal(cell1) = minMaxLTSRegion(2) + 2
               end if
            else if (LTSRegionLocal(cell1) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell2) == minMaxLTSRegion(2) + 2) then
               if ( cell1 .le. nCells) then
                  LTSRegionLocal(cell1) = minMaxLTSRegion(1) + 2
               end if
            else if (LTSRegionLocal(cell2) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell1) == minMaxLTSRegion(2) + 2) then
               if (cell2 .le. nCells) then
                  LTSRegionLocal(cell2) = minMaxLTSRegion(1) + 2
               end if
            end if
         end do

      end if   

      block => block % next
   end do

   call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')

   ! this below is to see if a processor that started with only fine or only coarse needs to fill out some extra LTS halo layers  
   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
      call mpas_pool_get_array(LTSPool, 'minMaxAfterHaloUpdate', minMaxAfterHaloUpdate)
      call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
      call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)

      minMaxAfterHaloUpdate(1) = minval(LTSRegionLocal(1:nCells))  
      minMaxAfterHaloUpdate(2) = maxval(LTSRegionLocal(1:nCells)) 

      if (minMaxLTSRegion(2) - minMaxLTSRegion(1) == 0) then

         if (minMaxAfterHaloUpdate(1) .ne. minMaxAfterHaloUpdate(2)) then
            doExtraLayers = 1
         end if

         if (minMaxAfterHaloUpdate(1) .ne. minMaxLTSRegion(1)) then
            call mpas_log_write('minMaxAfterHaloUpdate(1) not equal to minMaxLTSRegion(1), something is wrong, see init_LTS subroutine')
         end if
      end if 

      block => block % next
   end do

   ! this is a loop to change LTSRegionLocal for the extra halo regions that might be
   ! needed for higher order LTS schemes (for instance the 3rd order LTS (scheme))
   do iRegion = 1,2
      do iHalo = 1, nLTSHalosCopy - 1

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

            call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
            call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
            call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)

            if (doExtraLayers == 1) then 

               do iCell = 1, nCells
                  if (LTSRegionLocal(iCell) == (minMaxLTSRegion(iRegion) + 2 * iHalo)) then
                     do i = 1, nEdgesOnCell(iCell) 
                        iEdge = edgesOnCell(i, iCell)  
                        cell1 = cellsOnEdge(1,iEdge)
                        cell2 = cellsOnEdge(2,iEdge)
                        if (LTSRegionLocal(cell1) == minMaxLTSRegion(iRegion) ) then
                           if (cell1 .le. nCells) then
                              LTSRegionLocal(cell1) = LTSRegionLocal(cell1) +  2 * (iHalo+1)
                           end if
                        else if (LTSRegionLocal(cell2) == minMaxLTSRegion(iRegion)) then
                           if (cell2 .le. nCells) then
                              LTSRegionLocal(cell2) = LTSRegionLocal(cell2) +  2 * (iHalo+1)
                           end if
                        end if
                     end do
                  end if              
               end do   

            end if

            block => block % next
         end do


         call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')


         ! this below is to see if a processor that started with only fine or only coarse needs to fill out some extra LTS halo layers  
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

            call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

            call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
            call mpas_pool_get_array(LTSPool, 'minMaxAfterHaloUpdate', minMaxAfterHaloUpdate)
            call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
            call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)

            minMaxAfterHaloUpdate(1) = minval(LTSRegionLocal(1:nCells))
            minMaxAfterHaloUpdate(2) = maxval(LTSRegionLocal(1:nCells))

            if (doExtraLayers == 0) then

               if (minMaxAfterHaloUpdate(1) .ne. minMaxAfterHaloUpdate(2)) then
                  doExtraLayers = 1
               end if

               if (minMaxAfterHaloUpdate(1) .ne. minMaxLTSRegion(1)) then
                  call mpas_log_write('minMaxAfterHaloUpdate(1) not equal to minMaxLTSRegion(1), something is wrong, see init_LTS subroutine')
                  stop
               end if
            end if

            block => block % next
         end do

      end do
   end do


   call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')

   if (moreCellsOnInterface == 1) then
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(LTSPool, 'nLTSHalos', nLTSHalos)
         call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
    
         if (nLTSHalos == 1) then !this means we are using LTS2 
            do iCell = 1, nCells
               if( (mod(LTSRegionLocal(iCell),2) == 1) .and. (LTSRegionLocal(iCell) .ne. 1) ) then !if it's odd and not 1(fine) then must be 3(interface layer 1)
                  LTSRegionLocal(iCell) = 3
               else if ( (mod(LTSRegionLocal(iCell),2) == 0) .and. (LTSRegionLocal(iCell) .ne. 2) ) then !if it's even and not 2(coarse) then must be 4(interface layer 2)
                  LTSRegionLocal(iCell) = 4
               end if
            end do
         else if (nLTSHalos == 3) then !this means we are using LTS3
            do iCell = 1, nCells
               if( (mod(LTSRegionLocal(iCell),2) == 1) .and. (LTSRegionLocal(iCell) .ne. 1) ) then 
                  ! if we are here it could be either interface 1 or those two layers of fine we need for the third order         
                  if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy + 1) then
                     LTSRegionLocal(iCell) = 7 !this goes in cellsInLTSHalo(1,4)
                  else if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy - 1) then
                     LTSRegionLocal(iCell) = 5 !this goes in cellsInLTSHalo(1,3)
                  else
                      LTSRegionLocal(iCell) = 3 !this goes in cellsInLTSHalo(1,2)
                  end if
               else if ( (mod(LTSRegionLocal(iCell),2) == 0) .and. (LTSRegionLocal(iCell) .ne. 2) ) then
                  ! if we are here it could be either interface 2 or those two layers of coarse
                  if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy + 2) then
                     LTSRegionLocal(iCell) = 8 !this goes in cellsInLTSHalo(2,4)
                  else if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy) then
                     LTSRegionLocal(iCell) = 6 !this goes in cellsInLTSHalo(2,3)
                  else
                      LTSRegionLocal(iCell) = 4 !this goes in cellsInLTSHalo(2,2)
                  end if
               end if
            end do
         end if

         block => block % next
      end do

   call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')
   end if

   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(LTSPool, 'nLTSHalos', nLTSHalos)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
      call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
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
               do iHalo = 1, nLTSHalos
                  if(LTSRegionLocal(iCell) == (minMaxLTSRegion(iRegion) + 2 * iHalo) ) then
                     nCellsInLTSHalo(iRegion,iHalo+1) = nCellsInLTSHalo(iRegion,iHalo+1) + 1
                     cellsInLTSHalo(iRegion,iHalo+1,nCellsInLTSHalo(iRegion,iHalo+1)) = iCell
                  end if
               end do
            end if
         end do
      end do

      haveFineCounter = 0
      do iHalo = 1, nLTSHalos+1
         if ((nCellsInLTSHalo(1,iHalo) == 0) .and. (iHalo .ne. 2)  ) then
            haveFineCounter = haveFineCounter + 1
         end if
      end do
      if (haveFineCounter == nLTSHalos) then
         haveFine = 0
      end if

      haveCoarseCounter = 0
      do iHalo = 1, nLTSHalos+1
         if ((nCellsInLTSHalo(2,iHalo) == 0) .and. (iHalo .ne. 2)  ) then
            haveCoarseCounter = haveCoarseCounter + 1
         end if
      end do
      if (haveCoarseCounter == nLTSHalos) then
         haveCoarse = 0
      end if

      if ((nCellsInLTSHalo(1,2) == 0) .and. (nCellsInLTSHalo(2,2) == 0)) then
         haveInterface = 0
      end if

      ! below we fill out the lists for the edges, according to the LTSRegionLocal that have been assigned to the cells.
      ! we move from the fine to the coarse (i.e. from the fine to the nearest LTS region in the direction of the coarse). 
      ! Note that edges shared between cells of different LTS regions are owned by the cell in the LTS region closest to the fine region, 
      ! see Figure 3 in "Conservative explicit local time-stepping schemes for the shallow water equations" by Hoang et al.
    
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
   
      do iHalo = 1, nLTSHalos
         do iCell = 1, nCellsInLTSHalo(1, nLTSHalos - iHalo + 2)
            do i = 1, nEdgesOnCell(cellsInLTSHalo(1, nLTSHalos- iHalo + 2, iCell))
               iEdge = edgesOnCell(i,cellsInLTSHalo(1, nLTSHalos - iHalo + 2, iCell))
               if (iEdge .le. nEdgesSolve) then
                  if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                     nEdgesInLTSHalo(1, nLTSHalos - iHalo + 2) = nEdgesInLTSHalo(1, nLTSHalos - iHalo + 2) + 1
                     edgesInLTSHalo(1, nLTSHalos - iHalo + 2,  nEdgesInLTSHalo(1, nLTSHalos - iHalo + 2)) = iEdge
                     isLTSRegionEdgeAssigned(iEdge) = 1
                  end if
               end if
            end do
         end do
      end do

      do iHalo = 1, nLTSHalos
         do iCell = 1, nCellsInLTSHalo(2, iHalo+1)
            do i = 1, nEdgesOnCell(cellsInLTSHalo(2,iHalo+1,iCell))
               iEdge = edgesOnCell(i,cellsInLTSHalo(2,iHalo+1,iCell))
               if (iEdge .le. nEdgesSolve) then
                  if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                     nEdgesInLTSHalo(2, iHalo+1) = nEdgesInLTSHalo(2, iHalo+1) + 1
                     edgesInLTSHalo(2, iHalo+1,  nEdgesInLTSHalo(2, iHalo+1)) = iEdge
                     isLTSRegionEdgeAssigned(iEdge) = 1
                  end if
               end if
            end do
         end do
      end do
 
      do iCell = 1, nCellsInLTSHalo(2,1)
         do i = 1, nEdgesOnCell(cellsInLTSHalo(2,1,iCell))
            iEdge = edgesOnCell(i,cellsInLTSHalo(2,1,iCell))
            if (iEdge .le. nEdgesSolve) then
               if(isLTSRegionEdgeAssigned(iEdge) == 0) then
                  nEdgesInLTSHalo(2,1) = nEdgesInLTSHalo(2,1) + 1
                  edgesInLTSHalo(2,1,nEdgesInLTSHalo(2,1)) = iEdge
                  isLTSRegionEdgeAssigned(iEdge) = 1
               end if
            end if
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
   
      !if(block % blockID == 5) then 
      !if(haveInterface  == 1) then 
      !if (domain % dminfo % my_proc_id == 0) then
      !   print*, block % blockID
      !   print*, 'haveFine=', haveFine, 'haveCoarse=', haveCoarse, 'haveInterface=', haveInterface
      !   print*, 'nCellsInLTSHalo(1,1)=', nCellsInLTSHalo(1,1), 'nCellsInLTSHalo(1,3)=', nCellsInLTSHalo(1,3), 'nCellsInLTSHalo(1,4)=', nCellsInLTSHalo(1,4)
      !   print*, 'nCellsInLTSHalo(2,1)=', nCellsInLTSHalo(2,1), 'nCellsInLTSHalo(2,3)=', nCellsInLTSHalo(2,3), 'nCellsInLTSHalo(2,4)=', nCellsInLTSHalo(2,4)
      !   print*, 'nCellsInLTSHalo(1,2)=', nCellsInLTSHalo(1,2), 'nCellsInLTSHalo(2,2)=', nCellsInLTSHalo(2,2)
      !   print*, 'nEdgessInLTSHalo(1,1)=', nEdgesInLTSHalo(1,1), 'nEdgesInLTSHalo(1,3)=', nEdgesInLTSHalo(1,3), 'nEdgesInLTSHalo(1,4)=', nEdgesInLTSHalo(1,4)
      !   print*, 'nEdgessInLTSHalo(2,1)=', nEdgesInLTSHalo(2,1), 'nEdgesInLTSHalo(2,3)=', nEdgesInLTSHalo(2,3), 'nEdgesInLTSHalo(2,4)=', nEdgesInLTSHalo(2,4)
      !   print*, 'nEdgessInLTSHalo(1,2)=', nEdgesInLTSHalo(1,2), 'nEdgesInLTSHalo(2,2)=', nEdgesInLTSHalo(2,2)
      !end if

      !TO ERASE
      !if(haveFine == 1) then 
      !   if (nCellsInLTSHalo(1,1) + nCellsInLTSHalo(1,3) + nCellsInLTSHalo(1,4)  < 100) then
      !   !if (nCellsInLTSHalo(1,1) < 100) then
      !      print*, 'WARNING: fine block', block % blockID, 'has', nCellsInLTSHalo(1,1) + nCellsInLTSHalo(1,3) + nCellsInLTSHalo(1,4), 'cells'
      !      !print*, 'WARNING: fine block', block % blockID, 'has', nCellsInLTSHalo(1,1), 'cells'
      !   end if
      !end if
      !END

      if(haveInterface == 1) then 
         if (nCellsInLTSHalo(1,2) + nCellsInLTSHalo(2,2) < 100) then
            print*, 'WARNING: interface block', block % blockID, 'has', nCellsInLTSHalo(1,2) + nCellsInLTSHalo(2,2), 'cells'
         end if
      end if
   
      block => block % next
   end do
   ! END: CHECK TO SEE IF THE INTERFACE BLOCKS HAVE ENOUGH CELLS

   end subroutine mpas_init_LTS


   ! THIS ROUTINE WAS INTRODUCED ONLY FOR DEBUGGING REASONS AND SHOULD NOT BE PUSHED ULTIMATELY
   subroutine compute_weightsOnEdgeFortran(meshPool)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Computing the weights for the Coriolis term in the equation of u 
   !         
   !
   ! Reference: Pedro S. Peixoto, "Accuracy analysis of mimetic finite
   !            volume operators on geodesic grids and a consistent alternative" 
   !            J. of Comp. Phys., 310, pp. 127--160
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool 

      integer :: iEdge, jEdge, i, ii, jj, v2, sharedCell, v2Edge, jEdgeIsNextEdge, nextEdge, nextEdgeOld
      integer, pointer :: nEdges, vertexDegree

      integer, dimension(2) :: iEdgeCells, jEdgeCells, jEdgeVertices, nextEdgeCells
 
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      integer, dimension(:,:), pointer :: verticesOnEdge, edgesOnCell, cellsOnEdge, edgesOnEdge, edgesOnVertex

      real (kind=RKIND) :: areaSum, laCell, loCell, laVertex, loVertex, laEdge, loEdge, c, n, zMidpoint, yMidpoint, xMidpoint, deltaLon, Bx, By
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, xCell, yCell, zCell

      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge
 
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)

      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)

      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)

      weightsOnEdge(:,:) = 0.0
      do iEdge = 1, nEdges
         iEdgeCells(1) = cellsOnEdge(1,iEdge)
         iEdgeCells(2) = cellsOnEdge(2,iEdge)
         
         v2 = verticesOnEdge(1,iEdge) !this is v2 in the Thuburn paper  
         
         ! ---- BEGIN: compute c
         do ii = 1, vertexDegree
            v2Edge = edgesOnVertex(ii, v2)
            if (v2Edge == iEdge) then 
               ! Vector points from vertex 1 to vertex 2
               if(v2 == verticesOnEdge(1, iEdge)) then
                  c  = - 1.0 
               else
                  c  =  1.0
               end if
            end if
         end do       
         c = - c 
         ! ----- END: compute c

         do i = 1, nEdgesOnEdge(iEdge)
            areaSum = 0.0
            v2 = verticesOnEdge(1,iEdge)
            jEdge = edgesOnEdge(i,iEdge)
            jEdgeCells(1) = cellsOnEdge(1,jEdge)
            jEdgeCells(2) = cellsOnEdge(2,jEdge)
            do ii = 1,2
               do jj = 1,2
                  if (iEdgeCells(ii)  == jEdgeCells(jj)) then
                     sharedCell = iEdgeCells(ii) !this is to find what is the cell in common (in what cell jEdge is)
                  end if
               end do
            end do         

            ! ---- BEGIN: compute n (note we are using the MPAS convention, opposite to the 2010 Ringler paper)
            if (sharedCell < jEdgeCells(1) .or. sharedCell < jEdgeCells(2)) then
               n = - 1.0 
            else if (sharedCell > jEdgeCells(1) .or. sharedCell > jEdgeCells(2)) then
               n = 1.0
            end if
            ! --- END: compute n

            ! BEGIN: compute areaSum 
            ! here we compute the area of the only triangle associated with
            ! iEdge
            laCell = latCell(sharedCell)
            loCell = lonCell(sharedCell) 
            laVertex = latVertex(v2)
            loVertex = lonVertex(v2)
               
            deltaLon = lonCell(iEdgeCells(2)) - lonCell(iEdgeCells(1))
            Bx = cos(latCell(iEdgeCells(2))) * cos(deltaLon)
            By = cos(latCell(iEdgeCells(2))) * sin(deltaLon)
            laEdge = atan2(sin(latCell(iEdgeCells(1))) + sin(latCell(iEdgeCells(2))) , sqrt( (cos(latCell(iEdgeCells(1))) + & 
                     Bx) * (cos(latCell(iEdgeCells(1))) + Bx) + By * By ))
            loEdge = lonCell(iEdgeCells(1)) + atan2(By, cos(latCell(iEdgeCells(1))) + Bx)

            !BEGIN: COMMENT TO DO CELL AREA TEST
            areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell) / areaCell(sharedCell)
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex,
            !loVertex, laCell, loCell)
            !END: UNCOMMENT TO DO CELL AREA TEST

            jEdgeVertices(1)  = verticesOnEdge(1, jEdge)
            jEdgeVertices(2)  = verticesOnEdge(2, jEdge)
            
            !BEGIN: TO REMOVE
            !print*, 'v2 =', v2, 'jEdgeVertices= ', jEdgeVertices(1),
            !jEdgeVertices(2)
            !END: TO REMOVE
            
            !BEGIN: COMMENT TO DO CELL AREA TEST
            if ((jEdgeVertices(1) == v2) .or. (jEdgeVertices(2) == v2)) then 
               jEdgeIsNextEdge = 1
            else
               jEdgeIsNextEdge = 0
            end if
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !jEdgeIsNextEdge = 0
            !END: UNCOMMENT TO DO CELL AREA TEST

            if (jEdgeIsNextEdge == 0) then
               do ii = 1, vertexDegree
                  nextEdge = edgesOnVertex(ii,v2)
                  if ((cellsOnEdge(1,nextEdge) == sharedCell .or. cellsOnEdge(2,nextEdge) == sharedCell) .and. (nextEdge .ne. iEdge )) then
                     exit
                  end if
               end do
            end if

            !BEGIN: TO REMOVE
            !print*, 'BEFORE THE WHILE LOOP'
            !END: TO REMOVE

            nextEdgeOld = iEdge
            do while (jEdgeIsNextEdge == 0)

               !BEGIN: TO REMOVE
               !print*, 'nextEdge =', nextEdge, 'nextEdgeVertices=',
               !verticesOnEdge(1,nextEdge), verticesOnEdge(2,nextEdge) 
               !END: TO REMOVE

               ! this computes the areas of the two triangles associated with
               ! nextEdge
               nextEdgeCells(1) = cellsOnEdge(1,nextEdge)
               nextEdgeCells(2) = cellsOnEdge(2,nextEdge)

               deltaLon = lonCell(nextEdgeCells(2)) - lonCell(nextEdgeCells(1))
               Bx = cos(latCell(nextEdgeCells(2))) * cos(deltaLon)
               By = cos(latCell(nextEdgeCells(2))) * sin(deltaLon)
               laEdge = atan2(sin(latCell(nextEdgeCells(1))) + sin(latCell(nextEdgeCells(2))) , sqrt( (cos(latCell(nextEdgeCells(1))) + &
                        Bx) * (cos(latCell(nextEdgeCells(1))) + Bx) + By * By ))
               loEdge = lonCell(nextEdgeCells(1)) + atan2(By, cos(latCell(nextEdgeCells(1))) + Bx)

               do ii=1,2 
                  laVertex = latVertex(verticesOnEdge(ii, nextEdge))
                  loVertex = lonVertex(verticesOnEdge(ii, nextEdge))

                  !BEGIN: COMMENT TO DO CELL AREA TEST
                  areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell) / areaCell(sharedCell)
                  !END: COMMENT TO DO CELL AREA TEST

                  !BEGIN: UNCOMMENT TO DO CELL AREA TEST
                  !areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge,
                  !laVertex, loVertex, laCell, loCell)
                  !END: UNCOMMENT TO DO CELL AREA TEST
             
               end do

               ! new we have to find the nextEdge
               do ii=1,2
                  if (verticesOnEdge(ii,nextEdge) .ne. v2) then
                     v2 = verticesOnEdge(ii,nextEdge) !this is necessary to find the nextEdge
                     !BEGIN: TO REMOVE
                     !print*, 'v2=', v2
                     !END: TO REMOVE
                     exit
                  end if
               end do

               !BEGIN: TO REMOVE
               !print*, 'iEdgeVertices =', v2, verticesOnEdge(2,iEdge),
               !'jEdgeVertices= ', jEdgeVertices(1), jEdgeVertices(2)
               !END: TO REMOVE

               !BEGIN: COMMENT TO DO CELL AREA TEST
               if (jEdgeVertices(1) == v2 .or. jEdgeVertices(2) == v2) then
                  jEdgeIsNextEdge = 1
               end if
               !END: COMMENT TO DO CELL AREA TEST

               !BEGIN: UNCOMMENT TO DO CELL AREA TEST
               !if (v2 == verticesOnEdge(2,iEdge)) then
               !   jEdgeIsNextEdge = 1
               !end if
               !END: UNCOMMENT TO DO CELL AREA TEST

               if (jEdgeIsNextEdge == 0) then
                  nextEdgeOld = nextEdge
                  do ii = 1, vertexDegree
                     nextEdge = edgesOnVertex(ii,v2)
                     if ((cellsOnEdge(1,nextEdge) == sharedCell .or. cellsOnEdge(2,nextEdge) == sharedCell) .and. (nextEdge .ne. nextEdgeOld )) then
                        exit
                     end if
                  end do
               end if

            end do
 
            !BEGIN: TO REMOVE
            !print*, 'AFTER THE WHILE LOOP'
            !END: TO REMOVE
            
            ! here we compute the area of the only triangle associated with
            ! jEdge
         
            laVertex = latVertex(v2)
            loVertex = lonVertex(v2)

            !BEGIN: COMMENT TO DO CELL AREA TEST
            deltaLon = lonCell(jEdgeCells(2)) - lonCell(jEdgeCells(1))
            Bx = cos(latCell(jEdgeCells(2))) * cos(deltaLon)
            By = cos(latCell(jEdgeCells(2))) * sin(deltaLon)
            laEdge = atan2(sin(latCell(jEdgeCells(1))) + sin(latCell(jEdgeCells(2))) , sqrt( (cos(latCell(jEdgeCells(1))) + &
                     Bx) * (cos(latCell(jEdgeCells(1))) + Bx) + By * By ))
            loEdge = lonCell(jEdgeCells(1)) + atan2(By, cos(latCell(jEdgeCells(1))) + Bx)            
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !deltaLon = lonCell(iEdgeCells(2)) - lonCell(iEdgeCells(1))
            !Bx = cos(latCell(iEdgeCells(2))) * cos(deltaLon)
            !By = cos(latCell(iEdgeCells(2))) * sin(deltaLon)
            !laEdge = atan2(sin(latCell(iEdgeCells(1))) +
            !sin(latCell(iEdgeCells(2))) , sqrt( (cos(latCell(iEdgeCells(1))) +
            !Bx) * (cos(latCell(iEdgeCells(1))) + Bx) + By * By ))
            !loEdge = lonCell(iEdgeCells(1)) + atan2(By,
            !cos(latCell(iEdgeCells(1))) + Bx)            
            !END: UNCOMMENT TO DO CELL AREA TEST

            !BEGIN: COMMENT TO DO CELL AREA TEST
            areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell) / areaCell(sharedCell)
            !END: COMMENT TO DO CELL AREA TEST

            !BEGIN: UNCOMMENT TO DO CELL AREA TEST
            !areaSum = areaSum + sphericalTriangleArea(laEdge, loEdge, laVertex,
            !loVertex, laCell, loCell)
            !print*, 'areaSum =', areaSum, 'areaCell =', areaCell(sharedCell),
            !'difference =', abs(areaSum - areaCell(sharedCell)), 'edge =',
            !iEdge, 'cell =', sharedCell
            !END: UNCOMMENT TO DO CELL AREA TEST

            ! --- END: compute areaSum

            weightsOnEdge(i, iEdge) = c  * dvEdge(jEdge)/dcEdge(iEdge) * (0.5 - areaSum) * n 

         end do
      end do
  
   end subroutine compute_weightsOnEdgeFortran

   real function sphericalTriangleArea(laEdge, loEdge, laVertex, loVertex, laCell, loCell)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! sphericalTriangleArea uses the spherical analog of Heron's formula to
   ! compute
   ! the area of a triangle on the surface of a sphere
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      real (kind=RKIND), intent(in) :: laEdge, loEdge, laVertex, loVertex, laCell, loCell
      real (kind=RKIND) :: tanqe, s, a, b, c, radius
      radius = 1.0 

      a = sphere_distance(laCell, loCell, laVertex, loVertex, radius)
      b = sphere_distance(laCell, loCell, laEdge, loEdge, radius)
      c = sphere_distance(laEdge, loEdge, laVertex, loVertex, radius)
      s = 0.5*(a+b+c)

      tanqe = sqrt(tan(0.5*s)*tan(0.5*(s-a))*tan(0.5*(s-b))*tan(0.5*(s-c)))
      sphericalTriangleArea = 4.*atan(tanqe)

   end function sphericalTriangleArea 

   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
   !   sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance

   real function aa(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! A, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      aa = 0.5 * w * (2.0 * omega + w) * cos(theta)**2.0 + &
          0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 + 2.0*R**2.0 - R - 2.0 - 2.0*R**2.0 * cos(theta)**(-2.0))

   end function aa

   
   real function bb(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! B, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      bb = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * cos(theta)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*cos(theta))**2.0)

   end function bb


   real function cc(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      cc = 0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 - R - 2.0)

   end function cc

end module sw_test_cases 
