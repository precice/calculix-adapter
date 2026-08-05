!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!     This function is a copy of the function for calculation and printout of the lift and drag forces
!
!     Computes the physical coordinates of the Gauss/integration points
!     of a set of elements. Supports either C3D4 (tetrahedra, 1 GP) or
!     C3D8 (hexahedra, 8 GP) elements, selected via the "nope" argument
!     (number of nodes per element: 4 or 8).
!
      subroutine getelementgausspointcoords(nelem, elem_ids, nope,
     &    co, kon, ipkon, elem_gp_coord)

         implicit none

         !Input/Output variables
         integer     :: nelem            ! Number of elements
         integer     :: elem_ids(*)      ! Element IDs, need to be adjusted for Fortran
         integer     :: nope             ! Nodes per element: 4 (C3D4) or 8 (C3D8)
         real(8)     :: co(3,*)          ! Nodal coordinates of all nodes
         integer     :: kon(*)
         integer     :: ipkon(*)

         real(8)     :: elem_gp_coord(*) ! GP coords

         ! Internal variables
         integer(4)     :: iel, el_id, indexe, i, j, k, l
         integer(4)     :: iflag, idx, konl(8), gp_id, ngp
         real(8)        :: xl(3,8), xi, et, ze, xsj, shp(4,20)

         include "gauss.f"

         data iflag /1/

         ! Number of Gauss points per element: 1 for tetrahedra, 8 for hexahedra
         if (nope .eq. 4) then
            ngp = 1
         else
            ngp = 8
         end if

         ! Initialize gauss point coordinates to zero
         do l = 1, nelem*ngp*3
            elem_gp_coord(l) = 0.0
         end do

         do iel = 1, nelem
            el_id = elem_ids(iel)

            indexe=ipkon(el_id)

            ! connectivity + local nodal coordinates
            do i = 1, nope
               konl(i)=kon(indexe+i)
               do j = 1, 3
                  xl(j,i)=co(j,konl(i))
               end do
            end do

            ! Loop through gauss points of each element
            do gp_id = 1, ngp

               if (nope .eq. 4) then
                  ! gauss3d4: tet, 1 integration point
                  xi = gauss3d4(1,gp_id)
                  et = gauss3d4(2,gp_id)
                  ze = gauss3d4(3,gp_id)

                  call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
               else
                  ! gauss3d2: hex, 2-point integration (8 integration points)
                  xi = gauss3d2(1,gp_id)
                  et = gauss3d2(2,gp_id)
                  ze = gauss3d2(3,gp_id)

                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               end if

               ! Calculate the Gauss point coordinates
               ! Note: indexed by the loop counter (iel), not the element ID
               ! (el_id), since elem_gp_coord is a compact buffer holding only
               ! the "nelem" elements of this interface, in the order given by
               ! elem_ids.
               do k = 1, 3
                  idx = (iel - 1)*ngp*3 + (gp_id - 1)*3 + k
                  do l = 1, nope
                     elem_gp_coord(idx)=elem_gp_coord(idx)
     &                     +xl(k,l)*shp(4,l)
                  end do
               end do

            end do

         end do

         return

      end subroutine getelementgausspointcoords
