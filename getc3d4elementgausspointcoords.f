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
      subroutine getc3d4elementgausspointcoords(nelem, elem_ids,
     &    co, kon, ipkon, elem_gp_coord)

         implicit none

         !Input/Output variables
         integer     :: nelem            ! Number of elements
         integer     :: elem_ids(*)      ! Element IDs, need to be adjusted for Fortran
         real(8)     :: co(3,*)          ! Nodal coordinates of all nodes
         integer     :: kon(*)
         integer     :: ipkon(*)

         real(8)     :: elem_gp_coord(*) ! Gauss point coords

         ! Internal variables
         integer(4)     :: iel, el_id, indexe, i, j, k, l
         integer(4)     :: iflag, idx, konl(4), gp_id
         REAL(8)        :: xl(3,4), xi, et, ze, xsj, shp(4,20)

         include "gauss.f"

         data iflag /1/

         ! Initialize gauss point coordinates to zero
         do l = 1, nelem*3
            elem_gp_coord(l) = 0.0
         end do

         do iel = 1, nelem
            el_id = elem_ids(iel)

            indexe=ipkon(el_id)

            ! Local nodal coordinates
            do i = 1, 4
               konl(i) = kon(indexe+i)
               do j = 1, 3
                  xl(j,i)=co(j,konl(i))
               end do
            end do

            ! One Gauss point per element
            gp_id = 1

            ! gauss3d4: tet, 1 integration point
            xi = gauss3d4(1,gp_id)
            et = gauss3d4(2,gp_id)
            ze = gauss3d4(3,gp_id)

            ! Get the shape function
            call shape4tet(xi,et,ze,xl,xsj,shp,iflag)

            ! Calculate the Gauss point coordinates
            do k = 1, 3
               idx = (iel - 1) * 3 + k
               do l = 1, 4
                  elem_gp_coord(idx)=elem_gp_coord(idx)
     &                  +xl(k,l)*shp(4,l)
               end do
            end do

         end do

         return

      end subroutine getc3d4elementgausspointcoords
