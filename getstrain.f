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
      subroutine getstrain(co,ntmat_,vold,
        &  cocon,ncocon,iset,istartset,iendset,ipkon,lakon,kon,
        &  ialset,ielmat,mi,kdelta,reftemp)
   !
   !     returns the strain vector for given set of elements and corresponding gauss points. 
   !     Each set consists of the following: [element_id, gp_id, mat_id, eps_11, eps_22, eps_33, eps_12, eps_13, eps_23]
         implicit none
   !
!          character*8 lakonl,lakon(*)
!    !
!          integer konl(20),ifaceq(8,6),nelem,ii,i,j,i1,i2,j1,
!         &  ncocon(2,*),k1,jj,ig,ntmat_,nope,nopes,imat,
!         &  mint2d,ifacet(6,4),ifacew(8,5),iflag,indexe,jface,istartset(*),
!         &  iendset(*),ipkon(*),kon(*),iset,ialset(*),nset,
!         &  mi(*),ielmat(mi(3),*),fidx
!    !
!          real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),dvi,f(0:3),
!         &  vkl(0:3,3),t(3,3),div,
!         &  voldl(0:mi(2),20),cocon(0:6,ntmat_,*),xl2(3,8),xsj2(3),
!         &  shp2(7,8),vold(0:mi(2),*),xi,et,xsj,temp,xi3d,et3d,ze3d,weight,
!         &  xlocal20(3,9,6),xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),
!         &  xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),pres,
!         &  tf(0:3),tn,tt,dd,coords(3),cond, kdelta(*), reftemp(*), 
!         &  avgkdelta, avgreftemp
!    !
!          include "gauss.f"
!          include "xlocal.f"
!    !
!          data ifaceq /4,3,2,1,11,10,9,12,
!         &            5,6,7,8,13,14,15,16,
!         &            1,2,6,5,9,18,13,17,
!         &            2,3,7,6,10,19,14,18,
!         &            3,4,8,7,11,20,15,19,
!         &            4,1,5,8,12,17,16,20/
!          data ifacet /1,3,2,7,6,5,
!         &             1,2,4,5,9,8,
!         &             2,3,4,6,10,9,
!         &             1,4,3,8,10,7/
!          data ifacew /1,3,2,9,8,7,0,0,
!         &             4,5,6,10,11,12,0,0,
!         &             1,2,5,4,7,14,10,13,
!         &             2,3,6,5,8,15,11,14,
!         &             4,6,3,1,12,15,9,13/
!          data iflag /3/
!    !
!    !
!    !
!    !           initialisierung of the flux
!    !     
!                f(0)=0.d0
               
!    !			index for the output array
!                fidx=1
               
!    !     
!    !
!                do jj=istartset(iset),iendset(iset)
!    !     
!                   jface=ialset(jj)
!    !     
!                   nelem=int(jface/10.d0)
!                   ig=jface-10*nelem
!                   lakonl=lakon(nelem)
!                   indexe=ipkon(nelem)
!                   imat=ielmat(1,nelem) ! what is this?
!    !     
!                   if(lakonl(4:4).eq.'2') then
!                      nope=20
!                      nopes=8
!                   elseif(lakonl(4:4).eq.'8') then
!                      nope=8
!                      nopes=4
!                   elseif(lakonl(4:5).eq.'10') then
!                      nope=10
!                      nopes=6
!                   elseif(lakonl(4:4).eq.'4') then
!                      nope=4
!                      nopes=3
!                   elseif(lakonl(4:5).eq.'15') then
!                      nope=15
!                   elseif(lakonl(4:4).eq.'6') then
!                      nope=6
!                   endif
!    !     
!                   if(lakonl(4:5).eq.'8R') then
!                      mint2d=1
!                   elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
!         &             then
!                      if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
!         &                 (lakonl(7:7).eq.'E')) then
!                         mint2d=2
!                      else
!                         mint2d=4
!                      endif
!                   elseif(lakonl(4:4).eq.'2') then
!                      mint2d=9
!                   elseif(lakonl(4:5).eq.'10') then
!                      mint2d=3
!                   elseif(lakonl(4:4).eq.'4') then
!                      mint2d=1
!                   endif
!    !     
!    !     local topology
!    !     
!                   do i=1,nope
!                      konl(i)=kon(indexe+i)
!                   enddo
!    !     
!    !     computation of the coordinates of the local nodes
!    !     
!                   do i=1,nope
!                      do j=1,3
!                         xl(j,i)=co(j,konl(i))
!                      enddo
!                   enddo
!    !     
!    !     temperature, velocity and auxiliary variables
!    !     (rho*energy density, rho*velocity and rho)
!    !     
!                   do i1=1,nope
!                      do i2=0,mi(2)
!                         voldl(i2,i1)=vold(i2,konl(i1))
!                      enddo
!                   enddo
!    !     
!    !     treatment of wedge faces
!    !     
!                   if(lakonl(4:4).eq.'6') then
!                      mint2d=1
!                      if(ig.le.2) then
!                         nopes=3
!                      else
!                         nopes=4
!                      endif
!                   endif
!                   if(lakonl(4:5).eq.'15') then
!                      if(ig.le.2) then
!                         mint2d=3
!                         nopes=6
!                      else
!                         mint2d=4
!                         nopes=8
!                      endif
!                   endif
!    !     
!                   if((nope.eq.20).or.(nope.eq.8)) then
!                      do i=1,nopes
!                         do j=1,3
!                            xl2(j,i)=co(j,konl(ifaceq(i,ig)))
!                         enddo
!                      enddo
!                   elseif((nope.eq.10).or.(nope.eq.4)) then
!                      do i=1,nopes
!                         do j=1,3
!                            xl2(j,i)=co(j,konl(ifacet(i,ig)))
!                         enddo
!                      enddo
!                   else
!                      do i=1,nopes
!                         do j=1,3
!                            xl2(j,i)=co(j,konl(ifacew(i,ig)))
!                         enddo
!                      enddo
!                   endif
!    !     
   
   
   
!    !			very simple averaging over the integration points!!
!    !			TODO: check whether this needs to be improved!
!                   avgkdelta=0
!                   avgreftemp=0
   
   
!                   do i=1,mint2d
!    !     
!    !     local coordinates of the surface integration
!    !     point within the surface local coordinate system
!    !     
!                      if((lakonl(4:5).eq.'8R').or.
!         &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
!                         xi=gauss2d1(1,i)
!                         et=gauss2d1(2,i)
!                         weight=weight2d1(i)
!                      elseif((lakonl(4:4).eq.'8').or.
!         &                    (lakonl(4:6).eq.'20R').or.
!         &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
!                         xi=gauss2d2(1,i)
!                         et=gauss2d2(2,i)
!                         weight=weight2d2(i)
!                      elseif(lakonl(4:4).eq.'2') then
!                         xi=gauss2d3(1,i)
!                         et=gauss2d3(2,i)
!                         weight=weight2d3(i)
!                      elseif((lakonl(4:5).eq.'10').or.
!         &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
!                         xi=gauss2d5(1,i)
!                         et=gauss2d5(2,i)
!                         weight=weight2d5(i)
!                      elseif((lakonl(4:4).eq.'4').or.
!         &                    ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
!                         xi=gauss2d4(1,i)
!                         et=gauss2d4(2,i)
!                         weight=weight2d4(i)
!                      endif
!    !     
!    !     local surface normal
!    !     
!                      if(nopes.eq.8) then
!                         call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
!                      elseif(nopes.eq.4) then
!                         call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
!                      elseif(nopes.eq.6) then
!                         call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
!                      else
!                         call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
!                      endif
!    !
!    !                 global coordinates of the integration point
!    !
!                      do j1=1,3
!                         coords(j1)=0.d0
!                         do i1=1,nopes
!                            coords(j1)=coords(j1)+shp2(4,i1)*xl2(j1,i1)
!                         enddo
!                      enddo
!    !     
!    !     local coordinates of the surface integration
!    !     point within the element local coordinate system
!    !     
!                      if(lakonl(4:5).eq.'8R') then
!                         xi3d=xlocal8r(1,i,ig)
!                         et3d=xlocal8r(2,i,ig)
!                         ze3d=xlocal8r(3,i,ig)
!                         call shape8h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:4).eq.'8') then
!                         xi3d=xlocal8(1,i,ig)
!                         et3d=xlocal8(2,i,ig)
!                         ze3d=xlocal8(3,i,ig)
!                         call shape8h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:6).eq.'20R') then
!                         xi3d=xlocal8(1,i,ig)
!                         et3d=xlocal8(2,i,ig)
!                         ze3d=xlocal8(3,i,ig)
!                         call shape20h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:4).eq.'2') then
!                         xi3d=xlocal20(1,i,ig)
!                         et3d=xlocal20(2,i,ig)
!                         ze3d=xlocal20(3,i,ig)
!                         call shape20h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:5).eq.'10') then
!                         xi3d=xlocal10(1,i,ig)
!                         et3d=xlocal10(2,i,ig)
!                         ze3d=xlocal10(3,i,ig)
!                         call shape10tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:4).eq.'4') then
!                         xi3d=xlocal4(1,i,ig)
!                         et3d=xlocal4(2,i,ig)
!                         ze3d=xlocal4(3,i,ig)
!                         call shape4tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:5).eq.'15') then
!                         xi3d=xlocal15(1,i,ig)
!                         et3d=xlocal15(2,i,ig)
!                         ze3d=xlocal15(3,i,ig)
!                         call shape15w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      elseif(lakonl(4:4).eq.'6') then
!                         xi3d=xlocal6(1,i,ig)
!                         et3d=xlocal6(2,i,ig)
!                         ze3d=xlocal6(3,i,ig)
!                         call shape6w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
!                      endif
!    !     
!    !     calculating of
!    !     the temperature temp
!    !     in the integration point
!    !     
   
!                      temp=0.d0
!                      do j1=1,3
!                         vkl(0,j1)=0.d0
!                      enddo
!                      do i1=1,nope
!                         temp=temp+shp(4,i1)*voldl(0,i1)
!                         do k1=1,3
!                            vkl(0,k1)=vkl(0,k1)+shp(k1,i1)*voldl(0,i1)
!                         enddo
!                      enddo
!    !     
!    !     material data (conductivity)
!    !     
!                      call materialdata_cond(imat,ntmat_,temp,cocon,
!         &                    ncocon,cond)
!    !     
!    !     determining the stress 
!    !     
!                      dd=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
!         &                    xsj2(3)*xsj2(3))
!    !
!                      tf(0)=-cond*(vkl(0,1)*xsj2(1)+
!         &                            vkl(0,2)*xsj2(2)+
!         &                            vkl(0,3)*xsj2(3))
!                      f(0)=f(0)+tf(0)*weight
!                      tf(0)=tf(0)/dd
!    !                  dd = .1;
!    !				  print *, temp, tf(0), dd, temp+tf(0)/cond
!                      avgreftemp = avgreftemp + (temp + tf(0) / cond * dd)
!                      avgkdelta = avgkdelta + cond / dd
!    !     
!    !                  print *, fidx, avgreftemp, avgkdelta, tf(0), cond
!                   enddo
                  
!                   reftemp(fidx) = avgreftemp / mint2d
!                   kdelta(fidx) = avgkdelta / mint2d
                  
!                   fidx=fidx+1
                  
!                enddo
   !
   !     
         return
         end
         
         
   
   