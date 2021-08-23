c$Id: presld.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine presld(d,ul,xl,ix,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  Nodal force and tangent array for pressure loading

c      INPUT variables
c        d(60)      Value of constant pressure on face
c        xl(4,*)    Nodal coordinates
c        ul(4,*)    Nodal displacements
c        ndf        Number of DOF / node
c        ndm        Space dimension
c        nst        Dimension of residual vector

c      OUTPUT variables
c        p(ndf,4)   Contribution to residual
c        s(nst,nst) Contribution to striffness matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_tanfl.h'
      include  'eldata.h'
      include  'eltran.h'
      include  'iofile.h'
      include  'pglob1.h'
      include  'prld1.h'

      logical   pcomp, errck, tinput, mread, prsfl
      character texti*15
      integer   l, lint, ndf, ndm, nst, isw, ii, jj, kk, i1, j1, ix(*)
      real*8    pn, pp, x0,v0, xx, rr
      real*8    xii(4), eti(4), td(7), xu(3,4),dx(3,2)
      real*8    sg2(2,5), sg3(3,4), shp2(2,4), shp3(3,4)
      real*8    d(*), xl(ndm,*), ul(ndf,*), p(ndf,*), s(nst,*)

      save

      data      xii / -0.5d0, 0.5d0, 0.5d0,-0.5d0 /
      data      eti / -0.5d0,-0.5d0, 0.5d0, 0.5d0 /

      if(isw.eq.1) then

c       Axisymmetric type indicator

        d(81) = g2type

c       Deactivate dof in element for dof > ndm

        do ii = ndm+1,ndf
          ix(ii) = 0
        end do ! ii

        mread = .true.
        do while (mread)
          errck = tinput(texti,1,td,7)
          if(pcomp(texti,'load',4) .or. pcomp(texti,'pres',4)) then
            prsfl = .true.
            d(60) = td(1)
            d(61) = nint(td(2))
          elseif(pcomp(texti,'grad',4) .or. pcomp(texti,'bend',4)) then
            d(62) = td(1)
            d(63) = nint(td(2))
          elseif(pcomp(texti,'hydr',4)) then
            prsfl = .false.
            d(64) = td(1) ! j   = Direction of pressure gradient
            d(65) = td(2) ! x_j = Zero pressure elevation
            d(66) = td(3) ! v_j = Elevation velocity
            d(67) = td(4) ! rho = Density
            d(68) = td(5) ! g_j = Gravity value
            d(69) = td(6) ! p_n = Proportional load number pressure
            d(70) = td(7) ! p_h = Proportional load number height
          elseif(pcomp(texti,'axis',4)) then
            d(81) = 3.0d0
          elseif(pcomp(texti,'plan',4)) then
            d(81) = 0.0d0
          elseif(pcomp(texti,'fini',4) .or. pcomp(texti,'foll',4)) then
            d(82) = 1.d0
          else
            mread = .false.
          endif
        end do ! while

c       Output definitions

        if(prsfl) then
          write(iow,2000) d(60),nint(d(61)),d(62),nint(d(63))
          if(d(82).gt.0.0d0) then
            write(iow,2001)
          endif
          d(83) = 0.0d0
        else
          write(iow,2002) nint(d(64)),d(65),d(66),d(67),d(68),
     &                    nint(d(69)),nint(d(70))
          d(83) = 1.0d0
        endif

c       Set plot for point

        call pltpt1(iel)

      elseif(isw.eq.3 .or. (isw.eq.6 .and. .not.ddfl)) then

c       Compute nodal coordinates in correct reference frame

        if(d(82).eq.0.0d0) then ! Fixed force
          do ii = 1,nel
            do jj = 1,ndm
              xu(jj,ii) = xl(jj,ii)
            end do ! jj
          end do ! ii
        else                    ! Follower force
          do ii = 1,nel
            do jj = 1,ndm
              xu(jj,ii) = xl(jj,ii) + ul(jj,ii)
            end do ! ii
          end do ! ii
        endif

c       2-D Problems

        if(ndm.eq.2) then

c         Get quadrature information

          lint = nel
          if(nint(d(182)).gt.0) then
            call int1d (lint, sg2)
          else
            call int1d(lint, sg2)
          endif

c         Loop over quadrature points

          do l = 1,lint

c           Compute shape functions and geometric factors

            call shap1d( sg2(1,l), nel, shp2 )

            do ii = 1,ndm
              dx(ii,1) = shp2(1,1)*xu(ii,1)
              do jj = 2,nel
                dx(ii,1) = dx(ii,1) + shp2(1,jj)*xu(ii,jj)
              end do ! jj
            end do ! ii

c           Axisymmetric problem

            if(nint(d(81)).eq.3) then
              rr = 0.0d0
              do jj = 1,nel
                rr = rr + shp2(2,jj)*xl(1,jj)
              end do ! jj
            else
              rr = 1.0d0
            endif

c           Compute nodal loads for pressures

            if(nint(d(83)).eq.0) then
              kk = nint(d(61))
              if(kk.eq.0) then
                pn = dm*d(60)*sg2(2,l)
              else
                pn = prldv(kk)*d(60)*sg2(2,l)
              endif
              kk = nint(d(63))
              if(kk.eq.0) then
                pn = pn + dm*d(62)*sg2(1,l)*sg2(2,l)
              else
                pn = pn + prldv(kk)*d(62)*sg2(1,l)*sg2(2,l)
              endif
            elseif(nint(d(83)).eq.1) then
              ii = nint(d(64))
              x0 = d(65)
              v0 = d(66)
              kk = nint(d(70))
              if(kk.eq.0) then
                x0 = x0 + v0*dm
              else
                x0 = x0 + v0*prldv(kk)
              endif
              xx = 0.0
              do jj = 1,nel
                xx = xx + shp2(2,jj)*xl(ii,jj)
              end do ! jj
              if(xx.le.x0) then
                pn = d(67)*d(68)*(xx - x0)*sg2(2,l)
                kk = nint(d(69))
                if(kk.eq.0) then
                  pn = pn*dm
                else
                  pn = pn*prldv(kk)
                endif
              else
                pn = 0.0d0
              endif
            endif
            pn = pn*rr                   ! Axisymmetric correction
            do ii = 1,nel
              pp      = shp2(2,ii)*pn
              p(1,ii) = p(1,ii) + pp*dx(2,1)
              p(2,ii) = p(2,ii) - pp*dx(1,1)
            end do ! ii

c           Compute tangent if necessary

            if(d(82).gt.0.0d0) then
              i1 = 0
              do ii = 1,nel
                pp = shp2(2,ii)*pn*ctan(1)
                j1 = 0
                do jj = 1,nel
                  s(i1+1,j1+2) = s(i1+1,j1+2) - pp*shp2(1,jj)
                  s(i1+2,j1+1) = s(i1+2,j1+1) + pp*shp2(1,jj)
                  j1 = j1 + ndf
                end do ! jj
                i1 = i1 + ndf
              end do ! ii
            endif

          end do ! l

c       3-D Problems

        elseif(ndm.eq.3) then

c         Get quadrature information

          l  = 2
          call int2d (l, lint, sg3)

c         Loop over quadrature points

          do l = 1,lint

c           Compute shape functions and geometric factors

            do ii = 1,nel
              pp         = xii(ii)*sg3(1,l) + 0.5d0
              pn         = eti(ii)*sg3(2,l) + 0.5d0
              shp3(1,ii) = xii(ii)*pn
              shp3(2,ii) = eti(ii)*pp
              shp3(3,ii) = pn*pp
            end do ! ii

            do ii = 1,3
              dx(ii,1) = shp3(1,1)*xu(ii,1) + shp3(1,2)*xu(ii,2)
     &                 + shp3(1,3)*xu(ii,3) + shp3(1,4)*xu(ii,4)
              dx(ii,2) = shp3(2,1)*xu(ii,1) + shp3(2,2)*xu(ii,2)
     &                 + shp3(2,3)*xu(ii,3) + shp3(2,4)*xu(ii,4)
            end do ! ii

c           Compute nodal loads for pressures

            if(nint(d(83)).eq.0) then
              kk = nint(d(61))
              pn = d(60)*sg3(3,l)
              if(kk.eq.0) then
                pn = pn*dm
              elseif(kk.gt.0) then
                pn = pn*prldv(kk)
              endif
            elseif(nint(d(83)).eq.1) then
              ii = nint(d(64))
              x0 = d(65)
              v0 = d(66)
              kk = nint(d(70))
              if(kk.eq.0) then
                x0 = x0 + v0*dm
              else
                x0 = x0 + v0*prldv(kk)
              endif
              xx = 0.0
              do jj = 1,nel
                xx = xx + shp2(2,jj)*xl(ii,jj)
              end do ! jj
              if(xx.le.x0) then
                pn = d(67)*d(68)*(xx - x0)*sg3(3,l)
                kk = nint(d(69))
                if(kk.eq.0) then
                  pn = pn*dm
                else
                  pn = pn*prldv(kk)
                endif
              else
                pn = 0.0d0
              endif
            endif
            do ii = 1,4
              pp     = shp3(3,ii)*pn
              p(1,ii) = p(1,ii) + pp*(dx(2,1)*dx(3,2) - dx(3,1)*dx(2,2))
              p(2,ii) = p(2,ii) + pp*(dx(3,1)*dx(1,2) - dx(1,1)*dx(3,2))
              p(3,ii) = p(3,ii) + pp*(dx(1,1)*dx(2,2) - dx(2,1)*dx(1,2))
            end do ! ii

c           Compute tangent if necessary

            if(d(82).gt.0.0d0) then
              i1 = 0
              do ii = 1,4
                pp = shp3(3,ii)*pn*ctan(1)
                j1 = 0
                do jj = 1,4
                  s(i1+1,j1+2) = s(i1+1,j1+2) - pp*(shp3(1,jj)*dx(3,2)
     &                                        -     shp3(2,jj)*dx(3,1))
                  s(i1+2,j1+1) = s(i1+2,j1+1) + pp*(shp3(1,jj)*dx(3,2)
     &                                        -     shp3(2,jj)*dx(3,1))

                  s(i1+2,j1+3) = s(i1+2,j1+3) - pp*(shp3(1,jj)*dx(1,2)
     &                                        -     shp3(2,jj)*dx(1,1))
                  s(i1+3,j1+2) = s(i1+3,j1+2) + pp*(shp3(1,jj)*dx(1,2)
     &                                        -     shp3(2,jj)*dx(1,1))

                  s(i1+3,j1+1) = s(i1+3,j1+1) - pp*(shp3(1,jj)*dx(2,2)
     &                                        -     shp3(2,jj)*dx(2,1))
                  s(i1+1,j1+3) = s(i1+1,j1+3) + pp*(shp3(1,jj)*dx(2,2)
     &                                        -     shp3(2,jj)*dx(2,1))
                  j1 = j1 + ndf
                end do ! jj
                i1 = i1 + ndf
              end do ! ii
            endif

          end do ! l

        endif ! ndm test

c     External node check

      elseif(isw.eq.26) then

      endif ! isw test

c     Formats

2000  format(/'  P r e s s u r e   L o a d i n g'//
     &       10x,'Loading Intensity    ',1p,1e12.4/
     &       10x,'Proportional Load No.',i12/
     &       10x,'Gradient Intensity   ',1p,1e12.4/
     &       10x,'Proportional Load No.',i12/1x)

2001  format(10x,'Follower Loading'/1x)

2002  format(/'  P r e s s u r e   L o a d i n g'//
     &       10x,'Loading direction         =',i12/
     &       10x,'Zero pressure elevation   =',1p,1e12.4/
     &       10x,'Velocity of surface       =',1p,1e12.4/
     &       10x,'Density/specific weight   =',1p,1e12.4/
     &       10x,'Gravity acceleration      =',1p,1e12.4/
     &       10x,'Proportional Load Pressure=',i12/
     &       10x,'Proportional Load Height  =',i12/1x)

      end
