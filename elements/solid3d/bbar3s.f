c$Id: bbar3s.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine bbar3s(phi,shp,dvol,lint,nel,npm,hh,theta,bbar)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute mixed formulation for the volumetric response.
c              for small deformation problem.

c     Inputs:
c        phi(4,*)      - Stress functions
c        shp(4,27,*)   - Shape functions and derivatives
c        vol(*)        - Volume elements
c        lint          - Number of quadrature points
c        nel           - Number of nodes on element (should be 4 or 9)
c        npm           - Number of pressure modes
c        theta(3,*)    - Volumetric strain from displacements

c     Outputs:
c        hh(4,4)       - Volume/pressure shape integrals (inverse)
c        theta(3,*)    - Mixed volumetric strain
c        bbar(3,27,*)  - Mixed volumetric derivatives of shape functions.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   lint,   nel,  npm,   i,  j,  k, l
      real*8    shp(4,27,*),  dvol(*)   ,  ht(4,3),  h1
      real*8    gg(4,3,27) ,  hh(4,4)   ,  hv(4,3),  hg(4,3,16)
      real*8    phi(4,64)  ,  theta(3,*),  bbar(3,27,*)

      save

c     Constant pressure elements

      if(npm.eq.1) then

        do j = 1,nel
          do i = 1,3
            bbar(i,j,1) = 0.0d0
          end do ! i
        end do ! j
        hh(1,1) = 0.d0
        ht(1,1) = 0.0d0
        ht(1,2) = 0.0d0

        do l = 1,lint

c         H-array and G-array

          hh(1,1) = hh(1,1) + dvol(l)
          ht(1,1) = ht(1,1) + theta(1,l)*dvol(l)
          ht(1,2) = ht(1,2) + theta(2,l)*dvol(l)

c         G-array

          do j = 1,nel
            do i = 1,3
              bbar(i,j,1) = bbar(i,j,1) + shp(i,j,l) * dvol(l)
            end do ! i
          end do ! j
        end do ! l

c       Average Jacobian

        hh(1,1)    = 1.d0 / hh(1,1)

c       Small deformation case

        theta(1,1) = hh(1,1)*ht(1,1)
        theta(2,1) = hh(1,1)*ht(1,2)
        theta(3,1) = theta(1,1) - theta(2,1)

c       Modify bbar for B-bar type computations

        do j = 1,nel
          do i = 1,3
            bbar(i,j,1) = bbar(i,j,1)*hh(1,1)
          end do ! i
        end do ! j

c       Copy for other quadrature points

        do l = 2,lint
          theta(1,l) = theta(1,1)
          theta(2,l) = theta(2,1)
          theta(3,l) = theta(3,1)
          do j = 1,nel
            do i = 1,3
              bbar(i,j,l) = bbar(i,j,1)
            end do ! i
          end do ! j
        end do ! l

c     Higher order elements: npm = 4

      else

        do i = 1,npm
          do j = 1,nel
            gg(i,1,j) = 0.0d0
            gg(i,2,j) = 0.0d0
          end do ! j
          do j = 1,npm
            hh(i,j) = 0.0d0
          end do ! j
          ht(i,1) = 0.0d0
          ht(i,2) = 0.0d0
        end do ! i

c       Quadrature loop

        do l = 1,lint
          do j = 1,npm

c           H-array

            h1      = phi(j,l) * dvol(l)
            ht(j,1) = ht(j,1)  + theta(1,l)*h1
            ht(j,2) = ht(j,2)  + theta(2,l)*h1
            do i = 1,npm
              hh(i,j)   = hh(i,j)   + phi(i,l)*h1
            end do ! i

c           G-array

            do i = 1,nel
              gg(j,1,i) = gg(j,1,i) + shp(1,i,l)*h1
              gg(j,2,i) = gg(j,2,i) + shp(2,i,l)*h1
              gg(j,3,i) = gg(j,3,i) + shp(3,i,l)*h1
            end do ! i
          end do ! j

        end do ! l

c       Invert H-array

        call invert(hh,npm,6)

        do j = 1,nel
          do i = 1,npm
            hg(i,1,j) = 0.0d0
            hg(i,2,j) = 0.0d0
            hg(i,3,j) = 0.0d0
            do k = 1,npm
              hg(i,1,j) = hg(i,1,j) + hh(i,k)*gg(k,1,j)
              hg(i,2,j) = hg(i,2,j) + hh(i,k)*gg(k,2,j)
              hg(i,3,j) = hg(i,3,j) + hh(i,k)*gg(k,3,j)
            end do ! k
          end do ! i
        end do ! j

        do j = 1,2
          do i = 1,npm
            hv(i,j) = 0.0d0
            do k = 1,npm
              hv(i,j) = hv(i,j) + hh(i,k)*ht(k,j)
            end do ! k
          end do ! i
        end do ! j

        do l = 1,lint
          theta(1,l) = hv(1,1)
          theta(2,l) = hv(1,2)
          do k = 2,npm
            theta(1,l) = theta(1,l) + phi(k,l)*hv(k,1)
            theta(2,l) = theta(2,l) + phi(k,l)*hv(k,2)
          end do ! k
        end do ! l

        do l = 1,lint
          theta(3,l) = theta(1,l) - theta(2,l)
          do j = 1,nel
            do i = 1,3
              bbar(i,j,l) = hg(1,i,j)
            end do ! i
            do k = 2,npm
              do i = 1,3
                bbar(i,j,l) = bbar(i,j,l) + phi(k,l)*hg(k,i,j)
              end do ! i
            end do ! k
          end do ! j
        end do ! l

      endif

      end
