c$Id: tint3d.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine tint3d(ll,lint,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Gauss quadrature for 3-d tetrahedral element

c      Inputs:
c         ll       - Type of quadrature

c      Outputs:
c         lint     - Number of quadrature points
c         s(5,*)   - Values of volume coordinates and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i, j,k,l,m, ll, lint
      real*8    s(5,*), alp(2)

      integer  ind1(6),ind2(6),ind3(6),ind4(6)

      save

      data     ind1/ 1,2,2,1,1,2/
      data     ind2/ 1,1,2,2,2,1/
      data     ind3/ 2,1,1,2,1,2/
      data     ind4/ 2,2,1,1,2,1/

      data     alp /0.399403576166799219d+00,0.100596423833200785d+00/

c     1 pt. quadrature O(h^2)

      if(ll.eq.1) then
        lint = 1
        do i = 1,4
          s(i,1) = 0.25d0
        end do ! i
        s(5,1) = 1.0d0/6.d0

c     4 pt. quadrature O(h^2) - nodes of linear element

      elseif(ll.eq.-1) then
        lint = 4
        s(5,4) = 0.25d0/6.d0
        do i = 1,4
          do j = 1,4
            s(i,j) = 0.0d0
          end do ! j
          s(i,i) = 1.d0
          s(5,i) = s(5,4)
        end do ! i

c     4 pt. quadrature O(h^3)

      elseif(ll.eq.2) then
        lint = 4
        s(5,4) = 0.25d0/6.d0
        do i = 1,4
          do j = 1,4
            s(i,j) = 0.1381966011250105d+00
          end do ! j
          s(i,i) = 0.5854101966249658d+00
          s(5,i) = s(5,4)
        end do ! i

c     5  pt. quadrature O(h^4) -- has negative weight

      elseif(ll.eq.3) then
        lint = 5
        do i = 1,4
          s(i,1) = 1.d0/6.0d0
          s(i,2) = s(i,1)
          s(i,3) = s(i,1)
          s(i,4) = s(i,1)
          s(i,5) = 0.25d0
        end do ! i
        do i = 1,4
          s(i,i) = 0.5d0
        end do ! i
        s(5,1) =  0.45d0
        s(5,2) =  s(5,1)
        s(5,3) =  s(5,1)
        s(5,4) =  s(5,1)
        s(5,5) = -0.80d0

c     11 pt. quadrature O(h^4) -- has negative weight

      elseif(ll.eq.4) then
        lint = 11
        do i = 1,4
          do j = 1,4
            s(j,i) = 0.714285714285714285d-01
          end do ! j
        end do ! i
        do i = 1,4
          s(i,i) = 0.785714285714285714d+00
          s(5,i) = 0.762222222222222222d-02
        end do ! i

        do i = 5,10
          j      = ind1(i-4)
          k      = ind2(i-4)
          l      = ind3(i-4)
          m      = ind4(i-4)
          s(1,i) = alp(j)
          s(2,i) = alp(k)
          s(3,i) = alp(l)
          s(4,i) = alp(m)
          s(5,i) = 0.248888888888888880d-01
        end do ! i
        s(1,11) =  0.25d+00
        s(2,11) =  0.25d+00
        s(3,11) =  0.25d+00
        s(4,11) =  0.25d+00
        s(5,11) = -0.131555555555555550d-01
         
c     11 pt. quadrature O(h^4) -- has no negative weight

      elseif(ll.eq.-4) then
        lint = 11
        do i = 1,4
          do j = 1,10
            s(i,j) = 0.0d0
          end do ! j
          s(i, i  ) = 1.00d0
          s(i, i+4) = 0.50d0
          s(i, i+7) = 0.50d0
          s(i, 11 ) = 0.25d0
        end do ! i
        s(2, 5) = 0.50d0
        s(3, 6) = 0.50d0
        s(1, 7) = 0.50d0
        do j = 1,4
          s(5,j) = 1.d0/360.d0
        end do ! j
        do j = 5,10
          s(5,j) = 1.d0/90.d0
        end do ! j
        s(5,11) = 4.d0/45.d0

c     16 pt. quadrature O(h^5)

      else
        lint = 16
        s(5,4) = 0.8395632516687135d-02
        do i = 1,3
          do j = 1,4
            s(i,j) = 0.7611903264425430d-01
          end do ! j
          s(i,i) = 0.7716429020672371d+00
          s(5,i) = s(5,4)
        end do ! i
        do i = 5,16
          s(5,i) = 0.1109034477221540d-01
        end do ! i

        s(1, 5) = 0.1197005277978019d+00
        s(2, 5) = 0.7183164526766925d-01
        s(3, 5) = 0.4042339134672644d+00
        s(1, 6) = 0.4042339134672644d+00
        s(2, 6) = 0.1197005277978019d+00
        s(3, 6) = 0.7183164526766925d-01
        s(1, 7) = 0.4042339134672644d+00
        s(2, 7) = 0.4042339134672644d+00
        s(3, 7) = 0.1197005277978019d+00
        s(1, 8) = 0.7183164526766925d-01
        s(2, 8) = 0.4042339134672644d+00
        s(3, 8) = 0.4042339134672644d+00

        s(1, 9) = 0.1197005277978019d+00
        s(2, 9) = 0.4042339134672644d+00
        s(3, 9) = 0.7183164526766925d-01
        s(1,10) = 0.4042339134672644d+00
        s(2,10) = 0.1197005277978019d+00
        s(3,10) = 0.4042339134672644d+00
        s(1,11) = 0.7183164526766925d-01
        s(2,11) = 0.4042339134672644d+00
        s(3,11) = 0.1197005277978019d+00
        s(1,12) = 0.4042339134672644d+00
        s(2,12) = 0.7183164526766925d-01
        s(3,12) = 0.4042339134672644d+00

        s(1,13) = 0.1197005277978019d+00
        s(2,13) = 0.4042339134672644d+00
        s(3,13) = 0.4042339134672644d+00
        s(1,14) = 0.7183164526766925d-01
        s(2,14) = 0.1197005277978019d+00
        s(3,14) = 0.4042339134672644d+00
        s(1,15) = 0.4042339134672644d+00
        s(2,15) = 0.7183164526766925d-01
        s(3,15) = 0.1197005277978019d+00
        s(1,16) = 0.4042339134672644d+00
        s(2,16) = 0.4042339134672644d+00
        s(3,16) = 0.7183164526766925d-01

      endif

c     Compute fourth points

      do j = 1,lint
        s(4,j) = 1.d0 - (s(1,j) + s(2,j) + s(3,j))
      end do ! j

      end
