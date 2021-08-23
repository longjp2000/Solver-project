c$Id: tint3dn.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine tint3dn(ll,lint,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Nodal quadrature for 3-d tetrahedral element

c      Inputs:
c         ll       - Number of nodes on element

c      Outputs:
c         lint     - Number of quadrature points
c         s(5,*)   - Values of volume coordinates and weights
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   i, j, ll, lint
      real*8    s(5,*)

      save

c     4 pt. quadrature O(h^2) - nodes of linear element

      if(ll.eq.4) then
        lint = 4
        s(5,4) = 0.25d0/6.d0
        do i = 1,4
          do j = 1,4
            s(i,j) = 0.0d0
          end do ! j
          s(i,i) = 1.d0
          s(5,i) = s(5,4)
        end do ! i

c     11 pt. quadrature O(h^4) -- has no negative weight

      elseif(ll.eq.11) then
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

c     Error: No quadrature available

      else

        write(iow,2000) ll
        write(ilg,2000) ll
        if(ior.lt.0) then
          write(*,2000) ll
        endif
        call plstop()

      endif

c     Compute fourth points

      do j = 1,lint
        s(4,j) = 1.d0 - (s(1,j) + s(2,j) + s(3,j))
      end do ! j

c     Format

2000  format(' *ERROR* TINT3dN: No quadrature available - nel =',i3)

      end
