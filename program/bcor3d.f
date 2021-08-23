c$Id: bcor3d.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine bcor3d(ixl,xl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute missing coordinate values of 27-node element

c      Inputs:
c         ixl(*)    - Nodal connection list
c         xl(3,*)   - Unadjusted coordinate array

c      Outputs:
c         xl(3,*)   - Adjusted coordinate array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    ixl(27), imid(12),amid(12),bmid(12)
      real*8     xl(3,27)

      integer    i,j

      save

      data       imid/9,10,11,12, 13,14,15,16, 18,19,20,21/
      data       amid/1, 2, 3, 4,  1, 2, 3, 4,  5, 6, 7, 8/
      data       bmid/5, 6, 7, 8,  2, 3, 4, 1,  6, 7, 8, 5/

c     Mid edge coordinates

      do i = 1,12
        if(ixl(imid(i)).eq.0) then
          do j = 1,3
            xl(j,imid(i)) = 0.5d0*(xl(j,amid(i)) + xl(j,bmid(i)))
          end do ! j
          ixl(i) = i
        endif
      end do ! i

c     Bottom and top

      if(ixl(17).eq.0) then
        do j = 1,3
          xl(j,17) = 0.50d0*(xl(j,13) + xl(j,14) + xl(j,15) + xl(j,16))
     &             - 0.25d0*(xl(j, 1) + xl(j, 2) + xl(j, 3) + xl(j, 4))
        end do ! j
        ixl(17) = 17
      endif

      if(ixl(22).eq.0) then
        do j = 1,3
          xl(j,22) = 0.50d0*(xl(j,18) + xl(j,19) + xl(j,20) + xl(j,21))
     &             - 0.25d0*(xl(j, 5) + xl(j, 6) + xl(j, 7) + xl(j, 8))
        end do ! j
        ixl(22) = 22
      endif

c     Mid-face

      if(ixl(23).eq.0) then
        do j = 1,3
          xl(j,23) = 0.50d0*(xl(j,13) + xl(j, 9) + xl(j,10) + xl(j,18))
     &             - 0.25d0*(xl(j, 1) + xl(j, 2) + xl(j, 5) + xl(j, 6))
        end do ! j
        ixl(23) = 23
      endif

      if(ixl(24).eq.0) then
        do j = 1,3
          xl(j,24) = 0.50d0*(xl(j,14) + xl(j,10) + xl(j,11) + xl(j,19))
     &             - 0.25d0*(xl(j, 2) + xl(j, 3) + xl(j, 6) + xl(j, 7))
        end do ! j
        ixl(24) = 24
      endif

      if(ixl(25).eq.0) then
        do j = 1,3
          xl(j,25) = 0.50d0*(xl(j,15) + xl(j,11) + xl(j,12) + xl(j,20))
     &             - 0.25d0*(xl(j, 3) + xl(j, 4) + xl(j, 7) + xl(j, 8))
        end do ! j
        ixl(25) = 25
      endif

      if(ixl(26).eq.0) then
        do j = 1,3
          xl(j,26) = 0.50d0*(xl(j,16) + xl(j,12) + xl(j, 9) + xl(j,21))
     &             - 0.25d0*(xl(j, 4) + xl(j, 1) + xl(j, 8) + xl(j, 5))
        end do ! j
        ixl(26) = 26
      endif

c     Center node

      if(ixl(27).eq.0) then
        do j = 1,3
          xl(j,27) = 0.25d0*(xl(j,13) + xl(j,14) + xl(j,15) + xl(j,16)
     &                     + xl(j,18) + xl(j,19) + xl(j,20) + xl(j,21)
     &                     + xl(j,23) + xl(j,24) + xl(j,25) + xl(j,26)
     &                     - xl(j, 1) - xl(j, 2) - xl(j, 3) - xl(j, 4)
     &                     - xl(j, 5) - xl(j, 6) - xl(j, 7) - xl(j, 8))

        end do ! j
        ixl(27) = 27
      endif

      end
