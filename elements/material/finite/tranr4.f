c$Id: tranr4.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine tranr4(fl,fr,t)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:   Transformation array for 4th rank tensor
c                 t(a,b) = fl(i,I)*fr(j,J) : a -> I,J ; b -> i,j
c                   a,b  |  1    2    3    4    5    6
c                  ------+-----------------------------
c                  (I,J) | 1,1  2,2  3,3  1,2  2,3  3,1
c               or (i,j) |                2,1  3,2  1,3

c     Input:
c       fl(3,3) - left  deformation gradient
c       fr(3,3) - right deformation gradient
c     Output:
c       t(6,6) - transformation array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,j, i1(6),i2(6)
      real*8    fl(3,3),fr(3,3), t(6,6)

      data      i1 /1,2,3,1,2,3/
      data      i2 /1,2,3,2,3,1/

c     Form transformation array for a 4th rank tensor in matrix form

      do i = 1,3
        do j = 1,3
          t(i,j) =  fl(i1(j),i1(i))*fr(i2(j),i2(i))
        end do ! j
        do j = 4,6
          t(i,j) = (fl(i1(j),i1(i))*fr(i2(j),i2(i))
     &           +  fl(i2(j),i2(i))*fr(i1(j),i1(i)))*0.5d0
        end do ! j
      end do ! i

      do i = 4,6
        do j = 1,3
          t(i,j) =  fl(i1(j),i1(i))*fr(i2(j),i2(i))
     &           +  fl(i2(j),i2(i))*fr(i1(j),i1(i))
        end do ! j
        do j = 4,6
          t(i,j) = (fl(i1(j),i1(i))*fr(i2(j),i2(i))
     &           +  fl(i2(j),i1(i))*fr(i1(j),i2(i))
     &           +  fl(i1(j),i2(i))*fr(i2(j),i1(i))
     &           +  fl(i2(j),i2(i))*fr(i1(j),i1(i)))*0.5d0
        end do ! j
      end do ! i

      end
