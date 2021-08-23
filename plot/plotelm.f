c$Id: plotelm.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine plotelm(ix, x,xl, ndm, nen, num, axs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Plot curved edges on individual elements

c     Inputs:
c       ix(nen)   - List of nodes on element
c       x(ndm,*)  - Nodal coordinates of element
c       ndm       - Spatial dimension of mesh
c       nen       - Number of nodes on element
c       num       - Plot if true
c       axs       - Axes at node axs


c     Output:
c       Graphical plot of element edges
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'pdata1.h'
      include   'ppers.h'
     
      logical    num
      integer    ndm,nen, axs, nel, ne, n,i,ii
      integer    ix(nen), lin(3), tri(3,3),qud(3,4), tet(3,6),brk(3,12)
      real*8     x(ndm,*), xl(ndm,*), xx(3,3), xx0(3)

      real*8     oscale,oscaleg,os0(2),odx(2),osx(2),ofact, sz

      save

      data       lin / 1,3,2 /
      data       tri / 1,4,2, 2,5,3, 3,6,1 /
      data       qud / 1,5,2, 2,6,3, 3,7,4, 4,8,1 /
      data       tet / 1,5,2, 2,6,3, 3,7,1, 1,8,4, 2,9,4, 3,10,4/
      data       brk / 1,13,2, 2,14,3, 3,15,4, 4,16,1,
     &                 5,17,6, 6,18,7, 7,19,8, 8,20,5,
     &                 1, 9,5, 2,10,6, 3,11,7, 4,12,8 /

c     Determine number of nodes on element

      do n = 1,nen
        if(ix(n).ne.0) then
          nel = n
        endif
      end do ! n

c     Set xl array to do scaling

      do n = 1,nel
        ii = ix(n)
        if(ii.gt.0) then
          do i = 1,ndm
            xl(i,n) = x(i,ii)
          end do ! i
        else
          do i = 1,ndm
            xl(i,n) = xl(i,1)
          end do ! i
        endif
      end do ! n

      do i = 1,ndm
        xx0(i) = xl(i,1)
      end do ! i
      sz = 0.0d0
      do n = nel,1,-1
        do i = 1,ndm
          xl(i,n) = xl(i,n) - xx0(i)
          sz      = max(sz,abs(xl(i,n)))
        end do ! i
      end do ! n

c     Save old scaling

      oscale  = scale
      oscaleg = scaleg
      ofact   = fact
      do i = 1,2
        os0(i) = s0(i)
        odx(i) = dx(i)
        osx(i) = sx(i)
      end do ! i

c     Compute new scaling

      if(kpers.ne.0) then
        call frame(xl,ndm,nel,-1)
      else
        call frame(xl,ndm,nel, 1)
      endif

c     Set extra xx to zero

      do i = ndm+1,3
        do n = 1,3
          xx(i,n) = 0.0d0
        end do ! n
      end do ! i

c     3-node line

      if(nel.eq.3) then

        do n = 1,3
          ii = ix(lin(n))
          do i = 1,ndm
            xx(i,n) = x(i,ii) - xx0(i)
          end do ! i
        end do ! n
        call pltqln(xx, 20)

c     6-node triangle

      elseif(nel.eq.6) then

        do ne = 1,3
          do n = 1,3
            ii = ix(tri(n,ne))
            do i = 1,ndm
              xx(i,n) = x(i,ii) - xx0(i)
            end do ! i
          end do ! n
          call pltqln(xx, 20)
        end do ! ne

c     8/9-node quadrilateral

      elseif(nel.eq.8 .or. nel.eq.9) then

        do ne = 1,4
          do n = 1,3
            ii = ix(qud(n,ne))
            do i = 1,ndm
              xx(i,n) = x(i,ii) - xx0(i)
            end do ! i
          end do ! n
          call pltqln(xx, 20)
        end do ! ne

c     10-node tetrahedron

      elseif(nel.eq.10) then

        do ne = 1,6
          do n = 1,3
            ii = ix(tet(n,ne))
            do i = 1,ndm
              xx(i,n) = x(i,ii) - xx0(i)
            end do ! i
          end do ! n
          call pltqln(xx, 20)
        end do ! ne

c     20/27-node brick

      elseif(nel.eq.20 .or. nel.eq.27) then

        do ne = 1,12
          do n = 1,3
            ii = ix(brk(n,ne))
            do i = 1,ndm
              xx(i,n) = x(i,ii) - xx0(i)
            end do ! i
          end do ! n
          call pltqln(xx, 20)
        end do ! ne

      endif

c     Add nodes and numbers

      call pppcol (5,0)
      call pltelnd(xl,ix,ndm,nel, num)

      if(axs.gt.0) then
        call pltaxs(xl(1,axs),ndm,0.1d0*sz)
      endif

c     Restore old scaling

      scale  = oscale
      scaleg = oscaleg
      fact   = ofact
      do i = 1,2
        s0(i) = os0(i)
        dx(i) = odx(i)
        sx(i) = osx(i)
      end do ! i

      end
