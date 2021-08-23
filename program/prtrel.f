c$Id: prtrel.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine prtrel(r,ttim,ndf,nlist,list,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dabs' to 'abs'                           17/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal reactions based on list

c      Inputs:
c         r(*)      - Current value of reactions
c         ttim      - Value of solution time
c         ndf       - Number dof/node
c         nlist     - Number items in list
c         list(*)   - List of node numbers to output
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'iofile.h'

      logical   prth
      integer   ndf,nlist, i,k,n,nn, count, list(*)
      real*8    ttim, r(ndf,*),rsum(6),asum(6),psum(6)

      save

      do k = 1,ndf
        psum(k) = 0.d0
        rsum(k) = 0.d0
        asum(k) = 0.d0
      end do ! k
      do i = 1,numnp
        do k = 1,ndf
          rsum(k) = rsum(k) - r(k,i)
          asum(k) = asum(k) + abs(r(k,i))
        end do ! k
      end do ! i

      count = 0

      do nn = 1,nlist
        n = list(nn)
        count = count - 1
        do k = 1,ndf
          psum(k) = psum(k) - r(k,n)
        end do ! k
        if(count.le.0) then
          call prtitl(prth)
          write(iow,2000) ttim,(k,k=1,ndf)
          if(ior.lt.0) then
            write(*,2000) ttim,(k,k=1,ndf)
          endif
          count = 50000000
        endif
        if(ior.lt.0) then
          write(*,2001) n,(-r(k,n),k=1,ndf)
        endif
        write(iow,2001) n,(-r(k,n),k=1,ndf)
      end do ! nn

c     Print sum checks

      write(iow,2002) (psum(k),k=1,ndf)
      write(iow,2003) (rsum(k),k=1,ndf)
      write(iow,2004) (asum(k),k=1,ndf)
      if(ior.lt.0) then
        write(*,2002) (psum(k),k=1,ndf)
        write(*,2003) (rsum(k),k=1,ndf)
        write(*,2004) (asum(k),k=1,ndf)
      endif

c     Formats

2000  format('  N o d a l    R e a c t i o n s',12x,
     &  'Time',e18.5//'   Node',6(i8,' dof'):/(7x,6(i8,' dof')))

2001  format(i7,1p,6e12.4:/(7x,1p,6e12.4))

2002  format(/' Pr.Sum',1p,6e12.4:/(7x,1p,6e12.4))

2003  format( '   Sum ',1p,6e12.4:/(7x,1p,6e12.4))

2004  format( '  |Sum|',1p,6e12.4:/(7x,1p,6e12.4))

      end
