c$Id: prtrea.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine prtrea(r,x,ndm,ndf,n1,n2,n3,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Change 'dabs' to 'abs'                           17/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output nodal reactions for current solution

c      Inputs:
c         r(*)      - Current value of reactions
c         x(ndm,*)  - Nodal coordinates of mesh
c         ndm       - Spatial dimension of mesh
c         ndf       - Number dof/node
c         n1        - First node to output
c         n2        - Last node to output
c         n3        - Increment to n1
c         prth      - Output title/header data if true

c      Outputs:
c         None      - Outputs to file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'xtout.h'
      include  'pointer.h'
      include  'comblk.h'

      include  'p_int.h'

      logical   prth,nzprt
      integer   ndm,ndf,n1,n2,n3, i,k,n,  count, nxt1
      real*8    x(ndm,*),r(ndf,*),rsum(6),asum(6),psum(6)

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
      nxt1  = max(1,nxt)
      fp(1) = np(190) - 1
      do n = n1,n2,n3
        if( (mr(fp(1)+n).ge.0) .and. (nxt.eq.0 .or.
     &      (abs(x(nxt1,n)-xt).le.xtol) )  ) then
          nzprt = .false.
          do k = 1,ndf
            psum(k) = psum(k) - r(k,n)
            if(r(k,n).ne.0.0d0) then
              nzprt = .true.
            endif
          end do ! k
          if(nzprt) then
            count = count - 1
            if(count.le.0) then
              call prtitl(prth)
              write(iow,2000) (k,k=1,ndf)
              if(ior.lt.0.and.pfr) then
                write(*,2000) (k,k=1,ndf)
              endif
              count = 50000000
            endif
            if(ior.lt.0.and.pfr) then
              write(*,2001) n,(-r(k,n),k=1,ndf)
            endif
            write(iow,2001) n,(-r(k,n),k=1,ndf)
          endif
        endif
      end do ! n

c     Print sum checks

      write(iow,2002) (psum(k),k=1,ndf)
      write(iow,2003) (rsum(k),k=1,ndf)
      write(iow,2004) (asum(k),k=1,ndf)
      if(ior.lt.0.and.pfr) then
        write(*,2002) (psum(k),k=1,ndf)
        write(*,2003) (rsum(k),k=1,ndf)
        write(*,2004) (asum(k),k=1,ndf)
      endif

c     Formats

2000  format('  N o d a l    R e a c t i o n s'//
     &  '   Node',6(i8,' dof'):/(7x,6(i8,' dof')))

2001  format(i7,1p,6e12.4:/(7x,1p,6e12.4))

2002  format(/' Pr.Sum',1p,6e12.4:/(7x,1p,6e12.4))

2003  format( '   Sum ',1p,6e12.4:/(7x,1p,6e12.4))

2004  format( '  |Sum|',1p,6e12.4:/(7x,1p,6e12.4))

      end
