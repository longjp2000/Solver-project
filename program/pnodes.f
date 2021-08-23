c$Id: pnodes.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pnodes(xs,ndm,prt,prth)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Input super-node coordineates

c      Inputs:
c         ndm     - Spatial dimension of mesh
c         prt     - Print input values if true
c         prt     - Print header if true

c      Outputs:
c         xs(3,*) - Super-node coordinates
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'iofile.h'

      logical   prt,prth, setvar,pinput
      integer   ndm, i,n,node
      real*8    xs(3,*), td(6)

      save

      if(prt) then
        call prtitl(prth)
        write(iow,2000) (n,n=1,ndm)
        if(ior.lt.0) then
          write(*,2000) (n,n=1,ndm)
        endif
      endif

c     Input coordinates of supernodes for blending inputs

1     if(ior.lt.0) then
        write(*,3000)
        call pprint('   >')
      endif

      setvar = pinput(td,ndm+1)

      node   = td(1)

      if(node.gt.0) then

        numsn = max(numsn,node)
        do n = 1,ndm
          xs(n,node) = td(n+1)
        end do ! n

        if(prt) then
          write(iow,2001) node,(xs(i,node),i=1,ndm)
          if(ior.lt.0) then
            write(*,2001) node,(xs(i,node),i=1,ndm)
          endif
        endif

      else
        return
      endif

      go to 1

3000  format('   Input: node, xs(i,node),i=1,ndm')

2000  format('   S u p e r N o d e   C o o r d i n a t e s'//
     &   '       Node  ',i1,'-Coordinate':,'   ',i1,'-Coordinate':,
     &   '   ',i1,'-Coordinate'/)

2001  format(i10,1p,3e15.5)

      end
