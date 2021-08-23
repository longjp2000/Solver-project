c$Id: seproj.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine seproj(p,s,se,dt,st,ser,ix,nel,nen,numnp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Assemble element projection quantities into global ones.

c      Inputs:
c        p(*)    - Element nodal weights
c        s(*)    - Element projected quantities
c        se(*)   - Element error quantities
c        ix*)    - Element node numbers
c        nel     - Number element nodes
c        nen     - Number element nodes max
c        numnp   - Number global nodes

c      Outputs:
c        dt(*)   - Global nodal weights
c        st(*)   - Global projected quanties
c        ser(*)  - Global error quanties
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'strnum.h'

      integer    nel,nen,numnp, n,nn,i, ix(*)
      real*8     p(nen),s(nen,*),se(*),dt(*),st(numnp,*),ser(*)

      save

      do n = 1,nel
        nn = ix(n)
        if(nn.gt.0) then
          dt(nn)  = dt(nn)  + p(n)
          ser(nn) = ser(nn) + se(n)

          do i = 1,iste
            st(nn,i) = st(nn,i) + s(n,i)
          end do ! i
        endif
      end do ! n

      end
