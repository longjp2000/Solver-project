c$Id: sumcnt.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine sumcnt(ic,nneq,kp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Sum counts to get storage of connections

c      Inputs:
c         ic(*)  - Pointer array for equation connection list 
c         nneq   - Number total degrees of freedom

c      Outputs:
c         kp     - Length of sparse array
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    i,kp,nneq,ic(*)

      save

c     Set up pointers.

      do i = 2, nneq
         ic(i) = ic(i) + ic(i-1)
      end do ! i

      kp = ic(nneq)

      end
