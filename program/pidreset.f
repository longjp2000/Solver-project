c$Id: pidreset.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pidreset(id)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Renumber equations when all equations are active.
c               Used in parallel version with blocked form

c      Inputs:
c        id(ndf,*)  - Equation numbers compressed

c      Outputs:
c        id(ndf,*)  - Equation numbers uncompressed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'cdata.h'
      include   'sdata.h'

      integer    id(ndf,*),  i,n

      do n = 1,numnp
        do i = 1,ndf
          id(i,n) = ndf*(n-1) + i
        end do ! i
      end do ! k

      end
