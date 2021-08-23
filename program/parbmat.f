c$Id: parbmat.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine parbmat(phib,mb,cmass,lmass,t, neq)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute base matrix coupling term

c     Inputs:

c     Outputs:

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer    neq
      real*8     t(*),cmass(*),lmass(*),phib(*),mb(*)

c     Dummy routine in serial version

      end
