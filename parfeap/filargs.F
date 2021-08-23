c$Id: filargs.F,v 1.1 2006/11/21 16:44:39 rlt Exp $
      subroutine filargs(nargs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set up filename extender for input/output, etc.

c      Inputs:
c         None

c      Outputs:
c         File names returned in common /comfil/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

#     include  "include/finclude/petsc.h"
#     include  "pfeapa.h"

#     include  "comfil.h"

      logical   pcomp
      integer   nargs, ibeg, iend, ires

      save

c     Parallel

      iend = index(finp,' ') - 1
      ibeg = max(1,iend - 4)
      ires = iend - ibeg + 1
      if(ntasks.gt.1) then
        if(.not.pcomp(finp(ibeg:iend),'_0001',ires)) then
          ibeg = iend + 1
          iend = iend + 5
          finp(ibeg:iend) = '_0001'
          fout            = finp
          fout(1:1)       = 'O'
          fres            = finp
          fres(1:1)       = 'R'
          fsav            = finp
          fsav(1:1)       = 'R'
          fplt            = finp
          fplt(1:1)       = 'P'
        endif

c     Scalar

      else
        if(pcomp(finp(ibeg:iend),'_0001',ires)) then
          finp(ibeg:iend) = '     '
          fout            = finp
          fout(1:1)       = 'O'
          fres            = finp
          fres(1:1)       = 'R'
          fsav            = finp
          fsav(1:1)       = 'R'
          fplt            = finp
          fplt(1:1)       = 'P'
        endif

c     N.B. May need to remove call to 'doargs' or convert for computer used.
c          Call to doargs not allowed in parallel version.

        call doargs(finp,fout,fres,fsav,fplt,nargs)

      endif

      end
