c$Id: filargs.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine filargs(nargs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Serial version to get command line arguments

c      Inputs:
c         None

c      Outputs:
c         File names returned in common /comfil/
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'comfil.h'

      integer   nargs

      save

c     N.B. May need to remove call to 'doargs' or convert for computer used.
c          Call to doargs not allowed in parallel version.

      call doargs(finp,fout,fres,fsav,fplt,nargs)

      end
