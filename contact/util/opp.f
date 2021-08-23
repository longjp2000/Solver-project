c$Id: opp.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      integer function opp (ncom,nfea,nopt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Giorgio Zavarise           April 10, 1996            1.0

c      Acronim: OPtions command Pointer

c      Purpose: Set the appropriate pointer for OPTIONS in CIS array

c      Inputs :
c         ncom    - # of the command
c         nfea    - # of the feature
c         nopt    - # of the option

c      Outputs:
c         opp     - absolute pointer
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'c_0.h'
      include  'c_comnd.h'
      include  'c_dicti.h'

      integer   ncom,nfea,nopt, ofs,k

      save

      ofs = ofsop
      do k = 1,ncom-1
        ofs = ofs + nop(k)*nfe(k)
      end do
      opp = ofs + nop(ncom)*(nfea-1) + nopt

      end
