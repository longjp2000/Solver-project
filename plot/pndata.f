c$Id: pndata.f,v 1.1 2006/11/20 20:33:12 rlt Exp $
      subroutine pndata(tx,ct,nix,nxd,nxn,labl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Collected nodal data for plots

c      Inputs:
c         tx(2)     - Text identifier data
c         ct(3)     - Plot command parameters

c      Outputs:
c         none      - Plot outputs to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      character tx(2)*15
      logical   labl
      integer   nix,nxd,nxn
      real*8    ct(3)

      save

      write(*,*) ' *ERROR*: NDATa is allowed only in parallel version'

      end
