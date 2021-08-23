c$Id: plclos.f,v 1.1 2006/11/20 20:34:56 rlt Exp $
      subroutine plclos()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Close plot device

c      Inputs:
c         none

c      Outputs:
c         none      - Returns command outputs to text device
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'print.h'

      integer   status,vtxwin

      save

c     Close plot device

      fopn = .false.

      status = vtxwin()

      end
