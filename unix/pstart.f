c$Id: pstart.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine pstart()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:

c     Inputs:

c     Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit      none

      include  'memuse.h'
      include  'pfeapb.h'
      include  'prmptd.h'
      include  'comblk.h'

      save

c     Set flags for serial execution

      pfeap_on   = .false.
      pfeap_gnod = .false.

c     Set maximum memory use

      maxuse = 0

c     Start for X11 graphics driver

      call pdriver()

c     Set for file checking at startup

      fileck = .true.

c     Initialize memory

      call pinitm(hr(1),mr(1))

c     Check user installation options

      call pinstall()

c     Set initial file names

      call filnam()

      end
