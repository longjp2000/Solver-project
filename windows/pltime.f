c$Id: pltime.f,v 1.1 2006/11/20 20:34:56 rlt Exp $
      subroutine pltime()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place time on plot window

c      Inputs:  None

c      Outputs: To plot window
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'tdata.h'
      include  'pdatxt.h'

      character yy*15

      save

      data      yy / ' ' /

c     Display time for current view

      dtext = 0.11d0
      call pppcol(-1,1)
      call tplot(1.02d0 , 0.135d0, yy, 15, 1)
      call pppcol(1,1)
      write(yy, '(6hTime =,1p,1e9.2)' ) ttim
      call tplot(1.02d0 , 0.135d0, yy, 15, 1)

      end
