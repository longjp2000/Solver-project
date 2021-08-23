c$Id: pleigt.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine pleigt(eval)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c     1. Use 'pi' from 'pconstant.h'                        14/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place eigenvalue on plot window

c      Inputs:  None

c      Outputs: To plot window
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'evdata.h'
      include  'pconstant.h'
      include  'pdatxt.h'

      character yy*20
      real*8    eval,dd

      save

c     Set value to plot

      dtext = 0.00d0
      call pppcol(-1,1)
      yy = '                    '
      call tplot(1.13d0 , 0.160d0, yy, 20, 1)

c     Clear value for current vector

      if(imtyp.eq.1) then
        dd = sqrt(abs(eval))*0.5d0/pi
        write(yy, '(a7,1p,1e9.2,a4)' ) 'Value =',dd,' Hz.'
      else
        write(yy, '(a7,1p,1e9.2,a4)' ) 'Value =',eval,'    '
      endif

c     Display Value for current vector

      call pppcol(1,1)
      call tplot(1.13d0 , 0.160d0, yy, 20, 1)

      end
