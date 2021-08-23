c$Id: plabl.f,v 1.1 2006/11/20 20:34:56 rlt Exp $
      subroutine plabl(m)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Place numerical labels on plot data
c               N.B. Must be preceded by a move to location where
c                    value is centered

c      Inputs:
c         m         - Number to place on plot

c      Outputs:
c         none      - Plot output to screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'pdatap.h'
      include  'pdatps.h'
      include  'pdatxt.h'
      include  'plflag.h'
      include  'psdat3.h'

      character yyy*6
      integer   m,n,nchar
      real*8    x1,y1

      save

c     Set number of characters

      n = abs(m)

      if    (n.ge.0   .and. n.lt.10   ) then
        write(yyy,'(i1)') n
        nchar = 1
      elseif(n.ge.10  .and. n.lt.100  ) then
        write(yyy,'(i2)') n
        nchar = 2
      elseif(n.ge.100 .and. n.lt.1000 ) then
        write(yyy,'(i3)') n
        nchar = 3
      elseif(n.ge.1000 .and. n.lt.10000 ) then
        write(yyy,'(i4)') n
        nchar = 4
      elseif(n.ge.10000 .and. n.lt.100000 ) then
        write(yyy,'(i5)') n
        nchar = 5
      else
        write(yyy,'(i6)') n
        nchar = 6
      endif

      dtext = 0.0d0
      x1    = jx1/22000.d0
      y1    = jy1/22000.d0
      if(screfl) call tplot(x1,y1,yyy,nchar,0)

      end
