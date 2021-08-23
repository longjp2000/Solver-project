c$Id: plabl.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
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

      include  'pdata2.h'
      include  'pdatap.h'
      include  'pdatps.h'
      include  'pdatxt.h'
      include  'plflag.h'
      include  'psdat3.h'
      include  'x11f.h'

      character yyy*6
      integer   j,jj,m,n,nchar
      real*4    xp0(2)

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

c     X11

      x11(1) = (xp(2)+0.0125*(4-nchar+nchar/3))*xx(2)
      y11(1) = yp(2)*xx(3)*1.28
      y11(2) = 0
      do j = 1,nchar
         jj       = ichar(yyy(j:j))
         x11(j+1) = jj
      end do ! j
      if(screfl) then
        call gdx11(-1024-nchar,x11,y11)
      endif

c     PostScript

      if (hdcpy) then
        xp0(2) = xp(2) + 0.0125*(4-nchar+nchar/3)
        call fptplt(xp0,yp,yyy,nchar,0)
      endif

      end
