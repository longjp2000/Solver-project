c$Id: pprtd.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pprtd

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Show current dictionary entries

c      Inputs:
c         none

c      Outputs:
c         none      - To screen/file
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'allotn.h'
      include  'allotd.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'memuse.h'
      include  'pdata2.h'

      include  'pointer.h'

      character c*1, ptype(2)*8
      logical   skip
      integer   i, j, ity, lines

      save

      data      lines /9/
      data      ptype /' Program',' User   '/

c     Output dictionary names

      totimem = 0
      totrmem = 0
      skip   = idev.ne.1 .and. ior.lt.0
      if(ior.lt.0) then
        write(*,2000)
      endif
      write(iow,2000)
      do j = 1,ndict,lines
        do i = j,min(j+lines-1,ndict)
          if(ddict(i).eq.pdict(i)) then
            ity = 1
          else
            ity = 2
          endif
          if(ior.lt.0) then
            write(*,2001) i,dict(i),ddict(i),iprec(i),ipoint(i)
     &                     ,np(pdict(i)),ptype(ity)
          endif
          write(iow,2001) i,dict(i),ddict(i),iprec(i),ipoint(i)
     &                     ,np(pdict(i)),ptype(ity)
          if(iprec(i).eq.1) then
            totimem = totimem + ipoint(i)
          else
            totrmem = totrmem + ipoint(i)
          endif
        end do ! i
        if(skip .and. min(j+lines,ndict).ne.ndict) then
          write(*,*) '   ** PRESS ENTER **'
          read(*,1000) c
          write(  *,2000)
        endif
      end do ! j
      if(maxuse.gt.0) then
        if(skip) then
          write(*,*) '   ** PRESS ENTER **'
          read(*,1000) c
        endif
        if(ior.lt.0) then
          write(*,2002) totimem,totrmem,totimem+totrmem,maxuse
        endif
        write(iow,2002) totimem,totrmem,totimem+totrmem,maxuse
      else
        if(ior.lt.0) then
          write(*,2002) totimem,totrmem
        endif
        write(iow,2002) totimem,totrmem
      endif

c     Formats

1000  format(a)
2000  format(5x,'D i c t i o n a r y    o f   A r r a y s'//
     & 10x,' Entry  Array   Array  Array    Array            Pointer'/
     & 10x,'Number  Names  Number  Precn   Length          Value',
     &     ' Type')

2001  format(10x,i5,3x,a5,2i7,1i9,1i16,a8)
2002  format(10x,'Total memory used by FEAP:'/
     &       20x,'Integer Arrays = ',1i9/
     &       20x,'Real    Arrays = ',1i9:/
     &       20x,'Total used     = ',1i9/
     &       20x,'Total allowed  = ',1i9)

      end
