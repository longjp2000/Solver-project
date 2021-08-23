c$Id: pintio.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pintio(y,ifld)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Character string input routine for data
c               N.B. This routine has largely been superceded by
c                    pinput and tinput functions.

c      Inputs:
c         lfld   - Field width to separate input data items into.

c      Outputs:
c         y      - Character string input, in field widths of ifld
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'chdata.h'
      include  'comfil.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'iosave.h'

      character y*(*)
      integer   ifld

      save

c     Read a record from input file

100   if (ior.gt.0)  then
        read (ior,1000,err=901,end=902) record
        irecrd(isf) = irecrd(isf) + 1
        xxx         = record
      else
        read (  *,1000,err=901) xxx
      endif
      if(lsave) write(lfile,1000) xxx

c     Adjust and move record to y-array

      call acheck(xxx,y,ifld,255,255)
      return

c     Read error encountered

901   call  errclr ('PINTIO')
      goto  100

c     Eof encountered

902   call  endclr ('PINTIO',xxx)
      goto  100

c     Format

1000  format(a)

      end
