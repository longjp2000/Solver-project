c$Id: readfl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      logical function readfl(tx)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Read data from specified file.
c               (teminates with: read,end).

c      Inputs:
c         tx        - Name of file for reads

c      Outputs:
c         readfl    - Status of mesh input
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'iosave.h'

      logical   pcomp, lopen, lexist
      character tx*(*), fnamr*21, fext*8
      integer   isfile

      save

c     Set the default for returning error or close

      readfl = .false.

c     Close file and reset logical unit number

      if(pcomp(tx,'end',3)) then
        inquire(file=fnamr,opened=lopen,exist=lexist)
        if(lexist.and.lopen) close(lfile)
        if(lread) then
          ior    = isfile
          lread  = .false.
          readfl = .false.
        endif

c     Open file and set new logical unit number

      else
        fnamr    = tx
        fext     = tx
        call opnfil(fext,fnamr,-2,lfile,lread)
        if(.not.lread) return
        isfile = ior
        ior    = lfile
        readfl = .true.
      endif

      end
