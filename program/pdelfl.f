c$Id: pdelfl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pdelfl()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Clean up old files by erasing temporary mesh input files

c      Inputs:
c         none

c      Outputs:
c         none
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'comfil.h'
      include   'compas.h'
      include   'cornum.h'
      include   'iodata.h'
      include   'part1.h'
      include   'part3.h'
      include   'part8.h'
      include   'setups.h'

      include   'pointer.h'
      include   'comblk.h'

      include   'p_int.h'

      character  fname*132,fext*4,ftyp(13)*4, flnk(4)*4
      logical    exst,isopen, flags(5)
      integer    i,m,mmx

      save

      data       fp   / 10*0 /

      data       ftyp /'an0','bn0','co0','ds0','ep0','fr0',
     &                 'gb0','in0','ld0','sl0','ud0','ww0',
     &                 'yp0'/
      data       flnk /'eln','lnk','jnt','tem'/

c     Check if 'ios' open

      inquire(unit = ios, opened = exst)
      if(exst) then
        close(unit = ios)
      endif

c     Check for maximum number of files active

      mmx = max(nsurf,nbouf,ndisf,nforf,nangf) - 1

c     Delete mesh generation files

      do m = 0,mmx

        do i = 1,12
          fname = fsav
          fext  = ftyp(i)
          if(m.le.9) then
            write(fext(3:3),'(i1)') m
          else
            write(fext(2:3),'(i2)') m
          endif
          call addext(fname,fext,128,4)
          inquire(file=fname,exist=exst,opened=isopen)
          if(exst) then
            if(.not.isopen) then
              open (unit = ios, file = fname, status = 'old')
            endif
            close(unit = ios, status = 'delete')
          endif
        end do ! i
      end do ! m

c     Delete mesh manipulation files

      do m = 1,4
        fname = fsav
        fext  = flnk(m)
        call addext(fname,fext,128,4)
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif
      end do ! m

c     Delete 'feaploop' mesh blocks

      do m = 0,9
        fname = '      '
        write(fname,'(a9,i1)') 'feaploop.',m
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif
      end do ! m

c     Delete 'Material' file

      inquire(file= fmtl ,exist=exst,opened = isopen)
      if(exst) then
        if(.not.isopen) then
          open (unit = ios, file =  fmtl , status = 'old')
        endif
        close(unit = ios, status = 'delete')
      endif

c     Delete out-of-core blocks (additional to look for contacts)

      do m = 1,min(999,maxbl+10)

c       Check upper

        fname = '      '
        if(m.lt.10) then
          write(fname,'(a6,i1)') 'Aupper',m
        elseif(m.lt.100) then
          write(fname,'(a6,i2)') 'Aupper',m
        else
          write(fname,'(a6,i3)') 'Aupper',m
        endif
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif

c       Check lower

        fname = '      '
        if(m.lt.10) then
          write(fname,'(a6,i1)') 'Alower',m
        elseif(m.lt.100) then
          write(fname,'(a6,i2)') 'Alower',m
        else
          write(fname,'(a6,i3)') 'Alower',m
        endif
        inquire(file=fname,exist=exst,opened=isopen)
        if(exst) then
          if(.not.isopen) then
            open (unit = ios, file = fname, status = 'old')
          endif
          close(unit = ios, status = 'delete')
        endif

      end do ! m

c     Purge solver storage if necessary

      do i = 1,5
        if(.not.tflp(i)) then
          if(nsolver(i)) then
            fp(1) = nap (i)
            fp(2) = naup(i)
            fp(3) = nalp(i)
            fp(4) = np(20+i)
            call psolve(nittyp(i),hr(np(26)),fp,
     &                 .false.,.false.,.false.,.false.)
          else
            flags(1) = .false.
            flags(2) = .false.
            flags(3) = .false.
            flags(4) = .false.
            flags(5) = .true.
            call usolve(flags,hr(np(26)))
          endif
        endif
      end do ! i

      end
