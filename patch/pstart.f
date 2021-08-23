c$Id: pstart.f,v 1.1 2006/11/20 20:44:57 rlt Exp $
      subroutine pstart()

c      * * F E A P * * A Finite Element Analysis Program
c                        -      -       -        -
c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Start graphical outputs: Windows version with 'open' but
c               no file list inputs from window

c      Inputs: none

c      Outputs: none
c-----[--+---------+---------+---------+---------+---------+---------+-]
      use        dflib
      use        dfwin

      implicit   none

      include   'pdata2.h'
      include   'prmptd.h'

      logical*4  ldrv

      interface
        logical(4)  function initialsettings
        end function
      end interface

c     Files read from filnam inputs

      fileck = .true.   ! File checking at startup is on

c     Graphics Driver number for PC version

      idev = 2

c     Initialize memory

      call pinitm()

c     Start Windows

      call pwopn()

c     Check user installation options

      call pinstall()

c     Set initial file names

      call filnam()

      end

      logical(4)   function initialsettings ( )

      use          msflib

      implicit     none

      type(qwinfo) winfo
      logical      status, pcomp, dflag
      character    record*80,str*80

      save

c     Maximize Frame

      winfo.type = qwin$max
      status     = setwsizeqq(qwin$framewindow,winfo)

      initialsettings = .true.

      end function
