c$Id: umesh.f,v 1.1 2006/11/20 20:33:27 rlt Exp $
      logical function umesh(cc,prt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: User mesh command interface

c      Inputs:
c         cc     - User command option
c         prt    - Output if true

c      Outputs:
c         umesh  - Flag to indicate if successful in matching command
c         N.B. Users must provide other output via common blocks, etc.

c     IOR is logical unit number for data inputs (if negative input from *)
c     IOW is logical unit number for output data
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      logical   prt,pcomp
      character cc*4

c     Match on 'USER':  Can add as many checks as desired with 'user'
c                       replaced by 4-character word for each command.

      if(pcomp(cc,'user',4)) then

        if(prt) write(  *,*) ' '
        if(prt) write(  *,*) '    USER COMMAND = ',cc
        if(prt) write(  *,*) ' '

c       Return TRUE to indicate command was successful

        umesh = .true.

      elseif(pcomp(cc,'use2',4)) then

c       Second user command location, etc. for others.

        umesh = .true.

      else

c       Return FALSE for no match on user mesh command

        umesh = .false.

      endif

      end
