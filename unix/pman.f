c$Id: pman.f,v 1.1 2006/11/20 20:33:21 rlt Exp $
      subroutine pman(name,nn)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Output manual page to screen.  Manual page location
c               set in main program.  Must have extender '.t'

c      Inputs:
c         name      - Name of command to display
c         nn        - Type of command: 1=mesh, 2=macro, 3=plot

c      Outputs:
c         To screen
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iodata.h'
      include  'pathn.h'

      character name*4,filel*64
      integer   nn, iposl, ipos

      save

c     Path name to where manual files is stored in 'file'
c     set ipos to last entry in path

      filel             = ' '
      filel(1:9)        = 'acroread '
      iposl             = ipos(file(nn),44)
      filel(10:iposl+9) = file(nn)

c     Set for Mesh

      if(    nn.eq.1) then
        filel(iposl+10:iposl+18) = 'appxa.pdf'

c     Set for Mesh

      elseif(nn.eq.2) then
        filel(iposl+10:iposl+18) = 'appxd.pdf'

c     Set for Command Language

      elseif(nn.eq.3) then
        filel(iposl+10:iposl+18) = 'appxe.pdf'

c     Set for Manipulation

      elseif(nn.eq.4) then
        filel(iposl+10:iposl+18) = 'appxb.pdf'

c     Set for Contact

      elseif(nn.eq.5) then
        filel(iposl+10:iposl+18) = 'appxc.pdf'

c     Error

      else
        write(*,*) ' NO HELP '
        return
      endif

      call system(filel)

      end
