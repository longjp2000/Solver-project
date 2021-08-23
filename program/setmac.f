c$Id: setmac.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine setmac(yyy,wd,nwd,ll,jct,lct,lzz)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Set macro commands for procedures use

c      Inputs:
c         yyy      - String with command name
c         wd(nwd)  - Active command names
c         nwd      - Number of active commands

c      Outputs:
c         ll       - Execution command numbers
c         jct(*)   - List of command numbers
c         lct(*)   - String for command option
c         lzz(*)   - String for command paramaters
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      logical   pcomp
      integer   ll, i,j, nwd, jct(*)
      character clab1*4,clab2*15,wd(nwd)*4,lct(*)*15,xxx*80,yyy*80
      character lzz(*)*80, temp*47

      save

c     Repack parameter string

      read(yyy,1000) clab1,clab2,xxx

      temp( 1:15) = xxx( 1:15)
      temp(16:16) = ','
      temp(17:31) = xxx(16:30)
      temp(32:32) = ','
      temp(33:47) = xxx(31:45)

c     Strip any blanck characters

      xxx = ' '
      j = 0
      do i = 1,47
        if(temp(i:i).ne.' ') then
          j = j + 1
          xxx(j:j) = temp(i:i)
        end if
      end do ! i

c     Insert data into command list

      do i = 1,nwd
        if(pcomp(clab1,wd(i),4)) then
          ll      = ll + 1
          jct(ll) = i
          lct(ll) = clab2
          lzz(ll) = xxx
          return
        end if
      end do ! i
      write(*,3000) clab1

c     Formats

1000  format(a4,11x,a15,a)

3000  format(' *WARNING* Illegal command ',a4,' in procedure.')

      end
