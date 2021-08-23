c$Id: parexp.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine parexp(x,xs,v,nex,error)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Identify parenthetical expressions and evaluate

c      Inputs:
c         x(*)     - String containing expression to evaluate

c      Scratch:
c         xs(*)    - Array used to temporarily store expression
c         v(*)     - Array to hold values

c      Outputs:
c         x(*)     - Expression replaced by upper case letter
c         nex      - Number of upper case letters used
c         error    - Flag, true if error occurs

c      Common returns:
c         www(*)   - Upper case letters with values assigned
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'conval.h'

      logical   error
      character x(*)*1,xs(*)*1
      integer   i,j,k,l, i1,i2,nex
      real*8    val, v(*)

      save

c     Find parenthetical expressions and remove

      do i = 1,75
        if(x(i).eq.'(') then
          i1 = i + 1
          do j = i1,75
            if(x(j).eq.'(') then
              call errclr('PAREXP')
              call plstop()
            elseif(x(j).eq.')') then
              do l = 1,j-i+1
                xs(l) = ' '
              end do ! l
              i2 = j - 1
              if(i2.lt.i1) then
                call errclr('PAREXP')
                call plstop()
              else
                k = 0
                do l = i1,i2
                  k = k + 1
                  xs(k) = x(l)
                  x(l)  = ' '
                end do ! l
                x(i2+1)  = ' '

c               Evaluate expression in parenthesis

                call evalex(xs,v,val,k,error)
                if(error) return
                nex = nex + 1
                www(nex) = val

c               Put upper case letter in expression and close up remainder

                x(i) = char(nex +64)
                i2 = i2 -i1 + 2
                do l = i1,75
                  x(l) = ' '
                  if(l+i2.le.75) then
                    x(l) = x(l+i2)
                  endif
                end do ! l
              endif
              go to 100
            endif
          end do ! j
100       continue
        endif
      end do ! i

      end
