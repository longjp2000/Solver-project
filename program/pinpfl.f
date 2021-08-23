c$Id: pinpfl.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pinpfl(filnm,fext, type, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose:  Open file for input of saved input data

c     Inputs:
c        filnam   - Character name of calling routine
c        fext     - Extender name of file

c     Outputs:
c        type     - Type of action: 'set' or 'add'
c        isw      - Switch: 1 = open; 2 = close
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'cdat2.h'
      include   'comfil.h'
      include   'conval.h'
      include   'iodata.h'
      include   'iofile.h'
      include   'ioincl.h'
      include   'print.h'
      include   'trdata.h'

      logical    lsave
      character  filnm*(*), fnam*132,fext*4, type*4
      integer    isw, iosave, i,n, ipos

      save

c     Open unit for data saves

      if(isw.eq.1) then

        fnam = fsav
        call addext(fnam,fext,128,4)
        call opnfil(fext,fnam,-2,ios,lsave)
        if(lsave) then
          iosave = ior
          ior    = ios

          do i = 0,36
            do n = 1,26
              vvsave(n,i) = vvv(n,i)
            end do ! n
          end do ! i
          do i = 1,3
            x0sav(i) = x0(i)
          end do ! i

          read(ios,1000) type,fincld(isf),irecrd(isf),prt,prth
          read(ios,1001) vvv
          read(ios,1001) tr,xr,trdet,x0
        else
          i = ipos(filnm,128)
          n = ipos(fnam ,128)
          write(iow,3000) filnm(1:i),fnam(1:n)
          write(ilg,3000) filnm(1:i),fnam(1:n)
          call plstop()
        endif

c     Close file and restore parameters

      else

        close(ior,status = 'delete')
        ior = iosave

        do i = 0,36
          do n = 1,26
            vvv(n,i) = vvsave(n,i)
          end do ! n
        end do ! i
        do i = 1,3
          x0(i) = x0sav(i)
        end do ! i

      endif

c     Formats

1000  format(a4,2x,a12,i8,2l5)
1001  format(1p,4e20.12)

3000  format(' *ERROR* ',a,': Edge data file',a,' does not exist')

      end
