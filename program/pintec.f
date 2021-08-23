c$Id: pintec.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pintec()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Body force and reaction inputs from coordinates

c      Inputs:
c         none      - Data retrieved through common blocks

c      Outputs:
c         none      - Data stored in blank common locations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'cdat2.h'
      include  'comfil.h'
      include  'conval.h'
      include  'corset.h'
      include  'corfil.h'
      include  'cornum.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'print.h'
      include  'sdata.h'
      include  'trdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prtt, oprt,oprth
      character fnamr*132, fext*4, type*4
      integer   i,l1, iorsv,n

      save

c     Set body force internal loads

      if(intfl) then

        do l1 = 0,nintf-1
          fext = 'in0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          fnamr =  fsav
          call addext(fnamr,fext,128,4)
          call opnfil(fext,fnamr,-1,ios,prtt)

c         Read data from file

          iorsv = ior
          ior   = ios

          oprt  = prt
          oprth = prth
          do i = 0,36
            do n = 1,26
              vvsave(n,i) = vvv(n,i)
            end do ! n
          end do ! i
          do i = 1,3
            x0sav(i) = x0(i)
          end do ! i

          read(ior,1000) type,fincld(isf),irecrd(isf),prt,prth
          read(ior,1001) vvv
          read(ior,1001) tr,xr,trdet,x0

          prtt = prth .and. l1.eq.0

          call pbodyf(mr(np(33)),ndf,nen1,numel,prt,prtt)

          close(ior)
          ior  = iorsv

          prt  = oprt
          prth = oprth
          do i = 0,36
            do n = 1,26
              vvv(n,i) = vvsave(n,i)
            end do ! n
          end do ! i
          do i = 1,3
            x0(i) = x0sav(i)
          end do ! i

        end do ! l1
      endif

c     Set reaction data

      if(reafl) then
        call preain(mr(np(31)+ndf*numnp),hr(np(27)),hr(np(41)),
     &              ndf,numnp,reafi,prt,prth)
      endif

c     Formats

1000  format(a4,2x,a12,i8,2l5)
1001  format(1p,4e20.12)

      end
