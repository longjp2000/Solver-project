c$Id: pedgin.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pedgin()

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Control routine for data inputs based on edge coordinate

c      Inputs:
c         none      - Data retrieved through common blocks

c      Outputs:
c         none      - Data stored in blank common locations
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'conval.h'
      include  'edgdat.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'pointer.h'
      include  'print.h'
      include  'sdata.h'
      include  'trdata.h'
      include  'comblk.h'

      character fext*4, type*4
      logical   oprt,oprth
      integer   l1

      save

c     Set edge angle values

      oprt  = prt
      oprth = prth
      if(eanfl) then
        do l1 = 0,neang-1
          fext = 'ww0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          call pinpfl('PEDGIN',fext, type, 1)
          call peforc(hr(np(43)),hr(np(45)),mr(np(190)),
     &                ndm,1,numnp,type,prt,prth,'Angle')
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge boundary conditions

      if(ebcfl) then
        do l1 = 0,nebcs-1
          fext = 'co0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          call pinpfl('PEDGIN',fext, type, 1)
          call pedges(hr(np(43)),mr(np(31)+ndf*numnp),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'B.C.')
          call pinpfl('PEDGIN',fext, type, 2)

        end do ! l1
      endif

c     Set edge displacement values

      if(edifl) then
        do l1 = 0,nedis-1
          fext = 'ud0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          call pinpfl('PEDGIN',fext, type, 1)
          call peforc(hr(np(43)),hr(np(27)+ndf*numnp),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Displ')
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge force values

      if(efcfl) then
        do l1 = 0,nefrc-1
          fext = 'ld0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          call pinpfl('PEDGIN',fext, type, 1)
          call peforc(hr(np(43)),hr(np(27)),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Force')
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge proportional load numbers

      if(eprfl) then
        do l1 = 0,nepro-1
          fext = 'ep0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          call pinpfl('PEDGIN',fext, type, 1)
          type = 'set'
          call pedges(hr(np(43)),mr(np(29)),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Prop')
          call pinpfl('PEDGIN',fext, type, 2)
        end do ! l1
      endif

c     Set edge base conditions

      if(ebsfl) then
        do l1 = 0,nebas-1
          fext = 'gb0'
          if(l1.le.9) then
            write(fext(3:3),'(i1)') l1
          elseif(l1.le.99) then
            write(fext(2:3),'(i2)') l1
          endif
          call pinpfl('PEDGIN',fext, type, 1)
          call pedges(hr(np(43)),mr(np(125)),mr(np(190)),
     &                ndm,ndf,numnp,type,prt,prth,'Base')
          call pinpfl('PEDGIN',fext, type, 2)

        end do ! l1
      endif

      prt  = oprt
      prth = oprth

      end
