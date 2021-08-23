c$Id: pmacr6.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine pmacr6 (lct,ct,prt,j)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Command language instruction subprogram: Part 6

c      Inputs:
c         lct      - Command option for current command
c         ct(3)    - Command parameters for current command
c         prt      - Flag, print data if true
c         j        - Number of command to execute

c      Outputs:
c         Depends on value of command j
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'codat.h'
      include  'endata.h'
      include  'fdata.h'
      include  'iodata.h'
      include  'iofile.h'
      include  'ioincl.h'
      include  'iosave.h'
      include  'part0.h'
      include  'part3.h'
      include  'part7.h'
      include  'pbody.h'
      include  'plcapt.h'
      include  'plist.h'
      include  'prflag.h'
      include  'prlod.h'
      include  'ptdat1.h'
      include  'ptdat2.h'
      include  'ptdat3.h'
      include  'ptdat4.h'
      include  'ptdat5.h'
      include  'ptdat6.h'
      include  'ptdat7.h'
      include  'ptdat8.h'
      include  'ptdat9.h'
      include  'ptdata.h'
      include  'ptdatb.h'
      include  'sdata.h'
      include  'tdata.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   prt, pcomp, errck
      character lct*15, yyy*15
      integer   i, j, n, msopr
      real*8    propld, ct(3)

      save

c     Macro instruction subprogram - part 6.

c     [list,,i] i = list number (maximum is 3)

      if(j.eq.1) then

        i = ct(1)
        i = max(1,min(i,3))
        write(iow,2002) i
        if(ior.lt.0) write(*,2002) i
        call prlist(iolist(1,i), niols(i))

c     [tplo]t,,<interval> - time history plots

c      Options : disp n1 n2 x y z
c                velo n1 n2 x y z
c                acce n1 n2 x y z
c                reac n1 n2 x y z
c                stre n1 n2 x y z
c                elem n1 n2 x y z
c                user n1 n2 x y z
c                cont n1 n2 x y z
c                arcl n1 n2
c                rsum n1 n2
c                sums n1 n2
c                ener
c                show

      elseif(j.eq.2) then

        call ptplot (hr(np(43)),mr(np(33)),ct,prt)

c     [para] - input parameters for data use

      elseif(j.eq.3) then

        coflg = .true.
        call pconst(prt)

c     [func,name] - execute function 'name'

      elseif(j.eq.4) then

        yyy = lct
        n = index(yyy,' ')
        if(n.gt.1) then
          yyy(n:n+3) = '.fcn'
          inquire(file=yyy,exist=errck)
          if(errck) then
            icf = max(icl,abs(ior))
 41         inquire(unit=icf,opened=errck)
            if(errck) then
              icf = icf + 1
              go to 41
            endif
            call opnfil('fcn',yyy,-2,icf,lread)
            if(.not.lread) then
              if(ior.lt.0) return
              call plstop()
            endif
            isfile(isf) = ior
            ior         = icf
            isf         = isf + 1
            call pconst(prt)
            close(icf)
            isf         = max(1,isf-1)
            ior         = isfile(isf)
            icf         = max(icl,abs(ior))
            lread       =.false.
          else
            write(*,*) 'FILE ',yyy,' DOES NOT EXIST'
          end if
        else
          write(*,*) 'BLANK FILE NAME'
        end if

c     [dync,<off>] - dynamic contact updates

      elseif(j.eq.5) then

        if(pcomp(lct,'off',3)) then
          dyncon = .false.
        else
          dyncon = .true.
        endif

c     [part,,i] i=1 to npart (defined by mesh)

      elseif(j.eq.6) then

c       If monolithic done then recompute profile for first partition

        if(monofl) then
          monofl = .false.
          do n = 1,ndf
            ndfp(n) = ndfg(n)
          end do ! n

c         Set active partition 1 data

          do n = 1,4
            if(.not.tflp(n)) then
              npart = n
              call partpt(n,tflp(n),.false.)

c             Reset profile

              call profpart(mr(np(31)))
            endif
          end do ! n
        endif

c       Set active partition or report active partition

        i = nint(ct(1))
        if    (i.eq.0) then
          write(*,*) ' Current Partition =',npart
        elseif(i.lt.1 .or. i.gt.4) then
          if(ior.lt.0) then
            write(*,*) ' Bad partition number, reinput'
          else
            write(iow,*) ' Bad partition number, stop'
            call plstop()
          endif
        else
          errck = .false.
          call partpt(i,errck,.true.)
          if(npld.gt.0) then
            prop = propld(ttim,0)
          else
            prop = 1.0
          endif
          if(prt) then
            if(ior.lt.0) write(*,*) '   Partition =',i,
     &                              ': Prop. ld',prop
            write(iow,*) '   Partition =',i,': Prop. ld',prop
          endif
        endif

c     [mate,,i] i=0 to nummat (defined by mesh)

      elseif(j.eq.7) then

        msopr = msplt
        msplt = max(0,min(nint(ct(1)),nummat))
        if(msplt.ne.msopr) then
          if(ior.lt.0) then
            write(*,2003) msplt
          else
            write(iow,2003) msplt
          endif
          fl(11) = .false.
          rfl    = .false.
        endif

c     [capt,label]  Caption for next contour plot

      elseif(j.eq.8) then

        ncapt   = 1
        caption = lct

c     [mono]lithic -- no partitions

      elseif(j.eq.9) then

        monofl = .true.

c       Set active monolithic partition data

        call partpt(5,tflp(5),.false.)

        do i = 1,ndf
          ndfp(i) = 1
        end do ! i
        npart  =  1
        call profpart(mr(np(31)))

      endif

c     Formats

2002  format('   O u t p u t    L i s t',i5)

2003  format('     Plot results: Material set number    ',i3/
     &       '                   (0 - for all sets)')

      end
