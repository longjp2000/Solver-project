c$Id: dicont.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine dicont(id,numnp,ndf,lflag)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Provides information for displacement control
c               in arc length.

c      Programmed: 04/23/87 P. Wriggers

c      Inputs:
c         id(ndf,*) - Equation numbers for each dof
c         numnp     - number of nodal points in mesh
c         ndf       - Number dof/node
c         lflag     - If true, changes arc length

c      Outputs:
c         Equation number of assigned displacement
c         Factors to scale arc length control
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'arclel.h'
      include  'arclei.h'
      include  'arcler.h'
      include  'comfil.h'
      include  'ioincl.h'
      include  'iofile.h'

      character ch*1
      logical   errck, pinput
      integer   numnp,ndf,lflag
      integer   id(ndf,*)
      real*8    td(3)

      save

      if (lflag .ne. 0) go to 100

c     Read in if numerical damping desired or not

      if(ior.lt.0) then
        call pprint('   Input: numerical damping (1 = no damping)->')
      endif
      errck = pinput(td,1)
      ndamp = td(1)
      write (iow,2003) ndamp

c     Restart flag

      if (refl) go to 100
50    if (kflag.eq.4.or.kflag.eq.5) then
        if(ior.lt.0) then
          call pprint(
     &             '   Input: node nr., dof.nr.,  prescribed displ.-> ')
        endif
        errck = pinput(td,3)
        nodis = td(1)
        nddis = td(2)
        alfa0 = td(3)
        if(ior.lt.0) then
          if(nodis.le.0 .or. nodis.gt.numnp) go to 50
          if(nddis.le.0 .or. nddis.gt.ndf  ) go to 50
        endif

        ndis = id(nddis,nodis)
        write (iow,3000) nodis,nddis,alfa0
        if(ndis.le.0) then
          if(ior.lt.0) then
            write(*,2001)
            go to 50
          else
            write(iow,2001)
            call plstop()
          endif
        endif
      endif
      return

c     For restart only
c     Any method (displacement control stiff.param. just for chance)

 100  write(iow,2005) rlnew,c0,cs01,cs02

c     Arc length method (any)

      if (kflag.lt.4.or.kflag.eq.6) then
        if(ior.lt.0) then
          write(*,2006) ds0,r
          call pprint('   Keep arc-length and load-direction (y or n):')
          read (*,1000) ch
        else
          read (ior,1000,end=900) ch
          irecrd(isf) = irecrd(isf) + 1
          ch          = record(1:1)
        endif
        if(ch.eq.'n' .or. ch.eq.'N') then
          if(ior.lt.0) then
            call pprint(
     &         '   Input: new arc-length, new load direction->')
          endif
          errck = pinput(td,2)
          ds0 = td(1)
          r   = td(2)
          write(iow,2006) ds0,r
        endif
      endif

c     Displacement control

      if (kflag.eq.4.or.kflag.eq.5) then
        if(ior.lt.0) then
          write(*,2007) nodis,nddis,alfa0
          call pprint(
     &             '   Keep displacement control parameters (y or n): ')
          read (*,1000) ch
        else
          read (ior,1000,end=900) record
          irecrd(isf) = irecrd(isf) + 1
          ch          = record(1:1)
        endif
        if(ch.eq.'n' .or. ch.eq.'N') then
          go to 50
        endif
      endif
      return

c     Eof encountered

900   call  endclr ('DICONT',ch)

c     Formats

1000  format(a1)

2001  format('   Displacement control specified on restrained node')

2003  format('   Numerical damping = ',i3,3x,'(1 = no damping)')

2005  format('   v a l u e s  for  r e s t a r t:',/,
     & '     Current load level      = ',g12.4,/,
     & '     S t i f f n e s s  parameter values ',/,
     & '     Stiff.param first step  = ',g12.4,/,
     & '     Stiff.param 1.prev.step = ',g12.4,/,
     & '     Stiff.param 2.prev.step = ',g12.4,/)

2006  format('   Given arc length         = ',g12.4/
     *       '   Load direction           = ',f12.4)

2007  format('   Node number      = ',i3,/,
     &       '   Ndof number      = ',i3,/,
     &       '   Prescribed disp. = ',f10.3,/)

3000  format('   S i n g l e   D i s p l a c e m e n t   C o n t r o l '
     &       ,/,  '     Node     DOF   Displacement'
     &       ,/, i8,i8,1p,1e16.5)

      end
