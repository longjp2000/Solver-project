c$Id: crfunc.f,v 1.1 2006/11/20 20:32:46 rlt Exp $
      subroutine crfunc (nsopt,td,emax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Coded/modified by:                           Date:          rel.:
c               Robert L. Taylor        10 Novermber 2001            1.0

c      Acronym: Contact Read FUNCtions

c      Purpose: Input of contact surface descriptions as functions

c      Inputs:
c         nsopt   - Function number
c         td(14)  - Parameters for function

c      Outputs:
c         emax    - Element max number found
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'
      include  'print.h'

      integer   nsopt,emax, i,ii
      real*8    td(14)

      save

c     Check sub-command option and subcommand data not needed

      call cdebug0 ('      crfunc',-1)

c     Set defaults

      emax  = 0

c     Cylinder
      if    (nsopt.eq.1) then
        write(iow,2001) nint(td(1)),td(2),td(3),nint(td(4))
c     Sphere
      elseif(nsopt.eq.2) then
        write(iow,2002) nint(td(1)),td(2),td(3),nint(td(4))
c     Cartesian plane
      elseif(nsopt.eq.3) then
        write(iow,2003) nint(td(1)),td(2),td(3),nint(td(4))
c     Normal plane
      elseif(nsopt.eq.4) then
        write(iow,2004) nint(td(1)),td(2),td(3),td(4),td(5),
     &                  td(6),nint(td(7))
c     Polynomial surface
      elseif(nsopt.eq.5) then
        do i = 7,14
          if(td(i).ne.0.0d0) ii = i
        end do ! i
        write(iow,2005) (td(i),i=1,3),(nint(td(i)),i=4,6),
     &                  (i-6,td(i),i=7,ii)
      endif

c     Formats

2001  format(/7x,'Cylindrical Surface: Radial Expansion'/
     &       10x,'Motion Direction            =',i8/
     &       10x,'Initial Radius              =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8/)

2002  format(/7x,'Spherical Surface: Radial Expansion'/
     &       10x,'Motion Direction            =',i8/
     &       10x,'Initial Radius              =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8/)

2003  format(/7x,'Cartesian Surface: Radial Expansion'/
     &       10x,'Coordinate Direction        =',i8/
     &       10x,'Initial Coordinate          =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8/)

2004  format(/7x,'Plane Surface: Normal Motion'/
     &       10x,'Motion Direction            =',i8/
     &       10x,'Normal Vector (n_1)         =',1p,1e12.4/
     &       10x,'Normal Vector (n_2)         =',1p,1e12.4/
     &       10x,'Normal Vector (n_3)         =',1p,1e12.4/
     &       10x,'Initial Coordinate          =',1p,1e12.4/
     &       10x,'Displacement value (0: fix) =',1p,1e12.4/
     &       10x,'Expansion proportional load =',i8/)

2005  format(/7x,'Polynomial Surface:'/
     &       10x,'Displacement Vector (u_1)   =',1p,1e12.4/
     &       10x,'Displacement Vector (u_2)   =',1p,1e12.4/
     &       10x,'Displacement Vector (u_3)   =',1p,1e12.4/
     &       10x,'Proportional load   (p_1)   =',i8/
     &       10x,'Proportional load   (p_2)   =',i8/
     &       10x,'Proportional load   (p_3)   =',i8/
     &      (10x,'Polynomial term ',i3,'       =',1p,1e12.4))

      end
