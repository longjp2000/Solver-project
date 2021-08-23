c$Id: iterat.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine iterat(pu,prsd,oldrsd,d,t,
     &                  accrcy,v,w,prt,id,nbfgs,stol, etol)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: BFGS algorithm  for solution to equation

c          R(u) = F - P(u) = 0
c          Variable number of vectors
c          Max. no. of vectors = nbfgs ( < = 15 ) (input data)

c     CALL sequence
c      pmacr3.f -> bfgs.f (iterat) -> formfe.f (isw = 6)
c      -> pform.f -> elmlib.f ->  elmtxx.f

c    Description
c      routine will return solution and residual to equation:
c          K(u)*du - b = 0
c      Routine evaluates RHS of equation through routine operat.

c      Inputs:
c       pu          - Pointer for solution vectors
c       prsd        - Pointer for current residual
c       oldrsd(*)   - Old residual
c       prt         - Print if true
c       id(*)       - Equation numbers for each dof
c       stol        - Solution tolerance for line search
c       etol        - Solution tolerance for energy norm

c      Outputs:
c       d(*)        - Incremental solution vector for step
c       v(*)        - BFGS vectors
c       w(*)        - BFGS vectors
c       nbfgs       - Number BFGS steps

c      Scratch:
c       t(*)        - Temporary storage vectors
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'counts.h'
      include  'ddata.h'
      include  'endata.h'
      include  'fdata.h'
      include  'iofile.h'
      include  'hdatam.h'
      include  'ndata.h'
      include  'pointer.h'
      include  'rdata.h'
      include  'rdat1.h'
      include  'sdata.h'
      include  'comblk.h'

      logical   accrcy,prt
      integer   i,nupd,nbfgs, pu,prsd, id(*)
      real*8    g,g0,s, onorm,dnorm,stol,etol
      real*8    dot, gamma1, oldrsd(*),d(*),v(*),w(*),t(*)

      save

c     Set Starting values

      iform = 0
      nupd  = 0
      s     = 1.0d0
      g0    = 0.0d0

c     Initialization of residual R(u+s*d) [set to rsd]

      call pzero(oldrsd,neq)
      call pzero(hr(prsd),neq)
      call pzero(     d,neq)

c     Compute residual for du = 0

      g = gamma1(id,pu,prsd,d,t,s)

c     Loop for momentum balance iteration

      do i = 1,nbfgs

        if (prt) then
          rnorm = sqrt ( dot (hr(prsd),hr(prsd),  neq))
          onorm = sqrt ( dot (oldrsd,oldrsd,  neq))
          dnorm = sqrt ( dot (     d,     d,  neq))
          write(iow,2000) i, rnorm, onorm, dnorm
          if (ior.lt.0) then
            write(iow,2000) i, rnorm, onorm, dnorm
          endif
        endif

c       Compute search direction (d) by factorized form (u,w)
c                          d : not including time step

        call dfind(d,hr(prsd),oldrsd,nupd,g0,g,s,neq,v,w,nbfgs)

c       Do line search if necessary [new residual in rsd]

        s  = 1.d0
        g  = gamma1(id,pu,prsd,d,t,s)
        g0 = dot(d,oldrsd,neq)
        write(iow,2003) i,g0,g
        if(ior.lt.0) then
          write(*,2003) i,g0,g
        endif

c       Compute step size using scalar line search and update solution

        if(abs(g).gt.stol*abs(g0)) then
          call serchl(g0,id,prsd,pu,d,stol,t,neq, s)

c       Call to update histories

        else
          hflgu  = .true.
          h3flgu = .true.
          g = gamma1(id,pu,prsd,d,t,s)
        endif

        call update(id,hr(np(30)),hr(pu),hr(np(42)),d,fl(9),2)

c       Convergence checks: Norm of residual vector and energy incr.

        rnorm = sqrt(dot(hr(prsd),hr(prsd),neq))
        aengy = abs(g0*s)
        if (rnmax.eq.0.0d0) then
          rnmax  = aengy
        endif

        write(iow,2001) rnorm,aengy,rnmax
        if(ior.lt.0) then
          write(*,2001) rnorm,aengy,rnmax
        endif

        accrcy = aengy.le.(etol*rnmax)

        if (accrcy) go to 100

      end do ! i

      i = nbfgs

100   write(iow,2002) i,iform
      if(ior.lt.0) then
        write(*,2002) i,iform
      endif
      nform = nform + iform
      iform = 0

c     Formats

2000  format(4x,'Start of iteration no. ',i5/
     &       6x,'|Residual|  |Old Residual|   |Search Vect|'/
     &       1p,3d16.5)

2001  format('   | Res | =',g13.5,' | Energy |=',g13.5,
     &         ' | Max Engy |=',g13.5)

2002  format(4x,i3,' BFGS iterations with',i3,' RHS forms')

2003  format(4x,'Iteration',i3,': G0 =',1p,1e15.7,', G =',1p,1e15.7)

      end
