c$Id: umod1d.f,v 1.1 2006/11/20 20:33:27 rlt Exp $
      subroutine umod1d(umat,eps,td,d,ud,hn,h1,nh,ii,istrt, sig,dd, isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: User Constitutive Model

c     Input:
c          umat    -  User material type
c          eps     -  Current strain at point
c          td      -  Temperature change
c          d(*)    -  System material parameters
c          ud(*)   -  User material parameters
c          hn(nh)  -  History terms at point: t_n
c          h1(nh)  -  History terms at point: t_n+1
c          nh      -  Number of history terms
c          ii      -  Number of calls to constitution/element
c          istrt   -  Start state: 0 = elastic; 1 = last solution
c          isw     -  Element control parameter

c     Output:
c          sig     -  Stress at point.
c          dd      -  Current material tangent modulus

c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   umat,nh,ii,istrt, isw
      real*8    td, eps,sig,dd, d(*),ud(*),hn(nh),h1(nh)

      save

c     Material Model 1

      if(umat.eq.1) then

c       Dummy elastic model:  sig = E*eps

        dd  = d(1)
        sig = d(1)*eps

      endif

      end
