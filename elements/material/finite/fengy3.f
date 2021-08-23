c$Id: fengy3.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine fengy3(d,detf,u,up,upp, ha,hp,hpp, jsw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute pressure-entropy function

c     Inputs:
c       d(*)        - Material parameter array
c       detf        - Deformation gradient
c       jsw         - Function type: 1 - K*[ 0.25*( J^2 - 1 )-0.5*ln(J)]
c                                    2 - K*[ 0.5*( J - 1 )^2 ]
c                                    3 - K*[ 0.5*( ln(J) )^2 ]

c     Outputs:
c       u           - Internal energy
c       up          - First derivative of internal energy
c       upp         - Second derivative of internal energy
c       ha          - Augmentation function
c       hp          - First derivative of augmentation function
c       hpp         - Second derivative of augmentation function
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'augdat.h'

      integer   jsw
      real*8    d1, detf, u, up, upp, ha, hp, hpp, d(*)

c     Perform (augmented) iteration on penalty

      d1    = augf*d(21)

c     Free energy function
c     Psi   = U(J) -3*alpha*J*U'(J)*[T - T_ref]

c     up    = partial_J ( Psi )
c     upp   = [ partial_J ( J partial_J Psi ) - up ]/J
c           =   partial^2_J ( Psi )/J

c     Current volumetric functions are:

c     Model 1.) U(J) = K*0.25*(J^2 - 1 - 2*(log J))

      if    (jsw.eq.1) then

        u    = d1*( detf**2 - 1.d0  - 2.d0*log(abs(detf)))*0.25d0
        up   = d1*( detf - 1.d0/detf    )*0.5d0
        upp  = d1*( 1.d0 + 1.d0/detf**2 )*0.5d0

c     Model 2.) U(J) = K*0.5*(J-1)^2

      elseif(jsw.eq.2) then

        u    = d1*(detf - 1.d0)**2*0.5d0
        up   = d1*(detf - 1.d0)
        upp  = d1

c     Model 3.) U(J) = K*0.5*(log J)^2

      elseif(jsw.eq.3) then

c       up   = ( d1*log( detf ) - 3*K*alpha*(T - Tref) )/detf

        u    = d1*log(abs(detf))**2*0.5d0
        up   = d1*log(abs(detf)) / detf
        upp  = ( d1/detf  - up )/detf

      endif

c     Augmented Lagrangian function and derivatives
c     Current augmented Lagrangian function is

c     Model 1.) h(J) = (J - 1)

      ha  = detf - 1.d0
      hp  = 1.0d0
      hpp = 0.0d0

      end
