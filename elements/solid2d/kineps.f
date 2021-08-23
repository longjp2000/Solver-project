      subroutine kineps(f,fi,df,detf,f33n,f33, lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute kinematic quantities for plane stress

c      Inputs:
c         f33n(*)       - Plane stress deformation gradient   (t_n)
c         f33 (*)       - Plane stress deformation gradient   (t_n+1)
c         lint          - Number of quadrature points
c         f(3,3,2,*)    - Deformation gradient                (no F_33)
c         fi(3,3,*)     - Inverse deformation gradient        (no F_33)
c         df(3,3,*)     - Incremental deformation gradient    (no F_33)
c         detf(2,*)     - Determinant of deformation gradient (no F_33)

c      Outputs:
c         f(3,3,2,*)    - Deformation gradient
c         fi(3,3,*)     - Inverse deformation gradient
c         df(3,3,*)     - Incremental deformation gradient
c         detf(2,*)     - Determinant of deformation gradient
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   lint, l
      real*8    f(3,3,2,*),fi(3,3,*),df(3,3,*),detf(2,*),f33n(*),f33(*)

      do l = 1,lint
        f(3,3,1,l) = f33 (l) + 1.0d0
        f(3,3,2,l) = f33n(l) + 1.0d0
        fi(3,3,l)  = 1.d0/f(3,3,1,l)
        df(3,3,l)  = f33 (l) - f33n(l)
        detf(1,l)  = detf(1,l)*f(3,3,1,l)
        detf(2,l)  = detf(2,l)*f(3,3,2,l)
      end do ! l

      end

