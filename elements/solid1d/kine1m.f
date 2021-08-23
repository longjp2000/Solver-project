c$Id: kine1m.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
      subroutine kine1m(shp,ul,f,finv,df,ndf,nel,nen,detf,lint)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: Compute deformation gradient and its inverse at tn+1

c     Inputs:
c        shp(2,4,*)  - Shape functions and derivatives at gauss points
c        ul(ndf,*)   - Nodal solution parameters
c        ul(2,*)     - Nodal stress free reference displacements
c        ndf         - Number dof/node
c        nel         - Number nodes/element
c        nen         - Maximum number nodes/element
c        lint        - Number of quadrature points

c     Outputs:
c        f(9,*)      - Deformation gradient at gauss points
c        finv(9,*)   - Inverse deformation gradient at points
c        df(9,*)     - Incremental deformation gradient at points
c        detf(*)     - Determinant of deformation gradient at points
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndf, nel,nen, lint, k, l
      real*8    shp(2,4,*),ul(ndf,nen,*)
      real*8    df(9,*),f(9,2,*),finv(9,*),detf(2,*),dfi(9)

c     Compute deformation gradient at t-n+1

c        F**-1 = I - grad u

      do l = 1,lint
        do k = 1,8
          f (k,1,l) = 0.0d0
          f (k,2,l) = 0.0d0
          finv(k,l) = 0.0d0
          df  (k,l) = 0.0d0
        end do ! k
        f (5,1,l) = 1.0d0
        f (9,1,l) = 1.0d0
        f (5,2,l) = 1.0d0
        f (9,2,l) = 1.0d0
        finv(5,l) = 1.0d0
        finv(9,l) = 1.0d0
      end do ! l

      do l = 1,lint
        do k = 1,nel
          finv(1,l) = finv(1,l) - ul(1,k,1)*shp(1,k,l)
        end do ! k
        finv(1,l) = finv(1,l) + 1.0d0
      end do ! l

c     F = ( F^-1)^-1

      do l = 1,lint
        detf(1,l) = 1.d0/finv(1,l)
        f(1,1,l) = detf(1,l)
      end do ! l

c     Compute incremental deformation gradient

      do l = 1,lint
        dfi(1) = 0.0d0
        do k = 1,nel
          dfi(1)   = dfi(1) + ul(1,k,2)*shp(1,k,l)
        end do ! k
        df(1,l) = dfi(1)*f(1,1,l)

c       Compute deformation gradient F_n

        f(1,2,l)  = f(1,1,l) - df(1,l)
        detf(2,l) = f(1,2,l)

      end do ! l

      end
