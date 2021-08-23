c$Id: expo44.f,v 1.1 2006/11/20 20:32:56 rlt Exp $
        subroutine sh3flmda ( v , xlm )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c        Description:    SH3FLMDA is a subroutine which computes a unique
c                        orthogonal transformation matrix which rotates
c                        the 3-rd canonical basis vector, E3 into the
c                        vector v --- without drill.
c                        The singularity of the formula at (or near) the
c                        case: v = -E3 is avoided by computing the matrix
c                        which transforms -E3 into v when v.E3 < 0, and
c                        then applying a 180-degree rotation to it.

c        Authors:        M.S. Rifai, J.C. Simo, & D.D. Fox.

c        Date:           January 1991.

c       Version:        This routine was tested in FEAP version 6.3
c-----[--.----+----.----+----.-----------------------------------------]
c        Routine Input:
c        --------------
c        v ............. Arbitrary vector in R-3.

c        Routine Output:
c        ---------------
c        xlm ........... Orthogonal transformation matrix, which
c                        rotates E3 into v without drill.
c-----[--.----+----.----+----.-----------------------------------------]
        implicit  none

        real*8    fac, v(3) , xlm(3,2)

c       Compute [XLm], when v.E3 > 0

        if (v(3).gt.0.d0) then
           fac      =   1.d0 / ( 1.d0 + v(3) )
           xlm(1,1) =   v(3) + fac * v(2) * v(2)
           xlm(1,2) =        - fac * v(2) * v(1)
           xlm(2,1) =        - fac * v(1) * v(2)
           xlm(2,2) =   v(3) + fac * v(1) * v(1)
           xlm(3,1) = - v(1)
           xlm(3,2) = - v(2)

c       Compute [XLm], when v.E3 < 0

        else
           fac      =   1.d0 / ( 1.d0 - v(3) )
           xlm(1,1) = - v(3) + fac * v(2) * v(2)
           xlm(1,2) =          fac * v(2) * v(1)
           xlm(2,1) =        - fac * v(1) * v(2)
           xlm(2,2) =   v(3) - fac * v(1) * v(1)
           xlm(3,1) =   v(1)
           xlm(3,2) = - v(2)
        endif

        end
