c$Id: uformfe.f,v 1.1 2006/11/20 20:33:27 rlt Exp $
      subroutine uformfe(pnu,pna,pnl,pnb,aufl,alfl,bfl,dfl,
     &                   isw,nl1,nl2,nl3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:  Routine for constructing user contributions

c     Inputs:
c        pnu   - Pointer to current solution vector
c        pna   - Pointer to LHS diagonal (and upper) array
c        pnl   - Pointer to LHS lower array
c        pnb   - Pointer to RHS vector
c        aufl  - Logical flag to assemble LHS upper array
c        alfl  - Logical flag to assemble LHS lower array
c        bfl   - Logical flag to assembel RHS (active equations only)
c        dfl   - Logical flag to form RHS in full nodal form
c        isw   - Switch parameter
c        nl1   - First element to process
c        nl2   - Last  element to process
c        nl3   - Increment

c     Outputs:
c        Arrays in pointer locations 
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      logical   aufl,bfl,alfl,dfl
      integer   pnu,pnb,pna,pnl,isw,nl1,nl2,nl3

      save

      end
