c$Id: urotrd.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine urotrd(urotyp, id, ndf )

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: 5/6 Dof constraint check for shells (dummy)

c      Inputs:
c         urotyp - User rotation type (nn matches urotnn)
c         ndf    - Number dof/node

c      Outputs:
c         id(*)  - Constraint for 5/6 dof

c      Remark:
c         rotyp map:  5 = user defined (not checked correctly)
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   urotyp, ndf, id(ndf)

      save

c     id(1) = 0

      end
