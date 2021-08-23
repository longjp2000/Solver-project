c$Id: psetid.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine psetid(id,ix,ip,ie,iedof,nie,ndf,nen,nen1,
     &                  numel,numnp,nummat)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Remove unused dof's using ie(nie,*) array

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   i,ii,j,ma,n,nen,nen1,nie,ndf,numnp,numel,nummat
      integer   id(ndf,numnp,*),ie(nie,*),iedof(ndf,nen,*)
      integer   ip(ndf,*),ix(nen1,*)

      save

      do n = 1,numnp
        do j = 1,ndf
          ip(j,n)   = 0
          id(j,n,1) = id(j,n,2)
        end do ! j
      end do ! n

c     Check nodes on each element for active dof's

      do n = 1,numel

        if(ix(nen1-1,n).ge.0) then

c         Loop over the material sets

          do ma = 1,nummat
            if(ie(nie-2,ma).eq.ix(nen1,n)) then
              do i = 1,nen
                ii = ix(i,n)
                if(ii.gt.0) then
                  do j = 1,ndf
                    if(iedof(j,i,ma).gt.0) then
                      ip(iedof(j,i,ma),ii) = 1
                    endif
                  end do ! j
                endif
              end do ! i
            endif
          end do ! ma
        endif
      end do ! n

c     Set b.c. restraints for unused dof's

      do n = 1,numnp
        do j = 1,ndf
          if(ip(j,n).eq.0) then
            id(j,n,1) = 1
          end if
        end do ! j
      end do ! n

      end
