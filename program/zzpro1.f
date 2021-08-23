c$Id: zzpro1.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine zzpro1(ix,ib,ip,
     &                  ma,ndm,nen,nen1,numnp,numel,isw, ipmax)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Driver routine for Zienkiewicz-Zhu stress projection

c      Inputs:

c      Outputs:
c-----[--+---------+---------+---------+---------+---------+---------+-]

c     Zienkiewicz-Zhu patch determination

      implicit   none

      include   'pointer.h'
      include   'comblk.h'

      integer    ndm, ma, nen,nen1,numnp,numel,isw, ipmax
      integer    ij, j, n, nel,nnl

      integer    ix(nen1,numel), ib(numnp), ip(numnp)

      save

c     1.) Identify active vertex nodes

c     Set boundary node indicator

      do n = 1,numnp
        ib(n) =  0
      end do ! n

c     Loop through elements to set up list

      do n = 1,numel
        if(ma.eq.0 .or. ix(nen1,n).eq.ma) then

c         Compute number of nodes on element

          nel = 0
          do ij = 1,nen
            if(ix(ij,n).gt.0) then
              nel = ij
            endif
          end do ! ij

c         Set the number of vertex nodes for each element

          if(ndm.eq.1) then                      ! Line/pt
            nnl = min(2,nel)
          elseif(ndm.eq.2) then                  ! Line/pt
            if(nel.le.2) then
              nnl = min(2,nel)
            elseif(nel.eq.3 .or. nel.eq.6) then  ! Tri
              nnl = 3
            else                                 ! Quad
              nnl = 4
            endif
          elseif(ndm.eq.3) then                  ! Line/pt
            if(nel.le.3) then
              nnl = min(2,nel)
            elseif(nel.eq.4 .or. nel.eq.10) then ! Tets
              nnl = 4
            else                                 ! Hex
              nnl = 8
            endif
          endif

c         Look up element nodes

          do ij = 1,nnl
            j = ix(ij,n)
            if(j.gt.0) then
              if(ib(j).eq.0) then
                ib(j) = 1
              else
                ib(j) = ib(j) + 1
              endif
            endif
          end do ! ij

          do ij = nnl+1,nel
            j = ix(ij,n)
            if(j.gt.0) then
              if(ib(j).eq.0) then
                ib(j) = -1
              else
                ib(j) = ib(j) - 1
              endif
            endif
          end do ! ij

        endif
      end do ! n

c     2.) Check for active nodes

      do n = 1,numnp
        ip(n) = 0
      end do ! n

      do n = 1,numnp
        if(ib(n).gt.0) then
          if(isw.eq.1) ip(n) = ib(n)
          ib(n) = 0
        else
          ib(n) = 1
        endif
      end do ! n

c     3.) Set pointers for patch elements

      do n = 2,numnp
        ip(n) = ip(n) + ip(n-1)
      end do ! n

      ipmax = ip(numnp)

      end
