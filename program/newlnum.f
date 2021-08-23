c$Id: newlnum.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine newlnum(lagre,lagrn,id,ix,ren,irb,jnt)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Compute equation structure for lagrange multiplier unknowns

c      Inputs:
c        id(ndf,*)     - Solid element equation numbers
c        ix(nen1,*)    - Element connection array
c        ren(*)        - Nodal reorder list
c        irb(nrbdof,*) - Rigid body equation numbers
c        jnt(6,numjts) - Joint equation numbers

c      Outputs
c        lagre(numel)    - First element equation number
c        lagrn(numel)    - Number equations after node number
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cdata.h'
      include  'part0.h'
      include  'part1.h'
      include  'rigid1.h'
      include  'rjoint.h'
      include  'sdata.h'

      integer   eqad,neqad,i,m,mm,n,nn
      integer   lagre(numel),lagrn(numnp),id(ndf,numnp),ix(nen1,*)
      integer   ren(numnp),irb(nrbdof,nrbody),jnt(6,numjts)

      save

c     Zero nodal array

      do nn = 1,numnp
        lagrn(nn) = 0
      end do ! nn

c     Search elements

      neqad = 0
      eqad  = 0
      do n = 1,numel
        if(lagre(n).gt.0) then

c         Find maximum reordered node number on element

          nn = 0
          do i = 1,nen
            if(ix(i,n).gt.0) then
              nn = max(nn,ren(ix(i,n)))
            endif
          end do ! i

c         Adjust all other equations

          do m = 1,numnp
            mm = ren(m)
            if(mm.gt.nn) then
              do i = 1,ndf
                if(id(i,mm).gt.0 .and. ndfp(i).eq.npart) then
                  id(i,mm) = id(i,mm) + lagre(n)
                endif
              end do ! i
            else
              do i = 1,ndf
                if(id(i,mm).gt.0 .and. ndfp(i).eq.npart) then
                  eqad     = max(eqad,id(i,mm))
                endif
              end do ! i
            endif
          end do ! m

c         Check rigid body equation numbers

          if(rbody .and. nrbprt.eq.npart) then
            do m = 1,numjts
              if(jnt(6,m).gt.eqad) then
                jnt(6,m) = jnt(6,m) + lagre(n)
              end if
            end do ! m
            do m = 1,nrbody
              do i = 1,nrbdof
                if(irb(i,m) .gt. eqad) then
                  irb(i,m) = irb(i,m) + lagre(n)
                  neqr     = max(neqr,irb(i,m))
                endif
              end do ! i
            end do ! m
          endif

c         Set numbers for element

          mm        = lagre(n)
          if(lagrn(nn).eq.0) then
            lagre(n)  = eqad  + 1
          else
            lagre(n)  = lagrn(nn) + 1
          endif
          neqad     = neqad + mm
          eqad      = eqad  + mm
          lagrn(nn) = eqad

        endif
      end do ! n

c     Set final active equation number

      neq = neq + neqad
      if(.not.rbody .or. nrbprt.ne.npart) then
        neqr = neq
      endif

      end
