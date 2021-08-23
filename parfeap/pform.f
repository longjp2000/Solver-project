c$Id: pform.f,v 1.2 2006/11/25 02:06:29 rlt Exp $
      subroutine pform(ul,xl,tl,ld,p,s,ie,d,id,x,ix,rben,f,t,jp,
     &  u,ud,b,a,al,ndd,nie,ndf,ndm,nen1,nst,aufl,bfl,alfl,dfl,
     &  isw,nn1,nn2,nn3)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Remove interupt option                           11/24/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute element arrays and assemble global arrays

c      Inputs:
c         ie(nie,*)   - Element information for material set
c         d(ndd,*)    - Material set parameters
c         id(ndf,*)   - Equation numbers for each active dof
c         x(ndm,*)    - Nodal coordinates of mesh
c         ix(nen1,*)  - Element nodal connections of mesh
c         rben(*)     - Rigid/modal element indicator
c         f(ndf,*,2)  - Nodal force and displacement values
c         t(*)        - Nodal temperature values
c         jp(*)       - Pointer array for row/columns of tangent
c         u(*)        - Nodal solution values
c         ud(*)       - Nodal rate values
c         ndd         - Dimension for d array
c         nie         - Dimension for ie array
c         ndf         - Number dof/node
c         ndm         - Spatial dimension of mesh
c         nen1        - Dimension for ix array
c         nst         - Dimension for element array
c         aufl        - Flag, assemble coefficient array if true
c         bfl         - Flag, assemble vector if true
c         alfl        - Flag, coefficient array unsymmetric if true
c         dfl         - Flag, assemble reactions if true
c         isw         - Switch to control quantity computed
c         nn1         - First element number to process
c         nn2         - Last element number to process
c         nn3         - Increment to nn1

c      Local element arrays:
c         ul(*)       - Element solution and rate values
c         xl(*)       - Element nodal coordinates
c         tl(*)       - Element nodal temperatures
c         ld(nst,*)   - Element local/global equation numbers
c         p(nst,*)    - Element vector
c         s(nst,nst,*)- Element array

c      Outputs:
c         b(*)        - Global vector
c         a(*)        - Global matrix, diagonal and upper part
c         al(*)       - Global matrix, lower part
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'adapt1.h'
      include   'adapt2.h'
      include   'arcler.h'
      include   'cdata.h'
      include   'complx.h'
      include   'ddata.h'
      include   'elcount.h'
      include   'eldata.h'
      include   'elplot.h'
      include   'eluser.h'
      include   'erotas.h'
      include   'iofile.h'
      include   'hdata.h'
      include   'hdatam.h'
      include   'lmdata.h'
      include   'part0.h'
      include   'pbody.h'
      include   'pdata3.h'
      include   'pfeapb.h'
      include   'pointer.h'
      include   'prld1.h'
      include   'prlod.h'
      include   'prstrs.h'
      include   'ptdat1.h'
      include   'ptdat2.h'
      include   'ptdat8.h'
      include   'region.h'
      include   'rigid1.h'
      include   'rigid2.h'
      include   'setups.h'
      include   'tdatb.h'
      include   'comblk.h'

      include   'p_int.h'


      logical    aufl,bfl,alfl,dfl,efl,gfl,rel,rfl,erotflg,msflg,arotflg
      logical    mdfl, solsave
      character  y*1
      integer    isw, jsw, ksw
      integer    i, jj, nn1, nn2, nn3, nst, nneq, nov
      integer    ndf, ndm, nrot(2), ndd, nie, nen1
      integer    ld(nst,*),ie(nie,*),id(ndf,*),ix(nen1,*),rben(*),jp(*)
      real*8     xl(ndm,*),p(nst,*),s(nst,nst,*),d(ndd,*),ul(ndf,nen,*)
      real*8     x(ndm,*) ,f(ndf,numnp),u(ndf,*),ud(*),t(*),tl(*)
      real*8     un(20), dun(20), temp, prope, b(*), a(*), al(*)

      save

c     Save current equation solver type to permit stress projections

      solsave = solver

c     Set element proportional loading value

      prope = (theta(3)*(prop - propo) + propo)
      if(isw.ne.23) then
        prope = prope*rlnew
      endif

c     Recover nh1, nh2, nh3 pointers

      nh1 = np(50)
      nh2 = np(51)
      nh3 = np(52)

c     Set program and user material count parameters

      do i = 1,10
        nomats(1,i) = 0
        nomats(2,i) = 0
        unmats(1,i) = 0
        unmats(2,i) = 0
      end do ! i

c     Set parameters for solution

      iel = 0
      efl = .false.
      if(.not.dfl.and.isw.eq.6) efl = .true.
      if(bfl.and.isw.eq.3)      efl = .true.
      arotflg = aufl .or. bfl

      if(isw.eq.19) then
        if(bfl) efl = .true.
        jsw = 5
        ksw = 5
        gfl = .false.
      else
        jsw = isw
        ksw = 3
        gfl = .true.
      endif

      nneq   = numnp*ndf
      nrkn   = nrk*nneq - nneq
      nrcn   = nrc*nneq - nneq
      nrmn   = nrm*nneq - nneq

c     Check for a rigid body

      rfl = rbody .and. (nrbprt.eq.npart)

c     Loop over active elements

      do n = nn1,nn2,nn3

c       Check for active regions

        if((nreg.lt.0 .and. ix(nen1-1,n).ge.0)
     &                .or. (abs(ix(nen1-1,n)).eq.nreg)) then

c        Set up local arrays

         nov = 0
         do ma = 1, nummat

c         Check for rigid bodies

          if(rben(n).gt.0) then
            if(rbtype(rben(n)).eq.1) then
              if(jsw.eq.4 .or. jsw.eq.8) then
                erotflg = .true.
              elseif(jsw.eq.6) then
                do i = 1,nsplts
                  if(ispl(1,i).eq.n) then
                    erotflg = .true.
                  endif
                end do ! i
                do i = 1,nuplts
                  if(iupl(1,i).eq.n) then
                    erotflg = .true.
                  endif
                end do ! i
              else
                erotflg = .false.
              end if
            endif
          else
            erotflg = .false.
          endif

          if((.not.rbody .or. erotflg).or.(rfl.and.rben(n).le.0) ) then

            msflg = msplt.eq.0 .or. msplt.eq.ma

            if(ie(nie-2,ma).eq.ix(nen1,n)            .and.
     &        (jsw.ne.6 .or. (jsw.eq.6 .and. msflg)) .and.
     &        (jsw.ne.8 .or. (jsw.eq.8 .and. msflg))) then

c             Compute address and offset for history variables

              ht1 = np(49) + ix(nen+1,n) + ie(nie-3,ma)
              ht2 = np(49) + ix(nen+2,n) + ie(nie-3,ma)
              ht3 = np(49) + ix(nen+3,n) + ie(nie-4,ma)

c             If history variables exist move into nh1,nh2

              if(ie(nie,ma).gt.0) then
                do i = 0,ie(nie,ma)-1
                  hr(nh1+i) = hr(ht1+i)
                  hr(nh2+i) = hr(ht2+i)
                end do ! i
              endif

c             If Element variables exist move into nh3

              if(ie(nie-5,ma).gt.0) then
                do i = 0,ie(nie-5,ma)-1
                  hr(nh3+i) = hr(ht3+i)
                end do ! i
              endif

              if(ie(nie-1,ma).ne.iel) mct = 0
              iel   = ie(nie-1,ma)
              rotyp = ie(nie-6,ma)

c             Set local arrays for element

              fp(1) = ndf*nen*(ma-1) + np(240)        ! iedof
              call plocal(ld,id,ix(1,n),ie(1,ma),mr(fp(1)),
     &                    xl,ul,ule,tl,p(1,3),x,f,u,ud,t, un,dun,
     &                    nrot, rfl,rel, jsw)

              if(jsw.ne.8 .and. np(245).ne.0) then
                call pparlo(ld,mr(np(245)),ix(1,n),mr(np(31)+ndf*numnp))
              endif

c             Compute modal displacements for rigid body

              if(erotflg) then
                call moddis(xl,ul,ndm,ndf,nel,nen,hr(np(95)),
     &                      hr(np(104)),rben(n))
              end if

c             Set for projections

              if(jsw.eq.8) then
                erav = hr(np(207)+n-1)
                call pzero(hr(np(36)),nen*npstr)
                solver = .true.
              elseif(jsw.eq.11) then
                erav = 0.0d0
                call seserr(hr(np(36)),hr(np(58)+numnp),ix(1,n),
     &                      nen,npstr,numnp)
              else
                erav = 0.0d0
              endif

c             Call element library

              dm = prope
              call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,
     &                    ndf,ndm,nst,iel,jsw)

c             Numerical differentiation for tangent

              if(ndflg.or.ie(nie-7,ma).gt.0) then
                call ptdiff(n,0.0d0,.false.)
              endif

c             Store time history plot data from element

              if(jsw.eq.6) then

c               Standard element values

                do i = 1,nsplts
                  if(ispl(1,i).eq.n) then
                    jj = max(ispl(2,i),1)
                    spl(i) = tt(jj)
                  endif
                end do ! i

c               Standard user element values

                do i = 1,nuplts
                  if(iupl(1,i).eq.n) then
                    jj = max(iupl(2,i),1)
                    upl(i) = ut(jj)
                  endif
                end do ! i

              endif

c             Modify for rotated dof's

              if(nrot(1)+nrot(2).gt.0 .and. arotflg) then
                call ptlocal(ul,p,s, cplxfl, ndf,nen,nel,nst,nrot, 2)
              endif

c             Adaptive refinement quantity

              if(jsw.eq.11.and.em1.ne.0.d0) then
                if(ierr.eq.1) then
                  erav = sqrt(abs(verror))/em1
                elseif(em2.ne.0.0d0) then
                  erav = sqrt(abs(venere))/em2
                end if
                if(erav.ne.0.0d0) erav = heta/erav
                hr(np(207)+n-1) = erav
              end if

c             Position update terms 'ht1,ht2' from 'nh1,nh2' to save

              if(hflgu .and. ie(nie,ma).gt.0) then
                do i = 0,ie(nie,ma)-1
                  temp      = hr(ht1+i)
                  hr(ht1+i) = hr(nh1+i)
                  hr(nh1+i) = temp
                  temp      = hr(ht2+i)
                  hr(ht2+i) = hr(nh2+i)
                  hr(nh2+i) = temp
                end do ! i
              endif

c             Position update terms 'ht3' from 'nh3' to save

              if(h3flgu .and. ie(nie-5,ma).gt.0) then
                do i = 0,ie(nie-5,ma)-1
                  hr(ht3+i) = hr(nh3+i)
                end do ! i
              endif

c             Modify for non-zero displacement boundary conditions

              if(efl) then
                mdfl = .false.
                do i = 1,ndf
                  if(dun(i).gt.1.0d-10*un(i)) then
                    mdfl = .true.
                  endif
                end do ! i

                if(mdfl) then

c                 Get current element tangent matrix

                  if (.not.aufl) then
                    dm = prop
                    call elmlib(d(1,ma),ul,xl,ix(1,n),tl,s,p,
     &                          ndf,ndm,nst,iel,ksw)
                    if(nrot(1)+nrot(2).gt.0) then
                      call ptlocal(ul,p,s, .false., ndf,nen,nel,nst,
     &                             nrot, 2)
                    endif
                  end if

c                 Modify for displacements (non-zero Real part only)

                  do i = 1,nst
                    p(i,3) = p(i,3)*cc3
                  end do ! i
                  call pmodify(ld,p,s,p(1,3),nst)
                elseif(jsw.eq.3 .and. pfeap_blk) then
                  call smodify(ld,p,s,nst, .true.)
                elseif(jsw.eq.6 .and. pfeap_blk) then
                  call smodify(ld,p,s,nst, .false.)
                end if
              elseif(jsw.eq.3 .and. pfeap_blk) then
                call smodify(ld,p,s,nst, .true.)
              end if

c             Assemble arrays as necessary

              if(jsw.ne.14) then
                if (jsw.ne.8 .and. np(245).ne.0) then
                  call pcompress(ld, ld(1,7),p,s,nst)
                endif
                if(pfeap_on) then
                  if(jsw.eq.5) then               ! Mass assembly
                    call uasblem(s,p,ld,nst,aufl,bfl, 1)
                  else                            ! Other assembly
                    call passble(s,p,ld,ix(1,n), jp,a,al,b,
     &                           alfl,aufl,bfl,dfl,gfl,rel,nst,nov,jsw)
                  endif
                else
                  call passble(s,p,ld,ix(1,n), jp,a,al,b,
     &                         alfl,aufl,bfl,dfl,gfl,rel,nst,nov,jsw)
                endif
              endif
            end if

          end if

         end do ! ma

        end if ! regions

      end do ! n

c     Recover active solver type

      solver = solsave

c     Assemble mass for parallel form

      if(pfeap_on .and. jsw.eq.5) then
        call uasblem(s,p,ld,nst,aufl,bfl, 2)
      endif

c     Format statements

 1000 format(a)
 3000 format(' Interupt: Press "y" to interrupt or return to continue')

      end
