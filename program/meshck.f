c$Id: meshck.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      subroutine meshck(ip,ie,iedof,id,nty,ix,nie,nen,nen1,ndf,
     &                  numnp,numel,nummat,errs)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform check on mesh data to ensure nodes/elements input

c      Inputs:
c         ie(nie,*)      - Material set assembly information
c         iedof(ndf,*,*) - Material set nodal assembly information
c         id(ndf,*)      - Boundary condition and equation number array
c         nty(*)         - Nodal type
c         ix(nen1,*)     - Element nodal connection lists
c         nie            - Dimension of ie array
c         nen            - Maximum number of nodes/element
c         nen1           - Dimension for ix array
c         ndf            - Number dof/node
c         numnp          - Number of nodes in mesh
c         numel          - Number of elemenst in mesh
c         nummat         - Number of material sets

c      Outputs:
c         ip(ndf,*)      - List of active nodes, used for graphics
c         errs           - Flag, true if errors detected
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'cblend.h'
      include  'corset.h'
      include  'elflag.h'
      include  'iofile.h'
      include  'mdata.h'
      include  'prflag.h'
      include  'pointer.h'
      include  'comblk.h'

      logical   errs
      integer   i,ii,iel,j,ma,mg,n,nen,nen1,nie,ndf,numnp,numel,nummat
      integer   ip(ndf,*),ie(nie,*),iedof(ndf,nen,*),id(ndf,*)
      integer   ix(nen1,*),nty(*)

      save

c     Perform mesh checks to ensure nodes/elements input

      errs = .false.
      do n = 1,numel
        if (ix(nen1,n).le.0 .or. ix(nen1,n).gt.nummat) then
          write(ilg,2000) n
          write(iow,2000) n
          if(ior.lt.0) write(*,2000) n
          errs = .true.
        else
          do i = 1,nen
            ii = ix(i,n)
            if(ii.gt.numnp .or. ii.lt.0) then
              write(ilg,2001) ii,n
              write(iow,2001) ii,n
              if(ior.lt.0) write(*,2001) ii,n
              errs = .true.
            elseif(ii.ne.0 .and. nty(ii).lt.0) then
              write(ilg,2002) ii,n
              write(iow,2002) ii,n
              if(ior.lt.0) write(*,2002) ii,n
              errs = .true.
            endif
          end do ! i
        endif
      end do ! n

c     Remove unused dof's using ie(nie,*) array

      do n = 1,numnp
        do j = 1,ndf
          ip(j,n) = 0
        end do ! j
      end do ! n

c     Check nodes on each element for active dof's

      do n = 1,numel
        mg = ix(nen1,n)

c       Loop over material sets

        do ma = 1,nummat
          if(ie(nie-2,ma).eq.mg) then
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
      end do ! n

c     Check for point mass, dampers, stiffness

      if(nmfl) then
        call mshckpt(hr(np(86)),hr(np(87)),hr(np(88)),ip,ndf,numnp)
      endif

c     Set b.c. restraints for unused dof's

      do n = 1,numnp
        do j = 1,ndf
          if(ip(j,n).eq.0) then
            id(j,n) = -1000
          else
            id(j,n) = abs(id(j,n))
          end if
        end do ! j
      end do ! n

c     Remove unused nodes - for graphics

      do n = 1,numnp
        ip(1,n) = 0
      end do ! n

      do n = 1,numel
        do i = 1,nen
          ii = ix(i,n)
          if(ii.gt.0) ip(1,ii) = 1
        end do ! i
      end do ! n

c     Set flat to indicate node is not used

      do n = 1,numnp
        if(ip(1,n) .eq. 0) then
          nty(n) = -1
        end if
      end do ! n

c     Fix all unspecified coordinate dof's

      do n = 1,numnp
        if(nty(n).lt.0) then
          do i = 1,ndf
            id(i,n) = -1001
          end do ! i
        endif
      end do ! n

c     If supernodes used then

      if(numbd.gt.0) then
        if(numsn.gt.0) then
          if(numsd.gt.0) then
            call mshcksn(mr(np(162)),mr(np(164)),numsd,numsn,numbd,errs)
          else
            write(ilg,2003)
            write(iow,2003)
            if(ior.lt.0) write(*,2003)
            errs = .true.
          endif
        else
          write(ilg,2004)
          write(iow,2004)
          if(ior.lt.0) write(*,2004)
          errs = .true.
        endif
      endif

c     Set first and last element numbers for each material type
c     N.B. Limit set by dimension in include file elflag.h

      do ma = 1,min(80,nummat)
        do n = 1,numel
          if(ix(nen1,n).eq.ma) then
            elstart(ma) = n
            go to 100
          endif
        end do ! n
100     do n = numel,1,-1
          if(ix(nen1,n).eq.ma) then
            ellast(ma) = n
            go to 200
          endif
        end do ! n
200     continue
      end do ! ma

c     Check for rotation values

      do ma = 1,3
        dal(ma) = 0
        ral(ma) = 0
      end do ! ma

c     Check for angle condtions

      do ma = 1,nummat
        iel = ie(nie-1,ma)
        if(anglefl) then
          if(iel.gt.0) then
            do i = 1,2
              if(dal(i).eq.ia(i,iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ia(i,iel)
              else
                write(iow,2005) iel,i,ia(i,iel),dal(i)
                write(ilg,2005) iel,i,ia(i,iel),dal(i)
              endif
            end do ! i
            do i = 1,2
              if(ral(i).eq.ir(i,iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = ir(i,iel)
              else
                write(iow,2005) iel,i,ir(i,iel),ral(i)
                write(ilg,2005) iel,i,ir(i,iel),ral(i)
              endif
            end do ! i
          elseif(iel.lt.0) then
            do i = 1,2
              if(dal(i).eq.ea(i,-iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ea(i,-iel)
              else
                write(iow,2005) iel,i,ea(i,-iel),dal(i)
                write(ilg,2005) iel,i,ea(i,-iel),dal(i)
              endif
            end do ! i
            do i = 1,2
              if(ral(i).eq.er(i,-iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = er(i,-iel)
              else
                write(iow,2005) iel,i,er(i,-iel),ral(i)
                write(ilg,2005) iel,i,er(i,-iel),ral(i)
              endif
            end do ! i
          endif
        endif
        if(eulerfl) then
          if(iel.gt.0) then
            do i = 1,3
              if(dal(i).eq.ia3(i,iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ia3(i,iel)
              else
                write(iow,2005) i,iel,ia3(i,iel),dal(i)
                write(ilg,2005) i,iel,ia3(i,iel),dal(i)
              endif
            end do ! i
            do i = 1,3
              if(ral(i).eq.ir3(i,iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = ir3(i,iel)
              else
                write(iow,2005) iel,i,ir3(i,iel),ral(i)
                write(ilg,2005) iel,i,ir3(i,iel),ral(i)
              endif
            end do ! i
          elseif(iel.lt.0) then
            do i = 1,3
              if(dal(i).eq.ea3(i,-iel)) then

              elseif(dal(i).eq.0) then
                dal(i) = ea3(i,-iel)
              else
                write(iow,2005) iel,i,ea3(i,-iel),dal(i)
                write(ilg,2005) iel,i,ea3(i,-iel),dal(i)
              endif
            end do ! i
            do i = 1,3
              if(ral(i).eq.er3(i,-iel)) then

              elseif(ral(i).eq.0) then
                ral(i) = er3(i,-iel)
              else
                write(iow,2005) iel,i,er3(i,-iel),ral(i)
                write(ilg,2005) iel,i,er3(i,-iel),ral(i)
              endif
            end do ! i
          endif
        endif
      end do ! ma

c     Formats

2000  format(10x,' *ERROR* MESHCK: Data for element ',i10,' not input')

2001  format(10x,' *ERROR* MESHCK: Data for node ',i10,' on element',
     &       i10,' greater than maximum or negative')

2002  format(10x,' *ERROR* MESHCK: Data for node ',i10,' on element',
     &       i10,' not input')

2003  format(10x,' *ERROR* MESHCK: Blending functions used but no',
     &           ' SIDEs exist')

2004  format(10x,' *ERROR* MESHCK: Blending functions used but no',
     &           ' SNODes exist')

2005  format(10x,' *ERROR* MESHCK: Rotation of dof incompatible'/
     &       10x,'         ELEMENT TYPE =',i3,' DOF =',i3,' Current =',
     &           i3,' Previous =',i3) 
      end
