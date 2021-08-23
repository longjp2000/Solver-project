
                   !  CONTACT GEOMETRY COMMON (for surfaces selected in pair)

      integer         nsurf1,ifsty1,nope1,dnope1,neps1,ofs1,
     &                nsurf2,ifsty2,nope2,dnope2,neps2,ofs2,cndm
      common /c_geom/
     &  nsurf1,    !  # of SURFace 1 (slave)
     &  ifsty1,    !  Selection of Surface TYpe for surf. 1
     &  nope1,     !  # of NOdes Per Elements   for surf. 1
     &  dnope1,    !  Dimension NOdes Per Elemt for surf. 1
     &  neps1,     !  # of Elements Per Surface for surf. 1
     &  ofs1,      !  offset for data in ICS0   for surf. 1
     &  nsurf2,    !  # of SURFace 2 (master)
     &  ifsty2,    !  Selection of Surface TYpe for surf. 2
     &  nope2,     !  # of NOdes Per Elements   for surf. 2
     &  dnope2,    !  Dimension NOdes Per Elemt for surf. 2
     &  neps2,     !  # of Elements Per Surface for surf. 2
     &  ofs2,      !  offset for data in ICS0   for surf. 2
     &  cndm       !  Contact surface dimension
