
                             ! MAIN CONTACT COMMON

      integer  cck,numcs,numcm,numcp,ofssurf,ofsmate,isgp1,isgp3,indb
      logical  ifct,ifdb,ifistgn,ifchist,iffron,lagrm
      common /c_contac/
     &  cck(c_ncc),          ! Contact Commands Kounter (# of commands read)
     &  numcs,               ! # Contact Surfaces  (= cck(1))
     &  numcm,               ! # Contact Materials (= cck(2))
     &  numcp,               ! # Contact Pairs     (= cck(3))
     &  ofssurf,             ! counting offset for surface data
     &  ofsmate,             ! counting offset for material data
     &  isgp1,               ! CH1 & CH2 counting offset
     &  isgp3,               ! CH3       counting offset
     &  indb,                ! Debug level
     &  ifct,                ! ConTact Flag, if true contact is active
                             ! set by CONT,ON - CONT,OFF
     &  ifdb,                ! # DeBug Flag, if true debug is active
                             ! set by CONT,DEBU, output on fort.98, fort.99
     &  ifistgn,             ! Geometrical check flag
     &  ifchist,             ! CHange of geometry status flag
     &  iffron,              ! FRiction ON flag
     &  lagrm                ! LAGRange Multiplier flag
