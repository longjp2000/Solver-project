c$Id: palloc.f,v 1.1 2006/11/20 20:32:37 rlt Exp $
      logical function palloc(num,name,length,precis)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2006: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c       1. Add Parallel arrays                              02/11/2006
c       2. Move 'NFORC' from position 244 to 256            02/11/2006
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Define, delete, or resize a dictionary entry.
c               Pointer defined for integer (single) and real
c               (double precision arrays).

c      Inputs:
c         num        - Entry number for array (see below)
c         name       - Name of array          (see below)
c         length     - Length of array defined: =0 for delete
c         precis     - Precision of array: 1 = integers; 2 = reals

c      Outputs:
c         np(num)    - Pointer to first word of array in memory.
c-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      integer    list
      parameter (list = 256)

      include   'allotd.h'
      include   'allotn.h'
      include   'pointer.h'

      logical    ualloc
      character  name*(*)
      integer    i, num,length,precis

      save

c     Set active arrays for FEAP into list

      data   (nlist(i),i=1,list)/

     &         'TANG1', 'TANG2', 'TANG3', 'TANG4', 'UTAN1', 'UTAN2',
     &         'UTAN3', 'UTAN4', 'CMAS1', 'CMAS2', 'CMAS3', 'CMAS4',

     &         'LMAS1', 'LMAS2', 'LMAS3', 'LMAS4', 'DAMP1', 'DAMP2',
     &         'DAMP3', 'DAMP4', 'JP1  ', 'JP2  ', 'JP3  ', 'JP4  ',

c     Solution arrays

c              'TANG1',     !     1: Symm. tangent, partition 1        (na)
c              'TANG2',     !     2: Symm. tangent, partition 2
c              'TANG3',     !     3: Symm. tangent, partition 3
c              'TANG4',     !     4: Symm. tangent, partition 4

c              'UTAN1',     !     5: Unsym tangent, partition 1        (nul)
c              'UTAN2',     !     6: Unsym tangent, partition 2
c              'UTAN3',     !     7: Unsym tangent, partition 3
c              'UTAN4',     !     8: Unsym tangent, partition 4

c              'CMAS1',     !     9: Consist. mass, partition 1        (nm)
c              'CMAS2',     !    10: Consist. mass, partition 2
c              'CMAS3',     !    11: Consist. mass, partition 3
c              'CMAS4',     !    12: Consist. mass, partition 4

c              'LMAS1',     !    13: Diagonal mass, partition 1        (nl)
c              'LMAS2',     !    14: Diagonal mass, partition 2
c              'LMAS3',     !    15: Diagonal mass, partition 3
c              'LMAS4',     !    16: Diagonal mass, partition 4

c              'DAMP1',     !    17: Symm. damping, partition 1        (nc)
c              'DAMP2',     !    18: Symm. damping, partition 2
c              'DAMP3',     !    19: Symm. damping, partition 3
c              'DAMP4',     !    20: Symm. damping, partition 4

c              'JP1  ',     !    21: Profile pointer, partition 1      (n12p)
c              'JP2  ',     !    22: Profile pointer, partition 2      (n12p)
c              'JP3  ',     !    23: Profile pointer, partition 3      (n12p)
c              'JP4  ',     !    24: Profile pointer, partition 4      (n12p)

c     Mesh arrays

     &         'D    ', 'DR   ', 'F    ', 'F0   ', 'FPRO ', 'FTN  ',
     &         'ID   ', 'IE   ', 'IX   ', 'LD   ', 'P    ', 'S    ',
     &         'SLODI', 'T    ', 'TL   ', 'U    ', 'UL   ', 'VEL  ',
     &         'X    ', 'XL   ', 'ANG  ', 'ANGL ',


c              'D    ',     !    25: Material parameters               (n6)
c              'DR   ',     !    26: Residual/reactions                (1)

c              'F    ',     !    27: Nodal load/displacement, current  (n10)
c              'F0   ',     !    28: Nodal load/displace pattern, base (n13)
c              'FPRO ',     !    29: DOF proportional load numbers     (mpro)
c              'FTN  ',     !    30: Nodal load/displacement, current  (nt)

c              'ID   ',     !    31: Equation numbers/boundary conds   (n7)
c              'IE   ',     !    32: Element assembly information      (n5)
c              'IX   ',     !    33: Element connection data           (n9)

c              'LD   ',     !    34: Element local/global eq numbers   (n2)

c              'P    ',     !    35: Element vector                    (n3)

c              'S    ',     !    36: Element array                     (n4)
c              'SLODI',     !    37: Surface load data

c              'T    ',     !    38: Nodal temperatures                (n11)
c              'TL   ',     !    39: Nodal temperaturese, element      (n1)

c              'U    ',     !    40: Nodal solutions/increments        (n14)
c              'UL   ',     !    41: Nodal solutions/increments,element(nn)

c              'VEL'  ,     !    42: Nodal transient solution values   (nv)

c              'X    ',     !    43: Nodal coordinates (n8)
c              'XL   ',     !    44: Nodal coordinates, element (n0)

c              'ANG  ',     !    45: Nodal angles                      (n11b)
c              'ANGL ',     !    46: Nodal angles, element             (n11a)

c     Sparse solution data

     &         'PNTER', 'INVPT',

c              'PNTER',     !    47: Pointer array
c              'INVPT',     !    48: Inverse pointer array

c     History data

     &         'H    ', 'NH1  ', 'NH2  ', 'NH3  ',

c              'H    ',     !    49: Element history parameters
c              'NH1  ',     !    50: Element history data at t_n       (nh1)
c              'NH2  ',     !    51: Element history data at t_n+1     (nh2)
c              'NH3  ',     !    52: Element history data, time ind    (nh3)

c     Plot data

     &         'CT   ', 'FCIX ', 'FCZM ', 'FIDP ', 'NDER ', 'NDNP ',
     &         'NPAX ', 'NDNS ', 'OUTL ', 'SYMM ', 'TECON', 'TEFAC',
     &         'TENRM', 'VISN ',


c              'CT   ',     !    53: Plot deformed coordinate storage.
c              'FCIX ',     !    54: Face nodal connection list.
c              'FCZM ',     !    55: Face z-sort coordinates.
c              'FIDP ',     !    56: Faces for each node.
c              'NDER ',     !    57: Error estimator projections.      (nder)
c              'NDNP ',     !    58: Stress projection array.          (nph)
c              'NPAX ',     !    59: Principal axis projections.
c              'NDNS ',     !    60: Local stress projection array.
c              'OUTL ',     !    61: Outline construction.
c              'SYMM ',     !    62: Symmetric reflections table.
c              'TECON',     !    63: 3D edge definitions.
c              'TEFAC',     !    64: 3D edge factors.
c              'TENRM',     !    65: 3D edge normals
c              'VISN ',     !    66: Visible face list

c     Other data

     &         'SPTAN', 'AUR  ', 'BFGD ', 'BFGO ', 'BFGS ', 'BFGT ',
     &         'BFGV ', 'BFGW ', 'EIGE ', 'EVAL ', 'EVEC ', 'EXTND',
     &         'IPOS ', 'JPR  ', 'MO   ', 'MR   ', 'MT   ', 'MU1  ',
     &         'MU2  ', 'NDAM ', 'NMAS ', 'NSTI ', 'NREN ', 'OINMC',
     &         'OINMO', 'OINB ', 'OINC ', 'OINO ',


c              'SPTAN',     !    67: Sparse Tangent Array
c              'AUR  ',     !    68: Preconditioner Array

c              'BFGD ',     !    69: BFGS working vector               (mbd)
c              'BFGO ',     !    70: BFGS working vector               (mbo)
c              'BFGS ',     !    71: BFGS working vector               (mst)
c              'BFGT ',     !    72: BFGS working vector               (mbt)
c              'BFGV ',     !    73: BFGS vectors, U                   (mbv)
c              'BFGW ',     !    74: BFGS vectors, W                   (mbw)

c              'EIGE ',     !    75: Element eigenpairs                (ide)
c              'EVAL ',     !    76: Subspace eigenvalues              (md)
c              'EVEC ',     !    77: Subspace eigenvectors             (mv)
c              'EXTND',     !    78: External nodes

c              'IPOS ',     !    79: Tie node list

c              'JPR  ',     !    80: Preconditioner column pointers

c              'MO   ',     !    81: Rotational DOF option             (mropt)
c              'MR   ',     !    82: Rotational DOF lambda             (mrotas)
c              'MT   ',     !    83: Rotational DOF thickness          (mthick)
c              'MU1  ',     !    84:                                   (mu1)
c              'MU2  ',     !    85:                                   (mu2)

c              'NDAM ',     !    86: Nodal damping                     (ncc)
c              'NMAS ',     !    87: Nodal mass                        (ncm)
c              'NSTI ',     !    88: Nodal stiffness                   (nck)
c              'NREN ',     !    89: Nodal renumbering list            (nren)

c              'OINMC',     !    90: Consistent mass equation pointers (inmc)
c              'OINMO',     !    91: Consistent mass entries/equation  (inmo)

c              'OINB ',     !    92: Block information for out-of-core (inb)
c              'OINC ',     !    93: Sparse equation pointers          (inc)
c              'OINO ',     !    94: Sparse entries/equation           (inr)

c     Rigid Body data

     &         'RCG  ', 'REQRB', 'REVO ', 'RINER', 'RIRB ', 'RIXT ',
     &         'RJNT ', 'RJNX ', 'RJTU ', 'RLAMB', 'RLIST', 'RLOAD',
     &         'RMASS', 'RUROT', 'REXMS', 'REXIN',


c              'RCG  ',     !    95: Rigid body: center of mass        (nrcg)
c              'REQRB',     !    96: Rigid body: equation numbers      (nreqrb)
c              'REVO ',     !    97: Rigid body: revolutes             (nrev)
c              'RINER',     !    98: Rigid body: inertia dyadic        (niner)
c              'RIRB ',     !    99: Rigid body: rigid body nodes      (nrirb)
c              'RIXT ',     !   100: Rigid body: rigid node indicators (nrixt)
c              'RJNT ',     !   101: Rigid body: joint data            (nrjt)
c              'RJNX ',     !   102: Rigid body: joint coordinates     (nrjx)
c              'RJTU ',     !   103: Rigid body: joint solution values (nrju)
c              'RLAMB',     !   104: Rigid body: rotation lambda values(nrlamb)
c              'RLIST',     !   105: Rigid body:                       (nrlist)
c              'RLOAD',     !   106: Rigid body: loads                 (nrld)
c              'RMASS',     !   107: Rigid body: mass                  (nrmass)
c              'RUROT',     !   108: Rigid body: rotational solutions  (nrurot)

c              'REXMS',     !   109: Rigid body: Explicit mass coupling terms
c              'REXIN',     !   110: Rigid body: Explicit 6x6 inertias

c     Temporay arrays

     &         'TEMP1', 'TEMP2', 'TEMP3', 'TEMP4', 'TEMP5', 'TEMP6',
     &         'TEMP7', 'TEMP8', 'TEMP9', 'TEMP0',

c              'TEMP1',     !   111:  Temporary array
c              'TEMP2',     !   112:  Temporary array
c              'TEMP3',     !   113:  Temporary array
c              'TEMP4',     !   114:  Temporary array
c              'TEMP5',     !   115:  Temporary array
c              'TEMP6',     !   116:  Temporary array
c              'TEMP7',     !   117:  Temporary array
c              'TEMP8',     !   118:  Temporary array
c              'TEMP9',     !   119:  Temporary array
c              'TEMP0',     !   120:  Temporary array

c     Proportional loading table arrays (Type 2 prop. loads)

     &         'PROP0', 'PROP1', 'PROP2',

c              'PROP0',     !   121:  Prop. load offset table
c              'PROP1',     !   122:  Prop. load values table
c              'PROP2',     !   123:  Prop. load temporary use

c     Multiple support base excitations

     &         'PROBS', 'NUBAS', 'MASBS', 'PHIBS',

c              'PROBS',     !   124:  Base proportional factors
c              'NUBAS',     !   125:  Base pattern indicators
c              'MASBS',     !   126:  Static Mass projections
c              'PHIBS',     !   127:  Static modes

c     Follower nodal loads

     &         'APLOT', 'FOLLI', 'FOLLR',

c              'APLOT',     !   128:  Tag active plot elements
c              'FOLLI ',    !   129:  Follower force nodes
c              'FOLLR ',    !   130:  Follower force values

c     Contact arrays

     &         'C0   ', 'CM   ', 'ICS  ', 'HIC  ', 'CH   ',


c              'C0   ',     !   131:  Command control table            (ncc0)
c              'CM   ',     !   132:  Material table                   (ncm)
c              'ICS  ',     !   133:  List of nodes for geometry       (nics)
c              'HIC  ',     !   134:  History correspondence vector    (nhic)
c              'CH   ',     !   135:  History values            (ch1,ch2,ch3)

c     Contact

     &         'CTEM1', 'CTEM2', 'CTEM3', 'CTEM4', 'CTEM5','CTEM6',
     &         'CTEM7', 'CTEM8', 'CTEM9', 'CTE10', 'CTE11','CTE12',
     &         'CTE13', 'CTE14', 'CTE15',

c              'CTEM1',     !   136:  Contact temporary
c              'CTEM2',     !   137:  Contact temporary
c              'CTEM3',     !   138:  Contact temporary
c              'CTEM4',     !   139:  Contact temporary
c              'CTEM5',     !   140:  Contact temporary
c              'CTEM6',     !   141:  Contact temporary
c              'CTEM7',     !   142:  Contact temporary
c              'CTEM8',     !   143:  Contact temporary
c              'CTEM9',     !   144:  Contact temporary
c              'CTE10',     !   145:  Contact temporary
c              'CTE11',     !   146:  Contact temporary
c              'CTE12',     !   147:  Contact temporary
c              'CTE13',     !   148:  Contact temporary
c              'CTE14',     !   149:  Contact temporary
c              'CTE15',     !   150:  Contact temporary

c     User arrays

     &         'USER1', 'USER2', 'USER3', 'USER4', 'USER5', 'USER6',
     &         'USER7', 'USER8', 'USER9', 'USER0',

c              'USER1',     !   151:  Temporary array
c              'USER2',     !   152:  Temporary array
c              'USER3',     !   153:  Temporary array
c              'USER4',     !   154:  Temporary array
c              'USER5',     !   155:  Temporary array
c              'USER6',     !   156:  Temporary array
c              'USER7',     !   157:  Temporary array
c              'USER8',     !   158:  Temporary array
c              'USER9',     !   159:  Temporary array
c              'USER0',     !   160:  Temporary array

c     Blending arrays

     &         'BNODE', 'BSIDE', 'BTRAN', 'BLEND', 'BFACE', 'BNILR',

c              'BNODE',     !   161:  Super nodes for blending functions
c              'BSIDE',     !   162:  Sides for blending functions
c              'BTRAN',     !   163:  Transformations for blends
c              'BLEND',     !   164:  Blending function storage
c              'BFACE',     !   165:  Blending function storage
c              'BNILR',     !   166:  Blending layer storage

c     Rigid array

     &         'RLINK',

c              'RLINK',     !   167:  Rigid body link definitions.

c     Contact element connection array (total active)

     &         'IXC  ',

c              'IXC  ',     !   168:  Contact connection array

c     Control arrays for modal analyses

     &         'MCTRL', 'CCTRL', 'KCTRL',

c              'MCTRL',     !   169:  Mass      control array
c              'CCTRL',     !   170:  Damping   control array
c              'KCTRL',     !   171:  Stiffness control array

     &         'SVDA ', 'SVDV ', 'SVDW ', 'SVDX ',

c              'SVDA ',     !   172:  Singular valued decomp: A
c              'SVDV ',     !   173:  Singular valued decomp: V
c              'SVDW ',     !   174:  Singular valued decomp: W
c              'SVDX ',     !   175:  Singular valued decomp: X

c     Modal/Rigid node number array

     &         'IMODF', 'AFD  ', 'AFL  ', 'AFU  ', 'BFORC',

c              'IMODF',     |   176:  Modal equation numbers
c              'AFD  ',     |   177:  Modal stiffness diagonal
c              'AFL  ',     |   178:  Modal lower stiffness
c              'AFU  ',     |   179:  Modal upper stiffness
c              'BFORC',     |   180:  Modal force vector

c     Additional rigid body and modal arrays

     &         'RBEN','RBOU','UMOD',

c              'RBEN',      |   181:  Rigid body number associated with element
c              'RBOU',      |   182:  Rigid body boundary restraints
c              'UMOD',      |   183:  Modal displacement value

c     Modal data

     &         'CTEMP', 'KTEMP', 'MTEMP', 'FSMOD', 'YYMOD', 'WBASE',

c              'CTEMP',     !   184: Modal damping parameters
c              'KTEMP',     !   185: Modal stiffness parameters
c              'MTEMP',     !   186: Modal mass parameters
c              'FSMOD',     !   187: Modal forces                      (mfs)
c              'YYMOD',     !   188: Modal solution parameters         (my)
c              'WBASE',     !   189: Modal base solution parameters

c     Node type data

     &         'NDTYP',

c              'NDTYP'      !   190: Node type tags

c     Contact variables for surface descriptors

     &         'KNOTN', 'SURFP', 'INSEG', 'CNSEG', 'PNSEG', 'XISEG',

c              'KNOTN',     !   191: Node - surface listing
c              'SURFP',     !   192: Surface points
c              'INSEG',     !   193: In segment indicator
c              'CNSEG',     !   194: Contact flag
c              'PNSEG',     !   195: Point   flag
c              'XISEG',     !   196: Surface coordinates

c     Stress projection arrays

     &         'NS1  ', 'NS2  ',

c              'NS1  ',     !   197:                                   (ns1)
c              'NS2  ',     !   198:                                   (ns2)

c     Beam surface plot arrans

     &         'MXBEA', 'XBEAM', 'SBEAM', 'WBEAM',

c              'MXBEA',     !   199: Mesh for surface mesh of beams
c              'XBEAM',     !   200: Coordinates for surface mesh
c              'SBEAM',     !   201: Stress for surface mesh
c              'WBEAM',     !   202: Weights for surface mesh

c     Consistent damping arrays

     &         'OINDC', 'OINDO',

c              'OINDC',     !   203: Consistent damping equation pointers
c              'OINDO',     !   204: Consistent damping entries/equation

c     Slave boundary array

     &         'NSLAV',

c              'NSLAV',     !   205: Slave node numbers

c     Normal vector

     &         'NORMV','NSCR ','VTILD','DELTX',

c              'NORMV',     !   206: Normal vector (shell use)
c              'NSCR ',     !   207: TRI2D size projection data.
c              'VTILD',     !   208: Broyden vector 1
c              'DELTX',     !   209: Broyden vector 2

c     Interface storage and Lagrange multiplier

     &         'MATCH','LAGRE','LAGRN','ULAGR','HINTE','HINT1',
     &         'HINT2','HINT3',

c              'MATCH',     !   210: Lagrange multiplier solutions
c              'LAGRE',     !   211: Lagrange multiplier equations
c              'LAGRN',     !   212: Lagrange multiplier equations
c              'ULAGR',     !   213: Lagrange multiplier solutions
c              'HINTE',     !   214: Interface history variables
c              'HINT1',     !   215: Interface history variables
c              'HINT2',     !   216: Interface history variables
c              'HINT3',     !   217: Interface history variables

c     Zienkiewicz-Zhu Projector arrays

     &         'ZZIB ','ZZIP ',

c              'ZZIB ',     !   218: Zienkiewicz-Zhu boundary nodes
c              'ZZIP ',     !   219: Zienkiewicz-Zhu active nodes

c     Auto contact pointers

     &         'ACON2','ASLD2','ACIQ2','ACON3',

c     Short description of variables

c              'ACON2',     !   220: Autocon array: length = numnp
c              'ASLD2',     !   221: Autocon slideline: lg = 2*numnp+8
c              'ACIQ2',     !   222: Autocon IP   :     lg = ip(numnp)
c              'ACON3',     !   223: Autocon array :    lg = ip(numnp)*2

c     Sparse solver pointer

     &         'IAD  ',

c              'IAD  ',     !   224: Contact lagrange multiplier array

c     User solver pointers

     &         'USOL1','USOL2','USOL3','USOL4','USOL5','USOL6',
     &         'USOL7','USOL8','USOL9','USOL0',

c              'USOL1',     !   225: User solver array
c              'USOL2',     !   226: User solver array
c              'USOL3',     !   227: User solver array
c              'USOL4',     !   228: User solver array
c              'USOL5',     !   229: User solver array
c              'USOL6',     !   230: User solver array
c              'USOL7',     !   231: User solver array
c              'USOL8',     !   232: User solver array
c              'USOL9',     !   233: User solver array
c              'USOL0',     !   234: User solver array

c     Diagonal scaling array

     &         'DSCA1','DSCA2','DSCA3','DSCA4',

c              'DSCA1',     !   235: Reciprocal to diagonal sqare root
c              'DSCA2',     !   236: Reciprocal to diagonal sqare root
c              'DSCA3',     !   237: Reciprocal to diagonal sqare root
c              'DSCA4',     !   238: Reciprocal to diagonal sqare root

c     Interface type array

     &         'ITYPE','IEDOF',

c              'ITYPE',     !   239: Interface types
c              'IEDOF',     !   240: Interface types

c     Surface load real data

     &         'SLODR','EULER','LEULR',

c              'SLODR',     !   241: Surface load real data
c              'EULER',     !   242: Euler angles
c              'LEULR',     !   243: Euler angles, element

c     Parallel solver arrays

     &         'GN   ','EQ   ','DNNZ ','ONNZ ','GETP ',
     &         'GETV ','SENP ','SENV ',

c              'GN   ',     !   244: Local to Global node mapping
c              'EQ   ',     !   245: Local node to Global eq mapping
c              'DNNZ ',     !   246: Number non-zero diag-block entries
c              'ONNZ ',     !   247: No. non-zero off-diag-block entries
c              'GETP ',     !   248: Get data partition pointers
c              'GETV ',     !   249: Get data node values
c              'SENP ',     !   250: Send data partition pointers
c              'SENV ',     !   251: Send data node values

c     Mesh partioner (METIS/PARMETIS) arrays

     &         'XADJN ','NODG ','NODPT','VXDST',

c              'XADJN',     !   252: Pointer for nodal adjacency list
c              'NODG ',     !   253: Nodal adjacency list
c              'NODPT',     !   254: Nodal partition numbers
c              'VXDST',     !   255: Distribution array

c     Nodal follower forces

     &         'NFORC'/

c              'NFORC',     !   256: Nodal follower forces

c     Set memory pointers in blank common and list length

      if(num.eq.0) then

c       Zero pointer arrays

        do i = 1,num_nps
          np(i) = 0
        end do ! i

        do i = 1,num_ups
          up(i) = 0
        end do ! i

        do i = 1,400
          ilist(1,i) = 0
          ilist(2,i) = 0
        end do ! i
        llist  =  num_nps

      endif

c     Check user allocations then do allocation operation

      palloc = ualloc(num-llist,name,length,precis)

      end
