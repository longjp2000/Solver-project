
 NPROC     = 2                                      # Set number of processors

 FEAPRUN   = $(FEAPHOME8_1)/parfeap/feap               # Set path for executable
 ARPACKLIB = $(FEAPHOME8_1)/packages/arpack/archive/arpacklib.a # arpack library
#LAPACKLIB = $(FEAPHOME8_1)/packages/lapack/lapacklib.a # Local lapack library
#BLASLIB   = $(FEAPHOME8_1)/packages/blas/blaslib.a     # Local blas library

include $(FEAPHOME8_1)/makefile.in

include ${PETSC_DIR}/bmake/common/base

PETSC_INCLUDE := $(PETSC_INCLUDE) -I$(FINCLUDE)    # N.B. := avoids recursion

OBJECTS = $(FEAPHOME8_1)/main/feap81.o \
	adomnam.o filargs.o fparsop.o parbmat.o pcompress.o \
	pdomain.o pfeapmi.o pfeapsr.o pform.o parstop.o parmacio.o \
	pmacr7.o unix/p_metis.o pmodify.o ppmodin.o pparlo.o \
	pndata.o prwext.o psetb.o psproja.o psprojb.o pstart.o \
	psubsp.o scalev.o smodify.o uasble.o uasblem.o upremas.o \
	usolve.o upc.o pextnd.o tinput.o pprint.o

feap: $(OBJECTS) $(ARFEAP)
	ranlib $(ARFEAP)
	-${FLINKER} ${FFLAGS} -o feap $(OBJECTS) \
	$(ARPACKLIB) $(ARFEAP) \
	${PETSC_FORTRAN_LIB} ${PETSC_LIB} $(LDOPTIONS)

checkout:
	co -q RCS/*.h,v
	co -q RCS/*.f,v
	co -q RCS/*.F,v

rcs:
	mkdir RCS
	ci -t-"" *.h
	ci -t-"" *.f
	ci -t-"" *.F
	ci -t-"" makefile

clean:
	rcsclean -q
	rm -f *.o

feaprun:
	-@${MPIRUN} -np $(NPROC) $(FEAPRUN) -ksp_type cg -ksp_monitor -pc_type jacobi -log_summary -on_error_attach_debugger -ksp_view -options_left

feaprun-gmres:
	-@${MPIRUN} -np $(NPROC) $(FEAPRUN) -ksp_type gmres -ksp_monitor -pc_type bjacobi -log_summary -on_error_attach_debugger

feaprun-mg:
	-@${MPIRUN} -np $(NPROC) $(FEAPRUN) -ksp_type cg -ksp_monitor -log_summary -out_verbose 2 -pc_type prometheus -pc_mg_type multiplicative -pc_mg_cycles 1-prometheus_preduce_base 500 -prometheus_top_grid_limit 2500 -aggmg_avoid_resmooth -aggmg_smooths 1  -prometheus_repartition -options_left

.f.o:
	$(FF) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o
#
#.c.o:
#	$(CC) -c $(CCOPTFLAG)  $*.c -o $*.o

#FEAPRUN   = valgrind -v --leak-check=full --show-reachable=yes $(FEAPHOME8_1)/parfeap/feap               # set path for executable

