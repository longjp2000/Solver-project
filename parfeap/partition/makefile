#  Compile the parallel parmetis partioner with 'make partition'
#  Run the partioner with 'make run'
#  See further notes below!

CFLAGS=${PARMETIS_INCLUDE}

include ${PETSC_DIR}/bmake/common/base

partition: partition.o sparse-graph.o
	-${CLINKER} -o partition partition.o sparse-graph.o \
		${PARMETIS_LIB} ${MPI_LIB}

NPROC=2
NPART=2
INAME=Iblock.rev

#
# -np number controls the number of processors you will be using to partition
# your problem.   The number before the inputfile name controls the number of
# processors FOR WHICH you are partitioning the mesh
#
# The output will be graph.block (in the case below, strips off the leading I)
# Note that Iblock (in the case below) must be a flat file
#
# Next, start feap in the parfeap directory and then produce the parallel input
# files as:
# 
# Macro> graph,file Macro> outd
# 
# Note 'file' is the word "file" not "graph.block".
#

partrun:
	-@${MPIRUN} -np $(NPROC) \
	$(FEAPHOME)/parfeap/partition/partition $(NPART) $(INAME)

feaprun:
	-@${MPIRUN} -np 1 $(FEAPHOME)/parfeap/feap

run:
	-@${MPIRUN} -np $(NPROC) \
	$(FEAPHOME)/parfeap/partition/partition $(NPART) $(INAME)
	rm feapname
	echo $(INAME) > feap.in1
	cat feap.in1  $(FEAPHOME)/parfeap/partition/feap.in.base > feap.in
	cat feap.in | ${MPIRUN} -np 1 $(FEAPHOME)/parfeap/feap
	rm feap.in feap.in1

rcs:
	mkdir RCS
	ci -t-"" *.h
	ci -t-"" *.c
	ci -t-"" makefile
	ci -t-"" README

clean:
	rcsclean -q
	rm -f *.oq

checkout:
	co -q RCS/*.h,v
	co -q RCS/*.c,v
	co -q RCS/README,v

