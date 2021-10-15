all:sq

CXXPPFLAGS = 

SOURCECPP = ./src/main.cc ./src/lsdft_initObjs.cc  ./src/lsdft_matvec.cc  ./src/lsdft_scf.cc ./src/lsdft_readfiles.cc ./src/lsdft_tools.cc ./src/lsdft_nl.cc

SOURCEH = ./src/lsdft.h ./src/ilsdft.h 

OBJSC = ./src/main.o ./src/lsdft_initObjs.o  ./src/lsdft_matvec.o  ./src/lsdft_scf.o  ./src/lsdft_readfiles.o ./src/lsdft_tools.o ./src/lsdft_nl.o

LIBBASE = ./bin/sq

CLEANFILES = ./bin/sq

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

sq: ${OBJSC} chkopts
	${CLINKER} -Wall -o ${LIBBASE} ${OBJSC} ${PETSC_LIB}
	${RM} $(SOURCECPP:%.cc=%.o)
