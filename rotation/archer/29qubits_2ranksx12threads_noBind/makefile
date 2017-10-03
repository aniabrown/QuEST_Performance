#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- COMMON CONFIG
#

# COMPILER options: GNU, INTEL
COMPILER = INTEL
CRAY_SYSTEM = 1
EXE = demo
MY_FILE_NAME = timingRotation
USE_MPI = 1
USE_OPENMP = 1
QUEST_DIR = ../../../QuEST_v0.9.0

#
# --- compiler
#

ifeq ($(COMPILER), GNU)
	CFLAGS_OMP=-fopenmp
else 
	CFLAGS_OMP=-openmp
endif

ifneq ($(USE_OPENMP), 1)
	# disable OPENMP
	CFLAGS_OMP=
endif

ifneq ($(USE_MPI), 1)
	ifeq ($(COMPILER), GNU)
		# GCC compilers
		CCOMPILER         = gcc
		CFLAGS     = -O2 -std=c99 -mavx -Wall
	else ifeq ($(COMPILER), INTEL)
		# Intel compilers
		CCOMPILER         = icc
		CFLAGS     = -O2 -fprotect-parens -std=c99 -Wall -xAVX -axCORE-AVX2 -restrict
	else 
		$(error " *** error: invalid compiler")
	endif
else 
	ifeq ($(COMPILER), GNU)
		# GCC compilers
		CCOMPILER         = mpicc
		CFLAGS     = -O2 -std=c99 -mavx -Wall
	else ifeq ($(COMPILER), INTEL)
		# Mvapich2
		CCOMPILER         = mpicc
		CFLAGS     = -O2 -fprotect-parens -std=c99
	else 
		$(error " *** error: invalid compiler")
	endif
endif

ifeq ($(CRAY_SYSTEM), 1)
	# cray has a wrapper around c compilers
	CCOMPILER = cc
endif

#
# --- libraries
#
LIBS = -lm


#
# --- targets
#
OBJ = $(MY_FILE_NAME).o qubits.o
ifneq ($(USE_MPI), 0)
	OBJ += qubits_env_mpi.o
else
	OBJ += qubits_env_local.o
endif

#
# --- rules
#
%.o: %.c
	$(CCOMPILER) $(CFLAGS) $(CFLAGS_OMP) -c $<

%.o: $(QUEST_DIR)/%.c
	$(CCOMPILER) $(CFLAGS) $(CFLAGS_OMP) -c $<


#
# --- build
#
default:	$(EXE)

$(EXE):		$(OBJ)
		$(CCOMPILER) $(CFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) $(LIBS)

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o $(EXE)
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
