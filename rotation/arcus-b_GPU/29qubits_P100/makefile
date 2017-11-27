#======================================================================#
#                                                                      #
#      Makefile -- build the qubit function library                    #
#                                                                      #
#======================================================================#

#
# --- COMMON CONFIG
#

CC = nvcc
EXE = demo
MY_FILE_NAME = timingRotation
QUEST_DIR = ../../../QuEST_GPU_v0.5.0

#
# --- compiler
#

CFLAGS	   = -dc
CLFLAGS	   = -O2 -arch=compute_60 -code=sm_60 -lineinfo
	
#
# --- targets
#
OBJ = $(MY_FILE_NAME).o qubits.o qubits_env_localGPU.o

#
# --- rules
#
%.o: %.cu
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<
%.o: %.cpp
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<

%.o: $(QUEST_DIR)/%.cu
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<
%.o: $(QUEST_DIR)/%.cpp
	$(CC) $(CFLAGS) $(CLFLAGS) $(CFLAGS_OMP) $<


#
# --- build
#
default:	demo

demo:		$(OBJ)
		$(CC) $(CLFLAGS) $(CFLAGS_OMP) -o $(EXE) $(OBJ) 

.PHONY:		clean veryclean
clean:
		/bin/rm -f *.o demo
veryclean:	clean
		/bin/rm -f *.h~ *.c~ makefile~
