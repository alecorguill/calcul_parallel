CC = gcc
MPI_CC = mpicc
CFLAGS = -W -Wall -g
MPI_FLAGS = -lm -std=c99
LDFLAGS = 
SRC = $(wildcard *.c)
OBJS = $(SRC:.c=.o)
EXEC = fem.exe
TRASH = *~
NB_PROC = 4	


all : $(EXEC) 

create_partition : Preprocess.exe
	./Preprocess.exe $(NB_PROC) mesh_for_progc.data

Preprocess.exe : Preprocess.c
	 $(CC) $(CFLAGS) -o $@ -c $<

fem.exe : FemPar.c
	$(MPI_CC) $(MPI_FLAGS) -o fem.exe FemPar.c
clean :
	rm -f *.exe $(TRASH)

