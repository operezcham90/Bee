CC      = g++
MPI		= mpic++
INC     = `pkg-config --cflags opencv`
LIB     = -L/usr/lib
OPT     = `pkg-config --libs opencv`
BIN     = abejas sobel avi abejas_video abejas_foto abejas_foto_mpi abejas_debug fractal homogeneidad gabor rugosidad
all:    $(BIN)

homogeneidad:homogeneidad.cxx algos.cxx Makefile
			$(CC) homogeneidad.cxx algos.cxx $(INC) $(LIB) $(OPT) -o $@
abejas_foto:abejas_foto.cxx algos.cxx Makefile
			$(CC) abejas_foto.cxx algos.cxx $(INC) $(LIB) $(OPT) -o $@
abejas_foto_mpi:abejas_foto_mpi.cxx algos.cxx paralelizar.cxx Makefile
			$(MPI) abejas_foto_mpi.cxx algos.cxx paralelizar.cxx $(INC) $(LIB) $(OPT) -o $@
abejas_video:abejas_dyn.cxx algos.cxx Makefile
			$(CC) abejas_dyn.cxx algos.cxx $(INC) $(LIB) $(OPT) -o $@
abejas_video_mpi: abejas_dyn.cxx algos.cxx paralelizar.cxx Makefile
			$(MPI) abejas_dyn.cxx algos.cxx paralelizar.cxx $(INC) $(LIB) $(OPT) -o $@
abejas: abejas_NC.cxx algos.cxx vrml.cxx Makefile
			$(CC) abejas_NC.cxx algos.cxx vrml.cxx $(INC) $(LIB) $(OPT) -o $@
avi: avi.c Makefile
			$(CC) avi.c $(INC) $(LIB) $(OPT) -o $@
clean: 
	rm -f *˜ *.o $(BIN) core

