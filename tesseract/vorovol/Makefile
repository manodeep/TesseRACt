#--------------------------------------- Select some defaults
CC       = gcc
OPTIMIZE = -g -O0 -march=native -mavx -Wall -Wextra -std=c99

#--------------------------------------- Select target computer
SYSTYPE = "Accre"

#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"Accre")
CC        = gcc
OPTIMIZE  = -O3 -march=native -mavx -Wall -Wextra -std=c99
MLIBS     = -lm
QLIB      = -L/home/langmm/code/python/mmlpy/trunk/tesseract/qhull2002.1/src -lqhull
QINC      = -I/home/langmm/code/python/mmlpy/trunk/tesseract/qhull2002.1/src
endif


OPTIONS =  $(OPTIMIZE) $(OPT)

OBJS   = main.o io.o vozutil.o allvars.o

EXEC   = vorovol

INCL   = allvars.h proto.h Makefile

CFLAGS = $(OPTIONS) $(QINC) -fopenmp

LIBS   = $(MLIBS) $(QLIB)


$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
