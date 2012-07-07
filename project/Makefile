IDIR =./
CC=gcc
CFLAGS = -Wall -pedantic -Werror -g -fopenmp
#CFLAGS = -funroll-loops -O3 -fopenmp

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = constants.h\
        helper.h\
		init.h\
		boundary_val.h\
		uvp.h\
		sor.h\
		visual.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	visual.o\
      	main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

sim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
