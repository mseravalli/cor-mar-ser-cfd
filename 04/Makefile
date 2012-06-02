IDIR =./
CC=mpicc
CFLAGS = -Wall -pedantic -g

ODIR=obj
LDIR =../lib

LIBS=-lm 

_DEPS = helper.h\
		init.h\
		boundary_val.h\
		uvp.h\
		sor.h\
		visual.h\
        parallel.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	visual.o\
        parallel.o\
      	main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

sim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
