IDIR =./
CC=gcc
CFLAGS = -Wall -pedantic -Werror -g

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = helper.h\
		init.h\
		boundary_val.h\
		uvp.h\
		sor.h\
		visual.h\
		disc.h\
		matrix_op.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
      	sor.o\
      	main.o\
      	visual.o\
        disc.o\
        matrix_op.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

sim: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
