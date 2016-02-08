## Please consider using either Makefile.gsl or Makefile.mkl
## The bundled eigensolver is really for demonstration only!
SUFFIX =
CC = clang

# uncomment to compile debug
DEBUG=-g #-DSPECTRAL_DEBUG
#DEBUG=-O3

######################################################################
## shouldn't have to edit below
######################################################################
TARGETS = libspectral.a spectral_hk$(SUFFIX) pi$(SUFFIX)
OBJS = b32.o sha1.o jacobi.o spectral.o periodic.o inchi.o features.o ring.o
CFLAGS= -Wall $(DEBUG)
LIBS = -lm 

.c.o: $(OBJS)
	$(CC) $(CFLAGS) -c $<

all: $(TARGETS)

libspectral.a: $(OBJS)
	$(AR) -r $@ $(OBJS) 

spectral_hk$(SUFFIX): libspectral.a spectral_hk.c
	$(CC) $(CFLAGS) -o $@  spectral_hk.c libspectral.a $(LIBS)

pi$(SUFFIX): libspectral.a pi.c
	$(CC) $(CFLAGS) -o $@ pi.c libspectral.a $(LIBS)

test: spectral_hk$(SUFFIX)
	./spectral_hk$(SUFFIX) examples.txt | sort

clean:
	$(RM) $(OBJS) $(TARGETS)
