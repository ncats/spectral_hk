##

CC = gcc
# uncomment to compile debug
#DEBUG=-DSPECTRAL_DEBUG

# uncomment to compile with gsl
GSLFLAGS=-DHAVE_GSL -I/opt/local/include

# uncomment to compile with gsl library; for macports, use -lcblas instead
GSLLIBS=-L/opt/local/lib -lgsl -lcblas ## macports
#GSLLIBS=-lgsl -lgslcblas ## linux

######################################################################
## shouldn't have to edit below
######################################################################
TARGETS = libspectral.a spectral_hk
OBJS = b32.o sha1.o jacobi.o spectral.o interval.o
CFLAGS= -Wall -O3 $(GSLFLAGS) $(DEBUG)
LIBS = -lm $(GSLLIBS)

.c.o: $(OBJS)
	$(CC) $(CFLAGS) -c $<

all: $(TARGETS)

libspectral.a: $(OBJS)
	$(AR) -r $@ $(OBJS) 

spectral_hk: libspectral.a spectral_hk.c
	$(CC) $(CFLAGS) -o $@  spectral_hk.c libspectral.a $(LIBS)

test: spectral_hk
	./spectral_hk examples.txt

clean:
	$(RM) $(OBJS) $(TARGETS)
