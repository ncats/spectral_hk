
CC = gcc
CFLAGS= -Wall -g

TARGETS = libspectral.a spectral_hk
OBJS = b32.o sha1.o jacobi.o spectral.o

LIBS = -lm

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
