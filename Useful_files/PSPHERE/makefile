all: psphere

# Which compiler
CC = g++ 

# Where are include files kept
INCLUDE = . 

# make sure to install libglib2.0-dev
# Options for development
CFLAGS = -g -O2 -pedantic -Wall -ansi `pkg-config --cflags glib-2.0`
#-I /usr/lib/glib/include  
LFLAGS = -g -lgsl -lgslcblas -lm -lpng `pkg-config --cflags --libs glib-2.0`

# Options for release
#CFLAGS = -O6 -pedantic -mcpu=pentium -Wall -ansi
#LFLAGS = -O6 -g -lrfftw -lfftw -lgsl -lgslcblas -lm

#object files
OBS =  nrutil.o matrix.o png_image.o

psphere: psphere.o $(OBS)
	$(CC) -o psphere psphere.o $(OBS) $(LFLAGS)

#gcc -MM *.c output
psphere.o: psphere.cc nrutil.h matrix.h
	$(CC) -c psphere.cc $(CFLAGS)
nrutil.o: nrutil.cc nrutil.h
	$(CC) -c nrutil.cc $(CFLAGS)
matrix.o: matrix.cc matrix.h
	$(CC) -c matrix.cc $(CFLAGS)
png_image.o: png_image.cc png_image.h
	$(CC) -c png_image.cc $(CFLAGS)


