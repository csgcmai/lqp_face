#
# Choose your compiler
#
# 
CC = g++ #gcc #-4.1
CPP = g++ #-4.1
CFLAGS =?  O3 -ffast-math -fomit-frame-pointer -fPIC -Werror
#
# MPI KMEANS FLAGS
#
VERBOSEFLAG=-DKMEANS_VERBOSE=1 #-DINPUT_TYPE=1 # 0: silent, 1:iteration counter, 2:everything
#
# MATLAB
#
MATLABDIR=~/matlab/
MATLAB_INCLUDE=-I$(MATLABDIR)/extern/include

#
# BOOST LIBS (for standalone only)
#
BOOST_LIB= -L /usr/lib64/ -lboost_filesystem -lboost_program_options

BOOST_INCLUDE= #-I/kyb/agbs/pgehler/include/boost-1_36/

#
# PYTHON
#
PYTHON_INCLUDE=-I/usr/include/python2.5
PYTHON_LIB=-lpython2.5
NUMPY_INCLUDE=-I/usr/lib/python2.5/site-packages/numpy/core/include

#
# ARCHITECURE
# 

# 32 bit
#SUFFIX=mexglx
#MATLAB_LIB=-L$(MATLABDIR)/bin/glnx86 -lmex

# 64 bit
SUFFIX=mexa64
MATLAB_LIB=-L$(MATLABDIR)/bin/glnxa64 -lmex

#LIBS= /usr/lib/gcc/x86_64-redhat-linux/4.1.2/libstdc++.a /usr/lib64/libm.a -L /usr/lib64/
LIBS =? -L /usr/lib64/

all: standalone #matlab libmpikmeans python
matlab: 	mpi_kmeans_mex.$(SUFFIX) mpi_assign_mex.$(SUFFIX)
standalone:	mpi_kmeans_main mpi_assign_main
python:	cython_wrapper

mpi_kmeans.o:	mpi_kmeans.cxx mpi_kmeans.h
	$(CC) $(CFLAGS) $(VERBOSEFLAG) $(PRECISION) -c -o $@ mpi_kmeans.cxx

libmpikmeans:	mpi_kmeans.o
	ar rc libmpikmeans.a mpi_kmeans.o
	ranlib libmpikmeans.a
	$(CC) -shared -Wl,-soname=libmpikmeans.so -fPIC $(CFLAGS) -o libmpikmeans.so $(VERBOSEFLAGS) $(PRECISION) mpi_kmeans.cxx

mpi_kmeans_main.o:	mpi_kmeans_main.cxx 
	$(CC) $(CFLAGS) $(BOOST_INCLUDE) -c -o mpi_kmeans_main.o mpi_kmeans_main.cxx 

mpi_assign_main.o:	mpi_assign_main.cxx
	$(CC) $(CFLAGS) $(BOOST_INCLUDE) -c -o mpi_assign_main.o mpi_assign_main.cxx 

mpi_kmeans_main:	mpi_kmeans_main.o mpi_kmeans.o
	$(CC) $^  $(CFLAGS)  -o mpi_kmeans \
		$(LIBS) $(BOOST_LIB) 

mpi_assign_main:	 mpi_assign_main.o mpi_kmeans.o
	$(CC) $^  $(CFLAGS)  -o mpi_assign \
		$(LIBS) $(BOOST_LIB) 

%_mex.o:	%_mex.cxx
	$(CC) $(CFLAGS) $(MATLAB_INCLUDE) $(VERBOSEFLAG) $(PRECISION) -c $^ -o $@

mpi_kmeans_mex.$(SUFFIX):	libmpikmeans mpi_kmeans_mex.o 
	$(CC) mpi_kmeans_mex.o -shared -o mpi_kmeans_mex.$(SUFFIX) libmpikmeans.a $(MATLAB_LIB) 

mpi_assign_mex.$(SUFFIX):	libmpikmeans mpi_assign_mex.o
	$(CC) mpi_assign_mex.o -shared -o mpi_assign_mex.$(SUFFIX) libmpikmeans.a $(MATLAB_LIB)

cython_wrapper:	py_kmeans.c mpi_kmeans.o
	$(CPP) $(CFLAGS) $(PYTHON_INCLUDE) $(NUMPY_INCLUDE) -c -o py_kmeans.o py_kmeans.c
	$(CPP) $(CFLAGS) $(PYTHON_LIB) -lm -pthread -shared py_kmeans.o mpi_kmeans.o  -o py_kmeans.so 

test:	
	matlab -nojvm -r "test_code;exit"

clean:	
	rm -f *.o
	rm -f *.mexglx
	rm -f *.mexa64
	rm -f libmpikmeans.so
	rm -f libmpikmeans.a
	rm -f mpi_assign mpi_kmeans

