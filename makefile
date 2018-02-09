# Makefile to build lbstime.a library and driver routines 
#
# Andrew J. Pounds, Ph.D.
# Departments of Chemistry and Computer Science
# Mercer University
# Spring 2007
#

# Use gfortran to avoid g77/gcc incompatibilities on Fedora
# g77 works fine on systems where gcc/g77 were linked against the same library

# with this makefile if you type just "make" it will compile the code with the default
# dimension of 50.  If you type "make size=100" it will compile with dimension 100 or whatever
# other size you want to use

F77 = gfortran    
F90 = gfortran
CC  = gcc 
TIMINGLIB_FLAGS = -O3
OPT_FLAGS = #-O0  -fcombine-stack-adjustments -fipa-reference -fshrink-wrap -fssa-backprop -ftree-coalesce-vars -ftree-copy-prop#-O0 -DOPT1 -DOPT2 -DOPT3 -DOPT4 -DOPT5 -DOPT6 -fcompare-elim -fipa-pure-const -fsplit-wide-types -ftree-ccp -ftree-ch#-O0  -fcombine-stack-adjustments -fipa-reference -fshrink-wrap -fssa-backprop -ftree-coalesce-vars -ftree-copy-prop #-DOPT1 -DOPT2 -DOPT3 -DOPT4 -DOPT5 -DOPT6 -fcompare-elim -fipa-pure-const -fsplit-wide-types -ftree-ccp -ftree-ch #-ftree-dce -ftree-dominator-opts -ftree-dse -ftree-forwprop -ftree-fre -ftree-phiprop#-fssa-phiopt -ftree-bit-ccp -ftree-ccp -ftree-ch -ftree-coalesce-vars -ftree-copy-prop#-fmove-loop-invariants -freorder-blocks -fshrink-wrap -fsplit-wide-types -fssa-backprop #-fipa-pure-const -fipa-profile -fipa-reference 
 #-O1 -fstrict-aliasing   # -00 -fmove-loop-invariants -finline-functions-called-once
DEF_FLAGS = -DDO_TIMING #-DOPT #-DDO_TIMING 

TIMINGLIBS =  -L./ -llbstime 
CLIBS = -lm

OBJS = cputime.o walltime.o  

ifeq ($(origin size), undefined)
	PROFLAGS = -D CPPFLAGS=200
else
	PROFLAGS = -D CPPFLAGS=$(size)
endif

all: lib dusty77 dusty90 dusty 

cputime.o : cputime.cc   
	$(CC) $(TIMINGLIB_FLAGS) -c cputime.cc  

walltime.o : walltime.cc   
	$(CC) $(TIMINGLIB_FLAGS) -c walltime.cc  

dusty77.o : dusty77.f   
	$(F77) -cpp $(PROFLAGS) $(DEF_FLAGS) $(OPT_FLAGS) -c dusty77.f   

# Don't forget the -lstdc++
dusty77 : dusty77.o  $(OBJS) 
	$(F77) -o dusty77 dusty77.o  $(TIMINGLIBS) -lstdc++  

dusty90.o : dusty90.f90  
	$(F90) -cpp $(PROFLAGS) $(DEF_FLAGS) $(OPT_FLAGS) -c dusty90.f90   

dusty90 : dusty90.o  $(OBJS) 
	$(F90) -o dusty90 dusty90.o  $(TIMINGLIBS) -lstdc++  

dusty.o : dusty.c   
	$(CC) $(PROFLAGS) $(DEF_FLAGS) $(OPT_FLAGS) -c dusty.c   

# Don't forget the -lstdc++
dusty : dusty.o  $(OBJS) 
	$(CC) -o dusty dusty.o  $(TIMINGLIBS) $(CLIBS) -lstdc++   

# Default Targets for Cleaning up the Environment
clean :
	rm *.o
	rm *.a

pristine :
	rm *.o
	rm *.a
	touch *.c *.f *.f90 

ctags :
	ctags  *.c *.f *.f90

# Target for making the library

lib: $(OBJS) 
	ar -rc liblbstime.a $(OBJS) 
	ranlib liblbstime.a

run :
	dusty77 | tee testdata.dat
	dusty77 | tee -a testdata.dat
	dusty77 | tee -a testdata.dat
	dusty77 | tee -a testdata.dat
	dusty77 | tee -a testdata.dat

average :
	echo "mydata = read.table(\"testdata.dat\", header=FALSE)" >mean.R
	echo "mean(mydata[[\"V2\"]])" >>mean.R
	R -q --vanilla -f mean.R
	rm mean.R

