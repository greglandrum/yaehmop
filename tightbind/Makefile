#
#  This is the makefile for the tightbinding programs
# 
# created by greg Landrum August 1993
#


# to change where the program looks for the extended hueckel parameters
#   change the following:
# This is the one for isolobal
PARM_FILE_LOC = \"/home/glandrum/yaehmop/tightbind/eht_parms.dat\"

# to use floats instead of doubles, insert -DUSE_FLOATS here.
#  to perform real symmetry analysis (with principle axis location), use
#   -DREAL_SYMM_ANALYSIS  here.  To use this you must have the relevant
#     portions of the meschach libraries installed.

#stripped down with optimization and LAPACK
CFLAGS =  -g  -I/usr/local/include/ -DEHT_PARM_FILE=$(PARM_FILE_LOC) -DUNDERSCORE_FORTRAN -DANAL_ABOUT_PROTOTYPES

# these are the flags that get passed to the fortran compiler
#F77FLAGS = -g
F77FLAGS = -O

#A linux box w/o complib (using GNU compilers)
#LIBS = -lm -lf2c 
LIBS = -lm 
LOCALLIBS = /usr/lib/liblapack.a /usr/lib/libblas.a -lg2c

COBJS = main.o genutil.o fileio.o globals.o memory.o distance_mat.o \
 R_overlap_mat.o mov.o R_hamil.o K_hamil.o kpoints.o K_overlap_mat.o \
 mulliken.o charge_mat.o zetas.o electrostat.o new3_fileio.o walsh.o \
 transforms.o symmetry.o princ_axes.o avg_props.o DOS_stuff.o COOP_stuff.o \
 Zmat.o bands.o FMO_stuff.o xtal_coords.o matrices.o chg_it.o \
 mod_mulliken.o postprocess.o muller.o geom_frags.o solid_symmetry.o \
 recip_space.o netCDF_support.o 


#F2COBJS = lovlap.f2c.o abfns.f2c.o cboris.f2c.o diag.f2c.o
FOBJS = lovlap.o abfns.o cboris.o diag.o

all: bind


bind: $(COBJS) $(FOBJS) $(LAPACKOBJS)
	cc -o bind $(CFLAGS) $(COBJS) $(FOBJS)  $(LAPACKOBJS) $(LIBS) $(LOCALLIBS)

bind.64: $(COBJS) $(FOBJS) $(LAPACKOBJS)
	cc -o bind.64 $(CFLAGS) $(COBJS) $(FOBJS) $(LAPACKOBJS) $(LIBS)

install: bind 
	cp bind ../bin
	cd utils;make install

install64: bind.64
	cp bind.64 ../bin

clean:
	rm -f $(COBJS) $(FOBJS)


veryclean:
	rm -f bind bind.64
	rm -f $(COBJS) $(FOBJS) $(LAPACKOBJS)
	cd utils;make clean

.f.o: 
	f77 $(F77FLAGS) -c $*.f

.c.o: 
	cc $(CFLAGS) -c $*.c


fit_progs: $(FIT_PROGS)

utils/fit_dos: utils/fit_dos.c
	cd utils; make fit_dos

utils/fit_coop: utils/fit_coop.c
	cd utils; make fit_coop

utils/fit_walsh: utils/fit_walsh.c
	cd utils; make fit_walsh

utils/sub_dos: utils/sub_dos.c
	cd utils; make sub_dos
