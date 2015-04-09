# Computer name

MACHINE = $(shell hostname -s)

# Code name

MAKEFILE = Makefile

# Object files

OBJS = types.o vectors.o kernels.o grid.o mesh.o vars.o

# Conditionals

FF = gfortran
FCOMPILEFLAGS = -Ofast -fopenmp #-floop-parallelize-all -ftree-parallelize-loops=4
LD = $(FF)
LINKFLAGS = $(FCOMPILEFLAGS)
LIBS = tecio64.a -lstdc++

# Flags

FCOMPILE = $(FF) $(FCOMPILEFLAGS)
CCOMPILE = $(CC) $(CCOMPILEFLAGS)
LINK = $(LD) $(LINKFLAGS)

# Build targets

sph: $(OBJS) main.o $(MAKEFILE)
	$(LINK) -o sph $(OBJS) main.o $(LIBS)

plot: $(OBJS) plot.o $(MAKEFILE)
	$(LINK) -o plot $(OBJS) plot.o $(LIBS)

grid2d: $(OBJS) grid2d.o $(MAKEFILE)
	$(LINK) -o grid2d $(OBJS) grid2d.o $(LIBS)

couette2d: $(OBJS) couette2d.o $(MAKEFILE)
	$(LINK) -o couette $(OBJS) couette2d.o $(LIBS)
	
taylorgreen2d: $(OBJS) taylorgreen2d.o $(MAKEFILE)
	$(LINK) -o taylorgreen2d $(OBJS) taylorgreen2d.o $(LIBS)
	
taylorgreen3d: $(OBJS) taylorgreen3d.o $(MAKEFILE)
	$(LINK) -o taylorgreen3d $(OBJS) taylorgreen3d.o $(LIBS)

shear2d: $(OBJS) shear2d.o $(MAKEFILE)
	$(LINK) -o shear $(OBJS) shear2d.o $(LIBS)

dam2d: $(OBJS) dam2d.o $(MAKEFILE)
	$(LINK) -o dam $(OBJS) dam2d.o $(LIBS)

spindown2d: $(OBJS) spindown2d.o $(MAKEFILE)
	$(LINK) -o spindown2d $(OBJS) spindown2d.o $(LIBS)

dam3d: $(OBJS) dam3d.o $(MAKEFILE)
	$(LINK) -o dam $(OBJS) dam3d.o $(LIBS)

weir3d: $(OBJS) weir3d.o $(MAKEFILE)
	$(LINK) -o weir3d $(OBJS) weir3d.o $(LIBS)

clean:
	-rm *.o
	-rm *.mod
	-rm *.f90~
	-rm *.sh~
	-rm Makefile~
	-rm -rf core*
	-rm sph
	-rm plot
	-rm grid2d
	-rm couette
	-rm spindown2d
	-rm taylorgreen2d
	-rm taylorgreen3d
	-rm shear
	-rm dam
	-rm weir3d
	-rm sph.tar

tar:
	tar -cvf sph.tar *

# Dependencies

$(OBJS): $(MAKEFILE)

types.o: types.f90

vectors.o: types.f90 vectors.f90

kernels.o: types.f90 kernels.f90

grid.o: types.f90 vectors.f90 grid.f90

vars.o: types.f90 grid.f90 vars.f90

mesh.o: types.f90 vectors.f90 kernels.f90 grid.f90 mesh.f90

main.o: types.f90 kernels.f90 vars.f90 grid.f90 main.f90

main.o: types.f90 kernels.f90 vars.f90 grid.f90 mesh.f90 plot.f90

grid2d.o: types.f90 kernels.f90 vars.f90 grid.f90 grid2d.f90

couette2d.o: types.f90 kernels.f90 vars.f90 grid.f90 couette2d.f90

taylorgreen2d.o: types.f90 kernels.f90 vars.f90 grid.f90 taylorgreen2d.f90

taylorgreen3d.o: types.f90 kernels.f90 vars.f90 grid.f90 taylorgreen3d.f90

shear2d.o: types.f90 kernels.f90 vars.f90 grid.f90 shear2d.f90

dam2d.o: types.f90 kernels.f90 vars.f90 grid.f90 dam2d.f90

dam3d.o: types.f90 kernels.f90 vars.f90 grid.f90 dam3d.f90

weir3d.o: types.f90 kernels.f90 vars.f90 grid.f90 weir3d.f90

# Default suffix rules

%.o : %.f90
	$(FCOMPILE) -c $<

%.o : %.c
	$(CCOMPILE) -c $<
