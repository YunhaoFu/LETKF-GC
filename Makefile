FC=mpifort
FCFLAGS=-I$(NETCDF)/include
LDFLAGS=-L$(NETCDF)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lblas -llapack #-L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -ldl # module load mkl/v4.5

objs = precision.o ini_fnl.o io_control.o prepost.o kdtree.o localization.o letkf.o letkf_run.o

letkf_run.exe : $(objs)
	$(FC) $^ -o $@ $(LDFLAGS)
%.o : %.F90
	$(FC) $(FCFLAGS) -c $^
clean:
	rm -f *.o *.mod letkf_run.exe
