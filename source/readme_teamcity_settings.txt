Build SFINCS CPU EL7 Intel
Original:
	------------------------------
	Build step (1 of 2): Configure
	---
	mkdir -p $(pwd)/opt/sfincs
	module load intel/18.0.3
	module load netcdf/v4.7.4_v4.5.3_intel18.0.3
	
	autoreconf -ivf
	./autogen.sh &&
	CPPFLAGS=-I/opt/apps/netcdf/v4.7.4_v4.5.3_intel18.0.3/include FC=ifort F77=ifort ./configure --prefix $(pwd)/opt/sfincs
	
	------------------------------
	Build step (2 of 2): Make
	---
	module load intel/18.0.3
	module load netcdf/v4.7.4_v4.5.3_intel18.0.3
	make 
	make install
	rm -f opt.zip
	zip -r opt.zip  opt
	
	
Possibility to add -O3 optimization as in VS windows build settings:
	------------------------------
	Build step (1 of 2): Configure
	---
	mkdir -p $(pwd)/opt/sfincs
	module load intel/18.0.3
	module load netcdf/v4.7.4_v4.5.3_intel18.0.3
	
	autoreconf -ivf
	./autogen.sh &&
	CPPFLAGS=-I/opt/apps/netcdf/v4.7.4_v4.5.3_intel18.0.3/include FC=ifort F77=ifort FCFLAGS=-O3 ./configure --prefix $(pwd)/opt/sfincs
	
	------------------------------
	Build step (2 of 2): Make
	---
	module load intel/18.0.3
	module load netcdf/v4.7.4_v4.5.3_intel18.0.3
	make 
	make install
	rm -f opt.zip
	zip -r opt.zip  opt

