all:
	gfortran -c -cpp libanalysis.f90 
	gfortran -I./ -O2 RDF_tio2.f90 -o rdf_tio2.x ./libanalysis.o 
	gfortran -I./ -O2 ZDF.f90 -o zdf.x ./libanalysis.o
