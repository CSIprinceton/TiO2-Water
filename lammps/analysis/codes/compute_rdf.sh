#Atom types
ti5c=15
o2c=32
ow=3
h=2

#Compile Fortran code
cd codes
make
cd ..

#Generate input file for Fortrain code
./codes/make_input_files.sh ../pos.lammpstrj

#Compute RDFs
echo $ti5c $ow | ./codes/rdf_tio2.x ../pos.lammpstrj index_rdf
mv rdf.dat rdf_ti5c_ow.dat
echo $o2c $h | ./codes/rdf_tio2.x ../pos.lammpstrj index_rdf
mv rdf.dat rdf_o2c_h.dat
echo $o2c $ow | ./codes/rdf_tio2.x ../pos.lammpstrj index_rdf
mv rdf.dat rdf_o2c_ow.dat
