# Sample of analysis for TiO2/Water

Here I provide a simple [code](./codes/RDF_tio2.f90) to compute the RDF between surface TiO2 atoms and water atoms. I am assuming you have ran lammps and generated the trajectory file "pos.lammpstrj". If so, you can run the analysis in this directory with:

```
./codes/compute_rdf.sh
```

The code will generate 3 rdf files. I included my own output files in [reference](./reference)
