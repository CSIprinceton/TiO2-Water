## Sample analysis for the TiO2-water trajectory

Here I provide a simple [code](./codes/RDF_tio2.f90) to compute the RDF between surface TiO2 atoms and water atoms. I am assuming you have ran lammps and generated the trajectory file "pos.lammpstrj". If so, you can run the analysis in this directory with:

```
./codes/compute_rdf.sh
```

The code will generate 3 rdf files. I included my own output files in [reference](./reference)

I have also provided a [code](./codes/ZDF.f90) to compute the atomic density distribution as a function of the direction normal to the TiO2 surface. This code will output the number density of Ti (1), H (2) and O(3) atoms relative to the experimental equilibrium oxygen number density of bulk water at ambient conditions. You can run this analysis with:

```
./codes/compute_number_density_distribution.sh
```

You should run the command above in the same directory where the lammps trajectory (pos.lammpstrj) is located.
