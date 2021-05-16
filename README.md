# TiO2-Water DP Trainig Data

If you're using this data, please read and cite: Calegari Andrade, M. F., Ko, H.-Y., Zhang, L., Car, R. & Selloni, A. Free energy of proton transfer at the water–TiO2 interface from ab initio deep potential molecular dynamics. Chem. Sci. 11, 2335–2341 (2020).

## Description

This data contains atomic coordinates, energy and forces of TiO2 anatase, liquid water and the TiO2-water interface used to train a Deep Neural Network Potential. Data set was collected via a iterative training scheme described in [2]. In this repo you will find 3 main folders:

1. lammps: example of an input script (lammps.in) used to run Deep Potential Molecular Dynamics (DPMD) using the Lammps package. The file "pos.in" contain a initial configuration of the TiO2-water interface;

2. raw_data: in this root folder you will find 5 other folders containing the raw data of different systems. These systems are: 
   2.1. "tio2/": bulk anatase TiO2 (162 atoms); 
   2.2. "tio2_1water_vac/": anatase (101) surface with one adsorbed water molecule in vacuum;
   2.3. "tio2_2water_vac/": anatase (101) surface with two adsorbed water molecules in vaccum;
   2.4. "tio2_water/": TiO2-water interface (426 atoms); 
   2.5. "water/": bulk liquid water (192 atoms)

3. train: DeepMD-kit [3] input files used to train a Deep Neural Network potential (smooth edition) [4]. There are directories inside this folder, named from 1 to 4. These folders differ only in the random initialization of the initial NN parameters. We used the prediction deviation from models 1-4 as a coarse error estimatior of forces (or energy).   ***NOTE: This input script was used with DeepMD-Kit version 1.1. Different versions of the code have different input flags. Please check DeepMD-Kit documentation for more details.***

4. PW: example of PWscf [5] input file used to compute energy and forces of TiO2-water. ***Please note: if you want to expand this training data you HAVE to use the same pseudo-potential, wave-function cutoff and DFT functional (SCAN) described in the PW input file. It is also higly recommended to use PWscf as your force and energy evaluator.***

## Note on the format of raw files

   Please check the DeepMD-kit documentation for more details about the raw files format. Briefly, each raw file contains the full information of one snapshot in one line. So, if a "prefix.raw" contains 100 lines it means that this file has 100 different snapshots. Below we give a brief description of each raw file in our data:

   1. coord.raw: atomic coordinates in angstrom units. Format: C(1,x) C(1,y) C(1,z) ... C(N,x) C(N,y) C(N,z). C(i,j) is the cartesion component j of atom with index i.
   2. force.raw: atomic forces in eV/angstrom units. Format: F(1,x) F(1,y) F(1,z) ... F(N,x) F(N,y) F(N,z). F(i,j) is the cartesion component j of atom with index i.
   3. energy.raw: potential energy in eV units.
   4. box.raw: unit cell tensor in angstrom units.
   5. type.raw: index assigned for each atomic species. Format: I(1) I(2) ... I(N). I(i) is the label of atomic species of atom with index i. In this data set we use the following convention: Ti=1, H=2 and O=3.

##References

[1] Calegari Andrade, M. F., Ko, H.-Y., Zhang, L., Car, R. & Selloni, A. Free energy of proton transfer at the water–TiO2 interface from ab initio deep potential molecular dynamics. Chem. Sci. 11, 2335–2341 (2020).

[2] Zhang, L., Lin, D.-Y., Wang, H., Car, R. & E, W. Active learning of uniformly accurate interatomic potentials for materials simulation. Phys. Rev. Mater. 3, 023804 (2019).

[3] Wang, H., Zhang, L., Han, J. & E, W. DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics. Comput. Phys. Commun. 228, 178–184 (2018).

[4] Zhang, L. et al. End-to-end Symmetry Preserving Inter-atomic Potential Energy Model for Finite and Extended Systems. in Advances in Neural Information Processing Systems 31 (eds. Bengio, S. et al.) 4436–4446 (Curran Associates, Inc., 2018).

[5] Giannozzi, P. et al. Quantum ESPRESSO: a Modular and Open-Source Software Project for Quantum Simulations of Materials. J. Phys. Condens. matter 21, 395502 (2009). 
