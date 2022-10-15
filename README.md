# TiO2-Water DP Trainig Data #

This repo contains the data used to train a DPMD potential for the TiO2-water interface, as described in 

*Calegari Andrade, M. F., Ko, H.-Y., Zhang, L., Car, R. & Selloni, A. Free energy of proton transfer at the water–TiO2 interface from ab initio deep potential molecular dynamics. Chem. Sci. 11, 2335–2341 (2020).*

Before you start, please make sure you have DeepMD-Kit and Lammps installed. You can find good tutorials on how to compile DeepMD-kit and Lammps at the [DeepMD-Kit Documentation Page](https://docs.deepmodeling.com/projects/deepmd/en/master/)

## Downloading the dataset ##

Please make sure you have git-lfs running on your local machine. The raw data used to train the DNN can only be downloaded with git-lfs.

```
git clone https://github.com/CSIprinceton/TiO2-Water.git
cd TiO2-Water
git lfs fetch --all
get lfs pull
```

## Training your own DP models ##

I have provided trained DNN graphs that are ready to use. These graphs are located at "train/?/". In case you want to train your own DNN model, please do the following:

```
cd raw_data
./raw_to_set.sh
cd ../train/1
dp train tio2-water.json
```

## Running a simulation with Lammps ##

There is a simple lammps input example at [lammps/lammps.in](lammps/lammps.in) and a initial condition at [lammps/pos.in](lammps/pos.in). This example runs a NVT simulation of the TiO2-water interface at 330 K. 3 different DNN models are used to evaluate the interatomic interactions, and the deviation between the quatities predicted by these potentials will be outputed to the file model_devi.out.

Running Lammps is simple

```
lmp < lammps.in > lammps.out
```

**Note: there is no guarantee that the DNN potential will predict reasoable energy and forces if you start your dynamics from a unphysical configuration (or too far from the configurations sampled in the training data).**

## Detailed description of the dataset ##

This data contains atomic coordinates, energy and forces of TiO2 anatase, liquid water and the TiO2-water interface used to train a Deep Neural Network Potential. Data set was collected via a iterative training scheme described in [2]. In this repo you will find 3 main folders:

1. lammps: example of an input script (lammps.in) used to run Deep Potential Molecular Dynamics (DPMD) using the Lammps package. The file "pos.in" contain a initial configuration of the TiO2-water interface;

2. raw_data: in this root folder you will find 5 directories containing the raw data of different systems. These systems are: 
   * "tio2/": bulk anatase TiO2 (162 atoms); 
   * "tio2_1water_vac/": anatase (101) surface with one adsorbed water molecule in vacuum (219 atoms);
   * "tio2_2water_vac/": anatase (101) surface with two adsorbed water molecules in vaccum (222 atoms);
   * "tio2_water/": TiO2-water interface (426 atoms); 
   * "water/": bulk liquid water (192 atoms)

3. train: [DeepMD-kit](https://github.com/deepmodeling/deepmd-kit) [3] input files used to train a Deep Neural Network potential (smooth edition) [4]. There are directories inside this folder, named from 1 to 4. These folders differ only in the random initialization of the initial NN parameters. We used the prediction deviation from models 1-4 as a coarse error estimatior of forces (or energy).   *NOTE: This input script was used with DeepMD-Kit version 1.2. Different versions of the code have different input flags. Please check [DeepMD-kit documentation](https://deepmd.readthedocs.io/en/master/) for more details.*

4. PW: example of PWscf [5] input file used to compute energy and forces of TiO2-water. *Please note: if you want to expand this training data you HAVE to use the same pseudo-potential, wave-function cutoff and DFT functional (SCAN) described in the PW input file. It is also higly recommended to use PWscf as your force and energy evaluator.*

## Note on the format of raw files ##

   Please check the [DeepMD-kit documentation](https://deepmd.readthedocs.io/en/master/) for more details about the raw files format. Briefly, each raw file contains the full information of one snapshot in one line. So, if a "prefix.raw" contains 100 lines it means that this file has 100 different snapshots. Below we give a brief description of each raw file in our data:

   1. coord.raw: atomic coordinates in angstrom units. Format: C(1,x) C(1,y) C(1,z) ... C(N,x) C(N,y) C(N,z). C(i,j) is the cartesion component j of atom with index i.
   2. force.raw: atomic forces in eV/angstrom units. Format: F(1,x) F(1,y) F(1,z) ... F(N,x) F(N,y) F(N,z). F(i,j) is the cartesion component j of atom with index i.
   3. energy.raw: potential energy in eV units.
   4. box.raw: unit cell tensor in angstrom units.
   5. type.raw: index assigned for each atomic species. Format: I(1) I(2) ... I(N). I(i) is the label of atomic species of atom with index i. In this data set we use the following convention: Ti=0, H=1 and O=2.

## References ##

[1] Calegari Andrade, M. F., Ko, H.-Y., Zhang, L., Car, R. & Selloni, A. Free energy of proton transfer at the water–TiO2 interface from ab initio deep potential molecular dynamics. Chem. Sci. 11, 2335–2341 (2020).

[2] Zhang, L., Lin, D.-Y., Wang, H., Car, R. & E, W. Active learning of uniformly accurate interatomic potentials for materials simulation. Phys. Rev. Mater. 3, 023804 (2019).

[3] Wang, H., Zhang, L., Han, J. & E, W. DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics. Comput. Phys. Commun. 228, 178–184 (2018).

[4] Zhang, L. et al. End-to-end Symmetry Preserving Inter-atomic Potential Energy Model for Finite and Extended Systems. in Advances in Neural Information Processing Systems 31 (eds. Bengio, S. et al.) 4436–4446 (Curran Associates, Inc., 2018).

[5] Giannozzi, P. et al. Quantum ESPRESSO: a Modular and Open-Source Software Project for Quantum Simulations of Materials. J. Phys. Condens. matter 21, 395502 (2009). 
