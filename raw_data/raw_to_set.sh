#!/bin/bash

for i in tio2 tio2_?water_vac tio2_water water
do

  cd $i
  python -c 'import numpy as np; np.save ("coord",np.loadtxt("coord.raw" ).astype(np.float32))'
  python -c 'import numpy as np; np.save ("box",np.loadtxt("box.raw" ).astype(np.float32))'
  python -c 'import numpy as np; np.save ("energy",np.loadtxt("energy.raw" ).astype(np.float32))'
  python -c 'import numpy as np; np.save ("force",np.loadtxt("force.raw" ).astype(np.float32))'
  
  mkdir set.000
  mv *npy set.000
  cd ..
done
