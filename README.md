# LatPack: lattice protein software suite

LatPack is software originally written at the University of Freiburg
(http://www.bioinf.uni-freiburg.de/Software/LatPack/). This version
has received a few modifications for my research purposes:

- Random number generator uses Random123 library
  (https://www.deshawresearch.com/resources_random123.html)
- Simulation data (`latFold`, `latFoldVec`) can be output to HDF5 file
  format
- New program `latMapTraj` based on `latMap` to analyze simulation
  trajectories
- Ribosome feature in `latFoldVec` to imitate tethered C-terminus
- Conformation counting for more precise native conformation occupancy
  measurement

# Installation
One needs to install and compile BIU (included here) and ELL (Energy
Landscape Library) first. HDF5 also needs to be installed.

```
./configure --with-BIU=${PATHTOBIU} \
  --with-ELL=${PATHTOELL} \
  --with-HDF5=${HDF5PATH} \
  --enable-latFold CXXFLAGS="-std=c++0x"
make
make install

```

# Copyright
copyright 2020 by Victor Zhao
copyright 2008 by Martin Mann  (http://www.bioinf.uni-freiburg.de/~mmann/)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
