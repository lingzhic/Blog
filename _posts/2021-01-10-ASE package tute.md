---
layout: article
title: ASE package tute
tags: MD
# cover: assets/images/fuji400.jpg
article_header:
  type: overlay
  theme: dark
  background_color: '#203028'
  background_image:
    src: assets/images/fuji400.jpg
---

# ase package
https://wiki.fysik.dtu.dk/ase/ase/ase.html
# ase.lattice.cubic
https://wiki.fysik.dtu.dk/ase/ase/lattice.html?highlight=facecenteredcubic \
xtal -> abbreviation for crystal


```python
import numpy as np
import ase.calculators.lammpsrun as aselammps
from ase.io import *
from ase.lattice.cubic import *
```

## Create Lattice_object (supercell)
To set up a supercell of FCC copper with the [1,-1,0] direction along the x-axis, [1,1,-2] along the y-axis and [1,1,1] along the z-axis, use:


```python
from ase.lattice.cubic import FaceCenteredCubic
atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]], 
                          size=(2,2,3), 
                          symbol='Cu', 
                          pbc=(1,1,0), 
                          latticeconstant=4.0)
print(type(atoms))
```

    <class 'ase.lattice.bravais.Lattice'>


The minimal unit cell is repeated `2*2*3` times.\
The lattice constant `can be specified` or default taken from the database of lattice constants in `ase.data` module. \
There are periodic boundary conditions along the x and y axis, but free boundary conditions along the z axis. \
Since the three directions are perpendicular, a (111) surface is created.

## input arguments
### directions
Specify 3 crystal directions that along x, y, and z axis
### size
A tuple of three numbers, defining how many times the fundamental repeat unit is repeated. \
Default: (1,1,1). Be aware that if high-index directions are specified, the fundamental repeat unit may be large.
### symbol
The element, specified by the atomic number (an integer) or by the atomic symbol (i.e. ‘Au’). \
For compounds, a tuple or list of elements should be given. This argument is mandatory.
### pbc
Periodic boundary conditions\
pbc=(a, b, c)
* 0 -> free boundary condition
* 1 -> periodic boundary condition


```python
from ase import Atoms
d = 1.1
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])
print(co)
```

    Atoms(symbols='CO', pbc=False)



```python
# np.ceil(variable) -> to get the ceiling of a float number or a list of float
num = 1.2
np.ceil(num)
```




    2.0



# Lattice_object methods
https://wiki.fysik.dtu.dk/ase/ase/atoms.html
## Lattice_object.get_cell()
https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.get_cell \
Returns the three `unit cell vectors`, as a `ase.cell.Cell` object.\
The Cell object resembles a 3x3 ndarray, and cell[i, j] is the jth Cartesian coordinate of the ith cell vector.


```python
from ase.lattice.cubic import FaceCenteredCubic
atoms = FaceCenteredCubic(directions=[[1,0,0], [0,1,0], [0,0,1]], 
                          size=(3,2,1), 
                          symbol='Cu', 
                          pbc=(0,0,0), 
                          latticeconstant=4.0)
xm = atoms.get_cell() / 2.0
print(xm)
print(xm[2, 2])
```

    [[6. 0. 0.]
     [0. 4. 0.]
     [0. 0. 2.]]
    2.0


## Lattice_object.get_positions()
Returns an array of atom positions.\
https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.get_positions


```python
print(atoms.get_positions())
```

    [[-7.34846923 -2.82842712 -3.46410162]
     [-6.53197265 -1.41421356 -1.15470054]
     [-5.71547607 -2.82842712  1.15470054]
     [-4.89897949 -1.41421356 -3.46410162]
     [-4.0824829  -2.82842712 -1.15470054]
     [-3.26598632 -1.41421356  1.15470054]
     [-2.44948974 -2.82842712 -3.46410162]
     [-1.63299316 -1.41421356 -1.15470054]
     [-0.81649658 -2.82842712  1.15470054]
     [ 0.         -1.41421356 -3.46410162]
     [ 0.81649658 -2.82842712 -1.15470054]
     [ 1.63299316 -1.41421356  1.15470054]
     [ 2.44948974 -2.82842712 -3.46410162]
     [ 3.26598632 -1.41421356 -1.15470054]
     [ 4.0824829  -2.82842712  1.15470054]
     [ 4.89897949 -1.41421356 -3.46410162]
     [ 5.71547607 -2.82842712 -1.15470054]
     [ 6.53197265 -1.41421356  1.15470054]
     [-7.34846923  0.         -3.46410162]
     [-6.53197265  1.41421356 -1.15470054]
     [-5.71547607  0.          1.15470054]
     [-4.89897949  1.41421356 -3.46410162]
     [-4.0824829   0.         -1.15470054]
     [-3.26598632  1.41421356  1.15470054]
     [-2.44948974  0.         -3.46410162]
     [-1.63299316  1.41421356 -1.15470054]
     [-0.81649658  0.          1.15470054]
     [ 0.          1.41421356 -3.46410162]
     [ 0.81649658  0.         -1.15470054]
     [ 1.63299316  1.41421356  1.15470054]
     [ 2.44948974  0.         -3.46410162]
     [ 3.26598632  1.41421356 -1.15470054]
     [ 4.0824829   0.          1.15470054]
     [ 4.89897949  1.41421356 -3.46410162]
     [ 5.71547607  0.         -1.15470054]
     [ 6.53197265  1.41421356  1.15470054]]


## Lattice_object.position()
Edit position of the whole object directly\
https://wiki.fysik.dtu.dk/ase/ase/atoms.html?highlight=position#ase.Atoms.positions

## Lattice_object.translate((x, y, z))
Translate atomic positions. Non-return method, modify Lattice_object directly\
The displacement argument can be a float an xyz vector or an nx3 array (where n is the number of atoms).\
https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.translate


```python
from ase.lattice.cubic import FaceCenteredCubic
atoms = FaceCenteredCubic(directions=[[2,-1,-1], [0,1,-1], [1,1,1]], 
                          size=(3,2,1), 
                          symbol='Au', 
                          pbc=(0,0,0), 
                          latticeconstant=4.0)
xm = atoms.get_cell() / 2.0
print(xm)
print(xm[2, 2])
```

    [[7.34846923 0.         0.        ]
     [0.         2.82842712 0.        ]
     [0.         0.         3.46410162]]
    3.464101615137755



```python
print('Before using translate method\n', atoms)
atoms1 = atoms
atoms.translate((-xm[0, 0], -xm[1, 1], -xm[2, 2]))
print('After using translate method\n', atoms)
atoms2 = atoms
print(atoms1 == atoms2)
```

    Before using translate method
     Lattice(symbols='Au36', pbc=False, cell=[14.696938456699069, 5.656854249492381, 6.92820323027551])
    After using translate method
     Lattice(symbols='Au36', pbc=False, cell=[14.696938456699069, 5.656854249492381, 6.92820323027551])
    True


## Lattice_object.get_center_of_mass()
get the mass center position of a structure


```python
print(atoms.get_center_of_mass())
```

    [1. 1. 1.]



```python
atoms.get_cell()[2, 2]
```




    8.0



## Lattice_object.get_global_number_of_atoms()
Returns the total number of atoms in the lattice object\
https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.get_global_number_of_atoms

# Read and Write data file
https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html?highlight=write_lammps_data#ase.io.lammpsdata.write_lammps_data

## ase.io.lammpsdata.read_lammps_data(filename)

## ase.calculators.lammpsrun.write_lammps_data(filename, atoms_object)
Just use this one, ase.io.lammpsdata.write_lammps_data does not work for some reasons


```python
from ase.calculators.lammpsrun import *
help(ase.calculators.lammpsrun.write_lammps_data)
```

    Help on function write_lammps_data in module ase.io.lammpsdata:
    
    write_lammps_data(fileobj, atoms, specorder=None, force_skew=False, prismobj=None, velocities=False, units='metal', atom_style='atomic')
        Write atomic structure data to a LAMMPS data file.
    



```python
import ase.io.lammpsdata
dir(ase.io.lammpsdata)
ase.io.lammpsdata.write_lammps_data(fil)
```




    ['Atoms',
     'Prism',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__spec__',
     'convert',
     'np',
     'paropen',
     're',
     'read_lammps_data',
     'write_lammps_data']