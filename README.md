# CyMorph
Non-parametric Morphology Pipeline

## Compiling

 python compile.py build_ext --inplace
 cd ..
 python compile.py build_ext --inplace
 
## Running example
Single-core and single object:
    time python main.py config.ini
Multi-core:
    mpirun -np 3 python PCyMorph.py spirals30.csv
## Downloading from SDSS


