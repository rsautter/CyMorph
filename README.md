# CyMorph
### Non-parametric Morphology Pipeline
This software

## Compiling

    cd maskMaker/
    python compile.py build_ext --inplace
    cd ..
    python compile.py build_ext --inplace
 
## Running example
Single-core and single object:


    time python main.py config.ini
MPI run:


    mpirun -np 3 python PCyMorph.py spirals30.csv
    
## Configure File
In order to run, a config file is required (in the example the config.ini where used). To run with MPI support, the default configuration file is ParallelConfig.ini.
This configuration file contain information about files, the imput/output, and the parameters

## Downloading from SDSS
To download images from SDSS, a

