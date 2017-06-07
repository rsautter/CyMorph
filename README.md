# CyMorph
### Non-parametric Galaxy Morphology Pipeline
This pipeline  

## Requirements
 - Python 2.7
 - MPI / MPI4py
 - Astropy
 - Numpy
 - SciPy
 - Cython

## Compiling

    cd maskMaker/
    python compile.py build_ext --inplace
    cd ..
    python compile.py build_ext --inplace
    
### Paths
There should be specified in 'cfg/paths.ini', where is/how to call R, python, and SExtractor
 
## Running example
Single-core and single object:


    time python main.py config.ini
MPI run:


    mpirun -np 3 python PCyMorph.py spirals30.csv
    
## Configure File
In order to run, a config file is required (in the example the config.ini where used). To run with MPI support, the default configuration file is ParallelConfig.ini.
This configuration file contain information about files, the input/output, and the parameters.
An example is:

    [File_Configuration]
    indexes: C, H, A3, S3, Ga, OGa
    cleanit: False
    download: True
    stamp_size: 5

    [Output_Configuration]
    verbose: False
    savefigure: False

    [Indexes_Configuration]
    Entropy_Bins: 180
    Ga_Tolerance: 0.02
    Ga_Angular_Tolerance: 0.02
    Ga_Position_Tolerance: 0.00
    Concentration_Density: 100
    Concentration_Distances: 0.65, 0.25
    butterworth_order: 2
    smooth_degree: 0.2

## Downloading from SDSS
To download images from SDSS, a module in Download was created. Is required a csv file.
