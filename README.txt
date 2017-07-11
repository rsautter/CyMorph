To Compile:
	cd output
	python compile.py build_ext --inplace
	cd ..
	python compile.py build_ext --inplace

Running examples:
	time python main.py config.ini
	mpirun -np 3 python PCyMorph.py spirals30.csv


