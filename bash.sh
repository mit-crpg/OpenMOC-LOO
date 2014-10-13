#!/bin/bash          
./bin/openmoc -wc -lp >> cmfd-lp.dat
./bin/openmoc -wc -ep >> cmfd-ep.dat
./bin/openmoc -wl1 >> loo1-np.dat
./bin/openmoc -wl1 -lp >> loo1-lp.dat
./bin/openmoc -wl1 -ep >> loo1-ep.dat
./bin/openmoc -wl2 >> loo2-np.dat
./bin/openmoc -wl2 -lp >> loo2-lp.dat
./bin/openmoc -wl2 -ep >> loo2-ep.dat
