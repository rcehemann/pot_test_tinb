#!/bin/bash

#for setnum in 02 03 04 05 06 07 08 09; do
for setnum in 01 02 03 04 05 06 07 08 09; do
	for potnum in 1 2 3 4 5 6 7 8 9 10; do
	dir=`pwd`
	cd /n/jww-1/ehemann.2/pots/GMEAM/Ti-Nb/set.$setnum/pot.$potnum/
	mv lammps.pt old_lammps.pt
	python /n/jww-1/ehemann.2/testingScripts/rescale_gmeam_3body.py old_lammps.pt > lammps.pt
	cd $dir

	./masterBinary.sc Ti Nb $setnum $potnum
done
done
