% duct radiation

cd ~/palabos_acoustic/examples/codesByTopic/duct_radiation/
system('rm tmp/*');
system('make && time mpirun -n 6 duct_radiation');
