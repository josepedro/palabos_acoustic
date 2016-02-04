# Palabos Acoustic

Palabos Acoustic is a project derived by Palabos (http://palabos.org/) that aim to implement and simulate physical acoustic issues.

To run a example of acoustic implementation of Lattice Boltzmann in Palabos you can run:
$  make && time mpiexec anechoic_dynamics  &&  convert -delay 30 tmp/u*.gif tmp/animation.gif && animate tmp/animation.gif
