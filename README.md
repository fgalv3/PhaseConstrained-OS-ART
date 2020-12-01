# PhaseProjected-OS-ART

This code is posted here for completeness and to comply with reproducibility in science.
It refers to the paper "Model-driven reconstruction with phase-constrained highly-oversampled MRI", Galve F, Alonso J, Algarin JM, Benlloch JM.
about algebraic reconstruction technique applied to MRI when oversampled acquisition is used (timestep less than Nyquist) and projection of the 
image to real values is used. 

Two kind of MRI k-space trajectories are used: Echo-Planar Imaging (EPI) and Spiral trajectories, both with/without acceleration. 
Examples are given in program.c, which ouputs kspace, signal and reconstructed image.

The code needs armadillo-c++ library to work.


