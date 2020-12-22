# PhaseConstrained-OS-ART

This code is posted here for completeness and to comply with reproducibility in science.
It refers to the paper "Model-driven reconstruction with phase-constrained highly-oversampled MRI", Galve F, Alonso J, Algarin JM, Benlloch JM.
about algebraic reconstruction technique applied to MRI when oversampled acquisition (timestep less than Nyquist) and projection of the 
image to real values is used. 

Two kind of MRI k-space trajectories are used: Echo-Planar Imaging (EPI) and Spiral trajectories, both with/without acceleration. 
Four reconstruction methods are used: plane waves sum (equivalent to DFT for Nyquist grid), ART with/without phase projection, ART with l1-norm (approx.)
TV penalty. Examples are given in program.c, which ouputs kspace, signal and reconstructed image. A Mathematica program is given with plots of example 
reconstructions with a Shepp-Logan and Spiral trajectory.
One can now directly link with the static library libART.a

The code needs armadillo-c++ library to work.

Author: Fernando Galve, i3M, CSIC-UPV, Valencia, Spain.


