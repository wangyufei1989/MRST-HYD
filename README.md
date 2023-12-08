# MRST-HYD

This repository contains all the files for library: MRST-HYD, which is a toolbox for simulating offshore CO2 storage with hydrate formation. We offer four example case simulations for the purpose of introduction.  Users can learn MRST-HYD through these four given examples. Details about each exapmple are described inside the script.


![two_dimension_app](https://github.com/wangyufei1989/MRST-HYD/assets/97456379/17e76394-a813-421f-846b-ae150e6b0f7f)


**Installation**

MRST-HYD is handy, and no explicit operation is needed to install MRST-HYD. One only need to go into folder 'MRST-HYD/examples/', then go into one of the four subfolders each containing one example test, and choose the Matlab script inside this subfolder and run it.  


**Quick start (only for those who have zero knowledge of Matlab)**

* Launch Matlab Software
* Go to the folder 'MRST-HYD/examples/one_dimensional_Joule_Thomson_with_HYD_formation' 
* Open the Matlab script 'one_dimensional_Joule_Thomson_with_HYD_formation.m'
* Click on 'run' button

**Example simulations**



**Platform**

The code is developed based on the platform called MRST that can be dowloaded through https://www.sintef.no/projectweb/mrst [Lie,2019]. We develop the reactive transport module considering hydrate formation/dissociation,  enclose it in the folder called 'partial-equilibrium-reactive-three-phase', and test it with the test files in the folder called 'benchmarks and examples', while the files in the rest folders are from  MRST.

**Software requirement**

Only Matlab is needed. The current module is developed based on Matlab 2019b, and it should be compatible for new versions of Matlab. 

**Acknowledgement**

This work is fundded by the postsoc project of IFP Energies Nouvelles (France). We acknowledge the discussion with Dr. Maarten W. Saaltink (UPC- Barcelona Tech) during the early modeling stage.

**References**
Lie, K. A. (2019). An introduction to reservoir simulation using MATLAB/GNU Octave: User guide for the MATLAB Reservoir Simulation Toolbox (MRST). Cambridge University Press.
Yufei Wang, Eric Flauraud, Anthony Michel, Véronique Lachet and Clémentine Meiller (2023). Modeling Offshore Geological Carbon Storage (GCS) Undergoing Hydration with Matlab, submitted to Computational Geosciences.
