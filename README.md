# MRST-HYD

This repository contains the file for library: MRST-HYD, which is a toolbox for simulating offshore CO2 storage with hydrate formation.

![two_dimension_app](https://github.com/wangyufei1989/MRST-HYD/assets/97456379/a1e1d7bb-699f-482d-8e5e-e61a502cde4d)

**Installation**

MRST-HYD is handy, and no explicit operation is needed to install MRST-HYD. One only need to go into folder 'MRST-HYD/examples
/' and choose one example file and run it. We offer four example files for the purpose of introduction.  Users can learn MRST-HYD through these four given example files. Details about each exapmple are described inside the script. 


**Quick start**

* Launch Matlab Software
* Go to the folder 'MRST-HYD/examples/'
* Go to the folder 'Test_1', for instance, and open the Matlab script 'Test_1.m'
* Click on 'run' button





**Example simulations**

We offer four example simulations to help users to learn to use MRST-HYD. These four examples are listed in the folder called 'benchmarks and examples'. For instance, one goes to 'MRST-HYD/benchmarks and examples
/benchmark_against_Geoxim/' and run the file 'benchmark_against_geoxim.m'; the user can simulate their own case by modifying the parameters in the 'benchmark_against_geoxim.m'. For more details go to paper Wang et al [2023]

![plot]('https://github.com/wangyufei1989/MRST-HYD/params/BENCH_geoxim.JPG')

**Platform**

The code is developed based on the platform called MRST that can be dowloaded through https://www.sintef.no/projectweb/mrst [Lie,2019]. We develop the reactive transport module considering hydrate formation/dissociation,  enclose it in the folder called 'partial-equilibrium-reactive-three-phase', and test it with the test files in the folder called 'benchmarks and examples', while the files in the rest folders are from  MRST.

**Software requirement**

Only Matlab is needed. The current module is developed based on Matlab 2019b, and it should be compatible for new versions of Matlab. 

**Acknowledgement**

This work is fundded by the postsoc project of IFP Energies Nouvelles (France). We acknowledge the discussion with Dr. Maarten W. Saaltink (UPC- Barcelona Tech) during the early modeling stage.

**References**
Lie, K. A. (2019). An introduction to reservoir simulation using MATLAB/GNU Octave: User guide for the MATLAB Reservoir Simulation Toolbox (MRST). Cambridge University Press.
Yufei Wang, Eric Flauraud, Anthony Michel, Véronique Lachet and Clémentine Meiller (2023). Modeling Offshore Geological Carbon Storage (GCS) Undergoing Hydration with Matlab, submitted to Computational Geosciences.
