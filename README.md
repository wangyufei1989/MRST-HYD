# MRST-HYD

This repository contains the file for library: MRST-HYD, which is a toolbox for simulating offshore CO2 storage with hydrate formation


**Installation**

MRST-HYD is handy, and no explicit operation is needed to install MRST-HYD. One only need to go into folder 'MRST-HYD/benchmarks and examples
/' and find the example file and run it.

**Quick start**

* Go to the folder 'MRST-HYD/benchmarks and examples/benchmark_against_Geoxim/'
* Open the Matlab script 'benchmark_against_geoxim.m'
* Clic on run

Or in command window

* cd 'MRST-HYD/benchmarks and examples/benchmark_against_Geoxim/'
* benchmark_against_geoxim

Optional 
*  Modify the parameters in the 'benchmark_against_geoxim.m'

*  For more details go to paper Wang et al [2023]

**Example simulations**

We offer four example simulations to help users to learn to use MRST-HYD. These four examples are listed in the folder called 'benchmarks and examples'. For instance, one goes to 'MRST-HYD/benchmarks and examples
/benchmark_against_Geoxim/' and run the file 'benchmark_against_geoxim.m'; the user can simulate their own case by modifying the parameters in the 'benchmark_against_geoxim.m'. For more details go to paper Wang et al [2023]

![plot](./benchmarks and examples/benchmark_against_Geoxim/BENCH_geoxim.jpg?raw=true)

**Platform**

The code is developed based on the platform called MRST that can be dowloaded through https://www.sintef.no/projectweb/mrst [Lie,2019]. We develop the reactive transport module considering hydrate formation/dissociation,  enclose it in the folder called 'partial-equilibrium-reactive-three-phase', and test it with the test files in the folder called 'benchmarks and examples', while the files in the rest folders are from  MRST.

**Software requirement**

Only Matlab is needed. The current module is developed based on Matlab 2019b, and it should be compatible for new versions of Matlab. 

**Acknowledgement**

This work is fundded by the postsoc project of IFP Energies Nouvelles (France). We acknowledge the discussion with Dr. Maarten W. Saaltink (UPC- Barcelona Tech) during the early modeling stage.

**References**
Lie, K. A. (2019). An introduction to reservoir simulation using MATLAB/GNU Octave: User guide for the MATLAB Reservoir Simulation Toolbox (MRST). Cambridge University Press.
Yufei Wang, Eric Flauraud, Anthony Michel, Véronique Lachet and Clémentine Meiller (2023). Modeling Offshore Geological Carbon Storage (GCS) Undergoing Hydration with Matlab, submitted to Computational Geosciences.
