
ABOUT
This software is a kinetic Monte Carlo simulator. It extends the popular SPPARKS program and allows it to simulate more than one type of atom. 

SPPARKS        
This program is greatly aided by the work of code developed at Sandia National Laboratory. SPPARKS stands for Stochastic Parallel PARticle Kinetic Simulator.It performs kinetic Monte Carlo simulations and is distributed under the GPL license. Since this work is derivative, it is also subject to that license. To understand this code you will need to be familiar with the documentation of SPPARKS. The documentation can be found on Github.  

COMPILING     
To compile the code, type 'make serial' while in the src directory. This will create an executable called 'spk_serial'. An input script can be used by simply typing './spk_serial<in.spparks'. The only differences from the original src code are found in the app_diffusion.cpp file, however all the code is provided for ease of compiling. 

MODIFICATIONS TO SPPARKS
Modifications have been wrapped with comments //marshalltc to allow for easy searching. All changes are in the app_diffusion.cpp file. 

NEW COMMANDS
'deposition rate dirx diry dirz d0 lo hi pid'
The deposition command has been changed so that the last parameter is an integer that specifies what type of atom is deposited. 
'barrier hop atom1 atom2 ....atomN'
The barrier commmand has been changed to allow for the specification of different barriers for different atoms to diffuse. atom1 represents the barrier for an atom with pid = 1. You can specify as many barriers as you like. 
