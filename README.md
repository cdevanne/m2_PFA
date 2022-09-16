# pfa_timing

try to separates two neighboring hadronic showers. 

Every function is explain is headers. (xxx.h)

The main is PFA_SDHCAL_TIMING.C

## How to use : 

in file pfa_timing
	
> $ cd m2_pfa/pfa_timing

start root and launch your file with dayas (examples of input rootfile are in /rootFiles

> $ root yourRootFile.root

in root : 

> root [0] tree->Process("PFA_SDHCAL_TIMING.C")


A new rootFile with the reconstruction of particles is created in : 

> $ rootFile/outputPFA.root

## setting (in setting.C)

If you want to print each event performance, set : 

> print = true;

If you want to plot  each event (saved in graph file), set : 

> plot = true;
