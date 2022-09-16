# pfa_timing

try to separates two neighboring hadronic showers

## How to use : 

in file pfa_timing
	
> $ cd m2_pfa/pfa_timing

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
