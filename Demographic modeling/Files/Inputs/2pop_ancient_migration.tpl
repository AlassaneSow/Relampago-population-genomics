//Input parameters for the coalescence and recombination simulation program : simcoal2.exe
2
/Deme sizes (haploid number of genes)
NPOP1$
NPOP2$
//Sample Sizes - the HAPLOID sample size and generations since sampling (always 0)
68 0
68 0
//Growth rates
0
0
//Number of migration matrices : If 0 : No migration between demes: NoMig estimates the migration rates between populations
2
//NoMig matrix 0 :
0.000 0.000 
0.000 0.000
//NoMig matrix 1 :
0.000 MIG0-1$ 
MIG1-0$ 0.000
//Historical event: time (in gen), source, sink, proportion of migrants, new deme size, new growth rate, new migration matrix. MAKE SURE IN EST FILE T_STOP < T_DIV
2 historical events
T_STOP$ 0 0 1 1 0 1
T_DIV$ 0 1 1 1 0 keep
//Number of independent (unlinked) chromosomes, and "chromosome structure" flag:  0 for identical structure across chromosomes, and 1 for different structures on different chromosomes.
1 0
//Number of contiguous linkage blocks in chromosome 1:
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 1e-8 OUTEXP
