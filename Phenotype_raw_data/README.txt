- adult_hetero_Fisher_test.txt: count of YFP+ and YFP- of heterozygous adults. Data used to run Fisher exact test to know whether distribution of YFP+ and YFP- individuals are significantly different from distribution of the WT genotype ASK. 

- adult_homo_Fisher_test.txt: count of YFP+ and YFP- of homozygous adults. Data used to run Fisher exact test to know whether distribution of YFP+ and YFP- individuals are significantly different from distribution of the WT genotype ASK. 

- pupa_homo_Fisher_test.txt: count of YFP+ and YFP- of homozygous pupae. Data used to run Fisher exact test to know whether distribution of YFP+ and YFP- individuals are significantly different from distribution of the WT genotype ASK. 

- pupa_hetero_Fisher_test.txt: count of YFP+ and YFP- of heterozygous pupae. Data used to run Fisher exact test to know whether distribution of YFP+ and YFP- individuals are significantly different from distribution of the WT genotype ASK. 

- larva_hetero_Fisher_test.txt: count of YFP+ and YFP- of heterozygous larvae. Data used to run Fisher exact test to know whether distribution of YFP+ and YFP- individuals are significantly different from distribution of the WT genotype ASK. 

- larva_homo_Fisher_test.txt: count of YFP+ and YFP- of homozygous larvae. Data used to run Fisher exact test to know whether distribution of YFP+ and YFP- individuals are significantly different from distribution of the WT genotype ASK. 

- ratio_WT_homo.txt: ratio of YFP+ over YFP- normalized by expected ratio for homozygous larva, pupa, adult. Data used to plot Figure 2A. 

- ratio_WT_hetero.txt: ratio of YFP+ over YFP- normalized by expected ratio for heterozygous larva, pupa, adult. Data used to plot Figure 2A. 

- fertility.txt: number of progeny from individual crosses for ASK and VNR haplotypes ("number"). 

- climbing_assay.txt: raw data of climbing assay. We quantified locomotor performance by quantifying startle-induced negative geotaxis. Ten 1-day-old males were grouped with ten females into vials for 7-10 days (column "age"). Flies were let to recover for 40 to 60 min before the start of the experiment (column "recovery_time"). The assay tube was tapped 5 times on a fly pad so that all flies fell to the bottom of the tube, and the flies that reached the line in the following 10 sec were counted 3 times with 5 min recovery between each time (10s_1, 10s_2, 10s_3) and the mean (10s_mean) and the percent of the mean (mean_10s) were calcultated. The same assay tubes were then assayed the same way after 10 min recovery but counting flies that reached the line in the following 5 sec instead of 10. 


- Script_Fig: Contains scripts to make Figures and Stats: 
	- trio_climbing_assay.R: Plot and statistics for Fig 2B.
	- trio_fertilityR: plot and statistics for Fig S4
	- trio_viability.R: plot Fig 2A
	- viability_Fisher_tests: statistics for Fig 2A
