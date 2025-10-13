- Run_on_all_CDS.txt: pipeline used to run MASS-PRF

- file_category_v4.py: give counts for categories: 
1. Runs but lower gamma <4 for all sites 
2. Only one site with lower gamma >= 4
3. Several lower gamma >= 4 but >= codons apart
4. Several lower gamma >= 4 and at least two < 20 codons apart
and move category 4 into several_ge4_some_lt20 subfolder

- file_category_unique_OG: 
Adapt the script above to have one category per OG (in case OG is in multiple fragments). 

- list_adaptive_cluster.py: 
Create list of sequences that contains at least one cluster of adaptive substitutions (at least 2 substitutions less than 20 codons away)

- find_tractable_cluster_v2.py:
Will give a list of file and clusters that contain exactly 2 or 3 substitutions less than 20 codons away from each other

- make_plot_mass_prf_candidates_loop_pdf.R
go over a list of table file and plot MASS-PRF profile. in red, sites for which lower bond gamma>=4

- codeml: ctl command file used to run codeml

- treefile_v2.nwk: species tree used to run codeml

- create_consensus.py: python script used to create consensus sequences from Zambian and Dmel-Dsim ancestor alignment

- filter_gap.py: remove gaps from the Zambian and Dmel-Dsim ancestor alignment

- Split_in_900.py: split consensus sequences >900 codons into fragments of 900 pieces
