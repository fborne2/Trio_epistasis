## The folder contains: 

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



# Guideline (see Methods of the manuscript)

## Ancestral reconstruction with codeml (PAML 4.9)
codeml codeml.ctl
9310/9375 sequences gave ancestral sequences: either due to stop codon inferred in ancestral sequences or pb in the alignment

## Run multiple alignment on ZI + Ancestor sequences using macse (v. 2.07)
macse -prog alignSequences \
     -seq "$SEQFILE" \
     -max_refine_iter 1 \
     -out_NT "${BASENAME}_NT.fasta" \
     -out_AA "${BASENAME}_AA.fasta"
MACSE succesfully run on 9137 sequences


## Filter gaps (ancestor sequence contains bits that are not present in the Dmel sequences) 
Use python script filter_gap.py

## Create consensus sequences and scale down by 3 
Use python script create_consensus.py

## Split sequence > 900 by chunks of 900 (make sure about the header and save the full original sequence somewhere else) 
Use python script split_in_900.py
9137 proteins splitted into 10573

## Run MASS-PRF on all following https://github.com/Townsend-Lab-Yale/MASSPRF 
./massprf -ic 1 -sn 83 -p $transcript"_polymorphism_consensus.txt" -d $transcript"_divergence_consensus.txt" -o 1  -r 1 -ci_r 1 -ci_m 1 -s 1 -exact 0 -mn 30000 -t 2.5 >$transcript"_MASS-PRF_BIC.txt"

## Only keep table part of the output
for file in *.txt; do   awk '/^Position/{flag=1} /^Abbreviation:/{flag=0} flag' "$file" > "tables_only/table_$file"; done

## Create list of sequences that contains at least one cluster of adaptive substitutions 
use python script list_adaptive_cluster.py



