I run MASS-PRF on all transcripts using Dmel-Dsim ancestor as the outgroup.
See scripts/Run_on_all_CDS for more details 

-iso1_correspondence_table
Correspondence between iso1 transcript found in the orthogroup and the OG group

- species_alignment: CDS alignment between Dmel, Dsim, Dyak, Dsan, Dtei used for ancestral reconstruction

- Dmel_Anc_alignments: CDS alignment between Zambia population and Dmel/Dsim ancestor used to create consensus sequences

- Dmel_Anc_consensus_sequences: consensus sequences used to run MASS-PRF

- MASS-PRF_output: contains raw results of MASS-PRF that run succesfully 

- MASS-PRF_tables_only contains table of interest from the MASS-PRF results. 

- MASS-PRF_adaptive_clusters: contains table of sequences that contain at least one cluster (several sites lower gamma>4 and less than 20 codons apart).

- scripts: contain pipeline and scripts used for the analysis

- results: contains MASS-PRF profiles of protein with adaptive clusters and a list of experimentally tractable candidates

- baseml_output: output from BASEML run
