# primate (IV)SO analysis
 
 This repository hosts code and data files for replicating all statistical results reported in 
 
 Olivier, C.-A., Martin, J. S., Pilisi, C., Agnani, P., Kauffmann, C., Hayes, L., Jaeggi, A. V., & Schradin, C. (2024). Primate social organization evolved from a flexible pair-living ancestor. Proceedings of the National Academy of Sciences of the United States of America, 121, e2215401120. https://www.pnas.org/doi/full/10.1073/pnas.2215401120

Code will run most efficiently if all files are saved to a common directory. The `primate_ivso_prep.R` file can be used to see how the original database files (`db_R.csv`,`db_G.csv`,`db_GR.csv`) were modified for analysis, as well as how discrepancies between species names in the original database and phylogeny were identified and corrected. The `primate_ivso.R` file can then be used with the `dataR.csv`, `dataG.csv`, and `dataGR.csv` files for statistical analysis and visualization. 
