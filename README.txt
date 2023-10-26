RUNNING CELL RANGER CITE SEQ 




Step1

Copy over the ADT and GEX and any other libraries into a folder.

You will have folders like the following with fastq's in them,
ADT1, ADT2, ADT3 etc….
GEX1, GEX2, GEX3 etc….

Step2 
Create config files using:
create_config_file.pl
Make sure to change the options for you setup eg location of reference files

You will get config files like this:
multi_config_1.csv, multi_config_2.csv..


Step3
Run cellranger config files with the following command on HPC, make the necessary edits in the cellranger.sh for your setup.

qsub cellranger.sh 


Misc. Seurat functions wrappers for downstream processing
function_cite.R


For downstream processing see scripts:

Older Methods:
https://github.com/foocheung/H5N1_Single_cell_2022

New Methods:
https://github.com/foocheung/CHI-306/tree/main
