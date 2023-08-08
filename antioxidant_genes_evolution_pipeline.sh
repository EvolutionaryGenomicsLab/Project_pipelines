#This pipeline was created by Giovanna Selleghin Veiga during her IC project that aimed to study positive selection in antioxidant enzymes.

##Setting up the directories
mkdir -p project_directory/{src,data,results,sandbox}
mkdir data/raw_seqs/

##Get CDS sequences from NCBI
###The script for this stage is in the get_sequences repository
bash get_CDS_sequence_from_genes.sh genes_list.txt species_list.txt > data/raw_seqs/ #Not all genes are going to be found. Sometimes the acronym used by NCBI is different from the one you are using (like using cap letters) or the sequence is not annotated for the species. 

##Extracting the primary transcripts from NCBI files. Primary transcripts are the longest isoform and therefore used as the representative of the whole gene.
#Linearize the archives of CDS sequences
for files in $(cat data/raw_seqs/*.fasta);do awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${files} > data/linear_raw_seqs/${files}; done
# Remove the empty archives
grep -cH ">" data/linear_raw_seqs/*.fasta | awk -F":" '{ if($2 == 0) {print}}' | cut -d":" -f1 > to_remove.list
for file in $(cat to_remove.list) ; do rm ${file}; done
#Select the sequences in which the gene name is in the correct field ([gene=gene_name])
for file in $(ls data/linear_raw_seqs/*.fasta); do PATTERN=$(awk '{print $2}' ${file} | head -n 1) && grep -F ${PATTERN} ${file} > data/linear_raw_seqs/$(basename ${file}_filtered.fasta) ; done 
#Extract the biggest sequences from the filtered files
for file in $(ls data/linear_raw_seqs/*_filtered.fasta); do cat ${file} | awk '{ print length($(NF)), $0 }' | sort -n -r | cut -d" " -f2-| awk '/^>/{if(N)exit;++N;} {print;}' > data/longest_cds/$(basename ${file}); done
#Check if every file from the previous step has one sequence
for file in $(ls data/longest_cds/*.fasta); do grep -c ">" ${file}; done
#Trnasform into 60 characters in each line
for file in $(ls data/longest_cds/*.fasta); do cat ${file} | tr "\t" "\n" > data/CDS_sequences/$(basename ${file}); done
#Make the header simple 
for file in $(data/CDS_sequences/*.fasta); do sed "s/^>.*$/>${file%%.*}/; s/>[A-Z0-9]\{3,5\}_/>/" ${file}; done

##Alignments
mkdir data/seqs_to_align/
for gene in $(cat genes_list.txt); do cat data/CDS_sequences/${gene}* >> data/seqs_to_align/${gene}.fasta ; done
#Using MAFFT  
for file in $(ls data/seqs_to_align/*.fasta); do mafft --quiet ${file} > results/alignments/$(basename ${file%.fasta})_mafft.fasta; done
#Protein alignment (get protein sequences using Geneious)
for file in $(ls data/seqs_to_align/translated/*.fasta); do mafft --amino --quiet ${file} > results/alignments/$(basename ${file%.fasta})_mafft_protein.fasta; done
#Alignment using Muscle (protein and nucleotide)
for file in $(ls data/seqs_to_align/*.fasta); do ~/../software/muscle/muscle3.8.31_i86linux64  -quiet -in ${file} -out results/alignments/$(basename ${file%.fasta})_muscle.fasta; done

##Building phylogenetic trees
mkdir -p results/IQTree/{CAT,GPX3,GSR,PRDX1,PRDX3,SOD1,XDH}
#Copy alignments files into each directory above
for file in $(cat gene_list.txt); do cp results/alignments/${file}*fasta results/IQTree/${file}; done
#Run IQTree for each aligment
nohup iqtree -s results/IQTree/gene/gene_alignment.fasta -bb 1000 -wbt -nt AUTO &

##Condon alignment using PRANK http://wasabiapp.org/software/prank/ (input PAML)
#Run for each gene sequence
nohup ~/../software/prank/bin/prank -d=data/seqs_to_align/gene.fasta -o=result/alignments/gene_codon_alignment -codon -F &
#Convert to PAML input
~/../software/prank/bin/prank -convert -d=result/alignments/gene_codon_alignment.best.fas -f=paml -o=result/alignments/input_PAML_gene_codon_alignment.phy

##Running PAML
#Coping the input files to PAML directory
mkdir -p results/PAML/{CAT,GPX3,GSR,PRDX1,PRDX3,SOD1,XDH}
for file in $(cat gene_list.txt); do cp result/alignments/input_PAML_ ${file}_codon_alignment.phy results/PAML/${file}; done
#Run PAML
nohup codeml control_file.ctl &
#Script extract PAML infos from output (lnL, np, omegas and positively selected sites). Script in the repository extract_data
bash extract_info_PAML.sh
#cont sites found by BEB
for file in $(cat lista_genes.txt); do echo "${file}/${file}_modelA" && cat ${file}/${file}_modelA | grep -A50 "(BEB)" | g
rep -c -E "([0-9]\s[A-Z]\s[0-1]\.[0-9]|[0-9]\s(-)\s[0-1]\.[0-9])"; done
# Find sites with PP >= 0.900 in BEB
cat sitios.txt | awk ' $3 >= 0.900 {print $0}' 
#Calculate chi-square
chi2 "df" "LRT"
chi2 p

##Running HyPhY
#Activate conda environment
conda activate hyphy
#Basic command
hyphy method --alignment --options
#Infos about the methods 
hyphy method --help
#Running in an interactive way to get more parameter modifications
hyphy -i
#In the above mode you can select the method you are going to run and the available options to be changed. In the help command you can find the default values.
#Run the HyPhy inside the screen 
#show open screens
secreen -ls
#Attache to a screen
screen -r name.screen
#Detache from a screen
ctrl+A+D

##Run PartitionFinder
#Attached to a screen
screen -r nome.screen
#Activate python 2.7 environments used to run Partition Finder
conda activate python27
#Run the below command inside a directory with an alignment.phy and control file from Partition Finder. 
python partitionfinder-2.1.1/PartitionFinder.py -p 5 directory

##Run MrBayes
#The archive needs to have a specific block of command with MrBayes instructions (you get this from Partition Finder output)
nohup mb –i archive.nex &

##Install and run Godon
#Install (https://github.com/idavydov/godon-tutorial)
wget https://bitbucket.org/Davydov/godon/downloads/godon-master-linux-gnu-x86_64 -O godon
chmod +x godon
#Run
~/softwares/./godon test BS --m0-tree --ncat-codon-rate 4 codon_alignment.fasta tree.tre --procs 5 --json gene_result_GODON.json
# Use "#" to mark the trees  
#Extract the results
#In the results directory, run the script get_godon_results.py in the same directory. The script is in the repository extract_data
python get_godon_results.py > results_godon
