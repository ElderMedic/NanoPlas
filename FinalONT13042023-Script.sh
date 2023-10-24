# Create folder for the experimentation, the folder with the date will be selected in MinKnown for storing the RAW files 
# The program will create subfolders with the date, no_sample, random number and fast5 folder, we have to change this everytime in the input of the basecalling or open that Path before basecalling
mkdir -p ~/Documents/`date '+%y%m%d'`/{Basecalling,Barcoding,Results}

# Basecalling with Guppy from command line, this will generate two folders under the folder created "pass & fail" + basecaller logs. If we add the barcode kit it will generate folders called BarcodeXX with pass&fail subfolders
# The extension for the files created is .fastq.gz
guppy_basecaller -i ~/Documents/13042023/fast5 -s ~/Documents/Basecalling -c dna_r10.4.1_e8.2_260bps_hac.cfg -x auto -r --detect_adapter --detect_barcodes --detect_mid_strand_adapter --do_read_splitting 
--trim_adapters --min_qscore 10 --compress_fastq (--barcode_kits "SQK-RBK114-96" --detect_mid_strand_barcodes --enable_trim_barcodes --index --bam_out --sample_sheet) 

# Demultiplexing from the basecalled files 
guppy_barcoder -i ~/Documents/13042023/Basecalling/pass -s ~/Documents/13042023/Barcoding -c configuration.cfg --detect_barcodes --detect_mid_strand_barcodes -x auto -r --barcode_kits SQK-RBK114-96 --fastq_out --compress_fastq --enable_trim_barcodes 
(--index --bam_out)

# Concatenate the produced fastq.gz files  (Here we still have to modify where the directory will be)
for dir in ~/Documents/13042023/Barcoding/barcode*; do barcode_id=${dir: -2}; cd $dir && cat *.fastq.gz>b${barcode_id}.fastq.gz; done

# Plot the histogram from the concatenated files. The one that interest us has the name Non_weightedHistogramReadlength.png
for dir in ~/Documents/13042023/Barcoding/barcode*; do barcode_id=${dir: -2}; NanoPlot --fastq ~/Documents/13042023/Barcoding/barcode${barcode_id}/b${Pbarcode_id}.fastq.gz -o ~/Documents/13042023/Barcoding/barcode${barcode_id}; done

# De novo assembly using Flye + Medaka Polishing (maybe medaka is not needed if the quality from the sequencing is good)
# This step is difficult to do it in loop because of the genome size and maybe the coverage that depends on the sequencing experiment, maybe we can put similar size plasmids together and then it will be possible 
flye --nano-hq ~/Documents/13042023/Barcoding/barcode83/b83.fastq.gz --genome-size 3k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode83/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode83/b83.fastq.gz -d ~/Documents/13042023/Barcoding/barcode83/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode83/ -m r1041_e82_260bps_hac_g632
for P in 84 85 89; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 6.5k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632; done
for P in 81 82 86; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 9k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632; done
for P in 87 88 90; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 11.5k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632; done

for P in {84,85,89}; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode${P}/b${P}.fastq.gz --genome-size 6.5k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode${P}/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode${P}/b${P}.fastq.gz -d ~/Documents/13042023/Barcoding/barcode${P}/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode${P}/ -m r1041_e82_260bps_hac_g632; done

# Annotation using pLannotate (This requires to enter to the directory of the files and activate it)
for P in {81..90}; do cd ~/Documents/13042023/Barcoding/barcode${P} &&
conda activate plannotate && 
plannotate batch -i consensus.fasta --html && conda deactivate && cd ~/ ; done 

# Move the files to the specific Result folder created for each sample 
# cp copy the concatenated files and the files of interest, don't forget to rename them 
for P in {81..90}; do cd ~/Documents/Barcoding/barcode${P} && mkdir Results && cp b${P}.fastq.gz Non_weightedHistogramReadlength.png consensus.fasta consensus_pLann.html consensus_pLann.gbk ~/Documents/Barcoding/barcode${P}/Results 


for P in {81..90}; do cd ~/Documents/Barcoding/barcode${P} && mkdir Results && cp b${P}.fastq.gz ~/Documents/Barcoding/barcode${P}/Results && cp Non_weightedHistogramReadlength.png ~/Documents/Barcoding/barcode${P}/Results && 
cp consensus.fasta ~/Documents/Barcoding/barcode${P}/Results && cp consensus_pLann.html ~/Documents/Barcoding/barcode${P}/Results && cp consensus_pLann.gbk ~/Documents/Barcoding/barcode${P}/Results


for P in {81..90}; do cd ~/Documents/13042023/Barcoding/barcode${P} && cat *.fastq.gz>b${P}.fastq.gz; done
for P in {81..90}; do NanoPlot --fastq ~/Documents/13042023/Barcoding/barcode${P}/b${P}.fastq.gz -o ~/Documents/13042023/Barcoding/barcode${P}; done
flye --nano-hq ~/Documents/13042023/Barcoding/barcode83/b83.fastq.gz --genome-size 3k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode83/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode83/b83.fastq.gz -d ~/Documents/13042023/Barcoding/barcode83/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode83/ -m r1041_e82_260bps_hac_g632
for P in 84 85 89; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 6.5k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632; done
for P in 81 82 86; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 9k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632; done
for P in 87 88 90; do flye --nano-hq ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 11.5k --asm-coverage 30 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632; done
for P in {81..90}; do cd ~/Documents/13042023/Barcoding/barcode${P} &&
conda activate plannotate && 
plannotate batch -i consensus.fasta --html && conda deactivate && cd ~/ ; done 
for P in {81..90}; do cd ~/Documents/13042023/Barcoding/barcode${P} && mkdir Results && cp b${P}.fastq.gz Non_weightedHistogramReadlength.png consensus.fasta consensus_pLann.html consensus_pLann.gbk ~/Documents/13042023/Barcoding/barcode${P}/Results 

flye --nano-raw ~/Documents/13042023/Barcoding/barcode85/b85.fastq.gz --genome-size 6.5k --asm-coverage 20 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode85/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode85/b85.fastq.gz -d ~/Documents/13042023/Barcoding/barcode85/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode85/ -m r1041_e82_260bps_hac_g632
flye --nano-raw ~/Documents/13042023/Barcoding/barcode88/b88.fastq.gz --genome-size 11.8k --asm-coverage 10 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode88/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode88/b88.fastq.gz -d ~/Documents/13042023/Barcoding/barcode88/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode88/ -m r1041_e82_260bps_hac_g632
flye --nano-raw ~/Documents/13042023/Barcoding/barcode81/b81.fastq.gz --genome-size 9k --asm-coverage 20 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode81/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode81/b81.fastq.gz -d ~/Documents/13042023/Barcoding/barcode81/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode81/ -m r1041_e82_260bps_hac_g632





for P in 88 90; do flye --nano-raw ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz --genome-size 11.5k --asm-coverage 20 --iterations 5 -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632 && cd ~/Documents/13042023/Barcoding/barcode$P && conda activate plannotate && plannotate batch -i consensus.fasta --html && conda deactivate && cd ~/ && cd ~/Documents/13042023/Barcoding/barcode$P && mkdir Results && cp b$P.fastq.gz Non_weightedHistogramReadlength.png consensus.fasta consensus_pLann.html consensus_pLann.gbk ~/Documents/13042023/Barcoding/barcode$P/Results && cd ~/; done

for P in 84; do flye --nano-raw ~/Documents/13042023/Barcoding/barcode$P/fastq_runid_49ec5dea2a715fe3b188586c3891c7ce2d41174a_0.fastq.gz --genome-size 7k --asm-coverage 10 --iterations 10 --no-alt-contigs -o ~/Documents/13042023/Barcoding/barcode$P/ && 
medaka_consensus -i ~/Documents/13042023/Barcoding/barcode$P/b$P.fastq.gz -d ~/Documents/13042023/Barcoding/barcode$P/assembly.fasta -o ~/Documents/13042023/Barcoding/barcode$P/ -m r1041_e82_260bps_hac_g632 && cd ~/Documents/13042023/Barcoding/barcode$P && conda activate plannotate && plannotate batch -i consensus.fasta --html && conda deactivate && cd ~/ && cd ~/Documents/13042023/Barcoding/barcode$P && mkdir Results && cp b$P.fastq.gz Non_weightedHistogramReadlength.png consensus.fasta consensus_pLann.html consensus_pLann.gbk ~/Documents/13042023/Barcoding/barcode$P/Results && cd ~/; done
