# Scripts used for analyzing PAIso-Seq data

### Software Dependencies
```
# softwares
ccs (v5.0.0)
bam2fastq (v1.3.0)
pbbam (v1.0.6)
pbmerge (v1.0.6)
minimap2 (2.17-r941)
bedtools (v2.29.2)
samtools (v1.9)
subread (V2.0.0)
featureCounts (v2.0.0)
fastq-splitter.pl (http://kirill-kryukov.com/study/tools/fastq-splitter)
ParaFly (http://parafly.sourceforge.net/)
```
```
# python packages
pysam (https://github.com/pysam-developers/pysam)
regex (https://bitbucket.org/mrabarnett/mrab-regex)
parasail (https://github.com/jeffdaily/parasail-python)
itertools (https://github.com/more-itertools/more-itertools)
numpy (https://numpy.org/)
statistics (https://pypi.python.org/pypi/statistics)
```

### Installation
All the scripts we provided do not need installation and can be run directly in standard python3 enviroment.

### Generate CCS reads from sub reads</br>
**Input:** subreads.bam from a single SMRTcell run in PacBio BAM format. "SCGV.subreads.bam" for example. </br>
**Output:** Consensus reads in unaligned BAM (.bam) format.</br>

```
# Run ccs command
ccs SCGV.subreads.bam SCGV.ccs.bam -j 30 &> SCGV.ccs.log 

# Convert bam format to fastq format
bam2fastq -o SCGV.ccs SCGV.ccs.bam > SCGV.ccs.bam2fastq.log

# Extract the pass number of each CCS reads
perl GetCCSpass.pl SCGV.ccs.bam > SCGV.ccs.pass.txt
```

### Process of PAIso-seq data

#### Step1. Extract the clean reads containing the poly(A) inclusive full-length cDNA from the CCS reads based on the barcodes for each sample 
```
gunzip SCGV.ccs.fastq.gz
fastq-splitter.pl --n-parts 100 SCGV.ccs.fastq
for i in {001..100}; do echo 'python CCS_split_clean_end_extension_v1.py SCGV.ccs.part-'${i}'.fastq SCGV.ccs.pass.txt barcode.fa 2 1> SCGV. '${i}'.out.txt 2> SCGV.${i}.err.txt'; done > split_clean.sh
ParaFly -c split_clean.sh -CPU 25
cat *.out.txt > SCGV.all.out.txt
cat *.err.txt > SCGV.all.err.txt
awk '{print "@"$2"\n"$6"\n+\n"$7}' SCGV.all.out.txt | gzip -nc > SCGV.clean.fastq.gz
rm SCGV.*.err.txt SCGV.*.out.txt *fastq.log.txt SCGV.ccs.part-*.fastq
```
#### Step2. Align the clean-reads to the genome
```
paftools.js gff2bed ../genome/gencode.vM25.primary_assembly.annotation.gtf > gencode.vM25.primary_assembly.annotation.bed
minimap2 -x splice -t 20 -d gencode.vM25.mmi ../genome/GRCm38.primary_assembly.genome.fa.gz
minimap2 -ax splice -uf --secondary=no -t 40 -L --MD --cs --junc-bed gencode.vM25.primary_assembly.annotation.bed gencode.vM25.mmi SCGV.clean.fastq.gz 2> align.log | samtools view -F 3844 -bS > SCGV.clean.filter.bam
```
#### Step3. Extract and annotate poly(A) tails from mapped clean-reads
```
python PolyA_trim_V5.4.1.py SCGV.clean.filter.bam > SCGV.polyA_trim.out.txt
featureCounts -L -g gene_id -t exon -s 1 -R CORE -a gencode.vM25.primary_assembly.annotation.gtf -o SCGV.featureCounts SCGV.clean.filter.bam &> featureCounts.log
python PolyA_note_V2.1.py SCGV.polyA_trim.out.txt SCGV.clean.filter.bam.featureCounts 1> SCGV.polyA_note.txt 2> SCGV.polyA_note.err.txt
```
#### (Optional). Extract clean poly(A) sequences of poly(A) spike-ins from the CCS reads
```
python DNA_spikein_extract_2019_NC_V1.3.py SCGV.ccs.fastq SCGV.ccs.pass.txt PSI-barcode.fa 2 1> PSI.out.txt 2> PSI.err.txt
```
### 
#### Output file</br>
Now you have the SCGV.polyA_note.txt file containing the essential information for each of the sequenced poly(A) tail contains the following 13 columns of information for one read: contains the following 13 columns of information for one read: barcode, read_id, ensembl_id, pass_number, “1”, number_of_residue_A, number_of_residue_T, number_of_residue_C, number_of_residue_G, number_of_residue_T+C+G, “0”, poly(A)-tail_sequence and average_quality_value_of_poly(A)-tail_bases.
### 
