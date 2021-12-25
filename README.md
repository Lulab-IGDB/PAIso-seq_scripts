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
**Input:** subreads.bam from a single SMRTcell run in PacBio BAM format. "test.subreads.bam" for example. </br>
**Output:** Consensus reads in unaligned BAM (.bam) format.</br>

```
# Run ccs command
ccs test.subreads.bam test.ccs.bam -j 30 &> test.ccs.log 

# Convert bam format to fastq format
bam2fastq -o test.ccs test.ccs.bam > test.ccs.bam2fastq.log

# Extract the pass number of each CCS reads
perl GetCCSpass.pl test.ccs.bam > test.ccs.pass.txt
```

### Process of PAIso-seq data

#### Step1. Extract the clean reads containing the poly(A) inclusive full-length cDNA from the CCS reads based on the barcodes for each sample 
```
gunzip test.ccs.fastq.gz
fastq-splitter.pl --n-parts 100 test.ccs.fastq
for i in {001..100}; do echo 'python CCS_split_clean_end_extension_v1.py test.ccs.part-'${i}'.fastq test.ccs.pass.txt barcode.fa 2 1> test.'${i}'.out.txt 2> test.'${i}'.err.txt'; done > split_clean.sh
ParaFly -c split_clean.sh -CPU 25
cat *.out.txt > test.all.out.txt
cat *.err.txt > test.all.err.txt
awk '{print "@"$2"\n"$6"\n+\n"$7}' test.all.out.txt | gzip -nc > test.clean.fastq.gz
rm test.*.err.txt test.*.out.txt *fastq.log.txt test.ccs.part-*.fastq
```
Note: The parameter 2 in 'python CCS_split_clean_end_extension_v1.py test.ccs.part-*.fastq test.ccs.pass.txt barcode.fa 2 ' means number of mismatches allowed for matching barcodes in barcode.fa file with reads, which we recommend no more than 4.</br>
```
# demo for "CCS_split_clean_end_extension_v1.py" with the demo data in folder /demo_data
python CCS_split_clean_end_extension_v1.py test.fastq test.pass.txt test.barcode.fa 2 1> test.out.txt 2> test.err.txt
awk '{print "@"$2"\n"$6"\n+\n"$7}' test.out.txt | gzip -nc > test.clean.fastq.gz
## run time: 93s
## expected output: 6936 rows in test.out.txt and 221 rows in test.err.txt
```

#### Step2. Align the clean-reads to the genome
```
paftools.js gff2bed ../genome/gencode.vM25.primary_assembly.annotation.gtf > gencode.vM25.primary_assembly.annotation.bed
minimap2 -x splice -t 20 -d gencode.vM25.mmi ../genome/GRCm38.primary_assembly.genome.fa.gz
minimap2 -ax splice -uf --secondary=no -t 40 -L --MD --cs --junc-bed gencode.vM25.primary_assembly.annotation.bed gencode.vM25.mmi test.clean.fastq.gz 2> align.log | samtools view -F 3844 -bS > test.clean.filter.bam
```

#### Step3. Extract and annotate poly(A) tails from mapped clean-reads
```
python PolyA_trim_V5.4.1.py test.clean.filter.bam > test.polyA_trim.out.txt
featureCounts -L -g gene_id -t exon -s 1 -R CORE -a ../genome/gencode.vM25.primary_assembly.annotation.gtf -o test.featureCounts test.clean.filter.bam &> featureCounts.log
python PolyA_note_V2.1.py test.polyA_trim.out.txt test.clean.filter.bam.featureCounts 1> test.polyA_note.txt 2> test.polyA_note.err.txt
```
```
# demo for "PolyA_trim_V5.4.1.py" with the demo data in folder /demo_data
minimap2 -ax splice -uf --secondary=no -t 40 -L --MD --cs --junc-bed gencode.vM25.primary_assembly.annotation.bed gencode.vM25.mmi test.clean.fastq.gz 2> align.log | samtools view -F 3844 -bS > test.clean.filter.bam
python PolyA_trim_V5.4.1.py test.clean.filter.bam > test.polyA_trim.out.txt
## run time: 20s
## expected output: 5285 rows in test.polyA_trim.out.txt
```
```
# demo for "PolyA_note_V2.1.py" with the demo data in folder /demo_data
featureCounts -L -g gene_id -t exon -s 1 -R CORE -a ../genome/gencode.vM25.primary_assembly.annotation.gtf -o test.featureCounts test.clean.filter.bam &> test.featureCounts.log
python PolyA_note_V2.1.py test.polyA_trim.out.txt test.clean.filter.bam.featureCounts 1> test.polyA_note.txt 2> test.polyA_note.err.txt
## run time: 5s
## expected output: 5232 rows in test.polyA_note.txt and 0 row in test.polyA_note.err.txt
```

#### (Optional). Extract clean poly(A) sequences of poly(A) spike-ins from the CCS reads
```
python DNA_spikein_extract_2019_NC_V1.3.py test_spikein.ccs.fastq test_spikein.ccs.pass.txt PSI-barcode.fa 2 1 1> test.PSI.out.txt 2> test.PSI.err.txt
```
Note: The parameter 2 in 'python CCS_split_clean_end_extension_v1.py test.ccs.part-*.fastq test.ccs.pass.txt barcode.fa 2 ' means number of mismatches allowed for matching barcodes in barcode.fa file with reads, which we recommend no more than 4.</br>
```
# demo for "DNA_spikein_extract_2019_NC_V1.3.py" with the demo data in folder /demo_data
python DNA_spikein_extract_2019_NC_V1.3.py test_spikein.fastq test_spikein.pass.txt PSI-barcode.fa 2 1 1> test.PSI.out.txt 2> test.PSI.err.txt
## run time: 1s
## expected output: 14 rows in test.PSI.out.txt and 2 rows in test.PSI.err.txt
```

### 
#### Output file</br>
Now you have the test.polyA_note.txt file containing the essential information for each of the sequenced poly(A) tail contains the following 13 columns of information for one read: contains the following 13 columns of information for one read: barcode, read_id, ensembl_id, pass_number, “1”, number_of_residue_A, number_of_residue_T, number_of_residue_C, number_of_residue_G, number_of_residue_T+C+G, “0”, poly(A)_tail_sequence and average_quality_value_of_poly(A)_tail_bases.
### 
