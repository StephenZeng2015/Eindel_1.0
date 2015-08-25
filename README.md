#Eindel V 1.0 user manual
#####2015-08-01

A software use Split-read method to call indels on second generation sequencing data. Eindel can works on both exomic sequencing data and whole genome sequencing data. And paired-end reads are used in Eindel to validate candidate indels.


##Get Eindel Ready

###System Requirement:
***Linux, Unix(including Mac)***

Currently we only provide source file, and user can download it and run with python. We will provide standalone version in comming future.

###Prerequisite:

Please check and make sure following software and packages are installed on your computer/Server:

***Python (2.7.x)***: Eindel currently only support Python 2.7.x, Python can be downloaded from this website, https://www.python.org/download/releases/2.7/

	Prerequisite packages:
	1. numpy (version >= 1.9.1 http://www.numpy.org)
	2. Biopython (version >= 1.64, http://biopython.org/wiki/Main_Page)

***Samtools (version >=0.1.19)*** Details: http://samtools.sourceforge.net

***Bowtie (version >=1.0.0)*** Details: http://bowtie-bio.sourceforge.net/index.shtml

***Bowtie2 (version >=2-2.2.1)*** Details: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

###Build index for Bowtie and Bowtie2

Example command to build Bowtie reference index:
`bowtie-build Homo_sapiens_assembly19.fasta ref_seq/Homo_sapiens_assembly19_bowtie`

Example command to build Bowtie2 reference index:
`bowtie2-build Homo_sapiens_assembly19.fasta ref_seq/Homo_sapiens_assembly19_bowtie2`

###Get Eindel

Users can download either python source code or executable program:

1. Python Source Code (file can be downloaded here: LINK to Eindel file)

2. Executable file (Will be provided in comming future)

### Prepare Repetitive Regions files

Repetitive regions file is packed with python source code. 

###Prepare Configuration file

Eindel will call various softwares(bowtie and samtools) and read informations from varous files(for examples: reference genome fasta file)

Please fill path of Bowtie, Bowtie2, Bowtie2 reference, Bowtie reference, Genome reference fasta file, samtools
to configure files, 


Following is example in the configure file (Eindel_folder/Eindel.conf), input paths to:

1. Bowtie, for example:

		bowtie_path=/home/zengshuai/Bowtie/bowtie-1.0.1
	
2. Bowtie2, for example:

		bowtie2_path=/home/zengshuai/bowtie2/bowtie2-2.2.1
	
3. Bowtie reference, for example:

		BowRef_Index=/home/zengshuai/ReferenceGenome/Homo_sapiens_assembly19.bowtie
	
4. Bowtie2 reference, for example:

		Bow2Ref_Index=/home/zengshuai/ReferenceGenome/Homo_sapiens_assembly19.bowtie2
	
5. Genome reference fasta file, for example:

		GenomeRef_fasta=/home/zengshuai/ReferenceGenome/Homo_sapiens_assembly19.fasta
	
6. samtools, for example:

		bowtie_path=/home/zengshuai/Bowtie/bowtie-1.0.1
	
7. repetitive Region file.

		bowtie_path=/home/zengshuai/Bowtie/bowtie-1.0.1
	


##Run Eindel


Start Eindel from source file:

`Python Eindel.py	[options]	[PROCESS]	[FASTQ1 file]	[FASTQ2 file]	[OutputFolder&Prefix]`

	 |	Values | Detail
----------|----------|---------------------
 options	| 	(optional)	 |Maximum deletion size that Eindel can detect (unit: bp)
 PROCESS	|A / V / C			 |segment length (unit: bp)
 FASTQ1 file  | file path	 		 |path to fastq file of first paired read
 FASTQ2 file  | file path	 		 |path to fastq file of first paired read
 OutputFolder&Prefix  | path	 		 | path to results folder and prefix of results files
   

   

**for example:**

***Run All process by a command:***

`./Eindel.py A read1.fa read2.fa Outputfolder/sample1`

***Or run Candidate Detect and Candidate Validate process in two command:***

`./Eindel.py C read1.fa read2.fa Outputfolder/sample1`

`./Eindel.py V read1.fa read2.fa Outputfolder/sample1`

Both sets of commands can generate two results files:

`Outputfolder/sample1.csv` and `Outputfolder/sample1.vcf`

details about the results files are discussed below.

###Options:

####for candidate identifying step (C or A in command line)

	option |	Default | Value Detail
----------|----------|---------------------
 -MaxSize	| 1000		 |Maximum deletion size that Eindel can detect (unit: bp)
 -SegLen	|22			 |segment length (unit: bp)
 -Core	   |4	 		 |Max numbers of CPU-cores can be used (possible value: 1, 4)
 -MinSp	| 2	 		 |minimal number of reads that identify a candidate indel
 --IncRepReg	| NULL	 | Include Indels in repetitive Region. By default Eindel excludes indels in repetitive regions
 --KeepTmp|NULL| Keep temporal files, by default, Eindel only keeps result files
 
 
####for Candidate Indel Validation Step(V or A in command line)	
	option |	Default | Value Detail
----------|----------|---------------------
-numSRAlt	|2	|Minimal number of supporting reads on indel allele (alternative reference) required to validate an indel
-numPSAlt	|1	|Minimal number of paired supporting reads on alternative reference to validate an indel
-sprAlt	|10	|Minimal spread of supporting reads on alternative reference to validate an indel
- MulRatOri	|0.3|	Proportion of supporting reads on original reference that can map to other locations on genome
- MulRatAlt	|0.3|	Proportion of supporting reads on alternative reference that can map to other locations on genome
- MulRatAlt2	|0.9|	Proportion of supporting reads on alternative reference that can map to other alternative reference

###Example with Options

####Process Eindel using default parameters Example:

Process Candidate Indel Detecting step using default setting:
`python Eindel.py C fastq/Patient.read1.fastq fastq/Patient.read2.fastq Eindel_result/Patient`

Candidate Indel Validation Step using default setting:
`python Eindel.py V fastq/Patient.read1.fastq fastq/Patient.read2.fastq Eindel_result/Patient`

Process both Candidate Indel Detecting step and Candidate Indel Validation step in single command:
`python Eindel.py A fastq/Patient.read1.fastq fastq/Patient.read2.fastq Eindel_result/Patient`

####Process Eindel using Specific parameters Example:

Run Eindel using only 1 CPU core:
`python Eindel.py -Core 1 A fastq/Patient2.read1.fastq fastq/Patient2.read2.fastq Eindel_result/Patient2`

Include Indels in repetitive region:
`python Eindel.py –-IncRepReg A fastq/Patient.read1.fastq fastq/Patient.read2.fastq Eindel_result/Patient_all_region`

Validate candidate indels using user-customed setting:
`python Eindel.py -numSRAlt 1 -numPSAlt 0 - MulRatAlt2 0.6 –V fastq/Patient.read1.fastq fastq/Patient.read2.fastq Eindel_result/Patient`


##Result

###Path to Result
Eindel provide a standard vcf file, and to easy of data analysis we also provide a txt file. Results in both files are identical. Path to result file is:

	[OutputPath&Prefix].txt
	
	[OutputPath&Prefix].vcf


for example:

	Outputfolder/sample1/test.vcf
	
	Outputfolder/sample1/test.txt


###Contents and Details of Result (txt file)
 
Column id| Column name| Detail
---------|------------|---------
1	| CHROM 		|Chromosome ID
2	|POS  		|Location on chromosome
3	|COMT 		| comments
4	|TYPE 		|indel type ( I = insertion, D = deletion )
5	|SIZE 		|indel size (bps)
6	|SPO		|number of supporting reads on original reference
7	|PSO		|number of paired supporting reads on original reference
8	|MGO	|number of supporting reads that can mapped to other genome location
9	|SprO	|spread of supporting reads on original reference
10	|SPA		|number of supporting reads of alternative reference
11	|PSA|number of paired supporting reads on original reference
12	|MGA	|number of supporting reads that can mapped to other genomic location
13	|MAA	|number of supporting reads that can mapped to other alternative reference
14	|SprA	|spread of supporting reads on alternative reference
15	|GT	|genotype ( 0 = not an indel, 1 = heterozygous indel, 2 = homozygous indel )
16	|P1		|probability of genotype 0
17	|P2		|probability of genotype 1
18	|P3		|probability of genotype 2

###Contents and Details of Result (vcf file)
 
Column id| Column name| Detail
---------|------------|---------
1	| CHROM 		|Chromosome ID
2	|POS  		|Location on chromosome
3	|ID			|Indel ID
4	|REF		|Reference allele
5	|ALT		|Alternative allele
6	|FILTER	|Indel call results, Pass or Fail
7	|TYPE 		|indel type ( I = insertion, D = deletion )
8	|SIZE 		|indel size (bps)
9	|SPO		|number of supporting reads on original reference
10	|PSO		|number of paired supporting reads on original reference
11	|MGO	|number of supporting reads that can mapped to other genome location
12	|SprO	|spread of supporting reads on original reference
13	|SPA		|number of supporting reads of alternative reference
14	|PSA|number of paired supporting reads on original reference
15	|MGA	|number of supporting reads that can mapped to other genomic location
16	|MAA	|number of supporting reads that can mapped to other alternative reference
17	|SprA	|spread of supporting reads on alternative reference
18	|GT|genotype ( 0 = not an indel, 1 = heterozygous indel, 2 = homozygous indel )
19	|P1		|probability of genotype 0
20	|P2		|probability of genotype 1
21	|P3		|probability of genotype 2

##citation:



##frequency asked question:


####2. Can Eindel works on Whole genome data?

Yes, Eindel can works on second generation sequencing data including both exomic sequecing data and whole genome sequencing data. For whole genome data, we recommand use exclude indels in repetitive regions to reduce comuptation hours, which is a default setting of Eindel. 

####3. Can Eindel works on fastq file containing reads with differnt length?

Currently, Enidel can't works on fastq files with different read length.


##Problem Solving:

####1. screen shown following sentences when initial mapping:


	perl: warning: Setting locale failed.
	perl: warning: Please check that your locale settings:
		LANGUAGE = "en_HK:en",
		LC_ALL = (unset),
		LC_CTYPE = "UTF-8",
		LANG = "en_HK.UTF-8"
	    are supported and installed on your system.
	perl: warning: Falling back to the standard locale ("C").
	perl: warning: Setting locale failed.
	perl: warning: Please check that your locale settings:
		LANGUAGE = "en_HK:en",
		LC_ALL = (unset),
		LC_CTYPE = "UTF-8",
		LANG = "en_HK.UTF-8"
	    are supported and installed on your system.
	perl: warning: Falling back to the standard locale ("C").


This is warning output from bowtie2 software. Eindel can continue work on the data.


