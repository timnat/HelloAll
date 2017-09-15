DifCover (or DiRC – Difference in Read Coverage ?)

Description
The DifCover pipeline aims to identify regions in a reference genome for which the read coverage of a sample1 to the reference is significantly different from the read coverage of a sample2. “Significantly different” is determined by user defined threshold on a ration between average coverages of given samples. The pipeline allows to exclude from a consideration the under-represented fragments (with low coverage in sequencing of both samples) and/or the regions that carry repetitive sequences. Both cases can be misleading in the coverage analysis. The DifCover pipeline is specifically oriented to the analysis of large genomes and can handle very fragmented assemblies. 

Method
The alignment of short reads to a reference genome can be characterized by the depth of coverage computed for each genomic position as number of reads mapped over it. The large number of bases in a reference genome and “natural” fluctuations of coverage prevent or make unreasonable considering coverage for each individual base. Instead, an average coverage is computed over some intervals. Locations and sizes of the intervals can be defined differently (in a various ways). Traditional tools already offer solutions that allow to compute average coverage over intervals (windows) of a given fixed size (sambamba) or to split contig into separate intervals consisting only of bases with the same coverage (bedtools). However these tools are not suited on practice for the analysis of large complex genomes with gaps and repeats. DifCover addresses mentioned issues by introducing the notion of “stretching” windows. Essentially each of genomic scaffolds is scanned sequentially to form windows of variable size, but with predefined number of bases that have coverage within user defined limitations. These stretching windows allow to bridge over under- and over-represented fragments allowing more specific analysis (GIVE BETTER MOTIVATION HERE). For each window (interval) an average coverage is reported and compared to the coverage of another sample. If there are still too many windows (intervals) for downstream analysis, windows within a scaffold with similar coverage ratio can be combined to larger continues regions. Finally regions with significant difference in coverages can be filtered for downstream analyses.
Maybe picture of streching windows with different coverage

=====================================================

USAGE
Prerequisites (MUST be in your PATH)
	BEDTOOLS
	SAMTOOLS 
	AWK
	DNAcopy (for R) // optional (https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)

Quick start
The DifCover pipeline includes several bash scripts and one C++ program. They can be run separately stage by stage, to experiment with parameters, or run in a bulk from run_difcover.sh with predefined in it parameters. This section gives an example on how to run entire pipeline.

INPUT: two coordinate sorted BAM files presenting short read alignments from two samples to the same reference
OUTPUT: *.DNAcopyout.upp file with regions of significant coverage difference (p-fragments).  Format details can be found in the next section.

Download DifCover
Copy DifCover/dif_cover_scripts/run_difcover.sh to the directory with BAM files and replace parameters with your values

	FOLDER_PATH='path to dif_cover_scripts directory '
	BAM1='path to sample1.bam'
	BAM2='path to sample2.bam'
	a=10		# minimum coverage for sample1
	A=219		# maximum coverage for sample1
	b=10		# minimum coverage for sample2
	B=240		# maximum coverage for sample2
	v=1000	# target number of valid bases in streching windows
	l=500		# minimum size of window to output
	AC=1.095	# Adjustment Coefficient (set AC to 1, if modal coverage is equal) 
	p=2		# enrichment scores threshold (when p=2 will report regions with coverage in sample1 being roughly 4 times larger than coverage in sample2)
	bin=1		# for an auxiliary analitical stage (5); generates enrichment scores histogram with scores in bins with floating precision 1. For more detailed histogram use 10, 100.

Run entire pipeline
./run_difcover.sh

=====================================================
Pipeline overview and stage by stage usage example
The DifCover pipeline includes several bash scripts and one C++ program. They can be run separately stage by stage, to experiment with parameters, or run in a bulk from run_difcover.sh with predefined in it parameters.

INPUT: coordinate sorted bam files for two samples and mandatory parameters (explained for each stage below) 
OUTPUT:  *.DNAcopyout.upp file with regions of significant coverage difference (p-fragments)
                   Intermediate files (explained for each stage below)

 <<  sample1.bam, sample2.bam, a, A, b, B, v, l, AC, p >>
		|
		|   
   (1)  from_bams_to_unionbed.sh  (sample1.bam, sample2.bam)
		|
	    	|
   (2)  from_unionbed_to_ratio_per_window (a, A, b, B, v, l)
		|
	    	|
   (3)  from_ratio_per_window__to__DNAcopy_output.sh (AC)
		|
		|
   (4)  from_DNAcopyout_to_p_fragments.sh (p)
		|
		|
       << p-fragments >>

Stage by stage usage example
#prepear input data
cd DifCover
cp ./dif_cover_scripts/run_difcover.sh test_data/
cd test_data/

Open run_difcover.sh in text editor. Set FOLDER_PATH to a path to the directory  dif_cover_scripts/
FOLDER_PATH=../dif_cover_scripts

## run stage (1)
$FOLDER_PATH/from_bams_to_unionbed.sh sample1.bam sample2.bam
	OUPUT: sample1_sample2.unionbedcv	//
            	   ref.length.Vk1s_sorted //keep it for following stages
	NOTES: This script calls different functions from BEDTOOLS. File 	sample1_sample2.unionbedcv stores coverage information from both samples, allowing 	coverage comparisons between them.

## run stage (2)
$FOLDER_PATH/from_unionbed_to_ratio_per_window_fake0 sample1_sample2.unionbedcv 10 219 10 240 1000 500
	a=10 		minimum coverage for sample1
	A=219		maximum coverage for sample1
	b=10		minimum coverage for sample2
	B=240		maximum coverage for sample2
	v=1000 	target number of valid bases in the window
	 l=500		minimum size of window to output (window includes valid and non valid bases)

	NOTES:
	1. The program will merge bed intervals to larger intervals (windows), until number of valid bases 	won't exceed v.
	2. Valid bases satisfy following conditions  
		1) C1 < A and C2 < B; 	    and 	2) C1 > a or C2 > b.
	3.  Each window has approximately v valid bases, but because window is formed from bed intervals it can have
	- fewer than v bases – in a case if the window hits the end of the scaffold
	- more than v bases – to avoid breaking of the last added bed interval
      	4. For each window the program computes
	Q1 – average coverage of valid bases across all merged bed intervals for sample1  
	Q2 – average coverage of valid bases across all merged bed intervals for sample2
	W1 – is sum of coverages of merged bed interval for sample1
	W2 – is sum of coverages of merged bed interval for sample2
	R = W1/W2, if W2>0
	R = W1/FAKE0, if W2=0 .
	If coverage of sample2 is zero for some window, the program treats it as a case when only half of a single read is aligned (fake zero). It allows 1) avoid division by 0, and 2) reflect coverage information for sample1. FAKE0 is a predefined constant.
	5. ** The program from_unionbed_to_ratio_per_window_fake0_v2 calculates R differently: R = Q1/Q2 or Q1/FAKE0, if Q2=0.

OUTPUT: sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500
Columns are: scaffold, window_start, size_of_window, number_of_valid_bases_in_window, Q1, Q2, R

## run stage (3)
$FOLDER_PATH/from_ratio_per_window__to__DNAcopy_output.sh sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500 1.095
	NOTES:
	1. AC = 1.095 is an Adjustment Coefficient that allows to take in an account initial difference in the amount of sampling or produced coverage for each sample. We recommend compute AC as (modal coverage of sample2) : (modal coverage of sample1) . 
	2. First, the enrichment score is calculated for each window log2[AC*R]. Second, DNAcopy program runs to merge windows with similar enrichment scores (see details in DNAcopy description) to larger intervals and calculates final score for each of them. 
	OUTPUT: 
sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500.log2adj_1.095
sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500.log2adj_1.095.pdf
sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500.log2adj_1.095.DNAcopyout
	In *.DNAcopyout columns are: scaffold, start position of first window in the interval, start position of last window in the interval, number of merged windows, averaged enrichment score

	NOTES:
	1. DNAcopy may take long run time to generate results for large genomes. To reduce time one can decrease number of windows by increasing size of windows (v) on stage (2).
	2. Alternatively, input file can be split and processed in parallel, just note that it may change results slightly, because DNAcopy first does statistical analyses of all given windows.
	3. If number of windows is not large or if most scaffolds are covered just by several windows, you can directly proceed to filtering or analysis of *.log2adj_1.095 (scaffold, window_start, enrichment_score) without running DNAcopy.
Fig.1. Plot produced by DNAcopy. Each dot corresponds to a window.Green dots represent windows from first scaffold, black – from the second, green from the third, and so on. Red lines show an average enrichment scores for the regions computed from enrichments scores of merged windows.

## run stage (4)
Filter only genomic regions with enrichment scores > p.
FOLDER_PATH/from_DNAcopyout_to_p_fragments.sh sample1_sample2.unionbedcv.ratio_per_w_a2_A48_b2_B54_v1000_l500.log2adj_1.166.DNAcopyout 2
	OUTPUT:
sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500.log2adj_1.095.DNAcopyout.up2
sample1_sample2.ratio_per_w_fake0_a10_A219_b10_B240_v1000_l500.log2adj_1.095.DNAcopyout.down-2
	
	NOTES:
	1. Generated *pdf file provides visualization for the distribution of enrichment scores and can be helpful in the choice of threshold p. 
	2. The script filters from file *.DNAcopyout to *.DNAcopyout.upp fragments with enrichment scores ≥ p, i.e. fragments with coverage in a sample1 higher than in a sample2, and to *.DNAcopyout.down-p fragments with enrichment scores ≤-p, i.e. fragments with coverage in a sample2 higher than in a sample1.

Method Details
The DifCover works by comparing average depth of coverage across continuous intervals containing approximately v valid bases. The valid bases are determined by user defined lower and upper limits on depth of coverage for sample1 and sample2, defined respectively by a, A for sample1, and b, B for sample2. Some base with coverage C1 and C2 is defined to be valid if 1) C1 < A and C2 < B; and also 2) C1 > a or C2 > b. Upper limits allow to determination and skipping of fragments that contain repeats, while lower limits serve to exclude underrepresented fragments – gaps and fragments with too small number of reads in both samples. For identification of coverage differences only in single-copy regions we recommend lower limits to be assigned to one third of modal coverage and upper limits to 3X of modal coverage. 
The recruitment of valid bases of a given scaffold to the stretching windows (intervals) is done by traversing the scaffold from the beginning to the end. After the number of valid bases in a window reaches v, the window is closed and analysed. If the end of the scaffold was reached before v valid bases recruited to a current window, the window is discharged if its size is less than l. This approach provides flexibility for balancing between coarse granularity for large scaffolds and possibility to incorporate in the analysis short scaffolds of l < size < v, valuable feature for highly fragmented assemblies. 

phylosophy (about use of coverage analyses) may be an introduction in paper
The alignment of short read to a reference genome can be characterized by the depth of coverage computed for each genomic position as number of reads aligned to it. For high quality sufficient amount of sequencing and reasonably good reference assemblies, reads alignment usually reveals a normal distribution of depth of coverage, meaning that most genomic fragments are sampled and sequenced evenly, while other fragments can be sequenced in a different rate due to sequencing technologies, bringing in significant fluctuations in a read coverage. Presents of gaps, misassembled regions, repeats, sequencing errors and misalignments may further affect coverage distribution. Alternatively the difference in coverage across genomic segments can be caused by some genomic features like abundance of some subsequences in the genome (repeats), differences between reference genome and sequenced genome ( alignment to the genome of close related species, presents of different alleles, structural variations, SNPs, high level of polymorphisms, sex specific sequences…). Analysing read coverage of alignments from different samples we can learn about both underlying genome and sequenced samples. 1. If there are fragments in reference genome with zero coverage, those fragments could be either not sequenced properly, or absent in the sample. 2. Fragments with extremely high coverage tell us that underlying reference fragments are carrying some repeats. 3. If over some long continues fragments coverage is roughly half of modal coverage, that may be a sign of an allelic fragment, that is present only on one chromosome of the pair. 
Coverage analysis can be an informative instrument in a comparative study of several genomes. Here we present a pipeline for ...

