# de_novo_and_reference_based_short_read_genome
No consensus as of yet has been reached for which genome assembler is 'the' best. As a matter of fact, there probably isn't, with different assemblers performing better for specific data sets or objectives. Here I perform de novo assemblies using multiple algorithms as well as reference based assembly of two Anolis distichus females, and compare their outputs.
The first three steps that I present below have been worked out by Richard Glor (/rglor) and Alana Alexander (/laninsky).

Step 1: Preliminary QC & Quality Trimming
======
Preliminary QC of your sequences can be completed by applying [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to your fastq sequence files. [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular tool that "aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines." Carefully inspect the `html` output from [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) because this is the main way that you are going to make sure that nothing is seriously wrong with your data before you delve into the time consuming series of analyses discussed below.
```
fastqc -k 6 *
```

If your sequences look OK after preliminary QC, its time to get your sequences ready for downstream analyses by trimming adaptors and eliminating low quality sequences and basecalls. Here we do this operation on the reads generated from the short insert library. We are going to do this by using the function `cutadapt` to (1) trim Illumina adapter sequences, (2) discard reads <75 bp in length and (3) perform gentle trimming of low quality basecalls. This process should take around 12 hours to complete for a raw sequence file containing around 500 million 100-150 bp reads.
```
#PBS -N cutadapt_short.sh
#PBS -l nodes=1:ppn=1:avx,mem=16000m,walltime=48:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Short_Insert
#PBS -j oe
#PBS -o cutadapterror_short

fastqc -k 6 /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R1.fastq
fastqc -k 6 /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R2.fastq
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 25 -o Short_trimmed_R1.fastq.gz -p Short_trimmed_R2.fastq.gz /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R1.fastq /scratch/a499a400/anolis/deliv_Glor061715_raw/raw_fastq/Anolis_Genome_R2.fastq > cutadapt.log
```
After trimming is complete, use `fastqc` on the resulting files to check that adapters have been trimmed and that the newly generated `fastq` files look good. 

Step 2: Assembly Free Estimates of Genome Coverage, Size and Heterozygosity
======
We will use an assembly-free k-mer counting method to estimate genome size, genome coverage and heterozygosity. Such k-mer counting approaches are widely used in genome studies. Although numerous applications for interpreting k-mer counts have been introduced (e.g., [GenomeScope](http://qb.cshl.edu/genomescope/), [estimate Genomesize.pl](http://josephryan.github.io/estimate_genome_size.pl/)), most rely on the rapid k-mer counting program [jellyfish](http://www.genome.umd.edu/jellyfish.html) to generate k-mer counts. Selecting a k-mer size for these analyses is non-trivial. Advice on selecting a k-mer size for various types of applications can be found in a few places, including the FAQ for the error-correction function [Quake](http://www.cbcb.umd.edu/software/quake/faq.html), section 2 of the supplemental material for the program [KAT](https://doi.org/10.1093/bioinformatics/btw663), and the online documentation for [Genomescope](https://github.com/schatzlab/genomescope/blob/master/README.md). One approach that is widely-cited in recent genome papers follows the authors of the panda genome paper in relying on analyses of 17-mers ([Ruiqiang et al. 2010](http://dx.doi.org/10.1038/nature08696)), but other studies suggest that large k-mer sizes are appropriate for large complex genomes. However, large k-mer sizes can lead to longer runs and more intense memory demands. Running Jellyfish shouldn't take more than 12 hours when conducted on our short insert library including 183M paired reads.

```
#PBS -N jellyfish_short
#PBS -q bigm
#PBS -l nodes=1:ppn=16:avx,mem=20000m,walltime=12:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -j oe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Jellyfish
#PBS -o jellyfish_short_error

work_dir=$(mktemp -d)
cat /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R1.fastq /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R2.fastq > $work_dir/Short_trimmed_merged.fastq
jellyfish count -m 17 -o fastq.counts -C $work_dir/Short_trimmed_merged.fastq -s 10000000000 -U 500 -t 16
rm $work_dir/Short_trimmed_merged.fastq
mv $work_dir/* /scratch/glor_lab/rich/distichus_genome/Jellyfish/
rm -rf $work_dir
```
We can then make a histogram from this data as follows.
```
jellyplot.pl fastq.counts > fastq.counts.histo
```
The resulting histogram file can then be uploaded and evaluated by [GenomeScope](http://qb.cshl.edu/genomescope/), an online tool capable of providing assembly free estimates of genome size, coverage, heterozygosity and other statistics. For more information on Genomescope, please consult this [PDF version of a poster summarizing the platform](http://schatzlab.cshl.edu/publications/posters/2016/2016.AGBT.GenomeScope.pdf) and the [Genomescope github](https://github.com/schatzlab/genomescope/blob/master/README.md). Briefly, Genomescope estimates coverage by plotting the histogram of k-mer frequencies and searching for peaks on this plot. For organisms some degree of heterozygosity, we expect two peaks, a peak with lower coverage corresponding with heterozygous k-mers and a peak with higher coverage corresponding with homozygous k-mers. The relative height of the two peaks can be used to derive an overall estimate of heterozygosity.


Step 3: Correct Sequencing Errors
======
Although Illumina is generally regarded as a relatively error-free sequencing method, your sequences are still expected to include errors, most of which will be substitution errors. One way to fix these types of sequencing erros is to use a k-mer spectrum error correction framework such as EULER or Quake ([Kelley et al. 2010](http://genomebiology.com/2010/11/11/R116)). The basic idea of these approaches involves identification of particularly low frequency k-mers that are likely the result of sequencing errors. We have already seen these likely erroneous k-mers as the large peak near the Y-axis in the plots of k-mer coverage generated by GenomeScope at the end of Step 2. The [publication introducing Quake](http://genomebiology.com/2010/11/11/R116) has a nice introduction to the underlying algorithms. More recently Trowel has been introduced and has the advantage of incorporating information on quality scores.

Step 3a: Quake
-----
We use [Quake](http://www.cbcb.umd.edu/software/quake/index.html) here. Quake uses Jellyfish for k-mer counting. Running Quake is fairly straightforward, with simple instructions available via the [program's online manual](http://www.cbcb.umd.edu/software/quake/manual.html). If you are uncertain of what k-mer size to use, the [Quake FAQ](http://www.cbcb.umd.edu/software/quake/faq.html) suggests that you can determine an appropriate k-mer based on an estimate of genome size, where `optimal k-mer = ln(200 * Genome size in bases)/ln(4)`. We use a k-mer of 19 for Anolis distichus because we expect that anole genome is somewhere in the 1.5 Gb range based on our assembly-free analyses in GenomeScan. The `-q 33` tag denotes the manner in which your sequence quality scores are recorded and is correct for anole data. If you need to check for you data, you can use the simple solution offered during an [online discussion](https://www.biostars.org/p/63225/). Quake can't unzip zipped sequence files on the fly. For the second step of the Quake process 
```
#PBS -N quake_short
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=512000,walltime=12:00:00,file=200gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Quake
#PBS -j oe
#PBS -o quake_short_error

work_dir=$(mktemp -d)
cat /scratch/a499a400/anolis/contam_filtered/decontam_short_1.fastqc /scratch/a499a400/anolis/contam_filtered/decontam_short_2.fastqc | count-qmers -k 27 -q 33 > $work_dir/quake_27_counts
cov_model.py quake_27_counts
#cat /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R1.fastq /scratch/glor_lab/rich/distichus_genome/Short_Insert/Short_trimmed_R2.fastq | correct -k 27 -c 2 -m quake_19_counts -p 24 > $work_dir/
mv $work_dir/* /scratch/glor_lab/rich/distichus_genome/Quake
rm -rf $work_dir
```
