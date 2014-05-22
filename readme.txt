MPSequencing Readme v1.0
Xiaohu Shen  
xhshen@utexas.edu
1/2014

1 Installation
1.1 Requirements
ParticleCall should function on any standard Linux environment with required package
Gnu Scientific Library (GSL) (version >= 1.12)

1.2 Compiling instructions
From the directory containing the source files, simply type: 
> make

Two runnables are generated: FastqSimulate and SeqAssembly

2 Running instructions
After compilation, MPSequencing can be used to randomly generate short read DNA sequencing data file from a genome sequence file in FASTA format. The assemly sequence is obtained using alignment information of the reads to a reference sequence. 

ex. 
>./FastqSimulate

>./SeqAssembly


Output: FastqSimulate generates short reads data file with quality information in FASTQ format. SeqAssembly obtains the reference-guided assembled genome sequence from alignment information and outputs the number of mismatchs between assembled sequence and original sequence.

Note: In current implementation the input filenames are hard coded in the file.



