#include "ReadSimulator.h"

int main(){
	vector<double> err_profile;
	char* err_filename = "/home/local/calvine/xhshen/BPSeqAssembly/NMerr";
	int read_len = 100;
	int coverage = 30;
	
	char* seq_filename = "/home/local/calvine/xhshen/BPSeqAssembly/NMeningitidis.fasta";
	char* fastqfilename = "/home/local/calvine/xhshen/BPSeqAssembly/NMReads.fastq";
	
	
	err_profile = LoadErrProfile(err_filename, read_len);
	GenerateFastqFile(seq_filename, coverage, read_len, err_profile, fastqfilename);
	
	
	return 0;
}