/*
 * ReadSimulator.cpp
 *
 *  Created on: Jul 8, 2012
 *      Author: Xiaohu Shen
 */

#include "ReadSimulator.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


vector<int> GenerateSeq(int length){
	vector<int> seq;
	GSLRandom Random;
	int base;
	for(int i=0;i<length;i++){
		base = Random.discreteUniformRv(0,3);
		seq.push_back(base);
	}
	return seq;
}

vector<int> GenerateSeqVariation(vector<int> seq_orig, double vary_rate){
	int seq_len = seq_orig.size();
	vector<int> seq_vary = seq_orig;
	int length_vary = seq_len*vary_rate;
	int* ind_selected = new int[length_vary];
	GSLRandom Random;
	int num_selected = 0;
	double rand;
	for(int i=0;i<seq_len;i++){
		if(num_selected >= length_vary) {
			break;
		}
		else{
			rand = Random.uniformRv();
			double select_prob = (length_vary - num_selected)/(double)(seq_len-i);
			if(rand<select_prob){
				ind_selected[num_selected] = i;
				num_selected ++;
			}
		}
	}
//	for(int i=0;i<length_vary;i++){
//		cout<<ind_selected[i];
//	}
//	cout<<endl;
	for(int i=0;i<length_vary;i++){
		int base_orig = seq_orig[ind_selected[i]];
		if(base_orig > 3) continue;
		
		int bases_vary[3];
		int ind_tmp = 0;
		for(int j=0;j<4;j++){
			if(j!=base_orig){
				bases_vary[ind_tmp] = j;
				ind_tmp ++;
			}
		}
		int ind_rand = Random.discreteUniformRv(0,2);
		seq_vary[ind_selected[i]] = bases_vary[ind_rand];
	}

	delete [] ind_selected;
	return seq_vary;
}

vector<vector<int> > GenerateReads(vector<int> seq, int coverage, int read_len, vector<double> err_profile){
	int seq_len = seq.size();
	int num_reads = ceil(coverage*seq_len/read_len);
	vector<vector<int> > Allreads;
	vector<int> read;

	GSLRandom Random;
	for(int i=0;i<num_reads;i++){
		int st = Random.discreteUniformRv(0,seq_len-read_len);
		read.clear();
		for(int j=0;j<read_len;j++){
			if(seq[st+j]>3){
				read.push_back(4);
				continue;
			}
			double tmp_rnd = Random.uniformRv();
			int base;
			
			if(tmp_rnd>err_profile[j])
				base = seq[st+j]; //no error
			else{
				int baseset[4] = {0,1,2,3};
				baseset[seq[st+j]] = -1;
				int errorbases[3];
				int ind = 0;
				for(int k=0;k<4;k++){
					if(baseset[k]>=0){
						errorbases[ind] = baseset[k];
						ind++;
					}
				}
				int rnd = Random.discreteUniformRv(0,2);
				base = errorbases[rnd];
			}
			read.push_back(base);
		}
		Allreads.push_back(read);
	}
	return Allreads;
}

void GenerateFastqFile(char* seq_filename, int coverage, int read_len, vector<double> err_profile, char* fastqfilename){
	vector<char> sequence;	
	char base_set[4] = {'A','C','G','T'};
	ifstream data_file;
	data_file.open(seq_filename);
	string firstline;
	char singlebase;
	//int baseint;
	if (data_file.is_open()){
		getline(data_file,firstline);
		while (!data_file.eof()) {
			data_file>>singlebase;
			//cout<<signlebase;
			sequence.push_back(singlebase);
		}
	}
	else cout << "Unable to open file";
	data_file.close();

	int seq_len = sequence.size();
//	for(int i=0;i<seq_len;i++){
//		cout<<sequence[i];
//	}
//	cout<<endl;
	
	ofstream read_file;
	read_file.open(fastqfilename);
	
	int num_reads = coverage*seq_len/read_len;
	
	for(int i=0;i<num_reads;i++){
		//output 1 random read in fastq format
		read_file<<'@'<<endl;
		//cout<<'@'<<endl;
		int st = Random.discreteUniformRv(0,seq_len-read_len);
		for(int j=0;j<read_len;j++){
			double tmp_rnd = Random.uniformRv();
			double rand_rate = err_profile[j]*4/3;
			if(tmp_rnd>rand_rate){
				read_file<<sequence[st+j];
				//cout<<sequence[st+j];
			}
			else{
				read_file<<base_set[Random.discreteUniformRv(0,3)];
				//cout<<base_set[Random.discreteUniformRv(0,3)];
			}
		}
		read_file<<endl;
		//cout<<endl;
		read_file<<'+'<<endl;
		//cout<<'+'<<endl;
		for(int j=0;j<read_len;j++) {
			read_file<<'N';
			//cout<<'N';
		}
		read_file<<endl;
		//cout<<endl;
	}
	read_file.close();
}

vector<double> LoadErrProfile(char* err_filename, int read_len){
	ifstream err_file;
	err_file.open(err_filename);
	vector<double> err_profile;
	double err_rate;
	for(int i=0;i<read_len;i++){
		err_file>>err_rate;
		err_profile.push_back(err_rate);
	}
	return err_profile;
}