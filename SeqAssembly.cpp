#include <iostream>
#include <fstream>
#include <string>
#include "ReadSimulator.h"
#include "Matrix.h"
#include <math.h>


using namespace std;

//bowtie alignment?
//simple alignment

void displayvec(vector<int> vec){
	int len = vec.size();
	for(int i=0;i<len;i++){
		cout<<vec[i]<<" ";
	}
	cout<<endl;
}

vector<int> ReadDataFile(const char* filename){
	vector<int> sequence;	
	ifstream data_file;
	data_file.open(filename);
	string firstline;
	char singlebase;
	int baseint;
	if (data_file.is_open())
  {
  	getline(data_file,firstline);
    while (!data_file.eof()) {
	data_file>>singlebase;
	switch(singlebase){
	case 'A':
		baseint = 0;
		break;
	case 'C':
		baseint = 1;
		break;
	case 'G':
		baseint = 2;
		break;
	case 'T':
		baseint = 3;
		break;	
	default:
		baseint = 4;					
	}
	sequence.push_back(baseint);
 }
  }
  else cout << "Unable to open file";
  data_file.close();
  for(int i=0;i<20;i++){
  	cout<<sequence[i];
  	}
  	cout<<endl;
  return sequence;
}
	

vector<int> SimpleAlignment(vector<int> refseq, vector<vector<int> > AllReads){
	//return the aligned positions
	int seq_len = refseq.size();
	int num_reads = AllReads.size();
	int read_len = AllReads[0].size();
	vector<int> aligned_pos;
	for(int i=0;i<num_reads;i++){
		vector<int> read = AllReads[i];
		int max_match = 0;
		int match_pos = 0;
		for(int j=0;j<=seq_len-read_len;j++){
			int match = 0;
			for(int k=0;k<read_len;k++){
				if(AllReads[i][k]==refseq[j+k]) match++;
			}
			if(match>max_match){
				max_match = match;
				match_pos = j;
			}
		}
		aligned_pos.push_back(match_pos);
	}
	return aligned_pos;
}

vector<vector<int> > MultipleAlignment(vector<int> refseq, vector<vector<int> > AllReads, vector<int> &aligned_pos){
	int seq_len = refseq.size();
	int num_reads = AllReads.size();
	int read_len = AllReads[0].size();
	
	//vector<int> aligned_pos;
	vector<vector<int> > duplicate_reads;
	vector<int> num_align;  //store the number of alignments for each read
	int * match_score = new int[seq_len-read_len+1];
	for(int i=0;i<num_reads;i++){
		vector<int> read = AllReads[i];
		int max_match = 0;
		int num_multi_align = 0;
		for(int j=0;j<=seq_len-read_len;j++){
			int match = 0;
			match_score[j] = 0;
			for(int k=0;k<read_len;k++){
				if(AllReads[i][k]==refseq[j+k]) match++;
			}
			if(match>max_match){
				max_match = match;
			}
			match_score[j] = match;
		}
		//filter out bad reads
		int match_limit = 5; //set up a maximum number of possible matches
		if(max_match > 0.7*read_len){
			int num_match = 0;
		for(int j=0;j<=seq_len-read_len;j++){			
			if(match_score[j]==max_match){
				if(num_match<match_limit){
				duplicate_reads.push_back(AllReads[i]);
				aligned_pos.push_back(j);
				num_multi_align++;
				num_match++;
				}
			}
		}
		num_align.push_back(num_multi_align);
		}else{
			num_align.push_back(0);
		}
	}
	
//	int num_duplicate_reads = aligned_pos.size();
//	int* positions = new int[num_duplicate_reads];
//	for(int i=0;i<num_duplicate_reads;i++){
//		positions[i] = aligned_pos[i];
//	}
//	align_pos_ptr = &positions;
	//plot alignment histogram
	cout<<"alignment histogram"<<endl;
	int hist[6] = {0,0,0,0,0,0};
	for(int i=0;i<num_reads;i++){
		if(num_align[i]<=5){
			hist[num_align[i]]++;
		}else{
			cout<<"read alignment "<<num_align[i]<<endl;
		}
	}
	for(int i=0;i<6;i++){
		cout<<i<<": "<<hist[i]<<endl;
	}
	
	delete [] match_score;
	
	//check
	cout<<"size of duplicate reads: "<<duplicate_reads.size()<<" size of aligned positions: "<<aligned_pos.size()<<endl;
	return duplicate_reads;
}

vector<int> SequenceAssembly(vector<int> refseq, vector<int> actualseq, vector<vector<int> > AllReads, vector<int> positions, int k_max){
	// create connection tables
	int seq_len = refseq.size();
	int num_reads = AllReads.size();
	int read_len = AllReads[0].size();
	int num_bad_seq = 0;
	vector<int> base_covered;
	//cout<<"read_len is "<<read_len<<endl;
	int * reverseTable = new int [num_reads*read_len];  //the map of connected bases for each read
	for(int i=0;i<num_reads;i++){
		for(int j=0;j<read_len;j++){
			reverseTable[i*read_len+j] = positions[i]+j;
		}
	}
	vector<int> temp;
	vector<vector<int> > connectTable(seq_len, temp);  //the map of related reads for each base
	for(int i=0;i<num_reads;i++){
		for(int j=0;j<read_len;j++){
			connectTable[positions[i]+j].push_back(i);
		}
	}
	for(int i=0;i<seq_len;i++){
		//cout<<connectTable[i].size()<<endl;
		//displayvec(connectTable[i]);
		if(connectTable[i].size()>0) base_covered.push_back(i);
		if(actualseq[i]>3) num_bad_seq++;
	}
	
	//calculate coverage histogram
	int coverage[25];
	for(int i=0;i<25;i++) coverage[i]=0;
	for(int i=0;i<seq_len;i++){
		int cov = connectTable[i].size();
		if(cov<25){
			coverage[cov]++;
		}
	}
	cout<<"coverage histgram"<<endl;
	for(int i=0;i<25;i++){
		cout<<i<<": "<<coverage[i]<<endl;
	}
	cout<<"connect table finished"<<endl;
	//initialize read messages y
	double *y = new double [num_reads*read_len];
	for(int i=0;i<num_reads*read_len;i++) y[i] = 0.8; //should be 0.985

	vector<vector<Vector> > x;
	Vector tmp = Vector::zeros(4);
	vector<Vector> vec_tmp;
	for(int i=0;i<seq_len;i++){
		vec_tmp.clear();
		int con_len = connectTable[i].size();
		for(int j=0;j<con_len;j++){
			vec_tmp.push_back(tmp);
		}
		x.push_back(vec_tmp);
	}
	//cout<<x[4][0]<<endl;
	cout<<"initialization finished"<<endl;
	
	double init[4] = {-1.0, -1.0, -1.0, -1.0};
	//iterations
	for(int k=0;k<k_max;k++){
		// update base messages x
		vector<int>::iterator it;
		for(it=base_covered.begin();it<base_covered.end();it++){
			int i = *it;
			//cout<<"i="<<i<<endl;
			int con_len = connectTable[i].size();
			//if(con_len==0) {
				//cout<<"0 length occured"<<endl;
				//set the message to zeros
				//for(int l=0;l<con_len;l++){
				//x[i][l] = Vector::zeros(4);					
					//}
				//continue;
			//}
			//cout<<con_len<<endl;
			for(int l=0;l<con_len;l++){
				//cout<<"l "<<l<<endl;
				//int j=connectTable[i][l];
				x[i][l] = Vector::zeros(4);
				
				
				for(int j_prim=0;j_prim<con_len;j_prim++){
					if(j_prim!=l){
						int read_ind = connectTable[i][j_prim];
						int start_pos = positions[read_ind];
						int relative_pos = i-start_pos;
						//cout<<"read_ind "<<read_ind<<" start_pos "<<start_pos<<" relative_pos "<<relative_pos<<endl;
						//Vector edge = Vector::zeros(4);
						Vector edge(4);
						edge.setValues(init);
						int base = AllReads[read_ind][relative_pos];
						//cout<<j_prim<<relative_pos<<endl;
						if(base<4){
							edge[base] = 1;
						}
						//cout<<"1"<<endl;
						x[i][l]=x[i][l]+y[read_ind*read_len+relative_pos]*edge;
					}
				}
				//double sum = x[i][l][0]+x[i][l][1]+x[i][l][2]+x[i][l][3];
				double norm = x[i][l][0]*x[i][l][0]+x[i][l][1]*x[i][l][1]+x[i][l][2]*x[i][l][2]+x[i][l][3]*x[i][l][3];
				norm = sqrt(norm);
				//cout<<"sum "<<sum<<endl;
				if(norm>0){
					x[i][l] = x[i][l]/norm;
					}
				//cout<<"x "<<x[i][l]<<endl;
				//normalization according to degree
				//x[i][l] = x[i][l]/(con_len-1);
			}
		}
		
		cout<<"x message updated"<<endl;
		//update read messages y
		for(int j=0;j<num_reads;j++){
			int start_pos = positions[j];
			for(int l=0;l<read_len;l++){
				y[j*read_len+l] = 0;
				for(int i_prim=0;i_prim<read_len;i_prim++){
					if(i_prim!=l){
						int x_ind = start_pos+i_prim;
						int con_len = connectTable[x_ind].size();
						int base = AllReads[j][i_prim];
						Vector edge = Vector::zeros(4);
						if(base<4){
							edge[base] = 1;
						}
						for(int con_ind=0;con_ind<con_len;con_ind++){
							if(connectTable[x_ind][con_ind]==j){
								y[j*read_len+l] = y[j*read_len+l] + edge.dot(x[x_ind][con_ind]);
							}
						}

					}
				}
				y[j*read_len+l] = y[j*read_len+l]/(read_len-1);
			}
		}
	}
	
	cout<<"iteration finished"<<endl;
	//Estimation
	vector<int> base_est;
	for(int i=0;i<seq_len;i++){
		Vector xp = Vector::zeros(4);
		int con_len = connectTable[i].size();
		if(con_len==0) {
				//cout<<"0 length occured"<<endl;
				//set the message to zeros
				for(int l=0;l<con_len;l++){
				x[i][l] = Vector::zeros(4);					
					}
				base_est.push_back(-1);
				continue;
			}
		for(int j=0;j<con_len;j++){
			int read_ind = connectTable[i][j];
			int start_pos = positions[read_ind];
			int relative_pos = i-start_pos;
			//Vector edge = Vector::zeros(4);
			Vector edge(4);
			edge.setValues(init);
			int base = AllReads[read_ind][relative_pos];
			if(base<4){
				edge[base] = 1;
			}
			xp=xp+y[read_ind*read_len+relative_pos]*edge;
		}
		int maxbase = 0;
		double p_max = xp[0];
		for(int base_ind = 1;base_ind<4;base_ind++){
			if(xp[base_ind]>xp[maxbase]){
				p_max = xp[base_ind];
				maxbase = base_ind;
			}
		}
		base_est.push_back(maxbase);
	}
	
	//calculate error rate, snp FP, FN
	int err = 0;
	int FP = 0;
	int FN = 0;
	for(vector<int>::iterator it=base_covered.begin();it<base_covered.end();it++){
		if(base_est[*it]!=actualseq[*it] && actualseq[*it]<4 ){
			err++;
			//False positive
			if(refseq[*it]==actualseq[*it]){
				FP++;
			}else if(refseq[*it]==base_est[*it]){
				FN++;
			}
		}
	}
	cout<<"Error number is "<<err<<". Error rate is "<<double(err)/(base_covered.size()-num_bad_seq)<<endl;
	cout<<"number of uncovered bases is "<<seq_len - base_covered.size()<<endl;
	cout<<"False positive is "<<FP<<endl;
	cout<<"False negative is "<<FN<<endl;
	delete [] y;
	delete [] reverseTable;
	return base_est;
}



int main(){
	vector<int> refseq;
	vector<int> actualseq;
	vector<vector<int> > AllReads;
	vector<int> align_positions;

	int read_len = 76;
	int coverage = 9;
	int k_max = 11;
	//int seq_len = 200;
	//int seq_len = 20000;
	//refseq = GenerateSeq(seq_len);
	//concatanate vectors
//	vector<int> tmp1 = GenerateSeq(seq_len);
//	refseq = tmp1;
//	refseq.insert(refseq.end(),tmp1.begin(),tmp1.end());
	//actualseq = refseq;
	
	//load sequence from file
	char* datafilename = "Seq Data/ecoli.fasta";	
	vector<int> total_refseq = ReadDataFile(datafilename);
	int choose_length = 100000;
	for(int i=0;i<choose_length;i++){
		refseq.push_back(total_refseq[i]);
	}
	int seq_len = refseq.size();
	cout<<"sequence length is"<<seq_len<<endl;
	//return 0;
	actualseq = GenerateSeqVariation(refseq, 0.01);
	//displayvec(refseq);
	//displayvec(actualseq);
	vector<double> err_profile;
	double init_err= 0.008;
	for(int i=0;i<read_len;i++){
		err_profile.push_back(init_err);
		init_err += 0.00015;
	}
	AllReads = GenerateReads(actualseq, coverage, read_len, err_profile);
	cout<<"Read size"<<AllReads.size()<<endl;
	//cout<<"Reads generated"<<endl;
	//align_positions = SimpleAlignment(refseq, AllReads);
	vector<vector<int> > duplicate_reads;
	duplicate_reads = MultipleAlignment(refseq, AllReads, align_positions);
	
	cout<<"Alignment finished"<<endl;
	//cout<<align_positions[0]<<endl;
	vector<int> est_bases;
	//est_bases = SequenceAssembly(refseq, actualseq, AllReads, align_positions, k_max);
	est_bases = SequenceAssembly(refseq, actualseq, duplicate_reads, align_positions, k_max);
	cout<<"Assembly finished"<<endl;

	return 0;
}
