#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "util.h"
#include "Matrix.h"

std::map<char, int> ACGTMap = ACGTMapInitializer();
std::map<int, char> ACGTRevMap = ACGTRevMapInitializer();
GSLRandom Random;

std::map<char, int> ACGTMapInitializer(){
    std::map<char, int>* __map = new std::map<char, int>();
    std::map<char, int>& map = *__map;

    map['A'] = 0; map['a'] = 0;
    map['C'] = 1; map['c'] = 1;
    map['G'] = 2; map['g'] = 2;
    map['T'] = 3; map['t'] = 3;

    return map;
}

std::map<int, char> ACGTRevMapInitializer(){
    std::map<int, char>* __map = new std::map<int, char>();
    std::map<int, char>& map = *__map;

    map[-1] = 'N';
    map[0] = 'A';
    map[1] = 'C';
    map[2] = 'G';
    map[3] = 'T';
  
    return map;
}

    
std::vector<Vector> loadIntensityFile(const char* filename){
    std::vector<Vector> data;
    std::ifstream fin;
    
    fin.open(filename);
    double intensity[4];
    while( fin >> intensity[0] >> intensity[1] >> intensity[2] >> intensity[3]){
	Vector vec(4);
	vec.setValues(intensity);
	data.push_back(vec);
    }
    fin.close();

    return data;
}

std::vector<std::vector<Vector> > loadTrainingDataFile(const char* filename, int seqLen, int setSize){
	std::vector<std::vector<Vector> > data;
	std::vector<Vector> seq;
	double intensity[4];
	std::ifstream fin;
	fin.open(filename);
	for(int i=0; i<setSize; i++){
		seq.clear();
		for(int j=0; j<seqLen; j++){
			fin >> intensity[0] >> intensity[1] >> intensity[2] >> intensity[3];
			Vector vec(4);
			vec.setValues(intensity);
			seq.push_back(vec);
		}
		data.push_back(seq);
	}
	fin.close();

	return data;
}


void loadData(const char* filename, std::vector<std::pair<std::vector<int>, std::vector<Vector> > >& Data){
	std::ifstream fin;
    
	fin.open(filename);

	string buffer;
	while(getline(fin, buffer)){
		double intensity[4];
		vector<int> coor(4);
		std::vector<Vector> cluster;		
		istringstream fcluster;
		fcluster.str(buffer);

		fcluster >> coor[0] >> coor[1] >> coor[2] >> coor[3];
		while( fcluster >> intensity[0] >> intensity[1] >> intensity[2] >> intensity[3]){
			Vector vec(4);
			vec.setValues(intensity);
			cluster.push_back(vec);
		}
		Data.push_back(std::make_pair(coor, cluster));
	}
	fin.close();
}

void displayResult(const std::vector<std::pair<char, Vector> >& res, const char* ans_fname){
    std::ifstream fin;
    std::string MySeq, AnsSeq, SolexaSeq;

    for(std::vector<std::pair<char, Vector> >::const_iterator i=res.begin(); i<res.end(); i++)
	MySeq.append(1, i->first);
    
    fin.open(ans_fname);
    fin >> SolexaSeq >> AnsSeq;

    std::cout << "A: " << AnsSeq << std::endl;
    std::cout << "S: " << SolexaSeq << std::endl;
    std::cout << "P: " << MySeq << std::endl;

    std::cout << "E: ";
    for(unsigned int i=0; i<AnsSeq.size(); i++)
	if(AnsSeq[i]==MySeq[i] && AnsSeq[i]==SolexaSeq[i])
	    std::cout << ' ';
	else if(AnsSeq[i]==MySeq[i] && AnsSeq[i]!=SolexaSeq[i] )
	    std::cout << '*';
	else if(AnsSeq[i]!=MySeq[i] && AnsSeq[i]==SolexaSeq[i] )
	    std::cout << '|';
	else
	    std::cout << 'x';
    std::cout << std::endl;
    
    for(unsigned int i=0; i<res.size(); i++){
	for(int b=0; b<4; b++)
	    std::cout << res[i].second[b] << ' ';
	std::cout << std::endl;
    }    
}

double sum(const std::vector<double>& vec){
    double s = 0.0;
    for(std::vector<double>::const_iterator i=vec.begin(); i!=vec.end(); i++)
	s += *i;
    return s;
}
