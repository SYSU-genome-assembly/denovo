/*
 * ReadSimulator.h
 *
 *  Created on: Jun 21, 2012
 *      Author: Xiaohu Shen
 */
#include "util.h"
#include <vector>
#include <cmath>
using namespace std;

#ifndef READSIMULATOR_H_
#define READSIMULATOR_H_

vector<int> GenerateSeq(int length);

vector<int> GenerateSeqVariation(vector<int> seq_orig, double vary_rate);

vector<vector<int> > GenerateReads(vector<int> seq, int coverage, int read_len, vector<double> err_profile);

void GenerateFastqFile(char* seq_filename, int coverage, int read_len, vector<double> err_profile, char* fastqfilename);

vector<double> LoadErrProfile(char* err_filename, int read_len);

#endif /* READSIMULATOR_H_ */
