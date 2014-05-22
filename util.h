#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Matrix.h"

#ifndef __UTIL_H__
#define __UTIL_H__

#define PI 3.141592653589793


std::map<char, int> ACGTMapInitializer();
std::map<int, char> ACGTRevMapInitializer();
Vector* ACGT_XInitializer();


std::vector<Vector> loadIntensityFile(const char* filename);
std::vector<std::vector<Vector> > loadTrainingDataFile(const char* filename, int seqLen, int setSize);

void loadData(const char* filename, std::vector<std::pair<std::vector<int>, std::vector<Vector> > >& Data);
void displayResult(const std::vector<std::pair<char, Vector> >& res, const char* ans_fname);
double sum(const std::vector<double> & vec);

extern std::map<char, int> ACGTMap;
extern std::map<int, char> ACGTRevMap;



class GSLRandom{
    gsl_rng* r;
    const gsl_rng_type* T;    
 public:
    GSLRandom(){
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
    }

    inline double uniformRv(){
	return gsl_rng_uniform(r);
    }

    inline int discreteUniformRv(int st, int ed){
	return gsl_rng_uniform_int(r, ed-st+1) + st;
    }

    inline int MultinomialRv(int n, const double* prob){
	double cumulative_prob[n];

	cumulative_prob[0] = prob[0];
	for(int i=1; i<n; i++)
	    cumulative_prob[i] = cumulative_prob[i-1] + prob[i];
	double u = uniformRv();
	int selected = std::lower_bound(cumulative_prob, cumulative_prob+n, u) - cumulative_prob;

	return (selected>=n?n-1:selected);
    }

    inline double normalRv(){
	return gsl_ran_gaussian(r, 1.0);
    }

    inline double normalRv(double mu, double sigma) {
	return ( mu + sigma * normalRv() );
    }

    inline double lnMVariateNormalPdf(const Vector& y, const Vector& mu, const Matrix& sigma, const Matrix& inv_sigma, double det_sigma){
	double logProb = 0.0;
	double N = 1.0*mu.len();
    
	logProb += -0.5*(y-mu).dot(inv_sigma).dot(y-mu) - (N/2.0)*log(2.0*PI) - 0.5*log(det_sigma);
	return logProb;
    }
    
    inline double lnNormalPdf(double x, double mu, double sigma) {
	return -0.5 * std::log(2.0 * PI) - std::log(sigma) - 0.5 * (x - mu) * (x - mu) / (sigma * sigma);
    }
};

extern GSLRandom Random;

#endif


