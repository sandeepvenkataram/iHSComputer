#ifndef __iHSComputer_H__
#define __iHSComputer_H__

#include <iostream>
#include <iterator>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>
#include <map>
#include <algorithm>

using namespace std;

struct EHHReturn{
	double EHH;
	vector<string> haploGroups;
	vector<string> growingHaps;
};

void readParameters(char * filename);
void readSNPFile(char * filename);
void readRecombFile(char * filename);
void computeStandardizediHS(void);
void createMergedFile(void);
EHHReturn computeEHH(vector<string> strainSets, int position, vector<string> growingHaplos, int numTotalStrains);
void computeiHS(int threadNum, ofstream &outFileStream);
int computeHalfiHS(map<int, double> * EHHd,map<int, double> * EHHa, vector<string> present, vector<string> absent, int position, int indexIncrement,int numPresent, int numAbsent);
double nC2(int n);
void threadFunction(int threadNum, int startingPosition);
void printHelpStatement(void);
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
string join(vector<string> &strings, string delimiter);
bool fexists(const char *filename);

// all of the global variables are used during multithreading. MAKE SURE YOU DO NOT WRITE TO THE VARIABLES WITHIN THE THREADS OR YOU WILL GET RACE CONDITIONS!
extern double minEHH;
extern double minMinorAlleleFrequency;
extern double maxMissingData;
extern int minNumHaplotypes; //a hidden parameter to ensure that all SNPs have at least 4 haplotypes to compute iHS for.
extern char singleIndexType; //'G' signifies no single index, 'I' signifies position in haplotype, 'P' signifies physical position in genome
extern double singleIndexPosition;
extern bool randomize;
extern int gapCorrection;
extern bool unequalIntegration;
extern bool thresholdSubtraction;
extern int numIterations;
extern int numHaplos;
extern int haploLength;
extern string chromosome;
extern string outdir;
extern vector<string> haploArray; 
extern map<double, int> physPosToIndex;
extern vector<double> positions;
extern vector<double> physicalPositions;
extern vector<int> recombDistances;
extern vector<double> recombRates;
extern vector<int> singlePositions;

#endif