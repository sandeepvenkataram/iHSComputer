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
/*
 * Sandeep Venkataram and Yuan Zhu
 * Petrov Lab
 * Stanford University
 * svenkat1@stanford.edu
*/


using namespace std;

void computeStandardizediHS(void);
double nC2(int n);
void printHelpStatement(void);
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
string join(vector<string> &strings, string delimiter);

vector<vector<double> > freqBins;
vector<double> freqBinMeans;
vector<double> freqBinStdevs;
vector<int> singlePositions;
int numBins;
bool normalizeEqualSizedBins = false;


/*
 * Split/Join code taken from Evan Terna's code at http://stackoverflow.com/questions/236129/splitting-a-string-in-c
 * 
 */
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

string join(vector<string> &strings, string delimiter){
	string result = "";
	if(strings.size()==0){
		return result;
	}
	result=strings[0];
	for(int i=1; i<strings.size(); i++){
		result = result + delimiter + strings[i];
	}
	return result;
}

// 1 based coordinates for snp indexing
// fix rounding on genomic position output
void printHelpStatement(void){
	cout<<"This program standardizes the iHS values computed by iHS and iHS_multithreaded"<<endl;
	cout<<"arguments: fileListFile numBins equalSizedBins[optional]"<<endl;
	cout<<"\n\nThe fileListFile has the following tab-delimited format per line: \"path_to_unstandardized_iHS_file [tab] name_of_chromosome_in_file\""<<endl;
	cout<<"Standardization is done by first binning iHS values by the minor allele frequency of the core SNP, and then conducting Z-normalization within each bin. The number of bins parameter thus defines the width of each frequency bin for normalization. If the equalSizedBins parameter is set to \"true\", bins are defined to contain an equal number of values rather than an equal minor allele frequency range, which may be more useful in detecting significant values at positions with a large minor allele frequency (which tends to be rarer). More bins means fewer values in each bin, so extreme values may become less noticeable if there is insufficient values in each bin."<<endl;
}

/*
 * Computes a z-score standardized iHS value binned by derived allele frequency (subtract the mean, then divide by the standard deviation)
 */
//maybe use equal sized bins???
void computeStandardizediHS(vector<string> filenames, vector<string> chromosomes){
	
	vector<vector<double> > temp (numBins);
	freqBins=temp;
	
	if(normalizeEqualSizedBins){
		vector<double> iHSValues;
		for(int i=0; i<filenames.size(); i++){
			ifstream inFile(filenames[i].c_str());
			while(!inFile.eof()){
				string l;
				getline(inFile,l);
				if(l=="" || l.find("SNP") != std::string::npos){
					continue;
				}
				vector<string> sp;
				sp = split(l, '\t');
				double iHS = atof(sp[5].c_str());
				iHSValues.push_back(iHS);
				
			}
			inFile.close();
		}
		
		std::sort (iHSValues.begin(), iHSValues.end());
		double binSize = iHSValues.size()/numBins;
		for(int i=0; i<iHSValues.size(); i++){
			int binNumber = int(i/binSize);
			if(binNumber>=numBins){
				binNumber = numBins-1;
			}
			temp[binNumber].push_back(iHSValues[i]);
			
		}
		
	}else{
		for(int i=0; i<filenames.size(); i++){
			ifstream inFile(filenames[i].c_str());
			while(!inFile.eof()){
				string l;
				getline(inFile,l);
				if(l=="" || l.find("SNP") != std::string::npos){
					continue;
				}
				vector<string> sp;
				sp = split(l, '\t');
				int anc = atoi(sp[13].c_str());
				int der = atoi(sp[12].c_str());
				double iHS = atof(sp[5].c_str());
				int bin = int(numBins*der/(anc+der));
				freqBins[bin].push_back(iHS);
				
			}
			inFile.close();
		}
	}
	//go through the entire file to bin iHS values by derived allele frequency
	
	
	
	//compute the mean and stdev within each bin
	
	for(int i=0; i<freqBins.size(); i++){
		double sum = 0;
		for(int j=0; j<freqBins[i].size(); j++){
			double d = freqBins[i][j];
			sum+=d;
		}
		double m =  sum / freqBins[i].size();

		double accum = 0.0;
		for(int j=0; j<freqBins[i].size(); j++){
			double d = freqBins[i][j];
			accum += (d - m) * (d - m);
		};

		double stdev = sqrt(accum / (freqBins[i].size()-1));
		freqBinMeans.push_back(m);
		freqBinStdevs.push_back(stdev);
	}
	
	
	//compute standardized iHS and write to file
	
	for(int i=0; i<filenames.size(); i++){
		ifstream inFile2(filenames[i].c_str());
		cout<<"Chromosome\tSNP_Index\tRecombination_Position\tGenomic_Position\tiHHA\tiHHD\tUnstd_iHS\tMin_SNP_Index\tMax_SNP_Index\tMin_SNP_Recombination_Position\tMax_SNP_Recombination_Position\tMin_SNP_Physical_Position\tMax_SNP_Physical_Position\tNum_Derived_Haplos\tNum_Ancestral_Haplos\tNum_EHH_Ancestral_Integration\tNum_EHH_Derived_Integration\tError_Code\tStd_iHS"<<endl;
		while(!inFile2.eof()){
			string l;
			getline(inFile2,l);
			if(l=="" || l.find("SNP") != std::string::npos){
				continue;
			}
			vector<string> sp;
			sp=split(l, '\t');
			int anc = atoi(sp[13].c_str());
			int der = atoi(sp[12].c_str());
			double iHS = atof(sp[5].c_str());
			int bin = int(numBins*der/(anc+der));
			double stdiHS = (iHS-freqBinMeans[bin])/freqBinStdevs[bin];
			cout<<chromosomes[i]<<"\t"<<l<<"\t"<<stdiHS<<endl;
			
		}
		inFile2.close();
	}
}

double nC2(int n){
	return double(n)*double(n-1)/2.0;
}

bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

int main(int argc, char* argv[]){
	if(argc<3){
		printHelpStatement();
		exit(0);
	}
	char*   file1name = argv[1]; //parameter file
	vector<string> filenames;
	vector<string> chromosomes;
	
	if(!fexists(file1name)){
		cerr<<"File listing iHS files to be standardized does not exist!"<<endl;
		exit(1);
	}
	
	ifstream ifile(file1name);
	if(ifile.is_open()){
		while(!ifile.eof()){
			string line;
			getline(ifile, line);
			if(line.size()<1 || line.at(0) == '#'){
				continue;
			}
			vector<string> sp;
			sp= split(line,'\t');
			if(sp.size()!=2){
				cerr<<"File listing iHS files has a line without exactly 2 tab-delilmited elements!"<<endl;
				exit(1);
			}
			filenames.push_back(sp[0]);
			chromosomes.push_back(sp[1]);
		}
	}
	
	
	numBins = atoi(argv[2]); // number of frequency bins
	if(argc >=4){
		std::string boolString(argv[3]);
		std::transform(boolString.begin(), boolString.end(), boolString.begin(), ::tolower);
		std::istringstream is(boolString);
		is >> std::boolalpha >> normalizeEqualSizedBins;
	}
	computeStandardizediHS(filenames, chromosomes);
}