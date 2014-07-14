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
#include "iHSComputer.h"
/*
 * Sandeep Venkataram and Yuan Zhu
 * Petrov Lab
 * Stanford University
 * svenkat1@stanford.edu
*/


using namespace std;

double minEHH=0.05;
double minMinorAlleleFrequency=0.05;
double maxMissingData=0.1;
int minNumHaplotypes=4; //a hidden parameter to ensure that all SNPs have at least 4 haplotypes to compute iHS for.
char singleIndexType = 'G'; //'G' signifies no single index, 'I' signifies position in haplotype, 'P' signifies physical position in genome
double singleIndexPosition=0;
bool randomize=0;
int gapCorrection=100000;
bool unequalIntegration=1;
bool thresholdSubtraction=1;
int numIterations = 1;
int numHaplos;
int haploLength;
string chromosome;
string outdir;
vector<string> haploArray; 
map<double, int> physPosToIndex;
vector<double> positions;
vector<double> physicalPositions;
vector<int> recombDistances;
vector<double> recombRates;
vector<int> singlePositions;

/*
 * Split and join code taken from Evan Terna's code at http://stackoverflow.com/questions/236129/splitting-a-string-in-c
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
void printHelpStatement(void){
	cout<<"This program computes iHS from simple SNP files"<<endl;
	cout<<"arguments: parameterFile SNPfile recombinationFile outputDirectory outputPrefix"<<endl;
	cout<<"The multithreaded program has an additional parameter specifying the number of threads to be run."<<endl;
	cout<<"\n\nThe parameter file has the following format per line: \"parameterCharacter value\""<<endl;
	cout<<"The SNP file has format e.g \"124,ATTTAATATATAAATTATTTAAAAAAAAA\" where 124 defines the physical position in the chromosme, and the 2nd part defines the SNP state for each strain in the data set. \nNOTE: The first strain in the data set will be treated as an outgroup strain used to define the ancestral and derived alleles, and will not be used in the computation. There will be errors if this is not the case."<<endl;
	cout<<"The recombination file has format \"startPosition [tab] endPosition [tab] recombinationRate\". Undefined behavior will result if not all SNPs in the SNP file have recombination rates defined by the recombination file. If there are no recombination rate estimates for a particular region, set it to 0."<<endl;
	cout<<"The SNP and recombination files must cover at most a single chromosome."<<endl;
	cout<<"Error Codes: 0 == no error, 1 == ran out of valid sequence at 5' end of region, 2 == ran out of valid sequence at 3' end of region, 3 == ran out of valid sequence at both ends of region."<<endl;
	cout<<"Example execution: './iHSComputer exampleParameters.txt exampleSNPS.csv exampleRecombMap.txt /home/user/outputDirectory/ exampleRun &'"<<endl;
	cout<<"\n\nPossible Parameters for parameter file:"<<endl;
	cout<<"S - alternative use to calculate iHS at a single SNP instead of across the entire genome. "<<endl;
	cout<<"\tG indicates that the entire genome is computed,"<<endl;
	cout<<"\tP indicates the value given in the P parameter is the genomic position of the SNP (as defined in the SNP file)"<<endl;
	cout<<"\tI indicates the value given in the P parameter is the SNP index in the haplotype starting from 1 (e.g. the (i)th SNP in the chromosome)"<<endl;
	cout<<"P - if computing over entire genome, indicates starting index of computation, otherwise sets the position value if S indicates that a single position will be used to compute iHS. Default is 1 (beginning of the chromosome)"<<endl;
	cout<<"R- 0 indicates no randomization, 1 indicates use of randomization. Randomization conducted by randomly reassigning strains as containing the ancestral or derived haplotype, instead of using the SNP state"<<endl;
	cout<<"I- number of iterations to run code. Only used if R = 1"<<endl;
	cout<<"\n\nADVANCED PARAMETERS:\n\n"<<endl;
	cout<<"E - minimum EHH cutoff for computing iHS at a site. Default is 0.05"<<endl;
	cout<<"M - minimum minor allele frequency of SNPs for which iHS is computed. Default is 0.05"<<endl;
	cout<<"D - maxmimum allowable missing data in each haplotype to include it in the computations. Default is 0.1"<<endl;
	cout<<"G - 0 indicates no gap correction, Nonzero positive integer is gap correction threshold. Default is 100,000. Gap correction minimizes the impact of long genomic stretches without SNPs on the iHS computation. Is built into the Voight iHS software."<<endl;
	cout<<"U - 0 indicates no unequal integration, 1 indicates unequal integration (iHHAncestral and iHHDerived computed by integrating EHH over different lengths of genome, based on the minimum EHH cutoff). Default is 1. Is built into the Voight iHS software."<<endl;
	cout<<"T - 0 indicates no threshold subtraction, 1 indicates threshold subtraction from EHH values for iHH integration. Threshold subtraction subtracts the min EHH cutoff from the computed EHH values while integrating. Only possible if U==1. Default is 1. Is built into the Voight iHS software."<<endl;
	
	
}

/* 
 * Read in parameters from file filename and store 
*/
void readParameters(char* filename){
	string line;
	ifstream parameterFile(filename);
	if(parameterFile.is_open()){
		string outdirLocal = outdir;
		string filename = outdirLocal.append(chromosome).append(".parameters");
		ofstream ofstr (filename.c_str(), std::ofstream::out);
		while(! parameterFile.eof()){
			getline(parameterFile, line);
			ofstr<<line<<endl;
			if (line.size() < 1 || line.at(0) == '#'){
				continue; //skip comments and empty lines
			}
			stringstream sstr(line);
			char parameter;
			sstr>>parameter;
			switch (parameter){
				case 'E':
				sstr >> minEHH;
				break;
				case 'M':
				sstr >> minMinorAlleleFrequency;
				break;
				case 'S':
				sstr >> singleIndexType;
				break;
				case 'P':
				while(!sstr.eof()){
					int p;
					sstr >> p;
					singlePositions.push_back(p);
				}
				break;
				case 'D':
				sstr >> maxMissingData;
				break;
				case 'R':
				sstr >> randomize;
				break;
				case 'I':
				sstr >> numIterations;
				break;
				case 'G':
				sstr >> gapCorrection;
				break;
				case 'U':
				sstr >> unequalIntegration;
				break;
				case 'T':
				sstr >> thresholdSubtraction;
				break;
				default:
				{
					exit(1);
				}
			}
		}
		ofstr.close();
		parameterFile.close();
		if(singlePositions.size()==0){
			singlePositions.push_back(1);
		}
		if(!unequalIntegration){
			thresholdSubtraction=0;
		}
	}else{
		exit(1);
	}
	
}

/*
 * Read in tab-delimited recombination map file
 */ 
void readRecombFile(char * filename){
     ifstream recombFile(filename);
     if(recombFile.is_open()){
       double recombSum=0;
       while(! recombFile.eof()){
         string line;
         getline(recombFile,line);
         vector<string> recombVec;
         recombVec = split(line, '\t');
         if(recombVec.size()!=3){
           continue;
         }
         int begin = atoi(recombVec[0].c_str());
         int end = atoi(recombVec[1].c_str());
         double recombRate = atof(recombVec[2].c_str());
         recombDistances.push_back(end-begin+1);
         recombRates.push_back(recombRate);         
       }
     }
}

void readSNPFile(char* filename){
	ifstream SNPfile(filename);
	if(SNPfile.is_open()){
		while(! SNPfile.eof()){
			string snp;
			getline(SNPfile, snp);
			vector<string> snpVec;
			snpVec = split(snp, ',');
			if(snpVec.size()!=2){ //if the line is not properly defined, skip it
				continue;
			}
			
			// count the number of alleles at this site.
			double physPos = atof(snpVec[0].c_str());
			char ancestral = snpVec[1][0];
			int numAlleles = 0;
			if(snpVec[1].find('A')!=string::npos){
				numAlleles++;
			}
			if(snpVec[1].find('G')!=string::npos){
				numAlleles++;
			}
			if(snpVec[1].find('C')!=string::npos){
				numAlleles++;
			}
			if(snpVec[1].find('T')!=string::npos){
				numAlleles++;
			}
			if(numAlleles!=2){ // if the snp is not biallelic
				continue;
			}
			
			
			//compute the position of this site in terms of the linkage map
			double recombPos=0;
			double runningDistance=0;
			for(int j=0; j<recombDistances.size(); j++){
				if(runningDistance + recombDistances[j] > physPos){
					if(recombRates[j]<=0){
						recombPos=0;
					}else{
						recombPos+=recombRates[j]*(physPos - runningDistance);
					}
					break;
				}else{
					if(recombRates[j]<=0){
						recombPos=0;
					}else{
						recombPos+=recombRates[j]*recombDistances[j];
					}
					runningDistance+=recombDistances[j];
				}
			}
			
			
			//add in this snp to each haplotype's array
			for(int j=1; j<snpVec[1].size(); j++){
				char s = snpVec[1][j];
				string result;
				if(s == ancestral){ result = "0";}
				else if(s=='N'){ result = "2";}
				else{ result = "1";}
				if(haploArray.size()<(snpVec[1].size()-1)){ // if we still need to initialize the array
					haploArray.push_back(result);
				}else{
					haploArray[j-1].append(result);
				}
			}
			
			
			physicalPositions.push_back(physPos);
			positions.push_back(recombPos);	
			physPosToIndex[physPos] = physicalPositions.size()-1;
		}
		
		//done filling out haplotype arrays
		//if a haplotype has too much missing data, toss it
		vector<int> toRemove;
		for(int i=0; i<haploArray.size(); i++){
			int numN = std::count(haploArray[i].begin(),haploArray[i].end(),'2');
			if(double(numN) / double(haploArray[i].size()) > maxMissingData){
				toRemove.push_back(i);
			}
			haploLength=haploArray[i].size();
		}
		for(int i=toRemove.size()-1; i>0; i--){
			haploArray.erase(haploArray.begin()+toRemove[i]);
		}
		numHaplos = haploArray.size();
		
	}else{
		exit(1);
	}
	
	
	
}
double nC2(int n){ //computes N choose 2
	return double(n)*double(n-1)/2.0;
}

/* 
 * Computes the EHH value for a subset of haplotypes starting at SNP start and using the subsequent length SNPs
 * Ignores haplotypes with ambiguous haplogroup assignment, and reconsideres haplotypes with previously ambiguous assignments using the new data
 * 
 */
EHHReturn computeEHH(vector<string> strainSets, int position, vector<string> growingHaplos, int numTotalStrains){
	double EHH=0;
	int numUsed = 0;
	vector<string> haploGroups; //strains whose haplotypes are shared with at least one other strain, and thus should be considered for the next round of EHH computation
	map<string,int> allHaploHash; //hash for all haplotypes across all sets, including those from previously invalid haplotypes
	map<string,string> validHaploHash; //hash of valid haplotypes with the strains they are associated with across all strain sets
	for(int i=0; i<strainSets.size(); i++){//for each haplogroup
		stringstream ss(strainSets.at(i));
		istream_iterator<string> begin(ss);
		istream_iterator<string> end;
		vector<string> strainSet(begin, end);
		map<string,int> haploHash;
		for(int j=0; j<strainSet.size(); j++){ //for each strain in the haplogroup
			int idex = atoi(strainSet.at(j).c_str());
			string haplo = haploArray.at(idex).substr(position,1);
			growingHaplos[idex].append(haplo);
			haplo = growingHaplos[idex];
			vector<string> validHaplos;
			for(map<string, int>::iterator it = haploHash.begin(); it!=haploHash.end(); ++it){ //compute the number of haplogroups this strain belongs to
				if(haplo.compare(it->first)==0){
					validHaplos.push_back(it->first);
					break;
				}
			}
			if(validHaplos.size()==1){ //if the strain belongs to one haplogroup, add it to that haplogroup
				haploHash[haplo] +=1;
				allHaploHash[haplo] +=1;
				numUsed++;
				validHaploHash[haplo] += "\t";
				validHaploHash[haplo] += strainSet.at(j);
			}
			else if(validHaplos.size()==0){//if this strain defines a new haplogroup, create it
				haploHash[haplo]=1;
				allHaploHash[haplo]=1;
				numUsed++;
				validHaploHash[haplo]=strainSet.at(j);
			}else{// if there is ambiguity for haplogroup assignment, it is an invalid haplotype. This should never happen, so exit out.
				exit(1);
			}
		}
	}
	//compute EHH and return, along with the set of haplogroups containing more than 1 strain (no need to consider haplotype that are already unique)
	double sum = 0.0;
	for(map<string, int>::iterator it = allHaploHash.begin(); it!=allHaploHash.end(); ++it){
		sum+=((it->second)*(it->second));
		if(it->second>1){
			haploGroups.push_back(validHaploHash.at(it->first));
		}
	}
	
	if(numUsed==0){
		EHH=0;
	}else{
		numUsed = numTotalStrains;
		EHH = (sum/(numUsed*numUsed)-1./numUsed)/(1-1./numUsed);
	}
	EHHReturn rval;
	rval.EHH = EHH;
	rval.haploGroups = haploGroups;
	rval.growingHaps=growingHaplos;
	//delete invalidStrainValidHaplos;
	return rval;
}

/* 
 * compute iHS for a given core SNP index
 */
void computeiHS(int position, ofstream &outFileStream){
	vector<string> present; //strains with derived allele
	vector<string> absent; //strains with ancestral allele
	
	int leftError = 0, rightError = 0;
	for(int i=0; i<numHaplos;i++){ //determine which strains have the derived and ancestral allele
		string curSNPStr =haploArray[i].substr(position,1); 
		int curSNP = atoi(curSNPStr.c_str());
		ostringstream ss3;
		ss3 << i;
		if(curSNP==1){ //if this is the alternative allele, say that the SNP is in the strain. (derived state)
			present.push_back(ss3.str());
		}
		else if(curSNP==0){ //else, say the SNP is not in the strain (ancestral state)
			absent.push_back(ss3.str());
		}else{
		
		}
	}
	
	int numDerived = present.size();
	int numAncestral = absent.size();
	//if there are no strains with either the ancestral or the derived state, or the minor allele frequency is below the threshold, or there are less than 4 haplotypes with the minor allele, or the SNP is in a region of 0 recombination, ignore this SNP and move on
	if( (numDerived+numAncestral)==0 || (double(min(numDerived,numAncestral))/double(numDerived+numAncestral))<=minMinorAlleleFrequency || min(numDerived,numAncestral) < minNumHaplotypes || positions[position]<=0){
		return;
	}
	map<int, double> posEHHa;
	map<int, double> posEHHd;
	map<int, double> * posToEHHAncestral = &posEHHa;// maps a genomic position to an EHH value
	map<int, double> * posToEHHDerived = &posEHHd;
	if(randomize){ //if we are doing randomization, arbitrarily assign ancestral and derived state to haplotypes, keeping the number of each the same
		present.insert(present.end(),absent.begin(),absent.end());
		std::random_shuffle(present.begin(),present.end());
		vector<string> absent2 (&present[0],&present[absent.size()]);
		vector<string> present2 (&present[absent.size()],&present[present.size()]);
		present=present2;
		absent=absent2;
	}
	//compute iHH values walking 5' of the core SNP
	vector<string> presentCopy;
	vector<string> absentCopy;
	string presentString =join(present, "\t");
	string absentString =join(absent, "\t"); 
	presentCopy.push_back(presentString);
	absentCopy.push_back(absentString);
	int minIndex = computeHalfiHS(posToEHHDerived,posToEHHAncestral, presentCopy, absentCopy, position, -1,present.size(),absent.size());
	if(minIndex <= 0){
		leftError=1;
		minIndex = abs(minIndex);
	}
	
	
	//compute iHH values walking 3' of the core SNP
	presentCopy.clear();
	absentCopy.clear();
	presentString =join(present, "\t");
	absentString =join(absent, "\t"); 
	presentCopy.push_back(presentString);
	absentCopy.push_back(absentString);
	int maxIndex = computeHalfiHS(posToEHHDerived,posToEHHAncestral, presentCopy, absentCopy, position, 1,present.size(),absent.size());
	if(maxIndex <0){
		rightError=2;
		maxIndex = abs(maxIndex);
	}
	(*posToEHHAncestral)[position]=1;
	(*posToEHHDerived)[position]=1;
	
	double iHHA = 0;
	double iHHD = 0;
	int numGaps = 0;
	int numEHHIncreaseAncestral =0, numEHHIncreaseDerived=0;
	int numEHHA=0, numEHHD=0;
	//do integration equally, i.e. integrate over the same genomic region for both ancestral and derived
	if(!unequalIntegration){
		for(map<int,double>::iterator iter= posToEHHAncestral->begin(); iter!=posToEHHAncestral->end(); ++iter){
			map<int,double>::iterator iter2 = iter;
			++iter2;
			if(iter2 == posToEHHAncestral->end()){
				break;
			}
			int i = (*iter).first;
			int i2=(*(iter2)).first;
			double currentPos = positions[i];
			double nextPos = positions[i2];
			double ratio=1.;
			if(currentPos == 0){
				leftError=1;
				break;
			}
			if(nextPos==0){
				rightError=2;
				break;
			}
			if(gapCorrection>0){ // if there is a large gap between SNPs, minimize the impact on iHS
				double gap = physicalPositions[i2]-physicalPositions[i];
				if(gap>gapCorrection*3){
					ratio=0;
				}else if(gap >= gapCorrection){
					ratio = double(gapCorrection)/gap;
				}else{
					ratio=1;
				}
				
			}
			if(ratio<1){numGaps++;}
			if( (i<position && posToEHHAncestral->at(i) > posToEHHAncestral->at(i2)) || (i >= position && posToEHHAncestral->at(i) < posToEHHAncestral->at(i2))){
				numEHHIncreaseAncestral++;
			}
			if( (i<position && posToEHHDerived->at(i) > posToEHHDerived->at(i2)) || (i >= position && posToEHHDerived->at(i) < posToEHHDerived->at(i2))){
				numEHHIncreaseDerived++;
			}
			
			//trapezoid rule integration
			iHHA += ratio*(posToEHHAncestral->at(i)+posToEHHAncestral->at(i2))/2*(nextPos-currentPos);
			iHHD += ratio*(posToEHHDerived->at(i)+posToEHHDerived->at(i2))/2*(nextPos-currentPos);
			numEHHA++;
			numEHHD++;
		}
	}else{ // integrate unequally i.e. integrate iHH Ancestral and iHH Derived separately until EHH decays below minimum threshold.
		for(map<int,double>::iterator iter= posToEHHAncestral->begin(); iter!=posToEHHAncestral->end(); ++iter){ // integrate derived haplotypes
			map<int,double>::iterator iter2 = iter;
			++iter2;
			if(iter2 == posToEHHAncestral->end()){
				break;
			}
			int i = (*iter).first;
			int i2=(*(iter2)).first;
			double currentPos = positions[i];
			double nextPos = positions[i2];
			double ratio=1.;
			if(posToEHHAncestral->at(i)<minEHH && posToEHHAncestral->at(i2)<minEHH){
				continue;
			}
			if(currentPos == 0){
				leftError=1;
				break;
			}
			if(nextPos==0){
				rightError=2;
				break;
			}
			if(gapCorrection>0){
				double gap = physicalPositions[i2]-physicalPositions[i];
				if(gap>gapCorrection*3){
					ratio=0;
					continue;
				}else if(gap >= gapCorrection){
					ratio = double(gapCorrection)/gap;
				}else{
					ratio=1;
				}
				
			}
			if(ratio<1){numGaps++;}
			double currentVal=0.0;
			if(!thresholdSubtraction){
				currentVal= ratio*(posToEHHAncestral->at(i)+posToEHHAncestral->at(i2))/2*(nextPos-currentPos); //scale EHH by the distance it spans
			}else{ // subtract off minEHH cutoff from integration.
				if(posToEHHAncestral->at(i)<minEHH){
					currentVal= ratio*(nextPos-currentPos)*(1-(minEHH-posToEHHAncestral->at(i))/(posToEHHAncestral->at(i2)-posToEHHAncestral->at(i)))*.5*(posToEHHAncestral->at(i2)-minEHH);
				}else if(posToEHHAncestral->at(i2)<minEHH){
					currentVal= ratio*(nextPos-currentPos)*(1-(minEHH-posToEHHAncestral->at(i2))/(posToEHHAncestral->at(i)-posToEHHAncestral->at(i2)))*.5*(posToEHHAncestral->at(i)-minEHH);
				}else{
					currentVal= ratio*(posToEHHAncestral->at(i)+posToEHHAncestral->at(i2)-(minEHH*2))/2*(nextPos-currentPos);
				}
			}
			if( (i<position && posToEHHAncestral->at(i) > posToEHHAncestral->at(i2)) || (i >= position && posToEHHAncestral->at(i) < posToEHHAncestral->at(i2))){
				numEHHIncreaseAncestral++;
			}
			numEHHA++;
			iHHA+=currentVal;
		}
		//repeat integration for ancestral haplotypes
		for(map<int,double>::iterator iter= posToEHHDerived->begin(); iter!=posToEHHDerived->end(); ++iter){
			map<int,double>::iterator iter2 = iter;
			++iter2;
			if(iter2 == posToEHHDerived->end()){
				break;
			}
			int i = (*iter).first;
			int i2=(*(iter2)).first;
			double currentPos = positions[i];
			double nextPos = positions[i2];
			double ratio=1.;
			if(posToEHHDerived->at(i)<minEHH && posToEHHDerived->at(i2)<minEHH){
				continue;
			}
			if(currentPos == 0){
				//cout<<position<<" hit the end of the chromosome during integration!"<<endl;
				leftError=1;
				break;
			}
			if(nextPos==0){
				rightError=2;
				break;
			}
			if(gapCorrection>0){
				double gap = physicalPositions[i2]-physicalPositions[i];
				if(gap>gapCorrection*3){
					ratio=0;
					continue;
				}else if(gap >= gapCorrection){
					ratio = double(gapCorrection)/gap;
				}else{
					ratio=1;
				}
				
			}
			if(ratio<1){numGaps++;}
			double currentVal=0.0;
			if(!thresholdSubtraction){
				currentVal= ratio*(posToEHHDerived->at(i)+posToEHHDerived->at(i2)-.1)/2*(nextPos-currentPos);
			}else{
				if(posToEHHDerived->at(i)<minEHH){
					currentVal= ratio*(nextPos-currentPos)*(1-(minEHH-posToEHHDerived->at(i))/(posToEHHDerived->at(i2)-posToEHHDerived->at(i)))*.5*(posToEHHDerived->at(i2)-minEHH);
				}else if(posToEHHDerived->at(i2)<minEHH){
					currentVal= ratio*(nextPos-currentPos)*(1-(minEHH-posToEHHDerived->at(i2))/(posToEHHDerived->at(i)-posToEHHDerived->at(i2)))*.5*(posToEHHDerived->at(i)-minEHH);
				}else{
					currentVal= ratio*(posToEHHDerived->at(i)+posToEHHDerived->at(i2)-(minEHH*2))/2*(nextPos-currentPos);
				}
			}
			if( (i<position && posToEHHDerived->at(i) > posToEHHDerived->at(i2)) || (i >= position && posToEHHDerived->at(i) < posToEHHDerived->at(i2))){
				numEHHIncreaseDerived++;
			}
			numEHHD++;
			iHHD+=currentVal;
		}
	}
	if(iHHD>0){//write to file
		outFileStream<<fixed<<(position+1)<<"\t"<<positions[position]<<"\t"<<physicalPositions[position]<<"\t"<<iHHA<<"\t"<<iHHD<<"\t"<<log(iHHA/iHHD)<<"\t"<<minIndex<<"\t"<<maxIndex<<"\t"<<positions[minIndex]<<"\t"<<positions[maxIndex]<<"\t"<<physicalPositions[minIndex]<<"\t"<<physicalPositions[maxIndex]<<"\t"<<numDerived<<"\t"<<numAncestral<<"\t"<<numEHHA<<"\t"<<numEHHD<<"\t"<<(leftError+rightError)<<endl;//describes the output file
	}
}

/*
 * compute iHS in one direction (either walking 3' or 5' from core SNP
 */
int computeHalfiHS(map<int, double> * EHHd, map<int, double> * EHHa, vector<string> present, vector<string> absent, int position, int indexIncrement, int numPresent, int numAbsent){
	int curIndex = position+indexIncrement;
	int extremeIndex=-1;
	vector<string> presentLocal =present;
	vector<string> absentLocal = absent;
	vector<string> growingHaplos(numHaplos);
	int maxCountInvalid = 0;
	
	stringstream ssPres(present.at(0));
	istream_iterator<string> beginPres(ssPres);
	istream_iterator<string> endPres;
	vector<string> presSet(beginPres, endPres);
	for(int j=0; j<presSet.size(); j++){
		int idex = atoi(presSet.at(j).c_str());
		growingHaplos[idex]="1";		
	}
	stringstream ssAbs(absent.at(0));
	istream_iterator<string> beginAbs(ssAbs);
	istream_iterator<string> endAbs;
	vector<string> absSet(beginAbs, endAbs);
	
	for(int j=0; j<absSet.size(); j++){
		int idex = atoi(absSet.at(j).c_str());
		growingHaplos[idex]="0";		
	}
	while(curIndex>0 && curIndex < haploLength-1){ //work my way either 3' or 5' of the core SNP. The method is essentially a wierd kind of recursion, where the computeEHH call is repeated with values generated by the previous call of the method. however, this is done using a while loop and not true recursive calls
		if(positions[curIndex]==0){ //if the index is in a bad region of the genome (no available recombination map), then return out
			return -curIndex;
		}
		int start = 0; 
		int length = 0;
		//initialize iterators
		stringstream ss(present.at(0));
		istream_iterator<string> begin(ss);
		istream_iterator<string> end;
		vector<string> strainSet(begin, end);
		stringstream ss2(absent.at(0));
		istream_iterator<string> begin2(ss2);
		istream_iterator<string> end2;
		strainSet.insert(strainSet.end(),begin2,end2);
		bool done=0;
		string str;
		//add in one more site worth of information onto our haplotypes
		for(int j=0; j<strainSet.size(); j++){
			int idex = atoi(strainSet.at(j).c_str());
			string st2=haploArray.at(idex).substr(curIndex,1);
			str.append(st2);
			if(st2.compare("2")==0){
				done=1;
				break;
			}
		}
		if(done){
			curIndex+=indexIncrement;
			continue;
		}
		
		
		if(indexIncrement==1){ //index housekeeping depending on whether we are moving forwards or backwards in the genome
			start = position;
			length = curIndex + 1 - position;
		}else{
			start = curIndex;
			length = position + 1 - curIndex;
		}
		//actually compute EHH and update our arrays accordingly
		EHHReturn presentEHH = computeEHH(presentLocal,curIndex,growingHaplos, numPresent);
		growingHaplos=presentEHH.growingHaps;
		EHHReturn absentEHH =  computeEHH(absentLocal,curIndex,growingHaplos, numAbsent);
		growingHaplos=absentEHH.growingHaps;
		(*EHHd)[curIndex]=presentEHH.EHH;
		(*EHHa)[curIndex]=absentEHH.EHH;
		if(presentEHH.EHH < minEHH && absentEHH.EHH < minEHH){
			extremeIndex=curIndex;
			break;
		}
		presentLocal = presentEHH.haploGroups;
		absentLocal = absentEHH.haploGroups;
		extremeIndex = curIndex;
		curIndex+=indexIncrement;
	}
	if(curIndex <= 0 || curIndex >= haploLength-1){ //if the index is in a bad region of the genome (no available recombination map), then return out
		return -curIndex;
	}
	return extremeIndex;
}


bool fexists(const char *filename) //boolean to check if a file exists
{
  ifstream ifile(filename);
  return ifile;
}
