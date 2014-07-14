#include "iHSComputer.h"

using namespace std;

/*
 * Sandeep Venkataram and Yuan Zhu
 * Petrov Lab
 * Stanford University
 * svenkat1@stanford.edu
*/


void threadFunction(int threadNum, int startingPosition){ // function to compute iHS values within each thread
	string outdirLocal = outdir;
	string filename = outdirLocal.append(chromosome).append(".iHS");
	ofstream outFileStream(filename.c_str(), std::ofstream::out | std::ofstream::app);
	for(int k=threadNum+startingPosition; k<haploLength-1;k++){
		computeiHS(k,outFileStream);
	}
	outFileStream.close();
	
}

int main(int argc, char* argv[]){
	if(argc!=6){
		printHelpStatement();
		exit(0);
	}
	char*   file1name = argv[1]; //parameter file
	char*   file2name = argv[2]; //SNP file
	char*   file3name = argv[3]; //recomb file
	outdir = argv[4]; //out directory
	chromosome =argv[5]; //out directory
	readParameters(file1name);
	readRecombFile(file3name);
	readSNPFile(file2name);
	if(singleIndexType != 'G'){	//if we are not computing across the genome, compute iHS for the given site numIterations number of times.
		for(int k=0; k<singlePositions.size(); k++){
			singleIndexPosition=singlePositions[k];
			int pos;
			string type;
			switch (singleIndexType){
				case 'P':
					pos = physPosToIndex[singleIndexPosition];
					type = "_genomePosition";
					break;
				case 'I':
					pos = singleIndexPosition-1;
					type  ="_haplotypeIndex";
					break;			
				default:
				{
					exit(1);
				}
			}
			string filename = outdir;
			filename.append(chromosome).append(".iHS");
			bool fileExists = fexists(filename.c_str());
			ofstream outFileStream(filename.c_str(), std::ofstream::out | std::ofstream::app);
			if(!fileExists){
				outFileStream<<"SNP_Index\tRecombination_Position\tGenomic_Position\tiHHA\tiHHD\tUnstd_iHS\tMin_SNP_Index\tMax_SNP_Index\tMin_SNP_Recombination_Position\tMax_SNP_Recombination_Position\tMin_SNP_Physical_Position\tMax_SNP_Physical_Position\tNum_Derived_Haplos\tNum_Ancestral_Haplos\tNum_EHH_Ancestral_Integration\tNum_EHH_Derived_Integration\tError_Code"<<endl;
			}
			if(!randomize){
				numIterations=1;
			}
			for(int l=0;l<numIterations;l++){
				computeiHS(pos,outFileStream);
			}
			outFileStream.close();
		}
	}else{ // if we are computing across the genome, multithread, merge files and standardize
		singleIndexPosition=singlePositions[0];
		threadFunction(0, singleIndexPosition-1);
	}
	
}