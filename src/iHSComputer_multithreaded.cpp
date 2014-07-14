#include "iHSComputer.h"
#include <thread>



/*
 * Sandeep Venkataram and Yuan Zhu
 * Petrov Lab
 * Stanford University
 * svenkat1@stanford.edu
*/

using namespace std;

int numThreads;


/*
 * Merges the iHS parts computed by the separate threads into a single file for ease of use.
 */
void createMergedFile(void){
	ifstream streams[numThreads];
	
	string lines[numThreads];
	string outfiledir = outdir;	
	string outfile = outfiledir.append(chromosome).append(".iHS");
	ofstream outFile(outfile.c_str(), std::ofstream::out);
	outFile<<"SNP_Index\tRecombination_Position\tGenomic_Position\tiHHA\tiHHD\tUnstd_iHS\tMin_SNP_Index\tMax_SNP_Index\tMin_SNP_Recombination_Position\tMax_SNP_Recombination_Position\tMin_SNP_Physical_Position\tMax_SNP_Physical_Position\tNum_Derived_Haplos\tNum_Ancestral_Haplos\tNum_EHH_Ancestral_Integration\tNum_EHH_Derived_Integration\tError_Code"<<endl;
	for(int threadNum=0; threadNum<numThreads; threadNum++){
		ostringstream ss1;
		ss1 << (threadNum+1);
		ostringstream ss2;
		ss2 << numThreads;
		string outdirLocal = outdir;	
		string filename = outdirLocal.append(chromosome).append("_part_").append(ss1.str()).append("_of_").append(ss2.str()).append(".iHS");
		streams[threadNum].open(filename.c_str());
		if(!streams[threadNum].eof()){
			getline(streams[threadNum],lines[threadNum]);
		}
	}
	
	bool done = 0;
	while(!done){
		string line;
		int usedFile = 0;
		int lastIndex = std::numeric_limits<int>::max();
		done =1;
		for(int threadNum=0; threadNum<numThreads; threadNum++){
			string l = lines[threadNum];
			if(l==""){
				continue;
			}
			done=0;
			vector<string> sp;
			sp = split(l, '\t');
			int p = atoi(sp[0].c_str());
			if(p<lastIndex){
				lastIndex=p;
				usedFile=threadNum;
				line = l;
			}
		}
		if(done){
			break;
		}
		outFile<<line<<endl;
		lines[usedFile]="";
		if(!streams[usedFile].eof()){
			getline(streams[usedFile],lines[usedFile]);
		}
	}
	outFile.close();
	
	for(int threadNum=0; threadNum<numThreads; threadNum++){
		ostringstream ss1;
		ss1 << (threadNum+1);
		ostringstream ss2;
		ss2 << numThreads;
		string outdirLocal = outdir;	
		string filename = outdirLocal.append(chromosome).append("_part_").append(ss1.str()).append("_of_").append(ss2.str()).append(".iHS");
		remove(filename.c_str());
	}
}

void threadFunction(int threadNum, int startingPosition){ // function to compute iHS values within each thread
	ostringstream ss1;
	ss1 << (threadNum+1);
	ostringstream ss2;
	ss2 << numThreads;
	string outdirLocal = outdir;
	string filename = outdirLocal.append(chromosome).append("_part_").append(ss1.str()).append("_of_").append(ss2.str()).append(".iHS");
	ofstream outFileStream(filename.c_str(), std::ofstream::out | std::ofstream::app);
	for(int k=threadNum+startingPosition; k<haploLength-1;k+=numThreads){
		computeiHS(k,outFileStream);
	}	
	outFileStream.close();
	
}

int main(int argc, char* argv[]){
	if(argc!=7){
		printHelpStatement();
		exit(0);
	}
	char*   file1name = argv[1]; //parameter file
	char*   file2name = argv[2]; //SNP file
	char*   file3name = argv[3]; //recomb file
	outdir = argv[4]; //out directory
	chromosome =argv[5]; //out directory
	numThreads = atoi(argv[6]);
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
			outFileStream<<fixed;
			if(!fileExists){
				outFileStream<<"SNP_Index\tRecombination_Position\tGenomic_Position\tiHHA\tiHHD\tUnstd_iHS\tMin_SNP_Index\tMax_SNP_Index\tMin_SNP_Recombination_Position\tMax_SNP_Recombination_Position\tMin_SNP_Physical_Position\tMax_SNP_Physical_Position\tNum_Derived_Haplos\tNum_Ancestral_Haplos\tNum_EHH_Ancestral_Integration\tNum_EHH_Derived_Integration"<<endl;
			}
			for(int l=0;l<numIterations;l++){
				computeiHS(pos,outFileStream);
			}
			outFileStream.close();
		}
	}else{ // if we are computing across the genome, multithread, merge files and standardize
		singleIndexPosition=singlePositions[0];
		vector<thread> threads;
		for(int i=0; i < numThreads; i++ ){
			threads.push_back( thread(&threadFunction, i, singleIndexPosition-1));
		}
		for (auto& t : threads)
			t.join();

		createMergedFile();
	}
	
}
