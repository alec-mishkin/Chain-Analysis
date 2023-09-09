#include <fstream>
#include <iterator>
#include <iostream>
#include <string>
#include <cstring>
#include <time.h>
#include <sstream>
#include <stdio.h>
#include <stdio.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <ctime>
#include <sys/stat.h> 
#include <sys/types.h> 
using namespace std;

const int GLOBAL_N_MAX =3000; 
const int GLOBAL_N_NEIGHBORS = 50;
const int N_Bond_Types = 10;
const int N_EDGES = 10000;
float BIN_SIZE;
int HISTOGRAM_SIZE;
float DOI_BEG ;
float DOI_END;
int numberOfChains;


struct Bond //This structure contains the information for one bond type
{
        string startAtom; //The starting atom type in text 
        string endAtom; //The ending atom type in text
        float cutOffDist; //The maximum distance allowed between the two atom types for it to be considered a bond
};


struct ConfigData
{
	string inputFile;
	string outputFile;
	string allowedBond;
	string cutOffDistances;
	string binSize;
	string histogramSize;
	string doiBeg;
	string doiEnd;
	string nChains;
};

/*-----------------------------------------------------------------
 *The atom object with the necessary attributes for connecting edges between them
 *----------------------------------------------------------------
 */
struct XYZAtomData
{
	string atomType;
	int atomIndex;
	float xPos;
	float yPos;
	float zPos;

	void GetInfo()
	{
		cout << atomType << xPos << yPos << zPos << endl;
	}
};

struct Bonds //This structure contains all the information for allowed bonds
{
	int nBonds; //The number of bond types in the atomic structure being analyzed
	Bond bonds[N_Bond_Types]; //An array of the bond structures 	
	Bond tempBond; //A bond structure used to store specific bond date when checking conditions

	//Check if the two atom types can possibly make a bond
	//TODO REMOVE THIS
	bool CheckBondAtoms(string atom1, string atom2) 
	{
		//Loop through all possible defined bonds
		for (int bondIndex = 0; bondIndex < nBonds; bondIndex = bondIndex +1)
		{
			//If there exists a bond struture where atom1 and atom2 are the start atom and end atom or vice versa, than these atoms may make a bond
			if((atom1 == bonds[bondIndex].startAtom and atom2 == bonds[bondIndex].endAtom) or (atom1 == bonds[bondIndex].endAtom and atom2 == bonds[bondIndex].startAtom))	
			{
				return true;
			}
		}
		
		return false;
	}
	//TODO REMOVE THIS
	bool CheckBondDistances(float distance)
	{
		for (int bondIndex = 0; bondIndex < nBonds; bondIndex = bondIndex +1)
		{
			if((distance < bonds[bondIndex].cutOffDist))
			{
				return true;
			}
		}
		return false;
	}


	//Check the atom1 and atom2 types and the distance between them to see if they form a bond
	bool CheckIfBond(string atom1, string atom2, float distance)
	{
		//Loop through all possible defined bonds
		for (int bondIndex = 0; bondIndex < nBonds; bondIndex = bondIndex +1)
		{
			//If there exists a bond struture where atom1 and atom2 are the start atom and end atom or vice versa and the distance between these taoms is less than that same bond structures maximum distance, a bond is formed
			if((atom1 == bonds[bondIndex].startAtom and atom2 == bonds[bondIndex].endAtom) or (atom1 == bonds[bondIndex].endAtom and atom2 == bonds[bondIndex].startAtom) and distance < bonds[bondIndex].cutOffDist)
			{
				return true;
			}
		}
		return false;

	}
	//Output the starting atom type and ending atom type and maximum bond distance for all defined bonds
	void GetInfo()
	{
                for (int bondIndex = 0; bondIndex < nBonds; bondIndex = bondIndex +1)
                {
			cout << bonds[bondIndex].startAtom << " " << bonds[bondIndex].endAtom << " " << bonds[bondIndex].cutOffDist << endl;
		}

	}
};

/*--------------------------------------------------------------
 *The XYZ structure that contains all of the atoms
 * -------------------------------------------------------------
 */
struct XYZData
{
        float xLength;
	float yLength;
	float zLength;
	int nAtoms;
	XYZAtomData xYZAtomData[GLOBAL_N_MAX]; //You can change the size of GLOBAL_N_MAX but as long as it is bigger then the number of atoms in your structure it should be fine
};

/*------------------------------------------------------------------
 *The Neighbor structure is used for the Djikstra algorithm. 
 It allows the code to keep track of which atoms are connected to eachother
 * -----------------------------------------------------------------
 */
struct Neighbor
{
int atomId;
float distance;
};


/*-----------------------------------------------------------------
 *Vertices contain what atom id they are, a distance value they can updat  TODO????
 * Whether or not they have already been visited 
 * The number of neighbors they have been connected to
 * The array of neighbors they connected to
 * prevId TODO?????
 * ----------------------------------------------------------------
 */
struct Vertex
{
int atomId;
float distance;
bool visited;
int nNeighbors;
int prevId;
Neighbor neighbors[GLOBAL_N_NEIGHBORS]; //Make sure this is larger than the max coordination number of your atoms
};


/*-----------------------------------------
 *TODO: I am confused why there are two different vertices 
 TODO Next, look at what initialization function is being called
 Looks like this is really only being used by the MakeGraph Function
 * -------------------------------------------
 */
struct Vertices
{

	Vertices() //Empty Contrustor
	{
	}	

	Vertices(Vertex vertices[GLOBAL_N_MAX]) //The Vertices strucure can also be contstructed from an array of vertex
	{

		for (int i = 0; i < GLOBAL_N_MAX; i++)
		{	
			verts[i]= vertices[i];
		}
	}
		
	Vertex verts[GLOBAL_N_MAX]; //Create an array of vertex
	void Initialize() //Initialize the array
	{
		//Every vertex starts with an atom id, 0 neighbors, and a distance of infinity
		for (int i = 0; i < GLOBAL_N_MAX; i++)
		{
			verts[i].atomId = i;
			verts[i].distance = std::numeric_limits<float>::infinity();
			verts[i].nNeighbors = 0;
		}
	}
	
	//All vertices start with the a distance of infinity to every other atom except them seleves
	void Initialize(int nAtoms, int sourceAtomId)
	{
		for (int i = 0; i < nAtoms; i++)
		{
			verts[i].atomId = i;
			if(verts[i].atomId == sourceAtomId)
			{
				verts[i].distance = 0;
			}
			else
			{
				verts[i].distance = std::numeric_limits<float>::infinity();
			}
		
			verts[i].nNeighbors = 0; 
		}
	}
};

/*-------------------------------------------------------------------------
 *This structure represents a path of atoms or nodes connected by edges 
 * 
 *------------------------------------------------------------------------ 
 */
struct Path
{
	int pathIds[3000];
	int nPathAtoms;
	void GetInfo()
	{
		for (int pathIndex = 0; pathIndex < nPathAtoms; pathIndex = pathIndex +1)
		{
		cout << pathIds[pathIndex] << " " ;
		}
		cout << endl; 
		}
	//TODO What does this do????
	void AppendPath(Path path)
	{
		
		int oldPathIndex = 0 ;
		for(int j = nPathAtoms; j < nPathAtoms + path.nPathAtoms; j++)
		{
			
			pathIds[j] = path.pathIds[oldPathIndex];
			oldPathIndex = oldPathIndex + 1;
		}
		nPathAtoms = nPathAtoms + path.nPathAtoms;
		
	}
};

ConfigData ReadInConfigFile(const char* filename)
{
	
	ConfigData configData = {"empty","empty","empty","empty"};
	ifstream inFile;
        inFile.open(filename);
	int lineNumber = 0;
	char configValues[4][100];
	char valueName[100];
	char value[100];
        if (!inFile)
        {
            cerr << "Unable to open file configfile";
        }
	
       while (inFile >> valueName >> value)
        {
		if(strcmp(valueName, "INPUTFILE:") == 0)
		{
			configData.inputFile = value; 
		}
		else if(strcmp(valueName, "OUTPUTFILE:") == 0)
		{
			configData.outputFile = value;
		}
                else if(strcmp(valueName, "ALLOWEDATOMBONDS:") == 0)
                {
                       configData.allowedBond = value;
                }
                else if(strcmp(valueName, "CUTOFFDISTANCES:") == 0)
                {
                       configData.cutOffDistances = value;
                }

                else if(strcmp(valueName,"BIN_SIZE:") == 0)
                {
                        configData.binSize = value;
			BIN_SIZE = atof(configData.binSize.c_str());
                }
                else if(strcmp(valueName,"HISTOGRAM_SIZE:") == 0)
                {
                        configData.histogramSize = value;
			HISTOGRAM_SIZE = atoi(configData.histogramSize.c_str());
                }
                else if(strcmp(valueName,"DOI_BEG:") == 0)
                {
                        configData.doiBeg = value;
                        DOI_BEG = atof(configData.doiBeg.c_str());
                }
                else if(strcmp(valueName,"DOI_END:") == 0)
                {
                        configData.doiEnd = value;
			DOI_END = atof(configData.doiEnd.c_str());
                }
                else if(strcmp(valueName,"NCHAINS:") == 0)
                {
                        configData.nChains = value;
			numberOfChains = atoi(configData.nChains.c_str());
                }

		
		++lineNumber;
	}
	return configData; 
	

}

XYZData ReadInXYZFile(const char* filename)
{
	ifstream inFile;
	inFile.open(filename);
	int lineNumber = 0;
	char atomType[100];
	int atomIndex;
	int nAtoms;
	Bonds bonds;
	string line;
	XYZAtomData atomData;
	XYZData xYZData;
	while (inFile)
	{
			lineNumber++;
			if(lineNumber == 1)
			{
				inFile >> nAtoms;
			}
			else if(lineNumber == 2)
			{
				inFile>> xYZData.xLength >> xYZData.yLength >> xYZData.zLength;	
			}
			else if (lineNumber > 2 and nAtoms >= (lineNumber - 2))
			{
			
				atomData.atomIndex = lineNumber - 3;
				inFile >> atomData.atomType >> atomData.xPos >> atomData.yPos >> atomData.zPos;
				xYZData.xYZAtomData[lineNumber-3] = atomData;
			}
			else if (nAtoms < lineNumber  - 2)
			{
				//TODO need to make a better fix for the extra line that was added
				inFile >> atomData.atomType >> atomData.xPos >> atomData.yPos >> atomData.zPos;
				atomData.atomIndex = lineNumber - 2;
			}
	}
	xYZData.nAtoms = lineNumber - 3;//TODO need to make a better fix for the extra line that was added
	if(nAtoms != xYZData.nAtoms)
	{
		cout << "nAtoms" << nAtoms << endl;
		cout << "xYZData.nAtoms" << xYZData.nAtoms << endl;
		cout << "n Atoms ERROR" << endl; 
	}
	return xYZData;
}

Bonds LoadBonds(string allowBonds, string allowCutOffs)
{
	Bonds givenBonds;	
	cout << allowBonds << endl;
	string tempAtom="";
	int bondNumber = 0;
	int numStartAtoms;
	int numDists;
        
	for (int i = 0; i <allowBonds.length(); i++)
	{
		cout << allowBonds[i] << endl;
		if(allowBonds[i] == '-')
		{
			givenBonds.bonds[bondNumber].startAtom = tempAtom;		
			tempAtom = "";
		}
		else if(allowBonds[i] == ',')
		{
                        givenBonds.bonds[bondNumber].endAtom = tempAtom;
                        tempAtom = "";
			bondNumber = bondNumber + 1;

		}
		else
		{
			tempAtom += allowBonds[i];
		}
	}
	givenBonds.bonds[bondNumber].endAtom = tempAtom;
	numStartAtoms = bondNumber + 1;
        
	string tempDist = "";
	float dist;
	bondNumber = 0;
	for (int i = 0; i <allowCutOffs.length(); i++)
	{
                if(allowCutOffs[i] == ',')
                {
			givenBonds.bonds[bondNumber].cutOffDist = std::atof(tempDist.c_str());
			tempDist = "";
			bondNumber = bondNumber + 1;
		}
		else
		{
			cout << allowCutOffs[i] << endl;
			tempDist += allowCutOffs[i];
		}

	}
	givenBonds.bonds[bondNumber].cutOffDist = std::atof(tempDist.c_str());
	numDists = bondNumber + 1;
	if(numDists != numStartAtoms)
	{
		cout << "THE NUMBER OF DISTANCES IN THE CONFIG SHEET DID NOT EQUAL THE NUMBER OF STARTING ATOMS IN THE CONFIG SHEET" << endl;
		exit(EXIT_FAILURE);
	}
	givenBonds.nBonds = numDists;
	return givenBonds;
}

/*------------------------------------------------------------------------------------------------
 *Calculate the periodic distance between atoms TODO double check this
 * -----------------------------------------------------------------------------------------------
 */
float CalculateDistance(XYZAtomData atomData1, XYZAtomData atomData2, float xLength, float yLength, float zLength)
{
	float rVector1 = (1/xLength);
	float rVector2 = (1/yLength);
	float rVector3 = (1/zLength);
	float diff_x = atomData1.xPos*rVector1 - atomData2.xPos*rVector1;
        float box_dx = diff_x - round(diff_x);
        float diff_y = atomData1.yPos*rVector2 - atomData2.yPos*rVector2;
        float box_dy = diff_y - round(diff_y);
        float diff_z = atomData1.zPos*rVector3 - atomData2.zPos*rVector3;
        float box_dz = diff_z - round(diff_z);
	float xDist = box_dx * xLength;
	float yDist = box_dy * yLength;
	float zDist = box_dz * zLength;
	float dist = sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
	return dist;
}
Vertices MakeGraph(XYZData xYZData, Bonds bonds) //Takes in the XYZ data and all the possbile type of bonds being used
{
        Vertices vertices; //Create an object to hold the vertices
	vertices.Initialize(); //Give all vertices their own atomId, nNeigbhors variable set to 0, and a distance of infinity
	for (int atom1Index = 0; atom1Index < xYZData.nAtoms; atom1Index = atom1Index +1) //loop through all pairs of atoms
        {
                for (int atom2Index = 0; atom2Index < xYZData.nAtoms; atom2Index = atom2Index +1)
                {
                        
			
			if(bonds.CheckBondAtoms(xYZData.xYZAtomData[atom1Index].atomType, xYZData.xYZAtomData[atom2Index].atomType) and atom1Index != atom2Index) //Check if two atoms form a bond based on type and distance 
                        {
				float bondDist = CalculateDistance(xYZData.xYZAtomData[atom1Index], xYZData.xYZAtomData[atom2Index], xYZData.xLength, xYZData.yLength, xYZData.zLength); //Calculate the distance between the two atoms
				if (CheckIfBond(xYZData.xYZAtomData[atom1Index].atomType, xYZData.xYZAtomData[atom2Index].atomType, bondDist)) //Check if the distance between the two atom types is less than the cutoff bond distance
				{
					//Give the appropriate vertex its neighbor with its corresponding distance 
					vertices.verts[atom1Index].nNeighbors = vertices.verts[atom1Index].nNeighbors + 1; //Add 1 to the the total number of neighbors connected to the atom1 Vertex
					vertices.verts[atom1Index].neighbors[(vertices.verts[atom1Index].nNeighbors - 1)].atomId = atom2Index; //Give that neighbor the atom2Index
					vertices.verts[atom1Index].neighbors[(vertices.verts[atom1Index].nNeighbors - 1)].distance = bondDist;	//Give the vertex the distance to the newly added neighbor

				}
                        }
                }
        }
	//Return a graph with all bonded atoms
	return vertices;
}

int FindAtomPairs(float distance_of_interest, float binSize)
{
    int binIndex = round(distance_of_interest/binSize);
    return binIndex;
}

void WriteXYZFile(string outputDirectory, string startAtomType, string endAtomType, string xyzPath, int nChain, int chainNumber)
{
	string fileName = outputDirectory + startAtomType + endAtomType;
	std::ofstream outxyz(fileName.c_str(),std::ios_base::app);
	if(chainNumber < nChain)
	{
	outxyz << xyzPath << ":: ";
	}
	else
	{
	outxyz << xyzPath << endl;
	}
	outxyz.close();
}
string MakeAndWritePath(string outputFolder, int binIndex, Path path, XYZData xYZData, int nChain, int chainNumber)
{
	string stringPath = "";
	string xyzPath = "";
	std::ostringstream idss;
	idss << path.pathIds[0] << " " << path.pathIds[path.nPathAtoms-1] << " CN: " << chainNumber;
	xyzPath += idss.str();
	xyzPath  += " - "; //Put in the initial atom ids
	std::ofstream out("/ufrc/cheng/amishkin/BONDPATH_C++_TEST/output.xyz");
	string startAtomType = xYZData.xYZAtomData[path.pathIds[0]].atomType;	
        string endAtomType = xYZData.xYZAtomData[path.pathIds[path.nPathAtoms-1]].atomType;
	for (int pathIndex = 0; pathIndex < path.nPathAtoms; pathIndex++)
	{
		std::ostringstream ss;
		int id = path.pathIds[pathIndex];
		ss << id;
		string atomType = xYZData.xYZAtomData[id].atomType;
		stringPath += ss.str();
                std::ostringstream ssxyz;
                ssxyz << id  << " " << xYZData.xYZAtomData[id].atomType << " " <<  xYZData.xYZAtomData[id].xPos << " " << xYZData.xYZAtomData[id].yPos << " " << xYZData.xYZAtomData[id].zPos << " ";
                xyzPath += ssxyz.str();
		stringPath += ",";
	}
	float doi = binIndex * BIN_SIZE;
	std::ostringstream binString;
        binString << doi;

        string outputDirectory = outputFolder + "/" +binString.str() + "/";
        if (mkdir(outputDirectory.c_str(),0777) != 0)
        {
        mkdir(outputDirectory.c_str(),0777);
        }

        WriteXYZFile(outputDirectory, startAtomType, endAtomType, xyzPath, nChain, chainNumber);

	
	return stringPath;	
}



Vertex GetMinUnVisitedVertexId(Vertices vertices,XYZData xYZData)
{
        Vertex minVertex;
	minVertex.atomId = -1; //The minVertex atom id is instantiated at -1 so that is not a possible id 
        float minDistance = std::numeric_limits<float>::infinity();
        for (int i = 0; i <xYZData.nAtoms; i++)
        {
                if(vertices.verts[i].visited == false) //A visited atom can not be the minVisit atom
                {
                        if(vertices.verts[i].distance < minDistance) // If the distance is less than the min diststance that that can be the new minVertex
                        {
                                minVertex = vertices.verts[i];
				minDistance = vertices.verts[i].distance;  
                        }

                }
        }
//If the remainining not visited atoms have a distance of infinity that means they are not connected to any of the atoms in the path we already have
        return minVertex;
}


Path djikstra(Vertices vertices, int sourceId, int endId,  XYZData xYZData, Path oldPath)
{
	//To make a ring all I have to do is make a for loop that will go through all the oldPath and call it visited
	vertices.verts[endId].visited = false; //We want to be able to visit the ending atom so mark that as unvisited
	for (int pathIndex = 0; pathIndex < oldPath.nPathAtoms; pathIndex = pathIndex + 1)
	{
		if(oldPath.pathIds[pathIndex] != -1)
		{
		vertices.verts[oldPath.pathIds[pathIndex]].visited = true ; //If the atoms are already in the path then count them as visit so we don't visit them again
		}
	}
	

	//First initialize vertices
	vertices.verts[sourceId].visited = false; //We want to be able to visit the starting atom so mark that as unvisited
	vertices.verts[sourceId].distance = 0; // The distance between the source atom and istelf is 0 
	Vertex minVertex = GetMinUnVisitedVertexId(vertices, xYZData); //We want to get the atom that is closest to the source atom  (this should be the source atom if it we are just starting)
	while(minVertex.atomId != endId) //If the next nearest atom is the end atom then we have found the shorest chain
	{
		//Loop through all the neighbors
		vertices.verts[minVertex.atomId].visited = true; //We are now finishing this atom
		for (int neighborIndex = 0; neighborIndex < vertices.verts[minVertex.atomId].nNeighbors; neighborIndex++) //Lets now go through all the neigbors of the atom we are visiting
		{
			if(vertices.verts[minVertex.neighbors[neighborIndex].atomId].visited == false) //If we have already visited this atom ignore it
			{
				if(vertices.verts[minVertex.neighbors[neighborIndex].atomId].distance > minVertex.distance + minVertex.neighbors[neighborIndex].distance) //If the distanceof the neighbor we have visited is greater than the distance from the visited atom to the source plus the distance between the visited and the nuehibor than replace the distance with this new shorter distance. Also change the previous atom to the visited one 
				{
					vertices.verts[minVertex.neighbors[neighborIndex].atomId].distance = minVertex.distance + minVertex.neighbors[neighborIndex].distance;
					vertices.verts[minVertex.neighbors[neighborIndex].atomId].prevId = minVertex.atomId;
					
				}
			}

		}
		minVertex = GetMinUnVisitedVertexId(vertices, xYZData); //Now we want to get the next shortest atom
		vertices.verts[endId].visited = false; //If the min vertex id is the false one then the code will cancel anyway. If the endId was visited in the first path and is close enough that it would be a possible neighbor of the source id then we will still want it to be visited until a new minVertex has been picked
		if(minVertex.atomId == -1) //There is no path
		{
			break; //End the while loop
		}
	}	
	Path path;
	int prevId = vertices.verts[endId].prevId ;
	path.pathIds[0] = endId;
	path.pathIds[1] = vertices.verts[endId].prevId;
	path.nPathAtoms = 2;
	if(minVertex.atomId == -1) //The minVertex atom Id is negative 1 so that means the path doesn't exist so to make it clear we have the source id show up again 
	{
		cout << "PATH WAS NOT FOUND" << endl;
		path.pathIds[1] = -1; //This means that we have the end id twice which means we do not have a real chain
        	path.pathIds[2] = sourceId;
		path.nPathAtoms = 3;
		return path;
	}
	while(prevId != sourceId)
	{
		prevId = vertices.verts[prevId].prevId;
		path.pathIds[path.nPathAtoms] = prevId;
		path.nPathAtoms = path.nPathAtoms + 1;
	}
	return path;

}

void GatherAllPaths(string outputDirectory, XYZData xYZData, Vertices originalVertices, int nChains)
{
        int npaths = 0;
        std::ofstream out(outputDirectory.c_str()); //+"output.txt");
        for (int atom1Index = 0; atom1Index < xYZData.nAtoms; atom1Index = atom1Index +1)
        {
                cout << atom1Index << endl;
		for (int atom2Index = 0; atom2Index < xYZData.nAtoms; atom2Index = atom2Index +1)
                {
                        if(atom1Index < atom2Index)
                        {
                                float doi = CalculateDistance(xYZData.xYZAtomData[atom1Index], xYZData.xYZAtomData[atom2Index],xYZData.xLength, xYZData.yLength, xYZData.zLength);
				if(doi < DOI_END and doi > DOI_BEG)
				{
                                	int binIndex = FindAtomPairs(doi , BIN_SIZE);
					Path totalPath;
					totalPath.nPathAtoms = 0;
					for (int chainIndex = 0; chainIndex <  nChains; chainIndex++)
					{
						Vertices vertices(originalVertices.verts);
						Path path = djikstra(vertices, atom1Index, atom2Index, xYZData, totalPath);
						if(path.nPathAtoms == 0)
						{
							cout << "NO ATOMS FOUND IN PATH" << endl;
						}
						if(path.nPathAtoms == 2)
                                                {
							if(path.pathIds[1] == path.pathIds[0])
							{
								cout << "END PATH EQUALS START PATH" << endl;
							}
						}
						MakeAndWritePath(outputDirectory, binIndex, path, xYZData,nChains,(chainIndex+1));
						totalPath.AppendPath(path);
					}
				}
                        }
                }
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	clock_t begin = clock();
	char xyzFile[100];
	if (argc != 4)
	{
		cout << "You have entered the wrong number of arguments" << endl;
	}
	const char * inputFile = argv[1];
	const char * outputDirectory = argv[2];
	const char * configFile = argv[3];
	cout << inputFile << endl;
	cout << outputDirectory << endl;
	cout << configFile << endl;
	Bonds bonds;	
	ConfigData configData;
	XYZData xYZData;
	configData = ReadInConfigFile(configFile);
	cout << configData.inputFile << endl;
	bonds = LoadBonds(configData.allowedBond, configData.cutOffDistances);
	bonds.GetInfo();
	xYZData = ReadInXYZFile(configData.inputFile.c_str());
	Vertices vertices = MakeGraph(xYZData, bonds); //Set up all the vertices that we can use for the djikstra algrotihm
	GatherAllPaths(outputDirectory,xYZData,vertices, numberOfChains);
	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << elapsed_secs << endl;
	return 0;
}



