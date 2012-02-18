#include "Opls_Scan.h"

//constructor
Opls_Scan::Opls_Scan(string filename){
   fileName = filename;
   scanInOpls(filename);
}

Opls_Scan::~Opls_Scan(){
    oplsTable.clear();
}

//scans in the opls file and adds entries to the table
int Opls_Scan::scanInOpls(string filename){
    int numOfLines=0;
   ifstream oplsScanner(filename.c_str());
   if( !oplsScanner.is_open() )
      return -1;
   else {
      string line; 
      while( oplsScanner.good() )
      {
        numOfLines++;
        getline(oplsScanner,line);

        //check if it is a commented line,
		  //or if it is a title line
        if(line.at(0) != '#' && numOfLines > 1)
            addLineToTable(line,numOfLines);
      }
      oplsScanner.close();
   }
}

// adds a string line into the table
void Opls_Scan::addLineToTable(string line, int numOfLines){

    string hashNum;
    int secCol;
	 double charge,sigma,epsilon;
	 string name,extra;
	 
	 //check if line contains correct format
	 if(checkFormat(line)){
	     stringstream ss ; 
        ss << line;    	
        ss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon >> extra;
	 
	     Atom temp = createAtom(-1, -1, -1, -1, sigma, epsilon, charge);
	 
	     pair<map<string,Atom>::iterator,bool> ret;
	 	  ret = oplsTable.insert( pair<string,Atom>(hashNum,temp) );
	     if (ret.second==false)
        {
           cerr << "Err Opls Scanner: element "<< hashNum << "already existed" <<endl;
        }
	 }else
	     cerr << "Err Opls Scanner: line "<< numOfLines <<"contains bad format" << endl;

}

// check if line contains the right format...
bool Opls_Scan::checkFormat(string line){ 
    stringstream iss(line);
	 int count =0;
	 string word;
	 while (iss >> word)
	     count ++;
	 return (count > 5);
}

//return an Atom struct of the given hash value
Atom Opls_Scan::getAtom(string hashNum){
    if(oplsTable.count(hashNum)>0 ){
	     return oplsTable[hashNum];
	 }
	 else{
	    cerr << "Index does not exist" <<endl;
		 return createAtom(0, -1, -1, -1, -1, -1, -1);
	 }
}

//return the Sigma value associated with that hashNum/ atom number
double Opls_Scan::getSigma(string hashNum){
    if(oplsTable.count(hashNum)>0 ){
	     Atom temp = oplsTable[hashNum];
		  return temp.sigma;
	 }
	 else{
	    cerr << "Index does not exist" <<endl;
		return -1;
	 }
}

//return the Epsilon value associated with the  hashNum/ Atom number
double Opls_Scan::getEpsilon(string hashNum){
    if(oplsTable.count(hashNum)>0 ){
	     Atom temp = oplsTable[hashNum];
		  return temp.epsilon;
	 }
	 else{
	    cerr << "Index does not exist" <<endl;
		 return -1;
	 }
}

//return the Charge value associated with the  hashNum/ Atom number
double Opls_Scan::getCharge(string hashNum){
    if(oplsTable.count(hashNum)>0 ){
	     Atom temp = oplsTable[hashNum];
		  return temp.charge;
	 }
	 else{
	    cerr << "Index does not exist" <<endl;
		 return -1;
	 }
}
