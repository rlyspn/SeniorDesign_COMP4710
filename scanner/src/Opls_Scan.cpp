#include "opls_Scan.h"

//constructor
Opls_Scan::Opls_Scan(string filename){
   fileName = filename
   scanInOpls(filename);
}

Opls_Scan::~Opls_Scan(){
}

//scans in the opls file and adds entries to the table
int Opls_Scan::scanInOpls(string filename){
    int numOfLines=0;
   ifstrem oplsScanner(filename)a
   if( !oplsScanner.open() )
      return -1;
   else {
      string line; 
      while( oplsScanner.good() )
      {
        numOfLines++;
        getline(oplsScanner,line);

        //check if it is a commented line
        if(line.at(0) != '#')
            addLineToTable(line);
      }
      oplsScanner.close()
   }
}

// adds a string line into the table
void Opls_Scan::addLineToTable(string line){
}

