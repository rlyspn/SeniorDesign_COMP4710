// writing on a text file
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "../../parallelMetropolis/src/metroParallelUtil.h"
using namespace std;


class Opls_Scan{
   private:
      map<int,Atom> oplsTable; //"HashTable" the holds all opls ref.
      string fileName;
   public:
      Opls_Scan(string filename); // constructor
      ~Opls_Scan();
      int scanInOpls(string filename); // scans in the oplFile
      void addLineToTable(string line); //adds entry to map
};





