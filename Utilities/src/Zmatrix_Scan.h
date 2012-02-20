// writing on a text file
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "metroUtil.h"
#include "Opls_Scan.h"
using namespace std;


class Zmatrix_Scan{
   private:
      string fileName;
      Opls_Scan *oplsScanner;
      vector<Molecule> moleculePattern;
      vector<Atom> atomVector;
      vector<Bond> bondVector;
      vector<Angle> angleVector;
      vector<Dihedral> dihedralVector;
      vector<int>  moleculeLimits;

      bool startNewMolecule;
   public:
      Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef); // constructor
      ~Zmatrix_Scan();
		
		/**
		Scans in the z-matrix File calls sub-function parseLine
		@param filename - the name/path of the z-matrix file
		*/
      int scanInZmatrix(); 
		
		/**
		Parses out a line from the zmatrix file and gets the atom from the OPLS hash
		Creates structs for bond, angle, dihedral, if applicable
		calls sub-function checkFormat()
		@param line - a line from the zmatrix file
		@param numOflines - number of lines previously read in, used for error output
		*/
      void parseLine(string line, int numOfLines);
		
		/**
		Checks the format of the line being read in
		returns false if the format of the line is invalid
		@param line -  a line from the zmatrix file
		*/
        bool checkFormat(string line);

        /**
        Creates a molecule(s)  based on a starting unique ID and the pattern specified
        by the Z-matrix in the scan functions
        @param startingID - first ID for the molecule being built
        */
        vector<Molecule> buildMolecule(int startingID);
};
