// writing on a text file
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include "metroUtil.h"
using namespace std;


class Opls_Scan{
   private:
      map<int,Atom> oplsTable; //"HashTable" the holds all opls ref.
      string fileName;
   public:
      Opls_Scan(string filename); // constructor
      ~Opls_Scan();
		
		/**
		Scans in the opls File calls sub-function addLineToTable
		@param filename - the name/path of the opls file
		*/
      int scanInOpls(string filename); 
		
		/**
		Parses out a line from the opls file and gets (sigma, epsiolon, charge)
		adds an entry to the map at the hash number position, 
		calls sub-function checkFormat()
		@param line - a line from the opls file
		@param numOflines - number of lines previously read in, used for error output
		*/
      void addLineToTable(string line, int numOfLines);
		
		/**
		Checks the format of the line being read in
		returns false if the format of the line is invalid
		@param line -  a line from the opls file
		*/
		bool checkFormat(string line);
		
		/**
		Returns an Atom struct based on the hashNum (1st col) in Z matrix file
		The Atom struct has -1 for x,y,z and has the hashNum for an id. 
		sigma, epsilon, charges
		@param hashNum -  the hash number (1st col) in Z matrix file
		*/
		Atom getAtom(int hashNum);
		
		/**
		Returns the sigma value based on the hashNum (1st col) in Z matrix filee
		@param hashNum -  the hash number (1st col) in Z matrix file
		*/
		double getSigma(int hashNum);
		
		/**
		Returns the epsilon value based on the hashNum (1st col) in Z matrix filee
		@param hashNum -  the hash number (1st col) in Z matrix file
		*/
		double getEpsilon(int hashNum);
		
		/**
		Returns the charge value based on the hashNum (1st col) in Z matrix filee
		@param hashNum -  the hash number (1st col) in Z matrix file
		*/
		double getCharge(int hashNum);
};





