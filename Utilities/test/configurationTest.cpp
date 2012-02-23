#include <assert.h>
#include "../src/Config_Scan.h"

void testConfigScan(){
    string configPath = "configurationTest.txt"; 
    string oplsPath = "path/to/opls/file";
    string zMatrixPath = "path/to/zMatrix/file";
    string stateInputPath = "path/to/state/input";
    string stateOutputPath = "path/to/State/output";
    string pdbOutputPath = "path/to/pdb/output";
    Config_Scan cs (configPath);
    cs.readInConfig();
   
    //test configuration path
    assert(configPath.compare(cs.getConfigPath()) == 0);
    //test opls path
    cout << "Opls path: " << cs.getOplsusaparPath() << endl;
    assert(oplsPath.compare(cs.getOplsusaparPath()) == 0);
    //test zmatrix path
    cout << "Zmatrix path: " << cs.getZmatrixPath() << endl;
    assert(zMatrixPath.compare(cs.getZmatrixPath()) == 0);
    //test state input path
    cout << "State input path: " << cs.getStatePath() << endl;
    assert(stateInputPath.compare(cs.getStatePath()) == 0);
    //test state ouput path
    cout << "State output path: " << cs.getStateOutputPath() << endl;
    assert(stateOutputPath.compare(cs.getStateOutputPath()) == 0);
    //test pdb output
    cout << "PDB output path: " << cs.getPdbOutputPath() << endl;
    assert(pdbOutputPath.compare(cs.getPdbOutputPath()) == 0);

    //test box dimensions


}

int main(){
    testConfigScan();
}
