#ifndef CONFIG_SCAN_H
#define CONFIG_SCAN_H

#include <iostream>
#include <fstream>
#include "metroUtil.h"
#include <cstdlib>

using namespace std;

class Config_Scan{
    private:
        Environment enviro;
        string configpath;
        long numOfSteps;
        string oplsuaparPath;
        string zmatrixPath;
        string statePath;
        string stateOutputPath;
        string pdbOutputPath;

    public:

        /*
            Config file scanner instantiated with a path.
            @param configPath - path to paramaters.cfg
        */
        Config_Scan(string configPath);

        /*
            Reads in the config file located at config path
            given in teh constructor.
        */
        void readInConfig();

        /*
            getters
        */

        Environment getEnviro();

        string getConfigPath();

        long getSteps();
        
        string getOplsusaparPath();

        string getZmatrixPath();

        string getStatePath();

        string getStateOutputPath();

        string getPdbOutputPath();

};

#endif //CONFIG_SCAN_H
