/*
 
 paramreader.cpp
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <stdlib.h>

#include "rta.h"

using namespace std;

// this workhorse examines a key to see if it corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(const char *key, const char *value) {
    // integer params
    if (strcmp(key,"num")==0) num=atoi(value);
    if (strcmp(key,"maxiters")==0) maxiters=atoi(value);
    if (strcmp(key,"update")==0) update=atoi(value);
    if (strcmp(key,"snapupdate")==0) snapupdate=atoi(value);
    if (strcmp(key,"numPZ")==0) numPZ=atoi(value);
    if (strcmp(key,"numPT")==0) numPT=atoi(value);
    if (strcmp(key,"computeM")==0) computeMoments=atoi(value);
    if (strcmp(key,"computeF")==0) computeDist=atoi(value);
    if (strcmp(key,"fStep")==0) fStep=atoi(value);
    if (strcmp(key,"maxN")==0) maxN=atoi(value);
    if (strcmp(key,"maxM")==0) maxM=atoi(value);
    // double/float params
    if (strcmp(key,"t0")==0) t0=atof(value);
    if (strcmp(key,"tf")==0) tf=atof(value);
    if (strcmp(key,"T0")==0) T0=atof(value);
    if (strcmp(key,"a0")==0) a0=atof(value);
    if (strcmp(key,"fpieb")==0) fpieb=atof(value);
    if (strcmp(key,"maxPZ")==0) maxPZ=atof(value);
    if (strcmp(key,"maxPT")==0) maxPT=atof(value);
    return;
}

//
// This routine assumes that parameters are in text file with
// each parameter on a new line in the format
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
//
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParametersFromFile(string filename, int echo) {
    
    string commentmarker = "//";
    char space = ' ';
    char tab = '\t';
    
    int maxline = 128; // maximum line length used in the buffer for reading
    char buffer[maxline];
    ifstream paramFile(filename.c_str());
    
    while(!paramFile.eof()) {
        paramFile.getline(buffer,maxline,'\n');
        string line = buffer; int length = strlen(buffer);
        if (line.substr(0,commentmarker.length())!=commentmarker && line.length()>0) {
            char key[64]="",value[64]="";
            int founddelim=0;
            for (int i=0;i<length;i++) {
                if (buffer[i]==space || buffer[i]==tab) founddelim=1;
                else {
                    if (founddelim==0) key[strlen(key)] = buffer[i];
                    else value[strlen(value)] = buffer[i];
                }
            }
            if (strlen(key)>0 && strlen(value)>0) {
                setParameter(key,value);
                if (echo) cout << key << " = " << value << endl;
            }
        }
    }
    
    return;
}

//
// Read parameters from commandline
//
void readParametersFromCommandLine(int argc, char** argv, int echo) {
    int optind = 1;
    while (optind < argc)
    {
        if (argv[optind][0]=='-') {
            string key = argv[optind];
            key = key.substr(1,key.length()); // remove '-'
            string value = argv[optind+1]; // load value
            if (echo) cout << key << " = " << value << endl;
            setParameter(key.c_str(),value.c_str());
            optind++;
        }
        optind++;
    }
    return;
}

// parameter processing prior to run
void processParameters() {
    //VMIN = 2*log(PMIN);
    //VMAX = 2*log(PMAX);
}
