#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {
public:

    std::vector<std::string> fileNames;
    std::vector<std::string> intNames;

    std::string outFile;
    std::string metaFileList;

    size_t nExp;
    size_t nCov;
    int metaOpt;

    bool mb = false;
    bool rb = false;

    void processCommandLine(int argc, char* argv[]);
};


#endif
