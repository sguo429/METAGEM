#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {
public:

    std::vector<std::string> fileNames;
    std::vector<std::string> intNames;
    std::vector<std::string> lcIntNames;

    std::string outFile;
    std::string metaFileList;
    std::vector<std::string> additionalTest;

    size_t nInt;
    size_t nCov;
    int metaOpt;

    bool mb = false;
    bool rb = false;

    void processCommandLine(int argc, char* argv[]);
};


#endif
