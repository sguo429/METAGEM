#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

class CommandLine {
public:
	// Input files
	std::vector<std::string> fileNames;
	std::vector<std::string> intNames;

	size_t nExp;
	size_t nCov;
	std::string outFile;
	std::string metaFileList;

	int metaOpt;
	bool mb = false;
	bool rb = false;

	void processCommandLine(int argc, char* argv[]);
};


#endif
