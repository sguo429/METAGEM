/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/
#include "metagem.h"


void print_help();

// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {

    // GEM options. Details are printed from the print_help() function below.
    CLI::App app{""};

    // Defaults
    int metaOpt_in = 0;
    std::string outFile_in = "metagem.out";

    app.add_option("--input-files", fileNames, "")->expected(0, 1000000);
    app.add_option("--input-file-list", metaFileList, "")->expected(1);
    app.add_option("--exposure-names", intNames, "")->expected(1, 1000000)->required();
    app.add_option("--out", outFile_in, "")->expected(1);
    app.add_option("--meta-option", metaOpt_in, "");
    app.add_option("--additional-joint", additionalJointInfo, "") -> expected(0, 100);
    app.add_option("--additional-interaction", additionalInteractionInfo, "") -> expected(0, 100);
    app.add_option("--control-file", controlFilePath, "")->expected(1);
  
    try
    {
        app.parse( argc, argv);
      
        // Header-rename file
        size_t fhl = controlFilePath.length();
        if (fhl > 0) {
          renameHeaders = true;
        }

        if (renameHeaders) {
            if (fileNames.size() != 0 && metaFileList.length() != 0){
              cerr << "\nERROR: Only need to specify file names in the control file when changing file headers.\n\n";
            }
            auto result = loadHeaderRenaming(controlFilePath);
            fileNames = result.first;
            headerRenamings = result.second;
        }
      
        // Input files
        size_t fns = fileNames.size();
        size_t fls = metaFileList.length();
        if (fns == 0 && fls == 0) {
            cerr << "\nERROR: --input-files or --input-file-list is required.\n\n";
            exit(1);

        } else if (fns > 0 && fls > 0) {
            cerr << "\nERROR: Both --input-files and --input-file-list are specified.\n\n";
            exit(1);

        } else if (fns > 0) {
            if (fns < 2) {
                cerr << "\nERROR: METAGEM requires at least 2 input files.\n\n";
                exit(1);                
            }

            std::set<std::string> s(fileNames.begin(), fileNames.end());
            if (s.size() != fns) 
            {
                cerr << "\nERROR: There are duplicate input file names.\n\n";
                exit(1);
            }

        } else if (fls > 0) {
            std::ifstream file;
            std::string line;
            file.open(metaFileList);
            if (!file.is_open()) {
                cerr << "\nERROR: Cannot open the file: [" << metaFileList << "].\n\n";
                exit(1);
            }

            size_t nfiles = 0;
            while(getline(file, line)) {
                line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
                fileNames.push_back(line);
                nfiles++;
            }

            if (nfiles < 2)
            {
                cerr << "\nERROR: METAGEM requires at least 2 input files.\n\n";
                exit(1);
            }
        
            std::set<std::string> s(fileNames.begin(), fileNames.end());
            if (s.size() != nfiles) 
            {
                cerr << "\nERROR: There are duplicate input file names.\n\n";
                exit(1);
            }

        } else {
            cerr << "\nERROR: Unrecognized combination of --input-files and --input-file-list.\n\n";
            exit(1);
        }



        // Exposure names
        std::set<std::string> s(intNames.begin(), intNames.end());
        if (s.size() != intNames.size()) {
            cerr << "\nERROR: There are duplicate exposure names (--exposure-names).\n\n";
            exit(1);
        }
        nInt = intNames.size();

        lcIntNames = intNames;
        for(std::string &s : lcIntNames){
            std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
            s = "g-" + s;
        }
        lcIntNames.insert(lcIntNames.begin(), "g");



        // Output file
        outFile = outFile_in;
        std::ofstream results(outFile);
        if (!results) {
            printOpenFileError(outFile);
        }

        if (results.fail()) {
            printOpenFileError(outFile);
        }

        results << "test" << endl;
        if (results.fail()) {
            cerr << "\nERROR: Cannot write to output file.\n\n";
            results.close();
            
            if (std::remove(outFile.c_str()) != 0) {
                cerr << "\nERROR: Cannot delete output file.\n\n";
            }
            exit(1);
        }
        results.close();
        
        if (std::remove(outFile.c_str()) != 0) {
            cerr << "\nERROR: Cannot delete output file.\n\n";
            exit(1);
        }

        

        // Meta Option
        metaOpt = metaOpt_in;
        if (metaOpt < 0 || metaOpt > 2) {
            cerr << "\nERROR: The --meta-option integer value must be 0, 1, or 2.\n\n";
            exit(1);
        }

        if (metaOpt == 0) {
            mb = true;
            rb = true;
            cout << "Meta-analysis option: [model-based] and [robust].\n\n"; 
        } else if (metaOpt == 1) {
            mb = true;
            cout << "Meta-analysis option: [model-based].\n\n";
        } else if (metaOpt == 2) {
            rb = true;
            cout << "Meta-analysis option: [robust].\n\n";
        } else {
            cerr << "\nERROR: The --meta-option integer value must be 0, 1, or 2.\n\n";
            exit(1);
        }

        // Additional joint meta-analysis
        if (!additionalJointInfo.empty()) {
          additionalJoint = true;
          
          if (additionalJointInfo.size() > 1) {
            intNames2.assign(additionalJointInfo.begin(), additionalJointInfo.end() - 1);
          } else {
            intNames2.clear(); // Ensure it is empty if only the path is provided
          }

          outFile2 = additionalJointInfo.back();

          // Check if the file path was accidentally a variable name
          for (const auto& name : intNames) {
            if (outFile2 == name) {
              cerr << "ERROR: Please specify the full path of the additional output file at the end of '--additional-joint' flag.\n\n";
              exit(1);
            }
          }
          
          // Duplicate check (only if there are exposures to check)
          if (!intNames2.empty()) {
              std::set<std::string> s(intNames2.begin(), intNames2.end());
              if (s.size() != intNames2.size()) {
                  cerr << "\nERROR: There are duplicate exposure names in the additional joint meta-analysis.\n\n";
                  exit(1);
              }
          }
          
          nInt2 = intNames2.size() + 1;
          
          lcIntNames2 = intNames2;
          for(std::string &s : lcIntNames2){
              std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
              s = "g-" + s;
          }
          lcIntNames2.insert(lcIntNames2.begin(), "g");

          // Additional output file
          std::ofstream results2(outFile2);
          if (!results2) {
              printOpenFileError(outFile2);
          }

          if (results2.fail()) {
              printOpenFileError(outFile2);
          }

          results2 << "test" << endl;
          if (results2.fail()) {
              cerr << "\nERROR: Cannot write to the additional joint meta-analysis output file.\n\n";
              results2.close();
            
              if (std::remove(outFile2.c_str()) != 0) {
                  cerr << "\nERROR: Cannot delete the additional joint meta-analysis output file.\n\n";
              }
              exit(1);
          }
          results2.close();
        
          if (std::remove(outFile2.c_str()) != 0) {
              cerr << "\nERROR: Cannot delete the additional joint meta-analysis output file.\n\n";
              exit(1);
          }
        }

        // Additional interaction-only meta-analysis
        if (!additionalInteractionInfo.empty()) {
          additionalInteraction = true;

          outFile3 = additionalInteractionInfo.back();

          if (additionalInteractionInfo.size() > 1) {
              intNames3.assign(additionalInteractionInfo.begin(),
                               additionalInteractionInfo.end() - 1);
          } else {
              intNames3.clear();
              intNames3.push_back("G");   // genetic term only
          }

          // Check if the full path of the additional output file specified at the end
          const std::string& lastTestInfo = additionalInteractionInfo.back();
          for (const auto& name : intNames) {
            if (lastTestInfo == name) {
              cerr << "ERROR: Please specify the full path of the additional output file at the end of '--additional-interaction' flag.\n\n";
              exit(1);
            }
          }
          
          std::set<std::string> s(intNames3.begin(), intNames3.end());
          if (s.size() != intNames3.size()) {
              cerr << "\nERROR: There are duplicate names in --additional-interaction.\n\n";
              exit(1);
          }
      
          nInt3 = intNames3.size();
          if (nInt3 == 0) {
              cerr << "\nERROR: --additional-interaction produced 0 terms.\n\n";
              exit(1);
          }
      
          // Build lcIntNames3 for column matching in processFileHeader()
          lcIntNames3 = intNames3;
          for (auto &s : lcIntNames3) {
              std::transform(s.begin(), s.end(), s.begin(),
                             [](unsigned char c){ return std::tolower(c); });
      
              if (s == "g") {
                  // keep "g" for main effect column Beta_G
              } else {
                  s = "g-" + s;
              }
          }
          // Additional output file
          std::ofstream results3(outFile3);
          if (!results3) {
              printOpenFileError(outFile3);
          }

          if (results3.fail()) {
              printOpenFileError(outFile3);
          }

          results3 << "test" << endl;
          if (results3.fail()) {
              cerr << "\nERROR: Cannot write to the additional interaction-only meta-analysis output file.\n\n";
              results3.close();
            
              if (std::remove(outFile3.c_str()) != 0) {
                  cerr << "\nERROR: Cannot delete the additional interaction-only meta-analysis output file.\n\n";
              }
              exit(1);
          }
          results3.close();
        
          if (std::remove(outFile3.c_str()) != 0) {
              cerr << "\nERROR: Cannot delete the additional interaction-only meta-analysis output file.\n\n";
              exit(1);
          }
        }
        
    }
    catch( const CLI::CallForHelp &e )
    {
        print_help();
        exit(1);
    }
}





void print_help() {
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl;
    cout << endl << endl;

    cout << "Input File Options: " << endl
        << "   --input-files \t Output files from GEM 'meta' or 'full' option." << endl
        << "   --input-file-list \t A no header text file containing a single file name per line." << endl
        << "   --exposure-names \t The names of the exposure(s) to be included in the meta-analysis." << endl
        << "   --out \t\t Full path and extension to where METAGEM output results. \n \t\t\t    Default: metagem.out" << endl
        << "   --meta-option \t Integer value indicating which summary statistics should be used for meta-analysis. \n\t\t\t    0: Both model-based and robust summary statistics. \n \t\t\t    1: model-based summary statistics. \n \t\t\t    2: robust summary statistics. \n \t\t\t    Default: 0" << endl
        << "   --additional-joint \t The variable names and the full path of the output file for one additional joint meta-analysis." << endl
        << "   --additional-interaction \t The variable names and the full path of the output file for one additional interaction-only meta-analysis." << endl
        << "   --control-file \t A no header text file containing file names in seperate lines with a 'FILE' in front of the file name in each line, and containing both of the changed column name(s) and the original column name(s) following the line(s) of the file name(s) which need to do column name changing. This file should contain at least two file names." << endl;

    cout << endl << endl;
    cout << endl << endl;
}

std::pair<std::vector<std::string>, std::map<std::string, std::map<std::string, std::string>>> loadHeaderRenaming(std::string controlFilePath) {
    std::map<std::string, std::map<std::string, std::string>> fileColumnMappings;
    std::vector<std::string> filenames;

    std::ifstream titlesFile(controlFilePath);
    if (!titlesFile.is_open()) {
        std::cerr << "ERROR: Unable to open the control file: " << controlFilePath << std::endl;
        exit(1);
    }

    std::string line, currentFile;
    while (getline(titlesFile, line)) {
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) { return !std::isspace(ch); }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), line.end());

        if (line.empty()) continue;

        if (line.substr(0, 4) == "FILE") {
            currentFile = line.substr(5);
            filenames.push_back(currentFile);
            continue;
        }

        if (!currentFile.empty()) {
            std::istringstream iss(line);
            std::string newName, oldName;
            if (iss >> newName >> oldName) {
                fileColumnMappings[currentFile][oldName] = newName;
            } else {
                std::cerr << "ERROR: Invalid header renaming line: " << line << std::endl;
                exit(1);
            }
        } else {
            std::cerr << "ERROR: Header renaming found before any file declaration: " << line << std::endl;
            exit(1);
        }
    }
    titlesFile.close();

    return {filenames, fileColumnMappings};
}
