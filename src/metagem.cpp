// <METAGEM: META-analysis of GEM summary statistics>
// Copyright (C) <2021-2025> Duy T. Pham, Han Chen and Shuyi Guo 
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.


#include "metagem.h"
#include "time.h"

int main(int argc, char* argv[]) 
{
    printWelcome();
    
    // Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    // Process command line arguments
    CommandLine cmd;
    cmd.processCommandLine(argc, argv);

    metagem(cmd);

    // Stop timers
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    printTimeCompleted(wall0, wall1, cpu0, cpu1);
    
    return(0);
}


void metagem(CommandLine cmd)
{
    // Read the header of each file
    printProcessingFiles();
    FileInfo* fip = new FileInfo();
    processFileHeader(cmd.nInt + 1, cmd.nInt2, cmd.nInt3, cmd.mb, cmd.rb, cmd.additionalJoint, cmd.additionalInteraction, cmd.renameHeaders, cmd.lcIntNames, cmd.lcIntNames2, cmd.lcIntNames3, cmd.fileNames, cmd.headerRenamings, fip);

    // Initialize objects
    int nvars = 0;
    bool mb = cmd.mb;
    bool rb = cmd.rb;
    bool additionalJoint = cmd.additionalJoint;
    bool additionalInteraction = cmd.additionalInteraction;
    size_t nInt   = cmd.nInt;
    size_t nInt2   = cmd.nInt2;
    size_t nInt3   = cmd.nInt3;
    size_t nInt1  = nInt + 1; // plus G
    size_t nFiles = cmd.fileNames.size();

    int nInt1_sq = nInt1 * nInt1;
    int nInt2_sq = nInt2 * nInt2;
    int nInt3_sq = nInt3 * nInt3;
    std::vector<std::string> snpid;
    std::vector<std::string> effectAllele;
    std::vector<std::string> nonEffectAllele;
    std::vector<double> nSamples;
    std::vector<double> AF;
    std::vector<size_t> nSeen;

    std::vector<double> mb_MU;
    std::vector<double> mb_MV;
    std::vector<double> mb_U;
    std::vector<double> mb_V;        
    std::vector<double> rb_MU;
    std::vector<double> rb_MV;
    std::vector<double> rb_U;
    std::vector<double> rb_V;
    std::vector<double> mb_U2;
    std::vector<double> mb_V2;
    std::vector<double> rb_U2;
    std::vector<double> rb_V2;
    std::vector<double> mb_U3;
    std::vector<double> mb_V3;
    std::vector<double> rb_U3;
    std::vector<double> rb_V3;
    sparse_hash_map<std::pair<std::string, std::string>, int> snpid_idx;
    
    for (size_t f = 0; f < nFiles; f++)
    {
        // Initialize
        std::string fileName = fip->fileNames[f];
        int snpColumn        = fip->snpColumn[fileName];
        int chrColumn        = fip->chrColumn[fileName];
        int posColumn        = fip->posColumn[fileName];
        int effectColumn     = fip->effectColumn[fileName];
        int nonEffectColumn  = fip->nonEffectColumn[fileName];
        int nSampleColumn    = fip->nSampleColumn[fileName];
        int afColumn         = fip->freqColumn[fileName];
        int betaMargColumn   = fip->betaMargColumn[fileName];	
        int nheader          = fip->nheader[fileName];
        std::vector<int> betaIntColumn = fip->betaIntColumn[fileName];
        std::vector<int> betaIntColumn2 = fip->betaIntColumn2[fileName];
        std::vector<int> betaIntColumn3 = fip->betaIntColumn3[fileName];

        // Read input file
        std::ifstream file;
        file.open(fileName);
        if (!file.is_open()) {
            printOpenFileError(fileName);
            exit(1);
        }

        // Ignore comments in the text file (ex. #dispersion: )
        std::string line;
        getline(file, line);

        // Read in summary statistics
        int n = 0;
        std::string value;
        std::vector <std::string> values(nheader);
        while(line.rfind("#", 0) == 0) { getline(file, line); }
    
        
        std::vector<int> seen(nvars, 0);
        while(getline(file, line))
        {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            
            // Get values
            int z = 0;
            std::istringstream iss(line);
            while (getline(iss, value, '\t')) {
                values[z] = value;
                z++;
            }

            if (z != nheader) {
                printNColumnError(z, n + 1, fileName, nheader);
                exit(1);
            }

            std::string snpName = values[snpColumn];
            std::string chr = values[chrColumn];
            std::string pos = values[posColumn];
            std::string a1  = values[nonEffectColumn];
            std::string a2  = values[effectColumn];
            double file_ns  = std::stod(values[nSampleColumn]);
            double file_af  = std::stod(values[afColumn]);
            std::transform(a1.begin(), a1.end(), a1.begin(), [](char c){ return std::tolower(c); });
            std::transform(a2.begin(), a2.end(), a2.begin(), [](char c){ return std::tolower(c); });

            
            std::string identifier = chr + ":" + pos;
            std::string idA = identifier + ":" + a1 + ":" + a2;
            std::string idB = identifier + ":" + a2 + ":" + a1;
            std::pair<std::string, std::string> id = {idA, idB};
            sparse_hash_map<std::pair<std::string, std::string>, int>::iterator it = snpid_idx.find(id);
            if (it == snpid_idx.end()) {
                snpid_idx[id] = nvars;
                seen.push_back(1);

                snpid.push_back(snpName + "\t" + chr+ "\t" + pos + "\t" + a1 + "\t" + a2 + "\t");
                nonEffectAllele.push_back(a1);
                effectAllele.push_back(a2);
                AF.push_back(file_af * 2.0 * file_ns);
                nSamples.push_back(file_ns);
                nSeen.push_back(1);
                
                // Matrix indices
                int ui = nvars * nInt1;
                int vi = nvars * nInt1_sq;

                // Betas
                double fileBM = std::stod(values[betaMargColumn]);
                std::vector<double> fileBetaInt(nInt1);
                for (size_t i = 0; i < nInt1; i++) {
                    fileBetaInt[i] = std::stod(values[betaIntColumn[i]]);
                }

                // Model-based
                if (mb)
                {
                    // Marginal
                    int mb_seMargColumn = fip->mb_seMargColumn[fileName];
                    double mb_fileMV    = std::stod(values[mb_seMargColumn]);
                    double mb_margV     = 1 / (mb_fileMV * mb_fileMV);
                    mb_MV.push_back(1 / (mb_fileMV * mb_fileMV));
                    mb_MU.push_back(mb_margV * (fileBM));

                    // Interaction
                    mb_U.resize(mb_U.size() + nInt1);
                    std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t isq = (i * nInt1);
                        size_t ii  = vi + isq;
                        for (size_t j = 0; j < nInt1; j++) {
                            mb_V.push_back(std::stod(values[mb_covIntColumn[isq + j]]));
                        }
                        mb_V[ii + i] *= mb_V[ii + i];
                    }

                    // Add to main matrix
                    subMatInv(&mb_V[0], nInt1, vi);
                    subMatsubVecprod(&mb_V[0], &fileBetaInt[0], &mb_U[0], nInt1, nInt1, vi, 0, ui);
                }

                // Robust
                if (rb)
                {
                    // Marginal
                    int rb_seMargColumn = fip->rb_seMargColumn[fileName];
                    double rb_fileMV    = std::stod(values[rb_seMargColumn]);
                    double rb_margV     = 1 / (rb_fileMV * rb_fileMV);
                    rb_MV.push_back(1 / (rb_fileMV * rb_fileMV));
                    rb_MU.push_back(rb_margV * (fileBM));

                    // Interaction
                    rb_U.resize(rb_U.size() + nInt1);
                    std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t isq = (i * nInt1);
                        size_t ii  = vi + isq;
                        for (size_t j = 0; j < nInt1; j++) {
                            rb_V.push_back(std::stod(values[rb_covIntColumn[isq + j]]));
                        }
                        rb_V[ii + i] *= rb_V[ii + i];
                    }

                    // Add to main matrix
                    subMatInv(&rb_V[0], nInt1, vi);
                    subMatsubVecprod(&rb_V[0], &fileBetaInt[0], &rb_U[0], nInt1, nInt1, vi, 0, ui);
                }

                // Additional joint meta-analysis
                if(additionalJoint){
                    // Matrix indices
                    int ui = nvars * nInt2;
                    int vi = nvars * nInt2_sq;
    
                    // Betas
                    std::vector<double> fileBetaInt2(nInt2);
                    for (size_t i = 0; i < nInt2; i++) {
                        fileBetaInt2[i] = std::stod(values[betaIntColumn2[i]]);
                    }
    
                    // Model-based
                    if (mb)
                    {
                        // Interaction
                        mb_U2.resize(mb_U2.size() + nInt2);
                        std::vector<int> mb_covIntColumn2 = fip->mb_covIntColumn2[fileName];
                        for (size_t i = 0; i < nInt2; i++) {
                            size_t isq = (i * nInt2);
                            size_t ii  = vi + isq;
                            for (size_t j = 0; j < nInt2; j++) {
                                mb_V2.push_back(std::stod(values[mb_covIntColumn2[isq + j]]));
                            }
                            mb_V2[ii + i] *= mb_V2[ii + i];
                        }
    
                        // Add to main matrix
                        subMatInv(&mb_V2[0], nInt2, vi);
                        subMatsubVecprod(&mb_V2[0], &fileBetaInt2[0], &mb_U2[0], nInt2, nInt2, vi, 0, ui);
                    }
    
                    // Robust
                    if (rb)
                    {
                        // Interaction
                        rb_U2.resize(rb_U2.size() + nInt2);
                        std::vector<int> rb_covIntColumn2 = fip->rb_covIntColumn2[fileName];
                        for (size_t i = 0; i < nInt2; i++) {
                            size_t isq = (i * nInt2);
                            size_t ii  = vi + isq;
                            for (size_t j = 0; j < nInt2; j++) {
                                rb_V2.push_back(std::stod(values[rb_covIntColumn2[isq + j]]));
                            }
                            rb_V2[ii + i] *= rb_V2[ii + i];
                        }
    
                        // Add to main matrix
                        subMatInv(&rb_V2[0], nInt2, vi);
                        subMatsubVecprod(&rb_V2[0], &fileBetaInt2[0], &rb_U2[0], nInt2, nInt2, vi, 0, ui);
                    }
                }

                // Additional Interaction-only meta-analysis
                if(additionalInteraction){
                    // Matrix indices
                    int ui = nvars * nInt3;
                    int vi = nvars * nInt3_sq;
    
                    // Betas
                    std::vector<double> fileBetaInt3(nInt3);
                    for (size_t i = 0; i < nInt3; i++) {
                        fileBetaInt3[i] = std::stod(values[betaIntColumn3[i]]);
                    }
    
                    // Model-based
                    if (mb)
                    {
                        // Interaction
                        mb_U3.resize(mb_U3.size() + nInt3);
                        std::vector<int> mb_covIntColumn3 = fip->mb_covIntColumn3[fileName];
                        for (size_t i = 0; i < nInt3; i++) {
                            size_t isq = (i * nInt3);
                            size_t ii  = vi + isq;
                            for (size_t j = 0; j < nInt3; j++) {
                                mb_V3.push_back(std::stod(values[mb_covIntColumn3[isq + j]]));
                            }
                            mb_V3[ii + i] *= mb_V3[ii + i];
                        }
    
                        // Add to main matrix
                        subMatInv(&mb_V3[0], nInt3, vi);
                        subMatsubVecprod(&mb_V3[0], &fileBetaInt3[0], &mb_U3[0], nInt3, nInt3, vi, 0, ui);
                    }
    
                    // Robust
                    if (rb)
                    {
                        // Interaction
                        rb_U3.resize(rb_U3.size() + nInt3);
                        std::vector<int> rb_covIntColumn3 = fip->rb_covIntColumn3[fileName];
                        for (size_t i = 0; i < nInt3; i++) {
                            size_t isq = (i * nInt3);
                            size_t ii  = vi + isq;
                            for (size_t j = 0; j < nInt3; j++) {
                                rb_V3.push_back(std::stod(values[rb_covIntColumn3[isq + j]]));
                            }
                            rb_V3[ii + i] *= rb_V3[ii + i];
                        }
    
                        // Add to main matrix
                        subMatInv(&rb_V3[0], nInt3, vi);
                        subMatsubVecprod(&rb_V3[0], &fileBetaInt3[0], &rb_U3[0], nInt3, nInt3, vi, 0, ui);
                    }
                }
                nvars++;
            } else {
                
                int index = it->second;
                if (seen[index] == 0) {
                    seen[index] = 1;
                } else {
                    n++;
                    continue;
                }

                int dir = 1.0;
                std::string na = nonEffectAllele[index];
                std::string ea = effectAllele[index];
                if ((a1.compare(na) == 0) && (a2.compare(ea) == 0)) {
                    AF[index] += file_af * 2.0 * file_ns;

                } else if ((a1.compare(ea) == 0) && (a2.compare(na) == 0)){
                    dir = -1.0;
                    AF[index] += (1 - file_af) * 2.0 * file_ns;
                }

                nSeen[index] += 1;
                nSamples[index] += file_ns;


                // Matrix index
                int ui = index * nInt1;
                int vi = index * (nInt1 * nInt1);

                // Betas
                double fileBM = std::stod(values[betaMargColumn]) * dir;
                std::vector<double> fileBetaInt(nInt1);
                for (size_t i = 0; i < nInt1; i++) {
                    fileBetaInt[i] = std::stod(values[betaIntColumn[i]]) * dir;
                }

                // Model-based
                if (mb)
                {
                    // Marginal
                    int mb_seMargColumn = fip->mb_seMargColumn[fileName];
                    double mb_fileMV = std::stod(values[mb_seMargColumn]);

                    mb_fileMV = 1 / (mb_fileMV * mb_fileMV);
                    mb_MV[index] += mb_fileMV;
                    mb_MU[index] += (mb_fileMV * (fileBM));

                    // Interactions
                    std::vector<double> mb_fileU(nInt1);
                    std::vector<double> mb_fileV(nInt1 * nInt1);
                    std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];

                    for (size_t i = 0; i < nInt1; i++) {
                        size_t ii = i * nInt1;
                        for (size_t j = 0; j < nInt1; j++) {
                            mb_fileV[ii + j] = std::stod(values[mb_covIntColumn[ii + j]]);
                        }
                        mb_fileV[ii + i] *= mb_fileV[ii + i];
                    }

                    // Add to main matrix
                    matInv(&mb_fileV[0], nInt1);
                    matvecprod(&mb_fileV[0], &fileBetaInt[0], &mb_fileU[0], nInt1, nInt1);
                    for (size_t i = 0; i < nInt1; i++){
                        size_t ii = i * nInt1;
                        mb_U[ui + i] += mb_fileU[i];
                        for (size_t j = 0; j < nInt1; j++) {
                            mb_V[vi + ii + j] += mb_fileV[ii + j];
                        }
                    }
                }

                // Robust
                if (rb)
                {
                    // Marginal
                    int rb_seMargColumn = fip->rb_seMargColumn[fileName];
                    double rb_fileMV = std::stod(values[rb_seMargColumn]);

                    rb_fileMV = 1 / (rb_fileMV * rb_fileMV);
                    rb_MV[index] += rb_fileMV;
                    rb_MU[index] += (rb_fileMV * (fileBM));

                    // Interactions
                    std::vector<double> rb_fileU(nInt1);
                    std::vector<double> rb_fileV(nInt1 * nInt1);
                    std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];
                    
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t ii = i * nInt1;
                        for (size_t j = 0; j < nInt1; j++) {
                            rb_fileV[ii + j] = std::stod(values[rb_covIntColumn[ii + j]]);
                        }
                        rb_fileV[ii + i] *= rb_fileV[ii + i];
                    }

                    // Add to main matrix
                    matInv(&rb_fileV[0], nInt1);
                    matvecprod(&rb_fileV[0], &fileBetaInt[0], &rb_fileU[0], nInt1, nInt1);
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t ii = i * nInt1;
                        rb_U[ui + i] += rb_fileU[i];
                        for (size_t j = 0; j < nInt1; j++) {
                            rb_V[vi + ii + j] += rb_fileV[ii + j];
                        }
                    }
                }

                // Additional joint meta-analysis
                if(additionalJoint){
                    // Matrix index
                    int ui = index * nInt2;
                    int vi = index * (nInt2 * nInt2);
    
                    // Betas
                    std::vector<double> fileBetaInt2(nInt2);
                    for (size_t i = 0; i < nInt2; i++) {
                        fileBetaInt2[i] = std::stod(values[betaIntColumn2[i]]) * dir;
                    }
    
                    // Model-based
                    if (mb)
                    {
                        // Interactions
                        std::vector<double> mb_fileU2(nInt2);
                        std::vector<double> mb_fileV2(nInt2 * nInt2);
                        std::vector<int> mb_covIntColumn2 = fip->mb_covIntColumn2[fileName];
    
                        for (size_t i = 0; i < nInt2; i++) {
                            size_t ii = i * nInt2;
                            for (size_t j = 0; j < nInt2; j++) {
                                mb_fileV2[ii + j] = std::stod(values[mb_covIntColumn2[ii + j]]);
                            }
                            mb_fileV2[ii + i] *= mb_fileV2[ii + i];
                        }
    
                        // Add to main matrix
                        matInv(&mb_fileV2[0], nInt2);
                        matvecprod(&mb_fileV2[0], &fileBetaInt2[0], &mb_fileU2[0], nInt2, nInt2);
                        for (size_t i = 0; i < nInt2; i++){
                            size_t ii = i * nInt2;
                            mb_U2[ui + i] += mb_fileU2[i];
                            for (size_t j = 0; j < nInt2; j++) {
                                mb_V2[vi + ii + j] += mb_fileV2[ii + j];
                            }
                        }
                    }
    
                    // Robust
                    if (rb)
                    {
                        // Interactions
                        std::vector<double> rb_fileU2(nInt2);
                        std::vector<double> rb_fileV2(nInt2 * nInt2);
                        std::vector<int> rb_covIntColumn2 = fip->rb_covIntColumn2[fileName];
                        
                        for (size_t i = 0; i < nInt2; i++) {
                            size_t ii = i * nInt2;
                            for (size_t j = 0; j < nInt2; j++) {
                                rb_fileV2[ii + j] = std::stod(values[rb_covIntColumn2[ii + j]]);
                            }
                            rb_fileV2[ii + i] *= rb_fileV2[ii + i];
                        }
    
                        // Add to main matrix
                        matInv(&rb_fileV2[0], nInt2);
                        matvecprod(&rb_fileV2[0], &fileBetaInt2[0], &rb_fileU2[0], nInt2, nInt2);
                        for (size_t i = 0; i < nInt2; i++) {
                            size_t ii = i * nInt2;
                            rb_U2[ui + i] += rb_fileU2[i];
                            for (size_t j = 0; j < nInt2; j++) {
                                rb_V2[vi + ii + j] += rb_fileV2[ii + j];
                            }
                        }
                    }
                }

                // Additional interaction-only meta-analysis
                if(additionalInteraction){
                    // Matrix index
                    int ui = index * nInt3;
                    int vi = index * (nInt3 * nInt3);
    
                    // Betas
                    std::vector<double> fileBetaInt3(nInt3);
                    for (size_t i = 0; i < nInt3; i++) {
                        fileBetaInt3[i] = std::stod(values[betaIntColumn3[i]]) * dir;
                    }
    
                    // Model-based
                    if (mb)
                    {
                        // Interactions
                        std::vector<double> mb_fileU3(nInt3);
                        std::vector<double> mb_fileV3(nInt3 * nInt3);
                        std::vector<int> mb_covIntColumn3 = fip->mb_covIntColumn3[fileName];
    
                        for (size_t i = 0; i < nInt3; i++) {
                            size_t ii = i * nInt3;
                            for (size_t j = 0; j < nInt3; j++) {
                                mb_fileV3[ii + j] = std::stod(values[mb_covIntColumn3[ii + j]]);
                            }
                            mb_fileV3[ii + i] *= mb_fileV3[ii + i];
                        }
    
                        // Add to main matrix
                        matInv(&mb_fileV3[0], nInt3);
                        matvecprod(&mb_fileV3[0], &fileBetaInt3[0], &mb_fileU3[0], nInt3, nInt3);
                        for (size_t i = 0; i < nInt3; i++){
                            size_t ii = i * nInt3;
                            mb_U3[ui + i] += mb_fileU3[i];
                            for (size_t j = 0; j < nInt3; j++) {
                                mb_V3[vi + ii + j] += mb_fileV3[ii + j];
                            }
                        }
                    }
    
                    // Robust
                    if (rb)
                    {
                        // Interactions
                        std::vector<double> rb_fileU3(nInt3);
                        std::vector<double> rb_fileV3(nInt3 * nInt3);
                        std::vector<int> rb_covIntColumn3 = fip->rb_covIntColumn3[fileName];
                        
                        for (size_t i = 0; i < nInt3; i++) {
                            size_t ii = i * nInt3;
                            for (size_t j = 0; j < nInt3; j++) {
                                rb_fileV3[ii + j] = std::stod(values[rb_covIntColumn3[ii + j]]);
                            }
                            rb_fileV3[ii + i] *= rb_fileV3[ii + i];
                        }
    
                        // Add to main matrix
                        matInv(&rb_fileV3[0], nInt3);
                        matvecprod(&rb_fileV3[0], &fileBetaInt3[0], &rb_fileU3[0], nInt3, nInt3);
                        for (size_t i = 0; i < nInt3; i++) {
                            size_t ii = i * nInt3;
                            rb_U3[ui + i] += rb_fileU3[i];
                            for (size_t j = 0; j < nInt3; j++) {
                                rb_V3[vi + ii + j] += rb_fileV3[ii + j];
                            }
                        }
                    }
                }
            }
            n++;
        }
        file.close();

        if (n == 0) {
            printZeroVariantsError(fileName);
            exit(1);
        }

        printProcessedFiles(fileName, n);
    }
    printDone(2);

    // Create the output file and write column header names
    printOutputHeader(mb, rb, additionalJoint, additionalInteraction, cmd.outFile, cmd.outFile2, cmd.outFile3, nInt1, nInt2, nInt3, cmd.intNames, cmd.intNames2, cmd.intNames3);
    std::ofstream results(cmd.outFile, std::ios_base::app);
    std::ofstream results2(cmd.outFile2, std::ios_base::app);
    std::ofstream results3(cmd.outFile3, std::ios_base::app);
    std::ostringstream oss;
    std::ostringstream oss2;
    std::ostringstream oss3;




    double vMarg;
    double betaMarg;
    double varMarg;
    double statMarg;
    double statInt;
    double statJoint;
    double pvalMarg;
    double pvalInt;
    double pvalJoint;
    double* Ai = new double[nInt1 * nInt1];
    double* VE = new double[nInt * nInt];
    double* Ai2 = new double[nInt2 * nInt2];
    double* VE2 = new double[(nInt2-1) * (nInt2-1)];
    double* Ai3 = new double[nInt3 * nInt3];
    std::vector<double> StempE(nInt, 0.0);
    std::vector<double> StempGE(nInt1, 0.0);
    std::vector<double> betaInt(nInt1, 0.0);
    std::vector<double> StempE2(nInt2-1, 0.0);
    std::vector<double> StempGE2(nInt2, 0.0);
    std::vector<double> betaInt2(nInt2, 0.0);
    std::vector<double> StempE3(nInt3, 0.0);
    std::vector<double> betaInt3(nInt3, 0.0);
    boost::math::chi_squared chisq_dist_M(1);
    boost::math::chi_squared chisq_dist_Int(nInt);
    boost::math::chi_squared chisq_dist_Joint(nInt1);
    boost::math::chi_squared chisq_dist_Int2(nInt2-1);
    boost::math::chi_squared chisq_dist_Joint2(nInt2);
    boost::math::chi_squared chisq_dist_Int3(std::max((int)nInt3, 1));

    printMetaBegin(nFiles, nvars);
    for (int i = 0; i < nvars; i++)
    {
        int is  = (i * nInt1);
        int iss = (i * nInt1 * nInt1);

        oss << snpid[i] << nSeen[i] << "\t" << nSamples[i] << "\t" << AF[i] / 2.0 / nSamples[i] << "\t";
        oss2 << snpid[i] << nSeen[i] << "\t" << nSamples[i] << "\t" << AF[i] / 2.0 / nSamples[i] << "\t";
        oss3 << snpid[i] << nSeen[i] << "\t" << nSamples[i] << "\t" << AF[i] / 2.0 / nSamples[i] << "\t";
        if (mb)
        {
            subMatrix(&mb_V[0], Ai, nInt1, nInt1, nInt1, nInt1, iss);
            subMatInv(&mb_V[0], nInt1, iss);
            subMatrix(&mb_V[0], VE, nInt, nInt, nInt1, nInt, iss + nInt1 + 1);

            // Model-based Marginal Variance
            vMarg = mb_MV[i];
            varMarg  = 1 / vMarg;
            betaMarg = mb_MU[i] * varMarg;
            statMarg = (betaMarg * betaMarg) * vMarg;
            pvalMarg = (std::isnan(statMarg) || statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, statMarg));


            // Interaction effects
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    betaInt[j] += (mb_V[iss + (nInt1 * j) + k] * mb_U[is + k]);
                }
            }


            // Int P-value
            matInv(VE, nInt);
            for (size_t j = 0; j < nInt; j++) {
                for (size_t k = 0; k < nInt; k++) {
                    StempE[j] += (VE[(nInt * j) + k] * betaInt[k + 1]);
                }
            }
            
            statInt = 0.0;
            for (size_t j = 1; j < nInt1; j++) 
                statInt += betaInt[j] * StempE[j-1];
            pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));
    

            // Joint P-value
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    StempGE[j] += (Ai[(nInt1 * j) + k] * betaInt[k]);
                }
            }

            statJoint = 0.0;
            for (size_t k = 0; k < nInt1; k++)
                statJoint += betaInt[k] * StempGE[k];
            pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));


            // Print
            oss << betaMarg << "\t" << sqrt(varMarg) << "\t";
            for (size_t j = 0; j < nInt1; j++) {
                oss << betaInt[j] << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                oss << sqrt(mb_V[iss + (ii * nInt1) + ii]) << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                for (size_t jj = 0; jj < nInt1; jj++) {
                    if (ii < jj) {
                        oss << mb_V[iss + (ii * nInt1) + jj] << "\t";
                    }
                }
            }
            oss << pvalMarg << "\t" << pvalInt << "\t" << pvalJoint << ((rb) ? "\t" : "\n");

            std::fill(StempE.begin(), StempE.end(), 0.0);
            std::fill(StempGE.begin(), StempGE.end(), 0.0);
            std::fill(betaInt.begin(), betaInt.end(), 0.0);
        }

        if (rb)
        {
            subMatrix(&rb_V[0], Ai, nInt1, nInt1, nInt1, nInt1, iss);
            subMatInv(&rb_V[0], nInt1, iss);
            subMatrix(&rb_V[0], VE, nInt, nInt, nInt1, nInt, iss + nInt1 + 1);

            // Robust Marginal Variance
            vMarg = rb_MV[i];
            varMarg  = 1 / vMarg;
            betaMarg = rb_MU[i] * varMarg;
            statMarg = (betaMarg * betaMarg) * vMarg;
            pvalMarg = (std::isnan(statMarg) || statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, statMarg));


            // Interaction effects
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    betaInt[j] += (rb_V[iss + (nInt1 * j) + k] * rb_U[is + k]);
                }
            }


            // Int P-value
            matInv(VE, nInt);
            for (size_t j = 0; j < nInt; j++) {
                for (size_t k = 0; k < nInt; k++) {
                    StempE[j] += (VE[(nInt * j) + k] * betaInt[k + 1]);
                }
            }

            statInt = 0.0;
            for (size_t j = 1; j < nInt1; j++) 
                statInt += betaInt[j] * StempE[j-1];
            pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));
    

            // Joint P-value
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    StempGE[j] += (Ai[(nInt1 * j) + k] * betaInt[k]);
                }
            }

            statJoint = 0.0;
            for (size_t k = 0; k < nInt1; k++)
                statJoint += betaInt[k] * StempGE[k];
            pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));


            // Print
            oss << betaMarg << "\t" << sqrt(varMarg) << "\t";
            for (size_t j = 0; j < nInt1; j++) {
                oss << betaInt[j] << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                oss << sqrt(rb_V[iss + (ii * nInt1) + ii]) << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                for (size_t jj = 0; jj < nInt1; jj++) {
                    if (ii < jj) {
                        oss << rb_V[iss + (ii * nInt1) + jj] << "\t";
                    }
                }
            }
            oss << pvalMarg << "\t" << pvalInt << "\t" << pvalJoint << "\n";

            std::fill(StempE.begin(), StempE.end(), 0.0);
            std::fill(StempGE.begin(), StempGE.end(), 0.0);
            std::fill(betaInt.begin(), betaInt.end(), 0.0);
        }

        if (i % 100000 == 0)
        {
            results << oss.str();
            oss.str(std::string());
            oss.clear();
        }

        if(additionalJoint){
            int is  = (i * nInt2);
            int iss = (i * nInt2 * nInt2);
    
            if (mb)
            {
                subMatrix(&mb_V2[0], Ai2, nInt2, nInt2, nInt2, nInt2, iss);
                subMatInv(&mb_V2[0], nInt2, iss);
                subMatrix(&mb_V2[0], VE2, nInt2-1, nInt2-1, nInt2, nInt2-1, iss + nInt2 + 1);
    
                // Interaction effects
                for (size_t j = 0; j < nInt2; j++) {
                    for (size_t k = 0; k < nInt2; k++) {
                        betaInt2[j] += (mb_V2[iss + (nInt2 * j) + k] * mb_U2[is + k]);
                    }
                }
    
                // Int P-value
                matInv(VE2, nInt2-1);
                for (size_t j = 0; j < (nInt2 - 1); j++) {
                    for (size_t k = 0; k < (nInt2 - 1); k++) {
                        StempE2[j] += (VE2[((nInt2 - 1) * j) + k] * betaInt2[k + 1]);
                    }
                }
            
                statInt = 0.0;
                for (size_t j = 1; j < nInt2; j++) 
                    statInt += betaInt2[j] * StempE2[j-1];
                pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int2, statInt));
    
    
                // Joint P-value
                for (size_t j = 0; j < nInt2; j++) {
                    for (size_t k = 0; k < nInt2; k++) {
                        StempGE2[j] += (Ai2[(nInt2 * j) + k] * betaInt2[k]);
                    }
                }

                statJoint = 0.0;
                for (size_t k = 0; k < nInt2; k++)
                    statJoint += betaInt2[k] * StempGE2[k];
                pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint2, statJoint));
                
  
                
                // Print
                for (size_t j = 0; j < nInt2; j++) {
                    oss2 << betaInt2[j] << "\t";
                }
                for (size_t ii = 0; ii < nInt2; ii++) {
                    oss2 << sqrt(mb_V2[iss + (ii * nInt2) + ii]) << "\t";
                }
                for (size_t ii = 0; ii < nInt2; ii++) {
                    for (size_t jj = 0; jj < nInt2; jj++) {
                        if (ii < jj) {
                            oss2 << mb_V2[iss + (ii * nInt2) + jj] << "\t";
                        }
                    }
                }
                oss2 << pvalInt << "\t" << pvalJoint << ((rb) ? "\t" : "\n");
                
                std::fill(StempE2.begin(), StempE2.end(), 0.0);
                std::fill(StempGE2.begin(), StempGE2.end(), 0.0);
                std::fill(betaInt2.begin(), betaInt2.end(), 0.0);
            }
    
            if (rb)
            {
                subMatrix(&rb_V2[0], Ai2, nInt2, nInt2, nInt2, nInt2, iss);
                subMatInv(&rb_V2[0], nInt2, iss);
                subMatrix(&rb_V2[0], VE2, nInt2-1, nInt2-1, nInt2, nInt2-1, iss + nInt2 + 1);
    
                // Interaction effects
                for (size_t j = 0; j < nInt2; j++) {
                    for (size_t k = 0; k < nInt2; k++) {
                        betaInt2[j] += (rb_V2[iss + (nInt2 * j) + k] * rb_U2[is + k]);
                    }
                }
    
                // Int P-value
                matInv(VE2, nInt2-1);
                for (size_t j = 0; j < (nInt2 - 1); j++) {
                    for (size_t k = 0; k < (nInt2 - 1); k++) {
                        StempE2[j] += (VE2[((nInt2 - 1) * j) + k] * betaInt2[k + 1]);
                    }
                }
            
                statInt = 0.0;
                for (size_t j = 1; j < nInt2; j++) 
                    statInt += betaInt2[j] * StempE2[j-1];
                pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int2, statInt));
    
                // Joint P-value
                for (size_t j = 0; j < nInt2; j++) {
                    for (size_t k = 0; k < nInt2; k++) {
                        StempGE2[j] += (Ai2[(nInt2 * j) + k] * betaInt2[k]);
                    }
                }
    
                statJoint = 0.0;
                for (size_t k = 0; k < nInt2; k++)
                    statJoint += betaInt2[k] * StempGE2[k];
                pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint2, statJoint));

                
                // Print
                for (size_t j = 0; j < nInt2; j++) {
                    oss2 << betaInt2[j] << "\t";
                }
                for (size_t ii = 0; ii < nInt2; ii++) {
                    oss2 << sqrt(rb_V2[iss + (ii * nInt2) + ii]) << "\t";
                }
                for (size_t ii = 0; ii < nInt2; ii++) {
                    for (size_t jj = 0; jj < nInt2; jj++) {
                        if (ii < jj) {
                            oss2 << rb_V2[iss + (ii * nInt2) + jj] << "\t";
                        }
                    }
                }
                oss2 << pvalInt << "\t" << pvalJoint << "\n";

                std::fill(StempE2.begin(), StempE2.end(), 0.0);
                std::fill(StempGE2.begin(), StempGE2.end(), 0.0);
                std::fill(betaInt2.begin(), betaInt2.end(), 0.0);
            }
    
            if (i % 100000 == 0)
            {
                results2 << oss2.str();
                oss2.str(std::string());
                oss2.clear();
            }
        }
        
        if(additionalInteraction){
            int is  = (i * nInt3);
            int iss = (i * nInt3 * nInt3);
    
            if (mb)
            {
                subMatrix(&mb_V3[0], Ai3, nInt3, nInt3, nInt3, nInt3, iss);
                subMatInv(&mb_V3[0], nInt3, iss);
    
    
                // Interaction effects
                for (size_t j = 0; j < nInt3; j++) {
                    for (size_t k = 0; k < nInt3; k++) {
                        betaInt3[j] += (mb_V3[iss + (nInt3 * j) + k] * mb_U3[is + k]);
                    }
                }
    
                 
    
                // Int P-value
                for (size_t j = 0; j < nInt3; j++) {
                    for (size_t k = 0; k < nInt3; k++) {
                        StempE3[j] += (Ai3[(nInt3 * j) + k] * betaInt3[k]);
                    }
                }

                statInt = 0.0;
                for (size_t k = 0; k < nInt3; k++)
                    statInt += betaInt3[k] * StempE3[k];
                pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int3, statInt));
  
                
                // Print
                for (size_t j = 0; j < nInt3; j++) {
                    oss3 << betaInt3[j] << "\t";
                }
                for (size_t ii = 0; ii < nInt3; ii++) {
                    oss3 << sqrt(mb_V3[iss + (ii * nInt3) + ii]) << "\t";
                }
                for (size_t ii = 0; ii < nInt3; ii++) {
                    for (size_t jj = 0; jj < nInt3; jj++) {
                        if (ii < jj) {
                            oss3 << mb_V3[iss + (ii * nInt3) + jj] << "\t";
                        }
                    }
                }
                
                oss3 << pvalInt << ((rb) ? "\t" : "\n"); 
    
                std::fill(StempE3.begin(), StempE3.end(), 0.0);
                std::fill(betaInt3.begin(), betaInt3.end(), 0.0);
            }
    
            if (rb)
            {
                subMatrix(&rb_V3[0], Ai3, nInt3, nInt3, nInt3, nInt3, iss);
                subMatInv(&rb_V3[0], nInt3, iss);
    
    
                // Interaction effects
                for (size_t j = 0; j < nInt3; j++) {
                    for (size_t k = 0; k < nInt3; k++) {
                        betaInt3[j] += (rb_V3[iss + (nInt3 * j) + k] * rb_U3[is + k]);
                    }
                }
    
    
                // Int P-value
                for (size_t j = 0; j < nInt3; j++) {
                    for (size_t k = 0; k < nInt3; k++) {
                        StempE3[j] += (Ai3[(nInt3 * j) + k] * betaInt3[k]);
                    }
                }
    
                statInt = 0.0;
                for (size_t k = 0; k < nInt3; k++)
                    statInt += betaInt3[k] * StempE3[k];
                pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int3, statInt));

                
                // Print
                for (size_t j = 0; j < nInt3; j++) {
                    oss3 << betaInt3[j] << "\t";
                }
                for (size_t ii = 0; ii < nInt3; ii++) {
                    oss3 << sqrt(rb_V3[iss + (ii * nInt3) + ii]) << "\t";
                }
                for (size_t ii = 0; ii < nInt3; ii++) {
                    for (size_t jj = 0; jj < nInt3; jj++) {
                        if (ii < jj) {
                            oss3 << rb_V3[iss + (ii * nInt3) + jj] << "\t";
                        }
                    }
                }
                oss3 << pvalInt << "\n";
    
                std::fill(StempE3.begin(), StempE3.end(), 0.0);
                std::fill(betaInt3.begin(), betaInt3.end(), 0.0);
            }
    
            if (i % 100000 == 0)
            {
                results3 << oss3.str();
                oss3.str(std::string());
                oss3.clear();
            }
        }
    }

    results << oss.str();
    oss.str(std::string());
    oss.clear();
    results.close();

    results2 << oss2.str();
    oss2.str(std::string());
    oss2.clear();
    results2.close();

    results3 << oss3.str();
    oss3.str(std::string());
    oss3.clear();
    results3.close();

    delete[] Ai;
    delete[] VE;

    delete[] Ai2;
    delete[] VE2;

    delete[] Ai3;
    
    printDone(2);
    printOutputLocation1(cmd.outFile);
    if(additionalJoint){
        printOutputLocation2(cmd.outFile2);
    }
    if(additionalInteraction){
        printOutputLocation3(cmd.outFile3);
    }
}

void printWelcome() {
    cout << "\n*********************************************************\n";
    cout << "Welcome to METAGEM v" << VERSION << "\n";
    cout << "(C) 2021-2025 Duy T. Pham, Han Chen and Shuyi Guo \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";
}

void printMetaBegin(int nFiles, int nvars) {
  cout << "Performing meta-analysis with " << std::to_string(nFiles) << " files and " << std::to_string(nvars) << " variants...\n";
}

void printProcessedFiles(std::string fileName, int n) {
  cout << "Processed file [" << fileName << "] with " << n << " variants.\n";
}

void printProcessingFiles() {
  cout << "Processing files...\n";
}

void printDone(int nbs) {

  if (nbs == 1) {
    cout << "Done.\n";
  } else if (nbs == 2) {
    cout << "Done.\n\n";
  } else {
    cerr << "\n ERROR: Invalid number.\n\n";
  }
}

void printOutputLocation1(std::string outFile) {
  cout << "Main meta-analysis results are in [" << outFile << "].\n";
}

void printOutputLocation2(std::string outFile) {
  cout << "Additional joint meta-analysis results are in [" << outFile << "].\n";
}

void printOutputLocation3(std::string outFile) {
  cout << "Additional interaction-only meta-analysis results are in [" << outFile << "].\n";
}

void printZeroVariantsError(std::string fileName) {
  cerr << "\nERROR: No variants in file [" << fileName << "].\n\n";
}

void printNColumnError(int z, int n, std::string fileName, int nheader) {
  cerr << "\nERROR: Unexpect number of columns (" << z << ") for variant number " << n+1 << "in file [" << fileName << "]. Expected: " << nheader << ".\n\n";
}

void printOpenFileError(std::string fileName) {
    cerr << "\nERROR: Cannot open the file [" << fileName << "].\n\n";
}

void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1) {
    cout << "\n*********************************************************\n";
    cout << "Wall Time = " << wall1 - wall0 << " (sec)\n";
    cout << "CPU Time  = " << cpu1  - cpu0  << " (sec)\n\n";
}
