#include "metagem.h"

void processFileHeader(int nInt1, int nInt2, int nInt3, bool mb, bool rb, bool additionalJoint, bool additionalInteraction, bool renameHeaders, std::vector<std::string> lc_intNames, std::vector<std::string> lc_intNames2, std::vector<std::string> lc_intNames3, std::vector<std::string> fileNames, std::map<std::string, std::map<std::string, std::string>> headerRenamings, FileInfo* fip) 
{
    
    for (size_t f = 0; f < fileNames.size(); f++) {

        std::string fileName = fileNames[f];

        // Read input file
        std::ifstream file;
        file.open(fileName);
        if (!file.is_open()) {
            printOpenFileError(fileName);
            exit(1);
        }
        fip->fileNames.push_back(fileName);

        // Ignore comments in the text file (ex. #dispersion: )
        std::string line;
        getline(file, line);
        while(line.rfind("#", 0) == 0)
        {
            getline(file, line);
        }
        file.close();

        // Read the column headers
        std::unordered_map<std::string, int> columnNames;

        int header_i = 0;
        std::string header;
        std::istringstream iss(line);
        while (getline(iss, header, '\t')) 
        {
            header.erase(std::remove(header.begin(), header.end(), '\r'), header.end());

            if (columnNames.find(header) != columnNames.end()) {
                cerr << "\nERROR: There are duplicate header names (" << header << ") in file [" << fileName << "].\n\n";
                file.close();
                exit(1);
            }
            columnNames[header] = header_i;

            header_i++;
        }
        fip->nheader[fileName] = header_i;

        
        // Check if there are header renamings for this file
        if(renameHeaders){
            if (headerRenamings.find(fileName) != headerRenamings.end()) {
            // Apply header name changes
            const auto& mappings = headerRenamings[fileName];
                for (const auto& mapping : mappings) {
                    if (columnNames.find(mapping.first) != columnNames.end()) {
                        int colIndex = columnNames[mapping.first];
                        columnNames.erase(mapping.first);
                        columnNames[mapping.second] = colIndex;
                    } else {
                        std::cerr << "ERROR: Column '" << mapping.first << "' not found in file '" << fileName << "' for renaming to '" << mapping.second << "'.\n";
                        exit(1);
                    }
                }
            }
        }
        

        // Get variant information columns
        int snpid_col = 0;
        if (columnNames.find("SNPID") != columnNames.end()) 
        {
            snpid_col = columnNames["SNPID"];
            fip->snpColumn[fileName] = snpid_col;
        } 
        else 
        {
            printHeaderMissingError(fileName, "SNPID");
            exit(1);
        }

        if (columnNames.find("CHR") != columnNames.end())
        {
            fip->chrColumn[fileName] = columnNames["CHR"];
        }
        else
        {
            printHeaderMissingError(fileName, "CHR");
            exit(1);
        }

        if (columnNames.find("POS") != columnNames.end())
        {
            fip->posColumn[fileName] = columnNames["POS"];
        }
        else
        {
            printHeaderMissingError(fileName, "POS");
            exit(1);
        }

        if (columnNames.find("Non_Effect_Allele") != columnNames.end())
        {
            fip->nonEffectColumn[fileName] = columnNames["Non_Effect_Allele"];
        }
        else
        {
            printHeaderMissingError(fileName, "Non_Effect_Allele");
            exit(1);
        }

        if (columnNames.find("Effect_Allele") != columnNames.end())
        {
            fip->effectColumn[fileName] = columnNames["Effect_Allele"];
        }
        else
        {
            printHeaderMissingError(fileName, "Effect_Allele");
            exit(1);
        }

        if (columnNames.find("N_Samples") != columnNames.end())
        {
            fip->nSampleColumn[fileName] = columnNames["N_Samples"];
        }
        else
        {
            printHeaderMissingError(fileName, "N_Samples");
            exit(1);
        }

        if (columnNames.find("AF") != columnNames.end())
        {
            fip->freqColumn[fileName] = columnNames["AF"];
        }
        else
        {
            printHeaderMissingError(fileName, "AF");
            exit(1);
        }
        


        // Get the Beta Marginal column
        if (columnNames.find("Beta_Marginal") != columnNames.end())
        {
            fip->betaMargColumn[fileName] = columnNames["Beta_Marginal"];
        }
        else
        {
            printHeaderMissingError(fileName, "Beta_Marginal");
            exit(1);
        }

        if (mb)
        {
            // Get marginal model based SE columns if it exists
            if (columnNames.find("SE_Beta_Marginal") != columnNames.end())
            {
                fip->mb_seMargColumn[fileName] = columnNames["SE_Beta_Marginal"];
            }
            else 
            {
                printHeaderMissingError(fileName, "SE_Beta_Marginal");
                exit(1);
            }
        }

        if (rb)
        {
            // Get marginal robust SE columns if it exists
            if (columnNames.find("robust_SE_Beta_Marginal") != columnNames.end())
            {
                fip->rb_seMargColumn[fileName] = columnNames["robust_SE_Beta_Marginal"];
            }
            else 
            {
                printHeaderMissingError(fileName, "robust_SE_Beta_Marginal");
                exit(1);
            }
        }
            

        // Get Beta Interaction columns
        int nints = 0;
        std::vector<std::string> betaIntNames;
        if (columnNames.find("Beta_G") == columnNames.end()) {
            printHeaderMissingError(fileName, "Beta_G");
            exit(1);
        }
        else
        {
            betaIntNames.push_back("G");
            nints++;

            for (std::pair<std::string, int> element : columnNames)
            {
                std::string key = element.first;
                if (key.rfind("Beta_G-", 0) == 0) {
                    key.erase(0, 5);
                    betaIntNames.push_back(key);
                    nints++;
                }
            }
        }


        if(nints != nInt1)
        {
            cerr << "\nERROR: The file [" << fileName << "] does not contain the same number of interaction terms as --exposure-names.\n\n";
            exit(1);
        }


        // Convert to lowercase strings for matching
        std::vector<std::string> lc_betaIntNames = betaIntNames;
        for(std::string &s : lc_betaIntNames)
        {
            std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
        }
            
        // Order the interactions base on the first vector of interactions
        std::vector<int> betaIntColumn;
        std::vector<std::string> ord_betaIntNames;
        for (int i = 0; i < nInt1; i++)
        {
            auto it = std::find(lc_betaIntNames.begin(), lc_betaIntNames.end(), lc_intNames[i]);
            if (it != lc_betaIntNames.end())
            {
                auto idx = std::distance(lc_betaIntNames.begin(), it);
                ord_betaIntNames.push_back(betaIntNames[idx]);
                betaIntColumn.push_back(columnNames["Beta_" + betaIntNames[idx]]);
            }
            else 
            {
                cerr << "\nERROR: The file [" << fileName << "] does not contain the GxE term: " << lc_intNames[i] << ".\n\n";
                exit(1);
            }
        }
        fip->betaIntColumn[fileName] = betaIntColumn;

        


        // Get columns containing the model-based summary statistics for the interaction terms
        if (mb)
        {
            std::vector<int> mb_covIntColumn(nints * nints);

            for (int i = 0; i < nints; i++) 
            {
                std::string s1 = "SE_Beta_" + ord_betaIntNames[i];
                if (columnNames.find(s1) != columnNames.end()) 
                {
                    mb_covIntColumn[i*nints + i] = columnNames[s1];
                }
                else 
                {
                    printHeaderMissingError(fileName,  s1);
                    exit(1);
                }
            }

            for (int i = 0; i < nints; i++) 
            {
                for (int j = i+1; j < nints; j++) 
                {
                    int idx;
                    std::string s2 = "Cov_Beta_" + ord_betaIntNames[i] + "_" + ord_betaIntNames[j];
                    std::string s3 = "Cov_Beta_" + ord_betaIntNames[j] + "_" + ord_betaIntNames[i];
                    if (columnNames.find(s2) != columnNames.end()) 
                    {
                        idx = columnNames[s2];
                    } 
                    else if (columnNames.find(s3) != columnNames.end())
                    {
                        idx = columnNames[s3];
                    }
                    else 
                    {
                        cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
                        exit(1);
                    }
                    mb_covIntColumn[(i * nints) + j] = idx;
                    mb_covIntColumn[(j * nints) + i] = idx;
                }
            }
            fip->mb_covIntColumn[fileName] = mb_covIntColumn;
        }


        // Get columns containing the robust summary statistics for the interaction terms	
        if (rb)
        {
            std::vector<int> rb_covIntColumn(nints * nints);
            for (int i = 0; i < nints; i++) 
            {
                std::string s1 = "robust_SE_Beta_" + ord_betaIntNames[i];
                if (columnNames.find(s1) != columnNames.end()) 
                {
                    rb_covIntColumn[i*nints + i] = columnNames[s1];
                }
                else 
                {
                    printHeaderMissingError(fileName,  s1);
                    exit(1);
                }
            }

            for (int i = 0; i < nints; i++) 
            {
                for (int j = i+1; j < nints; j++) 
                {
                    int idx;
                    std::string s2 = "robust_Cov_Beta_" + ord_betaIntNames[i] + "_" + ord_betaIntNames[j];
                    std::string s3 = "robust_Cov_Beta_" + ord_betaIntNames[j] + "_" + ord_betaIntNames[i];

                    if (columnNames.find(s2) != columnNames.end()) 
                    {
                        idx = columnNames[s2];
                    } 
                    else if (columnNames.find(s3) != columnNames.end())
                    {
                        idx = columnNames[s3];
                    }
                    else 
                    {
                        cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
                        exit(1);
                    }
                    rb_covIntColumn[(i * nints) + j] = idx;
                    rb_covIntColumn[(j * nints) + i] = idx;
                }
            }
            fip->rb_covIntColumn[fileName] = rb_covIntColumn;
        }

        // Additional joint test
        if(additionalJoint){
            // Order the interactions base on the first vector of interactions
            std::vector<int> betaIntColumn2;
            std::vector<std::string> ord_betaIntNames2;
            for (int i = 0; i < nInt2; i++)
            {
                auto it = std::find(lc_betaIntNames.begin(), lc_betaIntNames.end(), lc_intNames2[i]);
                if (it != lc_betaIntNames.end())
                {
                    auto idx = std::distance(lc_betaIntNames.begin(), it);
                    ord_betaIntNames2.push_back(betaIntNames[idx]);
                    betaIntColumn2.push_back(columnNames["Beta_" + betaIntNames[idx]]);
                }
                else 
                {
                    cerr << "\nERROR: The file [" << fileName << "] does not contain the GxE term: " << lc_intNames2[i] << ".\n\n";
                    exit(1);
                }
            }
            fip->betaIntColumn2[fileName] = betaIntColumn2;

        


            // Get columns containing the model-based summary statistics for the interaction terms
            if (mb)
            {
                std::vector<int> mb_covIntColumn2(nInt2 * nInt2);

                for (int i = 0; i < nInt2; i++) 
                {
                    std::string s1 = "SE_Beta_" + ord_betaIntNames2[i];
                    if (columnNames.find(s1) != columnNames.end()) 
                    {
                        mb_covIntColumn2[i*nInt2 + i] = columnNames[s1];
                    }
                    else 
                    {
                        printHeaderMissingError(fileName,  s1);
                        exit(1);
                    }
                }

                for (int i = 0; i < nInt2; i++) 
                {
                    for (int j = i+1; j < nInt2; j++) 
                    {
                        int idx;
                        std::string s2 = "Cov_Beta_" + ord_betaIntNames2[i] + "_" + ord_betaIntNames2[j];
                        std::string s3 = "Cov_Beta_" + ord_betaIntNames2[j] + "_" + ord_betaIntNames2[i];
                        if (columnNames.find(s2) != columnNames.end()) 
                        {
                            idx = columnNames[s2];
                        } 
                        else if (columnNames.find(s3) != columnNames.end())
                        {
                            idx = columnNames[s3];
                        }
                        else 
                        {
                            cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
                            exit(1);
                        }
                        mb_covIntColumn2[(i * nInt2) + j] = idx;
                        mb_covIntColumn2[(j * nInt2) + i] = idx;
                    }
                }
                fip->mb_covIntColumn2[fileName] = mb_covIntColumn2;
            }


            // Get columns containing the robust summary statistics for the interaction terms	
            if (rb)
            {
                std::vector<int> rb_covIntColumn2(nInt2 * nInt2);
                for (int i = 0; i < nInt2; i++) 
                {
                    std::string s1 = "robust_SE_Beta_" + ord_betaIntNames2[i];
                    if (columnNames.find(s1) != columnNames.end()) 
                    {
                        rb_covIntColumn2[i*nInt2 + i] = columnNames[s1];
                    }
                    else 
                    {
                        printHeaderMissingError(fileName,  s1);
                        exit(1);
                    }
                }

                for (int i = 0; i < nInt2; i++) 
                {
                    for (int j = i+1; j < nInt2; j++) 
                    {
                        int idx;
                        std::string s2 = "robust_Cov_Beta_" + ord_betaIntNames2[i] + "_" + ord_betaIntNames2[j];
                        std::string s3 = "robust_Cov_Beta_" + ord_betaIntNames2[j] + "_" + ord_betaIntNames2[i];

                        if (columnNames.find(s2) != columnNames.end()) 
                        {
                            idx = columnNames[s2];
                        } 
                        else if (columnNames.find(s3) != columnNames.end())
                        {
                            idx = columnNames[s3];
                        }
                        else 
                        {
                            cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
                            exit(1);
                        }
                        rb_covIntColumn2[(i * nInt2) + j] = idx;
                        rb_covIntColumn2[(j * nInt2) + i] = idx;
                    }
                }
                fip->rb_covIntColumn2[fileName] = rb_covIntColumn2;
            }
        }

        // Additional interaction-only test
        if(additionalInteraction){
            // Order the interactions base on the first vector of interactions
            std::vector<int> betaIntColumn3;
            std::vector<std::string> ord_betaIntNames3;
            for (int i = 0; i < nInt3; i++)
            {
                auto it = std::find(lc_betaIntNames.begin(), lc_betaIntNames.end(), lc_intNames3[i]);
                if (it != lc_betaIntNames.end())
                {
                    auto idx = std::distance(lc_betaIntNames.begin(), it);
                    std::string actualName = betaIntNames[idx]; // This will be "G" or an exposure name
                    ord_betaIntNames3.push_back(actualName);
                    
                    // Construct the column name correctly
                    std::string colName = "Beta_" + actualName; 
                    if (columnNames.find(colName) != columnNames.end()) {
                        betaIntColumn3.push_back(columnNames[colName]);
                    } else {
                        cerr << "\nERROR: Could not find column " << colName << " in " << fileName << "\n";
                        exit(1);
                    }
                }
                else 
                {
                    cerr << "\nERROR: The file [" << fileName << "] does not contain the term: " << lc_intNames3[i] << ".\n\n";
                    exit(1);
                }
            }
            fip->betaIntColumn3[fileName] = betaIntColumn3;

        


            // Get columns containing the model-based summary statistics for the interaction terms
            if (mb)
            {
                std::vector<int> mb_covIntColumn3(nInt3 * nInt3);

                for (int i = 0; i < nInt3; i++) 
                {
                    std::string s1 = "SE_Beta_" + ord_betaIntNames3[i];
                    if (columnNames.find(s1) != columnNames.end()) 
                    {
                        mb_covIntColumn3[i*nInt3 + i] = columnNames[s1];
                    }
                    else 
                    {
                        printHeaderMissingError(fileName,  s1);
                        exit(1);
                    }
                }

                for (int i = 0; i < nInt3; i++) 
                {
                    for (int j = i+1; j < nInt3; j++) 
                    {
                        int idx;
                        std::string s2 = "Cov_Beta_" + ord_betaIntNames3[i] + "_" + ord_betaIntNames3[j];
                        std::string s3 = "Cov_Beta_" + ord_betaIntNames3[j] + "_" + ord_betaIntNames3[i];
                        if (columnNames.find(s2) != columnNames.end()) 
                        {
                            idx = columnNames[s2];
                        } 
                        else if (columnNames.find(s3) != columnNames.end())
                        {
                            idx = columnNames[s3];
                        }
                        else 
                        {
                            cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
                            exit(1);
                        }
                        mb_covIntColumn3[(i * nInt3) + j] = idx;
                        mb_covIntColumn3[(j * nInt3) + i] = idx;
                    }
                }
                fip->mb_covIntColumn3[fileName] = mb_covIntColumn3;
            }


            // Get columns containing the robust summary statistics for the interaction terms	
            if (rb)
            {
                std::vector<int> rb_covIntColumn3(nInt3 * nInt3);
                for (int i = 0; i < nInt3; i++) 
                {
                    std::string s1 = "robust_SE_Beta_" + ord_betaIntNames3[i];
                    if (columnNames.find(s1) != columnNames.end()) 
                    {
                        rb_covIntColumn3[i*nInt3 + i] = columnNames[s1];
                    }
                    else 
                    {
                        printHeaderMissingError(fileName,  s1);
                        exit(1);
                    }
                }

                for (int i = 0; i < nInt3; i++) 
                {
                    for (int j = i+1; j < nInt3; j++) 
                    {
                        int idx;
                        std::string s2 = "robust_Cov_Beta_" + ord_betaIntNames3[i] + "_" + ord_betaIntNames3[j];
                        std::string s3 = "robust_Cov_Beta_" + ord_betaIntNames3[j] + "_" + ord_betaIntNames3[i];

                        if (columnNames.find(s2) != columnNames.end()) 
                        {
                            idx = columnNames[s2];
                        } 
                        else if (columnNames.find(s3) != columnNames.end())
                        {
                            idx = columnNames[s3];
                        }
                        else 
                        {
                            cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
                            exit(1);
                        }
                        rb_covIntColumn3[(i * nInt3) + j] = idx;
                        rb_covIntColumn3[(j * nInt3) + i] = idx;
                    }
                }
                fip->rb_covIntColumn3[fileName] = rb_covIntColumn3;
            }
        }
    }
}



void printOutputHeader(bool mb, bool rb, bool additionalJoint, bool additionalInteraction, std::string output, std::string output2, std::string output3, size_t nInt1, size_t nInt2, size_t nInt3, std::vector<std::string> intNames, std::vector<std::string> intNames2, std::vector<std::string> intNames3) 
{
    for(std::string &s : intNames){
        s = "G-" + s;
    }
    intNames.insert(intNames.begin(), "G");


    std::ofstream results(output, std::ofstream::binary);

    results << "SNPID\t" << "CHR\t" << "POS\t" << "Non_Effect_Allele\t" << "Effect_Allele\t" << "N_files\t" << "N_Samples\t" << "AF\t";
    
    if (mb)
    {
        results << "Beta_Marginal\t" << "SE_Beta_Marginal\t";
       
        // Print Int beta header
        for (size_t i = 0; i < nInt1; i++) 
        {
            results << "Beta_" << intNames[i] << "\t";
        }

        // Print model-based covariance
        for (size_t i = 0; i < nInt1; i++) 
        {
            results << "SE_Beta_" << intNames[i] << "\t"; 
        }
        for (size_t i = 0; i < nInt1; i++) 
        {
            for (size_t j = 0; j < nInt1; j++) 
            {
                if (i < j) 
                {
                    results << "Cov_Beta_" << intNames[i] << "_" << intNames[j] << "\t"; 
                }
            }
        }

        // Print p-values
        results << "P_Value_Marginal\t" << "P_Value_Interaction\t" << "P_Value_Joint" << ((rb) ? "\t" : "\n");
    }

    if (rb)
    {
        results << "robust_Beta_Marginal\t" << "robust_SE_Beta_Marginal\t";

        // Print beta header
        for (size_t i = 0; i < nInt1; i++) 
        {
            results << "robust_Beta_" << intNames[i] << "\t";
        }

        // Print model-based covariance
        for (size_t i = 0; i < nInt1; i++) 
        {
            results << "robust_SE_Beta_" << intNames[i] << "\t"; 
        }
        for (size_t i = 0; i < nInt1; i++) 
        {
            for (size_t j = 0; j < nInt1; j++) 
            {
                if (i < j) 
                {
                    results << "robust_Cov_Beta_" << intNames[i] << "_" << intNames[j] << "\t"; 
                }
            }
        }

        results << "robust_P_Value_Marginal\t" << "robust_P_Value_Interaction\t" << "robust_P_Value_Joint\n";
    }
    
    results.close();

    if(additionalJoint){
        for(std::string &s : intNames2){
            s = "G-" + s;
        }
        intNames2.insert(intNames2.begin(), "G");

        std::ofstream results2(output2, std::ofstream::binary);

        results2 << "SNPID\t" << "CHR\t" << "POS\t" << "Non_Effect_Allele\t" << "Effect_Allele\t" << "N_files\t" << "N_Samples\t" << "AF\t";
    
        if (mb)
        {
            // Print Int beta header
            for (size_t i = 0; i < nInt2; i++) 
            {
                results2 << "Beta_" << intNames2[i] << "\t";
            }

            // Print model-based covariance
            for (size_t i = 0; i < nInt2; i++) 
            {
                results2 << "SE_Beta_" << intNames2[i] << "\t"; 
            }
            for (size_t i = 0; i < nInt2; i++) 
            {
                for (size_t j = 0; j < nInt2; j++) 
                {
                    if (i < j) 
                    {
                        results2 << "Cov_Beta_" << intNames2[i] << "_" << intNames2[j] << "\t"; 
                    }
                }
            }

            results2 << "P_Value_Interaction\t" << "P_Value_Joint" << ((rb) ? "\t" : "\n");
        }

        if (rb)
        {
            // Print beta header
            for (size_t i = 0; i < nInt2; i++) 
            {
                results2 << "robust_Beta_" << intNames2[i] << "\t";
            }

            // Print model-based covariance
            for (size_t i = 0; i < nInt2; i++) 
            {
                results2 << "robust_SE_Beta_" << intNames2[i] << "\t"; 
            }
            for (size_t i = 0; i < nInt2; i++) 
            {
                for (size_t j = 0; j < nInt2; j++) 
                {
                    if (i < j) 
                    {
                        results2 << "robust_Cov_Beta_" << intNames2[i] << "_" << intNames2[j] << "\t"; 
                    }
                }
            }
            
            results2 << "robust_P_Value_Interaction\t" << "robust_P_Value_Joint\n";
        }
    
        results2.close();
    }
    
    if(additionalInteraction){
        for(std::string &s : intNames3){
            if (std::tolower(s[0]) != 'g' || s.length() > 1) { // Only prefix if not already "G" or "g"
                 s = "G-" + s;
            }
        }


        std::ofstream results3(output3, std::ofstream::binary);

        results3 << "SNPID\t" << "CHR\t" << "POS\t" << "Non_Effect_Allele\t" << "Effect_Allele\t" << "N_files\t" << "N_Samples\t" << "AF\t";
    
        if (mb)
        {
            // Print Int beta header
            for (size_t i = 0; i < nInt3; i++) 
            {
                results3 << "Beta_" << intNames3[i] << "\t";
            }

            // Print model-based covariance
            for (size_t i = 0; i < nInt3; i++) 
            {
                results3 << "SE_Beta_" << intNames3[i] << "\t"; 
            }
            for (size_t i = 0; i < nInt3; i++) 
            {
                for (size_t j = 0; j < nInt3; j++) 
                {
                    if (i < j) 
                    {
                        results3 << "Cov_Beta_" << intNames3[i] << "_" << intNames3[j] << "\t"; 
                    }
                }
            }

            results3 << "P_Value_Interaction" << ((rb) ? "\t" : "\n");
        }

        if (rb)
        {
            // Print beta header
            for (size_t i = 0; i < nInt3; i++) 
            {
                results3 << "robust_Beta_" << intNames3[i] << "\t";
            }

            // Print model-based covariance
            for (size_t i = 0; i < nInt3; i++) 
            {
                results3 << "robust_SE_Beta_" << intNames3[i] << "\t"; 
            }
            for (size_t i = 0; i < nInt3; i++) 
            {
                for (size_t j = 0; j < nInt3; j++) 
                {
                    if (i < j) 
                    {
                        results3 << "robust_Cov_Beta_" << intNames3[i] << "_" << intNames3[j] << "\t"; 
                    }
                }
            }
            
            results3 << "robust_P_Value_Interaction\n";
        }
    
        results3.close();
    }

}


void printHeaderMissingError(std::string fileName, std::string column) {
    cerr << "\nERROR: The file [" << fileName << "] does not contain a " << column << " column.\n\n";
}
