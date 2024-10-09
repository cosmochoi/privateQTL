#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_file output_tsv_file" << std::endl;
        return 1;
    }

    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open input file." << std::endl;
        return 1;
    }

    std::ofstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output TSV file." << std::endl;
        return 1;
    }

    // Write the TSV header
    outputFile << "phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\ttstat\tpval_nom\tslope\tslope_se\tpval_perm\tpval_beta\n";

    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.find("Gene: ") != std::string::npos) {
            // Parse the relevant information from the line
            std::string gene;
            double beta_shape1, beta_shape2, dof, tstat, slope, slope_se, pval_perm, pval_beta,pval_true_df, pval_nom;
            int variantIdx;
            std::string temp;

            std::istringstream iss(line);
            iss >> gene >> gene;  // Extract the gene name

            std::getline(inputFile, line);  // Read the next line with Variant information
            iss.str(line);
            iss.clear();
            std::string variantId;
            iss >> temp>>variantId;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>>variantIdx;

            std::getline(inputFile, line);  // Skip the line with "std_ratio" information
            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>>dof;

            // Read the relevant information from the following lines
            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>>slope;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>>tstat;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> slope_se;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> pval_true_df;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> pval_nom;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> pval_perm;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> beta_shape1;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> beta_shape2;

            std::getline(inputFile, line);
            iss.str(line);
            iss.clear();
            iss >> temp>> pval_beta;

            // outputFile << "phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\ttstat\tslope\tslope_se\tpval_perm\tpval_beta\n";
            outputFile << gene<<"\t"<<variantId<<"\t"<<variantIdx<<"\t"<<beta_shape1<<"\t"<<beta_shape2
                       << "\t" << dof <<"\t"<<pval_true_df<<"\t"<<tstat<<"\t"<<pval_nom<<"\t"<<slope<<"\t"<<slope_se<<"\t"<<pval_perm
                       << "\t"<<pval_beta<<"\n";
            // Write the parsed information to the output TSV file

        }
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Conversion completed. Output written to " << outputFileName << std::endl;

    return 0;
}
