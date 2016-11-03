/**
 * Thomas Keck 2015
 * NeuroBayes crosstests
 *
 * Reimplementation of the crosstests previously written using dsa.
 * Performs classification and density training and compares expertise and results
 * with corresponding reference files,
 */

#include "NeuroBayesTeacher.hh"
#include "NeuroBayesExpert.hh"

#include <random>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <stdexcept>

const int RANDSEED = 4575;
const double EPSILON_PREDICTIONS = 1e-3;
const double EPSILON_EXPERTISE = 1e-2;

/* == User logging struct== */
typedef struct log_struct_t {
    char* str;
} log_struct_t;

/* == User logging function == */
void logf(char* msg, int debug, void* enclosed) {
    log_struct_t* ls = (log_struct_t*)(enclosed);
      if(!ls->str) ls->str = strdup(msg);
}

class NeuroBayesTester {

  public:
    NeuroBayesTester(std::string name) : name(name), ec(nullptr), ls(nullptr) {
      tests.push_back(this);
      test_directory = std::string("./test_dir_") + name;
    }

    ~NeuroBayesTester() {
      if(ec != nullptr)
        free_ec1(ec);
      if(ls != nullptr) {
        if(ls->str != nullptr) {
          free(ls->str);
          ls->str = nullptr;
        }
        delete ls;
        ls = nullptr;
      }
    }

    void create() {
      
      std::system((std::string("rm -rf ") + test_directory).c_str());
      std::system((std::string("mkdir ") + test_directory).c_str());
      
      ls = new log_struct_t;
      ls->str = nullptr;

      nb = NeuroBayesTeacher::Instance(&ec,  normal_dbg, (log_func_t)(logf), (ls));
      nb->SayHello();
      nb->NB_DEF(true);
      nb->NB_RANVIN(RANDSEED);
      nb->SetOutputFile((test_directory + std::string("/expertise.nb")).c_str());
      nb->SetCArrayFile((test_directory + std::string("/results.csv")).c_str());
      nb->SetHistosFile((test_directory + std::string("/ahist.txt")).c_str());
     
      // Ensure correct default values
      nb->NB_DEF_LOSS("ENT");
      nb->NB_DEF_RTRAIN(1.0);
      nb->NB_DEF_SPEED(1.0);
      nb->NB_DEF_MAXLEARN(1.0);
      nb->NB_DEF_RELIMPORTANCE(0.0);
      nb->NB_DEF_LOSSWGT(1.0);
      nb->NB_DEF_LEARNDIAG(0);
      nb->NB_DEF_SPLOT_MODE(0);

    }

    bool checkForError() {

      if(ec != nullptr) {
        std::cerr << "EC: " << ec->kind << " " << ec->reason << " " << ec->stacktrace << std::endl;
        if(ls != nullptr and ls->str != nullptr)
          std::cerr << "LOG: " << ls->str << std::endl;
        return true;
      }
      return false;

    }

    std::vector<int> loadPreProFlags(std::string filename) {
      
      std::string line;
      std::string column;
      std::map<std::string, int> prepro_flags;
      
      std::ifstream prepro_file(filename);
      if(not prepro_file)
        throw std::runtime_error("Failed to open crosstest preprocessing file");;
      
      // Skip header
      std::getline(prepro_file, line);

      // Read in individual prepro flags
      while(std::getline(prepro_file, line)) {
        std::vector<std::string> tokens;
        std::stringstream stream_line(line);
        while(std::getline(stream_line, column, ';')) {
          tokens.push_back(column);
        }
        if(tokens.size() != 2) {
          throw std::runtime_error("Error during read of prepro.csv");
        }
        prepro_flags[tokens[0]] = std::stoi(tokens[1]);
      }

      std::vector<int> flags;
      for(auto &feature : feature_columns)
        flags.push_back(prepro_flags[feature]);

      return flags;

    }

    std::vector<std::vector<float>> loadData(std::string filename) {

      std::string line;
      std::string column;
      
      std::vector<int> index;

      std::ifstream data_file(filename);
      if(not data_file)
        throw std::runtime_error("Failed to open crosstest data file");
      
      // Get header
      std::vector<std::string> header;
      std::getline(data_file, line);
      std::stringstream stream_line(line);
      while(std::getline(stream_line, column, ';')) {
        header.push_back(column);
      }

      std::vector<std::string> columns = feature_columns;
      columns.push_back(target_column);
      columns.push_back(weight1_column);
      columns.push_back(weight2_column);

      for(auto& column : columns) {
        int pos = std::distance(header.begin(), std::find(header.begin(), header.end(), column));
        if(pos < header.size()) {
          index.push_back(pos);
        }
      }
      
      std::vector<std::vector<float>> data;

      while(std::getline(data_file, line)) {
        std::stringstream stream_line(line);
        std::vector<float> row;
        while(std::getline(stream_line, column, ';')) {
          // Important! We MUST read in all data as doubles using std::stod, because the feature matrix is singular
          // and therefore is extremly sensitive to numeric uncertainties!
          // But we save and pass all data as float, as Neurobayes expects it.
          row.push_back(std::stod(column));
        }
        
        std::vector<float> features(columns.size(), 0);
        features[columns.size() - 3] = 0.0;
        features[columns.size() - 2] = 1.0;
        features[columns.size() - 1] = 1.0;

        unsigned int i = 0;
        for(auto &j: index) {
          features[i] = row[j];
          i++;
        }
        data.push_back(features);
      }

      return data;

    }
    
    bool compareExpertise(std::string reference_filename) {
      
      std::ifstream expert_file(test_directory + std::string("/expertise.nb"));
      if(not expert_file)
        throw std::runtime_error("Failed to open expert file");
      
      std::ifstream reference_file(std::string("../data/crosstest_data/reference_files/") + reference_filename);
      if(not reference_file)
        throw std::runtime_error("Failed to open reference file for expertise");

      std::string tmp;
      double dummy;
      // Throw away first line of reference, which just contains the string expertise
      std::getline(reference_file, tmp);

      double expert_value = 0;
      double reference_value = 0;
      int i = -1;
      while((reference_file >> reference_value) && (expert_file >> expert_value)) {
        ++i;
        if( std::abs(reference_value - expert_value) > EPSILON_EXPERTISE) {
          std::cerr << "Produced expertise differs from reference expertise in line " << i << std::endl;
          std::cerr << "Reference value: " << reference_value << " Expert Value: " << expert_value << std::endl;
          return false;
        }
      }

      if(reference_file) {
        std::cerr << "Total number of values in reference_file > expert_file." << std::endl;
        std::cerr << "However the shared " << i << " numbers were equal" << std::endl;
        return false;
      }

      if(expert_file) {
         while(expert_file >> expert_value) {
            if( std::abs(expert_value) > 0.0001) {
              std::cerr << "Encountered additional non-zero values at the end of the expert_file" << std::endl;
              return false;
            }
         }
      }
      return true;

    }
    
    bool compareOutputOfExpertise(std::string reference_filename, std::vector<std::vector<float>> data) {

      std::string line;
      std::string column;
      
      std::ifstream expert_file(test_directory + std::string("/expertise.nb"));
      if(not expert_file)
        throw std::runtime_error("Failed to open expert file");
      
      std::ifstream reference_file(std::string("../data/crosstest_data/reference_files/") + reference_filename);
      if(not reference_file)
        throw std::runtime_error("Failed to open reference file for expertise");
      
      // Create NeuroBayes Expert
      Expert* nb_expert = new Expert(test_directory + std::string("/expertise.nb"), quiet_dbg, false, &ec);
      
      if (checkForError())
        return false;

      // First line contains header
      std::vector<std::string> header;
      std::getline(reference_file, line);
      std::stringstream stream_line(line);
      while(std::getline(stream_line, column, ';')) {
        header.push_back(column);
      }

      unsigned int row_size = data[0].size();
      for(unsigned int i = 0; i < data.size(); ++i) {
        if (not std::getline(reference_file, line)) {
          std::cerr << "Failed to read next entry from reference data file" << std::endl;
          return false;
        }
        std::vector<float> reference_data;
        std::stringstream stream_line(line);
        while(std::getline(stream_line, column, ';')) {
          reference_data.push_back(std::stof(column));
        }
        
        for(unsigned int j = 0; j < header.size(); ++j) {
          float expert_value = 0;
          float reference_value = reference_data[j];
          
          if (header[j] == "prediction") {
            expert_value = (nb_expert->nb_expert(&data[i][0], Expert::BINCLASS) + 1.0)/2.0;
          } else if(header[j] == "median") {
            expert_value = nb_expert->nb_expert(&data[i][0], Expert::MEDIAN);
          } else if(header[j] == "mean") {
            expert_value = nb_expert->nb_expert(&data[i][0], Expert::MEAN);
          } else if(header[j] == "T") {
            expert_value = data[i][row_size-3];
          } else {
            std::cerr << "Unkown column in reference data file: " << header[j] << std::endl;
            return false;
          }

          if (checkForError())
            return false;

          if( std::abs(reference_value - expert_value) > EPSILON_PREDICTIONS) {
            std::cerr << "Produced expertise output differs from reference data in event " << i << std::endl;
            std::cerr << "Reference value: " << reference_value << " Expert Value: " << expert_value << std::endl;
            return false;
          }
        }
      }

      delete nb_expert;
      return true;
    }

    virtual void setup() { }
    virtual void finish() { }
    
    bool call() {

      // Set Weight mode
      if(weight1_column != "" and weight2_column != "") {
        // User weight mode 2, which means that the externally given weight_factor is used.
        nb->NB_DEF_WEIGHT_MODE(2);
      } else if(weight1_column != "") {
        // Use weight mode 1, which means that the weight_factor is set to 1 internally.
        nb->NB_DEF_WEIGHT_MODE(1);
      } else {
        nb->NB_DEF_WEIGHT_MODE(1);
        // Use weight mode 1 (default)
      }
  
      // Set individual prepro flags
      std::vector<int> prepro_flags = loadPreProFlags("../data/crosstest_data/prepro.csv");
      for(unsigned int i = 0; i < prepro_flags.size(); ++i) {
        nb->SetIndividualPreproFlag(i, prepro_flags[i]);
      }

      // Set data
      std::vector<std::vector<float>> data = loadData("../data/crosstest_data/covertype_crosstest.csv");
      unsigned int size = feature_columns.size() + 3;
      for(auto &row : data) {
        if(row.size() != size) {
          throw std::runtime_error("Received data row is too small.");
        }

        nb->SetWeight(row[size-2], row[size-1]);
        nb->SetTarget(row[size-3]);
        nb->SetNextInput(size-3, &row[0]);
      }
 
      // Train network and output expertise and other files
      nb->TrainNet(true);
      if (checkForError())
        return false;
      nb->nb_correl_signi((test_directory + std::string("/correl_signi.txt")).c_str(),
                          (test_directory + std::string("/correl_signi.html")).c_str()); 
      
      // Check if results are correct
      if (not compareExpertise(name + std::string("_expertise.csv"))) {
        return false;
      }
      if (not compareOutputOfExpertise(name + std::string(".csv"), data)) {
        return false;
      }

      return true;
    }

    void destroy() {
      
      //std::system((std::string("rm -rf ") + test_directory).c_str());
      if(ec != nullptr)
        free_ec1(ec);
      ec = nullptr;
      if(ls != nullptr) {
        if(ls->str != nullptr) {
          free(ls->str);
          ls->str = nullptr;
        }
        delete ls;
        ls = nullptr;
      }
    }

    static bool run() {
      bool ok = true;
      for(auto& test : tests) {
        test->create();
        test->setup();
        std::cerr << test->name << " ";
        try {
          if(test->call()) {
            std::cerr << "successfull";
            test->finish();
          } else {
            std::cerr << "failed";
            ok = false;
          }
        } catch(std::runtime_error &e) {
          std::cerr << e.what() << std::endl;
          return false;
        }
        std::cerr << std::endl;
        test->destroy();
      }
      return ok;
    }

  private:
    std::string test_directory;
    std::string name;
    ec_t* ec;
    log_struct_t* ls;

    static std::list<NeuroBayesTester*> tests;

  protected:
    NeuroBayesTeacher* nb;
    std::string target_column;
    std::string weight1_column;
    std::string weight2_column;
    std::vector<std::string> feature_columns;

};

std::list<NeuroBayesTester*> NeuroBayesTester::tests;

class ClassificationNeuroBayesTester : public NeuroBayesTester {
  public:
    ClassificationNeuroBayesTester(std::string name) : NeuroBayesTester(std::string("classify_") + name) { }
    virtual void setup() {
      nb->NB_DEF_TASK("CLA");
      nb->NB_DEF_NODE1(25);
      nb->NB_DEF_NODE2(13);
      nb->NB_DEF_NODE3(1);
      target_column = "Cover_Type";
      feature_columns = {"Elevation", "Aspect", "Slope", "Horizontal_Distance_To_Hydrology",
                         "Vertical_Distance_To_Hydrology", "Horizontal_Distance_To_Roadways", "Hillshade_9am",
                         "Hillshade_Noon", "Hillshade_3pm", "Horizontal_Distance_To_Fire_Points",
                         "Wilderness_Area_Rawah", "Wilderness_Area_Neota", "Wilderness_Area_Comanche_Peak",
                         "Wilderness_Area_Cache_la_Poudre", "Soil_Type0", "Soil_Type1", "Soil_Type2", "Soil_Type3",
                         "Soil_Type4", "Soil_Type5", "Soil_Type6", "Soil_Type7", "Soil_Type8", "Soil_Type9"};
    }
};

class ClassificationInternalBoost : public ClassificationNeuroBayesTester {
  public:
    ClassificationInternalBoost() : ClassificationNeuroBayesTester("internal_boost") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(422);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(2);
    }
} _ClassificationInternalBoost;

class ClassificationWithWeights : public ClassificationNeuroBayesTester {
  public:
    ClassificationWithWeights() : ClassificationNeuroBayesTester("weights") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(422);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(0);
      weight1_column = "W1";
    }
} _ClassificationWithWeights;

class ClassificationWithWeights2 : public ClassificationNeuroBayesTester {
  public:
    ClassificationWithWeights2() : ClassificationNeuroBayesTester("weights2") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(422);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(0);
      weight1_column = "W1";
      weight2_column = "W2";
    }
} _ClassificationWithWeights2;

class ClassificationRobust : public ClassificationNeuroBayesTester {
  public:
    ClassificationRobust() : ClassificationNeuroBayesTester("robust") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(422);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(0);
    }
} _ClassificationRobust;

class ClassificationSimpleNet : public ClassificationNeuroBayesTester {
  public:
    ClassificationSimpleNet() : ClassificationNeuroBayesTester("simple_net") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("OFF");
      nb->NB_DEF_ITER(2);
    }
} _ClassificationSimpleNet;

class ClassificationNetDia2 : public ClassificationNeuroBayesTester {
  public:
    ClassificationNetDia2() : ClassificationNeuroBayesTester("net_dia2") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(2);
    }
} _ClassificationNetDia2;

class ClassificationNetDiag : public ClassificationNeuroBayesTester {
  public:
    ClassificationNetDiag() : ClassificationNeuroBayesTester("net_diag") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("DIAG");
      nb->NB_DEF_ITER(2);
    }
} _ClassificationNetDiag;

class ClassificationNetBFGS : public ClassificationNeuroBayesTester {
  public:
    ClassificationNetBFGS() : ClassificationNeuroBayesTester("net_bfgs") { }
    virtual void setup() {
      ClassificationNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(100);
      nb->NB_DEF_METHOD("BFGS"); 
    }
} _ClassificationNetBFGS;

class DensityNeuroBayesTester : public NeuroBayesTester {
  public:
    DensityNeuroBayesTester(std::string name) : NeuroBayesTester(std::string("density_") + name) { }
    virtual void setup() {
      nb->NB_DEF_TASK("DEN"); 
      nb->NB_DEF_NODE1(25);
      nb->NB_DEF_NODE2(22);
      nb->NB_DEF_NODE3(20);
      target_column = "Elevation";
      feature_columns = {"Aspect", "Slope", "Horizontal_Distance_To_Hydrology", "Vertical_Distance_To_Hydrology",
                         "Horizontal_Distance_To_Roadways", "Hillshade_9am", "Hillshade_Noon", "Hillshade_3pm",
                         "Horizontal_Distance_To_Fire_Points", "Wilderness_Area_Rawah", "Wilderness_Area_Neota",
                         "Wilderness_Area_Comanche_Peak", "Wilderness_Area_Cache_la_Poudre", "Soil_Type0",
                         "Soil_Type1", "Soil_Type2", "Soil_Type3", "Soil_Type4", "Soil_Type5", "Soil_Type6",
                         "Soil_Type7", "Soil_Type8", "Soil_Type9", "Cover_Type"};
    }
};

class DensityWithWeights : public DensityNeuroBayesTester {
  public:
    DensityWithWeights() : DensityNeuroBayesTester("weights") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(432);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(0);
      weight1_column = "W1";
    }
} _DensityWithWeights;

class DensityWithWeights2 : public DensityNeuroBayesTester {
  public:
    DensityWithWeights2() : DensityNeuroBayesTester("weights2") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(432);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(0);
      weight1_column = "W1";
      weight2_column = "W2";
    }
} _DensityWithWeights2;

class DensityRobust : public DensityNeuroBayesTester {
  public:
    DensityRobust() : DensityNeuroBayesTester("robust") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(432);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(0);
    }
} _DensityRobust;

class DensityInternalBoost : public DensityNeuroBayesTester {
  public:
    DensityInternalBoost() : DensityNeuroBayesTester("internal_boost") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(432);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(2);
    }
} _DensityInternalBoost;

class DensitySimpleNet : public DensityNeuroBayesTester {
  public:
    DensitySimpleNet() : DensityNeuroBayesTester("simple_net") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("OFF");
      nb->NB_DEF_ITER(2);
    }
} _DensitySimpleNet;

class DensityNetDia2 : public DensityNeuroBayesTester {
  public:
    DensityNetDia2() : DensityNeuroBayesTester("net_dia2") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(2);
    }
} _DensityNetDia2;

class DensityNetDiag : public DensityNeuroBayesTester {
  public:
    DensityNetDiag() : DensityNeuroBayesTester("net_diag") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("DIAG");
      nb->NB_DEF_ITER(2);
    }
} _DensityNetDiag;

class DensityNetBFGS : public DensityNeuroBayesTester {
  public:
    DensityNetBFGS() : DensityNeuroBayesTester("net_bfgs") { }
    virtual void setup() {
      DensityNeuroBayesTester::setup();
      nb->NB_DEF_PRE(412);
      nb->NB_DEF_SHAPE("DIA2");
      nb->NB_DEF_ITER(100);
      nb->NB_DEF_METHOD("BFGS"); 
    }
} _DensityNetBFGS;


int main(int argc, char** argv) {
 
  if(NeuroBayesTester::run())
    return 0;
  else
    return -1;

}

