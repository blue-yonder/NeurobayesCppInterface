/**
 * Thomas Keck 2015
 * NeuroBayes C++ Interface tests
 * Test the correct failure of NeuroBayes for certain 
 * combinations of settings
 */

#include "NeuroBayesTeacher.hh"
#include "NeuroBayesExpert.hh"

#include <random>
#include <cstdlib>
#include <string.h>

const int nvar = 50;
const int nevents = 201;
float input_matrix[nevents][nvar];

/* == User logging struct== */
typedef struct log_struct_t {
    char* str;
} log_struct_t;

/* == User logging function == */
void logf(char* msg, int debug, void* enclosed) {
    log_struct_t* ls = (log_struct_t*)(enclosed);
      if(!ls->str) ls->str = strdup(msg);
}

bool training() {

  bool testok = true;
  // Create random distributions
  std::default_random_engine generator;
  std::normal_distribution<float> gaussian_noise;
  std::uniform_real_distribution<float> uniform_noise(0.0,1.0);


  // Do 5 NeuroBayes Trainings
  for(int k = 3; k < 5; k++) {
    ec_t* ec = nullptr;
    log_struct_t* ls = new log_struct_t;
    ls->str = nullptr;

    std::system("rm -rf test_dir");
    std::system("mkdir test_dir");

    NeuroBayesTeacher* nb;
    if(k ==0) {
      // First one with logging
      nb = NeuroBayesTeacher::Instance(&ec, normal_dbg, (log_func_t)(logf), (ls));
    } else {
      // Rest without
      nb =  NeuroBayesTeacher::Instance(&ec, quiet_dbg, NULL, NULL, NULL);
    }

    nb->SayHello();
    nb->NB_DEF(1);
    nb->NB_DEF_NODE1(nvar+1);
    nb->NB_DEF_NODE2(13);   
    nb->NB_DEF_NODE3(20);

    if(k == 1) {
      // Try invalid regularisation
      nb->NB_DEF_REG("INVALID_REGULARISATION");
    } else {
      nb->NB_DEF_REG("REG");
    }

    if(k == 2) {
      // Do classification
      nb->NB_DEF_TASK("CLA"); 
    } else {
      // Do density training
      nb->NB_DEF_TASK("DEN"); 
    }

    // Set common parameters
    nb->NB_DEF_PRE(32);
    nb->NB_DEF_SHAPE("DIA2");
    nb->NB_DEF_ITER( 0 );
    nb->SetIndividualPreproFlag(0, 15);
    for(int i = 1; i < nvar; ++i) {
      nb->SetIndividualPreproFlag(i, 14);
    }
    
    // Create and set input events
    nb->SetNEvents(nevents);
    float nb_input[nvar];
    for ( int i = 0; i< nevents; i++ ) {
      float target = k==0 ? 1./0. : gaussian_noise(generator); 
      for( int j=0; j<nvar; j++) {
        float u = uniform_noise(generator);
        float g = gaussian_noise(generator);
        if(j==0)
          nb_input[j] = target*u + g*(1-u);
        else 
          nb_input[j] = nb_input[j-1]*u + g*(1-u);
        if(k==4) input_matrix[i][j] = nb_input[j];
      }
      nb->SetTarget(k == 3 ? 0 : target);
      nb->SetNextInput(nvar, nb_input);
    }

    // Do Training and write out expertise and control information
    nb->SetOutputFile("./test_dir/density_test.nb");
    nb->TrainNet();
    nb->nb_correl_signi("./test_dir/correl_signi.txt", "./test_dir/correl_signi.html"); 

    // Check result of training
    bool ok = false;
    if (k == 0) {
      ok = (ec != nullptr) and (ec->kind == invalid_data_exc) and (std::string("NeuroBayes Teacher(R) \n") == std::string(ls->str));
    } else if(k == 1) {
      ok = (ec != nullptr) and (ec->kind == invalid_arg_exc);
    } else if(k == 2) {
      ok = (ec != nullptr) and (ec->kind == invalid_arg_exc);
    } else if(k == 3) {
      ok = (ec != nullptr) and (ec->kind == invalid_data_exc);
    } else {
      ok = (ec == nullptr);
    }
    
    if(not ok) {
      std::cerr << std::endl << "Run number " << k << " failed" << std::endl; 
      if(ec != nullptr)
        std::cerr << "EC: " << ec->kind << " " << ec->reason << " " << ec->stacktrace << std::endl;
      if(ls != nullptr and ls->str != nullptr)
        std::cerr << "LOG: " << ls->str << std::endl;
    }
      
    if(ec != nullptr)
      free_ec1(ec);
    if(ls != nullptr) {
      if(ls->str != nullptr)
        free(ls->str);
      delete ls;
    }

    testok = testok and ok;
  }
  return testok;
}

int expert() {
  int testok = true;
  for(int k=0;k<4;k++) {
    
    float tmean, trim;
    ec_t* ec = nullptr;
    Expert* nb;
    if(k == 0) {
      nb = new Expert(static_cast<const char*>(nullptr), quiet_dbg, false, &ec);
    } else {
      nb = new Expert("./test_dir/density_test.nb", quiet_dbg, false, &ec);
    }
    
    float nb_input[nvar];
    for ( int i = 0; i< nevents; i++ ) {
      for( int j=0; j<nvar; j++) {
	      nb_input[j] = input_matrix[i][j];
      }
      if (k==1) nb_input[0]=1./0. ;
      trim = k == 2 ? 15 : 0.05; 
      tmean = nb->nb_expert(nb_input, Expert::TRIM, trim);
    }
    if (k < 3) {
      testok &= (ec != nullptr);
      switch (k) {
        case 0: testok &= (ec->kind == file_open_exc);
                break;
        case 1:
                testok &= (ec->kind == invalid_data_exc); break;
        case 2:
                testok &= (ec->kind == invalid_arg_exc); break;
      }
      free_ec1(ec);
    } else {
      testok &= (ec == nullptr);
    }
  }
  return testok;
}

int main(int argc, char **argv) {
  
  std::cerr << argv[0] << std::endl;
  std::cerr << "Training ";
  
  if(training()) {
    std::cerr << " successfull" << std::endl;
    std::cerr << "Expert ";
  
    if(expert()) {
      std::cerr << " successfull" << std::endl;
      return 0;
    } else {
      std::cerr << " failed" << std::endl;
    }

  } else {
    std::cerr << " failed" << std::endl;
  }

  return 1;
}
