#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "NeuroBayesTeacher.hh" //NeuroBayes Header
#include "NeuroBayesExpert.hh"  //NeuroBayes Header

using namespace std;

float signal_distribution(int) {
    float r = (float)rand() / (float)RAND_MAX;
    return r;
}

float bckgrd_distribution(int) {
    float r = (float)rand() / (float)RAND_MAX + 0.1;
    return r;
}

int main() {
    const int  nvar = 2;    //number of input variables
    float input[nvar];
    int generator = 0;
    NeuroBayesTeacher *nb = NeuroBayesTeacher::Instance();
    nb->NB_DEF(true);

    nb->SetOutputFile("weightfile.nb");     // expert file

    nb->NB_DEF_NODE1(nvar+1);       // nodes in input layer 
    nb->NB_DEF_NODE2(nvar+2);           // nodes in hidden layer 
    nb->NB_DEF_NODE3(1);        // nodes in output layer

    nb->NB_DEF_TASK("CLA");     // binominal classification

    int i= 4701;
    int j=21; 
                // i has to be an odd number, the third argument is a debugging flag
    nb->NB_RANVIN(i,j,2);       // random number seed initialisation,
 
    nb->NB_DEF_PRE(112);                        // Global Preprocessing Flag 
    nb->NB_DEF_REG("ARD");              // 'OFF','REG' (def) ,'ARD','ASR','ALL'
    nb->NB_DEF_LOSS("ENTROPY");                 // 'ENTROPY'(def),'QUADRATIC'
    nb->NB_DEF_SHAPE("OFF");                    // 'OFF', 'INCL', 'TOTL'

    nb->NB_DEF_EPOCH(10);               // weight update after n events
    nb->NB_DEF_MOM(0.01);                       // Momentum 

    nb->NB_DEF_SPEED(1.0);              // multiplicative factor to enhance global learning speed
    nb->NB_DEF_MAXLEARN(0.9);           // multiplicative factor to limit the global learning speed
    nb->NB_DEF_ITER(100);               // number of training iteration
    nb->NB_DEF_METHOD("BFGS");                  // Training Method
                //nb->NB_DEF_INITIALPRUNE(fPruning);
                //nb->NB_DEF_RTRAIN(fTrainTestRatio);           // Ratio of Events to use for Trainig, Rest is used for Testing

    for(unsigned int j = 0; j < nvar; ++j) {
      nb->SetIndividualPreproFlag(j,34); //Set Individual Preprocessing Flag
    }

    std::cout << "Start NeuroBayes Training" << std::endl;
    for(unsigned int i = 0; i < 1000; ++i) {
      bool target = i % 2 == 0;
      for(unsigned int j = 0; j < nvar; ++j) {
        if(target)
          input[j] = signal_distribution(generator);
        else
          input[j] = bckgrd_distribution(generator);
      }
      nb->SetWeight(1.0);
      nb->SetTarget(target); 
      nb->SetNextInput(nvar, input);
    }
    char *c_varnames[nvar];
    for(unsigned int j = 0; j < nvar; ++j) {
      c_varnames[j] = new char[10];
      sprintf(c_varnames[j],"%d",j);
    }
    nb->TrainNet();
    nb->nb_correl_signi(c_varnames,"./correl_signi.txt","./correl_signi.html");

    std::cout << "Start evaluation using Expert" << std::endl;
    Expert* net = new Expert("weightfile.nb");

    double loss = 0;
    for(unsigned int i = 0; i < 1000; ++i) {
      bool target = i % 2 == 0;
      for(unsigned int j = 0; j < nvar; ++j) {
        if(target)
          input[j] = signal_distribution(generator);
        else
          input[j] = bckgrd_distribution(generator);
      }
            float prob = net->nb_expert(input);
      float t = (target) ? 1 : -1;
      loss += (t - prob) * (t - prob);
    }
    std::cout << "L2 " << loss / 1000 << std::endl;

  std::cout << "Execute this to produce analysis.pdf" << std::endl;
  std::cout << "root -b -q $NEUROBAYES/external/analysis.C'(\"ahist.txt\",\"analysis.pdf\",1,\"correl_signi.txt\")" << std::endl;;

}

