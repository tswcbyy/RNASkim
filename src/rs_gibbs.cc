#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "gflags/gflags.h"

#include "proto/rnasigs.pb.h"

using std::string;
using std::vector;
using std::map;

DEFINE_int32(num_replicates, 4,
             "The number of replicates in this condition.");
DEFINE_string(em_file_prefix, "",
              "The prefix of the path to the file that contains the EM results.");

namespace rs {
    class GibbsSampler {
    public:
        GibbsSampler(const int num_repicates, const string em_file_prefix){
            num_replicates_ = num_repicates;
            em_file_prefix_ = em_file_prefix;
            for (int i = 0; i < num_replicates_; ++i) {
                string em_file = em_file_prefix_ + std::to_string(i) + "_em";
                std::cout << em_file << std::endl;
            }
        }
    private:
        int num_replicates_; //R_c
        string em_file_prefix_;
        int num_trans; //T
        vector<int> num_sigmers; //L
        
        vector<double> mean; //m
        vector<vector<double>> theta;
        map<SelectedKey_Key, int> parent; //G
        
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file);
    gs.run();
}