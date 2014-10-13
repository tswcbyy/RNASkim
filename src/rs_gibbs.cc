#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

#include "gflags/gflags.h"

#include "proto/rnasigs.pb.h"

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::stringstream;

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
                string em_file = em_file_prefix_ + std::to_string(i+1) + "_em";
                fstream istream(em_file, ios::in);
                string line;
                int t = 0;
                while (getline(istream, line)) {
                    vector<string> tokens = split(line, '\t');
                    int occ = atoi(tokens[2].c_str());
                    theta[i][t] = occ;
                    ++t;
                }
                num_trans = t;
            }
            for (int i = 0; i < num_replicates_; ++i) {
                cout << i << endl;
                for (int t = 0; i < num_trans; ++t) {
                    cout << theta[i][t] << " ";
                }
                cout << endl;
            }
        }
        
        vector<string> split(const string &s, char delim) {
            vector<string> elems;
            stringstream ss(s);
            string item;
            while (getline(ss, item, delim)) {
                elems.push_back(item);
            }
            return elems;
        }
        
        void run(){
            cout << "done" << endl;
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
    rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}