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
        GibbsSampler(const int num_replicates, const string em_file_prefix){
            num_replicates_ = num_replicates;
            em_file_prefix_ = em_file_prefix;
            
            // For each replicate:
            // Create theta matrix: num_replicates x number of transcripts
            // Create M matrix
            // Create L vector
            
            // G per n sigmer occurrences:
            // Initialize n values of G for max t: P(G = t | ...) = theta_t * M_sigmer,t / L_t
            
            printf("Current theta size: %i x %i\n", theta.size(), theta[0].size());
            theta.resize(num_replicates);
            for (int i = 1; i < num_replicates_; i++) {
                
            }
            
            /*
            //read theta from EM results
            theta.resize(num_replicates_);
            vector<int> elem;
            string em_file = em_file_prefix_ + "1_em";
            fstream istream(em_file, ios::in);
            string line;
            while (getline(istream, line)) {
                vector<string> tokens = split(line, '\t');
                int occ = atoi(tokens[2].c_str());
                elem.push_back(occ);
            }
            num_trans = elem.size();
            theta[0] = elem;
            
            //read theta from EM results (cont)
            for (int i = 1; i < num_replicates_; ++i) {
                vector<int> elem(num_trans);
                string em_file = em_file_prefix_ + std::to_string(i+1) + "_em";
                fstream istream(em_file, ios::in);
                string line;
                int t = 0;
                while (getline(istream, line)) {
                    vector<string> tokens = split(line, '\t');
                    int occ = atoi(tokens[2].c_str());
                    elem[t] = occ;
                    ++t;
                }
                theta[i] = elem;
            }
            
            //map from the order of transcripts to their tids
            order2tid.resize(num_trans);
            fstream is(em_file, ios::in);
            int t = 0;
            while (getline(is, line)) {
                vector<string> tokens = split(line, '\t');
                string tid = tokens[0];
                order2tid[t] = tid;
                ++t;
            }
            
            //check the correntness of order2tid mapping
            for (int t = 0; t < num_trans; ++t) {
                cout << order2tid[t] << " ";
            }
            cout << endl;
            
            // check the correctness of reading theta from EM results
            for (int i = 0; i < num_replicates_; ++i) {
                cout << i << endl;
                cout << theta[i].size() << endl;
                for (int t = 0; t < num_trans; ++t) {
                    cout << theta[i][t] << " ";
                }
                cout << endl;
            }
             */
            
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
        vector<string> order2tid;
        vector<double> mean; //m
        vector<vector<int>> theta;
        map<SelectedKey_Key, int> parent; //G
        
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}
