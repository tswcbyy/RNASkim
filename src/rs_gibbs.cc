#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <assert.h>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "rs_common.h"
#include "rs_estimate_lib.h"
#include "proto/rnasigs.pb.h"
#include "proto_data.h"

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
#define DEBUG 1

namespace rs {
    class GibbsSampler {
    public:
        
        GibbsSampler(const int num_replicates, const string em_file_prefix){
            if (DEBUG) printf("File prefix: %s\n", em_file_prefix.c_str());
            
            string em_file, line, cf_file;
            int num_transcripts, num_clusters;
            
            // Create mapping for each cluster -> transcript IDs
            // Map transcript IDs to specific columns in the theta and F matrices(?)
            // TODO: discard less useful map(?)
            map<string, vector<string>> cluster_transcripts_map;
            map<string, string> transcript_cluster_map;
            cf_file = em_file_prefix + std::to_string(1) + ".cf";
            if (DEBUG) printf("Reading file: %s\n", cf_file.c_str());
            fstream istream(cf_file, ios::in | ios::binary);
            
            SelectedKey sk;
            int buffer_size = 200000000;
            ::google::protobuf::uint8 * buffer =
            new ::google::protobuf::uint8[buffer_size];
            // each sk is a cluster of sigmers
            // sk.tids are all the transcripts associated with the cluster
            while(load_protobuf_data(&istream, &sk, buffer, buffer_size)) {
                vector<string> transcript_list;
                for (int i = 0; i < sk.tids_size(); i++) {
                    transcript_list.push_back(sk.tids(i));
                    transcript_cluster_map[sk.tids(i)] = sk.gid();
                }
                cluster_transcripts_map[sk.gid()] = transcript_list;
            }
            num_transcripts = (int)transcript_cluster_map.size();
            num_clusters = (int)cluster_transcripts_map.size();
            if (DEBUG) printf("%i clusters, %i transcripts\n", num_clusters, num_transcripts);
            
            // For each replicate:
            
            // TODO: separate by cluster
            map<string, vector<vector<int>>> cluster_theta_map;
            // Each theta matrix: num_replicates x number of transcripts
            for (int i = 0; i < num_replicates; i++) {
                vector<int> elem;
                em_file = em_file_prefix + std::to_string(i + 1) + "_em";
                if (DEBUG) printf("Reading file: %s\t", em_file.c_str());
                fstream istream(em_file, ios::in);
                while (getline(istream, line)) {
                    vector<string> tokens = split(line, '\t');
                    int occ = atoi(tokens[2].c_str());
                    elem.push_back(occ);
                }
                if (DEBUG) printf("Row %i size: %zu\n", i, elem.size());
                theta.push_back(elem);
            }
            
            // Theta error check: each replicate has same number of transcripts
            for (int i = 0; i < (int)theta.size(); i++) {
                assert((int)theta[i].size() == num_transcripts);
            }
            if (DEBUG) printf("Theta size: %zu replicates x %zu transcripts\n", theta.size(), theta[0].size());
            
            // Create F matrix: number of sigmers x number of transcripts
            for (int i = 0; i < num_replicates; i++) {
                vector<int> transcript_list;
                cf_file = em_file_prefix + std::to_string(i + 1) + ".cf";
                if (DEBUG) printf("Reading file: %s\t", cf_file.c_str());
                fstream istream(cf_file, ios::in | ios::binary);
                
                SelectedKey sk;
                int buffer_size = 200000000;
                ::google::protobuf::uint8 * buffer =
                new ::google::protobuf::uint8[buffer_size];
                // each sk is a cluster of sigmers
                // sk.keys is the group of sigmers associated with each cluster
                // sk.keys(j).transcript_infos is the group of transcripts each sk.keys(j) sigmer is assoc with
                // sk.keys(j).transcript_infos(k).tidx is the transcript identifier
                // sk.keys(j).transcript_infos(k).positions_size() is the number of times a sigmer appears in the transcript
                while(load_protobuf_data(&istream, &sk, buffer, buffer_size)) {
                    //printf("Number of sigmers: %i for cluster %s\n", sk.keys_size(), sk.gid().c_str());
                }
            }
            
            
            // Create L vector
            
            // G per n sigmer occurrences:
            // Initialize n values of G for max t: P(G = t | ...) = theta_t * M_sigmer,t / L_t
            
            /*
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
        vector<vector<int>> F;
        map<SelectedKey_Key, int> parent; //G
        
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}
