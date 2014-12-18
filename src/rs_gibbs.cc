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
                    if (transcript_cluster_map.find(sk.tids(i)) != transcript_cluster_map.end()) {
                        printf("Transcript-cluster (%s,%s) already exists in map (new cluster %s with %i transcripts); original cluster with %zu transcripts\n", sk.tids(i).c_str(), transcript_cluster_map[sk.tids(i)].c_str(), sk.gid().c_str(), sk.tids_size(), cluster_transcripts_map[transcript_cluster_map[sk.tids(i)]].size());
                    } else {
                        transcript_list.push_back(sk.tids(i));
                        transcript_cluster_map[sk.tids(i)] = sk.gid();
                    }
                }
                
                if (transcript_list.size() > 0) {
                    if (cluster_transcripts_map.find(sk.gid()) != cluster_transcripts_map.end()) {
                        printf("Cluster %s already in map\n", sk.gid().c_str());
                    }
                    cluster_transcripts_map[sk.gid()] = transcript_list;
                }
            }
            istream.close();
            num_transcripts = (int)transcript_cluster_map.size();
            num_clusters = (int)cluster_transcripts_map.size();
            if (DEBUG) printf("%i clusters, %i transcripts\n", num_clusters, num_transcripts);
            
            // For each replicate:
            
            // cluster_theta_map: each cluster mapped to a theta matrix
            for (map<string, vector<string>>::iterator it= cluster_transcripts_map.begin(); it!=cluster_transcripts_map.end(); it++) {
                vector<vector<int>> theta_matrix;
                theta_matrix.resize(num_replicates);
                for (int i = 0; i < (int)theta_matrix.size(); i++) {
                    theta_matrix[i].resize(it->second.size());
                }
                cluster_theta_map[it->first] = theta_matrix;
            }
            
            // Each theta matrix: num_replicates x number of transcripts
            vector<int> initial_G;
            for (int i = 0; i < num_replicates; i++) {
                em_file = em_file_prefix + std::to_string(i + 1) + "_em";
                if (DEBUG) printf("Reading file: %s\n", em_file.c_str());
                fstream istream(em_file, ios::in);
                int max_theta = 0, max_transcript = 0;
                while (getline(istream, line)) {
                    vector<string> tokens = split(line, '\t');
                    // if (DEBUG) printf("Prepping line: %s\n", line.c_str());
                    int occ = atoi(tokens[2].c_str());
                    string transcript = tokens[0].c_str();
                    string cluster = transcript_cluster_map[transcript].c_str();
                
                    vector<string> transcripts = cluster_transcripts_map[cluster];
                    // Find the index of the transcript in the cluster
                    int j = find(transcripts.begin(), transcripts.end(), transcript) - transcripts.begin();
                    if (j >= (int)transcripts.size()) {
                        fprintf(stderr, "ERROR: failed to find index for %s in cluster %s\n", transcript.c_str(), cluster.c_str());
                    }
                    
                    vector<vector<int>> theta_matrix = cluster_theta_map[cluster];
                    theta_matrix[i][j] = occ;
                    cluster_theta_map[cluster] = theta_matrix;
                    
                    // if (DEBUG) printf("Set (%i, %i) on cluster %s for transcript %s to %i\n", i, j, cluster.c_str(), transcript.c_str(), occ);
                    
                    if (occ > max_theta) {
                        max_theta = occ;
                        max_transcript = j;
                    }
                }
                initial_G.push_back(max_transcript); // each replicate has same max theta value for transcripts
            }
            
            // Create set of F matrices per cluster: number of sigmers x number of transcripts
            // F matrix should be same for all replicates -- only analyze 1 replicate file
            cf_file = em_file_prefix + std::to_string(1) + ".cf";
            if (DEBUG) printf("Reading file for F matrix: %s\n", cf_file.c_str());
            fstream istream2(cf_file, ios::in | ios::binary);
            // each sk is a cluster of sigmers
            // sk.keys is the group of sigmers associated with each cluster
            // sk.keys(j).transcript_infos is the group of transcripts each sk.keys(j) sigmer is assoc with
            // sk.keys(j).transcript_infos(k).tidx is the transcript identifier
            // sk.keys(j).transcript_infos(k).positions_size() is the number of times a sigmer appears in the transcript
            while(load_protobuf_data(&istream2, &sk, buffer, buffer_size)) {
                vector<vector<int>> Fmatrix;
                string cluster = sk.gid().c_str();
                
                // iterate over sigmers
                for (int i = 0; i < sk.keys_size(); i++) {
                    // create vector -- size = num_transcripts per cluster, each entry is number of times a sigmer appears in the transcript
                    vector<int> transcript_vector;
                    transcript_vector.resize(cluster_transcripts_map[cluster].size());
                    
                    string sigmer = sk.keys(i).key();
                    // iterate over the transcripts associated with the sigmers
                    for (int j = 0; j < sk.keys(i).transcript_infos_size(); j++) {
                        int tidx = sk.keys(i).transcript_infos(j).tidx();
                        string transcript = sk.tids(tidx); // TODO: Verify this?
                        
                        if (DEBUG) { // assertion: every cluster assoc with transcripts should be same
                            string test_cluster = transcript_cluster_map[transcript];
                            assert(strcmp(cluster.c_str(), test_cluster.c_str()) == 0);
                        }
                        
                        // Find the index of the transcript in the cluster
                        vector<string> transcripts = cluster_transcripts_map[cluster];
                        int k = find(transcripts.begin(), transcripts.end(), transcript) - transcripts.begin();
                        if (k >= (int)transcripts.size()) {
                            fprintf(stderr, "ERROR: failed to find index for %s in cluster %s\n", transcript.c_str(), cluster.c_str());
                        }
                        
                        transcript_vector[k] = sk.keys(i).transcript_infos(j).positions_size();
                    }
                    Fmatrix.push_back(transcript_vector);
                    
                    // TODO: set the cluster_G_map[cluster] = vector of G maps
                }
                cluster_Fmatrix_map[cluster] = Fmatrix;
                //if (DEBUG) printf("Fmatrix for cluster %s: %zu sigmers x %zu transcripts\n", cluster.c_str(), Fmatrix.size(), Fmatrix[0].size());
                
                // Create r vectors for sigmer -> transcript id from initial_G
                vector<vector<int>> Gvector;
                for (int i = 0; i < num_replicates; i++) {
                    vector<int> G;
                    G.resize(sk.keys_size(), initial_G[i]);
                    Gvector.push_back(G);
                }
                cluster_G_map[cluster] = Gvector;
            }
            if (DEBUG) printf("cluster_Fmatrix_map number of clusters: %zu\n", cluster_Fmatrix_map.size());
            
            // Create L vector per cluster
            for (map<string, vector<vector<int>>>::iterator it= cluster_Fmatrix_map.begin(); it!=cluster_Fmatrix_map.end(); it++) {
                string cluster = it->first;
                vector<vector<int>> F = it->second;
                vector<int> L;
                L.resize(F[0].size(), 0); // length = number of transcripts (row size)

                for (int i = 0; i < (int)F.size(); i++) {
                    for (int j = 0; j < (int)F[i].size(); j++) {
                        L[j] += F[i][j];
                    }
                }
                
                cluster_L_map[cluster] = L;
            }
            
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
        
        map<string, vector<string>> cluster_transcripts_map; // cluster -> list of associated transcripts
        map<string, string> transcript_cluster_map; // transcript -> corresponding cluster
        map<string, vector<vector<int>>> cluster_theta_map; // cluster -> theta matrix of occurrences of sigmers per transcript
        map<string, vector<vector<int>>> cluster_Fmatrix_map; // cluster -> F matrix: sigmers x transcripts (num of sigmers per transcript)
        map<string, vector<int>> cluster_L_map; // cluster -> L vector
        map<string, vector<vector<int>>> cluster_G_map; // cluster -> #replicates x #sigmers -> index of transcript
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}
