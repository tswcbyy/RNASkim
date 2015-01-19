#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <assert.h>
#include <random>
#include <math.h>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "rs_common.h"
#include "rs_estimate_lib.h"
#include "proto/rnasigs.pb.h"
#include "proto_data.h"

#include <boost/random/discrete_distribution.hpp>

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::fstream;
using std::ifstream;
using std::ios;
using std::stringstream;
using std::find;

DEFINE_int32(num_replicates, 4,
             "The number of replicates in this condition.");
DEFINE_string(em_file_prefix, "",
              "The prefix of the path to the file that contains the EM results.");
DEFINE_string(input_file, "", "Text file containing all input files (tab-delimited) with 3 columns: condition name, .cf file, _em file");

#define DEBUG 1

namespace rs {
    class GibbsSampler {
    public:
        
        void GibbsSamplerOld(const int num_replicates, const string em_file_prefix){
            /*
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
            for (auto it : cluster_transcripts_map) {
                vector<vector<int>> theta_matrix;
                theta_matrix.resize(num_replicates);
                for (int i = 0; i < (int)theta_matrix.size(); i++) {
                    theta_matrix[i].resize(it.second.size());
                }
                cluster_theta_map[it.first] = theta_matrix;
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
                vector<string> transcripts = cluster_transcripts_map[cluster];
                //if (DEBUG) printf("Analyzing cluster %s for %i sigmers, with %zu transcripts\n", cluster.c_str(), sk.keys_size(), transcripts.size());
                if (transcripts.size() == 0) {
                    if (DEBUG) printf("SKIP: cluster %s\n", cluster.c_str());
                    continue;
                }

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
                        int k = find(transcripts.begin(), transcripts.end(), transcript) - transcripts.begin(); // TODO: fix efficiency
                        if (k >= (int)transcripts.size()) {
                            fprintf(stderr, "ERROR: failed to find index for %s in cluster %s\n", transcript.c_str(), cluster.c_str());
                        }

                        transcript_vector[k] = sk.keys(i).transcript_infos(j).positions_size();
                    }
                    Fmatrix.push_back(transcript_vector);
                }
                cluster_Fmatrix_map[cluster] = Fmatrix;
                //if (DEBUG) printf("Fmatrix for cluster %s: %zu sigmers x %zu transcripts\n", cluster.c_str(), Fmatrix.size(), Fmatrix[0].size());

                // Create r vectors for sigmer -> transcript id from initial_G
                // G per n sigmer occurrences:
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
            for (auto it : cluster_Fmatrix_map) {
                string cluster = it.first;
                vector<vector<int>> F = it.second;
                vector<string> transcripts = cluster_transcripts_map[cluster];
                // if (DEBUG) printf("F matrix %s: %zu\n", cluster.c_str(), F.size());

                vector<int> L;
                L.resize(F[0].size(), 0); // length = number of transcripts (row size)

                for (int i = 0; i < (int)F.size(); i++) {
                    for (int j = 0; j < (int)F[i].size(); j++) {
                        L[j] += F[i][j];
                    }
                }
                for (int i = 0; i < (int)L.size(); i++) {
                    if (L[i] == 0) {
                        printf("WARNING: transcript %s : %s has 0 sigmers\n", cluster.c_str(), transcripts[i].c_str());
                    }
                }

                cluster_L_map[cluster] = L;
            }*/

            /*

             map<string, vector<string>> cluster_transcripts_map; // cluster -> list of associated transcripts
             map<string, string> transcript_cluster_map; // transcript -> corresponding cluster
             map<string, vector<vector<int>>> cluster_theta_map; // cluster -> theta matrix of occurrences of sigmers per transcript
             map<string, vector<vector<int>>> cluster_Fmatrix_map; // cluster -> F matrix: sigmers x transcripts (num of sigmers per transcript)
             map<string, vector<int>> cluster_L_map; // cluster -> L vector, length = number of transcripts (#sigmers per transcript)
             map<string, vector<vector<int>>> cluster_G_map; // cluster -> #replicates x #sigmers -> index of transcript

            */


            /* Gibbs:
             1) For each cluster c, for each replicate r, for each sigmer occurrence s: re-sample G (transcript)
             P(G = t | ...) = theta_t * F_s,t / L_t

             2) theta matrix:
             For each replicate cluster c, for each replicate r:
             P(theta_rt | ...) = P(theta_r,t | m_t) * (theta_r,t / (sum_(i != t) theta_r,i + theta_r,t))^(#G == t)

             3) For each cluster c, produce m vector (length = num transcripts)
             For each transcript t:
             P(m_t| ...) = product_r (P(theta_r,t | m_t))
             P(theta_r,t | m_t) = choose(theta_r,t + q - 1, theta_r,t) (1-p)^q * p^(theta_rt)
             Normalize m vector; output

             */

            std::default_random_engine generator;
            // for each cluster:
            for (auto it : cluster_transcripts_map) {
                string cluster = it.first;
                vector<string> transcripts = it.second;
                int num_transcripts_cluster = (int) transcripts.size();
                if (DEBUG) printf("Gibbs sampling for cluster: %s [%i]\n", cluster.c_str(), num_transcripts_cluster);

                // initializations
                vector<vector<int>> G;// = cluster_G_map[cluster]; // #replicates x #sigmers
                vector<vector<int>> theta;// = cluster_theta_map[cluster]; // #replicates x #transcripts
                vector<vector<int>> F;// = cluster_Fmatrix_map[cluster]; // #sigmers x #transcripts
                vector<int> L;// = cluster_L_map[cluster]; // #transcripts
                vector<vector<double>> m; // 1000 x normalized vector #transcripts

                if (DEBUG) printf("G: %zu rows\n", G.size());

                for (int i = 0; i < 1000; i++) {
                    // Resample G
                    if (DEBUG) printf("G iteration %i\n", i);
                    for (int r = 0; r < (int)G.size(); r++) {
                        for (int s = 0; s < (int)G[0].size(); s++) {
                            vector<double> G_probs;
                            // choose random t from distribution over all t: theta[r][t] * F[s][t] / L[t]
                            for (int t = 0; t < num_transcripts_cluster; t++) {
                                if (L[t] == 0)
                                    G_probs.push_back(0);
                                else
                                    G_probs.push_back(theta[r][t] * F[s][t] / L[t]);
                            }
                            boost::random::discrete_distribution<> distribution(G_probs.begin(), G_probs.end());
                            G[r][s] = distribution(generator);
                        }
                    }

                    // Resample theta
                    for (int r = 0; r < (int)theta.size(); r++) {
                        for (int t = 0; t < (int)theta[0].size(); t++) {
                        }
                    }

                    // Resample m and output
                    vector<double> m_iter(num_transcripts_cluster);
                    m.push_back(m_iter);
                }

                // TODO:
                vector<double> m_mean(num_transcripts_cluster); // mean m across 1000 iterations x #transcripts
                vector<double> m_25(num_transcripts_cluster); // 25th m x #transcripts
                vector<double> m_975(num_transcripts_cluster); // 975th m x #transcripts
                vector<double> m_var(num_transcripts_cluster); // variance of m x #transcripts
                
                /*for (int i = 0; i < num_transcripts_cluster; i++) {
                    fprintf(output, "%s\t%f\t%f\t%f\t%f\n", transcripts[i].c_str(), m_mean[i] * size_condition, m_var[i] * size_condition * size_condition, m_25[i], m_975[i]);
                }*/
            }
            
        }

        GibbsSampler(const string input_file) {
            // conditions, replicates, cf_files, em_files
            parseInputFile(input_file);
            if (DEBUG) {
                printf("Checking input file parsing ...");
                for (auto it : conditions) {
                    string c = it.first;
                    int idx = it.second;
                    printf("Condition: %s [r = %i]\n", c.c_str(), replicates[idx]);
                    for (string f : cf_files[idx]) {
                        printf("\tcf file: %s\n", f.c_str());
                    }
                    for (string f : em_files[idx]) {
                        printf("\tem file: %s\n", f.c_str());
                    }
                }

                size_t num_conditions = conditions.size();
                assert(num_conditions > 0);
                assert(num_conditions == replicates.size());
                assert(num_conditions == cf_files.size());
                assert(num_conditions == em_files.size());
                for (auto it : conditions) {
                    int c = it.second;
                    int r = replicates[c];
                    assert(r > 0);
                    assert(r == (int) cf_files[c].size());
                    assert(r == (int) em_files[c].size());
                }
                printf(" complete\n");
            }

            // cluster_transcripts_map, transcript_cluster_map, clusters, transcripts
            parseClusterTranscripts(cf_files[0][0]);
            if (DEBUG) {
                printf("%zu clusters, %zu transcripts\n", clusters.size(), transcripts.size());
                printf("Checking cluster-transcript parsing ...");
                assert(clusters.size() == cluster_transcripts_map.size());
                assert(transcripts.size() == transcript_cluster_map.size());
                int i = 0;
                for (auto it : transcripts) {
                    string t = it.first;
                    int idx = it.second;
                    string c = transcript_cluster_map[t];
                    vector<string> ts = cluster_transcripts_map[c];
                    assert(strcmp(t.c_str(), ts[idx].c_str()) == 0);
                    i++;
                    if (i%1000 == 0) printf(".");
                }
                printf(" complete\n");
            }

            // theta_data
            parseEMData();
            if (DEBUG) {
                printf("EM data: complete\n");
            }

            // F_data, G_data, L_data
            parseCFData(cf_files[0][0]);
            if (DEBUG) {
                printf("CF data: complete\n");
            }

            // size_replicates, size_conditions
            calcSizeFactors();
            if (DEBUG) {
                printf("Size factors: complete\n");
            }

            /* Gibbs:
             1) For each cluster c, for each replicate r, for each sigmer occurrence s: re-sample G (transcript)
             P(G = t | ...) = theta_t * F_s,t / L_t

             2) theta matrix:
             For each replicate cluster c, for each replicate r:
             P(theta_rt | ...) = P(theta_r,t | m_t) * (theta_r,t / (sum_(i != t) theta_r,i + theta_r,t))^(#G == t)

             3) For each cluster c, produce m vector (length = num transcripts)
             For each transcript t:
             P(m_t| ...) = product_r (P(theta_r,t | m_t))
             P(theta_r,t | m_t) = choose(theta_r,t + q - 1, theta_r,t) (1-p)^q * p^(theta_rt)
             Normalize m vector; output
             
             */

            if (DEBUG) printf("Gibbs sampling:\n");
            for (auto cond_it : conditions) {
                int cond_idx = cond_it.second;
                string condition = cond_it.first;
                if (DEBUG) printf("\tCondition: %s\n", condition.c_str());

                for (int i = 0; i < 1000; i++) {
                    if (DEBUG && i%20 == 0) printf(".");
                    for (auto cl_it : clusters) {
                        int cl_idx = cl_it.second;
                        for (int r = 0; r < replicates[cond_idx]; r++) {
                            gibbsG(cond_idx, r, cl_idx);
                            gibbsTheta(cond_idx, r, cl_idx);
                        }
                        gibbsM(cond_idx, cl_idx);
                    }

                    output();
                }

                if (DEBUG) printf("\tcomplete\n");
            }
        }

        // Initializes: conditions (map), replicates, cf_files, em_files
        void parseInputFile(const string input_file) {
            string line;
            ifstream ifile(input_file);
            if (ifile.is_open()) {
                while (getline(ifile, line)) {
                    stringstream linestream(line);
                    string condition;
                    string cf_file;
                    string em_file;

                    getline(linestream, condition, '\t');
                    linestream >> cf_file >> em_file;

                    auto c = conditions.find(condition);
                    if (c == conditions.end()) { // new condition
                        conditions[condition] = (int)conditions.size() - 1;
                        replicates.push_back(1);

                        vector<string> cf;
                        cf.push_back(cf_file);
                        cf_files.push_back(cf);

                        vector<string> em;
                        em.push_back(em_file);
                        em_files.push_back(em);
                    } else {
                        int c_idx = conditions[condition];
                        replicates[c_idx]++;

                        cf_files[c_idx].push_back(cf_file);
                        em_files[c_idx].push_back(em_file);
                    }
                }
                ifile.close();
            } else {
                fprintf(stderr, "Unable to open %s\n", input_file.c_str());
            }
        }

        void parseClusterTranscripts(const string cf) {
            if (DEBUG) printf("parseClusterTranscript: %s\n", cf.c_str());

            fstream istream(cf, ios::in | ios::binary);
            SelectedKey sk;
            int buffer_size = 200000000;
            ::google::protobuf::uint8 * buffer =
            new ::google::protobuf::uint8[buffer_size];
            // each sk is a cluster of sigmers
            // sk.tids are all the transcripts associated with the cluster
            while(load_protobuf_data(&istream, &sk, buffer, buffer_size)) {
                vector<string> transcript_list;
                string cluster = sk.gid();
                for (int i = 0; i < sk.tids_size(); i++) {
                    string transcript = sk.tids(i);
                    if (transcript_cluster_map.find(transcript) != transcript_cluster_map.end()) {
                        printf("WARNING: Transcript-cluster (%s,%s) already exists in map (new cluster %s with %i transcripts); original cluster with %zu transcripts\n", sk.tids(i).c_str(), transcript_cluster_map[transcript].c_str(), cluster.c_str(), sk.tids_size(), cluster_transcripts_map[transcript_cluster_map[transcript]].size());
                    } else {
                        transcript_list.push_back(transcript);
                        transcript_cluster_map[transcript] = cluster;
                        transcripts[transcript] = transcript_list.size() - 1;
                    }
                }

                if (transcript_list.size() > 0) {
                    if (cluster_transcripts_map.find(sk.gid()) != cluster_transcripts_map.end()) {
                        printf("Cluster %s already in map\n", cluster.c_str());
                    } else {
                        clusters[cluster] = clusters.size() - 1;
                    }
                    cluster_transcripts_map[sk.gid()] = transcript_list;
                }
            }
            istream.close();
        }

        void parseEMData() {
            // theta_data: #conditions x #clusters x #replicates x #transcripts
            theta_data.resize(conditions.size());
            for (auto cond_it : conditions) {
                int cond_idx = cond_it.second;
                vector<vector<vector<int>>> theta_condition_data(clusters.size());
                for (auto cluster_it : clusters) {
                    string cluster = cluster_it.first;
                    int cluster_idx = cluster_it.second;
                    int num_replicates = replicates[cond_idx];
                    vector<vector<int>> theta_cluster_data(num_replicates);
                    for (int r = 0; r < num_replicates; r++) {
                        size_t num_transcripts = cluster_transcripts_map[cluster].size();
                        vector<int> theta_replicate_data(num_transcripts, 0);
                        theta_cluster_data[r] = theta_replicate_data;
                    }
                    theta_condition_data[cluster_idx] = theta_cluster_data;
                }
                theta_data[cond_idx] = theta_condition_data;
            }

            for (int cond_idx = 0; cond_idx < (int)em_files.size(); cond_idx++) {
                vector<string> files = em_files[cond_idx];

                for (int r_idx = 0; r_idx < (int)files.size(); r_idx++) {
                    string line, file = files[r_idx];
                    if (DEBUG) printf("Reading file: %s [condition %i, replicate %i]\n", file.c_str(), cond_idx, r_idx);
                    fstream istream(file, ios::in);

                    while (getline(istream, line)) {
                        vector<string> tokens = split(line, '\t');
                        int occ = atoi(tokens[2].c_str());
                        string transcript = tokens[0].c_str();
                        string cluster = transcript_cluster_map[transcript].c_str();

                        int t_idx = transcripts[transcript];
                        int cl_idx = clusters[cluster];

                        theta_data[cond_idx][cl_idx][r_idx][t_idx] = occ;
                    }
                    istream.close();
                }
            }
        }

        void parseCFData(const string cf) {
            /* F_data; #clusters x #sigmers x #transcripts: num (specific) sigmer per transcript
             * L_data; #clusters x #transcripts: num (general) sigmers per transcript
             * G_data; #conditions x #clusters x #replicates x #sigmers: transcript id (from transcripts) */

            G_data.resize(conditions.size());
            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<vector<vector<int>>> G_condition_data(clusters.size());
                for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                    vector<vector<int>> G_cluster_data(replicates[cond_idx]);
                    G_condition_data[cl_idx] = G_cluster_data;
                }
                G_data[cond_idx] = G_condition_data;
            }

            L_data.resize(clusters.size());
            F_data.resize(clusters.size());

            if (DEBUG) printf("Reading file for F matrix: %s\n", cf.c_str());
            fstream istream(cf, ios::in | ios::binary);
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
                string cluster = sk.gid().c_str();

                auto c = clusters.find(cluster);
                if (c == clusters.end()) {
                    if (DEBUG) printf("SKIP: cluster %s\n", cluster.c_str());
                    continue;
                }

                size_t num_transcripts = cluster_transcripts_map[cluster].size();
                int cl_idx = clusters[cluster];

                vector<vector<int>> F_cluster_data;
                vector<int> L_cluster_data(num_transcripts, 0);

                // iterate over sigmers
                for (int s = 0; s < sk.keys_size(); s++) {
                    // create vector -- size = num_transcripts per cluster, each entry is number of times a sigmer appears in the transcript
                    vector<int> F_sigmer_data(num_transcripts, 0);

                    string sigmer = sk.keys(s).key();
                    // iterate over the transcripts associated with the sigmers
                    for (int t = 0; t < sk.keys(s).transcript_infos_size(); t++) {
                        int tidx = sk.keys(s).transcript_infos(t).tidx();
                        string transcript = sk.tids(tidx); // TODO: Verify this?

                        int t_idx = transcripts[transcript];
                        int sigmer_count = sk.keys(s).transcript_infos(t).positions_size();
                        F_sigmer_data[t_idx] = sigmer_count;
                        L_cluster_data[t_idx] += sigmer_count;
                    }
                    F_cluster_data.push_back(F_sigmer_data);
                }

                F_data[cl_idx] = F_cluster_data;
                L_data[cl_idx] = L_cluster_data;

                // populate initial G_data
                for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        vector<int> theta_replicate_data = theta_data[cond_idx][cl_idx][r];
                        auto max_idx = std::max_element(theta_replicate_data.begin(), theta_replicate_data.end());
                        int g = std::distance(theta_replicate_data.begin(), max_idx);
                        vector<int> G_sigmer_data(sk.keys_size(), g);
                        G_data[cond_idx][cl_idx][r] = G_sigmer_data;
                    }
                }
            }
        }

        void calcSizeFactors() {
            int num_replicates = std::accumulate(replicates.begin(), replicates.end(), 0);
            int num_sigmers = 0;
            // size_replicate: median (over C clusters) for (#sigmers in cluster c for replicate r / sum #sigmers across all replicates in cluster c ^ 1/#replicates

            // calculate #sigmers per cluster per replicate
            vector<vector<vector<int>>> sigmer_data(conditions.size()); // #conditions x #replicates x #clusters: #sigmers (from theta_data)
            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<vector<int>> sigmer_condition(replicates[cond_idx]);

                for (int r = 0; r < replicates[cond_idx]; r++) {
                    // #clusters: #sigmers that occur per cluster in this replicate
                    vector<int> sigmer_replicate(clusters.size());

                    for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                        vector<int> sigmer_cluster = theta_data[cond_idx][cl_idx][r];
                        sigmer_replicate[cl_idx] = std::accumulate(sigmer_cluster.begin(), sigmer_cluster.end(), 0);
                        num_sigmers += sigmer_replicate[cl_idx];
                    }
                    sigmer_condition[r] = sigmer_replicate;
                }
                sigmer_data[cond_idx] = sigmer_condition;
            }

            // calculate size factors per replicate
            size_replicates.resize(conditions.size());
            double denom = pow(num_sigmers, 1 / (double)num_replicates);
            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<double> size_replicates_cond(replicates[cond_idx]);
                for (int r = 0; r < replicates[cond_idx]; r++) {
                    // find median of sigmer_replicate (divide by cross condition data first)
                    vector<double> rep_temp(sigmer_data[cond_idx][r].begin(), sigmer_data[cond_idx][r].end());
                    for (auto i = rep_temp.begin(); i != rep_temp.end(); ++i) {
                        (*i) /= denom;
                    }

                    size_replicates_cond[r] = median(rep_temp);
                }
                size_replicates[cond_idx] = size_replicates_cond;
            }

            // calculate size factors per condition
            size_conditions.resize(conditions.size());
            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                double numer = 0;
                for (int r = 0; r < replicates[cond_idx]; r++) {
                    vector<int> sigmers_replicate = sigmer_data[cond_idx][r];
                    double Nr = std::accumulate(sigmers_replicate.begin(), sigmers_replicate.end(), 0);
                    numer += (Nr / size_replicates[cond_idx][r]);
                }

                size_conditions[cond_idx] = numer / (double)replicates[cond_idx];
            }
        }

        void gibbsG(int cond_idx, int r, int cl_idx) {
            // initializations
            vector<int> G = G_data[cond_idx][cl_idx][r]; // #sigmers: transcript id (from transcripts)
            vector<int> theta = theta_data[cond_idx][cl_idx][r]; // #transcripts: occurrences of (general) sigmers per transcript
            vector<vector<int>> F = F_data[cl_idx]; // #sigmers x #transcripts
            vector<int> L = L_data[cl_idx]; // #transcripts

            int sigmer_count = (int)G.size();
            for (int s = 0; s < sigmer_count; s++) {
                vector<double> G_probs;
                // choose random t from distribution over all t: theta[r][t] * F[s][t] / L[t]
                for (int t = 0; t < (int)theta.size(); t++) {
                    if (L[t] == 0)
                        G_probs.push_back(0);
                    else
                        G_probs.push_back(theta[t] * F[s][t] / L[t]);
                }
                boost::random::discrete_distribution<> distribution(G_probs.begin(), G_probs.end());
                G[s] = distribution(generator);
            }

            G_data[cond_idx][cl_idx][r] = G;
        }

        void gibbsTheta(int cond_idx, int r, int cl_idx) {

        }

        void gibbsM(int cond_idx, int cl_idx) {

        }

        void output() {

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

        double median(vector<double> &v)
        {
            size_t n = v.size() / 2;
            nth_element(v.begin(), v.begin()+n, v.end());

            if (v.size()%2 == 1) {
                return v[n];
            } else {
                return 0.5 * (v[n] + v[n-1]);
            }
        }

        void run(){
            cout << "done" << endl;
        }
    private:
        map<string, int> conditions; // condition name -> index
        vector<int> replicates; // #conditions: number of replicates
        vector<vector<string>> cf_files; // #conditions x #replicates: cf file names
        vector<vector<string>> em_files; // #conditions x #replicates: em file names

        map<string, vector<string>> cluster_transcripts_map; // cluster -> list of associated transcripts
        map<string, string> transcript_cluster_map; // transcript -> corresponding cluster
        map<string, int> clusters; // cluster id -> index
        map<string, int> transcripts; // transcript id -> index (per cluster)

        vector<vector<vector<vector<int>>>> theta_data; // #conditions x #clusters x #replicates x #transcripts: occurrences of (general) sigmers per transcript
        vector<vector<vector<int>>> F_data; // #clusters x #sigmers x #transcripts: num (specific) sigmer per transcript
        vector<vector<int>> L_data; // #clusters x #transcripts: num (general) sigmers per transcript
        vector<vector<vector<vector<int>>>> G_data; // #conditions x #clusters x #replicates x #sigmers: transcript id (from transcripts)
        vector<vector<vector<double>>> m_data; // #conditions x #clusters x #transcripts: probability that general sigmer belongs to a transcript

        vector<vector<double>> size_replicates; // #conditions x #replicates
        vector<double> size_conditions; // #conditions

        std::default_random_engine generator;
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_input_file);
    //rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}
