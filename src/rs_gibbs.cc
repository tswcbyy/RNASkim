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
#include <errno.h>
#include <chrono>
#include <thread>
#include <mutex>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "rs_common.h"
#include "rs_estimate_lib.h"
#include "proto/rnasigs.pb.h"
#include "proto_data.h"
#include <boost/random/discrete_distribution.hpp>
#include <boost/lexical_cast.hpp>


/*extern "C" {
    #include "locfit/local.h"
}*/

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::fstream;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::stringstream;
using std::find;
using std::thread;

DEFINE_string(input_file, "", "Text file containing all input files (tab-delimited) with 3 columns: condition name, .cf file, _em file");

DEFINE_double(sampling_factor, 0, "(testing) Every nth sigmer to be included.");

#define DEBUG 1
#define NUM_THREADS 32
#define THRESHOLD 500

namespace rs {
    class GibbsSampler {
    public:
        GibbsSampler(const string input_file, const double sampling_factor) {
            cout << "sampling factor: " << sampling_factor << endl;
            // conditions, replicates, cf_files, em_files
            parseInputFile(input_file);
            if (DEBUG) {
                printf("%s: Checking input file parsing ...", currentDateTime().c_str());
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
                    assert(r > 1);
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
            parseCFData(cf_files[0][0], sampling_factor);
            if (DEBUG) {
                printf("CF data: complete\n");
            }

            // size_replicates, size_conditions
            calcSizeFactors();
            if (DEBUG) {
                printf("Checking size factors ...\n");
                for (auto cond_it : conditions) {
                    string condition = cond_it.first;
                    int cond_idx = cond_it.second;
                    printf("\t%s:\t[%f - %f]\n", condition.c_str(), *std::min_element(size_conditions[cond_idx].begin(), size_conditions[cond_idx].end()), *std::max_element(size_conditions[cond_idx].begin(), size_conditions[cond_idx].end()));
                    vector<double> sz = size_replicates[cond_idx];
                    for (double sz_factor : sz) {
                        printf("\t\t%f\n", sz_factor);
                    }
                }
            }

            // find max SrSc and calc +5 for ln_factorials
            size_t max_size = 0;
            for (size_t c = 0; c < size_conditions.size(); c++) {
                for (size_t r = 0; r < size_replicates[c].size(); r++) {
                    double max_Sc = *std::max_element(size_conditions[c].begin(), size_conditions[c].end());
                    max_size = std::max((size_t)round(max_Sc * size_replicates[c][r]), max_size);
                }
            }
            if (DEBUG) printf("Maximum SrSc = %zu\n", max_size);
            ln_factorials.resize(1, 0); // initialize ln(0!)
            calcLogFactorials(max_size + 5);
            if (DEBUG) printf("Initialized ln(n!) up to ln(%zu!) = %f\n", ln_factorials.size() - 1, ln_factorials[ln_factorials.size() - 1]);


            // m_data: first divide the theta vector for each replicate by its size factor. Then take average of all the replicates for one condition. Finally normalize it.
            initMData();
            filterClusters();

            cout << currentDateTime() << "\tClusters to be processed: " << unprocessed_clusters.size() << endl;

            for (auto cond_it : conditions) {
                output(ios::out | ios::trunc, 0, cond_it.first);
            }

            if (DEBUG) printf("initMData: complete\n");

            calcMVariance();
            if (DEBUG) printf("calcMVariance: complete\n");

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
            ofstream perf_file;
            perf_file.open ("rsgibbs_perf.dat", ios::out | ios::trunc);
            perf_file << "cluster\tsigmers\ttranscripts\tms\titeration" << endl;


            for (auto cond_it : conditions) {
                int cond_idx = cond_it.second;
                string condition = cond_it.first;
                if (DEBUG) printf("\tCondition [%i]: %s\n", cond_idx, condition.c_str());

                for (int i = 1; i < num_iterations; i++) {
                    if (DEBUG) printf("\t\t%s iteration %i: \n", currentDateTime().c_str(), i);

                    cluster_performance.clear();
                    cluster_performance.resize(clusters.size(), 0);

                    thread_unprocessed_clusters = unprocessed_clusters;
                    threads.resize(NUM_THREADS);
                    for (int thr = 0; thr < NUM_THREADS; thr++) {
                        threads.at(thr) = thread(&GibbsSampler::processCluster, this, cond_idx, thr);
                    }

                    for (int thr = 0; thr < NUM_THREADS; thr++) {
                        threads.at(thr).join();
                    }

                    output(ios::out | ios::app, i, condition);

                    for (int cl_idx = 0; cl_idx < (int)cluster_performance.size(); cl_idx++) {
                        if (cluster_performance[cl_idx] > 0) {
                            string cluster = cluster_vector[cl_idx];
                            perf_file << cluster << "\t" << F_data[cl_idx].size() << "\t" << cluster_transcripts_map[cluster].size() << "\t" << cluster_performance[cl_idx] << "\t" << i << endl;
                        }
                    }
                }

                if (DEBUG) printf("\tcomplete\n");
            }
            perf_file.close();

            finalOutput();
        }

        // order by #transcripts ascending to use threads efficiently
        // skip cluster processing with m[i] = {0,1}
        void filterClusters() {
            // TODO: fix me fix me temp stuff for figuring out timing issues
            unprocessed_clusters.clear();
            /*unprocessed_clusters.push_back(clusters["ENSMUSG00000000001"]); // 333618 sig, 8339 transcripts*/
            unprocessed_clusters.push_back(clusters["ENSMUSG00000074519"]); // 2134 sig, 121 transcripts
            unprocessed_clusters.push_back(clusters["ENSMUSG00000000182"]); // 3775	sig, 210 transcripts
            unprocessed_clusters.push_back(clusters["ENSMUSG00000030849"]);
            unprocessed_clusters.push_back(clusters["ENSMUSG00000020990"]);
            unprocessed_clusters.push_back(clusters["ENSMUSG00000051747"]);
            unprocessed_clusters.push_back(clusters["ENSMUSG00000028649"]);

            /****** end debuggy stuff *****/
            for (auto cl_it = unprocessed_clusters.begin(); cl_it != unprocessed_clusters.end(); ++cl_it) {
                string cluster = cluster_vector[*cl_it];

                bool to_process = false;
                for (auto cond_it : conditions) {
                    vector<double> m = m_data[cond_it.second][*cl_it];
                    for (double val : m) {
                        if (val != 0 && val != 1) {
                            to_process = true;
                            goto erase;
                        }
                    }
                }

                erase:
                if (!to_process) {
                    unprocessed_clusters.erase(cl_it);
                }
            }

            sortstruct s(this);
            std::sort(unprocessed_clusters.begin(), unprocessed_clusters.end(),  s);
        }

        bool clusterComparison(int cl1, int cl2) {
            string cluster1 = cluster_vector[cl1];
            string cluster2 = cluster_vector[cl2];
            size_t tr1 = cluster_transcripts_map[cluster1].size();
            size_t tr2 = cluster_transcripts_map[cluster2].size();

            return (tr1 <= tr2);
        }

        void finalOutput() {
            double m_mean, m_5 = 0, m_95 = 0, m_var = 0;

            for (auto cond_it : conditions) {
                string condition = cond_it.first;
                int cond_idx = cond_it.second;

                printf("Final output for condition %s\n", condition.c_str());

                ofstream mfile;
                mfile.open ("rsgibbs_" + condition + ".dat", ios::out | ios::trunc);
                mfile << "Transcript" << "\t" << "Mean x Sc" << "\t" << "Variance x Sc^2" << "\t" << "5%ile" << "\t" << "95%ile" << endl;

                printf("%zu iterations in m_Gibbs\n", m_Gibbs.size());
                for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                    double Sc = size_conditions[cond_idx][cl_idx];
                    string cluster = cluster_vector[cl_idx];
                    for (int t_idx = 0; t_idx < (int)cluster_transcripts_map[cluster].size(); t_idx++) {
                        string transcript = cluster_transcripts_map[cluster][t_idx];

                        vector<double> m;
                        for (int i = 0; i < (int)m_Gibbs.size(); i++) {
                            m.push_back(m_Gibbs[i][cond_idx][cl_idx][t_idx]);
                        }

                        m_mean = std::accumulate(m.begin(), m.end(), 0.0) / (double) m.size();
                        vector<double> diff(m.size());
                        std::transform(m.begin(), m.end(), diff.begin(), std::bind2nd(std::minus<double>(), m_mean));
                        m_var = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
                        m_var /= (double) m.size();

                        size_t fifth_percent = 0.05 * m.size();
                        std::nth_element(m.begin(), m.begin() + fifth_percent, m.end(), std::less<double>());
                        m_5 = m[fifth_percent];
                        std::nth_element(m.begin(), m.begin() + fifth_percent, m.end(), std::greater<double>());
                        m_95 = m[fifth_percent];

                        mfile << transcript << "\t" << m_mean * Sc << "\t" << m_var * Sc * Sc << "\t" << m_5 << "\t" << m_95 << endl;
                    }
                }

                mfile.close();
            }
        }

        void processCluster(int cond_idx, int thread_idx) {
            while(true) {
                cluster_mutex.lock();
                if (thread_unprocessed_clusters.size() == 0) {
                    cluster_mutex.unlock();
                    printf("[%i] exiting @%s\n", thread_idx, currentDateTime().c_str());
                    return;
                }

                int cl_idx = thread_unprocessed_clusters.back();
                string cluster = cluster_vector[cl_idx];
                thread_unprocessed_clusters.pop_back();
                cluster_mutex.unlock();

                if (DEBUG) printf("[%i]\tProcessing cluster [%i: %s]...\n", thread_idx, cl_idx, cluster.c_str());
                auto begin = std::chrono::high_resolution_clock::now();
                auto end = std::chrono::high_resolution_clock::now();
                auto dur = end - begin;

                // replicate threading
                vector<thread> rep_threads(replicates[cond_idx]);
                for (int r = 0; r < replicates[cond_idx]; r++) {
                    rep_threads.at(r) = thread(&GibbsSampler::processGandTheta, this, cond_idx, cl_idx, r, thread_idx);
                }
                for (int thr = 0; thr < (int)rep_threads.size(); thr++) {
                    rep_threads.at(thr).join();
                }

                gibbsM(cond_idx, cl_idx);

                end = std::chrono::high_resolution_clock::now();
                dur = end - begin;

                perf_mutex.lock();
                cluster_performance[cl_idx] = (double) std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
                perf_mutex.unlock();
            }
        }

        void processGandTheta(int cond_idx, int cl_idx, int r, int thread_idx) {
            if (DEBUG) printf("\t[%i:%i]\tProcessGandTheta...\n", thread_idx, r);
            auto begin = std::chrono::high_resolution_clock::now();
            auto end = std::chrono::high_resolution_clock::now();

            auto dur = end - begin;
            gibbsG(cond_idx, r, cl_idx);
            gibbsTheta(cond_idx, r, cl_idx);

            end = std::chrono::high_resolution_clock::now();
            dur = end - begin;
            if (DEBUG) printf("\t[%i:%i]\tCompleted processGandTheta [%f ms]\n", thread_idx, r, (double) std::chrono::duration_cast<std::chrono::milliseconds>(dur).count());
            return;
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
                        cluster_vector.push_back(cluster);
                        assert(std::strcmp(cluster_vector[clusters[cluster]].c_str(), cluster.c_str()) == 0);
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

        void parseCFData(const string cf, const double sampling_factor_) {
            /* F_data; #clusters x #sigmers x #transcripts: num (specific) sigmer per transcript
             * L_data; #clusters x #transcripts: num (general) sigmers per transcript
             * G_data; #conditions x #clusters x #replicates x #sigmers: transcript id (from transcripts) */

            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution(0.0,1.0);

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

                double sampling_factor = 0;
                if (sampling_factor_ != 0 && sk.keys_size() >= THRESHOLD) {
                    sampling_factor = sampling_factor_ * (double)num_transcripts / (double)sk.keys_size();
                    if (sampling_factor >= 1)
                        sampling_factor = 0;
                }

                // iterate over sigmers
                sigmers:
                for (int s = 0; s < sk.keys_size(); s++) {
                    // create vector -- size = num_transcripts per cluster, each entry is number of times a sigmer appears in the transcript
                    vector<int> F_sigmer_data(num_transcripts, 0);

                    if (sampling_factor != 0) {
                        double sample = distribution(generator);
                        if (sample > sampling_factor)
                            continue;
                    }

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

                if (sampling_factor != 0.) {
                    printf("Sampled %zu/%i sigmers for cluster %s with %zu transcripts [factor = %e]\n", F_cluster_data.size(), sk.keys_size(), cluster.c_str(), num_transcripts, sampling_factor);
                    if (F_cluster_data.size() == 0) {
                        sampling_factor *= 1.1;
                        printf("Re-sampling with sampling factor = %f\n", sampling_factor);
                        goto sigmers;
                    }
                }

                F_data[cl_idx] = F_cluster_data;
                L_data[cl_idx] = L_cluster_data;

                // populate initial G_data
                for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        vector<int> theta_replicate_data = theta_data[cond_idx][cl_idx][r];
                        auto max_idx = std::max_element(theta_replicate_data.begin(), theta_replicate_data.end());
                        int g = std::distance(theta_replicate_data.begin(), max_idx);
                        vector<int> G_sigmer_data(F_cluster_data.size(), g);
                        G_data[cond_idx][cl_idx][r] = G_sigmer_data;
                    }
                }
            }
        }

        void calcSizeFactors() {
            if (DEBUG) printf("Calculating size factors...\n");
            int num_replicates = std::accumulate(replicates.begin(), replicates.end(), 0);
            // size_replicate: median (over C clusters) for (#sigmers in cluster c for replicate r / product #sigmers across all replicates in cluster c ^ 1/#replicates

            // calculate #sigmers per cluster per replicate
            vector<vector<vector<int>>> sigmer_data(conditions.size()); // #conditions x #replicates x #clusters: #sigmers (from theta_data)
            vector<double> denom_data(clusters.size(), 0); // #clusters: (product_(all replicates) #sigmers) ^ 1/#replicates
            vector<int> zero_data(clusters.size(), 0); // #cluster ID : #0s

            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<vector<int>> sigmer_condition(replicates[cond_idx]);

                for (int r = 0; r < replicates[cond_idx]; r++) {
                    // #clusters: #sigmers that occur per cluster in this replicate
                    vector<int> sigmer_replicate(clusters.size());

                    for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                        vector<int> sigmer_cluster = theta_data[cond_idx][cl_idx][r];
                        sigmer_replicate[cl_idx] = std::accumulate(sigmer_cluster.begin(), sigmer_cluster.end(), 0);
                        if (sigmer_replicate[cl_idx] == 0) {
                            zero_data[cl_idx]++;
                        } else {
                            denom_data[cl_idx] += log(sigmer_replicate[cl_idx]);
                        }
                    }
                    sigmer_condition[r] = sigmer_replicate;
                }
                sigmer_data[cond_idx] = sigmer_condition;
            }

            // populate zero_cluster_map and denominator data
            for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                if (zero_data[cl_idx] == num_replicates)
                    zero_cluster_map[cl_idx] = 1;

                if (zero_data[cl_idx] != num_replicates)
                    unprocessed_clusters.push_back(cl_idx);

                errno = 0;
                if (zero_data[cl_idx] == 0)
                    denom_data[cl_idx] = exp(denom_data[cl_idx]/(double)num_replicates);
                else
                    denom_data[cl_idx] = 0;

                if (errno == ERANGE) {
                    printf("ERROR: operation overflowed: exp(%f / %i)\n", denom_data[cl_idx], num_replicates);
                    assert(false);
                }
            }
            if (DEBUG) printf("%zu/%zu clusters to be processed\n", unprocessed_clusters.size(), clusters.size());

            // calculate size factors per replicate
            size_replicates.resize(conditions.size());

            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<double> size_replicates_cond(replicates[cond_idx]);
                for (int r = 0; r < replicates[cond_idx]; r++) {

                    int zeros = 0;
                    vector<double> rep_temp;
                    for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                        if (zero_cluster_map.find(cl_idx) == zero_cluster_map.end()) {
                            if (denom_data[cl_idx] == 0 || sigmer_data[cond_idx][r][cl_idx] == 0) {
                                rep_temp.push_back(0);
                                zeros++;
                            } else {
                                rep_temp.push_back((double)sigmer_data[cond_idx][r][cl_idx] / denom_data[cl_idx]);
                            }
                        }
                    }

                    if (DEBUG) printf("\tPostfilter: %i/%zu clusters with 0 sigmers from theta data [condition:%i, replicate:%i]\n", zeros, rep_temp.size(), cond_idx, r);

                    double med = median(rep_temp);
                    size_replicates_cond[r] = med;
                }
                size_replicates[cond_idx] = size_replicates_cond;
            }

            // calculate size factors per condition & cluster
            size_conditions.resize(conditions.size());
            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<double> size_condition_cluster(clusters.size());

                for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                    double val = 0;
                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        val += sigmer_data[cond_idx][r][cl_idx] / size_replicates[cond_idx][r];
                    }
                    val /= (double) replicates[cond_idx];
                    size_condition_cluster[cl_idx] = val;
                }

                size_conditions[cond_idx] = size_condition_cluster;
            }
        }

        void calcMVariance() {
            // use LOCFIT to find wc(m) based on (mw, w)
            // v = m + wc(m) - z
            vector<vector<double>> w(conditions.size()); // #conditions x #transcripts
            vector<vector<double>> mw(conditions.size()); // #conditions x #transcripts
            vector<vector<vector<double>>> z(conditions.size()); // #conditions x #clusters x #transcripts: z bias correction term

            v_data.resize(conditions.size()); // #conditions x #clusters x #transcripts
            v_sampling_data.resize(conditions.size()); // #conditions x #clusters x Sc

            for (auto cond_it : conditions) {
                map<double, double> m_var; // LOCFIT m -> w
                int cond_idx = cond_it.second;
                string condition = cond_it.first;
                vector<double> w_condition;
                vector<double> mw_condition;
                vector<vector<double>> z_condition(clusters.size());
                vector<vector<double>> v_condition(clusters.size());
                vector<vector<double>> v_sampling_condition(clusters.size());

                for (auto cl_it : cluster_transcripts_map) {
                    size_t num_transcripts = cl_it.second.size();
                    string cluster = cl_it.first;
                    int cl_idx = clusters[cluster];

                    vector<double> z_cluster(num_transcripts, 0);
                    z_condition[cl_idx] = z_cluster;

                    vector<double> v_cluster(num_transcripts, 0);
                    v_condition[cl_idx] = v_cluster;

                    vector<double> v_sampling_cluster(round(size_conditions[cond_idx][cl_idx]) + 1, 0);
                    v_sampling_condition[cl_idx] = v_sampling_cluster;
                }

                for (auto t_it : transcripts) {
                    string transcript = t_it.first;
                    int t_idx = t_it.second;
                    string cluster = transcript_cluster_map[transcript];
                    int cl_idx = clusters[cluster];

                    if (size_conditions[cond_idx][cl_idx] == 0) continue;

                    double wt = 0;
                    double zt = 0;
                    double m = m_data[cond_idx][cl_idx][t_idx];
                    if (m == 0) {
                        w_condition.push_back(0);
                        mw_condition.push_back(0);
                        z_condition[cl_idx][t_idx] = 0;
                        continue;
                    }

                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        double theta = theta_data[cond_idx][cl_idx][r][t_idx];
                        double sz = size_conditions[cond_idx][cl_idx] * size_replicates[cond_idx][r];
                        wt += pow(theta/sz - m, 2);
                        zt += 1/sz;
                    }

                    wt *= 1 / ((double) replicates[cond_idx] - 1.);
                    zt *= m / (double) replicates[cond_idx];

                    w_condition.push_back(wt);
                    mw_condition.push_back(m);
                    z_condition[cl_idx][t_idx] = zt;
                }
                w[cond_idx] = w_condition;
                mw[cond_idx] = mw_condition;
                z[cond_idx] = z_condition;

                // create files for R consumption

                // write file for locfitting
                ofstream ofile;
                ofile.open (R_file, ios::out | ios::trunc);
                ofile << "m\tw" << endl;
                ofile.precision(20);

                for (size_t i = 0; i < mw_condition.size(); i++) {
                    ofile << mw_condition[i] << "\t" << w_condition[i] << endl;
                }

                ofile.close();

                // write file for prediction only (used in sampling m)
                ofile.open(R_predict_file, ios::out | ios::trunc);
                ofile << "m" << endl;

                for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                    if (size_conditions[cond_idx][cl_idx] == 0) continue;
                    for (int i = 0; i < (int)v_sampling_condition[cl_idx].size(); i++) {
                        double m_candidate = (double) i / size_conditions[cond_idx][cl_idx];
                        ofile << m_candidate << endl;
                    }
                }

                ofile.close();

                system("Rscript gibbs.R");

                string line;
                ifstream ifile;
                ifile.open(R_file, ios::in);
                if (ifile.is_open()) {
                    size_t l = 0;
                    while (getline(ifile, line)) {
                        vector<string> tokens = split(line, '\t');
                        // m, w, poisson
                        double w = boost::lexical_cast<double>(tokens[2].c_str()); // try poisson
                        double m = mw_condition[l];
                        m_var[m] = w;
                        if (m > 0) {
                            if (w <= 0)
                                printf("Line %zu: %s\n", l, line.c_str());
                            assert(w > 0);
                        }
                        l++;
                    }
                    assert(l == mw_condition.size());
                    ifile.close();
                } else {
                    printf("Failed to open %s\n", R_file.c_str());
                }

                if (DEBUG) { // keep copies of the locfit output files
                    ifstream src(R_file);
                    ofstream dst("rsgibbs_lfout_" + condition + ".dat", ios::trunc);

                    dst << src.rdbuf();
                }

                ifile.open(R_predict_file, ios::in);
                if (ifile.is_open()) {
                    size_t cl_idx = 0, sc = 0;
                    while (getline(ifile, line)) {
                        vector<string> tokens = split(line, '\t');
                        // m, fitted
                        double w = boost::lexical_cast<double>(tokens[1].c_str());
                        if (size_conditions[cond_idx][cl_idx] == 0 || sc >= v_sampling_condition[cl_idx].size()) {
                            cl_idx++;
                            sc = 0;
                        }
                        v_sampling_condition[cl_idx][sc] = w;
                        if (sc > 0) assert(w > 0);
                        sc++;
                    }
                    printf("Read up through cluster %zu, total = %zu\n", cl_idx, clusters.size());
                    ifile.close();
                } else {
                    printf("Failed to open %s\n", R_predict_file.c_str());
                }

                if (DEBUG) { // keep copies of the locfit output files
                    ifstream src(R_predict_file);
                    ofstream dst("rsgibbs_predict_" + condition + ".dat", ios::trunc);

                    dst << src.rdbuf();
                }

                // v = m + wc(m) - z
                for (auto t_it : transcripts) {
                    string transcript = t_it.first;
                    int t_idx = t_it.second;
                    string cluster = transcript_cluster_map[transcript];
                    int cl_idx = clusters[cluster];

                    double m = m_data[cond_idx][cl_idx][t_idx];
                    if (m != 0) {
                        assert(m_var.find(m) != m_var.end());
                        v_condition[cl_idx][t_idx] = m + m_var[m] - z_condition[cl_idx][t_idx];
                        if (v_condition[cl_idx][t_idx] <= 0) {
                            printf("WARNING: variance <= 0: m + w(%e) - z = %e + %e - %e = %e\n", m, m, m_var[m], z_condition[cl_idx][t_idx], v_condition[cl_idx][t_idx]);
                        }
                    } else
                        v_condition[cl_idx][t_idx] = 0;
                }

                v_data[cond_idx] = v_condition;
                v_sampling_data[cond_idx] = v_sampling_condition;
            }
        }

        void initMData() {
            /*
                1) divide each theta vector by size_replicate
                2) average over all replicates (per condition)
                3) normalize
             */

            // m_data; #conditions x #clusters x #transcripts: probability that general sigmer belongs to a transcript

            m_data.resize(conditions.size()); // #conditions x #clusters x #transcripts
            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<vector<double>> m_data_condition(clusters.size()); // #clusters x #transcripts

                for (auto cl_it : clusters) {
                    string cluster = cl_it.first;
                    int cl_idx = cl_it.second;

                    int num_transcripts = (int)cluster_transcripts_map[cluster].size();
                    vector<double> m_data_cluster(num_transcripts, 0);

                    // divide each theta vector by size_replicate
                    vector<vector<int>> temp = theta_data[cond_idx][cl_idx];
                    vector<vector<double>> theta(replicates[cond_idx]); // #replicates x #transcripts
                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        vector<double> theta_replicate(num_transcripts, 0);

                        double sz = size_replicates[cond_idx][r];
                        if (sz != 0 && zero_cluster_map.find(cl_idx) == zero_cluster_map.end()) {
                            for (int t = 0; t < num_transcripts; t++) {
                                theta_replicate[t] = (double) temp[r][t] / sz;
                            }
                        }

                        theta[r] = theta_replicate;
                    }

                    // avg theta over all replicates
                    for (int t = 0; t < num_transcripts; t++) {
                        double sum = 0;
                        for (int r = 0; r < replicates[cond_idx]; r++) {
                            sum += theta[r][t];
                        }
                        m_data_cluster[t] = sum / replicates[cond_idx];
                    }

                    // normalize
                    normalize(m_data_cluster);
                    m_data_condition[cl_idx] = m_data_cluster;
                }

                m_data[cond_idx] = m_data_condition;
            }

            //m_Gibbs; num_iterations x #conditions x #clusters x #transcripts: m data over all iterations
            m_Gibbs.resize(num_iterations);
            for (int i = 0; i < num_iterations; i++) {
                vector<vector<vector<double>>> m_cond(conditions.size());
                for (auto cond_it : conditions) {
                    int cond_idx = cond_it.second;
                    vector<vector<double>> m_clusters(clusters.size());
                    for (auto cl_it : clusters) {
                        int cl_idx = cl_it.second;
                        vector<double> m_transcripts;
                        m_clusters[cl_idx] = m_transcripts;
                    }
                    m_cond[cond_idx] = m_clusters;
                }
                m_Gibbs[i] = m_cond;
            }
        }

        // returns ln (Pr(theta | m))
        double negBinomialDist(int theta, double p, double q) {
            errno = 0;
            size_t a = theta + std::max(round(q), 1.) - 1;
            double r = lnFactorial(a) - lnFactorial((size_t)theta) - lnFactorial(a - theta);
            double l = r + q * log(1. - p) + (double)theta * log(p);

            if (errno != 0) {
                printf("ERROR: error occurred in negBinomialDist: %s (return value: ln C(a, th) + q ln(1-p) + th ln(p) = ln C(%zu, %i) + %e ln(1-%e) + %i ln(%e) = %e)\n", strerror(errno), a, theta, q, p, theta, p, l);
                assert(false);
            }

            return l;
        }

        void gibbsG(int cond_idx, int r, int cl_idx) {
            // 1) For each cluster c, for each replicate r, for each sigmer occurrence s: re-sample G (transcript)
            //      P(G = t | ...) = theta_t * F_s,t / L_t

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

            G_mutex.lock();
            G_data[cond_idx][cl_idx][r] = G;
            G_mutex.unlock();
        }

        void gibbsTheta(int cond_idx, int r, int cl_idx) {
            // For each replicate cluster c, for each replicate r:
            //      P(theta_rt | ...) = P(theta_r,t | m_t) * (theta_r,t / (sum_(i != t) theta_r,i + theta_r,t))^(#G == t)

            vector<int> theta = theta_data[cond_idx][cl_idx][r]; // #transcripts: occurrences of (general) sigmers per transcript
            vector<int> G = G_data[cond_idx][cl_idx][r]; // #sigmers: transcript id (from transcripts)
            if (size_replicates[cond_idx][r] == 0)
                return;

            int num_transcripts = (int)theta.size();
            int sz = round(size_replicates[cond_idx][r] * size_conditions[cond_idx][cl_idx]);

            double th_sum = 0;
            for (int t = 0; t < num_transcripts; t++) {
                th_sum += theta[t];
            }

            for (int t = 0; t < num_transcripts; t++) {
                vector<double> theta_probs;
                int g_transcript = 0; // #sigmers that match t
                for (int g_idx : G) {
                    if (g_idx == t)
                        g_transcript++;
                }

                double m = m_data[cond_idx][cl_idx][t];
                double v = v_data[cond_idx][cl_idx][t];

                double denom, factor, prob_theta;
                int erange_errors = 0, total_probs = 0;
                if (v != 0 && m != 0) {
                    double p = 1. - m / ((double)sz * v);
                    double q = ((double)sz * m * m) / ((double)sz * v - m);

                    if (p <= 0 || p >= 1) printf("ERROR: p = %.20e, q = %.20e\n", p, q);

                    if (size_conditions[cond_idx][cl_idx] * size_replicates[cond_idx][r] > THRESHOLD) {
                        int min_theta = floor(std::max(0., m - 3*pow(v, 0.5)) * size_conditions[cond_idx][cl_idx]) + 1;
                        int max_theta = ceil(std::min(1., m + 3*pow(v, 0.5)) * size_conditions[cond_idx][cl_idx]);

                        theta_probs.resize(max_theta + 1, 0);

                        for (int th = min_theta; th <= max_theta; th++) {
                            int theta_candidate = round(th * size_replicates[cond_idx][r]);
                            if (theta_candidate == 0) continue;

                            prob_theta = negBinomialDist(theta_candidate, p, q); // ln()
                            denom = th_sum - (double)theta[t] + (double)theta_candidate;
                            if (denom == 0) {
                                printf("ERROR: sum(theta) = 0 for cond[%i], cluster[%i], tr[%i], theta[%i]\n", cond_idx, cl_idx, t, theta[t]);
                            } else {
                                double frac = (double)theta_candidate / denom;
                                if (errno != 0) {
                                    printf("ERROR: error occurred in gibbsTheta: %s\n", strerror(errno));
                                    assert(false);
                                }
                                errno = 0;
                                factor = g_transcript * log(frac);
                                if (errno == ERANGE) {
                                    printf("ERROR: operation overflowed: g_transcript * log(factor) = %i * log(%e) = %e\n", g_transcript, frac, factor);
                                    assert(false);
                                }

                                theta_probs[th] = exp(prob_theta + factor);
                                if (errno == ERANGE) {
                                    if (prob_theta + factor > 0) {
                                        printf("ERROR: operation overflowed (positive) exp(%f)\n", prob_theta + factor);
                                        assert(false);
                                    }
                                    erange_errors++;
                                    errno = 0;
                                }
                                if (theta_probs[th] != 0) total_probs++;
                            }
                        }
                        
                        boost::random::discrete_distribution<> distribution(theta_probs.begin(), theta_probs.end());
                        theta[t] = distribution(generator) * size_replicates[cond_idx][r]; // result in {1, ..., sz}
                    } else {
                        theta_probs.resize(sz, 0);

                        for (int i = 0; i < sz; i++) {
                            int th = i + 1;
                            prob_theta = negBinomialDist(th, p, q); // ln()
                            denom = th_sum - (double)theta[t] + (double)th;
                            double frac = (double)th / denom;
                            factor = g_transcript * log(frac);
                            theta_probs[i] = exp(prob_theta + factor);
                        }
                        boost::random::discrete_distribution<> distribution(theta_probs.begin(), theta_probs.end());
                        theta[t] = distribution(generator) + 1; // result in {1, ..., sz}

                    }

                    /*if (erange_errors > 0)
                        printf("\t\t%s\t%i erange errors : %i nonzero sampling P(theta|...) for transcript %i/%i (v = %e, m = %e) (sz = %i)\n", currentDateTime().c_str(), erange_errors, total_probs, t, num_transcripts, v, m, sz);*/
                } else {
                    //printf("WARNING: v[%e], m[%e] for cond[%i], cluster[%i], tr[%i], theta[%i]\n", v, m, cond_idx, cl_idx, t, theta[t]);
                }

            }
            theta_mutex.lock();
            theta_data[cond_idx][cl_idx][r] = theta;
            theta_mutex.unlock();
        }

        void gibbsM(int cond_idx, int cl_idx) {
            vector<double> m = m_data[cond_idx][cl_idx]; // #transcripts
            double Sc = round(size_conditions[cond_idx][cl_idx]);

            for (int t = 0; t < (int)m.size(); t++) {
                if (m[t] == 0) continue;

                // limit probability vector to +-3 std dev
                double min_m = 0, max_m = 0;
                double v = v_data[cond_idx][cl_idx][t];

                for (int r = 0; r < replicates[cond_idx]; r++) {
                    min_m += (double) theta_data[cond_idx][cl_idx][r][t] / size_replicates[cond_idx][r];
                }
                min_m /= size_conditions[cond_idx][cl_idx];
                max_m = min_m;

                min_m -= (3. * size_conditions[cond_idx][cl_idx] * pow(v, 0.5));
                max_m += (3. * size_conditions[cond_idx][cl_idx] * pow(v, 0.5));

                min_m = floor(min_m);
                max_m = ceil(max_m);

                min_m = std::max(0., min_m);
                max_m = std::min(Sc, max_m); // TODO: what happens if max_m is 0? do we set this m = 0 or 1/Sc (current)

                vector<double> m_probs(max_m);
                int erange_errors = 0, total_probs = 0;
                for (int sz = 1; sz <= (int)m_probs.size(); sz++) {
                    double product = 1;
                    double m_candidate = (double)sz / Sc;
                    double v = v_sampling_data[cond_idx][cl_idx][sz];
                    if (v <= 0)
                        printf("ERROR: cond[%i], sz = %i, m = %e, v = %e\n", cond_idx, sz, m_candidate, v);
                    assert(v > 0);
                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        double s = round(size_conditions[cond_idx][cl_idx] * size_replicates[cond_idx][r]);
                        double p = 1. - m_candidate / (s * v);
                        double q = ((double)sz * m_candidate * m_candidate) / (s * v - m_candidate);

                        int theta = theta_data[cond_idx][cl_idx][r][t];

                        if (p <= 0 || p == 1) {
                            //printf("WARNING: p = 1 - m/sv = 1 - %e / (%e * %e) = %e\n", m_candidate, s, v, p);
                            product = 1;
                            break;
                        }
                        product += negBinomialDist(theta, p, q);
                    }

                    errno = 0;
                    if (errno != 0) {
                        printf("ERROR: error occurred in gibbsM: %s\n", strerror(errno));
                        assert(false);
                    }
                    m_probs[sz - 1] = exp(product);
                    if (errno == ERANGE) {
                        if (product > 0) {
                            printf("ERROR: operation overflowed (positive) exp(%f)\n", product);
                            assert(false);
                        }
                        errno = 0;
                        erange_errors++;
                    }
                    if (m_probs[sz - 1] > 0) total_probs++;
                }
                /*if (erange_errors > 0)
                    printf("\t\t%s\t%i erange errors : %i nonzero sampling P(m|...) for transcript %i/%zu\n", currentDateTime().c_str(), erange_errors, total_probs, t, m.size());*/

                boost::random::discrete_distribution<> distribution(m_probs.begin(), m_probs.end());
                m[t] = ((double) distribution(generator) + 1.) / Sc;
                //printf("gibbsM: m %f->%f from (%f, %f) out of possible %f\n", m_data[cond_idx][cl_idx][t], m[t], min_m, max_m, Sc);
            }

            normalize(m);
            m_mutex.lock();
            m_data[cond_idx][cl_idx] = m;
            m_mutex.unlock();
        }

        void normalize(vector<double>& v) {
            double norm_factor = std::accumulate(v.begin(), v.end(), 0.0);
            if (norm_factor != 0) {
                for (auto i = v.begin(); i != v.end(); ++i) {
                    (*i) /= norm_factor;
                }
            }
        }

        // for debugging: prints m and theta vector after every iteration in out/cluster.dat files
        void output(std::ios_base::openmode mode, int i, string condition) {
            int cond_idx = conditions[condition];
            for (int cl = 0; cl < (int)unprocessed_clusters.size(); cl++) {
                int cl_idx = unprocessed_clusters[cl];
                string cluster = cluster_vector[cl_idx];

                ofstream m_file, theta_file;
                string file_name = "out/m_" + condition + "." + cluster + ".dat";
                m_file.open (file_name.c_str(), mode);

                file_name = "out/theta_" + condition + "." + cluster + ".dat";
                theta_file.open (file_name.c_str(), mode);

                // print header line
                if (i == 0) {
                    theta_file << "R\t";
                    for (int t_idx = 0; t_idx < (int)m_data[cond_idx][cl_idx].size(); t_idx++) {
                        string transcript = cluster_transcripts_map[cluster][t_idx];
                        m_file << transcript << "\t";
                        theta_file << transcript << "\t";
                    }
                    m_file << endl;
                    theta_file << endl;
                }

                for (int t_idx = 0; t_idx < (int)m_data[cond_idx][cl_idx].size(); t_idx++) {
                    m_file << m_data[cond_idx][cl_idx][t_idx] << "\t";
                }

                for (int r = 0; r < replicates[cond_idx]; r++) {
                    theta_file << r << "\t";
                    for (int t_idx = 0; t_idx < (int)m_data[cond_idx][cl_idx].size(); t_idx++) {
                        theta_file << theta_data[cond_idx][cl_idx][r][t_idx] << "\t";
                    }
                    theta_file << endl;
                }

                m_file << endl;
                m_file.close();
                theta_file.close();
            }
            m_Gibbs[i][cond_idx] = m_data[cond_idx];
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

        void calcLogFactorials(size_t n) {
            if (ln_factorials.size() > n)
                return;

            size_t c = ln_factorials.size();
            ln_factorials.resize(n + 1, 0);
            for (; c < ln_factorials.size(); c++) {
                ln_factorials[c] = ln_factorials[c-1] + log((double)c);
            }
        }

        double lnFactorial(size_t n) {
            if (n < ln_factorials.size()) {
                return ln_factorials[n];
            }

            if (DEBUG) printf("Extending ln_factorials from %zu to %zu\n", ln_factorials.size() - 1, n);
            calcLogFactorials(n);
            return ln_factorials[n];
        }

        const string currentDateTime() {
            time_t     now = time(0);
            struct tm  tstruct;
            char       buf[80];
            tstruct = *localtime(&now);
            // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
            // for more information about date/time format
            strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
            
            return buf;
        }

        void run() {
            cout << "done" << endl;
        }
    private:
        int num_iterations = 100;
        map<string, int> conditions; // condition name -> index
        vector<int> replicates; // #conditions: number of replicates
        vector<vector<string>> cf_files; // #conditions x #replicates: cf file names
        vector<vector<string>> em_files; // #conditions x #replicates: em file names

        map<string, vector<string>> cluster_transcripts_map; // cluster -> list of associated transcripts
        map<string, string> transcript_cluster_map; // transcript -> corresponding cluster
        map<string, int> clusters; // cluster id -> index
        vector<string> cluster_vector; // list of clusters in order of id
        map<string, int> transcripts; // transcript id -> index (per cluster)
        vector<int> thread_unprocessed_clusters; // vector of clusters to be processed (thread)
        vector<int> unprocessed_clusters; // vector of clusters to be processed (original)

        std::mutex cluster_mutex;
        std::mutex G_mutex;
        std::mutex m_mutex;
        std::mutex theta_mutex;
        std::mutex perf_mutex;

        vector<vector<vector<vector<int>>>> theta_data; // #conditions x #clusters x #replicates x #transcripts: occurrences of (general) sigmers per transcript
        vector<vector<vector<int>>> F_data; // #clusters x #sigmers x #transcripts: num (specific) sigmer per transcript
        vector<vector<int>> L_data; // #clusters x #transcripts: num (general) sigmers per transcript
        vector<vector<vector<vector<int>>>> G_data; // #conditions x #clusters x #replicates x #sigmers: transcript id (from transcripts)
        vector<vector<vector<double>>> m_data; // #conditions x #clusters x #transcripts: probability that general sigmer belongs to a transcript
        vector<vector<vector<double>>> v_data; // #conditions x #clusters x #transcripts: variance of m
        vector<vector<vector<double>>> v_sampling_data; // #conditions x #clusters x (#Sc + 1): variance of m = i/Sc
        vector<vector<vector<vector<double>>>> m_Gibbs; // 1000 x #conditions x #clusters x #transcripts: m data over all 1000 iterations
        vector<vector<vector<vector<vector<double>>>>> theta_P_data; // P(theta | m) #conditions x #replicates x #clusters x #transcripts x #SrSc based on m_data

        vector<vector<double>> size_replicates; // #conditions x #replicates
        vector<vector<double>> size_conditions; // #conditions x #clusters

        vector<double> ln_factorials; // ln_factorials[i] = ln(i!) = ln(1) + ln(2) + ln(3) ... + ln(i)

        map<int, int> zero_cluster_map; // cluster ID -> no sigmers across all replicates

        std::default_random_engine generator;
        string R_file = "rsgibbs_locfit.dat"; // Note: if this value is edited, ensure that gibbs.R is also adjusted accordingly
        string R_predict_file = "rsgibbs_predict.dat";

        vector<double> cluster_performance;
        vector<thread> threads;

        struct sortstruct
        {
            // sortstruct needs to know its containing object
            rs::GibbsSampler* m;
            sortstruct(rs::GibbsSampler* p) : m(p) {};

            // this is our sort function, which makes use
            // of some non-static data (sortascending)
            bool operator() (int cl1, int cl2) {
                string cluster1 = m->cluster_vector[cl1];
                string cluster2 = m->cluster_vector[cl2];
                size_t tr1 = m->cluster_transcripts_map[cluster1].size();
                size_t tr2 = m->cluster_transcripts_map[cluster2].size();

                return (tr1 <= tr2);
            }
        };
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_input_file, FLAGS_sampling_factor);
    //rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}
