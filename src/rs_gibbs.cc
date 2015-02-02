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

extern "C" {
    #include "locfit/local.h"
}

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
/*
                // TODO:
                vector<double> m_mean(num_transcripts_cluster); // mean m across 1000 iterations x #transcripts
                vector<double> m_25(num_transcripts_cluster); // 25th m x #transcripts
                vector<double> m_975(num_transcripts_cluster); // 975th m x #transcripts
                vector<double> m_var(num_transcripts_cluster); // variance of m x #transcripts
                
                for (int i = 0; i < num_transcripts_cluster; i++) {
                    fprintf(output, "%s\t%f\t%f\t%f\t%f\n", transcripts[i].c_str(), m_mean[i] * size_condition, m_var[i] * size_condition * size_condition, m_25[i], m_975[i]);
                }
            }
            
        }*/

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
            parseCFData(cf_files[0][0]);
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
                    printf("\t%s:\t%f\n", condition.c_str(), size_conditions[cond_idx]);
                    vector<double> sz = size_replicates[cond_idx];
                    for (double sz_factor : sz) {
                        printf("\t\t%f\n", sz_factor);
                    }
                }
            }

            // m_data: first divide the theta vector for each replicate by its size factor. Then take average of all the replicates for one condition. Finally normalize it.
            initMData();
            if (DEBUG) printf("initMData: complete\n");


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
                if (DEBUG) printf("\tCondition [%i]: %s\n", cond_idx, condition.c_str());

                for (int i = 0; i < num_iterations; i++) {
                    if (DEBUG) printf("\t\tIteration %i: \n", i);
                    calcMVariance();
                    if (DEBUG) printf("\t\t\tcalcMVariance: complete\n");
                    for (auto cl_it : clusters) {
                        int cl_idx = cl_it.second;
                        if (zero_cluster_map.find(cl_idx) != zero_cluster_map.end()) {
                            printf("\t\t\tSKIPPING cluster %s\n", cl_it.first.c_str());
                            continue;
                        }

                        for (int r = 0; r < replicates[cond_idx]; r++) {
                            gibbsG(cond_idx, r, cl_idx);
                            gibbsTheta(cond_idx, r, cl_idx);
                        }
                        if (DEBUG) printf("\t\t\t\tCluster %i: %s sampling m\n", cl_idx, cl_it.first.c_str());
                        gibbsM(cond_idx, cl_idx, i);
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
            if (DEBUG) printf("Calculating size factors...\n");
            int num_replicates = std::accumulate(replicates.begin(), replicates.end(), 0);
            int num_sigmers = 0;
            // size_replicate: median (over C clusters) for (#sigmers in cluster c for replicate r / sum #sigmers across all replicates in cluster c ^ 1/#replicates

            // calculate #sigmers per cluster per replicate
            vector<vector<vector<int>>> sigmer_data(conditions.size()); // #conditions x #replicates x #clusters: #sigmers (from theta_data)

            // #cluster ID : #0s
            vector<int> zero_data(clusters.size(), 0);

            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<vector<int>> sigmer_condition(replicates[cond_idx]);

                for (int r = 0; r < replicates[cond_idx]; r++) {
                    // #clusters: #sigmers that occur per cluster in this replicate
                    vector<int> sigmer_replicate(clusters.size());

                    int zeros = 0;
                    for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                        vector<int> sigmer_cluster = theta_data[cond_idx][cl_idx][r];
                        sigmer_replicate[cl_idx] = std::accumulate(sigmer_cluster.begin(), sigmer_cluster.end(), 0);
                        num_sigmers += sigmer_replicate[cl_idx];
                        if (sigmer_replicate[cl_idx] == 0) {
                            zeros++;
                            zero_data[cl_idx]++;
                        }
                    }
                    sigmer_condition[r] = sigmer_replicate;
                }
                sigmer_data[cond_idx] = sigmer_condition;
            }

            // populate zero_cluster_map
            for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                if (zero_data[cl_idx] == num_replicates)
                    zero_cluster_map[cl_idx] = 1;
            }
            if (DEBUG) printf("%zu/%zu clusters to be skipped\n", zero_cluster_map.size(), clusters.size());

            // calculate size factors per replicate
            size_replicates.resize(conditions.size());
            double denom = pow(num_sigmers, 1 / (double)num_replicates);
            if (DEBUG) printf("\t#sigmers: %i\n\t#replicates = %i\n\t#s^(1/#r) = %f\n", num_sigmers, num_replicates, denom);

            for (int cond_idx = 0; cond_idx < (int)conditions.size(); cond_idx++) {
                vector<double> size_replicates_cond(replicates[cond_idx]);
                for (int r = 0; r < replicates[cond_idx]; r++) {

                    int zeros = 0;
                    vector<double> rep_temp;
                    for (int cl_idx = 0; cl_idx < (int)clusters.size(); cl_idx++) {
                        if (zero_cluster_map.find(cl_idx) == zero_cluster_map.end()) {
                            rep_temp.push_back((double)sigmer_data[cond_idx][r][cl_idx]);
                            if (sigmer_data[cond_idx][r][cl_idx] == 0)
                                zeros++;
                        }
                    }

                    if (DEBUG) printf("\tPostfilter: %i/%zu clusters with 0 sigmers from theta data [condition:%i, replicate:%i]\n", zeros, rep_temp.size(), cond_idx, r);

                    // find median of sigmer_replicate (divide by cross condition data first)
                    for (auto i = rep_temp.begin(); i != rep_temp.end(); ++i) {
                        (*i) /= denom;
                    }

                    double med = median(rep_temp);
                    size_replicates_cond[r] = med;
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

        void calcMVariance() {
            // use LOCFIT to find wc(m) based on (mw, w)
            // v = m + wc(m) - z
            vector<vector<double>> w(conditions.size()); // #conditions x #transcripts
            vector<vector<double>> mw(conditions.size()); // #conditions x #transcripts
            vector<vector<vector<double>>> z(conditions.size()); // #conditions x #clusters x #transcripts: z bias correction term

            v_data.resize(conditions.size()); // #conditions x #clusters x #transcripts

            for (auto cond_it : conditions) {
                map<double, double> m_var; // LOCFIT m -> w
                int cond_idx = cond_it.second;
                string condition = cond_it.first;
                vector<double> w_condition;
                vector<double> mw_condition;
                vector<vector<double>> z_condition(clusters.size());
                vector<vector<double>> v_condition(clusters.size());

                for (auto cl_it : cluster_transcripts_map) {
                    size_t num_transcripts = cl_it.second.size();
                    string cluster = cl_it.first;
                    int cl_idx = clusters[cluster];

                    vector<double> z_cluster(num_transcripts, 0);
                    z_condition[cl_idx] = z_cluster;

                    vector<double> v_cluster(num_transcripts, 0.000001); // TODO: fixme
                    v_condition[cl_idx] = v_cluster;
                }

                for (auto t_it : transcripts) {
                    string transcript = t_it.first;
                    int t_idx = t_it.second;
                    string cluster = transcript_cluster_map[transcript];
                    int cl_idx = clusters[cluster];

                    double wt = 0;
                    double zt = 0;
                    double m = m_data[cond_idx][cl_idx][t_idx];
                    if (m == 0) continue;

                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        double theta = theta_data[cond_idx][cl_idx][r][t_idx];
                        double sz = round(size_conditions[cond_idx] * size_replicates[cond_idx][r]);
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

                // TODO: LOCFIT to find wc(m) based on (mw, w)
                setuplf();
                if (DEBUG) printf("calcMVariance: Creating LOCFIT vars mw and w...\n");

                char namebuf[256];
                sprintf(namebuf, "mw");
                vari* loc_mw = createvar(namebuf,STREGULAR,mw_condition.size(),VDOUBLE);
                for (size_t i = 0; i < mw_condition.size(); ++i) {
                    loc_mw->dpr[i] = mw_condition[i];
                }
                sprintf(namebuf, "w");
                vari* loc_w = createvar(namebuf,STREGULAR,w_condition.size(),VDOUBLE);
                for (size_t i = 0; i < w_condition.size(); ++i) {
                    loc_w->dpr[i] = w_condition[i];
                }
                char locfit_cmd[2048];
                sprintf(locfit_cmd, "locfit w~mw");

                locfit_dispatch(locfit_cmd);

                sprintf(namebuf, "m");
                vari* m_domain = createvar(namebuf,STREGULAR,mw_condition.size(),VDOUBLE);
                for (size_t i = 0; i < mw_condition.size(); ++i) {
                    m_domain->dpr[i] = mw_condition[i];
                }

                sprintf(locfit_cmd, "w_variance=predict m");
                locfit_dispatch(locfit_cmd);

                int n = 0;
                sprintf(namebuf, "w_variance");
                vari* w = findvar(namebuf, 1, &n);
                assert(w != NULL);

                for (size_t i = 0; i < w->n; ++i) {
                    m_var[mw_condition[i]] = w->dpr[i];
                }

                // v = m + wc(m) - z
                for (auto t_it : transcripts) {
                    string transcript = t_it.first;
                    int t_idx = t_it.second;
                    string cluster = transcript_cluster_map[transcript];
                    int cl_idx = clusters[cluster];

                    double m = m_data[cond_idx][cl_idx][t_idx];

                    v_condition[cl_idx][t_idx] = m_var[m] - z_condition[cl_idx][t_idx];
                }
                v_data[cond_idx] = v_condition;
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

            //m_Gibbs; num_iterations x #conditions x #transcripts: m data over all iterations
            m_Gibbs.resize(num_iterations);
            for (int i = 0; i < num_iterations; i++) {
                vector<vector<double>> m_cond(conditions.size());
                for (auto cond_it : conditions) {
                    int cond_idx = cond_it.second;
                    vector<double> m_transcripts;
                    m_cond[cond_idx] = m_transcripts;
                }
                m_Gibbs[i] = m_cond;
            }
        }

        double negBinomialDist(int theta, double m, double v, int sz) {
            double p, q, s = (double) sz;

            if (s * v - m == 0) {
                printf("WARNING: SrScv - m = 0\n");
                return 0;
            }

            p = 1 - m / (s * v);
            q = (s * m * m) / (s * v - m);
            if (round(q) <= 0) {
                q = std::max(round(q), 1.);
            }

            // q = 1; TODO: FIX ME WHEN LOCFIT IS INCLUDED

            int a = theta + round(q) - 1;
            unsigned long long int c = choose(a, theta);
            if (q > 10)
                printf("WARNING: C(%i, %i) = %llu\n", a, theta, c);
            double r = ((double)c * pow(1 - p, q) * pow(p, theta));

            return r;
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

            G_data[cond_idx][cl_idx][r] = G;
        }

        void gibbsTheta(int cond_idx, int r, int cl_idx) {
            printf("gibbsTheta:  cond[%i], r[%i], cluster[%i]\n", cond_idx, r, cl_idx);
            // For each replicate cluster c, for each replicate r:
            //      P(theta_rt | ...) = P(theta_r,t | m_t) * (theta_r,t / (sum_(i != t) theta_r,i + theta_r,t))^(#G == t)

            vector<int> theta = theta_data[cond_idx][cl_idx][r]; // #transcripts: occurrences of (general) sigmers per transcript
            vector<int> G = G_data[cond_idx][cl_idx][r]; // #sigmers: transcript id (from transcripts)

            int num_transcripts = (int)theta.size();
            int sz = round(size_replicates[cond_idx][r] * size_conditions[cond_idx]);

            double th_sum = 0;
            for (int t = 0; t < num_transcripts; t++) {
                th_sum += theta[t];
            }

            for (int t = 0; t < num_transcripts; t++) {
                vector<double> theta_probs(sz, 0);
                int g_transcript = 0; // #sigmers that match t
                for (int g_idx : G) {
                    if (g_idx == t)
                        g_transcript++;
                }

                double m = m_data[cond_idx][cl_idx][t];
                double v = v_data[cond_idx][cl_idx][t];

                double denom, factor, prob_theta;

                if (v != 0 && m != 0) {
                    for (int i = 0; i < sz; i++) {
                        int th = i + 1;

                        prob_theta = negBinomialDist(th, m, v, sz);
                        denom = th_sum - (double)theta[t] + (double)th;
                        if (denom == 0) {
                            printf("WARNING: sum(theta) = 0 for cond[%i], cluster[%i], tr[%i], theta[%i]\n", cond_idx, cl_idx, t, theta[t]);
                        } else {
                            factor = (double)th / denom;
                            factor = pow(factor, g_transcript);
                            theta_probs[i] = (prob_theta * factor);
                        }
                    }
                    boost::random::discrete_distribution<> distribution(theta_probs.begin(), theta_probs.end());
                    theta[t] = distribution(generator) + 1; // result in {1, ..., sz}
                } else {
                    //printf("WARNING: v[%e], m[%e] for cond[%i], cluster[%i], tr[%i], theta[%i]\n", v, m, cond_idx, cl_idx, t, theta[t]);
                }

            }

            theta_data[cond_idx][cl_idx][r] = theta;
        }

        void gibbsM(int cond_idx, int cl_idx, int i) {
            vector<double> m = m_data[cond_idx][cl_idx];
            double Sc = round(size_conditions[cond_idx]);

            for (int t = 0; t < (int)m.size(); t++) {
                vector<double> m_probs(Sc);
                for (int sz = 1; sz <= (int)m_probs.size(); sz++) {
                    double product = 1;
                    double m_candidate = (double)sz / Sc;
                    double v = v_data[cond_idx][cl_idx][t];

                    for (int r = 0; r < replicates[cond_idx]; r++) {
                        int theta = theta_data[cond_idx][cl_idx][r][t];
                        product *= negBinomialDist(theta, m_candidate, v, round(size_conditions[cond_idx] * size_replicates[cond_idx][r]));
                    }
                    m_probs[sz - 1] = product;
                }
                boost::random::discrete_distribution<> distribution(m_probs.begin(), m_probs.end());
                m[t] = (distribution(generator) + 1) / Sc;
            }
            normalize(m);
            m_data[cond_idx][cl_idx] = m;
            m_Gibbs[i][cond_idx].insert(m_Gibbs[i][cond_idx].end(), m.begin(), m.end());
        }

        void normalize(vector<double>& v) {
            double norm_factor = std::accumulate(v.begin(), v.end(), 0.0);
            if (norm_factor != 0) {
                for (auto i = v.begin(); i != v.end(); ++i) {
                    (*i) /= norm_factor;
                }
            }
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

        unsigned long long int choose(int a, int b) {
            unsigned long long int r = 1, numer = 1, denom = 1, temp;
            assert(a >= b);
            assert(a >= 0  && b >= 0);

            b = std::min(b, a - b);

            for (int i = 0; i < b; i++) {
                if (numer%(i + 1) == 0)
                    numer /= (i + 1);
                else
                    denom *= (i + 1);

                temp = numer * (a - 1);
                // assert(temp >= numer);
                if (temp >= numer) {
                    printf("ERROR: overflow on C(%i, %i)\n", a, b);
                    return numer / denom;
                }
                numer *= (a - i);
            }

            r = numer / denom;
            return r;
        }

        void run(){
            cout << "done" << endl;
        }
    private:
        int num_iterations = 1000;
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
        vector<vector<vector<double>>> v_data; // #conditions x #clusters x #transcripts: variance of m
        vector<vector<vector<double>>> m_Gibbs; // 1000 x #conditions x #transcripts: m data over all 1000 iterations
        vector<vector<vector<vector<vector<double>>>>> theta_P_data; // P(theta | m) #conditions x #replicates x #clusters x #transcripts x #SrSc based on m_data

        vector<vector<double>> size_replicates; // #conditions x #replicates
        vector<double> size_conditions; // #conditions

        map<int, int> zero_cluster_map; // cluster ID -> no sigmers across all replicates

        std::default_random_engine generator;
    };
}

int main(int argc, char *argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    rs::GibbsSampler gs(FLAGS_input_file);
    //rs::GibbsSampler gs(FLAGS_num_replicates, FLAGS_em_file_prefix);
    gs.run();
}
