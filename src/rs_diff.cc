#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "gflags/gflags.h"
#include "glog/logging.h"

using namespace std;

DEFINE_double(FDR, 0.05,
             "Expected false discovery rate.");
DEFINE_string(gibbs_file_conditionA, "",
              "The path to the file that contains Gibbs Sampling results for condition A.");
DEFINE_string(gibbs_file_conditionB, "",
              "The path to the file that contains Gibbs Sampling results for condition B.");

class transcript{
public:
    transcript(string tid, double meanA, double varianceA, double meanB, double varianceB){
        this->tid = tid;
        this->meanA = meanA;
        this->varianceA = varianceA;
        this->meanB = meanB;
        this->varianceB = varianceB;
        this->p_value = 1;
        this->adjusted_p = 1;
        this->is_significant = false;
    }
    
    string get_tid(){
        return tid;
    }
    
    double get_meanA(){
        return meanA;
    }
    
    double get_varianceA(){
        return varianceA;
    }
    
    double get_meanB(){
        return meanB;
    }
    
    double get_varianceB(){
        return varianceB;
    }
    
    double get_p_value(){
        return p_value;
    }
    
    double get_adjusted_p(){
        return adjusted_p;
    }
    
    bool get_significant(){
        return is_significant;
    }
    
    void set_p_value(double value){
        p_value = value;
    }
    
    void set_adjusted_p(double value){
        adjusted_p = value;
    }
    
    void set_significant(bool value){
        is_significant = value;
    }
    
private:
    string tid;
    double meanA;
    double varianceA;
    double meanB;
    double varianceB;
    double p_value;
    double adjusted_p;
    bool is_significant; // whether this transcript is considered as differentially expressed
};

// split a string s seperated by delim into a vector of tokens
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

// read results from gibbs sampling. 
// Gibbs sampling result for condition A is in gibbs_fileA and Gibbs sampling result for condition B is in gibbs_fileB.
void read_gibbs_result(vector<transcript*>& trans, string gibbs_fileA, string gibbs_fileB){
    fstream istreamA(gibbs_fileA, ios::in);
    fstream istreamB(gibbs_fileB, ios::in);
    //count is the number of lines being read
    int count = 0;
    string lineA, lineB;
    while (getline(istreamA, lineA)) {
        ++count;
        //line = transcript_id    mean_condition variance_condition
        vector<string> tokensA = split(lineA, '\t');
        getline(istreamB, lineB);
        vector<string> tokensB = split(lineB, '\t');
        LOG_IF(ERROR, tokensA[0] != tokensB[0])
            << "Transcript ids in line" << count << " of two conditions are different.";
        transcript* t = new transcript(tokensA[0], atof(tokensA[1].c_str()), atof(tokensA[2].c_str()),
                                       atof(tokensB[1].c_str()), atof(tokensB[2].c_str()));
        trans.push_back(t);
    }
}


// calculate p value for transcript t under two conditions
void calc_p_value(transcript* t){
    double meanA = t->get_meanA();
    double varianceA = t->get_varianceA();
    double meanB = t->get_meanB();
    double varianceB = t->get_varianceB();
    
    // if it's not expressed in both conditions, it's not differential expressed, p value = 1.
    if (meanA == 0 && meanB == 0){
        t->set_p_value(1);
    }
    // if it's expressed in both conditions, use t test.
    if (meanA != 0 && meanB != 0){
        double t_value = (log(meanA) - log(meanB))/sqrt(varianceA/(meanA*meanA) + varianceB/(meanB*meanB));
        t->set_p_value(erfc(abs(t_value)/sqrt(2)));
    }
    // if it's expressed in only one condition,
    // a one-sided test is performed using the posterior distribution of the expressed transcript.
    if (meanA == 0 && meanB != 0){
        double p_comp = meanB/varianceB;
        double r = meanB*meanB/(varianceB-meanB);
        t->set_p_value(pow(p_comp, r));
    }
    if (meanA != 0 && meanB == 0){
        double p_comp = meanA/varianceA;
        double r = meanA*meanA/(varianceA-meanA);
        t->set_p_value(pow(p_comp, r));
    }
}

// compare the p value of two transcripts
// return true if the p value of transcript a is smaller than p value of b
bool comp_trans(transcript* a, transcript* b){
    return (a->get_p_value() < b->get_p_value());
}

// multiple testing using Benjamini and Hochberg
// adjusted p value is calculated
// transcripts with adjusted p value < FDR is set to be differentially expressed.
void multiple_testing(vector<transcript*> trans, double FDR){
    sort(trans.begin(), trans.end(), comp_trans);
    int size = trans.size();
    double min = trans[size-1]->get_p_value();
    for (int i = size-1; i >= 0; --i){
        double new_q = (double)size/(i+1)*trans[i]->get_p_value();
        if (new_q < min) min = new_q;
        trans[i]->set_adjusted_p(min);
        if (min <= FDR) {
            trans[i]->set_significant(true);
        }
    }
}

int main(int argc, char *argv[]){
    google::ParseCommandLineFlags(&argc, &argv, true);
    
    double FDR = FLAGS_FDR;
    string gibbs_fileA = FLAGS_gibbs_file_conditionA;
    string gibbs_fileB = FLAGS_gibbs_file_conditionB;
    
    vector<transcript*> trans;
    read_gibbs_result(trans, gibbs_fileA, gibbs_fileB);
    
    for (transcript* t: trans){
        calc_p_value(t);
    }
    
    multiple_testing(trans, FDR);
    
    //output
    cout << "transcript_id" << '\t' << "mean_conditionA" << '\t' << "variance_conditionA" << '\t'
         << "mean_conditionB" << '\t' << "variance_conditionB" << '\t' << "p_value" << '\t'
         << "adjusted_p_value" << '\t' << "is_significant" << endl;
    for (int i = 0; i < trans.size(); ++i){
        transcript* t = trans[i];
        cout << t->get_tid() << '\t' << t->get_meanA() << '\t' << t->get_varianceA() << '\t'
             << t->get_meanB() << '\t' << t->get_varianceB() << '\t' << t->get_p_value() << '\t'
             << t->get_adjusted_p() << '\t';
        if (t->get_significant()) {
            cout << "yes";
        }
        else cout << "no";
        cout << endl;
    }
}
