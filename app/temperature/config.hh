#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <map>

struct Configuration {
public:
    int J; // number of kernels
    int T;
    double rho;
    double deps;
    double deta;
    int S; // number of station

    double outer_delta_lambda_mu;
    double outer_delta_lambda_gamma;
    double outer_delta_lambda_Almp;

    std::string data_path;
    
    std::map<std::string, bnc::Matrix> data; // store data loaded from disk
};


#endif /* _CONFIG_H_ */
