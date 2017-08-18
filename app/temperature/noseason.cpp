/**
 * @file   noseason.cpp
 * @author Honglin Yu <honglin.yu@anu.edu.au>
 * @date   Fri Aug  4 20:24:18 2017
 * 
 * @brief  implement the temperature file without seasonality part.
 * 
 * 
 */


#include <iostream>
#include <vector>
#include <map>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <exception>

#include <matrix/matrix.hh>
#include <dist/norm.hh>
#include <dist/gamma.hh>
#include <dist/invgamma.hh>
#include <dist/unif.hh>
#include <rng/rng.hh>
#include <dist/tmvnorm.hh>
#include <dist/wishart.hh>
#include <util/logger.hh>
#include <util/profiling.hh>
#include <block/ssm/dlm.hh>
#include <io/text.hh>
#include <util/callstack.hh>
#include <util/rembed.hh>
#include <parallel/ThreadPool.hh>

#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <Eigen/Dense>
#include <fstream>

using namespace std;

#include <Eigen/Dense>
#include <vector>
#include <fstream>

using namespace Eigen;

#include "config.hh"
#include "data.hh"
#include "param.hh"

#define COUT(X) std::cout << '\n' << #X << '\n' << X;
//#define COUT(X) bnc::io::text::dump(X, #X);

template<typename M>
M read_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    if(!indata.good()) {
	std::cerr << "Open file '" << path << "' failed." << std::endl;
	exit(1);
    }
    std::string line;
    std::vector<double> values;
    int rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}

bnc::Matrix crossprod(const bnc::Matrix& M) {
    return M.transpose()*M;
}


template<class RNGType>
void do_sample(const int& tid, const int& n, const Data &data, Param &param, std::vector<Param> &store, RNGType *rng) {
    bnc::Matrix ct(data.T,data.S);
    // observation matrix
    bnc::Matrix F(data.S,data.J*2); // constant during sampling
    F.leftCols(data.J) = data.scaled_weight;
    F.rightCols(data.J) = data.scaled_weight;
    // dynamic matrix
    bnc::Matrix G = bnc::Matrix::Zero(data.J*2, data.J*2);
    G.topLeftCorner(data.J, data.J).diagonal().array() = 1.;
    // dynamic variance matrix
    bnc::Matrix W = bnc::Matrix::Zero(data.J*2, data.J*2); // constant during sampling
    W.topLeftCorner(data.J, data.J) = data.Sigma_epsilon;
    W.bottomRightCorner(data.J, data.J) = data.Sigma_eta;
    // observation variance matrix
    bnc::Matrix V(data.S, data.S); V.array() = 0.;

    // inv.Sigma.eta, computed for speed
    bnc::Matrix inv_Sigma_eta = data.Sigma_eta.inverse();

    // temp1 and temp2
    bnc::Matrix temp1(data.J, data.J);
    bnc::Vector temp2(data.J);
    
    // theta_alpha: sample from dlm
    bnc::Matrix theta_alpha(data.T+1, data.J*2);
    bnc::DLM dlm;

    // C_Phi and m_Phi
    bnc::Matrix C_Phi(data.J, data.J);
    bnc::Vector m_Phi(data.J);

    // upper bound and lower bound
    bnc::Vector u_Phi(data.J); u_Phi.array() = 1;
    bnc::Vector l_Phi(data.J); l_Phi.array() = -1;

    // xt
    bnc::Matrix xt(data.yt.rows(), data.yt.cols());

    // T.t1
    bnc::Matrix T_t1(data.yt.rows(), data.yt.cols());

    // tt2
    bnc::Matrix tt2(1, data.yt.cols());

    // C.mu.s
    bnc::Matrix C_mu_s(data.S, data.S);
    bnc::Matrix cmustmp(data.S, data.S);
    bnc::Matrix invexpDmu(data.S, data.S);

    // m.mu.s
    bnc::Vector m_mu_s(data.S);

    // C.mu.beta
    bnc::Matrix C_mu_beta(7,7);

    // m.mu.beta
    bnc::Vector m_mu_beta(7);

    // SST.NArm
    bnc::Matrix SST(data.T, data.S);

    // sigma2.mu
    bnc::Matrix R_lambda_mu_inv(data.S, data.S);
    double shape;
    double rate;
    bnc::Vector s2mutmp(data.S);

    // for sigma2.v/psi.v
    bnc::Vector vt(data.S);
    
    for (int MC = 0; MC < n; MC++) {
	/*
	 *    Module 1: sample theta.t and alpha.t
	 */
	for (int i=0; i<data.T; i++) {
	    ct.row(i) = data.yt.row(i).array() - param.mu_s.array().transpose();
	}
	// setting the model
	// only G and V are needed
	G.bottomRightCorner(data.J, data.J).diagonal() = param.Phi;
	V.diagonal().array() = param.sigma2_v;

//	theta_alpha
//	    = dlm.sample(ct.transpose(), G, F, W, V, data.m0, data.C0, rng).transpose();

	theta_alpha
	    = dlm.sample(ct.transpose(), G, F, W, V, bnc::Vector::Zero(data.m0.size()), bnc::Matrix::Zero(data.C0.rows(), data.C0.cols()), rng)
	    .transpose();
	

	
	param.theta_t = theta_alpha.leftCols(data.J);
	param.alpha_t = theta_alpha.rightCols(data.J);

	// save(theta.t)
	// save(alpha.t)
	/*
	 *    Module 2: sampling Phi
	 */
	temp1.array() = 0.;
	temp2.array() = 0.;
	for (int i=0; i<data.T; i++) {
	    temp1 += param.alpha_t.row(i).asDiagonal() * inv_Sigma_eta *
		param.alpha_t.row(i).asDiagonal();
	    temp2 += param.alpha_t.row(i).asDiagonal() * inv_Sigma_eta *
		param.alpha_t.row(i+1).transpose();
	}

	C_Phi = (temp1 + data.C_Phi_star_inv).llt()
	    .solve(bnc::Matrix::Identity(data.J, data.J));
	m_Phi = C_Phi * temp2;

	param.Phi = bnc::rtmvnorm(1, m_Phi, C_Phi, l_Phi, u_Phi, rng);

	// save(param.Phi)

	/*
	 *    Module 5: sampling mu.s (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	xt = data.yt -
	    (param.theta_t.bottomRows(data.T)+param.alpha_t.bottomRows(data.T))
	    * data.scaled_weight.transpose();
	T_t1 = xt;

        // tt1 = diag(data.T*1')
	tt2 = T_t1.colwise().sum().transpose();
	invexpDmu = (-data.D.array() / param.lambda_mu).exp()
	    .matrix().llt().solve(bnc::Matrix::Identity(data.S, data.S));

	// C.mu.s
	cmustmp = (param.psi_mu * invexpDmu.array()).matrix();
	cmustmp.diagonal().array() += param.psi_v * data.T;

        C_mu_s = cmustmp.llt().solve(bnc::Matrix::Identity(data.S, data.S));

	// m.mu.s
	m_mu_s = C_mu_s * (
	    param.psi_v*tt2.array() +
	    param.psi_mu*(invexpDmu*data.X_mu*param.beta_mu).array()
	    ).matrix();
	
	// mu_s
	param.mu_s = rmvnorm(m_mu_s, C_mu_s, rng);

	// save(mu_s)

	/*
	 *    Module 6: sampling mu.beta (NEED TO CHECK ERRORS in EQUATIONS)
	 */
        C_mu_beta = ((param.psi_mu *
                      (data.X_mu.transpose() * invexpDmu * data.X_mu).array())
		     .matrix() +
                     data.Sigma_beta_mu_star_inv)
	    .llt()
	    .solve(bnc::Matrix::Identity(7,7));
        m_mu_beta = C_mu_beta * 
	    (param.psi_mu *
	     (data.X_mu.transpose() * invexpDmu * param.mu_s).array()).matrix();

	param.beta_mu = rmvnorm(m_mu_beta, C_mu_beta, rng);

	// save(beta_mu)

	/*
	 *    Module 11: sampling sigma2.mu (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	// R.lambda.mu.inv == invexpDmu
	shape = data.a_psi_mu + data.S/2.0;
	s2mutmp = param.mu_s - data.X_mu*param.beta_mu;
	rate  = data.b_psi_mu + 0.5*(s2mutmp.transpose()*invexpDmu*s2mutmp)(0);

	param.psi_mu = bnc::rgamma(shape, 1/rate, rng);
	param.sigma2_mu = 1/param.psi_mu;

	// save(psi.mu)
	// save(sigma2.mu)

	// /*
	//  *    Module 14: sampling lambda.mu
	//  *               (NEED TO CHECK ERRORS in EQUATIONS)
	//  */
	// /// lambda.mu
	// // generate candidate
	// double log_cand_lambda_mu = bnc::rnorm(std::log(param.lambda_mu), data.delta_lambda_mu, rng);
	// double cand_lambda_mu     = std::exp(log_cand_lambda_mu);
	// // accept probability: alpha.lambda.mu
	// double numerator =
	//     bnc::dmvnorm(param.mu_s, data.X_mu*param.beta_mu,
	// 		 param.sigma2_mu*(-data.D.array()/cand_lambda_mu).exp().matrix(), bnc::LOG)
	//     + bnc::dinvgamma(cand_lambda_mu, 2, 1/data.d, bnc::LOG)
	//     + std::log(cand_lambda_mu);
	// double denominator =
	//     bnc::dmvnorm(param.mu_s, data.X_mu*param.beta_mu,
	// 		 param.sigma2_mu*(-data.D.array()/param.lambda_mu).exp().matrix(), bnc::LOG)
	//     + bnc::dinvgamma(param.lambda_mu, 2, 1/data.d, bnc::LOG)
	//     + std::log(param.lambda_mu);

	// double ratio = exp(numerator - denominator);
	
	// double alpha_lambda_mu = std::min(1., ratio);

	// // generate a random number from U(0,1)
	// double U = runif(0, 1, rng);

	// // updating lambda.mu
	// if (U < alpha_lambda_mu) {
	//     param.lambda_mu = cand_lambda_mu;
	//     param.count_lambda_mu++;
	// }

	param.lambda_mu = 0.1;
	
	// save(lambda.mu)
	
	/*
	 *    Module 15: sampling sigma2.v/psi.v
	 *               (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	double temp1_v = 0;
	for (int i = 0; i < data.T; i++) {
	    // Mt is just an identity matrix
	    vt.transpose() =
		data.yt.row(i)
		- (data.scaled_weight * (param.theta_t.row(i+1)+param.alpha_t.row(i+1)).transpose()).transpose()
		- param.mu_s.transpose();
	    temp1_v += vt.array().square().sum();
	}

	shape = data.a_psi_v + (data.S*data.T)/2.;
	rate  = data.b_psi_v + 0.5*temp1_v;

	param.psi_v = bnc::rgamma(shape, 1/rate, rng);
	param.sigma2_v = 1/param.psi_v;

	// save(psi.v)
	// save(sigma2.v)
	store.push_back(param);
    } // MCMC loop
}

Param * global_param_ptr;
std::vector<Param> * global_store_ptr;

#define ADD_NUM_TO_MAP(store, data, themap)		\
    themap[#data] = bnc::Vector(store.size());		\
    for (int i=0; i<store.size(); i++)			\
	themap[#data](i) = store[i].data;

#define ADD_VECTOR_TO_MAP(store, data, themap)		\
    for (int j=0; j < store[0].data.size(); j++) {	\
        std::ostringstream nms;				\
	nms << #data << "[" << j << "]";		\
	themap[nms.str()] = bnc::Vector(store.size());	\
	for (int i=0; i<store.size(); i++)		\
	    themap[nms.str()](i) = store[i].data[j];	\
    }

#define ADD_MATRIX_TO_MAP(store, data, themap)			\
    for (int k=0; k < store[0].data.rows(); k++) {	        \
        for (int j=0; j < store[0].data.cols(); j++) {		\
	    std::ostringstream nms;				\
	    nms << #data << "[" << k << "," << j << "]";	\
	    themap[nms.str()] = bnc::Vector(store.size());	\
	    for (int i=0; i<store.size(); i++)			\
		themap[nms.str()](i) = store[i].data(k,j);	\
	}							\
    }								


void store_to_map(const std::vector<Param>& store,
		  std::map<std::string, bnc::Vector>& m) {
    ADD_NUM_TO_MAP(store, sigma2_v, m);
    ADD_NUM_TO_MAP(store, psi_mu, m);
    ADD_NUM_TO_MAP(store, psi_mu            , m);
    ADD_NUM_TO_MAP(store, sigma2_mu         , m);
    ADD_NUM_TO_MAP(store, lambda_mu         , m);
    ADD_NUM_TO_MAP(store, psi_v            , m);    
    ADD_NUM_TO_MAP(store, sigma2_v, m    );

    ADD_NUM_TO_MAP(store, count_lambda_mu, m);

    ADD_VECTOR_TO_MAP(store, mu_s      , m    );
    ADD_VECTOR_TO_MAP(store, beta_mu   , m    );
    ADD_VECTOR_TO_MAP(store, Phi, m    );    

    ADD_MATRIX_TO_MAP(store, alpha_t, m    );
    ADD_MATRIX_TO_MAP(store, theta_t, m    );    
}



int main(int argc, char *argv[])
{
    //
    auto R = bnc::REmbed::get_instance(argc, argv);

    // RNG
    bnc::MersenneTwisterRNG<bnc::AHRENS_DIETER> rng(10);
    
    // do configuration
    Configuration config;
    config.J = 4; // number of kernels
    config.rho = 0.0125;
    config.deps = 0.2;
    config.deta = 0.2;
    config.T = 528-12; // use first T observation (max:528)
    config.S = 180; // use first S station (max:180)
    config.outer_delta_lambda_mu = .1;
    config.data_path = "./data/";

    // from init.param(1)$data$yt
    config.data["yt"] = read_csv<bnc::Matrix>(config.data_path + "/yt.csv")
	.topLeftCorner(config.T+12, config.S);

    // differencing the data
    config.data["yt"] = config.data["yt"].bottomRows(config.T) -
	config.data["yt"].topRows(config.T);

    // from kernel.locat[,1]
    config.data["x_kernel"] = read_csv<bnc::Matrix>(config.data_path + "/x_kernel.csv");
    // from kernel.locat[,2]
    config.data["y_kernel"] = read_csv<bnc::Matrix>(config.data_path + "/y_kernel.csv");
    // from xy.model[,1]
    config.data["x_model"] = read_csv<bnc::Matrix>(config.data_path + "/x_model.csv")
	.topRows(config.S);
    // from alt.model
    config.data["alt_model"] = read_csv<bnc::Matrix>(config.data_path + "/alt_model.csv")
	.topRows(config.S);
    // from xy.model[,2]
    config.data["y_model"] = read_csv<bnc::Matrix>(config.data_path + "/y_model.csv")
	.topRows(config.S);

    // normalisation
    {
	double x_mean = config.data["x_model"].mean();
	double x_range = config.data["x_model"].maxCoeff() -
	    config.data["x_model"].minCoeff();
	config.data["x_model"] = (config.data["x_model"].array() - x_mean)/x_range;
	config.data["x_kernel"] = (config.data["x_kernel"].array() - x_mean)/x_range;
    }

    {
	double y_mean = config.data["y_model"].mean();
	double y_range = config.data["y_model"].maxCoeff()
	    - config.data["y_model"].minCoeff();
	config.data["y_model"] = (config.data["y_model"].array() - y_mean)/y_range;
	config.data["y_kernel"] = (config.data["y_kernel"].array() - y_mean)/y_range;
    }

    {
	double alt_mean = config.data["alt_model"].mean();
	double alt_range = config.data["alt_model"].maxCoeff() -
	    config.data["alt_model"].minCoeff();
	config.data["alt_model"] = (config.data["alt_model"].array() - alt_mean) /
	    alt_range;
    }

    int totalN = 0;
    Data  data(config);

    int numThread = 3;
    std::vector<std::vector<Param>> storeVector;
    for (int i=0; i<3; i++) {
	storeVector.push_back(std::vector<Param>());
	storeVector[i].push_back(Param(config, data, &rng));
    }

    std::vector<Param> paramVector;

    std::vector<int> seeds = {100, 110, 101};
    std::vector<bnc::MersenneTwisterRNG<bnc::AHRENS_DIETER>> rngVector;
    for (auto i : seeds) {
	rngVector.push_back(bnc::MersenneTwisterRNG<bnc::AHRENS_DIETER>(i));
    }
    
    // expose variables to R
    while(true) {
	string cmd;
	std::cout << "input number of samples to continue : ";
	cin >> cmd;

	if (cmd == "q" or cmd == "quit") {
	    break;
	}
	
	int N;
	try {
	    N = std::stoi(cmd);	    
	} catch(...) {
	    std::cout << "invalid input" << std::endl;
	    continue;
	}

	// collect last param
	paramVector.clear();
	for (int i = 0; i< numThread; i++) {
	    paramVector.push_back(*(--storeVector[i].end()));
	    storeVector[i].clear();
	    storeVector[i].reserve(N);
	}

	ThreadPool pool(numThread);
	std::vector< std::future<void> > results;
	for (int i = 0; i < numThread; i++) {
	    results.emplace_back(
		pool.enqueue([i,&data,&paramVector,&storeVector,N, &rngVector] {
			do_sample(i, N, data, paramVector[i], storeVector[i], &rngVector[i]);
		    })
		);
	}

	for (auto& r : results) r.get();

	std::vector<std::map<std::string, bnc::Vector>> mapVector;
	for (int i=0; i<3; i++) {
	    std::map<std::string, bnc::Vector> m;
	    store_to_map(storeVector[i], m);
	    mapVector.push_back(m);
	}

	R.eval("rm(list='res')");
	R.eval("rm(list='yt')");
	R.define_coda_mcmc_list("res", mapVector);
	R.define_var("yt", data.yt);

	totalN += N;
	std::cout << "totalN = " << totalN << std::endl;
	
	R.repl();
	
    }
    
    return 0;
}
