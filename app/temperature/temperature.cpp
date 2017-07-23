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
void do_sample(const int& n, const Data &data, Param &param, std::vector<Param> &store, RNGType *rng) {
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
    bnc::Matrix V(data.S, data.S);

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

    // invexpDAlmp
    bnc::Matrix invexpDAlmp(data.S, data.S);
    bnc::Matrix C_Almp_s(data.S, data.S);
    bnc::Vector m_Almp_s(data.S);

    // C.Almp.beta & m.Almp.beta
    bnc::Matrix C_Almp_beta(6, 6);
    bnc::Vector m_Almp_beta(6);

    // C.gamma.s & m.gamma.s
    bnc::Matrix invexpDgamma(data.S, data.S);
    bnc::Matrix C_gamma_s(data.S, data.S);
    bnc::Vector m_gamma_s(data.S);
    bnc::Matrix tt3(data.T, data.S);

    // C_gamma_beta & m_gamma_beta
    bnc::Matrix C_gamma_beta(3, 3);
    bnc::Vector m_gamma_beta(3);

    // sigma2.mu
    bnc::Matrix R_lambda_mu_inv(data.S, data.S);
    double shape;
    double rate;
    bnc::Vector s2mutmp(data.S);

    // for sigma2.v/psi.v
    bnc::Vector vt(data.S);
    
    for (int MC = 0; MC < n; MC++) {
	if (MC%1000==0)
	    std::cout << "MC : " << MC << std::endl;
	if (MC==1021) {
//	    ofstream opf("param.cereal");
//	    cereal::BinaryOutputArchive oarchive(opf);
//	    oarchive(param);
	    std::cout << "for break" << std::endl;
	}
	/*
	 *    Module 1: sample theta.t and alpha.t
	 */
	for (int i=0; i<data.T; i++) {
	    ct.row(i) = data.yt.row(i).array() -
		( param.mu_s.array() +
		  ( param.Almp_s.array() *
		    (data.cost(i) + data.sint(i)*param.gamma_s.array()) ) ).transpose();
	}
	// setting the model
	// only G and V are needed
	G.bottomRightCorner(data.J, data.J).diagonal() = param.Phi;
	V.diagonal().array() = param.sigma2_v;

	//COUT(ct.rightCols(3));
	
	theta_alpha
	    = dlm.sample(ct.transpose(), G, F, W, V, data.m0, data.C0, rng)
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
	T_t1 = xt - data.cost*param.Almp_s.transpose() -
	    data.sint*(param.Almp_s.array()*param.gamma_s.array()).matrix().transpose();

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
                     (data.X_mu.transpose() * invexpDmu * param.mu_s).array())
                        .matrix();

	param.beta_mu = rmvnorm(m_mu_beta, C_mu_beta, rng);

	// save(beta_mu)

	/*
	 *    Module 7: sampling Almp.s (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	// T.t2 = xt.rowwise().array() - mu.s.tranpose();
	SST = (data.sint*param.gamma_s.transpose()).colwise() + data.cost;
	invexpDAlmp = (-data.D.array() / param.lambda_Almp).exp()
	    .matrix().llt().solve(bnc::Matrix::Identity(data.S, data.S));
	C_Almp_s = param.psi_Almp*invexpDAlmp;
	C_Almp_s.diagonal() += param.psi_v *
	    SST.array().square().colwise().sum().matrix();
	C_Almp_s = C_Almp_s.llt().solve(bnc::Matrix::Identity(data.S, data.S));
	
	m_Almp_s = C_Almp_s*((param.psi_v *
			     (SST.array() *
			      (xt.rowwise()-param.mu_s.transpose()).array())
			      .colwise().sum()).matrix().transpose() +
			     param.psi_Almp*invexpDAlmp*data.X_Almp*param.beta_Almp);
	
	param.Almp_s = bnc::rmvnorm(m_Almp_s, C_Almp_s, rng);

	// save(Almp.s)

	/*
	 *    Module 8: sampling beta.Almp (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	C_Almp_beta = (param.psi_Almp *
		       data.X_Almp.transpose()*invexpDAlmp*data.X_Almp +
		       data.Sigma_beta_Almp_star_inv)
	    .llt().solve(bnc::Matrix::Identity(6, 6));
	m_Almp_beta = C_Almp_beta *
	    (param.psi_Almp * data.X_Almp.transpose()*invexpDAlmp*param.Almp_s);

	param.beta_Almp = bnc::rmvnorm(m_Almp_beta, C_Almp_beta, rng);

	// saave(beta.Almp)
	
	/*
	 *    Module 9: sampling gamma.s (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	tt3 = data.sint * param.Almp_s.transpose();
	invexpDgamma = (-data.D.array() / param.lambda_gamma).exp()
	    .matrix().llt().solve(bnc::Matrix::Identity(data.S, data.S));
	C_gamma_s = param.psi_gamma*invexpDgamma;
	C_gamma_s.diagonal() += param.psi_v *
	    tt3.array().square().matrix().colwise().sum().transpose(); //tt2
	C_gamma_s = C_gamma_s.llt().solve(bnc::Matrix::Identity(data.S, data.S));
	m_gamma_s = C_gamma_s *
	    (param.psi_v *
	     ( tt3.array() *
	       (xt.rowwise() - param.mu_s.transpose())  // T.t3
	       .array() )                               
	     .matrix().colwise().sum().transpose()      // tt2
	     + param.psi_gamma*invexpDgamma*data.X_gamma*param.beta_gamma );
	
	param.gamma_s = bnc::rmvnorm(m_gamma_s, C_gamma_s, rng);

	// save(gamma.s)
	/*
	 *    Module 10: sampling beta.gamma (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	C_gamma_beta = (param.psi_gamma *
			data.X_gamma.transpose()*invexpDgamma*data.X_gamma +
			data.Sigma_beta_gamma_star_inv)
	    .llt().solve(bnc::Matrix::Identity(3,3));
	m_gamma_beta = C_gamma_beta *
	    (param.psi_gamma * data.X_gamma.transpose()*invexpDgamma*param.gamma_s);

	param.beta_gamma = bnc::rmvnorm(m_gamma_beta, C_gamma_beta, rng);

	// save(beta.gamma)

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

	/*
	 *    Module 12: sampling sigma2.Almp (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	shape = data.a_psi_Almp + data.S/2.0;
	s2mutmp = param.Almp_s - data.X_Almp*param.beta_Almp;
	rate = data.b_psi_Almp + 0.5*(s2mutmp.transpose()*invexpDAlmp*s2mutmp)(0);
	
	param.psi_Almp = bnc::rgamma(shape, 1/rate, rng);
	param.sigma2_Almp = 1/param.psi_Almp;

	// save(psi.Almp)
	// save(sigma2.Almp)

	/*
	 *    Module 13: sampling sigma2.gamma (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	shape = data.a_psi_gamma + data.S/2.0;
	s2mutmp = param.gamma_s - data.X_gamma*param.beta_gamma;
	rate  = data.b_psi_gamma + 0.5*(s2mutmp.transpose()*invexpDgamma*s2mutmp)(0);

	param.psi_gamma = bnc::rgamma(shape, 1/rate, rng);
	param.sigma2_gamma = 1/param.psi_gamma;

	// save(psi.gamma)
	// save(sigma2.gamma)

	/*
	 *    Module 14: sampling lambda.mu, lambda.Almp and lambda.gamma
	 *               (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	/// lambda.mu
	// generate candidate
	double log_cand_lambda_mu = bnc::rnorm(std::log(param.lambda_mu), data.delta_lambda_mu, rng);
	double cand_lambda_mu     = std::exp(log_cand_lambda_mu);
	// accept probability: alpha.lambda.mu
	double numerator =
	    bnc::dmvnorm(param.mu_s, data.X_mu*param.beta_mu,
			 param.sigma2_mu*(-data.D.array()/cand_lambda_mu).exp().matrix(), bnc::LOG)
	    + bnc::dinvgamma(cand_lambda_mu, 2, 1/data.d, bnc::LOG)
	    + std::log(cand_lambda_mu);
	double denominator =
	    bnc::dmvnorm(param.mu_s, data.X_mu*param.beta_mu,
			 param.sigma2_mu*(-data.D.array()/param.lambda_mu).exp().matrix(), bnc::LOG)
	    + bnc::dinvgamma(param.lambda_mu, 2, 1/data.d, bnc::LOG)
	    + std::log(param.lambda_mu);

	double ratio = exp(numerator - denominator);
	
	double alpha_lambda_mu = std::min(1., ratio);

	// generate a random number from U(0,1)
	double U = runif(0, 1, rng);

	// updating lambda.mu
	if (U < alpha_lambda_mu) {
	    param.lambda_mu = cand_lambda_mu;
	    param.count_lambda_mu++;
	}

	// save(lambda.mu)
	
	/// lambda.Almp
	// generate candidate
	double log_cand_lambda_Almp = bnc::rnorm(std::log(param.lambda_Almp), data.delta_lambda_Almp, rng);
	double cand_lambda_Almp = std::exp(log_cand_lambda_Almp);

	// accept probability: alpha.lambda.Almp
	numerator =
	    bnc::dmvnorm(param.Almp_s,
			 data.X_Almp*param.beta_Almp,
			 param.sigma2_Almp*(-data.D.array()/cand_lambda_Almp).exp(), bnc::LOG)
	    + std::log(bnc::dinvgamma(cand_lambda_Almp, 2, 1/data.d, bnc::LOG))
	    + std::log(cand_lambda_Almp);
	denominator =
	    bnc::dmvnorm(param.Almp_s,
			 data.X_Almp*param.beta_Almp,
			 param.sigma2_Almp*(-data.D.array()/param.lambda_Almp).exp(), bnc::LOG)
	    + std::log(bnc::dinvgamma(param.lambda_Almp, 2, 1/data.d, bnc::LOG))
	    + std::log(param.lambda_Almp);
	
	ratio = std::exp(numerator - denominator);

	double alpha_lambda_Almp = std::min(1., ratio);

	// generate a random number from U(0,1)
	U = bnc::runif(0, 1, rng);

	// updating lambda.Almp
	if (U < alpha_lambda_Almp) {
	    param.lambda_Almp = cand_lambda_Almp;
	    param.count_lambda_Almp ++;
	}

	/// lambda.gamma
	// generate candidate
	double log_cand_lambda_gamma =
	    bnc::rnorm(std::log(param.lambda_gamma), data.delta_lambda_gamma, rng);
	double cand_lambda_gamma = std::exp(log_cand_lambda_gamma);
	
	numerator =
	    bnc::dmvnorm(param.gamma_s,
			 data.X_gamma*param.beta_gamma,
			 param.sigma2_gamma*(-data.D.array()/cand_lambda_gamma))
	    + bnc::dinvgamma(cand_lambda_gamma, 2, 1/data.d, bnc::LOG)
	    + std::log(cand_lambda_gamma);
	denominator =
	    bnc::dmvnorm(param.gamma_s,
			 data.X_gamma*param.beta_gamma,
			 param.sigma2_gamma*(-data.D.array()/param.lambda_gamma))
	    + bnc::dinvgamma(param.lambda_gamma, 2, 1/data.d, bnc::LOG)
	    + std::log(param.lambda_gamma);
	ratio = std::exp(numerator - denominator);

	double alpha_lambda_gamma = std::min(1., ratio);

	// generate a random number from U(0,1)
	U = bnc::runif(0, 1, rng);

	// updating lambda.gamma
	if (U < alpha_lambda_gamma) {
	    param.lambda_gamma = cand_lambda_gamma;
	    param.count_lambda_gamma ++;
	}

	// save(lambda.gamma)

	/*
	 *    Module 15: sampling sigma2.v/psi.v
	 *               (NEED TO CHECK ERRORS in EQUATIONS)
	 */
	double temp1_v = 0;
	for (int i = 0; i < data.T; i++) {
	    // Mt is just an identity matrix
	    vt.transpose() =
		data.yt.row(i)
		- (data.scaled_weight * (param.theta_t.row(i+1)+param.alpha_t.row(i+1)).transpose()).transpose();
		- param.mu_s.transpose()
		- (param.Almp_s.array() * (data.cost(i) + data.sint(i)*param.gamma_s.array())).matrix().transpose();
	    temp1_v += vt.array().square().sum();
	}

	shape = data.a_psi_v + (data.S*data.T)/2.;
	rate  = data.b_psi_v + 0.5*temp1_v;

	param.psi_v = bnc::rgamma(shape, 1/rate, rng);
	param.sigma2_v = 1/param.psi_v;

	// save(psi.v)
	// save(sigma2.v)
	
	store.push_back(param);
//	COUT(param.lambda_mu);
//	COUT(param.lambda_Almp);
//	COUT(param.lambda_gamma);
    } // MCMC loop
}

Param * global_param_ptr;
std::vector<Param> * global_store_ptr;

int main(int argc, char *argv[])
{
    // RNG
    bnc::MersenneTwisterRNG<bnc::AHRENS_DIETER> rng(10);
    
    // do configuration
    Configuration config;
    config.J = 4; // number of kernels
    config.rho = 0.0125;
    config.deps = 0.2;
    config.deta = 0.2;
    config.T = 528; // use first T observation (max:528)
    config.S = 180; // use first S station (max:180)
    config.outer_delta_lambda_mu = .1;
    config.outer_delta_lambda_gamma = .1;
    config.outer_delta_lambda_Almp = 0.1;
    config.data_path = "./data/";

    // from init.param(1)$data$yt
    config.data["yt"] = read_csv<bnc::Matrix>(config.data_path + "/yt.csv")
	.topLeftCorner(config.T, config.S);
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

    int N = 50000;
    Data  data(config);    
//    Param param(config, data, &rng);
    Param param;
    std::ifstream pif("./tailparam_error.cereal");    
    cereal::BinaryInputArchive iarchive(pif);     
    iarchive(param);

    
    std::vector<Param> store; store.reserve(N);

//    global_param_ptr = &param;
//    global_store_ptr = &store;
//    set logger level
//    make warning throw then enable debug
//    auto cbsave = [](const char *message) {
//	ofstream pof("./store_error.cereal");
//	cereal::BinaryOutputArchive oarchive(pof);
//	oarchive(*global_store_ptr);
//	pof.flush();
//	std::cout << message << std::endl;
//	throw runtime_error(message);
//    };

//    bnc::Logger::showFatalCallback = cbsave;    
//    bnc::Logger::showWarningCallback = cbsave;
//    bnc::Logger::showErrorCallback = cbsave;
//    bnc::Logger::showStatusCallback = cbsave;
//    bnc::Logger::showMessageCallback = cbsave;

    bnc::tic();
    do_sample(N, data, param, store, &rng);
    bnc::toc();
    
    ofstream pof("./store.cereal");
    cereal::BinaryOutputArchive oarchive(pof);
    oarchive(store);
    
    return 0;
}
