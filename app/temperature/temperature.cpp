#include <iostream>
#include <vector>
#include <map>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <sstream>

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

using namespace std;

//#include <dist/dist.hh>

#include <Eigen/Dense>
#include <vector>
#include <fstream>

using namespace Eigen;

#define COUT(X) std::cout << '\n' << #X << '\n' << X;

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
    uint rows = 0;
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


/*
 *    Configuration
 */
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

    string data_path;
    
    std::map<std::string, bnc::Matrix> data; // store data loaded from disk
};


class Param {
public:
    // change from fitting    
    bnc::Vector mu_s;
    bnc::Vector Almp_s;
    bnc::Vector gamma_s;
    bnc::Vector beta_mu;
    bnc::Vector beta_Almp;
    bnc::Vector beta_gamma;
    // sampled params (add random noise)
    double psi_mu;
    double sigma2_mu;
    double psi_Almp;
    double sigma2_Almp;
    double psi_gamma;
    double sigma2_gamma;
    double lambda_mu;
    double lambda_Almp;
    double lambda_gamma;
    double psi_v;
    double sigma2_v;
    bnc::Vector Phi;
    //
    bnc::Matrix alpha_t;
    bnc::Matrix theta_t;

    // counters in MH steps
    int count_lambda_mu;
    int count_lambda_Almp;
    int count_lambda_gamma;

    template<class RNGType>
    Param(const Configuration& config, RNGType* rng) {
	// mu_s
	mu_s = 20 + rnorm(config.S, 0., 1., rng).array();
	// Almp_s
	Almp_s = 6 + rnorm(config.S, 0., 1., rng).array();
	// gamma_s
	gamma_s = 0.6 + rnorm(config.S, 0., 1., rng).array();
	// beta_mu
	beta_mu = bnc::Vector(7);
	beta_mu << 22.3574034, -1.1339522, -0.6861593, 3.2300531, -0.4381030,
	           0.1067707, -2.1484094;
	beta_mu += rnorm(7, 0., 1., rng);
	// beta_Almp
	beta_Almp = bnc::Vector(6);
	beta_Almp << 7.4468708, -1.7941546, -1.0614744, 1.6414963, -0.6990989, 
	             0.4470797;
	beta_Almp += rnorm(6, 0., 1., rng);
	// beta_gamma
	beta_gamma = bnc::Vector(3);
	beta_gamma << 0.61238031, -0.09456205, 0.01961614;
	beta_gamma = beta_gamma + rnorm(3, 0., 1., rng);
	// psi_mu
	psi_mu = 1.398935 + rnorm(0., 1., rng);
	// sigma2_mu
	sigma2_mu = 1/psi_mu;
	// psi_Almp
	psi_Almp = std::max(2.144541 + rnorm(0., 1., rng), 0.3);
	// sigma2_Almp
	sigma2_Almp = 1/psi_Almp;
	// psi_gamma
	psi_gamma = std::max(2.096279 + rnorm(0., 1., rng), 0.3);
	// sigma2_gamma
	sigma2_gamma = 1/psi_gamma;
	// lambda_mu
	lambda_mu = std::max(1.59522 + rnorm(0., 1., rng), 0.3);
	// lambda_Almp
	lambda_Almp = std::max(0.7601071 + rnorm(0., 1., rng), 0.3);
	// lambda_gamma
	lambda_gamma = std::max(2.35 + rnorm(0., 1., rng), 0.3);
	// psi_v
	psi_v = std::max(4.928772 + rnorm(0., 1., rng), 0.3);
	// sigma2_v
	sigma2_v = 1/psi_v;

	Phi = bnc::Vector::Zero(config.J);
	alpha_t = Eigen::Map<bnc::Matrix>(rnorm((config.T+1)*config.J, 0, 1, rng).data(),
					  config.T+1, config.J);
	theta_t = Eigen::Map<bnc::Matrix>(rnorm((config.T+1)*config.J, 0, 1, rng).data(),
					  config.T+1, config.J);

	count_lambda_mu = 0;
	count_lambda_Almp = 0;
	count_lambda_gamma = 0;
	
    }
};


class Data {
public:
    // transformed data (NOT change)
    bnc::Matrix scaled_weight;
    int J;
    int S;
    int T;
    bnc::Matrix D;
    double d;
    bnc::Matrix X_Almp;
    bnc::Matrix X_gamma;
    bnc::Matrix X_mu;
    bnc::Matrix yt;
    bnc::Vector cost;
    bnc::Vector sint;
    // prior (NOT change)
    int degree_eta;
    bnc::Matrix V_eta_prior;
    bnc::Matrix m0;
    bnc::Matrix C0;
    int degree_epsilon;
    bnc::Matrix V_epsilon_prior;
    bnc::Matrix C_Phi_star_inv;
    bnc::Matrix Sigma_beta_mu_star_inv;
    bnc::Matrix Sigma_beta_Almp_star_inv;
    bnc::Matrix Sigma_beta_gamma_star_inv;
    double a_psi_v;
    double b_psi_v;
    double a_psi_mu;
    double b_psi_mu;
    double a_psi_Almp;
    double b_psi_Almp;
    double a_psi_gamma;
    double b_psi_gamma;
    double delta_lambda_mu;
    double delta_lambda_Almp;
    double delta_lambda_gamma;
    bnc::Matrix Sigma_epsilon;
    bnc::Matrix Sigma_eta;

    Data(const Configuration& config) {

	// number of kernels
	J = config.J;
	// number of stations
	S = config.S;
	// number of time steps
	T = config.T;
	
	// load x-coordinate of stations
	bnc::Vector x_model = config.data.at("x_model");
	// load y-coordinate of stations
	bnc::Vector y_model = config.data.at("y_model");
	// load altitude
	bnc::Vector alt_model = config.data.at("alt_model");

	// compute D matrix
	D = bnc::Matrix(S,S);
	for (int i=0; i<S; i++)
	    for (int j=0; j<S; j++)
		D(i,j) = std::sqrt( pow(x_model(i)-x_model(j),2) +
				    pow(y_model(i)-y_model(j),2) );
	
	// compute d
	d = D.maxCoeff()/(-2*std::log(0.05));
	
	// X_Almp
	X_Almp = bnc::Matrix(S, 6);
	X_Almp.col(0).array() = 1.;
	X_Almp.col(1) = x_model;
	X_Almp.col(2) = x_model.array().square();
	X_Almp.col(3) = y_model;
	X_Almp.col(4) = y_model.array().square();
	X_Almp.col(5) = x_model.array()*y_model.array();

	// X_gamma
	X_gamma = bnc::Matrix(S, 3);
	X_gamma.col(0).array() = 1.;
	X_gamma.col(1) = y_model;
	X_gamma.col(2) = y_model.array().square();

	// X_mu
	X_mu = bnc::Matrix(S, 7);
	X_mu.col(0).array() = 1.;
	X_mu.col(1) = x_model;
	X_mu.col(2) = x_model.array().square();
	X_mu.col(3) = y_model;
	X_mu.col(4) = y_model.array().square();
	X_mu.col(5) = x_model.array()*y_model.array();
	X_mu.col(6) = alt_model;

	// load yt
	yt = config.data.at("yt");

	// cost and sint
	bnc::Vector t_range(T); for (int i=1; i<=T; i++) t_range(i-1)=i;
	cost = (t_range.array()*2*3.14159265358979323/12).cos();
	sint = (t_range.array()*2*3.14159265358979323/12).sin();

	// degree_eta
	degree_eta = J+1;

	// V_eta_prior
	V_eta_prior = bnc::Matrix::Identity(J,J)*(J+1);

	// m0
	m0 = bnc::Vector::Zero(2*J);

	// C0
	C0 = bnc::Matrix::Identity(2*J,2*J);
	C0.diagonal().head(J).array() = 0.0000001;
	C0.diagonal().tail(J).array() = 100;

	// degree_epsilon
	degree_epsilon = J+1;

	// V_epsilon_prior
	V_epsilon_prior = bnc::Matrix::Identity(J,J)*(J+1);

	// C_Phi_star_inv
	C_Phi_star_inv = bnc::Matrix::Identity(J,J)*0.01;

	// Sigma_beta_mu_star
	Sigma_beta_mu_star_inv = (bnc::Matrix::Identity(7,7)*100).inverse();

	// Sigma_beta_Almp_star_inv
	Sigma_beta_Almp_star_inv = (bnc::Matrix::Identity(6,6)*100).inverse();

	//Sigma_beta_gamma_star_inv
	Sigma_beta_gamma_star_inv = (bnc::Matrix::Identity(3,3)*100).inverse();

	// scalars
	a_psi_v      = 0.01;
	b_psi_v      = 0.01;
	a_psi_mu     = 0.01;
	b_psi_mu     = 0.01;
	a_psi_Almp   = 0.01;
	b_psi_Almp   = 0.01;
	a_psi_gamma  = 0.01;
	b_psi_gamma  = 0.01;
	delta_lambda_mu    = config.outer_delta_lambda_mu;
	delta_lambda_Almp  = config.outer_delta_lambda_Almp;
	delta_lambda_gamma = config.outer_delta_lambda_gamma;

	// compute K matrix
	bnc::Matrix x_kernel = config.data.at("x_kernel");
	bnc::Matrix y_kernel = config.data.at("y_kernel");
	bnc::Matrix K(J,J);
	for (int i=0; i<J; i++)
	    for (int j=0; j<J; j++)
		K(i,j) = std::sqrt(std::pow(x_kernel(i)-x_kernel(j),2) +
				   std::pow(y_kernel(i)-y_kernel(j),2));
	// Sigma_epsilon
	Sigma_epsilon = (-K.array()/config.deps).exp();

	// Sigma_eta
	Sigma_eta = (-K.array()/config.deta).exp();

	// scaled_weight
	bnc::Matrix model_weight(config.S, config.J);
	double kernel_dist;
	switch (config.J) {
	case 4:
	    kernel_dist = 600;
	    break;
	case 7:
	    kernel_dist = 500;
	    break;
	case 12:
	    kernel_dist = 400;
	    break;
	case 18:
	    kernel_dist = 275;
	    break;
	default:
	    cerr << "Unsupported number of kernels: " << config.J << endl;
	    exit(1);
	}
	for (int j=0; j<config.J; j++) {
	    bnc::Matrix kernel_shape = bnc::Matrix::Identity(2,2).array() *
		pow(1/(kernel_dist*config.rho),2);
	    for (int i=0; i<config.S; i++) {
		bnc::Vector tmp(2);
		tmp(0) = x_model(i) - x_kernel(j);
		tmp(1) = y_model(i) - y_kernel(j);
		model_weight(i,j) = std::exp(-0.5*(tmp.transpose()*kernel_shape*tmp)(0));
	    }
	}
	scaled_weight = bnc::Matrix(config.S, config.J);
	for (int i=0; i<config.S; i++) {
	    scaled_weight.row(i) = model_weight.row(i).array() /
		model_weight.row(i).sum();
	}
    }
};


bnc::Matrix crossprod(const bnc::Matrix& M) {
    return M.transpose()*M;
}


template<class RNGType>
void do_sample(const int& n, const Data &data, Param &param, RNGType *rng) {
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
        // tt1 = data.T*1'
	tt2 = T_t1.colwise().sum().transpose().array()*data.T;
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

    } // MCMC loop
}


int main(int argc, char *argv[])
{
    // RNG
    bnc::MersenneTwisterRNG<bnc::AHRENS_DIETER> rng(10);
    
    // do configuration
    Configuration config;
    config.J = 4; // number of kernels
    config.rho = 0.125;
    config.deps = 2.0;
    config.deta = 2.0;
    config.T = 528; // use first T observation (max:528)
    config.S = 180; // use first S station (max:180)
    config.outer_delta_lambda_mu = .9;
    config.outer_delta_lambda_gamma = 1.;
    config.outer_delta_lambda_Almp = 0.8;
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

    Param param(config, &rng);
    Data  data(config);

    bnc::tic();
    do_sample(10, data, param, &rng);
    bnc::toc();
	
    return 0;
}
