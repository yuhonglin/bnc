#ifndef DATA_H
#define DATA_H

#include <matrix/matrix.hh>

#include "config.hh"

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
	    std::cerr << "Unsupported number of kernels: " << config.J << std::endl;
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
    } // function Data
}; // class Data


#endif /* DATA_H */
