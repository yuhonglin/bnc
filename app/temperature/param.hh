#ifndef PARAM_H
#define PARAM_H

#include <cmath>
#include <dist/norm.hh>

#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>

#include <matrix/matrix.hh>

#include "data.hh"
#include "config.hh"

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
    Param(const Configuration& config, const Data& data, RNGType* rng) {

      bnc::Matrix yt = config.data.at("yt");

      bnc::Vector t_range(config.T); for (int i=1; i<=config.T; i++) t_range(i-1)=i;
      bnc::Vector cost = (t_range.array()*2*3.14159265358979323/12).cos();
      bnc::Vector sint = (t_range.array()*2*3.14159265358979323/12).sin();
      
      bnc::Matrix detrend(config.S, 4);
      bnc::Matrix X(config.T,4);
      X.col(0).array() = 1.;
      X.col(1)         = t_range;
      X.col(2)         = cost;
      X.col(3)         = sint;
      
      for (int i=0; i<config.S; i++) {
	detrend.row(i) =
	  X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
	  .solve(yt.col(i)).transpose();
      }
      
      bnc::Vector tr_coef =
	data.X_mu.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
	.solve(detrend.col(0));
      
      bnc::Vector Almp_coef =
	data.X_Almp.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
	.solve(detrend.col(2));

      bnc::Vector gamma_coef =
	data.X_gamma.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
	.solve( (detrend.col(3).array()/detrend.col(2).array()).matrix() );

      mu_s    = detrend.col(0);
      Almp_s  = detrend.col(2);
      gamma_s = detrend.col(3).array()/detrend.col(2).array();

      beta_mu    = tr_coef;
      beta_Almp  = Almp_coef;
      beta_gamma = gamma_coef;

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

      
	// // mu_s
	// mu_s = 20 + rnorm(config.S, 0., 1., rng).array();
	// // Almp_s
	// Almp_s = 6 + rnorm(config.S, 0., 1., rng).array();
	// // gamma_s
	// gamma_s = 0.6 + rnorm(config.S, 0., 1., rng).array();

	// // beta_mu
	// beta_mu = bnc::Vector(7);
	// beta_mu << 22.3574034, -1.1339522, -0.6861593, 3.2300531, -0.4381030,
	//            0.1067707, -2.1484094;
	// beta_mu += rnorm(7, 0., 1., rng);
	// // beta_Almp
	// beta_Almp = bnc::Vector(6);
	// beta_Almp << 7.4468708, -1.7941546, -1.0614744, 1.6414963, -0.6990989, 
	//              0.4470797;
	// beta_Almp += rnorm(6, 0., 1., rng);
	// // beta_gamma
	// beta_gamma = bnc::Vector(3);
	// beta_gamma << 0.61238031, -0.09456205, 0.01961614;
	// beta_gamma = beta_gamma + rnorm(3, 0., 1., rng);
	// // psi_mu
	// psi_mu = 1.398935 + rnorm(0., 1., rng);
	// // sigma2_mu
	// sigma2_mu = 1/psi_mu;
	// // psi_Almp
	// psi_Almp = std::max(2.144541 + rnorm(0., 1., rng), 0.3);
	// // sigma2_Almp
	// sigma2_Almp = 1/psi_Almp;
	// // psi_gamma
	// psi_gamma = std::max(2.096279 + rnorm(0., 1., rng), 0.3);
	// // sigma2_gamma
	// sigma2_gamma = 1/psi_gamma;
	// // lambda_mu
	// lambda_mu = std::max(1.59522 + rnorm(0., 1., rng), 0.3);
	// // lambda_Almp
	// lambda_Almp = std::max(0.7601071 + rnorm(0., 1., rng), 0.3);
	// // lambda_gamma
	// lambda_gamma = std::max(2.35 + rnorm(0., 1., rng), 0.3);
	// // psi_v
	// psi_v = std::max(4.928772 + rnorm(0., 1., rng), 0.3);
	// // sigma2_v
	// sigma2_v = 1/psi_v;

	// Phi = bnc::Vector::Zero(config.J);
	// alpha_t = Eigen::Map<bnc::Matrix>(rnorm((config.T+1)*config.J, 0, 1, rng).data(),
	// 				  config.T+1, config.J);
	// theta_t = Eigen::Map<bnc::Matrix>(rnorm((config.T+1)*config.J, 0, 1, rng).data(),
	// 				  config.T+1, config.J);

	// count_lambda_mu = 0;
	// count_lambda_Almp = 0;
	// count_lambda_gamma = 0;
	
    } // function Param

    template <class Archive>
    void serialize(Archive & ar) {
      ar(mu_s, Almp_s, gamma_s, beta_mu, beta_Almp, beta_gamma, psi_mu,
         sigma2_mu, psi_Almp, sigma2_Almp, psi_gamma, sigma2_gamma, lambda_mu,
         lambda_Almp, lambda_gamma, psi_v, sigma2_v, Phi, alpha_t, theta_t,
         count_lambda_mu, count_lambda_Almp, count_lambda_gamma);
    }

    Param(const Param& p) {
      mu_s = p.mu_s;
      Almp_s = p.Almp_s;
      gamma_s = p.gamma_s;
      beta_mu = p.beta_mu;
      beta_Almp = p.beta_Almp;
      beta_gamma = p.beta_gamma;
      // sampled params (noise)
      psi_mu = p.psi_mu;
      sigma2_mu = p.sigma2_mu;
      psi_Almp = p.psi_Almp;
      sigma2_Almp = p.sigma2_Almp;
      psi_gamma = p.psi_gamma;
      sigma2_gamma = p.sigma2_gamma;
      lambda_mu = p.lambda_mu;
      lambda_Almp = p.lambda_Almp;
      lambda_gamma = p.lambda_gamma;
      psi_v = p.psi_v;
      sigma2_v = p.sigma2_v;
      Phi = p.Phi;
      //
      alpha_t = p.alpha_t;
      theta_t = p.theta_t;

      // counters in MH steps
      count_lambda_mu = p.count_lambda_mu;
      count_lambda_Almp = p.count_lambda_Almp;
      count_lambda_gamma = p.count_lambda_gamma;
    }

}; // class param

#endif /* PARAM_H */
