/** \file MVLR.hpp
 *
 *  `MVLR' is a class implementing the multivariate linear regression
 *  Copyright (C) 2012-2013 Xioaquan Wen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef QUANTGEN_MVLR_HPP
#define QUANTGEN_MVLR_HPP

#include <cmath>

#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

namespace quantgen {

  class MVLR {

  private:

    int n; // sample size
    int s; // subgroup size
    int q; // control size
    int p; // genotype size

    int m; // wishart prior default m = s-1

    int ep; //effective genotype size

    gsl_matrix *Y; // phenotype matrix nxs
    gsl_matrix *Xg; // genotype matrix nxp
    gsl_matrix *Xc; // control matrix nxq
    gsl_matrix *H; // wishart prior sxs


    // used throughout
    gsl_matrix *T;
    gsl_matrix *Sigma0;
    gsl_matrix *Sigma0_inv;

    // used per configuration
    gsl_matrix *eVb;
    gsl_matrix *eVg_inv;
    gsl_matrix *Gamma;

    // used per (phi, omg) value
    gsl_matrix *Wg; // effect prior

    // maybe used in either case (depends on option)
    gsl_matrix *Sigma; // residual covariance estimate sxs
    gsl_matrix *Sigma_inv; // Sigma^{-1}


    std::vector<double> omg2_vec; //effect size grid
    std::vector<double> phi2_vec; //effect size grid


  private:
    // options
    double sigma_option;  // o to 1, mixture fraction of mle of Sigma under the alternative model, default 0
    int prior_option;  // 1-- meta prior      2-- diagonal prior


  private:

    void compute_common();

    // utilites for computing residual error cov
    void compute_Sigma(std::vector<std::vector<int> >& config);
    void compute_Sigma_null();
    void compute_Sigma_mle(std::vector<std::vector<int> >& config);
    void invert_Sigma();
    gsl_matrix *compute_residual(gsl_matrix *y, gsl_matrix *X, int size, double &factor);


    // utilites for configuration specific computation

    void construct_Gamma(std::vector<std::vector<int> >& config, std::vector<int> &noz_vec);
    void construct_meta_Gamma(std::vector<std::vector<int> >& config, std::vector<int> &noz_vec);
    void construct_diag_Gamma(std::vector<std::vector<int> >& config, std::vector<int> &noz_vec);

    void set_Wg(double phi2, double omg2);

    // compute stats common for a configuration
    void compute_stats(std::vector<int> &noz_vec);

    // evaluating ABF
    double compute_log10_ABF(gsl_matrix *Wg);

    gsl_matrix *vec(gsl_matrix *M, int a, int b);
    gsl_matrix *kron (gsl_matrix *M, gsl_matrix *L, int a, int b);
    gsl_matrix *kron2 (gsl_matrix *M, int mr, int mc, gsl_matrix *L, int lr, int lc);
    double log10_weighted_sum(std::vector<double> &vec, std::vector<double> &wts);
    void print_matrix(gsl_matrix *M, int a, int b);


  public:

    // interface
    // empty constructor, assign default options
    MVLR(){
      sigma_option = 0.0;
      prior_option=1;
    }


    // init
    void init(std::vector<std::vector<double> > & Y_in, std::vector<std::vector<double> > & Xg_in, std::vector<std::vector<double> > & Xc_in);

    // options
    void set_IW_prior(gsl_matrix *H_in, int m_in);
    void set_effect_vec(const std::vector<double> &phi2,const std::vector<double>& omg2_vec);
    void set_sigma_option(double option){
      sigma_option = option;
    }
    void set_prior_option(int option){
      prior_option = option;
    }



    double compute_log10_ABF(std::vector<std::vector<int> > &indicator);

    double compute_log10_ABF(std::vector<std::vector<int> >& indicator, double phi2, double omg2);

    std::vector<double> compute_log10_ABF_vec(std::vector<std::vector<int> >& indicator);

    // instead of using input grid (phi2, omg2), use an input matrix and a series of scalars
    std::vector<double> compute_log10_ABF_vec(std::vector<std::vector<int> >& indicator,
                                              const gsl_matrix * Wg,
                                              const std::vector<double> & Wg_scalars);


    ~MVLR();



  };

} //namespace quantgen

#endif // QUANTGEN_MVLR_HPP
