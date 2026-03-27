#include <RcppArmadillo.h>

#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

inline void shuffle_uvec_inplace(arma::uvec& x) {
  if (x.n_elem <= 1) {
    return;
  }
  for (arma::uword i = x.n_elem - 1; i > 0; --i) {
    arma::uword j = static_cast<arma::uword>(std::floor(R::unif_rand() * (i + 1)));
    std::swap(x(i), x(j));
  }
}

// [[Rcpp::export(name = ".permutation_ssq_kernel_cpp")]]
Rcpp::NumericVector permutation_ssq_kernel_cpp(const arma::mat& DD,
                                               const arma::mat& DR,
                                               const Rcpp::List& groups,
                                               const int n_perm) {
  Rcpp::RNGScope scope;

  const int n_groups = groups.size();
  std::vector<arma::uvec> idx_store;
  idx_store.reserve(n_groups);

  for (int g = 0; g < n_groups; ++g) {
    Rcpp::IntegerVector block = groups[g];
    if (block.size() <= 1) {
      continue;
    }
    arma::uvec idx(block.size());
    for (int i = 0; i < block.size(); ++i) {
      idx(i) = static_cast<unsigned int>(block[i] - 1);
    }
    idx_store.push_back(idx);
  }

  const arma::uword N = DR.n_rows;
  arma::uvec base_order = arma::regspace<arma::uvec>(0, N - 1);
  arma::uvec order = base_order;
  arma::uvec perm_idx;

  Rcpp::NumericVector vals(n_perm);
  for (int p = 0; p < n_perm; ++p) {
    order = base_order;

    for (std::size_t g = 0; g < idx_store.size(); ++g) {
      perm_idx = idx_store[g];
      shuffle_uvec_inplace(perm_idx);
      order.elem(idx_store[g]) = perm_idx;
    }

    arma::mat score = DD * DR.rows(order);
    vals[p] = arma::accu(score % score);
  }

  return vals;
}

// [[Rcpp::export(name = ".rotation_ssq_kernel_cpp")]]
Rcpp::NumericVector rotation_ssq_kernel_cpp(const arma::mat& DD,
                                            const arma::mat& DR,
                                            const Rcpp::List& groups,
                                            const int n_rot) {
  Rcpp::RNGScope scope;

  const int n_groups = groups.size();
  std::vector<arma::uvec> idx_store;
  idx_store.reserve(n_groups);

  for (int g = 0; g < n_groups; ++g) {
    Rcpp::IntegerVector block = groups[g];
    if (block.size() <= 1) {
      continue;
    }
    arma::uvec idx(block.size());
    for (int i = 0; i < block.size(); ++i) {
      idx(i) = static_cast<unsigned int>(block[i] - 1);
    }
    idx_store.push_back(idx);
  }

  Rcpp::NumericVector vals(n_rot);
  arma::mat DR_rot;
  arma::mat Q;
  arma::mat R;

  for (int r = 0; r < n_rot; ++r) {
    DR_rot = DR;

    for (std::size_t g = 0; g < idx_store.size(); ++g) {
      const arma::uvec& idx = idx_store[g];
      const arma::uword m = idx.n_elem;

      arma::mat Z = arma::randn<arma::mat>(m, m);
      arma::qr_econ(Q, R, Z);

      arma::vec d = arma::sign(R.diag());
      d.elem(arma::find(d == 0)).ones();
      Q.each_row() %= d.t();

      DR_rot.rows(idx) = Q * DR.rows(idx);
    }

    arma::mat score = DD * DR_rot;
    vals[r] = arma::accu(score % score);
  }

  return vals;
}
