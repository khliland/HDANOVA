#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

/**
 * Fisher-Yates shuffle using R's internal RNG for consistency with set.seed()
 */
inline void shuffle_uvec_inplace(arma::uvec& x) {
    const arma::uword n = x.n_elem;
    if (n <= 1) return;
    
    for (arma::uword i = n - 1; i > 0; --i) {
        // Use R::unif_rand() and ensure proper scaling to [0, i]
        // static_cast to double for precise multiplication before floor
        arma::uword j = static_cast<arma::uword>(std::floor(R::unif_rand() * static_cast<double>(i + 1)));
        if (j > i) j = i; // Bounds safety
        std::swap(x(i), x(j));
    }
}

// [[Rcpp::export(name = ".permutation_ssq_kernel_cpp")]]
Rcpp::NumericVector permutation_ssq_kernel_cpp(const arma::mat& DD,
                                               const arma::mat& DR,
                                               const Rcpp::List& groups,
                                               const int n_perm) {
    Rcpp::RNGScope scope;

    const arma::uword n_groups = static_cast<arma::uword>(groups.size());
    std::vector<arma::uvec> idx_store;
    idx_store.reserve(n_groups);

    for (arma::uword g = 0; g < n_groups; ++g) {
        Rcpp::IntegerVector block = groups[g];
        if (block.size() <= 1) continue;
        
        arma::uvec idx(block.size());
        for (int i = 0; i < block.size(); ++i) {
            // Convert R 1-based index to C++ 0-based index
            idx(i) = static_cast<arma::uword>(block[i] - 1);
        }
        idx_store.push_back(idx);
    }

    const arma::uword N = DR.n_rows;
    const arma::uword K = DR.n_cols;
    
    // base_order stays constant, order is modified per iteration
    const arma::uvec base_order = arma::regspace<arma::uvec>(0, N - 1);
    arma::uvec order = base_order;
    arma::uvec perm_idx;
    
    // Pre-allocate the score matrix to avoid allocation inside the loop
    // Dimension: (Rows of DD) x (Cols of DR)
    arma::mat score(DD.n_rows, K);
    Rcpp::NumericVector vals(n_perm);

    for (int p = 0; p < n_perm; ++p) {
        // Reset order to baseline
        order = base_order;

        for (std::size_t g = 0; g < idx_store.size(); ++g) {
            perm_idx = idx_store[g];
            shuffle_uvec_inplace(perm_idx);
            // Re-map the block indices
            order.elem(idx_store[g]) = perm_idx;
        }

        // Calculation: score = DD * DR[order, ]
        // Using square() + accu() is generally faster than % (element-wise mult)
        score = DD * DR.rows(order);
        vals[p] = arma::accu(arma::square(score));
    }

    return vals;
}

// [[Rcpp::export(name = ".rotation_ssq_kernel_cpp")]]
Rcpp::NumericVector rotation_ssq_kernel_cpp(const arma::mat& DD,
                                            const arma::mat& DR,
                                            const Rcpp::List& groups,
                                            const int n_rot) {
    Rcpp::RNGScope scope;

    const arma::uword n_groups = static_cast<arma::uword>(groups.size());
    std::vector<arma::uvec> idx_store;
    idx_store.reserve(n_groups);

    for (arma::uword g = 0; g < n_groups; ++g) {
        Rcpp::IntegerVector block = groups[g];
        if (block.size() <= 1) continue;
        
        arma::uvec idx(block.size());
        for (int i = 0; i < block.size(); ++i) {
            idx(i) = static_cast<arma::uword>(block[i] - 1);
        }
        idx_store.push_back(idx);
    }

    Rcpp::NumericVector vals(n_rot);
    arma::mat DR_rot = DR;
    arma::mat Q, R, Z, score;

    for (int r = 0; r < n_rot; ++r) {
        DR_rot = DR;

        for (std::size_t g = 0; g < idx_store.size(); ++g) {
            const arma::uvec& idx = idx_store[g];
            const arma::uword m = idx.n_elem;

            Z = arma::randn<arma::mat>(m, m);
            arma::qr_econ(Q, R, Z);

            arma::vec d = arma::sign(R.diag());
            d.transform( [](double val) { return (val == 0.0) ? 1.0 : val; } );
            
            Q.each_row() %= d.t();

            DR_rot.rows(idx) = Q * DR.rows(idx);
        }

        score = DD * DR_rot;
        vals[r] = arma::accu(arma::square(score));
    }

    return vals;
}