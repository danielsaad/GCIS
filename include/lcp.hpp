#include "util.hpp"

template <class text_t, class sa_t, class lcp_t>
int lcp_PHI(const text_t &T, const sa_t &SA, lcp_t &LCP, uint_t n) {

    int_t *PLCP = new int_t[n];
    std::fill_n(PLCP, n, 0);

    // PHI is stored in PLCP array
    int_t i, j;
    for (i = 0, j = 0; i < n; ++i) {
        PLCP[SA[i]] = j;
        j = SA[i];
    }

    int_t l;
    for (i = 0, l = 0; i < n - 1; ++i) {
        uint_t ii = PLCP[i];

        while (T[i + l] == T[ii + l]) {
            ++l;
        }
        PLCP[i] = l;
        if (l)
            --l;
    }

    for (i = 1; i < n; ++i) {
        LCP[i] = PLCP[SA[i]];
    }

    delete[] PLCP;

    return 0;
}