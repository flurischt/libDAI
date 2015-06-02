#include <dai/bp.h>

namespace dai {

#define OLDPROD0(j)    (_oldProd[j][0])
#define OLDPROD1(j)    (_oldProd[j][1])

void BP::calcIncomingMessageProduct_0101_0011(__m256d& prod_vec, size_t I, size_t i) const {
        const Neighbors& n = nbF(I);
        DAI_ASSERT(nbF(I).size() == 2);


        const size_t n0  = n[0].node;
        if (i!=n0) {
            const size_t _I0 = n[0].dual;
            DAI_ASSERT(_edges[n0][_I0].index == INDEX_0101);

            double temp1 = OLDPROD0(n0) * _edges[n0][_I0].reciprocals[0];
            double temp2 = OLDPROD1(n0) * _edges[n0][_I0].reciprocals[1];
            __m256d temp_vec = _mm256_set_pd(temp2, temp1, temp2, temp1);
            prod_vec = _mm256_mul_pd(prod_vec, temp_vec);
        }

        const size_t n1  = n[1].node;
        if (i!=n1) {
            const size_t _I1 = n[1].dual;
            DAI_ASSERT(_edges[n1][_I1].index == INDEX_0011);

            double temp1 = OLDPROD0(n1) * _edges[n1][_I1].reciprocals[0];
            double temp2 = OLDPROD1(n1) * _edges[n1][_I1].reciprocals[1];
//            __m128d hi = _mm_set1_pd(temp2);
//            __m128d lo = _mm_set1_pd(temp1);
//            __m256d temp_vec = _mm256_setr_m128d(hi, lo);
            __m256d temp_vec = _mm256_set_pd(temp2, temp2, temp1, temp1);
            prod_vec = _mm256_mul_pd(prod_vec, temp_vec);
        }
    }



void BP::calcIncomingMessageProduct_0101_0011(double* prod, size_t I, size_t i) const {
    const Neighbors& n = nbF(I);
    DAI_ASSERT(nbF(I).size() == 2);

    const size_t n0  = n[0].node;
    if (i!=n0) {
        const size_t _I0 = n[0].dual;
        DAI_ASSERT(_edges[n0][_I0].index == INDEX_0101);

        // prod._p[0] *= _oldProd[n0][0] / _edges[n0][_I0].message._p[0];
        // prod._p[1] *= _oldProd[n0][1] / _edges[n0][_I0].message._p[1];
        // prod._p[2] *= _oldProd[n0][0] / _edges[n0][_I0].message._p[0];
        // prod._p[3] *= _oldProd[n0][1] / _edges[n0][_I0].message._p[1];

        double temp1 = OLDPROD0(n0) * _edges[n0][_I0].reciprocals[0];
        double temp2 = OLDPROD1(n0) * _edges[n0][_I0].reciprocals[1];
        prod[0] *= temp1;
        prod[1] *= temp2;
        prod[2] *= temp1;
        prod[3] *= temp2;
    }

    const size_t n1  = n[1].node;
    if (i!=n1) {
        const size_t _I1 = n[1].dual;
        DAI_ASSERT(_edges[n1][_I1].index == INDEX_0011);

        // prod._p[0] *= _oldProd[n1][0] / _edges[n1][_I1].message._p[0];
        // prod._p[1] *= _oldProd[n1][0] / _edges[n1][_I1].message._p[0];
        // prod._p[2] *= _oldProd[n1][1] / _edges[n1][_I1].message._p[1];
        // prod._p[3] *= _oldProd[n1][1] / _edges[n1][_I1].message._p[1];

        double temp1 = OLDPROD0(n1) * _edges[n1][_I1].reciprocals[0];
        double temp2 = OLDPROD1(n1) * _edges[n1][_I1].reciprocals[1];
        prod[0] *= temp1;
        prod[1] *= temp1;
        prod[2] *= temp2;
        prod[3] *= temp2;
    }
}
} // namespace dai