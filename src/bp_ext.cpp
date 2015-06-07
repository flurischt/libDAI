#include <dai/bp.h>

namespace dai {


#ifdef DAI_VECTORIZATION
void BP::calcIncomingMessageProduct_0101_0011(__m256d& prod_vec, size_t I, size_t i) const {
        const Neighbors& n = nbF(I);
        DAI_ASSERT(nbF(I).size() == 2);
        const size_t n0  = n[0].node;
        if (i!=n0) {
            const size_t _I0 = n[0].dual;
            DAI_ASSERT(_edges[n0][_I0].index == INDEX_0101);

            __m256d temp_vec = _mm256_mul_pd(_oldProd[n0], _mm256_cvtps_pd(_edges[n0][_I0].reciprocals));
            // a a b b into a b a b
            //std::cout << "tempvec1" << ((double*)&temp_vec)[0] << "   "<< ((double*)&temp_vec)[1] << "   "<< ((double*)&temp_vec)[2] << "   "<< ((double*)&temp_vec)[3] << std::endl;
            __m256d xswap = _mm256_permute2f128_pd(temp_vec, temp_vec, 0x01);
            temp_vec = _mm256_blend_pd(xswap, temp_vec, 0b1001);
            //std::cout << "tempvec2" << ((double*)&temp_vec)[0] << "   "<< ((double*)&temp_vec)[1] << "   "<< ((double*)&temp_vec)[2] << "   "<< ((double*)&temp_vec)[3] << std::endl;
            prod_vec = _mm256_mul_pd(prod_vec, temp_vec);
        }

        const size_t n1  = n[1].node;
        if (i!=n1) {
            const size_t _I1 = n[1].dual;
            DAI_ASSERT(_edges[n1][_I1].index == INDEX_0011);
            __m256d temp_vec = _mm256_mul_pd(_oldProd[n1], _mm256_cvtps_pd(_edges[n1][_I1].reciprocals));
            prod_vec = _mm256_mul_pd(prod_vec, temp_vec);
        }
    }
#else
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

        double temp1 = _oldProd[n0][0] * _edges[n0][_I0].reciprocals[0];
        double temp2 = _oldProd[n0][1] * _edges[n0][_I0].reciprocals[1];
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

        double temp1 = _oldProd[n1][0] * _edges[n1][_I1].reciprocals[0];
        double temp2 = _oldProd[n1][1] * _edges[n1][_I1].reciprocals[1];
        prod[0] *= temp1;
        prod[1] *= temp1;
        prod[2] *= temp2;
        prod[3] *= temp2;
    }
}
#endif

} // namespace dai