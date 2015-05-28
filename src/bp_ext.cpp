#include <dai/bp.h>

namespace dai {

#define MESSAGE0(j,_I) (_edges[j][_I].message)
#define MESSAGE1(j,_I) (1. - _edges[j][_I].message)
#define OLDPROD0(j)    (_oldProd[j][0])
#define OLDPROD1(j)    (_oldProd[j][1])

void BP::calcIncomingMessageProduct_0011(BP::ProbProduct &prod, size_t I) const
{
    for(const Neighbor &j: nbF(I)) {
        size_t _I = j.dual;
        DAI_ASSERT(_edges[j][_I].index == INDEX_0011);

        // prod._p[0] *= _oldProd[j.node][0] / _edges[j][_I].message._p[0];
        // prod._p[1] *= _oldProd[j.node][0] / _edges[j][_I].message._p[0];
        // prod._p[2] *= _oldProd[j.node][1] / _edges[j][_I].message._p[1];
        // prod._p[3] *= _oldProd[j.node][1] / _edges[j][_I].message._p[1];

        // Pre-calculate reziprocals.
        const ProbProduct::value_type r_m0 = 1. / MESSAGE0(j,_I);
        const ProbProduct::value_type r_m1 = 1. / MESSAGE1(j,_I);

        prod._p[0] *= (OLDPROD0(j.node) * r_m0);
        prod._p[1] *= (OLDPROD0(j.node) * r_m0);
        prod._p[2] *= (OLDPROD1(j.node) * r_m1);
        prod._p[3] *= (OLDPROD1(j.node) * r_m1);
    }
}

void BP::calcIncomingMessageProduct_0101(BP::ProbProduct &prod, size_t I) const
{
    for(const Neighbor &j: nbF(I)) {
        size_t _I = j.dual;
        DAI_ASSERT(_edges[j][_I].index == INDEX_0101);

        // prod._p[0] *= _oldProd[j.node][0] / _edges[j][_I].message._p[0];
        // prod._p[1] *= _oldProd[j.node][1] / _edges[j][_I].message._p[1];
        // prod._p[2] *= _oldProd[j.node][0] / _edges[j][_I].message._p[0];
        // prod._p[3] *= _oldProd[j.node][1] / _edges[j][_I].message._p[1];

        // Pre-calculate reziprocals.
        const ProbProduct::value_type r_m0 = 1. / MESSAGE0(j,_I);
        const ProbProduct::value_type r_m1 = 1. / MESSAGE1(j,_I);

        prod._p[0] *= (OLDPROD0(j.node) * r_m0);
        prod._p[1] *= (OLDPROD1(j.node) * r_m1);
        prod._p[2] *= (OLDPROD0(j.node) * r_m0);
        prod._p[3] *= (OLDPROD1(j.node) * r_m1);
    }
}

void BP::calcIncomingMessageProduct_01(BP::ProbProduct &prod, size_t I) const
{
    for(const Neighbor &j: nbF(I)) {
        size_t _I = j.dual;
        DAI_ASSERT(_edges[j][_I].index == INDEX_01);

        //prod._p[0] *= _oldProd[j.node][0] / _edges[j][_I].message._p[0];
        //prod._p[1] *= _oldProd[j.node][1] / _edges[j][_I].message._p[1];

        // Pre-calculate reziprocals.
        const ProbProduct::value_type r_m0 = 1. / MESSAGE0(j,_I);
        const ProbProduct::value_type r_m1 = 1. / MESSAGE1(j,_I);

        prod._p[0] *= (OLDPROD0(j.node) * r_m0);
        prod._p[1] *= (OLDPROD1(j.node) * r_m1);
    }
}

void BP::calcIncomingMessageProduct_0101_0011(double* prod, size_t I, size_t i) const
{
    const Neighbors& n = nbF(I);
    DAI_ASSERT(nbF(I).size() == 2);

    const size_t _I0 = n[0].dual;
    const size_t n0  = n[0].node;


    if (i!=n0) {
        DAI_ASSERT(_edges[n0][_I0].index == INDEX_0101);

        // prod._p[0] *= _oldProd[n0][0] / _edges[n0][_I0].message._p[0];
        // prod._p[1] *= _oldProd[n0][1] / _edges[n0][_I0].message._p[1];
        // prod._p[2] *= _oldProd[n0][0] / _edges[n0][_I0].message._p[0];
        // prod._p[3] *= _oldProd[n0][1] / _edges[n0][_I0].message._p[1];

        prod[0] *= OLDPROD0(n0) * _edges[n0][_I0].reciprocals[0];
        prod[1] *= OLDPROD1(n0) * _edges[n0][_I0].reciprocals[1];
        prod[2] *= OLDPROD0(n0) * _edges[n0][_I0].reciprocals[0];
        prod[3] *= OLDPROD1(n0) * _edges[n0][_I0].reciprocals[1];
    }

    const size_t _I1 = n[1].dual;
    const size_t n1  = n[1].node;

    if (i!=n1) {
        DAI_ASSERT(_edges[n1][_I1].index == INDEX_0011);

        // prod._p[0] *= _oldProd[n1][0] / _edges[n1][_I1].message._p[0];
        // prod._p[1] *= _oldProd[n1][0] / _edges[n1][_I1].message._p[0];
        // prod._p[2] *= _oldProd[n1][1] / _edges[n1][_I1].message._p[1];
        // prod._p[3] *= _oldProd[n1][1] / _edges[n1][_I1].message._p[1];

        prod[0] *= OLDPROD0(n1) * _edges[n1][_I1].reciprocals[0];
        prod[1] *= OLDPROD0(n1) * _edges[n1][_I1].reciprocals[0];
        prod[2] *= OLDPROD1(n1) * _edges[n1][_I1].reciprocals[1];
        prod[3] *= OLDPROD1(n1) * _edges[n1][_I1].reciprocals[1];
    }
}


} // namespace dai
