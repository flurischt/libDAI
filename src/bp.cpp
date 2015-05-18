/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <dai/dai_config.h>
#ifdef DAI_WITH_BP


#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <dai/bp.h>
#include <dai/util.h>
#include <dai/properties.h>

//#define DAI_VERBOSE
#ifdef DAI_VERBOSE
#   define DAI_LOG(MESSAGE) do { std::cout << MESSAGE << std::endl; } while(0)
#else
#   define DAI_LOG(MESSAGE)
#endif

namespace dai {

using namespace std;

void BP::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("tol") );

    props.tol = opts.getStringAs<Real>("tol");

    if( opts.hasKey("maxiter") )
        props.maxiter = opts.getStringAs<size_t>("maxiter");
    else
        props.maxiter = 10000;
    if( opts.hasKey("maxtime") )
        props.maxtime = opts.getStringAs<Real>("maxtime");
    else
        props.maxtime = INFINITY;
    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<Real>("damping");
    else
        props.damping = 0.0;
}


PropertySet BP::getProperties() const {
    PropertySet opts;
    opts.set( "tol", props.tol );
    opts.set( "maxiter", props.maxiter );
    opts.set( "maxtime", props.maxtime );
    opts.set( "verbose", props.verbose );
    opts.set( "damping", props.damping );
    return opts;
}


string BP::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "maxtime=" << props.maxtime << ",";
    s << "verbose=" << props.verbose << ",";
    s << "damping=" << props.damping << "]";
    return s.str();
}


void BP::construct() {
    // create edge properties
    _edges.clear();
    _edges.reserve( nrVars() );
    _indices.clear();
    _oldProd.clear();
    _edge2lutNew.clear();
    _edge2lutNew.reserve( nrVars() );

    // Add the predefined indices first.
    // The order matters and will be used later!
    _indices.push_back({0,0,1,1}); DAI_DEBASSERT(_indices.size()-1 == INDEX_0011);
    _indices.push_back({0,1,0,1}); DAI_DEBASSERT(_indices.size()-1 == INDEX_0101);
    _indices.push_back({0,1});     DAI_DEBASSERT(_indices.size()-1 == INDEX_01);

    for( size_t i = 0; i < nrVars(); ++i ) {
        _edges.push_back( vector<EdgeProp>() );
        _edges[i].reserve( nbV(i).size() );
        _oldProd.push_back(vector<double>(var(i).states(), 1));
        _edge2lutNew.push_back( vector<heap_data_handle>() );
        _edge2lutNew[i].reserve( nbV(i).size() );
        for( const Neighbor &I : nbV(i) ) {
            EdgeProp newEP;
#ifndef DAI_RECOMMENDER_BOOST
            newEP.message = Prob( var(i).states() );
            newEP.newMessage = Prob( var(i).states() );
#else
            newEP.message = 0;
            newEP.newMessage = 0;
#endif

            ind_t index;
            for( IndexFor k( var(i), factor(I).vars() ); k.valid(); ++k )
                index.push_back( k );
            auto it = std::find(_indices.begin(), _indices.end(), index);
            if ( it == _indices.end() ) {
                _indices.push_back(index);
                newEP.index = _indices.size() - 1;
                //cout << "Added index: " << index << endl;
            }
            else {
                newEP.index = it - _indices.begin();
                DAI_DEBASSERT(newEP.index < _indices.size());
            }
            DAI_DEBASSERT(_indices[newEP.index] == index);

            newEP.residual = 0.0;
            _edges[i].push_back( newEP );
            _edge2lutNew[i].push_back( _lutNew.push( make_pair( newEP.residual, make_pair( i, _edges[i].size() - 1 ))));
        }
    }

    // create old beliefs
    _oldBeliefsV.clear();
    _oldBeliefsV.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i )
        _oldBeliefsV.push_back( Factor( var(i) ) );
    _oldBeliefsF.clear();
    _oldBeliefsF.reserve( nrFactors() );
    for( size_t I = 0; I < nrFactors(); ++I )
        _oldBeliefsF.push_back( Factor( factor(I).vars() ) );

    // create update sequence
    _updateSeq.clear();
    _updateSeq.reserve( nrEdges() );
    for( size_t I = 0; I < nrFactors(); I++ )
        for( const Neighbor &i : nbF(I) )
            _updateSeq.push_back( Edge( i, i.dual ) );
}


void BP::init() {

    // The initialization of the messages has a strong impact on the
    // convergence. The value originally suggested by libDAI was c = 1.0.
    // I consider this choice unnatural because it will lead to messages
    // that don't sum up to 1. This seems wrong, because throughout
    // the algorithm, messages are "normalized".
    // Changing c to (e.g.) 0.5 leads to faster convergence: 40% less
    // messages are required for example for the u1 dataset!
    // For DAI_RECOMMENDER_BOOST, it is required to choose a value
    // different from 1.0 - otherwise we run into problems with
    // division by zero (as p[1] = 1-p[0] = 0)...
    // The most appropriate choice seems to me the uniform distribution:
    // c = (Real)1./message.size();
    // Addendum: the observed faster convergence was observed only for
    // single precision. Double precision behaves identical (in terms
    // of number of messages processed).

    //Real c = 1.0;
    Real c = 0.5;
    for( size_t i = 0; i < nrVars(); ++i ) {
        for( const Neighbor &I : nbV(i) ) {
#ifndef DAI_RECOMMENDER_BOOST
            message( i, I.iter ).fill( c );
            newMessage( i, I.iter ).fill( c );
            updateResidual( i, I.iter, 0.0 );
#else
            message( i, I.iter ) = c;
            newMessage( i, I.iter ) = c;
            updateResidual( i, I.iter, 0.0 );
#endif
        }
    }
    _iters = 0;
    messageCount = 0;
}


bool BP::findMaxResidual( size_t &i, size_t &_I ) {
    DAI_ASSERT( !_lutNew.empty() );
    i  = _lutNew.top().second.first;
    _I = _lutNew.top().second.second;

#if 0
    static int count = 0; count ++;
    static Real sum = 0.f; sum += _lutNew.top().first;
    if (isnan(_lutNew.top().first))
        cout << "Warning: invalid residual occured for " << i << " <-- " << _I << endl;
    if (count % 1000 == 0) {
        cout << "Moving avgerage of residuals: "<< sum / 1000 << " " << _lutNew.top().first << endl;
        sum = 0;
    }
#endif

    return _lutNew.top().first > 0;
}


// TODO: Optimize (in progress)
void BP::calcIncomingMessageProduct(ProbProduct &prod, size_t I, bool without_i, size_t i) const {

    // Calculate product of incoming messages and factor I
    for(const Neighbor &j: nbF(I)) {
        if( !(without_i && (j == i)) ) {

            // TODO: the calculation in this loop got very cryptic.
            // One might want to have some explanations at some point...

            // The message that should not go into the product is the one from that
            // node that that message will be sent to. Conveniently, the value is
            // already available: j.dual.
            size_t _I = j.dual;

            // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
            const ind_t &ind = index(j, _I);
#ifdef DAI_RECOMMENDER_BOOST
            DAI_DEBASSERT(var(j).states() == 2);
            Real message_0 = _edges[j][_I].message;
            Real message_1 = (Real)1 - message_0;
#endif

            for(size_t r = 0; r < prod.size(); ++r) {

                // Let's divide by that message that should not go into the product.
                // Calculate with double precision!
#ifndef DAI_RECOMMENDER_BOOST
                double prod_jk = _oldProd[j.node][ind[r]] / _edges[j][_I].message._p[ind[r]];
#else
                double prod_jk = (ind[r] == 0)
                        ? _oldProd[j.node][0] / message_0
                        : _oldProd[j.node][1] / message_1;
#endif

                // And multiply it with the target.
                prod._p[r] *= prod_jk;

                // This is a hack to handle issues with precision.
                //if (normalize) prod.normalize();
            }
        }
    }
}

#ifndef DAI_RECOMMENDER_BOOST
void BP::marginalizeProductOntoMessage(const ProbProduct &prod, size_t i, size_t _I)
{
#ifdef DAI_SINGLE_PRECISION
        MessageType &marg = newMessage(i,_I);
        if (_marg.size() != marg.size())
            _marg.resize(marg.size());

        // Shortcut, to avoid code duplication. We are making use of the
        // fact that m is always of type ProbD == ProbProduct!
        // m <-- _marg = ProbProduct(newMessage(i, _I))
        ProbProduct &m = _marg;
#else
        MessageType &marg = newMessage(i,_I);

        // Shortcut, to avoid code duplication. We are making use of the
        // fact that m is always of type ProbD == ProbProduct!
        // m <-- marg = newMessage(i, _I)
        ProbProduct &m = marg;
#endif

        // Calculate marginal AND normalize probability.
        // Avoid the indirect lookup via ind_t if possible.
        switch (_edges[i][_I].index) {
            case INDEX_0011: {
                const ProbProduct::value_type a = (prod._p[0]+prod._p[1]);
                const ProbProduct::value_type s = a + (prod._p[2]+prod._p[3]);
                m._p[0] = a/s;
                m._p[1] = 1. - m._p[0];
            } break;
            case INDEX_0101: {
                const ProbProduct::value_type a = (prod._p[0]+prod._p[2]);
                const ProbProduct::value_type s = a + (prod._p[1]+prod._p[3]);
                m._p[0] = a/s;
                m._p[1] = 1. - m._p[0];
            } break;
            default: {
                std::fill(m._p.begin(), m._p.end(), 0.0);
                // ind is the precalculated IndexFor(i,I) i.e. to x_I == k
                // corresponds x_i == ind[k]
                const ind_t& ind = index(i,_I);
                for( size_t r = 0; r < prod.size(); ++r )
                    m._p[ind[r]] += prod[r];
                m.normalizeFast();
            }
        }

#ifdef DAI_SINGLE_PRECISION
        // Copy (and cast) from m to newMessage(i,_I);
        std::copy(m._p.begin(), m._p.end(), marg._p.begin());
#endif
}
#else
void BP::marginalizeProductOntoMessage(const ProbProduct &prod, size_t i, size_t _I)
{
    MessageType &marg = newMessage(i,_I);

    // Calculate marginal AND normalize probability.
    // Avoid the indirect lookup via ind_t if possible.
    switch (_edges[i][_I].index) {
        case INDEX_0011: {
            const ProbProduct::value_type a = (prod._p[0]+prod._p[1]);
            const ProbProduct::value_type s = a + (prod._p[2]+prod._p[3]);
            marg = a/s;
        } break;
        case INDEX_0101: {
            const ProbProduct::value_type a = (prod._p[0]+prod._p[2]);
            const ProbProduct::value_type s = a + (prod._p[1]+prod._p[3]);
            marg = a/s;
        } break;
        default: {
            marg = 0;
            // ind is the precalculated IndexFor(i,I) i.e. to x_I == k
            // corresponds x_i == ind[k]
            const ind_t& ind = index(i,_I);
            ProbProduct::value_type a = 0.;
            ProbProduct::value_type s = 0.;
            for( size_t r = 0; r < prod.size(); ++r )
            {
                if (ind[r] == 0)
                    a += prod[r];
                s += prod[r];
                DAI_ASSERT(ind[r] == 0 || ind[r] == 1);
            }
            marg = a/s;
        }
    }
}
#endif


void BP::calcNewMessage( size_t i, size_t _I) {

    // load
    size_t I = _G.nb1(i)[_I].node;

    // The following applies only rarely  (uV2New1: 50x, uNew1 and u1: 135x)
    // TODO: investigate further, can this still be useful?
    // UPDATE: image segmentation example doesn't converge if this "optimization"
    // is removed. I don't fully get it though. NJU
    if( _factors[I].vars().size() == 1 ) {    // optimization
#ifndef DAI_RECOMMENDER_BOOST
        std::copy(_factors[I].p().begin(),
                  _factors[I].p().end(),
                  newMessage(i,_I)._p.begin());
#else
        newMessage(i,_I) = _factors[I].p()[0];
#endif
    }
    else {
        // calculate updated message I->i

        // The capacity of _prod is not changed here. malloc/free will be called
        // very rarely. However, this can be further improved, because
        //  _factors[I].p().size() cleanly toggles between 2 and 4:
        // TODO: create two containers _prod4 and _prod2 and cleverly call
        // calcNewMessage() with either one as argument.
        if (_prod.size() != _factors[I].p().size())
            _prod.resize(_factors[I].p().size());
        std::copy(_factors[I].p().begin(), _factors[I].p().end(), _prod.begin());

        // Calc the message product.
        DAI_LOG("calcNewMessage " << I << " <-> " << i);
#ifdef DAI_RECOMMENDER_BOOST
        calcIncomingMessageProduct_0101_0011(_prod, I, i);
#else
        calcIncomingMessageProduct(_prod, I, true, i);
#endif

        // Marginalize onto i
        marginalizeProductOntoMessage(_prod, i, _I);
    }

    // Update the residual if necessary
#ifdef DAI_RECOMMENDER_BOOST
    // Make use of the fact that message.size() == 2 and that
    // the messages are normalized to 1.
    Real r = std::abs( newMessage( i, _I ) - message( i, _I ) );
#else
    Real r = distFast( newMessage( i, _I ), message( i, _I ) );
#endif
    updateResidual( i, _I , r );
}


// BP::run does not check for NANs for performance reasons
// Somehow NaNs do not often occur in BP...
Real BP::run() {
    if( props.verbose >= 1 )
        cerr << "Starting " << identify() << "...";
    if( props.verbose >= 3)
        cerr << endl;

    double tic = toc();

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    Real maxDiff = INFINITY;
    bool hasValidResidual = true;
    for( ; (_iters < props.maxiter)
           && (maxDiff > props.tol)
           && ((toc() - tic) < props.maxtime)
           && hasValidResidual; _iters++ ) {
        if( _iters == 0 ) {
            // do the first pass
            for( size_t i = 0; i < nrVars(); ++i )
                for( const Neighbor &I : nbV(i) ) {
                    calcNewMessage( i, I.iter);
                }
        }
        // Maximum-Residual BP [\ref EMK06]
        for( size_t t = 0; t < _updateSeq.size(); ++t ) {
            // update the message with the largest residual
            size_t i, _I;
            hasValidResidual = findMaxResidual( i, _I );
            if (!hasValidResidual) break;

            DAI_LOG("updating message from " << i << " to " << _I);
            updateMessage( i, _I );

            // I->i has been updated, which means that residuals for all
            // J->j with J in nb[i]\I and j in nb[J]\i have to be updated
            for( const Neighbor &J: nbV(i) ) {
                if( J.iter != _I ) {
                    for( const Neighbor &j: nbF(J) ) {
                        size_t _J = j.dual;
                        if( j != i )
                            calcNewMessage( j, _J);
                    }
                }
            }
        }

        // calculate new beliefs and compare with old ones
        maxDiff = -INFINITY;
        for( size_t i = 0; i < nrVars(); ++i ) {
            Factor b( beliefV(i) );
            maxDiff = std::max( maxDiff, distFast( b.p(), _oldBeliefsV[i].p() ) );
            _oldBeliefsV[i] = b;
        }
        for( size_t I = 0; I < nrFactors(); ++I ) {
            Factor b( beliefF(I) );
            maxDiff = std::max( maxDiff, distFast( b.p(), _oldBeliefsF[I].p() ) );
            _oldBeliefsF[I] = b;
        }

        if( props.verbose >= 3 )
            cerr << name() << "::run:  maxdiff " << maxDiff << " after " << _iters+1 << " passes" << endl;
    }

    if( maxDiff > _maxdiff )
        _maxdiff = maxDiff;

    if( props.verbose >= 1 ) {
        if( maxDiff > props.tol ) {
            if( props.verbose == 1 )
                cerr << endl;
                cerr << name() << "::run:  WARNING: not converged after " << _iters << " passes (" << toc() - tic << " seconds)...final maxdiff:" << maxDiff << endl;
        } else {
            if( props.verbose >= 3 )
                cerr << name() << "::run:  ";
                cerr << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

#if 0
    // Print how the messages look like at the end of the calculation.
    for( size_t i = 0; i < nrVars(); ++i )
        for ( const Neighbor &I : nbV(i))
            cout << message(i,I.iter) << " ";
#endif

    return maxDiff;
}


void BP::calcBeliefV( size_t i, ProbProduct &p ) const {
    p.resize(var(i).states());
    std::fill(p._p.begin(), p._p.end(), 1.0);
    for ( const Neighbor &I : nbV(i) )
    {
#ifndef DAI_RECOMMENDER_BOOST
        for (size_t j=0; j<p.size(); ++j)
            p._p[j] *= newMessage( i, I.iter )._p[j];
#else
        p._p[0] *= newMessage( i, I.iter );
        p._p[1] *= ((Real)1-newMessage( i, I.iter ));
#endif
    }
}


Factor BP::beliefV( size_t i ) const {
    calcBeliefV( i, _probTemp );
    _probTemp.normalize();

    // Factor is created each time. Could be avoided...
    // Currently not a bottleneck, so no need to change InfAlg interface.
    return( Factor( var(i), Prob(_probTemp.begin(), _probTemp.end(), _probTemp.size()) ) );
}


Factor BP::beliefF( size_t I ) const {
    Factor fac( factor(I) );
    Prob &p = fac.p();
#ifdef DAI_SINGLE_PRECISION
    ProbProduct pd(p.begin(), p.end(), p.size());
    calcBeliefF( I, pd );
    pd.normalize();
    std::copy(pd.begin(), pd.end(), p.begin());
#else
    calcBeliefF( I, p );
    p.normalize();
#endif

    return fac;
}


vector<Factor> BP::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); ++i )
        result.push_back( beliefV(i) );
    for( size_t I = 0; I < nrFactors(); ++I )
        result.push_back( beliefF(I) );
    return result;
}


Factor BP::belief( const VarSet &ns ) const {
    if( ns.size() == 0 )
        return Factor();
    else if( ns.size() == 1 )
        return beliefV( findVar( *(ns.begin() ) ) );
    else {
        size_t I;
        for( I = 0; I < nrFactors(); I++ )
            if( factor(I).vars() >> ns )
                break;
        if( I == nrFactors() )
            DAI_THROW(BELIEF_NOT_AVAILABLE);
        return beliefF(I).marginal(ns);
    }
}


Real BP::logZ() const {
    Real sum = 0.0;
    for( size_t i = 0; i < nrVars(); ++i )
        sum += (1.0 - nbV(i).size()) * beliefV(i).entropy();
    for( size_t I = 0; I < nrFactors(); ++I )
        sum -= dist( beliefF(I), factor(I), DISTKL );
    return sum;
}


void BP::init( const VarSet &ns ) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); ++n ) {
        size_t ni = findVar( *n );
        for( const Neighbor &I : nbV( ni ) ) {
            Real val = 1.0;
#ifndef DAI_RECOMMENDER_BOOST
            message( ni, I.iter ).fill( val );
            newMessage( ni, I.iter ).fill( val );
            updateResidual( ni, I.iter, 0.0 );
#else
            message( ni, I.iter ) = val;
            newMessage( ni, I.iter ) = val;
            updateResidual( ni, I.iter, 0.0 );
#endif
        }
    }
    _iters = 0;
    messageCount = 0;
}

#ifndef DAI_RECOMMENDER_BOOST
void BP::updateMessage( size_t i, size_t _I ) {
    for (size_t j=0; j<_oldProd[i].size(); ++j) {
        _oldProd[i][j] =  _oldProd[i][j] / _edges[i][_I].message._p[j] * _edges[i][_I].newMessage._p[j];
    }

    // Count message.
    messageCount++;
    if( recordSentMessages )
        _sentMessages.push_back(make_pair(i,_I));
    if( props.damping == 0.0 ) {
        message(i,_I) = newMessage(i,_I);
        updateResidual( i, _I, 0.0 );
    } else {
        message(i,_I) = (message(i,_I) ^ props.damping) * (newMessage(i,_I) ^ (1.0 - props.damping));
        updateResidual( i, _I, distFast( newMessage(i,_I), message(i,_I) ) );
    }
}
#else
void BP::updateMessage( size_t i, size_t _I ) {

    // Damping is not supported here.
    DAI_DEBASSERT(props.damping == false);
    DAI_DEBASSERT(_oldProd[i].size() == 2);

    _oldProd[i][0] =  _oldProd[i][0] /       _edges[i][_I].message  *       _edges[i][_I].newMessage;
    _oldProd[i][1] =  _oldProd[i][1] / (1. - _edges[i][_I].message) * (1. - _edges[i][_I].newMessage);

    // Count message.
    messageCount++;
    if( recordSentMessages )
        _sentMessages.push_back(make_pair(i,_I));
    message(i,_I) = newMessage(i,_I);
    updateResidual( i, _I, 0.0 );
}
#endif

// TODO: Optimize: We are using a heap now but this is not faster then the
// multimap solution. So we might have to revert to it.
void BP::updateResidual( size_t i, size_t _I, Real r ) {
    EdgeProp* pEdge = &_edges[i][_I];
    pEdge->residual = r;

    // rearrange look-up table (delete and reinsert new key)
    _lutNew.update(_edge2lutNew[i][_I], make_pair( r, make_pair(i, _I) ));
}


} // end of namespace dai


#endif
