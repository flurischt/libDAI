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
    DAI_ASSERT( opts.hasKey("updates") );

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
    _oldProd.clear();
    _edge2lutNew.clear();
    _edge2lutNew.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i ) {
        _edges.push_back( vector<EdgeProp>() );
        _edges[i].reserve( nbV(i).size() );
        _oldProd.push_back(vector<double>(var(i).states(), 1));
        _edge2lutNew.push_back( vector<heap_data_handle>() );
        _edge2lutNew[i].reserve( nbV(i).size() );
        for( const Neighbor &I : nbV(i) ) {
            EdgeProp newEP;
            newEP.message = Prob( var(i).states() );
            newEP.newMessage = Prob( var(i).states() );

            newEP.index.reserve( factor(I).nrStates() );
            for( IndexFor k( var(i), factor(I).vars() ); k.valid(); ++k )
                newEP.index.push_back( k );

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
    Real c = 1.0;
    for( size_t i = 0; i < nrVars(); ++i ) {
        for( const Neighbor &I : nbV(i) ) {
            message( i, I.iter ).fill( c );
            newMessage( i, I.iter ).fill( c );
            updateResidual( i, I.iter, 0.0 );
        }
    }
    _iters = 0;
}


void BP::findMaxResidual( size_t &i, size_t &_I ) {
    DAI_ASSERT( !_lutNew.empty() );
    i  = _lutNew.top().second.first;
    _I = _lutNew.top().second.second;
}


// TODO: Optimize (in progress)
void BP::calcIncomingMessageProduct(Prob &prod, size_t I, bool without_i, size_t i) const {
    // Calculate product of incoming messages and factor I
    for(const Neighbor &j: nbF(I)) {
        if( !(without_i && (j == i)) ) {
            // prod_j will be the product of messages coming into j
            vector<double> prod_j;
            prod_j.reserve(_oldProd[j.node].size());
            // We need to find out which message we should not take into the product.
            // TODO: In the future we can try to pass this value as a parameter or compute it in a faster way.
            size_t Iiter = -1;
            for(const Neighbor &J: nbV(j)) {
                if( J == I ) {
                    Iiter = J.iter;
                    break;
                }
            }
            // Now let us divide by that message.
            for (size_t k=0; k<_oldProd[j.node].size(); ++k) {
                prod_j.push_back(_oldProd[j.node][k] / _edges[j][Iiter].message._p[k]);
            }

            DAI_LOG("Product of incoming messages into " << j << " is " << prod_j);

            // TODO: If we understand this we might be able to get rid of this whole function call and use _oldProd directly.
            // multiply prod with prod_j
            size_t _I = j.dual;
            // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
            const ind_t &ind = index(j, _I);
            for(size_t r = 0; r < prod.size(); ++r) {
                prod._p[r] *= prod_j[ind[r]];
            }
        }
    }
}


void BP::calcNewMessage( size_t i, size_t _I) {

    // load
    size_t I = _G.nb1(i)[_I].node;
    Prob marg;

    // The following applies only rarely  (uV2New1: 50x, uNew1 and u1: 135x)
    // TODO: investigate further, can this still be useful?
#if 0
    if( factor(I).vars().size() == 1 ) // optimization
        marg = factor(I).p();
    else {
#endif

    // calculate updated message I->i
    Prob prod = Prob(_factors[I].p());
    calcIncomingMessageProduct(prod, I, true, i);
    DAI_LOG("calcNewMessage " << I << " <-> " << i);

    // Marginalize onto i
    marg = Prob( var(i).states(), 0.0 );
    // ind is the precalculated IndexFor(i,I) i.e. to x_I == k corresponds x_i == ind[k]
    const ind_t ind = index(i,_I);
    for( size_t r = 0; r < prod.size(); ++r )
        marg.set( ind[r], marg[ind[r]] + prod[r] );
    marg.normalize();

    // Store result
    newMessage(i,_I) = marg;

    // Update the residual if necessary
    updateResidual( i, _I , dist( newMessage( i, _I ), message( i, _I ), DISTLINF ) );
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
    for( ; _iters < props.maxiter && maxDiff > props.tol && (toc() - tic) < props.maxtime; _iters++ ) {
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
            findMaxResidual( i, _I );
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
            maxDiff = std::max( maxDiff, dist( b, _oldBeliefsV[i], DISTLINF ) );
            _oldBeliefsV[i] = b;
        }
        for( size_t I = 0; I < nrFactors(); ++I ) {
            Factor b( beliefF(I) );
            maxDiff = std::max( maxDiff, dist( b, _oldBeliefsF[I], DISTLINF ) );
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

    return maxDiff;
}


void BP::calcBeliefV( size_t i, Prob &p ) const {
    p = Prob( var(i).states(), 1.0);
    for ( const Neighbor &I : nbV(i) )
            p *= newMessage( i, I.iter );
}


Factor BP::beliefV( size_t i ) const {
    Prob p;
    calcBeliefV( i, p );
    p.normalize();

    return( Factor( var(i), p ) );
}


Factor BP::beliefF( size_t I ) const {
    Factor Fprod( factor(I) );
    Prob &p = Fprod.p();
    calcBeliefF( I, p );
    p.normalize();

    return( Factor( factor(I).vars(), p ) );
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
            message( ni, I.iter ).fill( val );
            newMessage( ni, I.iter ).fill( val );
            updateResidual( ni, I.iter, 0.0 );
        }
    }
    _iters = 0;
}


void BP::updateMessage( size_t i, size_t _I ) {
    for (size_t j=0; j<_oldProd[i].size(); ++j) {
        _oldProd[i][j] =  _edges[i][_I].newMessage._p[j] * _oldProd[i][j] / _edges[i][_I].message._p[j];
    }


    if( recordSentMessages )
        _sentMessages.push_back(make_pair(i,_I));
    if( props.damping == 0.0 ) {
        message(i,_I) = newMessage(i,_I);
        updateResidual( i, _I, 0.0 );
    } else {
        message(i,_I) = (message(i,_I) ^ props.damping) * (newMessage(i,_I) ^ (1.0 - props.damping));
        updateResidual( i, _I, dist( newMessage(i,_I), message(i,_I), DISTLINF ) );
    }
}

// TODO: Optimize: We are using a heap now but this is not faster then the multimap solution. So we might have to revert to it.
void BP::updateResidual( size_t i, size_t _I, Real r ) {
    EdgeProp* pEdge = &_edges[i][_I];
    pEdge->residual = r;

    // rearrange look-up table (delete and reinsert new key)
    _lutNew.update(_edge2lutNew[i][_I], make_pair( r, make_pair(i, _I) ));
}


} // end of namespace dai


#endif
