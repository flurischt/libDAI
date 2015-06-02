/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


/// \file
/// \brief Defines class BP, which implements (Loopy) Belief Propagation
/// \todo Consider using a priority_queue for maximum residual schedule


#ifndef __defined_libdai_bp_h
#define __defined_libdai_bp_h


#include <dai/dai_config.h>
#ifdef DAI_WITH_BP


#include <string>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/enum.h>

#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/skew_heap.hpp>
#include <boost/heap/priority_queue.hpp>

//#include "immintrin.h"
//#include "avxintrin.h"

namespace dai {


/// Approximate inference algorithm "(Loopy) Belief Propagation"
/** The Loopy Belief Propagation algorithm uses message passing
 *  to approximate marginal probability distributions ("beliefs") for variables
 *  and factors (more precisely, for the subset of variables depending on the factor).
 *  There are two variants, the sum-product algorithm (corresponding to
 *  finite temperature) and the max-product algorithm (corresponding to
 *  zero temperature).
 *
 *  The messages \f$m_{I\to i}(x_i)\f$ are passed from factors \f$I\f$ to variables \f$i\f$.
 *  In case of the sum-product algorith, the update equation is:
 *    \f[ m_{I\to i}(x_i) \propto \sum_{x_{N_I\setminus\{i\}}} f_I(x_I) \prod_{j\in N_I\setminus\{i\}} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}\f]
 *  and in case of the max-product algorithm:
 *    \f[ m_{I\to i}(x_i) \propto \max_{x_{N_I\setminus\{i\}}} f_I(x_I) \prod_{j\in N_I\setminus\{i\}} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}\f]
 *  In order to improve convergence, the updates can be damped. For improved numerical stability,
 *  the updates can be done in the log-domain alternatively.
 *
 *  After convergence, the variable beliefs are calculated by:
 *    \f[ b_i(x_i) \propto \prod_{I\in N_i} m_{I\to i}(x_i)\f]
 *  and the factor beliefs are calculated by:
 *    \f[ b_I(x_I) \propto f_I(x_I) \prod_{j\in N_I} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}(x_j) \f]
 *  The logarithm of the partition sum is calculated by:
 *    \f[ \log Z = \sum_i (1 - |N_i|) \sum_{x_i} b_i(x_i) \log b_i(x_i) - \sum_I \sum_{x_I} b_I(x_I) \log \frac{b_I(x_I)}{f_I(x_I)} \f]
 *
 *  For the max-product algorithm, a heuristic way of finding the MAP state (the
 *  joint configuration of all variables which has maximum probability) is provided
 *  by the findMaximum() method, which can be called after convergence.
 *
 *  \note There are two implementations, an optimized one (the default) which caches IndexFor objects,
 *  and a slower, less complicated one which is easier to maintain/understand. The slower one can be
 *  enabled by defining DAI_BP_FAST as false in the source file.
 */
    class BP : public DAIAlgFG {
    protected:

        /// Type used for products of probabilities.
        /// Indifferent to presence of flag DAI_SINGLE_PRECISION.
        typedef ProbD ProbProduct;

        /// Type used for index cache
        typedef std::vector<size_t> ind_t;
        /// Type used for storing edge properties
        typedef Real MessageType;
        struct EdgeProp {
            /// Index cached for this edge
            size_t      index;
            /// Old message living on this edge
            MessageType message;
            /// New message living on this edge
            MessageType newMessage;
            /// Residual for this edge
            Real        residual;
            /// Precalculated reciprocals of the message for this edge:
            /// 1/message and 1/(1-message).
            Real reciprocals[2];
        };

        /// Stores all edge properties
        std::vector<std::vector<EdgeProp> > _edges;

        /// Stores the pre-calculated indices for the edges.
        std::vector<ind_t> _indices;

        // We store the product for each variable. Every time a message gets
        // updated we also update the corresponding product. We can then reuse
        // the result and make the algorithm much faster.
        // Use double precision!
        // TODO: use std::vector<ProbProd> (for consistent notation)
        std::vector< std::vector<double> > _oldProd;

        TProb<double> _factorsFixed;

        __m256d _prod_vec;
        double _prod_double[4];

#ifdef DAI_SINGLE_PRECISION
        ProbProduct _marg;
#endif

        // Buffer for simple calculations.
        mutable ProbProduct _probTemp;

        /// Type of lookup table (only used for maximum-residual BP)
        typedef std::multimap<Real, std::pair<size_t, size_t> > LutType;
        /// Lookup table (only used for maximum-residual BP)
        std::vector<std::vector<LutType::iterator> > _edge2lutOld;
        typedef std::pair<Real, std::pair<size_t, size_t>> heap_data;
        typedef boost::heap::pairing_heap<heap_data, boost::heap::mutable_<true>>::handle_type heap_data_handle;
        std::vector<std::vector<heap_data_handle> > _edge2lutNew;
        /// Lookup table (only used for maximum-residual BP)
        LutType _lut;
        boost::heap::pairing_heap<heap_data, boost::heap::mutable_<true>> _lutNew;
        /// Maximum difference between variable beliefs encountered so far
        Real _maxdiff;
        /// Number of iterations needed
        size_t _iters;
        /// The history of message updates (only recorded if \a recordSentMessages is \c true)
        std::vector<std::pair<size_t, size_t> > _sentMessages;
        /// Stores variable beliefs of previous iteration
        std::vector<Factor> _oldBeliefsV;
        /// Stores factor beliefs of previous iteration
        std::vector<Factor> _oldBeliefsF;
        /// Stores the update schedule
        std::vector<Edge> _updateSeq;

    public:
        /// Parameters for BP
        struct Properties {
            /// Enumeration of possible update schedules
            /** The following update schedules have been defined:
             *  - PARALL parallel updates
             *  - SEQFIX sequential updates using a fixed sequence
             *  - SEQRND sequential updates using a random sequence
             *  - SEQMAX maximum-residual updates [\ref EMK06]
             * \obsolete BPFast uses always SEQMAX.
             */
            DAI_ENUM(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL);

            /// Enumeration of inference variants
            /** There are two inference variants:
             *  - SUMPROD Sum-Product
             *  - MAXPROD Max-Product (equivalent to Min-Sum)
             */
            DAI_ENUM(InfType,SUMPROD,MAXPROD);

            /// Verbosity (amount of output sent to stderr)
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Maximum time (in seconds)
            double maxtime;

            /// Tolerance for convergence test
            Real tol;

            /// Damping constant (0.0 means no damping, 1.0 is maximum damping)
            Real damping;

            size_t specialFactors;

        } props;

        /// Specifies whether the history of message updates should be recorded
        bool recordSentMessages;
        int messageCount;

    public:
        /// \name Constructors/destructors
        //@{
        /// Default constructor
        BP() : DAIAlgFG()
                , _edges()
                , _indices()
                , _edge2lutOld()
                , _lut()
                , _maxdiff(0.0)
                , _iters(0U)
                , _sentMessages()
                , _oldBeliefsV()
                , _oldBeliefsF()
                , _updateSeq()
                , props()
                , recordSentMessages(false)
        {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        /** \param fg Factor graph.
         *  \param opts Parameters @see Properties
         */
        BP( const FactorGraph & fg, const PropertySet &opts ) : DAIAlgFG(fg)
                , _edges()
                , _indices()
                , _maxdiff(0.0)
                , _iters(0U)
                , _sentMessages()
                , _oldBeliefsV()
                , _oldBeliefsF()
                , _updateSeq()
                , props()
                , recordSentMessages(false) {
            setProperties( opts );
            construct();
        }

        /// Copy constructor
        BP( const BP &x ) : DAIAlgFG(x)
                , _edges(x._edges)
                , _indices(x._indices)
                , _oldProd(x._oldProd)
                , _edge2lutOld(x._edge2lutOld)
                , _lut(x._lut)
                , _maxdiff(x._maxdiff)
                , _iters(x._iters)
                , _sentMessages(x._sentMessages)
                , _oldBeliefsV(x._oldBeliefsV)
                , _oldBeliefsF(x._oldBeliefsF)
                , _updateSeq(x._updateSeq)
                , props(x.props)
                , recordSentMessages(x.recordSentMessages) {
            for( LutType::iterator l = _lut.begin(); l != _lut.end(); ++l )
                _edge2lutOld[l->second.first][l->second.second] = l;
        }

        /// Assignment operator
        BP& operator=( const BP &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _edges = x._edges;
                _indices = x._indices;
                _oldProd = x._oldProd;
                _lut = x._lut;
                for( LutType::iterator l = _lut.begin(); l != _lut.end(); ++l )
                    _edge2lutOld[l->second.first][l->second.second] = l;
                _maxdiff = x._maxdiff;
                _iters = x._iters;
                _sentMessages = x._sentMessages;
                _oldBeliefsV = x._oldBeliefsV;
                _oldBeliefsF = x._oldBeliefsF;
                _updateSeq = x._updateSeq;
                props = x.props;
                recordSentMessages = x.recordSentMessages;
            }
            return *this;
        }
        //@}

        /// \name General InfAlg interface
        //@{
        virtual BP* clone() const { return new BP(*this); }
        virtual BP* construct( const FactorGraph &fg, const PropertySet &opts ) const { return new BP( fg, opts ); }
        virtual std::string name() const { return "BP"; }
        virtual Factor belief( const Var &v ) const { return beliefV( findVar( v ) ); }
        virtual Factor belief( const VarSet &vs ) const;
        virtual Factor beliefV( size_t i ) const;
        virtual Factor beliefF( size_t I ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        /** \pre Assumes that run() has been called and that \a props.inference == \c MAXPROD
         */
        std::vector<size_t> findMaximum() const { return dai::findMaximum( *this ); }
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        virtual void setMaxIter( size_t maxiter ) { props.maxiter = maxiter; }
        virtual void setProperties( const PropertySet &opts );
        virtual PropertySet getProperties() const;
        virtual std::string printProperties() const;
        //@}

        /// \name Additional interface specific for BP
        //@{
        /// Returns history of which messages have been updated
        const std::vector<std::pair<size_t, size_t> >& getSentMessages() const {
            return _sentMessages;
        }

        /// Clears history of which messages have been updated
        void clearSentMessages() { _sentMessages.clear(); }
        //@}

    protected:
        /// Returns constant reference to message from the \a _I 'th neighbor of variable \a i to variable \a i
        const MessageType & message(size_t i, size_t _I) const { return _edges[i][_I].message; }
        /// Returns reference to message from the \a _I 'th neighbor of variable \a i to variable \a i
        MessageType & message(size_t i, size_t _I) { return _edges[i][_I].message; }
        /// Returns constant reference to updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        const MessageType & newMessage(size_t i, size_t _I) const { return _edges[i][_I].newMessage; }
        /// Returns reference to updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        MessageType & newMessage(size_t i, size_t _I) { return _edges[i][_I].newMessage; }
        /// Returns constant reference to cached index for the edge between variable \a i and its \a _I 'th neighbor
        const ind_t & index(size_t i, size_t _I) const { return _indices[_edges[i][_I].index]; }
        /// Returns reference to cached index for the edge between variable \a i and its \a _I 'th neighbor
        ind_t & index(size_t i, size_t _I) { return _indices[_edges[i][_I].index]; }
        /// Returns constant reference to residual for the edge between variable \a i and its \a _I 'th neighbor
        const Real & residual(size_t i, size_t _I) const { return _edges[i][_I].residual; }
        /// Returns reference to residual for the edge between variable \a i and its \a _I 'th neighbor
        Real & residual(size_t i, size_t _I) { return _edges[i][_I].residual; }

        /// Calculate the product of factor \a I and the incoming messages
        /** If \a without_i == \c true, the message coming from variable \a i is omitted from the product
         *  \note This function is used by calcNewMessage() and calcBeliefF()
         */
        virtual void calcIncomingMessageProduct(ProbProduct &prod, size_t I, bool without_i, size_t i) const;

        /// Specialised versions of calcIncomingMessageProduct for special patterns.
        /// Implementation in bp_ext.cpp.
        void calcIncomingMessageProduct_0101_0011(__m256d& prod_vec, size_t I, size_t i) const;
        void marginalizeProductOntoMessage(__m256d& prod_vec, size_t i, size_t _I, size_t prodsize);

        void calcIncomingMessageProduct_0101_0011(double* prod, size_t I, size_t i) const;
        void marginalizeProductOntoMessage(double* prod, size_t i, size_t _I, size_t prodsize);

        static const int INDEX_0011      = 0;
        static const int INDEX_0101      = 1;
        static const int INDEX_01        = 2;
        static const int INDEX_0101_0011 = 3;


        /// Calculate the updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        virtual void calcNewMessage( size_t i, size_t _I);
        /// Replace the "old" message from the \a _I 'th neighbor of variable \a i to variable \a i by the "new" (updated) message
        void updateMessage( size_t i, size_t _I );
        /// Set the residual (difference between new and old message) for the edge between variable \a i and its \a _I 'th neighbor to \a r
        void updateResidual( size_t i, size_t _I, Real r );
        /// Finds the edge which has the maximum residual (difference between new and old message)
        bool findMaxResidual( size_t &i, size_t &_I );
        /// Calculates unnormalized belief of variable \a i
        virtual void calcBeliefV(size_t i, ProbProduct &p ) const;
        /// Calculates unnormalized belief of factor \a I
        virtual void calcBeliefF( size_t I, ProbProduct &p ) const {
            calcIncomingMessageProduct(p, I, false, 0);
        }

        /// Helper function for constructors
        virtual void construct();
    };


} // end of namespace dai


#endif


#endif
