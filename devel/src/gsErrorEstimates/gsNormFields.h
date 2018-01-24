/** @file gsNormFields.h

    @brief Provides generic routines for computing norms of the difference of two fields.
    Structure of the class is very close to sturcture of the gsNorm class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#pragma once
#include <gsErrorEstimates/gsEstimator.h>


namespace gismo
{

/** @brief Provides generic routines for computing norms of the residual functionals
 * element-wise or globally
 *
 * \ingroup ErrorEstimates
*/
    template <class T>
    class gsNormFields : public gsEstimator<T>
    {

        friend class gsEstimator<T>;

    public:

        gsNormFields(const gsField<T> & _field1,
                     const gsField<T> & _field2,
                     bool _rhsFunctParam = false)
                : gsEstimator<T>(_field1, _field2),
                  m_f2param(_rhsFunctParam){
            m_bcInitialized = false;    // no bc is given
            m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
        }


        /** \brief Computes the norm:
         * \f$ \| u_1 - u_2 \|^2 := \sum_{K \in \Tau_K} \| u_1 - u_2 \|^2_{K}\f$, where
         *
         * \param[in] storeElWise: Boolean flag, indicating whether the element-wise contributions must be stored.
         * If <em>storeEletWise = true</em>, the gsVector of element-wise estimates \f$\underline{\rm m}^2_{K}\f$
         * can be accessed by calling elementNorms().
         *
         * \returns The total value of the estimate \f$ \underline{\rm m}_d \f$.
         *
         */
        inline T compute(bool storeElWise = false, int elemNum = 0)
        {
            m_elWiseFull.clear();   // clear the vector of element-wise values
            m_storeElWiseType = unsigned( storeElWise );    // set the bool tag for whether to store element-wise values
            this->apply2Func(*this, storeElWise);   // m_storeElWiseType instead of storeElWise ?
            return this->m_value;
        }

        inline T takeRoot(const T v) { return math::sqrt(v);}

    protected:

        /**
         * @brief Initializes the error estimate, namely
         *
         * based on the degree of \em basis,
         * sets up \em rule (the quadrature rule) and
         * \em evFlags (for the gsGeometryEvaluator) that are needed for this specific problem.
         *
         * \param[in] basis
         * \param[out] rule
         * \param[out] evFlags
         */
        void initialize(const gsBasis<T> & basis,
                        gsQuadRule<T> & rule,
                        unsigned      & evFlags) // replace with geoEval ?
        {
            // Access the dimension of the geometry
            m_parDim = basis.dim();

            // Check if the dimension is 2 or 3
            // TODO: check that it works for dim = 1
            GISMO_ASSERT(m_parDim == 2 || m_parDim == 3,
                         "Called error estimator with dimension other than 2 or 3.");

            // Initialize a number of quadrature nodes for each of the sub-dimensions
            // based on the degree of a given basis
            gsVector<index_t> numQuadNodes( m_parDim );
            for (unsigned i = 0; i < m_parDim; ++i)
                numQuadNodes[i] = basis.degree(i) + 1; // number of the QuadNodes is + 1 more then the degree of basis

            // Initialize the Gauss rule for the integration
            // base on the number of quadrature nodes numQuadNodes
            rule = gsGaussRule<T>(numQuadNodes); // harmless slicing occurs here

            // Set Geometry evaluation flags, which is used in evaluate()
            evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN | NEED_2ND_DER | NEED_GRAD_TRANSFORM;
        }

        /**
        * @brief Evaluate needed data on the given quadrature nodes
        *
        * Executes and stores needed function evaluations at \em quNodes.\n
        * The gsGeometryEvaluator \em geoEval is also evaluated at the nodes,
        * using evaluation flags (evFlags) specified in initialize().
        *
        * @param[in, out] geoEval
        * @param[in] u1
        * @param[in] u2
        * @param[in] quNodes
        */
        inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                             const gsFunction<T>    & u1,
                             const gsFunction<T>    & u2,
                             gsMatrix<T>            & quNodes)
        {
            // temporary matrices for the parallel omp
            gsMatrix<T> m_u1_tmp, m_u2_tmp;

            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            //gsInfo << "geoEval measures:\n" << geoEval.measures() << "\n";

            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    // Evaluate u1's value
                    u1.eval_into(quNodes, m_u1_tmp);
                }
                #pragma omp section
                {
                    // Evaluate u2's value
                    u2.eval_into(quNodes, m_u2_tmp);
                }
            }

            m_u1 = m_u1_tmp;
            m_u2 = m_u2_tmp;

        }

        /**
        * @brief Assembles (computes) the local contributions on the element.
        * See documentation of the class for the computed error estimate.
        *
        * @param[in] element: Specifies the element \f$ K \f$.
        * @param[in] geoEval: gsGeometryEvaluator as evaluated in evaluate().
        * @param[in] quWeights: Quadrature weights \em before transformation of the element,
        * i.e., the sum of the weights should be 1.
        * @return: The contribution \f$ \| u_1 - u_2 \|_{K}^2 \f$ on the element \f$ K \f$.
        */

        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & total_val) {
            // Access the index of the active patch
            // unsigned actPatch = geoEval.id();

            // Initialize the variable for the integral on the element
            T incr_val(0.0);

            Eigen::setNbThreads(1);

            //#pragma omp parallel num_threads(4)
            {

                //#pragma omp for reduction(+:incr_val)
                for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
                {
                    const T weight = quWeights[k] * geoEval.measure(k);

                    //gsInfo << "m_u1.col(k) :" << m_u1.col(k) << "\n";
                    //gsInfo << "m_u2.col(k) :" << m_u2.col(k) << "\n";
                    //gsInfo << " quWeights[k] = " << quWeights[k] << "\n";
                    //gsInfo << " geoEval.measure(k) = " << geoEval.measure(k) << "\n";
                    //gsInfo << "( m_u1.col(k) - m_u2.col(k) ).squaredNorm() " << ( m_u1.col(k) - m_u2.col(k) ).squaredNorm() << "\n";
                    //gsInfo << "( m_u1.col(k) - m_u2.col(k))^2 " << ( m_u1.col(k) - m_u2.col(k) ) * ( m_u1.col(k) - m_u2.col(k) ) << "\n";

                    //gsInfo << "incr = " << weight * ( m_u1.col(k) - m_u2.col(k) ).squaredNorm() << "\n";

                    if (weight != 0)
                        incr_val += weight * ( m_u1.col(k) - m_u2.col(k) ).squaredNorm();

                    //gsInfo << "incr_val" << incr_val << "\n";

                }
            }
            Eigen::setNbThreads(0);
            total_val += incr_val;

            return incr_val;
        }

    public:

        const std::vector< std::vector<T> > & elementNormsFullData()  const
        {
            return m_elWiseFull;
        }

    private:

        gsBoundaryConditions<T> m_bcInfo;

        gsMatrix<T> m_u1;    // u_1
        gsMatrix<T> m_u2;    // u_2

        unsigned m_parDim;              // dimension of the parameter domain

        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsNorm is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        using gsEstimator<T>::patchesPtr;    // pointer to all the patches

        unsigned m_storeElWiseType;                 // auxiliary flag, mainly for debugging
        std::vector< std::vector<T> > m_elWiseFull; // vector of the element-wise values
    };

} // namespace gismo

