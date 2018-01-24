/** @file gsErrEstDualSpaceTimeMajorantII.h

    @brief

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
    class gsErrEstDualSpaceTimeMajorantII : public gsEstimator<T>
    {
        friend class gsEstimator<T>;

    public:

        gsErrEstDualSpaceTimeMajorantII(const gsField<T> & _field1,
                                        const gsField<T> & _field2,
                                        const gsField<T> & _field3,
                     bool _rhsFunctParam = false)
                : gsEstimator<T>(_field1, _field2, _field3),
                  m_f2param(_rhsFunctParam){
            m_bcInitialized = false;    // no bc is given
            //m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
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
            this->apply3Fields(*this, storeElWise);   // m_storeElWiseType instead of storeElWise ?
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
                             const gsFunction<T>    & v,
                             const gsFunction<T>    & flux,
                             const gsFunction<T>    & w,
                             gsMatrix<T>            & quNodes)
        {
            // temporary matrices for the parallel omp
            gsMatrix<T> m_flux_tmp, m_gradv_tmp, m_gradw_tmp;

            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            #pragma omp parallel sections
            {
                #pragma omp section
                flux.eval_into(quNodes, m_flux_tmp);

                #pragma omp section
                v.deriv_into(quNodes, m_gradv_tmp);

                #pragma omp section
                w.deriv_into(quNodes, m_gradw_tmp);

            }

            m_fluxVals = m_flux_tmp;
            m_approxGradV = m_gradv_tmp;
            m_approxGradW = m_gradw_tmp;

            //gsInfo << "m_fluxVals : \n" << m_fluxVals << "\n";
            //gsInfo << "m_approxGradV : \n" << m_approxGradV << "\n";
            //gsInfo << "m_approxGradW : \n" << m_approxGradW << "\n";
            //gsInfo << "\n";
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

            const std::pair<int,int> dim = geoEval.geometry().function(0).dimensions();
            const int d = dim.first;

            gsMatrix<T> m_phApproxGradV;             // grad(u_h)
            gsMatrix<T> m_phApproxSpaceGradV;        // grad_x(u_h)
            gsMatrix<T> m_phApproxGradW;             // grad(w_h)
            gsMatrix<T> m_phApproxSpaceGradW;        // grad_x(w_h)
            //#pragma omp parallel num_threads(4)

            m_phApproxGradV.setZero(d, 1);
            m_phApproxGradW.setZero(d, 1);

            for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights[k] * geoEval.measure(k);

                // Map m_approxGrad from the param domain to the physical domain
                geoEval.transformGradients(k, m_approxGradV, m_phApproxGradV);
                geoEval.transformGradients(k, m_approxGradW, m_phApproxGradW);
                //gsInfo << "m_phApproxGradV : \n" << m_phApproxGradV << "\n";
                //gsInfo << "m_phApproxGradW : \n" << m_phApproxGradW << "\n";

                // get grad_x(v) from grad(v)
                m_phApproxSpaceGradV.setZero(d, 1);
                m_phApproxSpaceGradV.topRows(d - 1) = m_phApproxGradV.topRows(d - 1); // copy just top d-1 rows

                // get grad_x(w) from grad(w)
                m_phApproxSpaceGradW.setZero(d, 1);
                m_phApproxSpaceGradW.topRows(d - 1) = m_phApproxGradW.topRows(d - 1); // copy just top d-1 rows

                //gsInfo << "m_phApproxSpaceGradV : \n" << m_phApproxSpaceGradV << "\n";
                //gsInfo << "m_phApproxSpaceGradW : \n" << m_phApproxSpaceGradW << "\n";

                if (weight != 0)
                    incr_val += weight * (m_phApproxSpaceGradW + m_fluxVals.col(k) - 2 * m_phApproxSpaceGradV).squaredNorm();
            }
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

        gsMatrix<T> m_fluxVals;         // m_flux
        gsMatrix<T> m_approxGradV;      // m_gradv
        gsMatrix<T> m_approxGradW;      // m_gradw

        unsigned m_parDim;              // dimension of the parameter domain

        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsNorm is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        using gsEstimator<T>::patchesPtr;    // pointer to all the patches

        unsigned m_storeElWiseType;                 // auxiliary flag, mainly for debugging
        std::vector< std::vector<T> > m_elWiseFull; // vector of the element-wise values
    };

} // namespace gismo

