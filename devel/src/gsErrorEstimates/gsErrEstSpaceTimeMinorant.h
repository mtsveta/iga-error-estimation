/** @file gsErrEstSpaceTimeMinorant.h

@brief A term (measuring the residual in the equilibrium equation) of the functional-type error estimator
 for the space-time formulation of the heat problem.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): S. Matculevich
*/

#pragma once

#include <gsCore/gsEvalUtils.hpp>
#include <gsErrorEstimates/gsEstimator.h>

namespace gismo
{

/** \brief Provides computation of equilibrium term of fucntional-type element-wise error indicator for the Poisson problem.
*
* Let the Poisson problem on the domain \f$ Q := \Omega \times (0, T)\f$ be given by
* \f[
 u_t - \Delta u = f,  \quad \mathrm{\in\ } Q
 u = u_D,             \quad  \mathrm{\on\ } \partial \Omega \cup [0, T]
 u(0, x) = u_0,       \quad \mathrm{\on\ } \Omega \times {0}
\f]
* where\n
* \f$ f \f$ is a given right-hand-side,\n
* \f$ u_D \f$ is given Dirichlet data on the Dirichlet boundary * \f$ \partial \Omega \cup [0, T] \f$.
* \f$ u_0 \f$ is given initial data at the * \f$ \Omega \times {0} \f$.
*
* The error indicator \f$\overline{\rm m}^2_eq\f$ for a computed discrete solution \f$ v \f$ is given by
* \f[ \overline{\rm m}^2_eq := \sum_K \overline{\rm eq}^2_{d, K} \f]
* where the local estimate \f$ \overline{\rm eq}^2_{eq, K} \f$ on an element \f$ K \f$ is given by
 \f[
 \overline{\rm m}^2_{eq, K} := \int_K ( div y + f - v_t )^2 dx
 \f]
* \f$ f \f$ is the RHS of the heat equation,
* \f$ y \f$ is the approximate flux  on the element \f$ \K \f$, and
* \f$ v \f$ is the approximate solution on the element \f$ \K \f$.
*
* \ingroup Assembler
*/
    template <class T>
    class gsErrEstSpaceTimeMinorant : public gsEstimator<T>
    {
        friend class gsEstimator<T>;

    public:

        /**
        * \brief Constructor with
        * \param _approxSol: Approximate solution
        * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _approxW: Auxiliary function, which physically mimics the approximate solution but reconstructed on
        * reacher space.
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */
        gsErrEstSpaceTimeMinorant(const gsField<T> & _approxV,
                                       const gsField<T> & _approxW,
                                       const gsFunction<T> & _rhsFunction,
                                       bool _rhsFunctParam = false)
                : gsEstimator<T>(_approxV,_approxW, _rhsFunction),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = false;    // since no bc is given
            m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
        }


    public:

        /** \brief Computes the error indicator, namely
         * computes the term of the functional-type error estimate \f$\overline{\rm m}_eq\f$ that measures
         * the error in balance equation, i.e.,
         * \f$\overline{\rm m}^2_eq := \sum_{K \in \Tau_K} \overline{\rm m}^2_{eq, K}\f$, where
         * \f$\overline{\rm m}^2_{eq, K} := \int_K ( div y + f - v_t )^2 \dx\f$
         *
         * \param[in] storeElWise: Boolean flag, indicating whether the element-wise contributions must be stored.
         * If <em>storeEletWise = true</em>, the gsVector of element-wise estimates \f$\overline{\rm m}^2_{eq, K}\f$
         * can be accessed by calling elementNorms().
         *
         * \returns The total value of the estimate \f$ \overline{\rm m}_eq \f$.
         *
         */

        inline T compute(bool storeElWise, int elemNum)
        {
            m_elWiseFull.clear();
            m_storeElWiseType = unsigned( storeElWise );

            this->apply3Func(*this,storeElWise);
            return this->m_value;
        }

        inline T takeRoot(const T v) { return v;}

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
            GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Called error estimator with dimension other than 2 or 3.");
            // whey it cannot be 1?

            // Initialize a number of quadrature nodes for each of the sub-dimenisions based on the
            // degree of a given basis
            gsVector<index_t> numQuadNodes( m_parDim );
            for (unsigned i = 0; i < m_parDim; ++i)
                numQuadNodes[i] = basis.degree(i) + 1; // number of the QuadNodes is 1 more then the degree of basis

            // Initialize the Gauss rule for the integration
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
        * @param[in] approxSol
        * @param[in] approxFlux
        * @param[in] quNodes
        */
        inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                             const gsFunction<T>    & approxV,
                             const gsFunction<T>    & approxW,
                             const gsFunction<T>    & rhsFunction,
                             gsMatrix<T>            & quNodes)
        {

            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            // Evaluate approximate solution's gradient
            approxV.eval_into(quNodes, m_approxV);
            approxV.deriv_into(quNodes, m_approxGradV);

            // Evaluate approximate flux's values
            approxW.eval_into(quNodes, m_approxW);
            approxW.deriv_into(quNodes, m_approxGradW);

            // Evaluate right-hand-side function (defined of physical domain)
            rhsFunction.eval_into(geoEval.values(), m_rhsVals);

        }

        /**
        * @brief Assembles (computes) the local error estimate contributions on the element.
        * See documentation of the class for the computed error estimate.
        *
        * @param[in] element: Specifies the element \f$ K \f$.
        * @param[in] geoEval: gsGeometryEvaluator as evaluated in evaluate().
        * @param[in] quWeights: Quadrature weights \em before transformation of the element,
        * i.e., the sum of the weights should be 1.
        * @return: The \em squared indicator \f$ \overline{\rm m}_{eq, K}^2 \f$ of local error
        * on the element \f$ K \f$.
        */

        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & total)
        {
            // Initialize the variable for the integral on the element and it's facets
            T incrSq(0.0);

            const std::pair<int,int> dim = geoEval.geometry().function(0).dimensions();
            const int d = dim.first;

            gsMatrix<T> m_phApproxGradV;        // grad (v) on physical domain
            gsMatrix<T> m_phApproxSpaceGradV;   // grad_x (v) on physical domain
            gsMatrix<T> m_phApproxTimeGradV;    // grad_t (v) on physical domain

            gsMatrix<T> m_phApproxGradW;        // grad (v) on physical domain
            gsMatrix<T> m_phApproxSpaceGradW;   // grad_x (v) on physical domain
            gsMatrix<T> m_phApproxTimeGradW;    // grad_t (v) on physical domain

            //#pragma omp for reduction(+:meqSq)
            for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights[k] * geoEval.measure(k);

                geoEval.transformGradients(k, m_approxGradV, m_phApproxGradV);
                geoEval.transformGradients(k, m_approxGradW, m_phApproxGradW);

                m_phApproxSpaceGradV = m_phApproxGradV.topRows(d-1);
                m_phApproxTimeGradV = m_phApproxGradV.bottomRows(1);

                m_phApproxSpaceGradW = m_phApproxGradW.topRows(d-1);
                m_phApproxTimeGradW = m_phApproxGradW.bottomRows(1);

                T G11 = - m_phApproxSpaceGradV.col(0).transpose() * (m_phApproxSpaceGradW.col(0) - m_phApproxSpaceGradV.col(0));
                T G12 = - 0.5 * (m_phApproxSpaceGradW.col(0) - m_phApproxSpaceGradV.col(0)).transpose() * (m_phApproxSpaceGradW.col(0) - m_phApproxSpaceGradV.col(0));
                T G21 = m_approxV(0, k) * (m_approxW(0, 0) - m_approxV(0, 0));
                T G22 = - 0.5 * (m_approxW(0, 0) - m_approxV(0, 0)) * (m_approxW(0, 0) - m_approxV(0, 0));
                T F   = m_rhsVals(0, k) * (m_approxW(0, k) - m_approxV(0, k));

                //gsInfo << "G11 = " << G11 << "\n";
                //gsInfo << "G12 = " << G12 << "\n";
                //gsInfo << "G21 = " << G21 << "\n";
                //gsInfo << "G22 = " << G22 << "\n";
                //gsInfo << "F   = " << F << "\n";

                incrSq += weight * (G11 + G12 + G21 + G22  + F);
                // min_K(v, w) := { - grad_x v * grad_x (w - v) - 1/2 * |grad_x (w - v)|^2 + v * (w - v)_t - 1/2 * |(w - v)_t|^2 + f * (w - v) }_K
                /*incrSq += weight * ( - m_phApproxSpaceGradV.col(0).transpose() * (m_phApproxSpaceGradW.col(0) - m_phApproxSpaceGradV.col(0))
                                     - 0.5 * (m_phApproxSpaceGradW.col(0) - m_phApproxSpaceGradV.col(0)).transpose() * (m_phApproxSpaceGradW.col(0) - m_phApproxSpaceGradV.col(0))
                                     + m_approxV(0, k) * (m_phApproxTimeGradW(0, 0) - m_phApproxTimeGradV(0, 0))
                                     - 0.5 * (m_phApproxTimeGradW(0, 0) - m_phApproxTimeGradV(0, 0)) * (m_phApproxTimeGradW(0, 0) - m_phApproxTimeGradV(0, 0))
                                     + m_rhsVals(0, k) * (m_approxW(0, k) - m_approxV(0, k)));
                                     */
            }
            total += incrSq;
            return incrSq;
        }

    public:

        const std::vector< std::vector<T> > & elementNormsFullData()  const
        {
            return m_elWiseFull;
        }

    private:

        gsBoundaryConditions<T> m_bcInfo;

        gsMatrix<T> m_approxV, m_approxGradV;   // v, grad(v)
        gsMatrix<T> m_approxW, m_approxGradW;   // w, grad(w)
        gsMatrix<T> m_rhsVals;                  // f_h

        unsigned m_parDim;              // dimension

        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsNorm is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        //using gsNorm<T>::patchesPtr;    // pointer to all the patches

        // auxiliary flag, mainly for debugging
        unsigned m_storeElWiseType;
        std::vector< std::vector<T> > m_elWiseFull;

    };

} // namespace gismo
