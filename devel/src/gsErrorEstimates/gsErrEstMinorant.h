/** @file gsErrEstMinorant.h

@brief Functional-type error minorant for the Poisson problem.

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

/** \brief Provides computation of the element-wise error minorant of the approximate solution v constructed for
*          the Poisson-Dirichlet problem.
*
* Let the Poisson problem on the domain \f$ \Omega \f$ be given by
* \f[
-\Delta u = f, \quad
 u = u_D \mathrm{\ on\ } \Gamma
\f]
* where\n
* \f$ f \f$ is a given right-hand-side,\n
* \f$ u_D \f$ is given Dirichlet data on the Dirichlet boundary * \f$ \Gamma_D \f$.
*
* The error minorant \f$\underline{\rm m}^2\f$ for a computed discrete solution \f$ v \f$ is given by
* \f[ \underline{\rm m}^2 := \sum_K \underline{\rm m}^2_{K} \f]
* where the local estimate \f$ \underline{\rm m}^2_{K} \f$ on an element \f$ K \f$ is given by
 \f[
 \underline{\rm m}^2_{K} := \int_K ( \grad v - y )^2 dx
 \f]
* \f$ v \f$ denotes the approximate solution on the element \f$ \K \f$ and
* \f$ y \f$ denotes the approximate flux  on the element \f$ \K \f$.
*
* \ingroup Assembler
*/
    template <class T>
    class gsErrEstMinorant : public gsEstimator<T>
    {
        friend class gsEstimator<T>;

    public:

        /**
        * \brief Constructor with
        * \param _approxSol: Approximate solution.
        * \param _approxW: Auxiliary function, which physically mimics error of the approximate solution
        * and helps to reconstruct optimal error minorant.
        * \param _rhsFunction RHS-Source-function \f$ f \f$ of the Poisson problem.
        * \param bcInfo: Boundary conditions
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */

        /*
        gsErrEstMinorant(const gsField<T> & _approxSol,
                         const gsField<T> & _approxW,
                         const gsFunction<T> & _rhsFunction,
                         const gsBoundaryConditions<T> & bcInfo,
                         bool _rhsFunctParam = false)
                : gsEstimator<T>(_approxSol, _approxW, _rhsFunction),
                  m_bcInfo(bcInfo),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = true; // since bc are given
            m_storeElWiseType = 0;  // auxiliary flag, mainly for debugging
        }

         */

        /**
        * \brief Constructor with
        * \param _approxSol: Discrete solution
        * \param _approxW: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */
        gsErrEstMinorant(const gsField<T> & _approxSol,
                         const gsField<T> & _approxW,
                         const gsFunction<T> & _rhsFunction,
                         bool _rhsFunctParam = false)
                : gsEstimator<T>(_approxSol,_approxW, _rhsFunction),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = false;    // since no bc is given
            m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
        }

    public:

        /** \brief Computes the error minorant:
          * \f$\underline{\rm m}^2 := \sum_{K \in \Tau_K} \overline{\rm m}^2_{K}\f$, where
         * \f$\underline{\rm m}^2_{K} := \int_K ( 2 ( f * w - \nabla v * \nabla w   - \nabla w * \nabla w ) \dx\f$
         *
         * \param[in] storeElWise: Boolean flag, indicating whether the element-wise contributions must be stored.
         * If <em>storeEletWise = true</em>, the gsVector of element-wise estimates \f$\underline{\rm m}^2_{K}\f$
         * can be accessed by calling elementNorms().
         *
         * \returns The total value of the estimate \f$ \underline{\rm m}_d \f$.
         *
         */
        T compute(bool storeElWise = false, int elemNum = 0)
        {
            m_elWiseFull.clear();   // clear the vector of element-wise values
            m_storeElWiseType = unsigned( storeElWise );    // set the bool tag for whether to store element-wise values
            this->apply3Func(*this, storeElWise);   // m_storeElWiseType instead of storeElWise ?
            return this->m_value;
        }

        inline T takeRoot(const T v) {
            gsInfo << "minorant value before taking a square root = " << v << "\n";
            //return math::sqrt(math::abs(v));
            return math::sqrt(v);
        }

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
        * @param[in] approxSol
        * @param[in] approxFlux
        * @param[in] quNodes
        */
        inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                             const gsFunction<T>    & approxSol,
                             const gsFunction<T>    & approxW,
                             const gsFunction<T>    & rhsFunction,
                             gsMatrix<T>            & quNodes)
        {
            // temporary matrices for the parallel omp
            gsMatrix<T> m_approxSolVals_tmp, m_approxSolGrad_tmp,
                    m_approxWVals_tmp, m_approxWGrad_tmp,
                    m_rhsFunctVals_tmp;
            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            #pragma omp parallel sections num_threads(3)
            {
                #pragma omp section
                {
                    // Evaluate approximate solution's value
                    approxSol.eval_into(quNodes, m_approxSolVals_tmp);
                    // Evaluate approximate solution's gradient
                    approxSol.deriv_into(quNodes, m_approxSolGrad_tmp);
                }
                #pragma omp section
                {
                    // Evaluate approximate w's values
                    approxW.eval_into(quNodes, m_approxWVals_tmp);
                    // Evaluate approximate w's gradient
                    approxW.deriv_into(quNodes, m_approxWGrad_tmp);
                }
                #pragma omp section
                {
                    // Evaluate right-hand-side function (defined of physical domain)
                    rhsFunction.eval_into(geoEval.values(), m_rhsFunctVals_tmp);
                }
            }

            m_approxSolVals = m_approxSolVals_tmp;
            m_approxSolGrad = m_approxSolGrad_tmp;
            m_approxWVals = m_approxWVals_tmp;
            m_approxWGrad = m_approxWGrad_tmp;
            m_rhsFunctVals = m_rhsFunctVals_tmp;
        }

        /**
        * @brief Assembles (computes) the local error estimate contributions on the element.
        * See documentation of the class for the computed error estimate.
        *
        * @param[in] element: Specifies the element \f$ K \f$.
        * @param[in] geoEval: gsGeometryEvaluator as evaluated in evaluate().
        * @param[in] quWeights: Quadrature weights \em before transformation of the element,
        * i.e., the sum of the weights should be 1.
        * @return: The \em squared indicator \f$ \overline{\rm m}_{d, K}^2 \f$ of local error
        * on the element \f$ K \f$.
        */

        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & min) {
            // Access the index of the active patch
            // unsigned actPatch = geoEval.id();

            // Initialize the variable for the integral on the element
            T minSq(0.0);

            Eigen::setNbThreads(1);

            #pragma omp parallel num_threads(4)
            {
                gsMatrix<T> m_phApproxSolGrad, m_phApproxWGrad;    // grad(u_h)

                #pragma omp for reduction(+:minSq)
                for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
                {
                    const T weight = quWeights[k] * geoEval.measure(k);

                    // Map m_approxGrad from the param domain to the physical domain
                    geoEval.transformGradients(k, m_approxSolGrad, m_phApproxSolGrad);
                    // Map m_approxWGrad from the param domain to the physical domain
                    geoEval.transformGradients(k, m_approxWGrad, m_phApproxWGrad);

                    gsMatrix<T> grad_v_grad_v = m_phApproxSolGrad.transpose() * m_phApproxSolGrad;
                    gsMatrix<T> grad_w_grad_w = m_phApproxWGrad.transpose() * m_phApproxWGrad;

                    minSq += weight * ( 2.0 * m_rhsFunctVals(0, k) * (m_approxWVals(0, k) - m_approxSolVals(0, k))
                                        + ( grad_v_grad_v(0, 0) - grad_w_grad_w(0, 0)));

                }
            }
            Eigen::setNbThreads(0);
            min += minSq;

            return minSq;
        }

    public:

        const std::vector< std::vector<T> > & elementNormsFullData()  const
        {
            return m_elWiseFull;
        }

    private:

        gsBoundaryConditions<T> m_bcInfo;

        gsMatrix<T> m_approxSolVals;    // u_h
        gsMatrix<T> m_approxSolGrad;    // grad(u_h)
        gsMatrix<T> m_approxWVals;      // w_h
        gsMatrix<T> m_approxWGrad;      // grad(w_h)
        gsMatrix<T> m_rhsFunctVals;      // f_h

        unsigned m_parDim;              // dimension of the parameter domain

        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsNorm is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        using gsEstimator<T>::patchesPtr;    // pointer to all the patches

        unsigned m_storeElWiseType;                 // auxiliary flag, mainly for debugging
        std::vector< std::vector<T> > m_elWiseFull; // vector of the element-wise values
    };

} // namespace gismo
