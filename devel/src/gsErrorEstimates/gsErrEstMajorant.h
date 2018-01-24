/** @file gsErrEstDualMajorant.h

@brief Functional-type error estimator for the Poisson problem.

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

/** \brief Provides computation of dual part of the fucntional-type element-wise error indicator for the Poisson problem.
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
* The error estimate \f$\overline{\rm M}^2\f$ for a computed discrete solution \f$ v \f$ is given by
* \f[ \overline{\rm M}^2 := \sum_K \overline{\rm M}^2_{d, K} \f]
* where the local estimate \f$ \overline{\rm M}^2_{d, K} \f$ on an element \f$ K \f$ is given by
 \f[
 \overline{\rm M}^2_{d, K} := \int_K ((1 + \beta) \overline{\rm m}^2_{d, K}
         *                             + (1 + 1/\beta)\overline{\rm m}^2_{eq, K})
 \f]
* \f$ v \f$ denotes the approximate solution on the element \f$ \K \f$ and
* \f$ y \f$ denotes the approximate flux  on the element \f$ \K \f$.
*
* \ingroup Assembler
*/
    template <class T>
    class gsErrEstMajorant : public gsEstimator<T>
    {
        friend class gsEstimator<T>;

    public:

        /**
        * \brief Constructor with
        * \param _approxSol: Approximate solution.
        * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _rhsFunction RHS-Source-function \f$ f \f$ of the Poisson problem.
        * \param bcInfo: Boundary conditions
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */

        gsErrEstMajorant(const gsField<T> & _approxSol,
                         const gsField<T> & _approxFlux,
                         const gsFunction<T> & _rhsFunction,
                         const T & _cFriedrichs,
                         const T & _beta,
                         const gsBoundaryConditions<T> & bcInfo,
                         bool _rhsFunctParam = false)
                : gsEstimator<T>(_approxSol, _approxFlux, _rhsFunction),
                  m_cFriedrichs(_cFriedrichs),
                  m_beta(_beta),
                  m_bcInfo(bcInfo),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = true; // since bc are given
            m_storeElWiseType = 0;  // auxiliary flag, mainly for debugging
        }


        /**
        * \brief Constructor with
        * \param _approxSol: Approximate solution
        * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _rhsFunction RHS-Source-function \f$ f \f$ of the Poisson problem.
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */
        gsErrEstMajorant(const gsField<T> & _approxSol,
                         const gsField<T> & _approxFlux,
                         const gsFunction<T> & _rhsFunction,
                         const T & _cFriedrichs, const T & _beta,
                         bool _rhsFunctParam = false)
                : gsEstimator<T>(_approxSol, _approxFlux, _rhsFunction),
                  m_cFriedrichs(_cFriedrichs),
                  m_beta(_beta),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = false;    // since no bc is given
            m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
        }

    public:

        /** \brief Computes the error indicator, namely
         * computes the functional-type error estimate
         * \f$\overline{\rm M} := (1 + \beta) \overline{\rm m}^2_d + (1 + 1/\beta) \overline{\rm m}^2_eq
         *                     := (1 + \beta) \sum_{K \in \Tau_K} \overline{\rm m}^2_{d, K}
         *                        + (1 + 1/\beta) \sum_{K \in \Tau_K} \overline{\rm m}^2_{eq, K} \f$,
         * where
         * \f$\overline{\rm m}^2_{d, K} := \int_K ( \nabla v - y )^2 \dx\f$, and
         * \f$\overline{\rm m}^2_{eq, K} := \int_K ( \div y +f )^2 \dx\f$.

         *
         * \param[in] storeElWise: Boolean flag, indicating whether the element-wise contributions must be stored.
         * If <em>storeEletWise = true</em>, the gsVector of element-wise estimates \f$\overline{\rm m}^2_{d, K}\f$
         * can be accessed by calling elementNorms().
         *
         * \returns The total value of the estimate \f$ \overline{\rm M}\f$.
         *
         */
        T compute(bool storeElWise = true)
        {
            m_elWiseFull.clear();                           // clear the vector of element-wise values
            m_storeElWiseType = unsigned( storeElWise );    // set the bool tag for whether to store element-wise values
            this->apply3Func(*this,storeElWise);            // m_storeElWiseType instead of storeElWise ?
            return this->m_value;
        }

        T computeMajorant(bool storeElWise = true)
        {
            m_elWiseFull.clear();                           // clear the vector of element-wise values
            m_storeElWiseType = unsigned( storeElWise );    // set the bool tag for whether to store element-wise values
            this->apply3Func(*this,storeElWise);            // m_storeElWiseType instead of storeElWise ?
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
        * @param[in] approxSol
        * @param[in] approxFlux
        * @param[in] quNodes
        */
        inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                             const gsFunction<T>    & approxSol,
                             const gsFunction<T>    & approxFlux,
                             const gsFunction<T>    & rhsFunction,
                             gsMatrix<T>            & quNodes)
        {
            // Evaluate approximate solution's gradient
            approxSol.deriv_into(quNodes, m_approxSolGrad);
            // Evaluate approximate flux's values and divergence of the fluxes values
            approxFlux.eval_into(quNodes, m_approxFluxVals);
            approxFlux.div_into(quNodes, m_approxFluxDiv);

            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            // Evaluate right-hand-side function (defined of physical domain)
            rhsFunction.eval_into(geoEval.values(), m_rhsFunctVals);

            //gsInfo << "quNodes:\n"   << quNodes << "\n";
            //gsInfo << "geoEval.values():\n"   << geoEval.values() << "\n";
            //gsInfo << "geoEval.numPoints() in evaluate:\n"   << geoEval.numPoints() << "\n";

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
        /*
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & maj)
        {
            // Access the index of the active patch
            unsigned actPatch = geoEval.id();

            // Initialize the variable for the integral on the element
            T mdSq(0.0), meqSq(0.0), majSq(0.0);

            gsVector<index_t> numQuadNodesRef( m_parDim );  // reference-vector with default num. of quad. points

            // Run through all the dimensions collecting the amount of quadrature nodes
            for (unsigned i = 0; i < m_parDim; ++i)
                numQuadNodesRef[i] = patchesPtr->basis(actPatch).degree(i) + 1;

            //gsInfo << "m_approxSolGrad:\n" << m_approxSolGrad << "\n";

            for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights[k] * geoEval.measure(k);

                // TODO: should one use m_phApproxSolGrad? ask Stefan...
                //geoEval.transformGradients(k, m_approxSolGrad, m_phApproxSolGrad); // should one use m_phApproxSolGrad
                //gsInfo << "m_phApproxSolGrad:\n" << m_phApproxSolGrad << "\n";

                // |grad v - y|^2 := (v_{x1} - y_1)^2 + (v_{x1} - y_2)^2
                mdSq += weight * (math::pow(m_approxSolGrad(0, k) - m_approxFluxVals(0, k), 2.0 )
                                  + math::pow(m_approxSolGrad(1, k) - m_approxFluxVals(1, k), 2.0 ));
                // |\div y + f|^2_K
                meqSq += weight * math::pow( m_approxFluxDiv(0,k) + m_rhsFunctVals(0,k), 2.0);
            }

            // Add into the total value of the indicator

            majSq = (1.0 + m_beta) * mdSq + (1.0 + 1.0 / m_beta) * m_cFriedrichs * meqSq;
            maj     += majSq;

            return maj;
        }
        */
        inline T compute(gsDomainIterator<T>    & element,
                                 gsGeometryEvaluator<T> & geoEval,
                                 gsVector<T> const      & quWeights,
                                 T & maj, T & mD, T & mEq)
        {
            // Access the index of the active patch
            unsigned actPatch = geoEval.id();

            // Initialize the variable for the integral on the element
            T mdSq(0.0), meqSq(0.0), majSq(0.0);

            //gsInfo << "m_approxSolGrad:\n" << m_approxSolGrad << "\n";

            for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
            {
                const T weight = quWeights[k] * geoEval.measure(k);

                // TODO: should one use m_phApproxSolGrad? ask Stefan...
                //geoEval.transformGradients(k, m_approxSolGrad, m_phApproxSolGrad); // should one use m_phApproxSolGrad
                //gsInfo << "m_phApproxSolGrad:\n" << m_phApproxSolGrad << "\n";

                // |grad v - y|^2 := (v_{x1} - y_1)^2 + (v_{x1} - y_2)^2
                mdSq  += weight * (m_approxSolGrad.col(k) - m_approxFluxVals.col(k)).squaredNorm();
                //mdSq  += weight * (math::pow( m_approxSolGrad(0, k) - m_approxFluxVals(0, k), 2.0 )
                //                  + math::pow( m_approxSolGrad(1, k) - m_approxFluxVals(1, k), 2.0 ));
                // |\div y + f|^2_K
                meqSq += weight * math::pow( m_approxFluxDiv(0, k) + m_rhsFunctVals(0, k), 2.0 );
            }

            majSq = (1.0 + m_beta) * mdSq + (1.0 + 1.0 / m_beta) * math::pow(m_cFriedrichs, 2) * meqSq;

            // Add into the total value of the majorant and its components
            maj += majSq;
            mD  += mdSq;
            mEq += meqSq;

            return majSq;
        }

    public:

        const std::vector< std::vector<T> > & elementNormsFullData()  const
        {
            return m_elWiseFull;
        }

    private:

        gsBoundaryConditions<T> m_bcInfo;
        T m_cFriedrichs, m_beta;

        gsMatrix<T> m_rhsFunctVals;      // f_h

        gsMatrix<T> m_approxSolGrad;    // grad(u_h)
        gsMatrix<T> m_approxFluxVals;   // y_h
        gsMatrix<T> m_approxFluxDiv;    // div(y_h)

        // TODO: check the difference between m_approxFluxVals and m_phApproxFluxVals
        gsMatrix<T> m_phApproxFluxVals;   // y_h on physical domain
        gsMatrix<T> m_phApproxSolGrad;    // grad(u_h) on physical domain
        gsMatrix<T> m_phApproxFluxDiv;    // div(y_h) on physical domain

        unsigned m_parDim;              // dimension of the parameter domain
        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsEstimator is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        using gsEstimator<T>::patchesPtr;    // pointer to all the patches
        using gsEstimator<T>::field1;

        unsigned m_storeElWiseType;                 // auxiliary flag, mainly for debugging
        std::vector< std::vector<T> > m_elWiseFull; // vector of the element-wise values

    };

} // namespace gismo
