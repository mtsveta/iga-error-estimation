/** @file gsErrEstEquilMajorant.h

@brief Functional-type error estimator for the Poisson problem.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): S. Matculevich
*/

#pragma once

#include <gsCore/gsEvalUtils.hpp>

namespace gismo
{

/** \brief Provides computation of equilibrium term of fucntional-type element-wise error indicator for the Poisson problem.
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
* The error indicator \f$\overline{\rm m}^2_eq\f$ for a computed discrete solution \f$ v \f$ is given by
* \f[ \overline{\rm m}^2_eq := \sum_K \overline{\rm eq}^2_{d, K} \f]
* where the local estimate \f$ \overline{\rm eq}^2_{eq, K} \f$ on an element \f$ K \f$ is given by
 \f[
 \overline{\rm m}^2_{eq, K} := \int_K ( \div y + f )^2 dx
 \f]
* \f$ f \f$ denotes the RHS of the Poisson problem and
* \f$ y \f$ denotes the approximate flux  on the element \f$ \K \f$.
*
* \ingroup Assembler
*/
    template <class T>
    class gsErrEstEquilMajorant : public gsNorm<T>
    {
        friend class gsNorm<T>;

    public:

        /**
        * \brief Constructor with
        * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _rhsFunction RHS-Source-function \f$ f \f$ of the Poisson problem.
        * \param bcInfo: Boundary conditions
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */
        gsErrEstEquilMajorant(const gsField<T> & _approxFlux,
                              const gsFunction<T> & _rhsFunction,
                              const gsBoundaryConditions<T> & bcInfo,
                              bool _rhsFunctParam = false)
                : gsNorm<T>(_approxFlux, _rhsFunction),
                  m_bcInfo(bcInfo),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = true; // since bc are given
            m_storeElWiseType = 0;  // auxiliary flag, mainly for debugging
            m_parDim = 0;           // default initilization
            m_f2param = false;      // default initilization
        }


        /**
        * \brief Constructor with
         * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _rhsFunction RHS-Source-function \f$ f \f$ of the Poisson problem.
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */
        gsErrEstEquilMajorant(const gsField<T> & _approxFlux,
                              const gsFunction<T> & _rhsFunction,
                              bool _rhsFunctParam = false)
                : gsNorm<T>(_approxFlux,_rhsFunction),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = false;    // since no bc is given
            m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
            m_parDim = 0;           // default initilization
            m_f2param = false;      // default initilization
        }

    public:

        /** \brief Computes the error indicator, namely
         * computes the term of the functional-type error estimate \f$\overline{\rm m}_eq\f$ that measures
         * the error in balance equation, i.e.,
         * \f$\overline{\rm m}^2_eq := \sum_{K \in \Tau_K} \overline{\rm m}^2_{eq, K}\f$, where
         * \f$\overline{\rm m}^2_{eq, K} := \int_K ( \div y +f )^2 \dx\f$
         *
         * \param[in] storeElWise: Boolean flag, indicating whether the element-wise contributions must be stored.
         * If <em>storeEletWise = true</em>, the gsVector of element-wise estimates \f$\overline{\rm m}^2_{eq, K}\f$
         * can be accessed by calling elementNorms().
         *
         * \returns The total value of the estimate \f$ \overline{\rm m}_eq \f$.
         *
         */
        T compute(bool storeElWise, int elemNum)
        {
            m_elWiseFull.clear();
            m_storeElWiseType = unsigned( storeElWise );

            this->applyElem(*this,storeElWise, elemNum);
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
            GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Called error estimator with dimension other than 2 or 3.");
            // whey it cannot be 1?

            // Initialize a number of quadrature nodes for each of the sub-dimenisions based on the
            // degree of a given basis
            gsVector<index_t> numQuadNodes( m_parDim );
            for (int i = 0; i < m_parDim; ++i)
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
                             const gsFunction<T>    & approxFlux,
                             const gsFunction<T>    & rhsFunction,
                             gsMatrix<T>            & quNodes)
        {
            // Evaluate approximate flux's values
            approxFlux.deriv_into(quNodes, m_approxFluxGrad);
            // Compute geometry related values
            geoEval.evaluateAt(quNodes);
            // Evaluate right-hand-side function (defined of physical domain)
            rhsFunction.eval_into(geoEval.values(), m_rhsFunctVals);

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
                         T & maj)
        {
            // Access the index of the active patch
            // unsigned actPatch = geoEval.id();

            // Initialize the variable for the integral on the element and it's facets
            T meqSq(0.0);
            Eigen::setNbThreads(1);
            const std::pair<int,int> dim = geoEval.geometry().function(0).dimensions();

            //#pragma omp num_threads(4)
            {
                gsMatrix<T> m_phApproxFluxDiv;    // div(y_h) on physical domain
                gsMatrix<T> m_phApproxFluxGrad;

                m_phApproxFluxGrad.setZero(dim.first, dim.second);
                m_phApproxFluxDiv.setZero(1, 1);

                //#pragma omp for reduction(+:meqSq)
                for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
                {
                    const T weight = quWeights[k] * geoEval.measure(k);

                    // Transform the gradients from the parametric to physical domain
                    geoEval.transformGradients(k, m_approxFluxGrad, m_phApproxFluxGrad);

                    // Resize the the obtained m_phApproxFluxGrad:
                    // [G^(-1) * dx B^(1)_k, G^(-1) * dx B^(2)_k;
                    //  G^(-1) * dy B^(1)_k, G^(-1) * dy B^(2)_k]
                    //
                    // to the shape
                    //
                    // [G^(-1) * dx B^(1)_k, G^(-1) * dy B^(1)_k, G^(-1) * dx B^(2)_k, G^(-1) * dy B^(2)_k]^T
                    //
                    m_phApproxFluxGrad.resize(dim.first * dim.second, 1);
                    // phApproxFluxDiv = G^(-1) * dx B^(1)_k + G^(-1) * dy B^(2)_k
                    convertValue<T>::derivToDiv(m_phApproxFluxGrad, dim, m_phApproxFluxDiv);

                    // |\div y + f|^2_K
                    meqSq += weight * math::pow(m_phApproxFluxDiv(0, 0) + m_rhsFunctVals(0, k), 2.0);
                }
            }
            Eigen::setNbThreads(0);

            maj += meqSq;
            return meqSq;
        }

    public:

        const std::vector< std::vector<T> > & elementNormsFullData()  const
        {
            return m_elWiseFull;
        }

    private:

        gsBoundaryConditions<T> m_bcInfo;

        gsMatrix<T> m_approxFluxGrad;    // grad(y_h)
        gsMatrix<T> m_rhsFunctVals;      // f_h

        int m_parDim;              // dimension

        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsNorm is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        using gsNorm<T>::patchesPtr;    // pointer to all the patches
        
        // auxiliary flag, mainly for debugging
        unsigned m_storeElWiseType;
        std::vector< std::vector<T> > m_elWiseFull;

    };

} // namespace gismo
