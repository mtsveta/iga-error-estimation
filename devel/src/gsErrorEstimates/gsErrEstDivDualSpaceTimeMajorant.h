/** @file gsErrEstDivDualSpaceTimeMajorant.h

@brief A term (measuring the residual in the dual equation) of the functional-type error estimator
 for the space-time formulation of the heat problem.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): S. Matculevich
*/

#pragma once


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
    class gsErrEstDivDualSpaceTimeMajorant : public gsNorm<T>
    {
        friend class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;

    public:

        /**
        * \brief Constructor with
        * \param _approxSol: Approximate solution.
        * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param bcInfo: Boundary conditions
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */

        gsErrEstDivDualSpaceTimeMajorant(const gsField<T> & _approxSol,
                                      const gsMultiPatch<T> & _approxFlux,
                                      const gsBoundaryConditions<T> & bcInfo,
                                      bool _rhsFunctParam = false)
                : gsNorm<T>(_approxSol, _approxFlux),
                  m_bcInfo(bcInfo),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = true; // since bc are given
            m_storeElWiseType = 0;  // auxiliary flag, mainly for debugging
        }


        /**
        * \brief Constructor with
        * \param _approxSol: Discrete solution
        * \param _approxFlux: Auxiliary function, which physically mimics fluxes of the approximate solution
        * and helps to reconstruct optimal error estimate.
        * \param _rhsFunctParam: Boolean flag indicating whether the \em _rhsFunct is parameterized
        * (if true, the evaluation points must be given on the parameter domain)
        */
        gsErrEstDivDualSpaceTimeMajorant(const gsField<T> & _approxSol,
                                      const gsMultiPatch<T> & _approxFlux,
                                      bool _rhsFunctParam = false)
                : gsNorm<T>(_approxSol,_approxFlux),
                  m_f2param(_rhsFunctParam)
        {
            m_bcInitialized = false;    // since no bc is given
            m_storeElWiseType = 0;      // auxiliary flag, mainly for debugging
        }

    public:

        /** \brief Computes the error indicator, namely
         * computes the term of the functional-type error estimate \f$\overline{\rm m}_d\f$ that measures
         * the error in dual variable, i.e.,
         * \f$\overline{\rm m}^2_d := \sum_{K \in \Tau_K} \overline{\rm m}^2_{d, K}\f$, where
         * \f$\overline{\rm m}^2_{d, K} := \int_K ( \nabla v - y )^2 \dx\f$
         *
         * \param[in] storeElWise: Boolean flag, indicating whether the element-wise contributions must be stored.
         * If <em>storeEletWise = true</em>, the gsVector of element-wise estimates \f$\overline{\rm m}^2_{d, K}\f$
         * can be accessed by calling elementNorms().
         *
         * \returns The total value of the estimate \f$ \overline{\rm m}_d \f$.
         *
         */
        T compute(bool storeElWise = false, int elemNum = 0)
        {
            m_elWiseFull.clear();   // clear the vector of element-wise values
            m_storeElWiseType = unsigned( storeElWise );    // set the bool tag for whether to store element-wise values
            this->applyElem(*this,storeElWise, elemNum);                 // m_storeElWiseType instead of storeElWise ?
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
                             gsMatrix<T>            & quNodes)
        {
            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            // Evaluate approximate solution's gradient
            approxSol.deriv_into(quNodes, approxGrad); // Evaluate approximate solution's gradient
            approxSol.deriv2_into(quNodes, approxDeriv2);

            // Evaluate approximate flux's values
            approxFlux.deriv_into(quNodes, fluxGrad);
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
                         T & maj) {
            // Access the index of the active patch
            // unsigned actPatch = geoEval.id();

            // Initialize the variable for the integral on the element
            T mdSq(0.0);

            Eigen::setNbThreads(1);

            const std::pair<int,int> dim = geoEval.geometry().function(0).dimensions();
            const int d = dim.first;

            gsMatrix<> ones = gsMatrix<>::Identity(d-1, 1);
            ones.setOnes();
            gsMatrix<T> approxSpaceLaplace;

            gsMatrix<T> ph_approxDeriv2;       // grad_x(u_h)
            gsMatrix<T> ph_fluxSpaceDiv;       // div_x (y_h) on physical domain
            gsMatrix<T> ph_fluxGrad;           // grad (y_h) on physical domain
            gsMatrix<T> ph_fluxSpaceGrad;      // grad_x (y_h) on physical domain


            gsMatrix<T> resDual, resDualSpaceComponent;
            resDualSpaceComponent.setZero(approxGrad.rows(), 1);

            // divide the execution of a loop into multiple threads and
            // execute those loop slices in parallel using SIMD
            for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
            {
                ph_fluxSpaceGrad.setZero(dim.first, dim.second);

                const T weight = quWeights[k] * geoEval.measure(k);
                // Map approxGrad from the param domain to the physical domain
                geoEval.transformDeriv2Hgrad(k, approxGrad, approxDeriv2, ph_approxDeriv2);
                geoEval.transformGradients(k, fluxGrad, ph_fluxGrad);

                // |grad_x v - y|^2 := (v_{x_1} - y_1)^2 + ... + (v_{x_{d-1}} - y_1)^2
                ph_approxSpaceLaplace = ph_approxDeriv2.leftCols(d-1);
                approxSpaceLaplace = ph_approxSpaceLaplace * ones;

                //gsInfo << "ones : \n" << ones << "\n";

                geoEval.transformGradients(k, fluxGrad, ph_fluxGrad);
                ph_fluxSpaceGrad.topRows(d-1) = ph_fluxGrad.topRows(d-1);
                ph_fluxSpaceGrad.resize(dim.first * dim.second, 1);
                convertValue<T>::derivToDiv(ph_fluxSpaceGrad, dim, ph_fluxSpaceDiv);

                // |Delta_x v - div_x y|^2 := (Delta_x v - div_x y)^2
                mdSq += weight * std::pow(approxSpaceLaplace(0, 0) - ph_fluxSpaceDiv(0, 0), 2.0);
            }

            Eigen::setNbThreads(0);
            maj += mdSq;

            return mdSq;
        }

    public:

        const std::vector< std::vector<T> > & elementNormsFullData()  const
        {
            return m_elWiseFull;
        }

    private:

        gsBoundaryConditions<T> m_bcInfo;

        gsMatrix<T> ph_approxSpaceLaplace;
        gsMatrix<T> approxGrad, approxDeriv2;   // grad(u_h)
        gsMatrix<T> fluxGrad;     // y_h

        unsigned m_parDim;              // dimension of the parameter domain

        bool m_f2param;                 // boolean flag telling whether the f2 function in the gsNorm is parametric
        bool m_bcInitialized;           // boolean flag telling whether the BC initialized

        using gsNorm<T>::patchesPtr;    // pointer to all the patches

        unsigned m_storeElWiseType;                 // auxiliary flag, mainly for debugging
        std::vector< std::vector<T> > m_elWiseFull; // vector of the element-wise values

    };

} // namespace gismo
