/** @file gsVisitorLocalisedSpaceTime.h

    @brief Heat equation localised element visitor.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#pragma once

namespace gismo
{

/** \brief Visitor (localised) for the heat equation.
 *
 * Assembles the bilinear terms
 * \f[ (\partial_t u- \Delta u = f  \f]
 * For \f[ u = g \quad on  \quad \Sigma \f],
 *     \f[ u(x,0) = u_0(x) \quad on \quad \Sigma_0 \f],
 */

    template <class T, bool paramCoef = false>
    class gsVisitorLocalisedSpaceTime
    {
    public:

        /** \brief Constructor for gsVisitorLocalisedSpaceTime.
         *
         * \param[in] rhs Given right-hand-side function/source term that, for
         * \param[in] theta Given coefficient
         */
        /// Constructor with the right hand side function of the Poisson equation
        gsVisitorLocalisedSpaceTime(const gsPde<T> & pde)
        {
            const gsPoissonHeterogeneousPde<T>* pde_ptr
                    = static_cast<const gsPoissonHeterogeneousPde<T>*>(&pde);
            rhs_ptr = pde_ptr->rhs();
        }

        void initialize(const gsBasis<T> & basis,
                        const index_t patchIndex,
                        const gsOptionList & options,
                        gsQuadRule<T>    & rule,
                        unsigned         & evFlags )
        {
            // Setup Quadrature
            rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here

            // Set Geometry evaluation flags
            evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
        }

        // Evaluate on element.
        inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                             gsGeometryEvaluator<T> & geoEval,
                             gsMatrix<T> const      & quNodes)
        {

            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            basis.active_into(quNodes.col(0), actives);
            numActive = actives.rows();

            // Evaluate basis functions on element
            basis.evalAllDers_into( quNodes, 2, basisData );

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval.evaluateAt(quNodes);

            // Evaluate right-hand side at the geometry points paramCoef
            // specifies whether the right hand side function should be
            // evaluated in parametric(true) or physical (false)
            rhs_ptr->eval_into( (paramCoef ?  quNodes :  geoEval.values() ), rhsVals );
            // Evaluate the coefficient, we know that it is the expression

            // Initialize local matrix/rhs
            localMat.setZero(numActive, numActive      );
            localRhs.setZero(numActive, rhsVals.rows() );//multiple right-hand m_sides

        }

        inline void assemble(gsDomainIterator<T>    & element,
                             gsGeometryEvaluator<T> & geoEval,
                             gsVector<T> const      & quWeights)
        {
            const unsigned d = element.dim();
            gsMatrix<> ones = gsMatrix<>::Identity(d-1, 1);
            ones.setOnes();


            // basisVals =
            // B1(p1) B1(p2) .... B1(pk) ...
            // ...
            // BN(p1) BN(p2) .... BN(pk)
            const gsMatrix<T> & basisVals      = basisData[0];
            const gsMatrix<T> & basisGrads     = basisData[1];
            const gsMatrix<T> & basis2ndDerivs = basisData[2];

            gsMatrix<T> localMat_, localRhs_;
            gsMatrix<T> upwindFunc, upwindSpaceGradFunc;
            gsMatrix<T> approxSpaceLaplace;

            // look over the number of evaluation points
            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k] * geoEval.measure(k);


                // Compute physical gradients at k-th evaluation poins as a d x N matrix,
                // where d is the physical dimension and
                // N is the number of active basis-functions ar this point
                // [dx_1 B1 dx_1 B2 dx_1 B3 dx_1 B4 ]
                // ...
                // [dx_d B1 dx_d B2 dx_d B3 bx_d B4 ]
                geoEval.transformGradients(k, basisGrads, ph_basisGrad);
                // Compute the physical 2nd derivative at k-th evaluation poins as a d(d+1)/2 x N matrix
                // [d2x B1 d2y B1 dxdy B1]
                // ...
                // [d2x BN d2y BN dxdy BN]
                // TODO: to ph_basisDeriv2.transposeInPlace()
                // [d2x B1 d2x B2 dxdy B1]
                // [d2y B1 d2y
                // [dxy B1 BN d2y BN dxdy BN]
                geoEval.transformDeriv2Hgrad(k, basisGrads, basis2ndDerivs, ph_basisDeriv2);

                // Physical gradients for d-1 variables, matrix of a size (d-1) x N,
                // where N is NumActive matrix
                // [ dx_1 B1     ... dx_1 BN ]
                // ...
                // [ dx_{d-1} B1 ... dx_{d-1} BN]
                //
                ph_basisSpaceGrad = ph_basisGrad.topRows(d-1);                   // is a (Dim-1) x NumActive matrix

                // Physical gradients for dth variable (time variable), matrix of a size 1 x N,
                // where N is NumActive matrix
                // [ dt B1 ... dt BN ]

                ph_basisTimeGrad      = ph_basisGrad.bottomRows(1);
                ph_basisMixedDeriv  = ph_basisDeriv2.rightCols(d-1);
                ph_basisSpace2ndDeriv = ph_basisDeriv2.leftCols(d-1);
                ph_basisSpaceLaplace  = ph_basisSpace2ndDeriv * ones;


                //gsInfo << "ph_basisSpace2ndDeriv : \n" << ph_basisSpace2ndDeriv << "\n";
                //gsInfo << "ones : \n" << ones << "\n";
                //gsInfo << "ph_basisSpaceLaplace : \n" << ph_basisSpaceLaplace << "\n";

                T h_K     = element.getCellSize();
                T C_invK  = 1;
                T theta_K = h_K / ( d * math::pow(C_invK, 2));
                T delta_K = theta_K * h_K;

                // upwind vector = w = [wh_1 .... wh_N], where N is number of active basis function
                upwindFunc          = basisVals.col(k).transpose() + delta_K * ph_basisTimeGrad;
                // grad_x upwind vector = grad_x ([wh_1 .... wh_N])
                upwindSpaceGradFunc = ph_basisSpaceGrad + delta_K * ph_basisMixedDeriv.transpose();
                //upwindSpaceGradFunc = ph_basisSpaceGrad - delta_h * ph_basisMixedDeriv.transpose();

                //----------------------------------------------------------------------------------------------------//
                // Debug report:
                //----------------------------------------------------------------------------------------------------//
                /*
                gsInfo << "basisVals.col(k) = \n " << basisVals.col(k) << "\n";
                if (d == 2)
                    gsInfo << "ph_basisGrad_[d x N] = \n "
                           << "[dx B1 dx B2 dx B3 dx B4 ... dx BN] \n "
                           << "[dy B1 dy B2 dy B3 by B4 ... by BN] = \n"
                           << ph_basisGrad << "\n\n";
                if (d == 3)
                    gsInfo << "ph_basisGrad_[d x N] = \n "
                           << "[dx B1 dx B2 dx B3 dx B4 ... dx BN] \n "
                           << "[dy B1 dy B2 dy B3 by B4 ... by BN] \n "
                           << "[dz B1 dz B2 dz B3 bz B4 ... bz BN] = \n"
                           << ph_basisGrad << "\n\n";

                if (d == 2)
                    gsInfo << "ph_basisDeriv2_[N x d*(d-1)] = \n "
                           << "[d2x B1 d2y B1 dxdy B1] \n "
                           << "[d2x B2 d2y B2 dxdy B2] \n "
                           << "[d2x B3 d2y B3 dxdy B3] \n "
                           << "[d2x B4 d2y B4 dxdy B4] = \n"
                           << ph_basisDeriv2 << "\n\n";
                if (d == 3)
                    gsInfo << "ph_basisDeriv2_[N x d*(d-1)] = \n "
                           << "[d2x B1 d2y B1 d2z B1 dxdy B1 dxdz B1 dydz B1] \n "
                           << "[d2x B2 d2y B2 d2z B2 dxdy B2 dxdz B2 dydz B2] \n "
                           << "... \n "
                           << "[d2x BN d2y BN d2z BN dxdy BN dxdz BN dydz BN] = \n"
                           << ph_basisDeriv2 << "\n\n";
                if (d == 2)
                    gsInfo << "ph_basisSpaceGrad = \n "
                           << "[ dx B1 dx B2 dx B3 dx B4 ... dx BN] = \n "
                           << ph_basisSpaceGrad << "\n";
                else if (d == 3)
                    gsInfo << "ph_basisSpaceGrad = \n "
                           << "[dx B1 dx B2 dx B3 dx B4 ... dx BN] \n "
                           << "[dy B1 dy B2 dy B3 dy B4 ... dy BN] = \n "
                           << ph_basisSpaceGrad << "\n\n";
                if (d == 2)
                    gsInfo << "ph_basisTimeGrad = \n [ dy B1 dy B2 dy B3 by B4 ... dy BN] = \n " // is a 1 x NumActive matrix
                           << ph_basisTimeGrad << "\n\n";
                else if (d == 3)
                    gsInfo << "ph_basisTimeGrad = \n [ dz B1 dz B2 dz B3 bz B4 ... dz BN] = \n " // is a 1 x NumActive matrix
                           << ph_basisTimeGrad << "\n\n";
                if (d == 2)
                    gsInfo << "ph_basisMixedDeriv = \n [dxdy B1] \n [dxdy B2] \n [dxdy B3] \n [dxdy B4] = \n " // is a NumActive x (d - 1) matrix
                           << ph_basisMixedDeriv << "\n\n";
                else if (d == 3)
                    gsInfo << "ph_basisMixedDeriv = \n [dxdz B1 dydz B1] \n [dxdz B2 dydz B2] \n [dxdz B3 dydz B3] \n ... \n [dxdz BN dydz BN] = \n " // is a NumActive x (d - 1) matrix
                           << ph_basisMixedDeriv << "\n\n";
                gsInfo << "----------------------- \n\n";
                gsInfo << "h_K : " << h_K << "\n";
                gsInfo << "theta_K : " << theta_K << "\n";
                gsInfo << "----------------------- \n\n";
                gsInfo << "ph_basisMixedDeriv * ph_basisSpaceGrad : \n" << ph_basisMixedDeriv * ph_basisSpaceGrad << "\n\n";
                gsInfo << "rhsVals.col(k) : \n" << rhsVals.col(k) << "\n\n";
                gsInfo << "(v + delta_K * dt v)_[1 x N] : \n" << upwindFunc << "\n\n";
                gsInfo << "(grad_x (v + delta_K * dt v))_[2 x N] : \n" << upwindSpaceGradFunc << "\n\n";
                gsInfo << "(dt v * (v + delta_K * dt v))_[N x N] : \n" << upwindFunc.transpose() * ph_basisTimeGrad << "\n\n";
                gsInfo << "(grad_x v * grad_x (v + delta_K * dt v))_[N x N] : \n" << upwindSpaceGradFunc.transpose() * ph_basisSpaceGrad << "\n\n";
                gsInfo << "grad_x v * grad_x v: \n" <<  ph_basisSpaceGrad.transpose() * ph_basisSpaceGrad << "\n\n";
                gsInfo << "( - delta_K * v_t)' * Laplace_x v' : \n" << - delta_K * ph_basisTimeGrad.transpose() * ph_basisSpaceLaplace.transpose() << "\n\n";
                gsInfo << "Laplace_x v * ( - delta_K * v_t): \n" << - delta_K * ph_basisSpaceLaplace * ph_basisTimeGrad << "\n\n";
                */
                // u_t * v + grad_x u * grad_x v +
                // dt ([vh_1 ... vh_N]^T) * [vh_1 + theta * h * dt vh_1,  .... , vh_N + theta * h * dt vh_N]
                // + grad_x ([vh_1 ... vh_N]^T) * grad_x ([vh_1 + theta * h * dt vh_1,  .... , vh_N + theta * h * dt vh_N])
                // [N x 1] * [1 x N] + [N * (d-1)] * [(d-1) * N]
                localMat.noalias() += weight * ( upwindFunc.transpose() * ph_basisTimeGrad                              // (v + delta_K * v_t) * v_t
                                                 + ph_basisSpaceGrad.transpose() * ph_basisSpaceGrad                    // grad_x v * grad_x v
                                                 - delta_K * ph_basisTimeGrad.transpose() * ph_basisSpaceLaplace.transpose());      // ( - delta_K * v_t) * Laplace_x v
                // [f_1 .. f_N]^T * [vh_1 + theta * h * dt vh_1,  .... , vh_N + theta * h * dt vh_N]
                // ([1 x 1] * [N x 1])^T = [1 x N]
                localRhs.noalias() += weight * ( (rhsVals.col(k) * upwindFunc).transpose()) ;

                /*
                gsInfo << "--------------------------------------------------------------------------\n";
                gsInfo << "h : " << h << "\n";
                gsInfo << "--------------------------------------------------------------------------\n";

                gsInfo << "upwindFunc : \n" << upwindFunc << "\n";
                gsInfo << "upwindSpaceGradFunc : \n" << upwindSpaceGradFunc << "\n";


                gsInfo << "localMat increment : \n " << weight * ( upwindFunc.transpose() * ph_basisTimeGrad + ph_basisSpaceGrad.transpose() * upwindSpaceGradFunc) << "\n";
                gsInfo << "localMat : \n " << localMat << "\n";

                gsInfo << "localRhs increment : \n " << weight * ( (rhsVals.col(k) * upwindFunc).transpose()) << "\n";
                gsInfo << "localRhs : \n " << localRhs << "\n";
                gsInfo << "\n";
                */

            }
        }

        inline void localToGlobal(const int patchIndex,
                                  const std::vector<gsMatrix<T> > & eliminatedDofs,
                                  gsSparseSystem<T>     & system)
        {
            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives, patchIndex, actives);

            // Add contributions to the system matrix and right-hand side
            system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);
            /*
            gsInfo << "localMat to push: \n " << localMat << "\n";
            gsInfo << "localRhs to push: \n " << localRhs << "\n";
            gsInfo << "actives: \n " << actives << "\n";
            gsInfo << "eliminatedDofs.front(): \n " << eliminatedDofs.front() << "\n";
            gsInfo << "\n";
            */

        }

    protected:
        // Right hand side
        const gsFunction<T> * rhs_ptr;

    protected:
        // Basis values
        std::vector<gsMatrix<T> > basisData;
        gsMatrix<unsigned> actives;
        index_t numActive;

        gsMatrix<T> ph_basisGrad, ph_basisDeriv2;
        gsMatrix<T> ph_basisTimeGrad, ph_basisSpaceGrad, ph_basisMixedDeriv, ph_basisSpace2ndDeriv, ph_basisSpaceLaplace;


    protected:
        // Local values of the right hand side
        gsMatrix<T> rhsVals;
        boxSide side1;

    protected:
        // Local matrices
        gsMatrix<T> localMat, localMass;
        gsMatrix<T> localRhs;
    };

} // namespace gismo
