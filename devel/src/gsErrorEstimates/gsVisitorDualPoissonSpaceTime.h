/** @file gsVisitorDualPoissonSpaceTime.h

    @brief Visitor for the equation

    y = grad v, grad v \in L^2(\Omega, R^d) (*)

    where y \in H(Omega, div) := \{ y \in L^2(\Omega, R^d) :  div y \in L^2(\Omega)\}.

    We assume that y \in Y_div \subset H(Omega, div) and Y^M_div := span \{ phi_1, ..., phi_M\},
    where phi_i are vector functions in H(\Omega, div). Function y is approximated as

    y = \Sum_{i = 1}^{M} Yh_i phi_i, Yh \in R^M is the vector of unknown DOFs.

    We test (*) with phi_j and obtain local contributions

    \int_{\Omega} phi_i * phi_j \dx and

    \int_{\Omega} grad v * phi_j \dx.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#pragma once

#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

/** \brief Visitor for the DivDiv equation.
 *
 *
 */

    template <class T, bool paramCoef = false>
    class gsVisitorDualPoissonSpaceTime
    {
    public:

        /** \brief Constructor for gsVisitorDualPoisson.
         */
        gsVisitorDualPoissonSpaceTime(const gsPde<T> & pde)
        {
            pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde);
        }

        void initialize(const gsBasis<T> & basis,
                        const index_t patchIndex,
                        const gsOptionList & options,
                        gsQuadRule<T>    & rule,
                        unsigned         & evFlags )
        {
            // Grab right-hand side for current patch
            rhs_ptr = &pde_ptr->rhs()->piece(patchIndex);

            // Setup Quadrature
            rule = gsGaussRule<T>(basis, options);// harmless slicing occurs here

            // Set Geometry evaluation flags
            evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        inline void evaluate(gsBasis<T> const       & basis,
                             gsGeometryEvaluator<T> & geoEval,
                             gsMatrix<T> const      & quNodes)
        {

            // Get the dimension of parametric domain
            int d = basis.dim();

            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            basis.active_into(quNodes.col(0), actives);
            numActive = actives.rows();

            //gsInfo << "actives = " << actives << "\n";

            // Evaluate basis functions on element
            basis.evalAllDers_into( quNodes, 1, basisData);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval.evaluateAt(quNodes);// is this generic ??

            // Evaluate gradient of approximate solution v
            //rhs_ptr->deriv_into( (paramCoef ?  quNodes :  geoEval.values() ), m_rhsGradV );
            rhs_ptr->deriv_into( quNodes, m_rhsGradV );

            // Initialize local matrix/rhs for vector valued basis function
            localMat.setZero(d * numActive, d * numActive);
            localRhs.setZero(d * numActive, 1); //multiple right-hand m_sides
        }

        inline void assemble(gsDomainIterator<T>    & element,
                             gsGeometryEvaluator<T> & geoEval,
                             gsVector<T> const      & quWeights) {

            gsVector<T> elementCenter = element.centerPoint();
            //gsInfo << "MM: center elem. = " << elementCenter << "\n";

            // Get basis values from the basisData
            gsMatrix<T> bVals = basisData[0];

            int d = geoEval.geometry().targetDim();
            numActive = bVals.rows(); // number of the active basis function on the element

            // [! Create local matrix]
            // Create tmp matrix to stole actual mass matrix of the existing basis
            gsMatrix<T> localMatrix_(numActive, numActive);
            //gsInfo << "bVals :\n" << bVals << "\n";
            //gsInfo << "quWeights :\n" << quWeights << "\n";
            localMatrix_ = bVals * quWeights.asDiagonal() * geoEval.measures().asDiagonal() * bVals.transpose();
            // Create the local matrix that corresponds to mass matrix in d dimensions
            localMat.noalias() = localMatrix_.blockDiag(d);

            //gsInfo << "localMatrix :\n" << localMatrix_ << "\n";
            //gsInfo << "localMat :\n" << localMat << "\n";

            // [! Create local matrix]

            // [! Create local rhs]
            gsMatrix<T> locRhs(d * numActive, 1);
            locRhs.setZero(d * numActive, 1);

            // Initialize th k-th basis-funtion
            //#pragma omp parallel num_threads(4)
            {
                gsMatrix<T> m_phRhsGradV;
                gsMatrix<T> m_phRhsSpaceGradV;
                gsMatrix<T> Bk(numActive, 1);
                gsMatrix<T> dDimBk(d * numActive, d);

                gsMatrix<T> localRhs_;
                localRhs_.setZero(d * numActive, 1);

                //#pragma omp for nowait //fill result_private in parallel
                for (index_t k = 0; k < quWeights.rows(); ++k) { // loop over quadrature nodes

                    m_phRhsSpaceGradV.resize(d, 1);
                    m_phRhsSpaceGradV.setZero(d, 1);

                    // [d x 1] vector function
                    Bk = bVals.col(k);
                    dDimBk = Bk.blockDiag(d);

                    // Multiply weight by the geometry measure
                    const T weight = quWeights[k] * geoEval.measure(k);

                    // Transform the gradients of v from the parametric domain to the physical domain
                    geoEval.transformGradients(k, m_rhsGradV, m_phRhsGradV);

                    //gsInfo << "Bk :\n" << Bk << "\n";
                    /*
                     gsInfo << "dDimBk :\n" << dDimBk << "\n";
                     gsInfo << "m_phRhsGradV = \n" << m_phRhsGradV << "\n";
                     gsInfo << "physSpaceGrad (before assiging) = \n" << m_phRhsSpaceGradV << "\n";
                     */
                    m_phRhsSpaceGradV.topRows(d-1) = m_phRhsGradV.topRows(d-1);
                    //gsInfo << "physSpaceGrad (after assiging) = \n" << m_phRhsSpaceGradV << "\n";

                    // We should multiply here by the vector RHS [d * numVector, 1]
                    localRhs_ += weight * dDimBk * m_phRhsSpaceGradV;
                    //localRhs_ += weight * dMinus1DimBk * m_phRhsSpaceGradV;
                    //gsInfo << "localRhs_ :\n" << localRhs_ << "\n";


                }
                //#pragma omp critical
                {
                    //gsInfo << "creating localRhs_ :\n";
                    localRhs.noalias() += localRhs_;
                }
            }
            // [! Create local rhs]
        }

        inline void localToGlobal(const int patchIndex,
                                  const std::vector<gsMatrix<T> > & eliminatedDofs,
                                  gsSparseSystem<T>     & system)
        {

            // Add contributions to the system matrix and right-hand side
            int d = localMat.rows() / numActive;

            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives, patchIndex, actives);

            gsMatrix <T> localMatToPush(numActive, numActive);
            gsMatrix <T> localRhsToPush(numActive, 1);

            // since localMat is block-diagonal matrix of d localMatToPush blocks
            // we can just fetch it from localMat once
            localMatToPush = localMat.block(0, 0, numActive, numActive);
            //gsInfo << "localMatToPush :\n" << localRhsToPush << "\n";

            //gsInfo << "localToGlobal function :\n";

            for (index_t i = 0; i < d; ++i) {
                /*
                if (i == d - 1) {
                    // at i == d-1 (time direction) the localPHS is zero (due to the earlier grad_x v multiplication)
                    // no need to push only zeros
                    //gsInfo << "pushing just localMatToPush :\n";
                    system.pushToMatrix(localMatToPush, actives, eliminatedDofs.front(), i, i);
                }
                else {
                    //gsInfo << "i = " << i << "\n";
                    //gsInfo << "actives = " << actives << "\n";
                    //gsInfo << "numActive = " << numActive << "\n";
                    localRhsToPush = localRhs.middleRows(i * numActive, numActive);
                    // push localMatrix and updated localRhs
                    //gsInfo << "pushing localMatToPush and localRhsToPush:\n";
                    system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), i, i);
                }
                */
                ///*
                // there is no point in pushing localMatToPush in d-th block
                // since det(localMatToPush) = 0, it won't effect the flux approximation anyways
                localRhsToPush = localRhs.middleRows(i * numActive, numActive);
                system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), i, i);
                //*/
            }
        }

    protected:
        // Pointer to the pde data
        const gsPoissonPde<T> * pde_ptr;

    protected:
        // Basis values
        std::vector<gsMatrix<T> > basisData;
        gsMatrix<T>        physGrad;
        gsMatrix<unsigned> actives;
        index_t numActive;

    protected:
        // Right hand side ptr for current patch
        const gsFunction<T> * rhs_ptr;

        // Local values of the right hand side
        gsMatrix<T> m_rhsGradV;

    protected:
        // Local matrices
        gsMatrix<T> localMat;
        gsMatrix<T> localRhs;
    };


} // namespace gismo

