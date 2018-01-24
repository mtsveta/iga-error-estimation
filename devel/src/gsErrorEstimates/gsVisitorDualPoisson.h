/** @file gsVisitorDualPoisson.h

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
    class gsVisitorDualPoisson
    {
    public:

        /** \brief Constructor for gsVisitorDualPoisson.
         */
        explicit gsVisitorDualPoisson(const gsPde<T> & pde)
        {
            pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde);
            numActive = 0;
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


            // Evaluate basis functions on element
            basis.evalAllDers_into( quNodes, 1, basisData);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval.evaluateAt(quNodes);// is this generic ??

            // Evaluate gradient of approximate solution v
            //rhs_ptr->deriv_into( (paramCoef ?  quNodes :  geoEval.values() ), m_rhsGradV );
            rhs_ptr->deriv_into( quNodes, m_rhsGradV );

            // Initialize local matrix/rhs for vector valued basis function
            localMat.setZero(d * numActive, d * numActive);
            localRhs.setZero(d * numActive, 1); //multiple right-hand sides
           }

        inline void assemble(gsDomainIterator<T>    & element,
                             gsGeometryEvaluator<T> & geoEval,
                             gsVector<T> const      & quWeights) {
            // Get basis values from the basisData
            gsMatrix<T> bVals = basisData[0];

            int d = geoEval.geometry().targetDim();
            numActive = bVals.rows(); // number of the active basis function on the element

            // [! Create local matrix]
            // Create tmp matrix to stole actual mass matrix of the existing basis
            gsMatrix<T> localMatrix(numActive, numActive);
            localMatrix = bVals * quWeights.asDiagonal() * geoEval.measures().asDiagonal() * bVals.transpose();

            // Create the local matrix that corresponds to mass matrix in d dimensions
            localMat.noalias() = localMatrix.blockDiag(d);
            // [! Create local matrix]

            // [! Create local rhs]
            gsMatrix<T> locRhs(d * numActive, 1);
            locRhs.setZero(d * numActive, 1);

            // Initialize th k-th basis-funtion
            //#pragma omp parallel num_threads(4)
            {
                gsMatrix<T> m_phRhsGradV;
                gsMatrix<T> Bk(numActive, 1);
                gsMatrix<T> dDimBk(d * numActive, d);

                gsMatrix<T> localRhs_;
                localRhs_.setZero(d * numActive, 1);

                //#pragma omp for nowait //fill result_private in parallel
                for (index_t k = 0; k < quWeights.rows(); ++k) { // loop over quadrature nodes

                    // [d x 1] vector function
                    Bk = bVals.col(k);
                    dDimBk = Bk.blockDiag(d);

                    // Multiply weight by the geometry measure
                    const T weight = quWeights[k] * geoEval.measure(k);

                    // Transform the gradients of v from the parametric domain to the physical domain
                    geoEval.transformGradients(k, m_rhsGradV, m_phRhsGradV);

                    // We should multiply here by the vector RHS [d * numVector, 1]
                    localRhs_ += weight * dDimBk * m_phRhsGradV;
                }
                //#pragma omp critical
                {
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

            localMatToPush = localMat.block(0, 0, numActive, numActive); // it's the same for the lower-right block

            for (index_t i = 0; i < d; ++i) {
                localRhsToPush = localRhs.middleRows(i * numActive, numActive);
                system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), i, i);
            }

            /*
            switch (d) {
                case 2:
                    #pragma omp parallel sections
                    {
                        #pragma omp section
                        {
                            localRhsToPush = localRhs.middleRows(0, numActive);
                            system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), 0, 0);
                        }
                        #pragma omp section
                        {
                            localRhsToPush = localRhs.middleRows(numActive, numActive);
                            system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), 1, 1);
                        }
                    }
                    break;
                case 3:
                    #pragma omp parallel sections
                    {
                        #pragma omp section
                        {
                            localRhsToPush = localRhs.middleRows(0, numActive);
                            system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), 0, 0);
                        }
                        #pragma omp section
                        {
                            localRhsToPush = localRhs.middleRows(1 * numActive, numActive);
                            system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), 1, 1);
                        }
                        #pragma omp section
                        {
                            localRhsToPush = localRhs.middleRows(2 * numActive, numActive);
                            system.push(localMatToPush, localRhsToPush, actives, eliminatedDofs.front(), 2, 2);
                        }
                    }
                    break;
             }
             */

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

