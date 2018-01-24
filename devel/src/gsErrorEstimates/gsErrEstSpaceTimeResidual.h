/** @file gsErrEstSpaceTimeResidual.h

    @brief Residual-type error estimator for the Parabolic problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Adapted from gsErrEstPoissonResidual.h --S. Kleiss
    
    Author(s): S. Moore
*/

#pragma once


namespace gismo
{


/** \brief Provides a residual-type and element-wise error estimator
 * for the Poisson problem.
 *
 * Let the Poisson problem on the domain \f$ \Omega \f$ be given by
 * \f[ \partial_t u -\Delta u = f,\quad u = g \mathrm{\ on\ } \Sigma,
 *  \quad u_0(x) = 0 \mathrm{\ on\ } \Sigma_0, \f]
 * where\n
 * \f$ f   \f$ is a given right-hand-side,\n
 * \f$ g   \f$ is given Dirichlet data, \f$ \Sigma \f$ is the Dirichlet boundary,\n
 * \f$ u_0 \f$ is given initial data, \f$ \Sigma_0 \f$ is the initial COndition.
 *
 * The error estimate \f$\eta\f$ for a computed discrete solution
 * \f$ u_h \f$ is given by
 * \f[ \eta^2 = \sum_K \eta_K^2 \f]
 * where the local estimate \f$ \eta_K \f$ on an element \f$ K \f$ is given by
     \f[
     \eta_K^2 =
      \int_K (f -\partial_t u_h + \Delta u_h)^2 dx
     + \int_{\partial K \cap \Sigma_0} ( u_h - u_0 )^2 ds \f]
 * \f$ h \f$ denotes the size of element \f$ K \f$,\n.
 *
 * \ingroup Assembler
 */
template <class T>
class gsErrEstSpaceTimeResidual : public gsNorm<T>
{
    friend class gsNorm<T>;
    typedef typename gsMatrix<T>::RowsBlockXpr Rows;

public:

    // f1 in gsNorm corresponds to discrete Solution
    // f2 in gsNorm corresponds to right-hand-side


    /**
     * \brief Constructor
     * \param _discSolution Discrete solution
     * \param _rhsFunction Right-hand-side-/Source-function \f$ f \f$
     * of the Poisson problem.
     * \param bcInfo Boundary conditions
     * \param _rhsFunctionParam Flag indicating whether the \em _rhsFunction
     * is parameterized (if true, the evaluation points must be given
     * on the parameter domain)
     *
     */
    gsErrEstSpaceTimeResidual(const gsField<T> & _discSolution,
                            const gsFunction<T> & _rhsFunction,
                            const gsBoundaryConditions<T> & bcInfo,
                            bool _rhsFunctionParam = false)
    : gsNorm<T>(_discSolution,_rhsFunction), m_bcInfo(bcInfo),  m_f2param(_rhsFunctionParam)
    {
        m_bcInitialized = true;

        // auxiliary flag, mainly for debugging
        m_storeElWiseType = 0;
    }


     /** \brief Constructor
     * \param _discSolution Discrete solution
     * \param _rhsFunction Right-hand-side-/Source-function \f$ f\f$
     * of the Poisson problem.
     * \param _rhsFunctionParam Flag indicating whether the \em _rhsFunction
     * is parameterized (in this case, the evaluation points must be given
     * on the parameter domain
    */
    gsErrEstSpaceTimeResidual(const gsField<T> & _discSolution,
             const gsFunction<T> & _rhsFunction,
             bool _rhsFunctionParam = false)
    : gsNorm<T>(_discSolution,_rhsFunction), m_f2param(_rhsFunctionParam)
    {
        m_bcInitialized = false;

        // auxiliary flag, mainly for debugging
        m_storeElWiseType = 0;
    }     
    
public:

    /** \brief Computes the error estimate.
     *
     * Computes the residual-based error estimate \f$\eta\f$
     * (see class-documentation at the top).
     *
     *
     * \param storeElWise Bool indicating whether the element-wise
     * errors should be stored also. If <em>storeEletWise = true</em>,
     * the gsVector of element-wise estimates \f$\eta_K^2\f$
     * can be obtained by
     * calling elementNorms().
     *
     * \returns The total estimated error \f$ \eta \f$.
     *
     */
    T compute(bool storeElWise = true)
    {
        m_elWiseFull.clear();
        m_storeElWiseType = unsigned( storeElWise );

        this->apply(*this,storeElWise);
        return this->m_value;
    }

    // like compute, but allowing for m_storeElWiseType == 2
    // in which case more information on the local error estimates is stored
    T compute2(unsigned storeType = 1)
    {
        bool storeElWise = ( storeType == 0 ? false : true );

        m_elWiseFull.clear();
        m_storeElWiseType = storeType;

        this->apply(*this,storeElWise);
        return this->m_value;
    }

protected:

    /**
     * @brief Initializes the error estimator
     *
     * Sets up the quadrature rule (based on the degree of \em basis) and
     * the \em evFlags for the gsGeometryEvaluator that are needed
     * for this specific problem.
     *
     * \param[in] basis
     * \param[out] rule
     * \param[out] evFlags
     */
    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        m_parDim = basis.dim();

        GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Called error estimator with dimension other than 2 or 3.");

        // Setup Quadrature
        gsVector<index_t> numQuadNodes( m_parDim );
        for (unsigned i = 0; i < m_parDim; ++i)
            numQuadNodes[i] = 2*basis.degree(i) + 1;

        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        // is used in evaluate()
        evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN | NEED_2ND_DER | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    /**
     * @brief Evaluate some needed data on the given quadrature nodes
     *
     * Executes and stores needed function evaluations at \em quNodes.\n
     * The gsGeometryEvaluator \em geoEval is also evaluated at the nodes,
     * using evaluation flags specified in initialize().
     *
     * @param[in,out] geoEval
     * @param[in] discSolution
     * @param[in] rhsFunction
     * @param[in] quNodes
     */
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsFunction<T>    & discSolution,
                         const gsFunction<T>    & rhsFunction,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate discrete solution
        discSolution.deriv_into(quNodes, m_discSolDer);
        discSolution.deriv2_into(quNodes, m_discSol2ndDer);        

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate right-hand-side function (defined of physical domain)
        rhsFunction.eval_into(geoEval.values(), m_rhsFctVals);
    }

    // assemble on element

    /**
     * @brief Computes the local error estimate on an element.
     *
     * See documentation of the class for the computed error estimate.
     *
     * @param[in] element specifies the element \f$ K \f$.
     * @param[in] geoEval gsGeometryEvaluator as evaluated in evaluate().
     * @param[in] quWeights Quadrature weights \em before transformation the the
     * element, i.e., the sum of the weights should be 1.
     * @return The \em squared estimate \f$ \eta_K^2 \f$ of local error on
     * element \f$ K \f$.
     */
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {

        //  const int d = element.dim();
        T sumVolSq(0.0);
        //geoEval.transformGradients(k, basisGrads , m_phdiscSolDer);
        // will be used to set up quadrature points on the side
        for (index_t k = 0; k < quWeights.size(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval.measure(k);
            // Compute the APPROXIMATION of the transformation of the
            // Laplacian and Gradient to the physical domain.
            geoEval.transformGradients(k, m_discSolDer , m_phdiscSolDer);
            geoEval.transformLaplaceHgrad(k, m_discSolDer, m_discSol2ndDer , m_phLaplace);
            geoEval.transformDeriv2Hgrad(k, m_discSolDer, m_discSol2ndDer, m_phBasisLaplace);
            m_phBasisLaplace.transposeInPlace();
            
            // residual squared: Laplace of solution + RHS.
            sumVolSq += weight * ( m_phLaplace(0,0) + m_rhsFctVals(0,k) - m_phdiscSolDer.bottomRows(1)(0,0) ) \
                               * ( m_phLaplace(0,0) + m_rhsFctVals(0,k) - m_phdiscSolDer.bottomRows(1)(0,0) );
                    
        }

        accumulated += sumVolSq;
        return sumVolSq;
    }

    inline T takeRoot(const T v) { return math::sqrt(v);}

private:

    gsBoundaryConditions<T> m_bcInfo;
    
    gsMatrix<T> m_discSol2ndDer, m_discSolDer;;
    gsMatrix<T> m_rhsFctVals;
    gsMatrix<T> m_phHessVals,m_phdiscSolDer, m_phBasisLaplace,  m_phLaplace;
    gsMatrix<T> m_discSolDerTime,m_discSol2ndDerTime ;
    gsVector<T> unormal;
    unsigned m_parDim;

    bool m_f2param;

    bool m_bcInitialized;

    // auxiliary flag, mainly for debugging
    unsigned m_storeElWiseType;
    std::vector< std::vector<T> > m_elWiseFull;
    

};

} // namespace gismo
