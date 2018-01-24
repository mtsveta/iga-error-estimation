/** @file gsSpaceTimeSolOperNorm.h

    @brief Computes the Space-Time norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore
*/

#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

/** @brief The gsSpaceTimeNorm class provides the functionality
 * to calculate the Space-Time norm between a field and a function.
 *
 * \ingroup Assembler
*/
    template <class T>
    class gsSpaceTimeSolOperNorm : public gsNorm<T>
    {
        friend  class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;

    public:

        gsSpaceTimeSolOperNorm(const gsField<T> & _field1,
                        const gsFunction<T> & _func2,
                        bool _f2param = false)
                : gsNorm<T>(_field1,_func2), dfunc2(NULL), f2param(_f2param)
        {

        }

        gsSpaceTimeSolOperNorm(const gsField<T> & _field1,
                        const gsFunction<T> & _func2,
                        const gsFunction<T> & _dfunc2,
                        bool _f2param = false)
                : gsNorm<T>(_field1,_func2), dfunc2(&_dfunc2), f2param(_f2param)
        {

        }

        gsSpaceTimeSolOperNorm(const gsField<T> & _field1)
                : gsNorm<T>(_field1), dfunc2(NULL), f2param(false)
        { }

    public:

        T compute(bool storeElWise = false)
        {
            this->apply(*this,storeElWise);
            return m_value;
        }


    protected:

        void initialize(const gsBasis<T> & basis,
                        gsQuadRule<T> & rule,
                        unsigned      & evFlags) // replace with geoEval ?
        {
            // Setup Quadrature
            const unsigned d = basis.dim();
            gsVector<index_t> numQuadNodes( d );
            for (unsigned i = 0; i < d; ++i)
                numQuadNodes[i] = basis.degree(i) + 1;

            // Setup Quadrature
            rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

            // Set Geometry evaluation flags
            evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN | NEED_2ND_DER | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        void evaluate(gsGeometryEvaluator<T> & geoEval,
                      const gsFunction<T>    & func1,
                      const gsFunction<T>    & _func2,
                      gsMatrix<T>            & quNodes)
        {
            //#pragma omp parallel sections
            {

//#pragma omp section
                {
                    // Evaluate first function
                    func1.deriv_into(quNodes, func1Derivs);
                    func1.deriv2_into(quNodes, func1Deriv2s);
                    // get the gradients to columns
                }
//#pragma omp section
                {
                    // Evaluate second function (defined of physical domain)
                    geoEval.evaluateAt(quNodes);

                    if (dfunc2 == NULL) {
                        // get the gradients to columns
                        _func2.deriv_into(geoEval.values(), func2Derivs);
                        _func2.deriv2_into(geoEval.values(), func2Deriv2s);
                    } else {
                        dfunc2->eval_into(geoEval.values(), func2Derivs);
                        dfunc2->deriv_into(geoEval.values(), func2Deriv2s);
                    }
                }
            }

        }

        // assemble on element
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & accumulated)
        {
            T sum(0.0);
            const unsigned d = element.dim();
            gsMatrix<> ones = gsMatrix<>::Identity(d-1, 1);
            ones.setOnes();
            gsMatrix<T> func1SpaceLaplace, func2SpaceLaplace;

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Transform the gradients
                // f1pders : N X dim
                // f1pdersSpace, f1pdersTime : N X (dim-1)

                //#pragma omp parallel sections
                {
                    //#pragma omp section
                    {
                        geoEval.transformGradients(k, func1Derivs, ph_func1Derivs);
                    }
                    //#pragma omp section
                    {
                        geoEval.transformDeriv2Hgrad(k, func1Derivs, func1Deriv2s, ph_func1Deriv2s);
                        ph_func1SpaceLaplace = ph_func1Deriv2s.leftCols(d-1);
                        ph_func2SpaceLaplace = func2Deriv2s.topRows(d-1).col(k);
                        func1SpaceLaplace = ph_func1SpaceLaplace * ones;
                        func2SpaceLaplace = ph_func2SpaceLaplace.transpose() * ones;
                    }
                }
                Rows ph_func1TimeGrad = ph_func1Derivs.bottomRows(1);
                Rows ph_func2TimeGrad  = func2Derivs.bottomRows(1);

                //gsInfo << "ph_func1Derivs \n =" << ph_func1Derivs << "\n";
                //gsInfo << "ph_func1Deriv2s \n =" << ph_func1Deriv2s << "\n";
                //gsInfo << "func2Deriv2s.topRows(d-1) \n =" << func2Deriv2s.topRows(d-1) << "\n";

                //gsInfo << "ph_func1SpaceLaplace \n =" << ph_func1SpaceLaplace << "\n";
                //gsInfo << "ph_func2SpaceLaplace \n =" << ph_func2SpaceLaplace << "\n";
                //gsInfo << "ones \n =" << ones << "\n";

                //gsInfo << "ph_func1TimeGrad \n =" << ph_func1TimeGrad << "\n";
                //gsInfo << "ph_func2TimeGrad.col(k) \n =" << ph_func2TimeGrad.col(k) << "\n";

                //gsInfo << "f1SpaceLaplace \n =" << func1SpaceLaplace << "\n";
                //gsInfo << "f2SpaceLaplace \n =" << func2SpaceLaplace << "\n";

                const T weight = quWeights[k] *  geoEval.measure(k);

                sum += weight * ( (ph_func1TimeGrad - ph_func2TimeGrad.col(k)).squaredNorm()
                                  + (func1SpaceLaplace - func2SpaceLaplace).squaredNorm() );

            }

            accumulated += sum;
            return sum;
        }


        inline T takeRoot(const T v) { return math::sqrt(v);}

    private:
        // first derivative of func2:
        const gsFunction<T> * dfunc2; // If this is NULL a numerical approximation will be used

        using gsNorm<T>::m_value;
        using gsNorm<T>::m_elWise;

        gsMatrix<T> func1Derivs, func1Deriv2s, func2Derivs, func2Deriv2s;
        gsMatrix<T> ph_func1Derivs, ph_func1Deriv2s, ph_func1SpaceLaplace,
                ph_func2Derivs, ph_func2Deriv2s, ph_func2SpaceLaplace;
        gsVector<T> unormal;

        bool f2param;// not used yet
    };
} // namespace gismo
