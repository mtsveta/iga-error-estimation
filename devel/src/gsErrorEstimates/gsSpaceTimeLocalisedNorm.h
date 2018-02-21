/** @file gsSpaceTimeNorm.h

    @brief Computes the Space-Time norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore and S. Matculevich
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
    class gsSpaceTimeNorm : public gsNorm<T>
    {
        friend  class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;

    public:

        gsSpaceTimeNorm(const gsField<T> & _field1,
                        const gsFunction<T> & _func2,
                        bool _f2param = false)
                : gsNorm<T>(_field1, _func2), thetaFunc(NULL), f2param(_f2param)
        {

        }

        gsSpaceTimeNorm(const gsField<T> & _field1,
                        const gsFunction<T> & _func2,
                        const gsFunction<T> & _thetaFunc,
                        bool _f2param = false)
                : gsNorm<T>(_field1,_func2), thetaFunc(&_thetaFunc), f2param(_f2param)
        {

        }

        gsSpaceTimeNorm(const gsField<T> & _field1)
                : gsNorm<T>(_field1), thetaFunc(NULL), f2param(false)
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
            evFlags = NEED_MEASURE | NEED_VALUE | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        void evaluate(gsGeometryEvaluator<T> & geoEval,
                      const gsFunction<T>    & func1,
                      const gsFunction<T>    & func2,
                      gsMatrix<T>            & quNodes)
        {
            // Evaluate first function
            func1.deriv_into(quNodes, f1ders);
            thetaFunc->eval_into(quNodes, thetaVals);

            // Evaluate second function (defined of physical domain)
            geoEval.evaluateAt(quNodes);

            func2.deriv_into(f2param ? quNodes : geoEval.values(), f2ders);
        }

        // assemble on element
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & accumulated)
        {
            const int d = element.dim();
            T sum(0.0);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Transform the gradients
                // f1pders : N X dim
                // f1pdersSpace, f1pdersTime : N X (dim-1)
                geoEval.transformGradients(k, f1ders, f1pders);
                Rows _f1pdersSpace = f1pders.topRows(d-1);
                Rows _f1pdersTime =  f1pders.bottomRows(1);

                Rows f2dersSpace = f2ders.topRows(d-1);
                Rows f2dersTime =  f2ders.bottomRows(1);

                T weight = quWeights[k] *  geoEval.measure(k);
                T h_K = element.getCellSize();
                T C_invK = 1;
                T theta_K = h_K / ( d * math::pow(C_invK, 2));

                T delta_K = theta_K * h_K;

                // f2ders : N X 1
                sum += weight * ( (_f1pdersSpace - f2dersSpace.col(k)).squaredNorm()
                                  + delta_K * (_f1pdersTime - f2dersTime.col(k)).squaredNorm() );

            }

            accumulated += sum;
            return sum;
        }


        inline T takeRoot(const T v) { return math::sqrt(v);}

    private:
        const gsFunction<T>* thetaFunc; // If this is NULL a numerical approximation will be used

        using gsNorm<T>::m_value;
        using gsNorm<T>::m_elWise;

        gsMatrix<T> f1ders, f2ders, thetaVals;
        gsMatrix<T> f1pders;
        gsVector<T> unormal;

        bool f2param;// not used yet
    };
} // namespace gismo


