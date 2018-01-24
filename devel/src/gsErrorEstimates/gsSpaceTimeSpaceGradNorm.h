/** @file gsSpaceTimeSpaceGrad.h

    @brief Computes the Space-Time norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

/** @brief The gsSpaceTimeSpaceGrad class provides the functionality
 * to calculate the Space-Time seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
    template <class T>
    class gsSpaceTimeSpaceGrad : public gsNorm<T>
    {
        friend  class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;

    public:

        gsSpaceTimeSpaceGrad(const gsField<T> & _field1,
                        const gsFunction<T> & _func2,
                        bool _f2param = false)
                : gsNorm<T>(_field1,_func2), dfunc2(NULL), f2param(_f2param)
        {

        }

        gsSpaceTimeSpaceGrad(const gsField<T> & _field1,
                        const gsFunction<T> & _func2,
                        const gsFunction<T> & _dfunc2,
                        bool _f2param = false)
                : gsNorm<T>(_field1,_func2), dfunc2(&_dfunc2), f2param(_f2param)
        {

        }


        gsSpaceTimeSpaceGrad(const gsField<T> & _field1)
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
            evFlags = NEED_MEASURE | NEED_VALUE | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        void evaluate(gsGeometryEvaluator<T> & geoEval,
                      const gsFunction<T>    & func1,
                      const gsFunction<T>    & _func2,
                      gsMatrix<T>            & quNodes)
        {
            // Evaluate first function
            func1.deriv_into(quNodes, f1ders);

            // Evaluate second function (defined of physical domain)
            geoEval.evaluateAt(quNodes);

            if(dfunc2==NULL)
            {
                // get the gradients to columns
                _func2.deriv_into(f2param ? quNodes : geoEval.values(), f2ders);

            }
            else
                dfunc2->eval_into(f2param ? quNodes : geoEval.values(), f2ders);

            // ** Evaluate function v
            //gsMatrix<T> f2val = func2Param ? func2.deriv(quNodes)
            //: func2.eval( geoEval->values() );
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
                Rows f2dersSpace = f2ders.topRows(d-1);

                const T weight = quWeights[k] *  geoEval.measure(k);

                sum += weight *( (_f1pdersSpace - f2dersSpace.col(k)).squaredNorm() );

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

        gsMatrix<T> f1ders, f2ders;
        gsMatrix<T> f1pders;
        gsVector<T> unormal;

        bool f2param;// not used yet
    };
} // namespace gismo
