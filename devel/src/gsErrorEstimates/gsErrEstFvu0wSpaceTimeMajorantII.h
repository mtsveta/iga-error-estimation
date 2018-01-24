#include<gsAssembler/gsNorm.h>

#pragma once

#include <gsErrorEstimates/gsEstimator.h>

namespace gismo
{

    template <class T>
    class gsErrEstFvu0wSpaceTimeMajorantII : public gsEstimator<T>
    {
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;
        friend class gsEstimator<T>;

    public:



        gsErrEstFvu0wSpaceTimeMajorantII(const gsField<T> & _approxV,
                                       const gsField<T> & _approxW,
                                       const gsFunction<T> & _u0Function,
                                      const std::vector<boundaryInterface>& SliceInterfaces,
                                      const std::vector<patchSide>& TopSides,
                                      bool _u0FunctParam = false)
                : gsEstimator<T>(_approxV,_approxW, _u0Function), m_sliceInterfaces(SliceInterfaces), m_sides(TopSides), m_f2param(_u0FunctParam)
        {
        }

    public:

        T compute(bool storeElWise = false)
        {

            this->m_value = 0.0;

            /*
             * m_elWise.clear();   // clear the vector of element-wise values

            for ( std::vector<patchSide>::const_iterator it=m_sides.begin();it!=m_sides.end();++it) // *it ---> interface
            {
                this->apply3Func(*this, storeElWise, it->patch, it->side());
            }
            */

            side = * m_sides.begin(); // get the side from the vector of m_sides

            /*
            for (typename gsMultiPatch<T>::const_biterator bit = patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit) {
                if (bit->side() == side) {
                    gsInfo << "side : " << side << "\n";
                    side = bit->side();
                    this->apply3Func(*this, storeElWise, bit->patch, side);
                }

            }
            */
            for (typename gsMultiPatch<T>::const_biterator bit = patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit)
                if (bit->side() == side)
                    this->apply3Func(*this, storeElWise, bit->patch, bit->side());

            return this->m_value;
        }

        // we replace takeRoot function with the function that just returns the value
        inline T takeRoot(const T v) { return v;}


    protected:

        void initialize(const gsBasis<T> & basis,
                        gsQuadRule<T> & rule,
                        unsigned      & evFlags) // replace with geoEval ?
        {
            // Setup Quadrature
            const unsigned d = basis.dim();
            const int dir = side.direction();
            gsVector<index_t> numQuadNodes( d );
            for (unsigned i = 0; i < d; ++i)
                numQuadNodes[i] = basis.degree(i) + 1;
            numQuadNodes[dir] = 1;

            // Setup Quadrature
            rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

            // Set Geometry evaluation flags
            evFlags = NEED_MEASURE | NEED_VALUE | NEED_GRAD_TRANSFORM ;
        }

        // Evaluate on element.
        inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                             const gsFunction<T>    & _vApprox,
                             const gsFunction<T>    & _wApprox,
                             const gsFunction<T>    & _u0Func,
                             gsMatrix<T>            & quNodes)
        {
            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            // Evaluate first function
            _vApprox.eval_into(quNodes, vVals);

            // Evaluate second function (defined on physical domain)
            _wApprox.eval_into(quNodes, wVals);

            // Evaluate right-hand-side function (defined of physical domain)
            _u0Func.eval_into(geoEval.values(), m_funcVals);

        }

        // assemble on element
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & accumulated)
        {
            T sum(0.0);
            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                geoEval.outerNormal(k, side, unormal);

                const T weight = quWeights[k] * unormal.norm();
                /*
                gsInfo << "wVals : " << wVals << "\n";
                gsInfo << "vVals : " << vVals << "\n";
                gsInfo << "m_funcVals : " << m_funcVals << "\n";
                gsInfo << "weight : "   << weight << "\n";
                gsInfo << "unormal : "   << unormal << "\n";
                */
                if (weight != 0.0)  sum += weight * (math::pow(m_funcVals(0, k) - vVals(0, k), 2)
                                                     - 2 * (wVals(0, k) - vVals(0, k)) * (m_funcVals(0, k) - vVals(0, k)));
            }
            accumulated += sum;
            return sum;
        }

    private:

        const std::vector<boundaryInterface>& m_sliceInterfaces;
        const std::vector<patchSide>& m_sides;

        using gsEstimator<T>::m_value;
        using gsEstimator<T>::m_elWise;
        using gsEstimator<T>::patchesPtr;
        gsMatrix<T> vVals, wVals;
        gsMatrix<T> m_funcVals;

        gsVector<T> unormal;

        using gsEstimator<T>::field1;
        using gsEstimator<T>::func2;

        bool m_f2param;
        boxSide side;

    };



} // namespace gismo
