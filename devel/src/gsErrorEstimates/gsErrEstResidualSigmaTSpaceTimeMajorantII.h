#include<gsAssembler/gsNorm.h>

#pragma once

#include <gsErrorEstimates/gsEstimator.h>

namespace gismo
{

    template <class T>
    class gsErrEstResidualSigmaTSpaceTimeMajorantII : public gsEstimator<T>
    {
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;
        friend class gsEstimator<T>;

    public:

        gsErrEstResidualSigmaTSpaceTimeMajorantII(const gsField<T> & _field1,
                             const gsField<T> & _func2,
                             const std::vector<boundaryInterface>& SliceInterfaces,
                             const std::vector<patchSide>& TopSides,
                             bool _f2param = true)
                :gsEstimator<T>(_field1,_func2), m_sliceInterfaces(SliceInterfaces), m_sides(TopSides), f2param(_f2param)
        {}

    public:

        T compute(bool storeElWise = false)
        {

            this->m_value = 0.0;

            m_elWise.clear();   // clear the vector of element-wise values

            side = * m_sides.begin(); // get the side from the vector of m_sides
/*
            for (typename gsMultiPatch<T>::const_biterator bit = patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit){
                if (bit->side() == side) {
                    //gsInfo << "side : " << side << "\n";
                    side = bit->side();
                    this->apply2Func(*this, storeElWise, bit->patch, side);
                }
            }
*/
            for (typename gsMultiPatch<T>::const_biterator bit = patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit)
                if (bit->side() == side)
                    this->apply2Func(*this, storeElWise, bit->patch, bit->side());


            return this->m_value;
        }

        inline T takeRoot(const T v) { return math::sqrt(v);}


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
                             const gsFunction<T>    & _func1,
                             const gsFunction<T>    & _func2,
                             gsMatrix<T>            & quNodes)
        {
            // Compute geometry related values
            geoEval.evaluateAt(quNodes);

            // Evaluate first function
            _func1.eval_into(quNodes, f1vals);

            // Evaluate second function (defined on physical domain)
            _func2.eval_into(quNodes, f2vals);

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
                //gsInfo << "unormal : "   << unormal << "\n";

                if (weight != 0.0)  sum += weight * ((f1vals.col(k) - f2vals.col(k)).squaredNorm());
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

        gsMatrix<T> f1vals, f2vals;

        gsVector<T> unormal;

        using gsEstimator<T>::field1;
        using gsEstimator<T>::func2;

        bool f2param;
        boxSide side;

    };



} // namespace gismo
