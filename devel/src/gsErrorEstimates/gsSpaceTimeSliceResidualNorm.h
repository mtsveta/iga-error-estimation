#include<gsAssembler/gsNorm.h>

#pragma once

#include <gsErrorEstimates/gsEstimator.h>

namespace gismo
{

    template <class T>
    class gsSpaceTimeSliceResidualNorm : public gsEstimator<T>
    {
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;
        friend class gsNorm<T>;

    public:

        gsSpaceTimeSliceResidualNorm(const gsField<T> & _field1,
                                     const gsField<T> & _func2,
                                     const std::vector<boundaryInterface>& SliceInterfaces,
                                     const std::vector<patchSide>& TopSides,
                                     bool _f2param = true)
                :gsEstimator<T>(_field1,_func2), m_sliceInterfaces(SliceInterfaces), m_sides(TopSides), f2param(_f2param)
        {}

    public:

        inline T compute(bool storeElWise) {

            m_value = 0.0;

            for ( std::vector<patchSide>::const_iterator it=m_sides.begin();it!=m_sides.end();++it) // *it ---> interface
            {
                gsInfo << "it->patch " << it->patch << "\n";
                gsInfo << "it->side() " << it->side().m_index << "\n";
                this->apply2Func(*this, storeElWise, it->patch, it->side());
            }
            return this->m_value;
        }


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
            _func1.eval_into(f1param? quNodes : geoEval.values(), f1vals);

            // Evaluate second function (defined of physical domain)
            _func2.eval_into(f2param? quNodes : geoEval.values(), f2vals);

            gsInfo << "f1vals \n " << f1vals << "\n";
            gsInfo << "f2vals \n " << f2vals << "\n";

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

                gsInfo << "normVals : " << (f1vals.col(k) - f2vals.col(k)).squaredNorm() << "\n";
                gsInfo << "weight : "   << weight << "\n";
                gsInfo << "unormal : "   << unormal << "\n";
                gsInfo << "quWeights[k] : "   << quWeights[k] << "\n";

                if (weight != 0.0)  sum += weight * ((f1vals.col(k) - f2vals.col(k)).squaredNorm());
            }
            accumulated += sum;
            return sum;
        }

        inline T takeRoot(const T v) { return math::sqrt(v);}

    private:

        const std::vector<boundaryInterface>& m_sliceInterfaces;
        const std::vector<patchSide>& m_sides;

        using gsEstimator<T>::m_value;
        using gsEstimator<T>::m_elWise;

        gsMatrix<T> f1vals, f2vals;

        gsVector<T> unormal;

        using gsEstimator<T>::field1;
        using gsEstimator<T>::func2;

        bool f2param;
        boxSide side;

    };



} // namespace gismo
