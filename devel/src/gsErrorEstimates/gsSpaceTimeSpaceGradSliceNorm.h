#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

    template <class T>
    class gsSpaceTimeSpaceGradSliceNorm : public gsNorm<T>
    {
        friend class gsNorm<T>;
        typedef typename gsMatrix<T>::RowsBlockXpr Rows;
        typedef gsNorm<T> Base;

    public:

        gsSpaceTimeSpaceGradSliceNorm(const gsField<T> & _field1,
                                const gsFunction<T> & _func2,
                                const std::vector<boundaryInterface>& sliceInterfaces,
                                const std::vector<patchSide>& sides,
                                bool _f2param = false)
                :gsNorm<T>(_field1,_func2), m_sliceInterfaces(sliceInterfaces), m_sides(sides), f2param(_f2param)
        {}

        gsSpaceTimeSpaceGradSliceNorm(const gsField<T> & _field1,
                                const gsFunction<T> & _func2,
                                const gsFunction<T> & _dfunc2,
                                const std::vector<boundaryInterface>& sliceInterfaces,
                                const std::vector<patchSide>& sides,
                                bool _f2param = false)
                : gsNorm<T>(_field1,_func2), m_sliceInterfaces(sliceInterfaces), m_sides(sides), f2param(_f2param)
        {}


        gsSpaceTimeSpaceGradSliceNorm(const gsField<T> & _field1,const std::vector<boundaryInterface>& sliceInterfaces,
                                const std::vector<patchSide>& sides)
                : gsNorm<T>(_field1), m_sliceInterfaces(sliceInterfaces), m_sides(sides), f2param(false)
        {}

    public:

        T compute(bool storeElWise = false)
        {

            m_value = T(0.0);
            side = * m_sides.begin(); // get the side from the vector of m_sides

            for (typename gsMultiPatch<T>::const_biterator bit =
                    patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit){
                if (bit->side() == side) {

                    //gsInfo << "side : " << side << "\n";
                    side = bit->side();
                    this->apply1(*this, storeElWise, bit->patch, side);
                }

            }

            return takeRoot(this->m_value);
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
            evFlags = NEED_MEASURE| NEED_VALUE | NEED_GRAD_TRANSFORM | NEED_JACOBIAN;
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
            _func1.deriv_into(quNodes, f1grads);
            //_func1.deriv_into(geoEval.values(), f1grads);

            // Evaluate second function (defined of physical domain)
            _func2.deriv_into(f2param? quNodes : geoEval.values(), f2grads);

            /*
            gsInfo << "quNodes \n " << quNodes << "\n";
            gsInfo << "geoEval.values() \n " << geoEval.values() << "\n";
            gsInfo << "f1grads (quNodes)\n " << f1grads << "\n";
            gsInfo << "f2grads \n " << f2grads << "\n";
             */

        }

        // assemble on element
        inline T compute(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights,
                         T & accumulated)
        {
            T sum(0.0);
            const int d = element.dim();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                geoEval.transformGradients(k, f1grads, ph_f1grads);
                //geoEval.transformGradients(k, f2grads, ph_f2grads);

                Rows ph_f1gradsSpace = ph_f1grads.topRows(d-1);
                Rows ph_f2gradsSpace = f2grads.topRows(d-1);

                geoEval.outerNormal(k, side, unormal);
                const T weight = quWeights[k] * unormal.norm();

                /*
                gsInfo << "ph_f1gradsSpace \n " << ph_f1gradsSpace << "\n";
                gsInfo << "ph_f2gradsSpace \n " << ph_f2gradsSpace << "\n";

                gsInfo << "unormal \n " << unormal << "\n";
                gsInfo << "weight \n " << weight << "\n\n";
                */

                if (weight != 0.0)  sum +=weight * ((ph_f1gradsSpace - ph_f2gradsSpace.col(k)).squaredNorm());
            }
            accumulated += sum;
            return sum;
        }

    private:

        const std::vector<boundaryInterface>& m_sliceInterfaces;
        const std::vector<patchSide>& m_sides;

        using Base::m_value;
        using Base::m_elWise;
        using Base::patchesPtr;

        gsMatrix<T> f1grads, f2grads, ph_f1grads, ph_f2grads;

        gsVector<T> unormal;

        using gsNorm<T>::field1;
        using gsNorm<T>::func2;

        bool f2param;
        boxSide side;
    };



} // namespace gismo
