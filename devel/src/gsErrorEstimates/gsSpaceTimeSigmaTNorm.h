#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

template <class T>
class gsSpaceTimeSigmaTNorm : public gsNorm<T>
{
    typedef gsNorm<T> Base;
    typedef typename gsMatrix<T>::RowsBlockXpr Rows;
    friend class gsNorm<T>;

public:

    gsSpaceTimeSigmaTNorm(const gsField<T> & _field1,
                         const gsFunction<T> & _func2,
                         const std::vector<boundaryInterface>& SliceInterfaces,
                         const std::vector<patchSide>& TopSides,
                         const T & theta,
                         bool _f2param = false)
            :gsNorm<T>(_field1,_func2), m_sliceInterfaces(SliceInterfaces), m_sides(TopSides), m_theta(theta), f2param(_f2param)
    {}

    gsSpaceTimeSigmaTNorm(const gsField<T> & _field1,
                         const gsFunction<T> & _func2,
                         const gsFunction<T> & _dfunc2,
                         const std::vector<boundaryInterface>& SliceInterfaces,
                         const std::vector<patchSide>& TopSides,
                          const T & theta,
                          bool _f2param = false)
            : gsNorm<T>(_field1,_func2), m_sliceInterfaces(SliceInterfaces), m_sides(TopSides), m_theta(theta), f2param(_f2param)
    {}


    gsSpaceTimeSigmaTNorm(const gsField<T> & _field1,
                          const std::vector<boundaryInterface>& SliceInterfaces,
                          const std::vector<patchSide>& TopSides,
                          const T & theta,
                          bool _f2param = false)
            : gsNorm<T>(_field1), m_sliceInterfaces(SliceInterfaces), m_sides(TopSides), m_theta(theta), f2param(_f2param)
    {}

public:

    T compute(bool storeElWise = false)
    {

        m_value = 0.0;

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
        /*
        for ( std::vector<patchSide>::const_iterator it=m_sides.begin();it!=m_sides.end();++it) { // *it ---> interface
            gsInfo << "it->patch " << it->patch << "\n";
            gsInfo << "it->side() " << it->side().m_index << "\n";
            this->apply1(*this, storeElWise, it->patch, it->side());
        }

        m_value = takeRoot(m_value);
        return m_value;
         */
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

        //gsInfo << "quNodes \n " << quNodes << "\n";
        //gsInfo << "geoEval.values() \n " << geoEval.values() << "\n";

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        //int d = geoEval.parDim();
        //gsMatrix<real_t> zeros(d, quNodes.cols() / 2.0);
        //zeros.setZero();

        /*
        gsInfo << "d : " << d << "\n";
        gsInfo << "quNodes \n " << quNodes << "\n";
        gsInfo << "geoEval.values() \n " << geoEval.values() << "\n";
        gsInfo << "zero \n " << zeros << "\n";
        */
        //_func1.eval_into(geoEval.values(), f1vals);
        //_func1.deriv_into(geoEval.values(), f1grads);

        //gsInfo << "f1vals (geoEval.values())\n " << f1vals << "\n";
        //gsInfo << "f1grads (geoEval.values())\n " << f1grads << "\n";

        //quNodes.block(0, quNodes.cols() / 2.0, d, quNodes.cols() / 2.0) = zeros;
        //gsInfo << "quNodes \n " << quNodes << "\n";

        // Evaluate first function
        _func1.eval_into(quNodes, f1vals);
        _func1.deriv_into(quNodes, f1grads);

        //quNodes.block(0, quNodes.cols() / 2.0, d, quNodes.cols() / 2.0) = zeros;
        //gsInfo << "quNodes \n " << quNodes << "\n";

        //gsInfo << "f1vals (quNodes)\n " << f1vals << "\n";
        //gsInfo << "f1grads (quNodes)\n " << f1grads << "\n";

        // Evaluate second function (defined of physical domain)
        _func2.eval_into(f2param? quNodes : geoEval.values(), f2vals);
        _func2.deriv_into(f2param? quNodes : geoEval.values(), f2grads);

        //gsInfo << "f2vals \n " << f2vals << "\n";
        //gsInfo << "f2grads \n " << f2grads << "\n";

    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        const int d = element.dim();
        T sum(0.0);
        T normVals(0.0), normGrads(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            geoEval.outerNormal(k, side, unormal);

            geoEval.transformGradients(k, f1grads, ph_f1grads);
            geoEval.transformGradients(k, f2grads, ph_f2grads);

            Rows ph_f1gradsSpace = ph_f1grads.topRows(d-1);
            Rows ph_f2gradsSpace = ph_f2grads.topRows(d-1);

            const T weight = quWeights[k] * unormal.norm();
            //T theta = std::pow(element.getMinCellLength(), 2);
            //theta = 1.0;
            const T h = element.getCellSize();
            T delta_h = m_theta * h;

            //gsInfo << "unormal \n " << unormal << "\n";
            //gsInfo << "side : " << side.m_index << "\n";
            //gsInfo << "ph_f1gradsSpace \n " << ph_f1gradsSpace << "\n";
            //gsInfo << "ph_f2gradsSpace \n " << ph_f2gradsSpace << "\n";
            normVals  = (f1vals.col(k) - f2vals.col(k)).squaredNorm();
            normGrads = (ph_f1gradsSpace - ph_f2gradsSpace).squaredNorm();

            sum += weight * (normVals + delta_h * normGrads);
        }
        accumulated += sum;
        return sum;
    }

    inline T takeRoot(const T v) { return math::sqrt(v);}

private:

    const std::vector<boundaryInterface>& m_sliceInterfaces;
    const std::vector<patchSide>& m_sides;

    using Base::m_value;
    using Base::m_elWise;
    using Base::patchesPtr;


    gsMatrix<T> f1vals, f2vals;
    gsMatrix<T> f1grads, f2grads, ph_f1grads, ph_f2grads;

    gsVector<T> unormal;

    using gsNorm<T>::field1;
    using gsNorm<T>::func2;

    T m_theta;

    bool f2param;
    boxSide side;


};



} // namespace gismo
