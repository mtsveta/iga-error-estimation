/** @file gsEstimator.h

    @brief Provides generic routines for computing norms of the residual functionals.
    Structure of the class is very close to sturcture of the gsNorm class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#pragma once


namespace gismo
{

/** @brief Provides generic routines for computing norms of the residual functionals
 * element-wise or globally
 *
 * \ingroup ErrorEstimates
*/
    template <class T>
    class gsEstimator
    {
    public:
        gsEstimator(const gsField<T> & _field1,
                    const gsField<T> & _field2)
                : m_zeroFunction(T(0.0), _field1.parDim()),
                  patchesPtr(&_field1.patches()),
                  field1(&_field1),
                  field2(&_field2),
                  field3(nullptr),
                  func2(&m_zeroFunction)
        {}
        gsEstimator(const gsField<T> & _field1,
                    const gsFunctionSet<T> & _func2)
                : m_zeroFunction(T(0.0), _field1.parDim()),
                  patchesPtr(&_field1.patches()),
                  field1(&_field1),
                  field2(nullptr),
                  field3(nullptr),
                  func2(&_func2)
        { }
        gsEstimator(const gsField<T> & _field1,
                    const gsField<T> & _field2,
                    const gsFunctionSet<T> & _func2)
                    : m_zeroFunction(T(0.0), _field1.parDim()),
                      patchesPtr(&_field1.patches()),
                      field1(&_field1),
                      field2(&_field2),
                      field3(nullptr),
                      func2(&_func2)
        { }

        gsEstimator(const gsField<T> & _field1,
                    const gsField<T> & _field2,
                    const gsField<T> & _field3)
                : m_zeroFunction(T(0.0), _field1.parDim()),
                  patchesPtr(&_field1.patches()),
                  field1(&_field1),
                  field2(&_field2),
                  field3(&_field3),
                  func2(&m_zeroFunction)
        { }
        explicit gsEstimator(const gsField<T> & _field1)
                             : m_zeroFunction(gsVector<T>::Zero(_field1.dim()), _field1.parDim()),
                               patchesPtr(&_field1.patches()),
                               field1(&_field1),
                               field2(nullptr),
                               func2(&m_zeroFunction),
                               field3(nullptr)
        { }

        //virtual T compute(bool storeElWise = false);
        //virtual T compute(bool storeElWise = false, real_t elemNum);

        /*{
            this->apply2Func(*this,storeElWise);
            return m_value;
        }*/
        /*
        virtual T compute(bool storeElWise = false, real_t elemNumber)
        {
            this->apply3Func(*this,storeElWise, elemNumber);
            return m_value;
        }*/
        /** \brief Main function for the norm computation for 2 functions
         *
         * The computed value can be accessed by value().
         *
         * \param[in] visitor The Norm-visitor to be used.
         * \param[in] storeElWise Flag indicating whether the element-wise norms should be stored.
         * See also elementNorms().
         * \param[in] side To be used, if the visitor will iterate over a side of the domain.
         */
        template <class NormVisitor>
        void apply2Func(NormVisitor & visitor, bool storeElWise = false, int patchIndex = 0, boxSide side = boundary::none)
        {
            if ( storeElWise )  {   m_elWise.clear();}

            gsMatrix<T> quNodes; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            m_value = T(0.0);

            for (unsigned pn=0; pn < patchesPtr->nPatches(); ++pn )// for all patches
            {
                const gsFunction<T> & f1  = field1->function(pn);
                const gsFunction<T> & f2  = field2->function(pn);
                // Obtain an integration domain
                const gsBasis<T> & dom = field1->isParametrized() ?
                                         field1->igaFunction(pn).basis() : field1->patch(pn).basis();

                // Initialize visitor
                visitor.initialize(dom, QuRule, evFlags);

                // Initialize geometry evaluator
                typename gsGeometry<T>::Evaluator geoEval(
                        patchesPtr->patch(pn).evaluator(evFlags));

                typename gsBasis<T>::domainIter domIt = dom.makeDomainIterator(side);

                /*
               std::vector<gsBasis<T>::domainIter> iterators;
               for (; domIt->good(); domIt->next()){
                   iterators.emplace_back(domIt);
               }

               //#pragma omp parallel for num_threads(4)
               for (index_t i = 0; i < iterators.size(); i++){
                   typename gsBasis<T>::domainIter domIt = iterators[i];
                   // Map the Quadrature rule to the element
                   QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                   //gsInfo << "domIt->lowerCorner() = \n" << domIt->lowerCorner() << "\n";
                   //gsInfo << "domIt->upperCorner() = \n" << domIt->upperCorner() << "\n";
                   //gsInfo << "quNodes = \n" << quNodes << "\n";
                   //gsInfo << "quWeights = \n" << quWeights << "\n";

                   // Evaluate on quadrature points
                   visitor.evaluate(*geoEval, func1, func2p, quNodes);

                   // Accumulate value from the current element (squared)
                   const T result = visitor.compute(*domIt, *geoEval, quWeights, m_value);
                   if ( storeElWise )
                       m_elWise.emplace_back( visitor.takeRoot(result) );
                }
                */
                for (; domIt->good(); domIt->next())
                {
                    // Map the Quadrature rule to the element
                    QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                    // Evaluate on quadrature points
                    visitor.evaluate(*geoEval, f1, f2, quNodes);

                    // Accumulate value from the current element (squared)
                    const T result = visitor.compute(*domIt, *geoEval, quWeights, m_value);

                    if ( storeElWise )
                        m_elWise.push_back( visitor.takeRoot(result) );
                }
            }

            m_value = visitor.takeRoot(m_value);

        }

        /** \brief Main function for the norm computation for 3 functions
        *
        * The computed value can be accessed by value().
        *
        * \param[in] visitor The Norm-visitor to be used.
        * \param[in] storeElWise Flag indicating whether the element-wise norms should be stored.
        * See also elementNorms().
        * \param[in] side To be used, if the visitor will iterate over a side of the domain.
        */
        template <class NormVisitor>
        void apply3Func(NormVisitor & visitor, bool storeElWise = false, int patchIndex = 0, boxSide side = boundary::none)
        {
            if ( storeElWise )  {   m_elWise.clear();}

            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            m_value = T(0.0);

            for (unsigned pn=0; pn < patchesPtr->nPatches(); ++pn )// for all patches
            {

                const gsFunction<T> & f1  = field1->function(pn);
                const gsFunction<T> & f2  = field2->function(pn);
                const gsFunction<T> & func2p = func2->function(pn);

                // Obtain an integration domain
                const gsBasis<T> & dom = field1->isParametrized() ?
                                         field1->igaFunction(pn).basis() : field1->patch(pn).basis();

                // Initialize visitor
                visitor.initialize(dom, QuRule, evFlags);

                // Initialize geometry evaluator
                typename gsGeometry<T>::Evaluator geoEval(patchesPtr->patch(pn).evaluator(evFlags));

                //gsInfo << "dom.makeDomainIterator(side)" << dom.makeDomainIterator(side) << "\n\n";
                typename gsBasis<T>::domainIter domIt = dom.makeDomainIterator(side);

                for (; domIt->good(); domIt->next())
                {
                    // Map the Quadrature rule to the element
                    QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                    // Evaluate on quadrature points
                    visitor.evaluate(*geoEval, f1, f2, func2p, quNodes);

                    // Accumulate value from the current element (squared)
                    const T resultEl = visitor.compute(*domIt, *geoEval, quWeights, m_value);

                    if ( storeElWise )
                        m_elWise.push_back( visitor.takeRoot(resultEl) );
                }
            }

            m_value = visitor.takeRoot(m_value);
        }

        /** \brief Main function for the norm computation for 3 functions
        *
        * The computed value can be accessed by value().
        *
        * \param[in] visitor The Norm-visitor to be used.
        * \param[in] storeElWise Flag indicating whether the element-wise norms should be stored.
        * See also elementNorms().
        * \param[in] side To be used, if the visitor will iterate over a side of the domain.
        */
        template <class NormVisitor>
        void apply3Fields(NormVisitor & visitor, bool storeElWise = false, boxSide side = boundary::none)
        {
            if ( storeElWise )  {   m_elWise.clear();}

            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            m_value = T(0.0);

            for (unsigned pn=0; pn < patchesPtr->nPatches(); ++pn )// for all patches
            {

                const gsFunction<T> & f1 = field1->function(pn);
                const gsFunction<T> & f2 = field2->function(pn);
                const gsFunction<T> & f3 = field3->function(pn);

                // Obtain an integration domain
                const gsBasis<T> & dom = field1->isParametrized() ?
                                         field1->igaFunction(pn).basis() : field1->patch(pn).basis();

                // Initialize visitor
                visitor.initialize(dom, QuRule, evFlags);

                // Initialize geometry evaluator
                typename gsGeometry<T>::Evaluator geoEval(patchesPtr->patch(pn).evaluator(evFlags));

                //gsInfo << "dom.makeDomainIterator(side)" << dom.makeDomainIterator(side) << "\n\n";
                typename gsBasis<T>::domainIter domIt = dom.makeDomainIterator(side);

                for (; domIt->good(); domIt->next())
                {
                    // Map the Quadrature rule to the element
                    QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                    // Evaluate on quadrature points
                    visitor.evaluate(*geoEval, f1, f2, f3, quNodes);

                    // Accumulate value from the current element (squared)
                    const T resultEl = visitor.compute(*domIt, *geoEval, quWeights, m_value);

                    if ( storeElWise )
                        m_elWise.push_back( visitor.takeRoot(resultEl) );
                }
            }

            m_value = visitor.takeRoot(m_value);
        }

    public:

        /// @brief Return the multipatch
        const gsMultiPatch<T> & patches() const { return *patchesPtr; }

        /** @brief Function elementNorms() returns the computed norm values element-wise.
         *
         *
         * \returns The SQUARED element-norms in a std::vector.\n
         * It is assumed that they were actually computed by providing
         * the proper flag \em storeEltWise in the call of apply().
         *
         * \remarks The order of the element-norms is "defined"
         * firstly by the numbering of the patches, and secondly
         * by the order in which the gsDomainIterator iterates
         * over the elements of the mesh! As of now (02.Dec.2014),
         * there is no way to access a particularly numbered element
         * on a mesh directly; you have to run the corresponding
         * gsDomainIterator again to reach the respective element.
         *
         */
        const std::vector<T> & elementNorms() const { return m_elWise; }

        /// @brief Returns the computed norm value
        T value() const { return m_value; }


    private:
        // constant function zero that is used  when no function is specified
        // to calculate the norm and not the distance
        // maybe its possible to optimize this
        const gsConstantFunction<T> m_zeroFunction;

    protected:

        const gsMultiPatch<T> * patchesPtr;

        const gsField<T>    * field1;
        const gsField<T>    * field2;
        const gsField<T>    * field3;

        const gsFunctionSet<T> * func2;

    protected:

        std::vector<T> m_elWise;    // vector of the element-wise values of the norm
        T              m_value;     // the total value of the norm
    };

} // namespace gismo

