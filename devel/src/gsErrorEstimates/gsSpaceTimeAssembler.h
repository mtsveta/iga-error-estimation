/** @file gsSpaceTimeAssembler.h

    @brief Provides assembler for a homogeneous Heat equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  S. Matculevich
*/


#pragma once

#include <gsPde/gsPoissonHeterogeneousPde.h>
#include <gsErrorEstimates/gsVisitorSpaceTime_.h>
#include <gsAssembler/gsVisitorNeumann.h>

namespace gismo
{

/** @brief
    Implementation of a Space Time Assembler.

    It inherits from the gsPoissonAssembler
*/
template <class T>
class gsSpaceTimeAssembler: public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
/** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Sigma_D, \text{ and Initial Condition: } u(x,0) = u_0(x) \text{ on } \Sigma_0\f]
    \f$ v \f$ is the test function and \f$ \Sigma_D \f$ is the boundary side.
    \param[in] rhs is the right-hand side of the parabolic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsSpaceTimeAssembler( gsMultiPatch<T> const         & patches,
                          gsMultiBasis<T> const         & bases,
                          gsBoundaryConditions<T> const & bconditions,
                          const gsFunction<T>           & rhs,
                          const gsFunction<T>           & theta)
    {
        typename gsPde<T>::Ptr pde( new gsPoissonHeterogeneousPde<T>(patches,bconditions,rhs,theta) );
        Base::initialize(pde, bases, m_options);

    }

    gsSpaceTimeAssembler( gsMultiPatch<T> const         & patches,
                          gsMultiBasis<T> const         & bases,
                          gsBoundaryConditions<T> const & bconditions,
                          const gsFunction<T>           & rhs)
    {
        typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rhs) );
        Base::initialize(pde, bases, m_options);

    }

    gsSpaceTimeAssembler() : Base() {
    }

    virtual ~gsSpaceTimeAssembler() { }

    void refresh()
    {
        Base::scalarProblemGalerkinRefresh();
    }

    /// Main assembly routine
    void assemble()
    {
        GISMO_ASSERT(m_system.initialized(), 
                     "Sparse system is not initialized, call initialize() or refresh()");

        // Reserve sparse system
        m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

        Base::computeDirichletDofs();

        if ( 0 == this->numDofs() ) // Are there any interior dofs ?
        {
            gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
            return;
        }
       
        // Assemble modified volume stiffness and load vector integrals
        Base::template push<gsVisitorSpaceTime_<T> >();

        // Enforce Neumann boundary conditions (in the space-time case do-nothing conditions)
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

        // Assembly is done, compress the matrix
        Base::finalize();

    }
    void basisUpdate(const gsMultiBasis<T> & bases){
        m_bases.clear();
        m_bases.push_back(bases);
    }

    void thetaUpdate(const gsFunctionExpr<real_t> & f, const gsFunctionExpr<real_t> & theta){
        typename gsPde<T>::Ptr updatedPde( new gsPoissonHeterogeneousPde<T>(m_pde_ptr->patches(), m_pde_ptr->bc(), f, theta));
        m_pde_ptr.reset();
        m_pde_ptr = updatedPde;
  }

protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

};


} // namespace gismo
