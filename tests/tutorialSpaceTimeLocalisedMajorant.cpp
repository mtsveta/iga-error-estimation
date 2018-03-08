/** @file tutorialSpaceTimeLocalisedMajorant.cpp

    @brief Example for testing the gsSpaceTimeLocalisedSolver
    with adaptive refinement with THB-splines using
    functional error estimates as a refinement criteria.


    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/
#include <iostream>
#include <gismo.h>

#include <gsErrorEstimates/gsSpaceTimeLocalisedAssembler.h>
#include <gsErrorEstimates/gsSpaceTimeLocalisedNorm.h>
#include <gsErrorEstimates/gsSpaceTimeLocalisedSigmaTNorm.h>
#include <gsErrorEstimates/gsTestSpaceTimeLocalisedMajorant.h>


#include <gsErrorEstimates/gsSpaceTimeSliceNorm.h>
#include <gsErrorEstimates/gsSpaceTimeSpaceGradSliceNorm.h>
#include <gsErrorEstimates/gsErrEstSpaceTimeResidual.h>
#include <gsErrorEstimates/gsNormFields.h>
#include <gsErrorEstimates/gsSpaceTimeSpaceGradNorm.h>
#include <gsErrorEstimates/gsSpaceTimeSolOperNorm.h>
#include <gsErrorEstimates/gsSpaceTimeDeltaxNorm.h>
#include <gsErrorEstimates/gsSpaceTimeDtNorm.h>
#include <gsErrorEstimates/gsSpaceTimeDeltaxDeltaKNorm.h>
#include <gsErrorEstimates/gsSpaceTimeDtDeltaKNorm.h>
#include <gsErrorEstimates/gsSpaceTimeErrorIdentity.h>
#include <gsErrorEstimates/gsErrEstSpaceTimeMinorant.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

// S. Matculevich
//
// This is a test example for a illustrating the adaptive refinement procedure implemented for the gsPoissonAssembler
// Flags, parameters, geometry and prescribed exact solution are specified within the main() function

bool gsParseCommandLine(int argc, char **argv, bool plot)
{
    //! [Parse command line]
    gsCmdLine cmd("Tutorial on solving a heat eqaution with guaranteed error control using the funcional error estimate.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the gmresSolver", plot);
    cmd.getValues(argc,argv);
    return true;
}

int main(int argc, char *argv[])
{

    // Define stopwatch to measure the performance of the routines
    gsCPUStopwatch clock, clock_total;
    clock_total.restart();

    //! [Initialize Test Parameters]
    // -------------------------------------------------------------------------------------------------------------- //
    // Define constants and preliminaries
    const int NUM_PATCHES = 1;     // All the geometries are single-patch geometries

    real_t rho = 1.0 / 16.0;
    real_t TOL = 1e-3;

    // Define parameters of I/O
    bool plotToParaview = false; // Flag indicating whether objects must be plotted in ParaView
    bool saveToFile     = true;  // Flag indicating whether objects must be plotted in ParaView
    bool isAdaptive     = true;
    bool withMajorant   = true;
    bool withMajorantOptimization = true;
    bool withMajorantEquilibration = false;

    if ( !gsParseCommandLine(argc, argv, plotToParaview) ) return 0;

    // Define test-case parameters (number and dimension)
    const unsigned
    // exampleNumber(2), d(2);        // 2 example: 2d unit square, u = (1 - x)*x*x*(1 - t)*t
    // exampleNumber(3), d(2);        // 3 example: 2d unit square, u = sin(pi*x)*sin(pi*t)
    // exampleNumber(4), d(2);        // 4 example: 2d unit square, u = sin(6.0*pi*x)*sin(3.0*pi*t)
    // exampleNumber(5), d(2);        // 5 example: 2d unit square, u = cos(x)*exp(t)
    // exampleNumber(6), d(2);        // 6 example: 2d unit square, u = (x^2 - x) * (y^2 - y) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))
    // exampleNumber(7), d(2);        // 7 example: 2d rectangle $(0, 2) x (0, 1)$, u = (2 - x)*x*x*(1 - y)*y
    // exampleNumber(8), d(2);        // 8 example: 2d rectangle $(0, 2) x (0, 1)$, u = (x^2 - 2*x) * (y^2 - y) * exp(-100 * ((x - 1.4)^2 + (y - 0.95)^2))
    // exampleNumber(9), d(2);     // 9 example: 2d rectangle (0, 1) x (0, 2), u = (x^2 + t^2)^(1.0/3.0) * sin(2.0/3.0*atan2(t,x) + pi)
    // exampleNumber(10), d(2);     // 10 example: 2d rectangle (0, 1) x (0, 2), u = cos(x)*exp(t)
    // exampleNumber(11), d(2);     // 11 example: 2d unit square, u = sin(1 / (1/pi/10 + (x^2 + y^2)^(1.0/2.0)))
    // exampleNumber(12), d(2);        // 11 example: 2d unit square, u =
    // exampleNumber(16), d(2);       // 22 example: 2d unit square, u = (x^2 - x)*(y^2 - y)*exp(-100*((x - 0.25)^2 + (y - 0.25)^2))
    // exampleNumber(33), d(2);       // 33 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)
    // exampleNumber(34), d(2);       // 34 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^0.5
    // exampleNumber(35), d(2);        // 35 example: 2d [0, 1] x (0, 2), u = sin(pi*x)*abs(1 - y)^1.5
    // exampleNumber(37), d(2);
    // exampleNumber(38), d(2);

    // 3d examples:
    // exampleNumber(20), d(3);     // 20 example: unit cube, u = sin(pi*x)*sin(pi*y)*sin(pi*t),        non-homogeneous BC
    // exampleNumber(25), d(3);     // 25 example: unit cube, u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(21), d(3);       // 21 example: unit cube, u = cos(x)*exp(y)*sin(pi*t),              non-homogeneous BC
    exampleNumber(23), d(3);     // 23 example: 2d+1 quater annulus + [0, 1] in time, u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2
    // exampleNumber(24), d(3);     // 24 example: 2d+1 quater annulus + [0, 1] in time, u = sin(pi*x)*sin(pi*y)*sin(pi*z)
    // exampleNumber(26), d(3);       // 26 example: 2d+1 quater annulus + [0, 1] in time, u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2
    // exampleNumber(27), d(3);     // 27 example: square x [0, 2], u = (1 - x)*x^2*(1 - y)*y^2*(2 - z)*z^2,  homogeneous BC
    // exampleNumber(28), d(3);     // 28 example: [0, 2] x [0, 1] x [0, 1], u = (2 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(29), d(3);      // 29 example: [0, 2] x [0, 3] x [0, 1], u = (2 - x)*x^2*(3 - y)*y^2*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(30), d(3);      // 30 example: 2d+1 quater annulus + [0, 1] in time, u = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)*(1 - z)*z^2,  homogeneous BC
    // exampleNumber(31), d(3);      // 30 example: G-domain + [0, 1] in time, u = (x^2 + y^2 - 1)*(x^2 + y^2 - 4)*exp(-100 * ((z - 0.8)^2)),  homogeneous BC
    // exampleNumber(15), d(3);     // 15 example: 3d l-shape x (0, 2), u = if( y > 0, (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // exampleNumber(31), d(3);     // 31 example: 3d unit cube, (x^2 - x)*(y^2 - y)*(z^2 - z)*exp(-100*((x - 0.25)^2 + (y - 0.25)^2 + (z - 0.25)^2))
    // exampleNumber(32), d(3);     // 32 example: 3d, unit square + [0, 2] in time, u = if( y > 0, (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // exampleNumber(36), d(3);     // 32 example: 3d, unit square + [0, 2] in time, u = if( y > 0, (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (z^2 + z + 1) * (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )

    // Define test-case dimentions
    //const unsigned d(2);
    //const unsigned d(3);    // Dimension must be prescribed

    gsTestSpaceTimeLocalisedMajorant<d> testSpaceTime(exampleNumber,
                                             isAdaptive,
                                             withMajorant, withMajorantOptimization,
                                             withMajorantEquilibration);

    // Init the degree of basis S^{p, p}, where p = vDegree
    int vDegree(2), m(1), l(1);
    // Init the degree of the basis for the flux: y \in S^{p + k, p + k} + S^{p + k, p + k}, p = vDegree
    // yDegree must be >= 2 to ensure that S^{p + k, p + k} + S^{p + k, p + k} \subset H(\Omega, div§)
    int yDegree(vDegree + m);
    int wDegree(vDegree + l);

    // Setting up the refinement strategy
    // Number of initial uniform and total unif./adapt. refinement steps
    unsigned int numInitUniformRefV(1), numInitUniformRefY(1), numInitUniformRefW(1), numTotalAdaptRef(6);

    MarkingStrategy adaptRefCrit(BULK); // with alternatives GARU, PUCA, and BULK
    //MarkingStrategy adaptRefCrit(GARU);
    real_t markingParamTheta(0.4);  // parameter theta in marking strategy
    unsigned int yBasisRefDelay(3); // parameter for the delay in updating y_h basis for the refinement
    unsigned int wBasisRefDelay(3); // parameter for the delay in updating w-h basis for the refinement

    testSpaceTime.gsCreateResultsFolder(saveToFile, exampleNumber,
                                        vDegree, yDegree, wDegree, yBasisRefDelay, wBasisRefDelay,
                                        numTotalAdaptRef, adaptRefCrit, markingParamTheta);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Initialize Test Parameters]

    //! [Define Problem Data]
    // -------------------------------------------------------------------------------------------------------------- //
    gsFunctionExpr<> uDFunc, fFunc, uFunc;
    gsBoundaryConditions<> bcInfo;
    testSpaceTime.gsInitializeProblemData(uDFunc, fFunc, uFunc, bcInfo);
    testSpaceTime.gsLogProblemData();

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Problem Data]

    //! [Define Basis]
    // -------------------------------------------------------------------------------------------------------------- //
    gsMultiBasis<> basisV, basisY, basisW;
    testSpaceTime.gsGetInitialBasis(vDegree, yDegree, wDegree,
                                    basisV, basisY, basisW,
                                    numInitUniformRefV,
                                    numInitUniformRefY,
                                    numInitUniformRefW);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Basis]

    //! [Define Auxiliary Structures for Storing of the Results]
    // -------------------------------------------------------------------------------------------------------------- //
    // Initialize arrays to store the error and estimators on each refinement step
    gsVector<real_t> eL2Vector(numTotalAdaptRef),           // || e_i ||_Q          := || u - u_i ||_Q
            eH1Vector(numTotalAdaptRef),                    // || grad e_i ||_Q     := || grad(u - u_i) ||_Q
            eSpaceTimeVector(numTotalAdaptRef),             // ||| e_i |||_h        := (|| grad_x e_i ||^2_Q + detla_h * || (u - u_i)_t ||^2_Q)^1/2
            eFullSpaceTimeVector(numTotalAdaptRef),         // ||| e_i |||_h        := (|| grad_x e_i ||^2_Q + detla_h * || (e_i)_t ||^2_Q + || e_i ||^2_T + detla_h * || grad_x e_i ||^2_T)^1/2
            eFullSpaceTimeLocVector(numTotalAdaptRef),      // ||| e_i |||_{h, loc} := (|| grad_x e_i ||^2_Q + || e_i ||^2_T + 1/2 * \sum_{K \in K_h} detla_K * (|| (e_i)_t ||^2_K)^1/2
            eFullSpaceTimeFullLocVector(numTotalAdaptRef),  // ||| e_i |||_{h, loc} := (|| grad_x e_i ||^2_Q + || e_i ||^2_T + 1/2 * \sum_{K \in K_h} detla_K * (|| (e_i)_t ||^2_K - || laplas_x e_i ||^2_K)^1/2
            eSpaceTimeSpaceGradVector(numTotalAdaptRef),    // || grad_x e_i ||_Q
            eFullSpaceTimeSpaceGradVector(numTotalAdaptRef),// (|| grad_x e_i ||^2_Q + || e_i ||^2_T)^2
            eSpaceTimeSolOperVector(numTotalAdaptRef),      // (|| Delta_x e_i ||^2_Q + || (e_i)_t||^2_Q )^1/2
            eSpaceTimeDeltaxVector(numTotalAdaptRef),       // || Delta_x e_i ||_Q
            eSpaceTimeDeltaxDeltaKVector(numTotalAdaptRef), // (delta_K * || Delta_x e_i ||^_Q)^1/2
            eSpaceTimeDtVector(numTotalAdaptRef),           // || (e_i)_t||^2_Q
            eSpaceTimeDtDeltaKVector(numTotalAdaptRef),     // (delta_K * || (e_i)_t||^2_Q)^1/2
            eFullSpaceTimeSolOperVector(numTotalAdaptRef),  // || e_i ||_{L, Q}   := (|| Delta_x e_i ||^2_Q + || (e_i)_t||^2_Q + || grad_x e_i ||^2_T)^1/2
            eIdentVector(numTotalAdaptRef),                 // || Delta_x u_i + f - (u_i)_t ||_Q
            eFullIdentVector(numTotalAdaptRef),             // (|| Delta_x u_i + f - (u_i)_t ||_Q + || grad_x e_i ||_0)^1/2
            relErrorVector(numTotalAdaptRef-1),             // delta_i = || u_i - u_{i-1} ||, i > 0
            relError0Vector(numTotalAdaptRef-1),            // eps0_i  = || u_i - u_0 ||, i > 0
            thetaVector(numTotalAdaptRef-1),                // theta_i  = || u - u_i || / || u - u_{i-1}||, i > 0
            stopcritVector(numTotalAdaptRef-1),             // stop_crit_i = (1 - theta_i) * TOL * eps0_i * rho, i > 0
            majVector(numTotalAdaptRef),
            mdVector(numTotalAdaptRef),
            meqVector(numTotalAdaptRef),
            majhVector(numTotalAdaptRef),
            minVector(numTotalAdaptRef),
            majDeltaHVector(numTotalAdaptRef),
            etaVector(numTotalAdaptRef),
            hmaxVector(numTotalAdaptRef),
            hminVector(numTotalAdaptRef);
    majVector.setZero(numTotalAdaptRef);
    majhVector.setZero(numTotalAdaptRef);
    minVector.setZero(numTotalAdaptRef);

    gsVector<real_t> gradxeTVector(numTotalAdaptRef),
            gradxe0Vector(numTotalAdaptRef),
            eTVector(numTotalAdaptRef),
            e0Vector(numTotalAdaptRef),
            ehSigmaTVector(numTotalAdaptRef);

    gsVector<real_t> majIIVector, majIIGapVector;
    majIIVector.setZero(numTotalAdaptRef);
    majIIGapVector.setZero(numTotalAdaptRef);

    gsVector<real_t> eWL2Vector(numTotalAdaptRef), eWH1Vector(numTotalAdaptRef);
    eWL2Vector.setZero(numTotalAdaptRef);
    eWH1Vector.setZero(numTotalAdaptRef);

    // Initialize arrays to DOFs for v and y on each refinement step
    gsVector<index_t> vDOFs(numTotalAdaptRef), yDOFs(numTotalAdaptRef), wDOFs(numTotalAdaptRef);
    // Initialize vectors with assembling and computation times on each refinement step
    gsVector<double> timeAsmbV(numTotalAdaptRef),
            timeAsmbDivDivY(numTotalAdaptRef),
            timeAsmbMMY(numTotalAdaptRef),
            timeAsmbY(numTotalAdaptRef),
            timeAsmbW(numTotalAdaptRef),
            timeAsmbH1Error(numTotalAdaptRef),
            timeAsmbL2Error(numTotalAdaptRef),
            timeAsmbSpaceTimeError(numTotalAdaptRef),
            timeAsmbSpaceTimeGradSpaceError(numTotalAdaptRef),
            timeAsmbSpaceTimeSolOperError(numTotalAdaptRef),
            timeAsmbSpaceTimeDeltaxError(numTotalAdaptRef),
            timeAsmbSpaceTimeDtError(numTotalAdaptRef),
            timeAsmbSpaceTimeDeltaxDeltaKError(numTotalAdaptRef),
            timeAsmbSpaceTimeDtDeltaKError(numTotalAdaptRef),
            timeAsmbSpaceTimeErrorIdentity(numTotalAdaptRef),
            timeAsmbMajorant(numTotalAdaptRef),
            timeAsmbMajorantII(numTotalAdaptRef),
            timeAsmbDeltaHMajorant(numTotalAdaptRef),
            timeAsmbMinorant(numTotalAdaptRef),
            timeAsmbEtaIndicator(numTotalAdaptRef),
            timeAsmbEquilY(numTotalAdaptRef);


    int numOfSolvers = 2;
    // matrix are used since we compare the performance of different solvers
    // [0] is the direct solver
    // [1] is the iterative solver
    gsMatrix<double> timeSolvV(numTotalAdaptRef, numOfSolvers),
            timeSolvY(numTotalAdaptRef, numOfSolvers),
            timeSolvW(numTotalAdaptRef, numOfSolvers),
            timeSolvEquilY(numTotalAdaptRef, numOfSolvers);
    timeAsmbMajorant.setZero(numTotalAdaptRef);
    timeAsmbDeltaHMajorant.setZero(numTotalAdaptRef);
    timeAsmbMajorantII.setZero(numTotalAdaptRef);
    timeSolvY.setZero(numTotalAdaptRef, 2);
    timeSolvV.setZero(numTotalAdaptRef, 2);
    timeSolvW.setZero(numTotalAdaptRef, 2);
    etaVector.setZero(numTotalAdaptRef);

    // Initialize auxiliary vectors of all the basis' for y to store the history along the refinement
    std::vector< gsMultiBasis<> > basisYVector;
    std::vector< gsMultiBasis<> > basisWVector;

    // Initialize auxiliary vectors of with reconstructed solutions fields and multipatches
    std::vector< gsField<> > solutionFieldVector;
    std::vector< gsMultiPatch<> > solutionMPVector;

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Auxiliary Structures to Store the Results]
    gsPoissonPde<> heatPde(testSpaceTime.patches, bcInfo, fFunc);

    gsSpaceTimeLocalisedAssembler<real_t> spaceTimeAssemblerV;
    //if (!testSpaceTime.isAdaptive)
    spaceTimeAssemblerV = gsSpaceTimeLocalisedAssembler<real_t>(testSpaceTime.patches, basisV, bcInfo, *heatPde.rhs());
    //else
    //    spaceTimeAssemblerV = gsSpaceTimeAssembler<real_t>(testSpaceTime.patches, basisV, bcInfo, *heatPde.rhs());
    //spaceTimeAssemblerV.initialize(heatPde, basisV);
    spaceTimeAssemblerV.options().setInt("DirichletValues", dirichlet::l2Projection);
    spaceTimeAssemblerV.options().setInt("InterfaceStrategy", iFace::glue);

    gsSpaceTimeLocalisedAssembler<real_t> spaceTimeAssemblerW;
    //if (!testSpaceTime.isAdaptive)
    spaceTimeAssemblerW = gsSpaceTimeLocalisedAssembler<real_t>(testSpaceTime.patches, basisW, bcInfo, *heatPde.rhs());
    //else
    //    spaceTimeAssemblerW = gsSpaceTimeAssembler<real_t>(testSpaceTime.patches, basisW, *bcInfoP, *heaPdeRhsP);
    //spaceTimeAssemblerW.initialize(heatPde, basisW);
    spaceTimeAssemblerW.options().setInt("DirichletValues", dirichlet::l2Projection);
    spaceTimeAssemblerW.options().setInt("InterfaceStrategy", iFace::glue);

    gsMultiPatch<> mpV, mpY, mpW;
    gsMatrix<> vVector(1, 1), yVector(1, 1), wVector(1, 1);
    gsMatrix<> vRefVector(1, 1);
    gsField<> w;
    gsField<> v;

    //! [Refinement Iterations]
    // -------------------------------------------------------------------------------------------------------------- //
    for( unsigned int refCount = 0; refCount < numTotalAdaptRef ; refCount++ )
    {
        testSpaceTime.gsLogRefinementBasisInfo(refCount, NUM_PATCHES, numTotalAdaptRef,
                                               spaceTimeAssemblerV.multiBasis(NUM_PATCHES-1), basisY,
                                               spaceTimeAssemblerW.multiBasis(NUM_PATCHES-1));

        hmaxVector[refCount] = spaceTimeAssemblerV.multiBasis(0).basis(0).getMaxCellLength();
        hminVector[refCount] = spaceTimeAssemblerV.multiBasis(0).basis(0).getMinCellLength();

        testSpaceTime.gsRecontructV(refCount, spaceTimeAssemblerV, bcInfo, vVector, mpV, v, vDOFs,
                                    stopcritVector, timeAsmbV, timeSolvV);
        testSpaceTime.gsRecontructW(refCount, spaceTimeAssemblerW, bcInfo, wVector, mpW, w, wDOFs,
                                    stopcritVector, timeAsmbW, timeSolvW);

        solutionFieldVector.push_back(v);
        solutionMPVector.push_back(mpV);

        //! [Error, Majorant (Optimal Flux), and Residual Estimate Computation]
        // ---------------------------------------------------------------------------------------------------------- //
        // ---------------------------------------------------------------------------------------------------------- //
        int elemNumber = spaceTimeAssemblerV.multiBasis(0).basis(0).numElements();
        std::vector<real_t> mdDistr, mIIdDistr, eH1Distr, eL2Distr, eH2Distr,
                eSpaceTimeDistr, eSpaceTimeSpaceGradDistr,
                eSpaceTimeSolOperDistr, eSpaceTimeDeltaxDistr, eSpaceTimeDtDistr,
                etaDistr, minDistr, eIdentDistr;
        mdDistr.resize(elemNumber);
        eH1Distr.resize(elemNumber);
        etaDistr.resize(elemNumber);
        minDistr.resize(elemNumber);
        eIdentDistr.resize(elemNumber);
        mIIdDistr.resize(elemNumber);

        std::vector<gsFunctionExpr<> *> uForOMP;
        std::vector<gsFunctionExpr<> *> fForOMP;
        index_t numOfU(17);
        index_t numOfF(1);
        index_t numOfWatches(21);
        std::vector<gsCPUStopwatch> watches(numOfWatches);
        for (index_t i = 0; i < numOfU; i++)    uForOMP.emplace_back(new gsFunctionExpr<>(uFunc));
        for (index_t i = 0; i < numOfF; i++)    fForOMP.emplace_back(new gsFunctionExpr<>(fFunc));

        // add top and bottom m_sides into the collections
        std::vector<patchSide> topSides, bottomSides;
        topSides.push_back(spaceTimeAssemblerV.patches().boundaries()[d == 2 ? 0 : d*(d-1) - 1]); // spaceTimeAssemblerV.patches().boundaries()[5]:back (top of the cylinder)
        bottomSides.push_back(spaceTimeAssemblerV.patches().boundaries()[d == 2 ? 1 : d*(d-1) - 2]); // spaceTimeAssemblerV.patches().boundaries()[4]:front (bottom of the cylinder)

#pragma omp parallel sections
        {

#pragma omp section
            {
                gsSpaceTimeSolOperNorm<real_t> eSpaceTimeSolOper(v, * uForOMP[1]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeSolOper, eSpaceTimeSolOperDistr, elemNumber,
                                                      eSpaceTimeSolOperVector, timeAsmbSpaceTimeSolOperError, refCount);
                gsInfo << "t_{e/w} (||| e |||_L)        = " <<  timeAsmbSpaceTimeSolOperError[refCount] << " sec.\n";
            }
#pragma omp section
            {
                gsSpaceTimeDeltaxNorm<real_t> eSpaceTimeDeltax(v, * uForOMP[2]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeDeltax, eSpaceTimeDeltaxDistr, elemNumber,
                                                      eSpaceTimeDeltaxVector, timeAsmbSpaceTimeDeltaxError, refCount);
                gsInfo << "t_{e/w} (|| Delta_x e ||)    = " <<  timeAsmbSpaceTimeDeltaxError[refCount] << " sec.\n";
            }
#pragma omp section
            {
                gsSpaceTimeDeltaxDeltaKNorm<real_t> eSpaceTimeDeltaxDeltaK(v, * uForOMP[15]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeDeltaxDeltaK, eSpaceTimeDeltaxDistr, elemNumber,
                                                      eSpaceTimeDeltaxDeltaKVector, timeAsmbSpaceTimeDeltaxDeltaKError, refCount);
                gsInfo << "t_{e/w} (delta_K * || Delta_x e ||)    = " <<  timeAsmbSpaceTimeDeltaxDeltaKError[refCount] << " sec.\n";
                gsInfo << "delta_K * || Delta_x e || = " <<  eSpaceTimeDeltaxDeltaKVector[refCount] << " sec.\n";

            }

            #pragma omp section
            {
                gsSpaceTimeDtNorm<real_t> eSpaceTimeDt(v, * uForOMP[3]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeDt, eSpaceTimeDtDistr, elemNumber,
                                                      eSpaceTimeDtVector, timeAsmbSpaceTimeDtError, refCount);
                gsInfo << "t_{e/w} (|| Dt e ||)         = " <<  timeAsmbSpaceTimeDtError[refCount] << " sec.\n";
            }

#pragma omp section
            {
                gsSpaceTimeDtDeltaKNorm<real_t> eSpaceTimeDtDeltaK(v, * uForOMP[16]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeDtDeltaK, eSpaceTimeDtDistr, elemNumber,
                                                      eSpaceTimeDtDeltaKVector, timeAsmbSpaceTimeDtDeltaKError, refCount);
                gsInfo << "t_{e/w} (delta_K * || Dt e ||) = " <<  timeAsmbSpaceTimeDtDeltaKError[refCount] << " sec.\n";
                gsInfo << "delta_K * || Dt e || = " <<  eSpaceTimeDtDeltaKVector[refCount] << " sec.\n";
            }

#pragma omp section
            {
                gsSpaceTimeErrorIdentity<real_t> eIdentity(v, * fForOMP[0]);
                testSpaceTime.gsCalculateDistribution(eIdentity, eIdentDistr, elemNumber,
                                                      eIdentVector, timeAsmbSpaceTimeErrorIdentity, refCount);
                gsInfo << "t_{e/w} ( Id )               = " <<  timeAsmbSpaceTimeErrorIdentity[refCount] << " sec.\n";
            }
#pragma omp section
            {
                gsSeminormH1<real_t> eH1Seminorm(v, * uForOMP[5]);
                testSpaceTime.gsCalculateDistribution(eH1Seminorm, eH1Distr, elemNumber,
                                                      eH1Vector, timeAsmbH1Error, refCount);
                gsInfo << "t_{e/w} (| e |_H1)           = " <<  timeAsmbH1Error[refCount] << " sec.\n";
            }
#pragma omp section
            {
                gsNormL2<real_t> eL2Norm(v, * uForOMP[6]);
                testSpaceTime.gsCalculateDistribution(eL2Norm, eL2Distr, elemNumber,
                                                      eL2Vector, timeAsmbL2Error, refCount);
                gsInfo << "t_{e/w} (|| e ||_L2)         = " <<  timeAsmbL2Error[refCount] << " sec.\n";
            }
#pragma omp section
            {
                gsSpaceTimeLocalisedNorm<real_t> eSpaceTimeNorm(v, * uForOMP[7]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeNorm, eSpaceTimeDistr, elemNumber,
                                                      eSpaceTimeVector, timeAsmbSpaceTimeError, refCount);
                gsInfo << "t_{e/w} (|| e ||_{s, h, Q})  = " <<  timeAsmbSpaceTimeError[refCount] << " sec.\n";
            }
#pragma omp section
            {
                gsSpaceTimeSpaceGrad<real_t> eSpaceTimeSpaceGrad(v, * uForOMP[8]);
                testSpaceTime.gsCalculateDistribution(eSpaceTimeSpaceGrad, eSpaceTimeSpaceGradDistr, elemNumber,
                                                      eSpaceTimeSpaceGradVector, timeAsmbSpaceTimeGradSpaceError,
                                                      refCount);
                gsInfo << "t_{e/w} (|| grad_x e ||)     = " <<  timeAsmbSpaceTimeGradSpaceError[refCount] << " sec.\n";
            }
            //}
            //#pragma omp parallel sections
            //{
#pragma omp section
            {
                watches[0].restart();
                gsSpaceTimeSpaceGradSliceNorm<real_t> eSpaceTimeSpaceGradTopSliceNorm(solutionFieldVector[refCount], * uForOMP[10], spaceTimeAssemblerV.patches().interfaces(), topSides);
                gradxeTVector[refCount] = eSpaceTimeSpaceGradTopSliceNorm.compute();
                gsInfo << "t_{e/w} (|| grad_x e ||_T)   = " << watches[0].stop() << " sec.\n";
            }
#pragma omp section
            {
                watches[1].restart();
                gsSpaceTimeSpaceGradSliceNorm<real_t> eSpaceTimeSpaceGradBottomSliceNorm(solutionFieldVector[refCount], * uForOMP[11], spaceTimeAssemblerV.patches().interfaces(), bottomSides);
                gradxe0Vector[refCount] = eSpaceTimeSpaceGradBottomSliceNorm.compute();
                gsInfo << "t_{e/w} (|| grad_x e ||_0)   = " << watches[1].stop() << " sec.\n";
            }
#pragma omp section
            {
                watches[2].restart();
                gsSpaceTimeSliceNorm<real_t> eSpaceTimeTopSliceNorm(solutionFieldVector[refCount], * uForOMP[0], spaceTimeAssemblerV.patches().interfaces(), topSides);
                eTVector[refCount] = eSpaceTimeTopSliceNorm.compute();
                gsInfo << "t_{e/w} (|| e ||_T)          = " << watches[2].stop() << " sec.\n";
            }
#pragma omp section
            {
                watches[3].restart();
                gsSpaceTimeSliceNorm<real_t> eSpaceTimeBottomSliceNorm(solutionFieldVector[refCount], *uForOMP[12], spaceTimeAssemblerV.patches().interfaces(), bottomSides);
                e0Vector[refCount] = eSpaceTimeBottomSliceNorm.compute();
                gsInfo << "t_{e/w} (|| e ||_0)          = " << watches[3].stop() << " sec.\n";
            }
#pragma omp section
            {
                watches[20].restart();
                gsSpaceTimeLocalisedSigmaTNorm<real_t> eSigmaTNorm(v, * uForOMP[14], spaceTimeAssemblerV.patches().interfaces(), bottomSides);
                ehSigmaTVector[refCount] = eSigmaTNorm.compute();
                gsInfo << "t_{e/w} (|| e ||_{s,h, T})   = " << watches[20].stop() << " sec.\n";

            }
            ///*
#pragma omp section
            {
                watches[4].restart();
                gsSeminormH1<real_t> eWH1Seminorm(w, * uForOMP[12]);
                testSpaceTime.gsCalculateDistribution(eWH1Seminorm, eH1Distr, elemNumber, eWH1Vector, timeAsmbH1Error,
                                                      refCount);
                gsInfo << " t_{e/w} (| u - w |_H1)      = " << watches[4].stop() << " sec.\t || grad e_w ||_Q   = "<< eWH1Vector[refCount] << "\n";
            }
#pragma omp section
            {
                watches[5].restart();
                gsNormL2<real_t> eWL2Norm(w, * uForOMP[13]);
                testSpaceTime.gsCalculateDistribution(eWL2Norm, eL2Distr, elemNumber, eWL2Vector, timeAsmbL2Error, refCount);
                gsInfo << " t_{e/w} (|| u - w ||_L2)    = " << watches[5].stop() << " sec.\t" << "|| e_w ||_Q   = " << eWL2Vector[refCount] << "\n";

            }
            //*/
            //}
            //#pragma omp parallel sections
            //{

#pragma omp section
            {
                // Compute the residual between two successive iterations
                if (refCount > 0) {
                    gsField<> u_cur = solutionFieldVector[refCount];
                    gsMultiPatch<> u_cur_MP = solutionMPVector[refCount];
                    gsMultiPatch<> u_prev_MP = solutionMPVector[refCount-1];

                    gsMatrix<> u_curIntPoints = u_cur.igaFunction(0).basis().anchors();

                    // truncate the fixed Dirichlet BC node and leave only free coefficients
                    gsMatrix<> u_prevIntValsFree;
                    testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerV, u_prev_MP, u_curIntPoints,
                                                                        u_prevIntValsFree);
                    // reconstruct a multi-patch function from the obtained free coeffs.
                    gsMultiPatch<> u_prevIntMP;
                    spaceTimeAssemblerV.constructSolution(u_prevIntValsFree, u_prevIntMP);

                    //const gsGeometry<> &geom = u_cur_MP.piece(0);
                    //const gsMultiPatch<> mp(geom);
                    const gsMultiPatch<> mp(u_cur_MP.piece(0));

                    // reconstruct a field multi-patch function from the obtained free coeffs.
                    gsField<> u_prevInt(mp, u_prevIntMP);

                    watches[6].restart();
                    gsNormFields<real_t> ucur_minus_uprev(u_cur, u_prevInt);
                    relErrorVector[refCount - 1] = ucur_minus_uprev.compute(false, elemNumber);
                    gsInfo << " t_{e/w} (|| u_i - u_{i-1} ||_L2)    = " << watches[6].stop() << " sec.\n";
                }
            }
#pragma omp section
            {
                // Compute the residual between the current iteration and initial u_0
                if (refCount > 0) {
                    gsField<> u_cur = solutionFieldVector[refCount];

                    gsMultiPatch<> u_0_MP = solutionMPVector[0];
                    gsMultiPatch<> u_cur_MP = solutionMPVector[refCount];

                    gsMatrix<> u_curIntPoints = u_cur.igaFunction(0).basis().anchors();

                    // truncate the fixed Dirichlet BC node and leave only free coefficients
                    gsMatrix<> u_0IntValsFree;
                    testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerV, u_0_MP, u_curIntPoints,
                                                                        u_0IntValsFree);
                    // reconstruct a multi-patch function from the obtained free coeffs.
                    gsMultiPatch<> u_0IntMP;
                    spaceTimeAssemblerV.constructSolution(u_0IntValsFree, u_0IntMP);

                    //const gsGeometry<> & geom = u_cur_MP.piece(0);
                    //const gsMultiPatch<> mp(geom);
                    const gsMultiPatch<> mp(u_cur_MP.piece(0));

                    // reconstruct a field multi-patch function from the obtained free coeffs.
                    gsField<> u_0Int(mp, u_0IntMP);

                    watches[7].restart();
                    gsNormFields<real_t> ucur_minus_u0(u_cur, u_0Int);
                    relError0Vector[refCount-1] = ucur_minus_u0.compute(false, elemNumber);
                    gsInfo << " t_{e/w} (|| u_i - u_0 ||_L2)        = " << watches[7].stop() << " sec.\n";
                }
            }
        }

        // || e ||_{s, h}
        eFullSpaceTimeVector[refCount]          = math::sqrt(math::pow(eSpaceTimeVector[refCount], 2)
                                                             +  math::pow(ehSigmaTVector[refCount], 2));
        // || e ||_{L}
        eFullSpaceTimeSolOperVector[refCount]   = math::sqrt(math::pow(eSpaceTimeSolOperVector[refCount], 2)
                                                             + math::pow(gradxeTVector[refCount], 2));
        // Id
        eFullIdentVector[refCount]              = math::sqrt(math::pow(eIdentVector[refCount], 2)
                                                             + math::pow(gradxe0Vector[refCount], 2));
        // [e]
        eFullSpaceTimeSpaceGradVector[refCount] = math::sqrt(math::pow(eSpaceTimeSpaceGradVector[refCount], 2)
                                                             + math::pow(eTVector[refCount], 2));
        // || e ||_{s, loc}
        eFullSpaceTimeFullLocVector[refCount]       = (math::pow(eSpaceTimeSpaceGradVector[refCount], 2)
                                                             + math::pow(eTVector[refCount], 2)
                                                             + 0.5 * (math::pow(eSpaceTimeDtDeltaKVector[refCount], 2)
                                                                      - math::pow(eSpaceTimeDeltaxDeltaKVector[refCount], 2)));
        // || e ||_{s, loc} - 0.5 sum delta_K || Delta_x e||
        eFullSpaceTimeLocVector[refCount]       = math::sqrt(math::pow(eSpaceTimeSpaceGradVector[refCount], 2)
                                                             + math::pow(eTVector[refCount], 2)
                                                             + 0.5 * math::pow(eSpaceTimeDtDeltaKVector[refCount], 2));

        testSpaceTime.gsLogRefinementIterationErrorReport(refCount, hmaxVector, hminVector,
                                                          eL2Vector, eH1Vector,
                                                          eSpaceTimeSpaceGradVector, eFullSpaceTimeSpaceGradVector,
                                                          eSpaceTimeVector, eFullSpaceTimeVector, eFullSpaceTimeLocVector,
                                                          eSpaceTimeSolOperVector, eFullSpaceTimeSolOperVector,
                                                          eIdentVector, eFullIdentVector,
                                                          gradxe0Vector, gradxeTVector, e0Vector, eTVector, ehSigmaTVector);


        // Calculating norm || f || for some checks
        real_t fL2NormSq(0.0);
        gsFunctionExpr<> zeroFunc("0.0", (int) d);
        gsField<> zeroField(testSpaceTime.patches, zeroFunc);
        gsNormL2<real_t> fL2Norm(zeroField, fFunc);
        fL2NormSq = math::pow(fL2Norm.compute(), 2);

        testSpaceTime.gsRecontructMajorantBasedOnOptimalFlux(refCount, basisY, yDegree,
                                                             yVector, mpY, yDOFs,
                                                             mpV, v, mpW, w,
                                                             stopcritVector,
                                                             fFunc, uFunc, fL2NormSq,
                                                             hmaxVector, hminVector,
                                                             timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, timeSolvY,
                                                             timeAsmbMajorant, timeAsmbDeltaHMajorant, timeAsmbMajorantII,
                                                             majVector, mdVector, meqVector, majhVector, majIIVector, majIIGapVector,
                                                             mdDistr, mIIdDistr,
                                                             e0Vector,
                                                             elemNumber,
                                                             spaceTimeAssemblerV, topSides, bottomSides);

        if (refCount > 0) {
            thetaVector[refCount-1]    = eH1Vector[refCount] / eH1Vector[refCount-1];
            stopcritVector[refCount-1] = (1 - thetaVector[refCount-1]) * TOL * rho * relError0Vector[refCount-1];
        }

        // TODO: incapsulate it into a function with lots of outputs
        // testSpaceTime.gsLogRefinementIterationErrorReport(
        gsInfo << "\ntime for solving V-system \n";
        gsInfo << "\t with direct solver: " << timeSolvV(refCount, 0) << " sec.\n";
        gsInfo << "\t with linear solver: " << timeSolvV(refCount, 1) << " sec.\n\n";
        gsInfo << "\n%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << " Testing the majorant \n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "time for assembling DivDiv-system          : " << timeAsmbDivDivY[refCount] << " sec.\n";
        gsInfo << "time for assembling MM-system              : " << timeAsmbMMY[refCount] << " sec.\n";
        gsInfo << "time for assembling DivDiv- and MM-systems : " << timeAsmbY[refCount] << " sec.\n";
        gsInfo << "\ntime for solving Y-system: \n" <<
               "\t with direct solver : " << timeSolvY(refCount, 0) << " sec.\n" <<
               "\t with linear solver : " << timeSolvY(refCount, 1) << " sec.\n";


        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << " Testing the adanced majorant \n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "time for assembling W-system : " << timeAsmbW[refCount] << " sec.\n";
        gsInfo << "time for solving W-system    : \n" <<
               "\twith direct solver : " << timeSolvW(refCount, 0) << " sec.\n" <<
               "\twith linear solver : " << timeSolvW(refCount, 1) << " sec.\n";

        gsInfo << "time for element-wise evaluation of the majorant:           " << timeAsmbMajorant[refCount] << " sec.\n";
        gsInfo << "time for element-wise evaluation of the minorant:           " << timeAsmbMajorantII[refCount] << " sec.\n";
        gsInfo << "time for element-wise evaluation of the error:              " << timeAsmbSpaceTimeError[refCount] << " sec.\n";

        //! [Error and Residual Estimate Computation]

        // Log and plotToParaview the results
        testSpaceTime.gsLogRefinementIterationInfo(refCount, vDOFs, yDOFs, wDOFs,
                                                   eH1Vector, eL2Vector,
                                                   eFullSpaceTimeLocVector, eSpaceTimeSpaceGradVector, eFullSpaceTimeSpaceGradVector, eFullSpaceTimeSolOperVector,
                                                   relErrorVector, relError0Vector, thetaVector, stopcritVector,
                                                   majVector, majhVector, majIIVector, majIIGapVector, minVector,
                                                   etaVector, eFullIdentVector);

        if (refCount <= numTotalAdaptRef - 1) {
            testSpaceTime.gsSaveToFileRefinementIterationInfo(saveToFile,
                                                              v, spaceTimeAssemblerV.multiBasis(),
                                                              eSpaceTimeSpaceGradDistr, mdDistr,
                                                              eSpaceTimeSolOperDistr, eIdentDistr,
                                                              e0Vector, eTVector, gradxe0Vector, gradxeTVector,
                                                              refCount,
                                                              refCount, numTotalAdaptRef,
                                                              exampleNumber);
        }
        //! [Refine]
        if (refCount < numTotalAdaptRef - 1) {
            testSpaceTime.gsExecuteRefinement(spaceTimeAssemblerV,
                                              basisY, spaceTimeAssemblerW,
                                              basisYVector, basisWVector, mdDistr, //mdDistr, mdDistr, //mIIdDistr, // eIdentDistr, //eSpaceTimeSolOperDistr, //eH1Distr, //mdDistr, //eSpaceTimeDistr, //eSpaceTimeSolOperDistr, //eH1Distr, //eIdentDistr, //mdDistr,
                                              adaptRefCrit,
                                              markingParamTheta,
                                              refCount,
                                              yBasisRefDelay,
                                              wBasisRefDelay);
            spaceTimeAssemblerV.refresh();
            spaceTimeAssemblerW.refresh();

            // get new interpolation points for v, y, w from new reconstructed basises
            gsMatrix<> vInterpPoints = spaceTimeAssemblerV.multiBasis().basis(0).anchors();
            gsMatrix<> yInterpPoints = basisY.basis(0).anchors();
            gsMatrix<> wInterpPoints = spaceTimeAssemblerW.multiBasis().basis(0).anchors();
            gsMatrix<> vRefVector(1, 1), yRefVector(1, 1), wRefVector(1, 1);

            // evaluate new values of v, y, w based on the obtained interpolation points
            testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerV, mpV, vInterpPoints, vRefVector);
            testSpaceTime.gsSetVRefVector(vRefVector);

            yRefVector = (mpY.patch(0).eval(yInterpPoints)); // returns the [1 x N] vector that needed to be transposed
            yRefVector.transposeInPlace();
            testSpaceTime.gsSetYRefVector(yRefVector);

            testSpaceTime.interpolateToRefVectorWithDirichletBC(spaceTimeAssemblerW, mpW, wInterpPoints, wRefVector);
            testSpaceTime.gsSetWRefVector(wRefVector);
        }
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Refinement Iterations]

    testSpaceTime.gsLogTestResults(vDegree, yDegree, wDegree,
                                   m, l,
                                   yBasisRefDelay, wBasisRefDelay,
                                   markingParamTheta, numInitUniformRefV, numTotalAdaptRef,
                                   vDOFs, yDOFs, wDOFs,
                                   timeAsmbV, timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, timeAsmbW,
                                   timeSolvV, timeSolvY, timeSolvW,
                                   timeAsmbH1Error, timeAsmbSpaceTimeSolOperError, timeAsmbSpaceTimeDeltaxError, timeAsmbSpaceTimeDtError,
                                   timeAsmbMajorant, timeAsmbDeltaHMajorant, timeAsmbMajorantII, timeAsmbEtaIndicator, timeAsmbSpaceTimeErrorIdentity,
                                   eL2Vector, eH1Vector,
                                   eFullSpaceTimeLocVector, eSpaceTimeSpaceGradVector, eFullSpaceTimeSpaceGradVector, eFullSpaceTimeSolOperVector, eSpaceTimeDeltaxVector, eSpaceTimeDtVector,
                                   relErrorVector, relError0Vector,
                                   majVector, mdVector, meqVector, majhVector, majIIVector, majIIGapVector, minVector, etaVector, eIdentVector);


    testSpaceTime.gsSaveToFileTestResults(saveToFile, vDOFs, yDOFs, wDOFs,
                                          eL2Vector, eH1Vector, eSpaceTimeSpaceGradVector, eSpaceTimeVector, eSpaceTimeSolOperVector,
                                          majVector, majhVector, majIIVector, majIIGapVector, eIdentVector,
                                          numTotalAdaptRef, exampleNumber);

    gsInfo << "\nTotal execution time : " << clock_total.stop() << "\n";
    return 0;
}
