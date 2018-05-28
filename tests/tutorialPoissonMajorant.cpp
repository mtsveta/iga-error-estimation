/** @file tutorialPoissonMajorant.cpp

    @brief Example for using the gsPoissonSolver with adaptive refinement with THB-splines using
    functional error estimates as a refinimenent criterial.

    This file is poissonAssemblerrt of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#include <iostream>
#include <gismo.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsAssembler/gsErrEstPoissonResidual.h>
#include <gsErrorEstimates/gsTestMajorant.h>
#include <gsErrorEstimates/gsNormFields.h>

using namespace gismo;

// S. Matculevich
//
// This is a test example for a illustrating the adaptive refinement procedure implemented for the gsPoissonAssembler
// Flags, parameters, geometry and prescribed exact solution are specified within the main() function

bool gsParseCommandLine(int argc, char **argv, bool plot)
{
    //! [Parse command line]
    gsCmdLine cmd("Tutorial on solving a Poisson problem with guaranteed error control using the funcional error estimate.");
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
    real_t TOL = 1e-2;

    // Define parameters of I/O
    bool plotToParaview = false; // Flag indicating whether objects must be plotted in ParaView
    bool saveToFile     = false;  // Flag indicating whether objects must be plotted in ParaView
    bool isAdaptive     = true;

    if ( !gsParseCommandLine(argc, argv, plotToParaview) ) return 0;

    // Define test-case parameters
    ///*
    const unsigned exampleNumber(2);        // 2 example: 2d unit square, u = (1 - x)*x*x*(1 - y)*y
    // int exampleNumber(3);        // 3 example: 2d unit square, u = sin(pi*x)*sin(pi*y)
    // int exampleNumber(4);        // 4 example: 2d unit square, u = sin(6.0*pi*x)*sin(3.0*pi*y)
    // int exampleNumber(6);        // 6 example: 2d unit square, u = (x^2 - x) * (y^2 - y) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))
    // int exampleNumber(7);        // 7 example: 2d rectangle $(0, 2) \times (0, 1)$, u = (2 - x)*x*x*(1 - y)*y
    // int exampleNumber(8);        // 8 example: 2d rectangle $(0, 2) \times (0, 1)$, u = (x^2 - 2*x) * (y^2 - y) * exp(-100 * ((x - 1.4)^2 + (y - 0.95)^2))
    // int exampleNumber(9);        // 9 example: 2d l-shape, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    //int exampleNumber(10);       // 10 example: 2d quater annulus, u = cos(x)*exp(y)
    // int exampleNumber(11);       // 11 example: 2d l-shape, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // int exampleNumber(12);        // 12 example: 2d l-shape, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // int exampleNumber(16);       // 16 example: 2d mix of quater annulus and l-shape domains, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    const unsigned d(2);
    //*/

    /*
    const unsigned exampleNumber(13);  // 3d unit cube: u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2
    //const unsigned exampleNumber(14);  // 3d G-shape: u = cos(x)*exp(y)*10*z
    //const unsigned exampleNumber(18);    // 3d G-shape: u = tanh(1 - (x + 2*y + 4*z - 2))
    //const unsigned exampleNumber(19);    // 3d unit cube: u = tanh(1 - (x + 2*y + 4*z - 2))
    //const unsigned exampleNumber(20);    // 3d unit cube: u = sin(pi*x)*sin(pi*y)*sin(6*pi*x)
    //const unsigned exampleNumber(17);  // 3d unit cube: u = cos(x)*exp(y)
    const unsigned d(3);    // Dimension must be prescribed

    */
    gsTestMajorant<d> testMajorant(exampleNumber, isAdaptive);

    // Init the degree of basis S^{p, p}, where p = vDegree
    int vDegree(2), m(7), l(8);
    // Init the degree of the basis for the flux: y \in S^{p + k, p + k} + S^{p + k, p + k}, p = vDegree
    // yDegree must be >= 2 to ensure that S^{p + k, p + k} + S^{p + k, p + k} \subset H(\Omega, div§)
    int yDegree(vDegree + m);
    int wDegree(vDegree + l);

    // Refinement strategy
    int numInitUniformRefV(4), numInitUniformRefY(4), numInitUniformRefW(4), numTotalAdaptRef(6);    // Number of initial uniform  and total adaptive refinement steps
    //MarkingStrategy adaptRefCrit(GARU); // with alternatives GARU, PUCA, and BULK

    MarkingStrategy adaptRefCrit(BULK); // with alternatives GARU, PUCA, and BULK
    real_t markingParamTheta(0.4);
    int yBasisRefDelay(6);              // parameter for the delay in updating the flux for uniform and adaptive refinement
    int wBasisRefDelay(6);

    testMajorant.gsCreateResultsFolder(saveToFile,
                                       vDegree, yDegree, wDegree,
                                       yBasisRefDelay, wBasisRefDelay,
                                       numTotalAdaptRef, adaptRefCrit, markingParamTheta);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Initialize Test Parameters]

    //! [Define Problem Data]
    // -------------------------------------------------------------------------------------------------------------- //
    gsFunctionExpr<> uDFunc, fFunc, uFunc;
    gsBoundaryConditions<> bcInfo;
    testMajorant.gsInitializeProblemData(uDFunc, fFunc, uFunc, bcInfo);
    testMajorant.gsLogProblemData();
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Problem Data]

    //! [Define Basis]
    // -------------------------------------------------------------------------------------------------------------- //
    gsMultiBasis<> basisV, basisY, basisW;
    testMajorant.gsGetInitialBasis(vDegree, yDegree, wDegree,
                                   basisV, basisY, basisW,
                                   numInitUniformRefV,
                                   numInitUniformRefY,
                                   numInitUniformRefW);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Basis]

    //! [Define Auxiliary Structures for Storing of the Results]
    // -------------------------------------------------------------------------------------------------------------- //
    // Initialize arrays to store the error and estimators on each refinement step
    gsVector<real_t> eH1Vector(numTotalAdaptRef),
            relErrorVector(numTotalAdaptRef-1),     // delta_i = || u_i - u_{i-1} ||, i > 0
            relError0Vector(numTotalAdaptRef-1),    // eps0_i = || u_i - u_0 ||, i > 0
            thetaVector(numTotalAdaptRef-1),       // theta_i = || u_i - u || / || u_{i-1} - u ||, i > 0
            stopcritVector(numTotalAdaptRef-1),    // stop_crit_i = (1 - theta_i) * TOL * eps0_i * rho, i > 0
            majVector(numTotalAdaptRef), mdVector(numTotalAdaptRef), meqVector(numTotalAdaptRef),
            minVector(numTotalAdaptRef),
            etaVector(numTotalAdaptRef);
    // Initialize arrays to DOFs for v and y on each refinement step
    gsVector<index_t> vDOFs(numTotalAdaptRef), yDOFs(numTotalAdaptRef), wDOFs(numTotalAdaptRef);
    // Initialize vectors with assembling and computation times on each refinement step
    gsVector<real_t> timeAsmbV(numTotalAdaptRef),
            timeAsmbDivDivY(numTotalAdaptRef), timeAsmbMMY(numTotalAdaptRef), timeAsmbY(numTotalAdaptRef), timeAsmbW(numTotalAdaptRef),
            timeAsmbError(numTotalAdaptRef), timeAsmbMajorant(numTotalAdaptRef), timeAsmbMinorant(numTotalAdaptRef), timeAsmbEtaIndicator(numTotalAdaptRef),
            timeAsmbEquilY(numTotalAdaptRef);
    int numOfSolvers = 2;
    // matrix are used since we compare the performance of different solvers
    // [0] is the direct solver
    // [1] is the iterative solver
    gsMatrix<real_t> timeSolvV(numTotalAdaptRef, numOfSolvers),
            timeSolvY(numTotalAdaptRef, numOfSolvers),
            timeSolvW(numTotalAdaptRef, numOfSolvers),
            timeSolvEquilY(numTotalAdaptRef, numOfSolvers);
    timeAsmbMajorant.setZero(numTotalAdaptRef);
    timeSolvY.setZero(numTotalAdaptRef, 2);
    timeSolvV.setZero(numTotalAdaptRef, 2);
    timeSolvW.setZero(numTotalAdaptRef, 2);
    etaVector.setZero(numTotalAdaptRef);
    // Initialize auxiliary vector of all the basis' for y to store the history along the refinement
    std::vector< gsMultiBasis<> > basisYVector;
    std::vector< gsMultiBasis<> > basisWVector;

    std::vector< gsField<> > solutionFieldVector;
    std::vector< gsMultiPatch<> > solutionMPVector;

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Auxiliary Structures to Store the Results]
    gsPoissonPde<> poissonPde(testMajorant.patches, bcInfo, fFunc);
    gsPoissonAssembler<real_t> poissonAssembler(testMajorant.patches, basisV, bcInfo, fFunc);
    gsPoissonAssembler<real_t> poissonAssemblerW(testMajorant.patches, basisW, bcInfo, fFunc);
    poissonAssembler.options().setInt("DirichletValues", dirichlet::l2Projection);
    poissonAssemblerW.options().setInt("DirichletValues", dirichlet::l2Projection);

    gsMultiPatch<> mpV, mpY, mpW;
    gsMatrix<> vVector(1, 1), yVector(1, 1), wVector(1, 1);
    gsField<> w;
    gsField<> v;

    //! [Refinement Iterations]
    // -------------------------------------------------------------------------------------------------------------- //
    for( int refCount = 0; refCount < numTotalAdaptRef ; refCount++ )
    {
        testMajorant.gsLogRefinementBasisInfo(refCount, NUM_PATCHES, numTotalAdaptRef, poissonAssembler.multiBasis(NUM_PATCHES-1), basisY, poissonAssemblerW.multiBasis(NUM_PATCHES-1));

        poissonAssembler.refresh();

        //! [Assemble System of Discretized Poisson Equation]
        // ---------------------------------------------------------------------------------------------------------- //
        clock.restart();
        poissonAssembler.assemble();
        timeAsmbV[refCount] = clock.stop();
        gsInfo << "time for assembling V-system: " << timeAsmbV[refCount] << " sec.\n";
        clock.restart();
        //! [Assemble System of Discretized Poisson Equation]

        //! [Dof mapper]
        // ---------------------------------------------------------------------------------------------------------- //
        gsDofMapper poissonMapper; // Gets the indices mapped from Basis --> Matrix
        poissonAssembler.multiBasis(NUM_PATCHES-1).getMapper((dirichlet::strategy) poissonAssembler.options().getInt("DirichletStrategy"),
                                                             (iFace::strategy) poissonAssembler.options().getInt("InterfaceStrategy"),
                                                             bcInfo, poissonMapper, 0);
        vDOFs[refCount] = poissonMapper.size();
        //! [Dof mapper]

        //! [Solve System]
        // ---------------------------------------------------------------------------------------------------------- //
        testMajorant.gsSolveKhfhSystem(poissonAssembler, vVector, timeSolvV, refCount, stopcritVector);
        //! [SolvevSystem]

        //! [Recover the Approximation Field]
        // Construct the solution as gsMultiPatch and gsField
        //gsInfo << "vVector size: " << vVector.dim() << "\n";
        // ---------------------------------------------------------------------------------------------------------- //
        //gsInfo << "vVector : \n" << vVector << "\n";
        poissonAssembler.constructSolution(vVector, mpV);
        //gsInfo << "mpV.coefs():\n" <<  mpV.patch(0).coefs() << "\n";
        //gsInfo << "poissonAssembler mapper free dofs:\n" << poissonAssembler.system().colMapper(0).freeSize() << "\n";


        v = gsField<>(poissonAssembler.patches(), mpV);
        solutionFieldVector.push_back(v);
        solutionMPVector.push_back(mpV);

        //! [Recover the Approximation Field]

        //! [Error, Majorant (Optimal Flux), and Residual Estimate Computation]
        // ---------------------------------------------------------------------------------------------------------- //
        // ---------------------------------------------------------------------------------------------------------- //
        unsigned long elemNumber = poissonAssembler.multiBasis(0).basis(0).numElements();
        std::vector<real_t> mdDistr, majDistr, edDistr, etaDistr, minDistr;
        mdDistr.reserve(elemNumber);
        edDistr.reserve(elemNumber);
        etaDistr.reserve(elemNumber);
        minDistr.reserve(elemNumber);

        gsSeminormH1<real_t> eH1Seminorm(v, uFunc);
        gsErrEstPoissonResidual<real_t> etaIndicator(v, fFunc, bcInfo);

        // Calculating norm || f || for some checks
        real_t fL2NormSq(0.0);
        gsFunctionExpr<> zeroFunc("0.0", (int) d);
        gsField<> zeroField(testMajorant.patches, zeroFunc);
        gsNormL2<real_t> fL2Norm(zeroField, fFunc);
        fL2NormSq = math::pow(fL2Norm.compute(), 2);

        clock.restart();

        // commentede ut to avoid the error in cdash
        //omp_set_nested(2);

        //#pragma omp parallel sections
        {
            //#pragma omp section
            {
                // Majorant Estimate and Optimal Flux Computation
                testMajorant.gsRecontructMajorantBasedOnOptimalFlux(refCount, basisY, yDegree,
                                                                    yVector, mpY, yDOFs,
                                                                    mpV, v, stopcritVector,
                                                                    fFunc, fL2NormSq,
                                                                    timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, timeSolvY,
                                                                    timeAsmbMajorant,
                                                                    majVector, mdVector, meqVector, mdDistr, majDistr,
                                                                    elemNumber);
            }
            //#pragma omp section
            {
                // Minorant Estimate and Optimal W  Computation
                testMajorant.gsRecontructMinorantBasedOnAuxiliaryW(refCount, poissonAssemblerW, bcInfo,
                                                                   wVector, mpW, w, wDOFs, stopcritVector,
                                                                   v, fFunc,
                                                                   timeAsmbW, timeSolvW, timeAsmbMinorant,
                                                                   minVector, minDistr,
                                                                   elemNumber);

            }
        }
        #pragma omp parallel sections
        {
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
                    testMajorant.interpolateToRefVectorWithDirichletBC(poissonAssembler, u_prev_MP, u_curIntPoints,
                                                                       u_prevIntValsFree);
                    // reconstruct a multi-patch function from the obtained free coeffs.
                    gsMultiPatch<> u_prevIntMP;
                    poissonAssembler.constructSolution(u_prevIntValsFree, u_prevIntMP);

                    //const gsGeometry<> &geom = u_cur_MP.piece(0);
                    //const gsMultiPatch<> mp(geom);
                    const gsMultiPatch<> mp(u_cur_MP.piece(0));

                    // reconstruct a field multi-patch function from the obtained free coeffs.
                    gsField<> u_prevInt(mp, u_prevIntMP);

                    gsNormFields<real_t> ucur_minus_uprev(u_cur, u_prevInt);
                    relErrorVector[refCount - 1] = ucur_minus_uprev.compute(false, elemNumber);
                }
            }
            #pragma omp section
            {
                if (refCount > 0) {
                    gsField<> u_cur = solutionFieldVector[refCount];

                    gsMultiPatch<> u_0_MP = solutionMPVector[0];
                    gsMultiPatch<> u_cur_MP = solutionMPVector[refCount];

                    gsMatrix<> u_curIntPoints = u_cur.igaFunction(0).basis().anchors();

                    // truncate the fixed Dirichlet BC node and leave only free coefficients
                    gsMatrix<> u_0IntValsFree;
                    testMajorant.interpolateToRefVectorWithDirichletBC(poissonAssembler, u_0_MP, u_curIntPoints,
                                                                       u_0IntValsFree);
                    // reconstruct a multi-patch function from the obtained free coeffs.
                    gsMultiPatch<> u_0IntMP;
                    poissonAssembler.constructSolution(u_0IntValsFree, u_0IntMP);

                    //const gsGeometry<> & geom = u_cur_MP.piece(0);
                    //const gsMultiPatch<> mp(geom);
                    const gsMultiPatch<> mp(u_cur_MP.piece(0));

                    // reconstruct a field multi-patch function from the obtained free coeffs.
                    gsField<> u_0Int(mp, u_0IntMP);

                    gsNormFields<real_t> ucur_minus_u0(u_cur, u_0Int);
                    relError0Vector[refCount-1] = ucur_minus_u0.compute(false, elemNumber);
                }

            }
            #pragma omp section
            {
                // Compute the H1-seminorm of error
                testMajorant.gsCalculateDistribution(eH1Seminorm, edDistr, elemNumber,
                                                     eH1Vector, timeAsmbError, refCount);
            }
        }
        //#pragma omp parallel sections
        {
            //#pragma omp section
            /*
            {
                // Compute the residual error indicator
                testMajorant.gsCalculateDistribution(etaIndicator, etaDistr, elemNumber,
                                                     etaVector, timeAsmbEtaIndicator, refCount);
            }
            */
        }

        if (refCount > 0) {
            thetaVector[refCount-1]    = eH1Vector[refCount] / eH1Vector[refCount-1];
            stopcritVector[refCount-1] = (1 - thetaVector[refCount-1]) * TOL * rho * relError0Vector[refCount-1];
        }

                gsInfo << "\ntime for solving V-system \n";
        gsInfo << "\twith SLDLT solver: " << timeSolvV(refCount, 0) << " sec.\n";
        gsInfo << "\twith CG solver   : " << timeSolvV(refCount, 1) << " sec.\n\n";

        gsInfo << "\n%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << " Testing the majorant \n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "time for assembling DivDiv-system          : " << timeAsmbDivDivY[refCount] << " sec.\n";
        gsInfo << "time for assembling MM-system              : " << timeAsmbMMY[refCount] << " sec.\n";
        gsInfo << "time for assembling DivDiv- and MM-systems : " << timeAsmbY[refCount] << " sec.\n";
        gsInfo << "\ntime for solving Y-system: \n" <<
               "\twith SLDLT : " << timeSolvY(refCount, 0) << " sec.\n" <<
               "\twith CG    : " << timeSolvY(refCount, 1) << " sec.\n";


        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << " Testing the minorant \n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "time for assembling W-system : " << timeAsmbW[refCount] << " sec.\n";
        gsInfo << "time for solving W-system    : \n" <<
               "\twith SLDLT : " << timeSolvW(refCount, 0) << " sec.\n" <<
               "\twith CG    : " << timeSolvW(refCount, 1) << " sec.\n";

        gsInfo << "time for element-wise evaluation of the majorant:           " << timeAsmbMajorant[refCount] << " sec.\n";
        gsInfo << "time for element-wise evaluation of the minorant:           " << timeAsmbMinorant[refCount] << " sec.\n";
        gsInfo << "time for element-wise evaluation of the error:              " << timeAsmbError[refCount] << " sec.\n";
        gsInfo << "time for element-wise evaluation of the residual indicator: " << timeAsmbEtaIndicator[refCount] << " sec.\n";
        gsInfo << "time for element-wise evaluation of all residuals:          " << clock.stop() << "sec.\n";

        //! [Error and Residual Estimate Computation]


        // Log and plotToParaview the results
        testMajorant.gsLogRefinementIterationInfo(refCount, vDOFs, yDOFs, wDOFs,
                                                  eH1Vector, relErrorVector, relError0Vector, thetaVector, stopcritVector,
                                                  majVector, minVector, etaVector);

        if (refCount <= numTotalAdaptRef - 1) {
            testMajorant.gsSaveToFileRefinementIterationInfo(saveToFile,
                                                             v, poissonAssembler.multiBasis(),
                                                             edDistr, mdDistr, refCount);
        }
        //! [Refine]
        if (refCount < numTotalAdaptRef - 1) {
            testMajorant.gsExecuteRefinement(poissonAssembler,
                                             basisY, poissonAssemblerW,
                                             basisYVector, basisWVector,
                                             mdDistr,
                                             adaptRefCrit,
                                             markingParamTheta,
                                             refCount, yBasisRefDelay,
                                             wBasisRefDelay);
            poissonAssembler.refresh();
            poissonAssemblerW.refresh();

            gsMatrix<> vInterpPoints = poissonAssembler.multiBasis().basis(0).anchors(); // newCoeff
            gsMatrix<> yInterpPoints = basisY.basis(0).anchors();
            gsMatrix<> wInterpPoints = poissonAssemblerW.multiBasis().basis(0).anchors(); // newCoeff
            gsMatrix<> vRefVector(1, 1), yRefVector(1, 1), wRefVector(1, 1);

            testMajorant.interpolateToRefVectorWithDirichletBC(poissonAssembler, mpV, vInterpPoints, vRefVector);
            testMajorant.gsSetVRefVector(vRefVector);

            yRefVector = mpY.patch(0).eval(yInterpPoints);
            yRefVector.transposeInPlace();
            testMajorant.gsSetYRefVector(yRefVector);

            testMajorant.interpolateToRefVectorWithDirichletBC(poissonAssemblerW, mpW, wInterpPoints, wRefVector);
            testMajorant.gsSetWRefVector(wRefVector);
        }
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Refinement Iterations]

    testMajorant.gsLogTestResults(vDegree, yDegree, wDegree,
                                  m, l,
                                  yBasisRefDelay, wBasisRefDelay,
                                  markingParamTheta, numInitUniformRefV, numTotalAdaptRef,
                                  vDOFs, yDOFs, wDOFs,
                                  timeAsmbV, timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, timeAsmbW,
                                  timeSolvV, timeSolvY, timeSolvW,
                                  timeAsmbError, timeAsmbMajorant, timeAsmbMinorant, timeAsmbEtaIndicator,
                                  eH1Vector, relErrorVector, relError0Vector,
                                  majVector, mdVector, meqVector, minVector, etaVector);

    testMajorant.gsSaveToFileTestResults(saveToFile,
                                         vDOFs, yDOFs, wDOFs,
                                         eH1Vector, majVector, minVector, etaVector,
                                         numTotalAdaptRef);

    gsInfo << "\nTotal execution time : " << clock_total.stop() << "\n";
    return 0;
}
