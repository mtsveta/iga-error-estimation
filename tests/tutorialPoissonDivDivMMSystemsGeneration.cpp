/** @file tutorialPoissonDivDivMMSystemsGeneration.cpp

    @brief Tutorial for generation of the DivDiv and MM matrices typical for construction of
    functional error estimates on the example of the Poisson-Dirichlet problem

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Matculevich
*/

#include <iostream>
#include <gismo.h>

#include <gsAssembler/gsAdaptiveRefUtils.h>
//#include <gsAssembler/gsErrEstPoissonResidual.h>
#include <gsErrorEstimates/gsTestMajorant.h>

using namespace gismo;

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

    //! [Initialize Test Parameters]
    // -------------------------------------------------------------------------------------------------------------- //
    // Define constants and preliminaries
    const int NUM_PATCHES = 1;     // All the geometries are single-patch geometries
    // Define stopwatch to measure the performance of the routines
    gsCPUStopwatch clock;

    // Define parameters of I/O
    bool plotToParaview = false; // Flag indicating whether objects must be plotted in ParaView
    bool saveToFile     = false; // Flag indicating whether objects must be plotted in ParaView
    bool isAdaptive     = false;

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
    // int exampleNumber(10);       // 10 example: 2d quater annulus, u = cos(x)*exp(y)
    // int exampleNumber(11);       // 11 example: 2d l-shape, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    //int exampleNumber(12);        // 12 example: 2d l-shape, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    // int exampleNumber(16);       // 16 example: 2d mix of quater annulus and l-shape domains, u = if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), (x^2+y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )
    const unsigned d(2);
    //*/

    /*
    //const unsigned exampleNumber(13);  // 3d unit cube: u = (1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2
    //const unsigned exampleNumber(14);  // 3d G-shape: u = cos(x)*exp(y)
    //const unsigned exampleNumber(18);    // 3d G-shape: u = tanh(1 - (x + 2*y + 4*z - 4))
    //const unsigned exampleNumber(19);    // 3d unit cube: u = tanh(1 - (x + 2*y + 4*z - 2))
    //const unsigned exampleNumber(20);    // 3d unit cube: u = sin(pi*x)*sin(pi*y)*sin(6*pi*x)
    //const unsigned exampleNumber(17);  // 3d unit cube: u = cos(x)*exp(y)
    const unsigned d(3);    // Dimension must be prescribed

    */
    gsTestMajorant<d> testMajorant(exampleNumber, isAdaptive);

    // Init the degree of basis S^{p, p}, where p = vDegree
    int vDegree(2), k(1);
    // Init the degree of the basis for the flux: y \in S^{p + k, p + k} + S^{p + k, p + k}, p = vDegree
    // yDegree must be >= 2 to ensure that S^{p + k, p + k} + S^{p + k, p + k} \subset H(\Omega, div§)
    int yDegree(vDegree + k);
    //int yBasisRefDelay(6);              // parameter for the delay in updating the flux

    // Refinement strategy
    int numTotalAdaptRef(10);    // Number of initial uniform  and total adaptive refinement steps
    //int numInitUniformRefV(0), numInitUniformRefY(0);
    //MarkingStrategy adaptRefCrit(GARU); // with alternatives GARU, PUCA, and BULK
    //real_t markingParamTheta(0.1);
    //MarkingStrategy adaptRefCrit(BULK); // with alternatives GARU, PUCA, and BULK
    //real_t markingParamTheta(0.4);

    /*
    if (saveToFile) {
        testMajorant.gsCreateResultsFolder(saveToFile, vDegree, yDegree, yBasisRefDelay,
                                           numTotalAdaptRef, adaptRefCrit, markingParamTheta);
    }
    */

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
    gsMultiBasis<> basisV, basisY;
    //testMajorant.gsGetInitialBasis(vDegree, yDegree, basisV, basisY, numInitUniformRefV, numInitUniformRefY);
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Basis]

    //! [Define Auxiliary Structures for Storing of the Results]
    // -------------------------------------------------------------------------------------------------------------- //
    // Initialize arrays to DOFs for v and y on each refinement step

    gsVector<index_t> vDOFs(numTotalAdaptRef), yDOFs(numTotalAdaptRef);
    // Initialize vectors with assembling and computation times on each refinement step
    gsVector<double> timeAsmbV(numTotalAdaptRef),
            timeAsmbDivDivY(numTotalAdaptRef), timeAsmbMMY(numTotalAdaptRef), timeAsmbY(numTotalAdaptRef), timeSolvY(numTotalAdaptRef),
            timeAsmbError(numTotalAdaptRef), timeAsmbMajorant(numTotalAdaptRef), timeAsmbEtaIndicator(numTotalAdaptRef);
    int numOfSolvers = 2;
    // matrix are used since we compare the performance of different solvers
    // [0] is the direct solver
    // [1] is the iterative solver
    gsMatrix<double> timeSolvV(numTotalAdaptRef, numOfSolvers);
    timeSolvV.setZero(numTotalAdaptRef, 2);

    gsVector<real_t> stopcritVector(numTotalAdaptRef-1);

    // Initialize vector of all the basis' for y to store the history along the refinement
    std::vector< gsMultiBasis<> > basisYVector;
    // -------------------------------------------------------------------------------------------------------------- //
    //! [Define Auxiliary Structures for Storaging the Results]

    gsPoissonPde<> poissonPde(testMajorant.patches, bcInfo, fFunc);
    gsPoissonAssembler<real_t> poissonAssembler(testMajorant.patches, basisV, bcInfo, fFunc);
    poissonAssembler.options().setInt("DirichletValues", dirichlet::l2Projection);
    gsMultiPatch<> mpV;

    //! [Refinement Iterations]
    // -------------------------------------------------------------------------------------------------------------- //
    for( int refCounter = 0; refCounter < numTotalAdaptRef ; refCounter++ )
    {
        //testMajorant.gsLogRefinementBasisInfo(refCounter, NUM_PATCHES, numTotalAdaptRef, poissonAssembler.multiBasis(NUM_PATCHES-1), basisY);

        poissonAssembler.refresh();

        //! [Assemble System of Discretized Poisson Equation]
        // ---------------------------------------------------------------------------------------------------------- //
        clock.restart();
        poissonAssembler.assemble();
        timeAsmbV[refCounter] = clock.stop();
        gsInfo << "time for assembling PDE: " << timeAsmbV[refCounter] << " sec.\n";

        if (saveToFile)
            testMajorant.gsSaveToFileKMatrix(true, poissonAssembler, refCounter);
        //! [Assemble System of Discretized Poisson Equation]

        //! [Solve System]
        // ---------------------------------------------------------------------------------------------------------- //
        gsMatrix<> vVector(1, 1);
        stopcritVector[refCounter-1] = poissonAssembler.multiBasis(0).basis(0).getMaxCellLength();
        testMajorant.gsSolveKhfhSystem(poissonAssembler, vVector, timeSolvV, refCounter, stopcritVector);
         //! [SolvevSystem]

        //! [Recover the Approximation Field]
        // Construct the solution as gsMultiPatch and gsField
        // ---------------------------------------------------------------------------------------------------------- //
        poissonAssembler.constructSolution(vVector, mpV);
        gsField<> v(poissonAssembler.patches(), mpV);
        vDOFs[refCounter] = poissonAssembler.multiBasis(NUM_PATCHES-1).getMapper((dirichlet::strategy) poissonAssembler.options().getInt("DirichletStrategy"),
                                                                                 (iFace::strategy) poissonAssembler.options().getInt("InterfaceStrategy"),
                                                                                 bcInfo, 0).size();
        //! [Recover the Approximation Field]

        //! [Dof mapper]
        // ---------------------------------------------------------------------------------------------------------- //
        /*
        gsDofMapper poissonMapper; // Gets the indices mapped from Basis --> Matrix
        poissonAssembler.multiBasis(NUM_PATCHES-1).getMapper((dirichlet::strategy) poissonAssembler.options().getInt("DirichletStrategy"),
                                                             (iFace::strategy) poissonAssembler.options().getInt("InterfaceStrategy"),
                                                             bcInfo, poissonMapper, 0);
        vDOFs[refCounter] = poissonMapper.size();
        */
        //! [Dof mapper]

        //! [Error, Majorant (Optimal Flux), and Residual Estimate Computation]
        // ---------------------------------------------------------------------------------------------------------- //
        // ---------------------------------------------------------------------------------------------------------- //

        // Initialtize the assembler for div(y) = f with free BC based on new basis for the flux y
        gsAssembler<> divdivAssembler;
        // Initizlize assembler with y = grad(v) with free BC
        gsAssembler<> dualAssembler;

        gsBoundaryConditions<> freeBC;
        gsPoissonPde<> divPde(testMajorant.patches, freeBC, fFunc);
        divdivAssembler.initialize(divPde, basisY);

        gsPoissonPde<> dualPde(testMajorant.patches, freeBC, mpV.piece(0));
        dualAssembler.initialize(dualPde, basisY);

        testMajorant.gsGenerateDivDivMMMatrices(refCounter, divPde, dualPde,
                                                basisY, yDegree, yDOFs,
                                                mpV, v, fFunc,
                                                divdivAssembler, dualAssembler,
                                                timeAsmbDivDivY, timeAsmbMMY, timeAsmbY, saveToFile);

        std::vector<real_t> mdDistr, edDistr, etaDistr;

        //! [Refine]
        if (refCounter < numTotalAdaptRef - 1) {
            /*
            testMajorant.gsExecuteRefinement(poissonAssembler, basisY, basisYVector, mdDistr, adaptRefCrit,
                                             markingParamTheta, clock,
                                             refCounter, yBasisRefDelay);
                                             */
            poissonAssembler.refresh();
        }
        //! [Refine]
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //! [Refinement Iterations]

    return 0;
}
