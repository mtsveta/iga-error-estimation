//
// Created by Svetlana Matculevich on 26/02/2017.
//
#include <gsErrorEstimates/gsVisitorDivDivSpaceTime.h>
#include <gsErrorEstimates/gsVisitorDivDivSpaceTimeMh.h>
#include <gsErrorEstimates/gsVisitorDualPoissonSpaceTime.h>
#include <gsErrorEstimates/gsErrEstEquilSpaceTimeMajorant.h>
#include <gsErrorEstimates/gsErrEstDualSpaceTimeMajorant.h>
#include <gsErrorEstimates/gsErrEstDualSpaceTimeMajorantII.h>
#include <gsErrorEstimates/gsErrEstFvwSpaceTimeMajorantII.h>
#include <gsErrorEstimates/gsErrEstFvu0wSpaceTimeMajorantII.h>
#include <gsErrorEstimates/gsErrEstResidualSigmaTSpaceTimeMajorantII.h>
#include <gsErrorEstimates/gsErrEstDivDualSpaceTimeMajorant.h>
//#include <gsErrorEstimates/gsSpaceGradSpaceTimePde.h>
#include <gsErrorEstimates/gsErrEstMinorant.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <sys/stat.h>
#include <fstream>
#include <functional>
#include <gsSolver/gsIterativeSolver.h>
#include <gsSolver/gsSolverUtils.h>
#include <gsErrorEstimates/gsDivPde.h>

namespace gismo {

    template<unsigned d>
    class gsTestSpaceTimeMajorant {

    private:
        void initGivenProblemData();


    public:
        // Constructor
        gsTestSpaceTimeMajorant(const unsigned _exampleNumber,
                                bool _isAdaptive,
                                bool _withMajorant,
                                bool _withMajorantOptimization,
                                bool _withMajorantEquilibration) :
                exampleNumber(_exampleNumber),
                isAdaptive(_isAdaptive),
                withMajorant(_withMajorant),
                withMajorantOptimization(_withMajorantOptimization),
                withMajorantEquilibration(_withMajorantEquilibration)
        {
            dim = (int) d;
            initGivenProblemData();
        }

        void gsGetInitialBasis(int, int, int,
                               gsMultiBasis<> &, gsMultiBasis<> &, gsMultiBasis<> &,
                               int, int, int);

        void gsLogProblemData();

        void gsCreateResultsFolder(bool, int, int, int, int, int, int, int, MarkingStrategy, real_t);

        void gsLogRefinementBasisInfo(int, const int, int, gsMultiBasis<> &, gsMultiBasis<> & , gsMultiBasis<> & );

        void gsInitializeProblemData(gsFunctionExpr<> &, gsFunctionExpr<> &, gsFunctionExpr<> &, gsBoundaryConditions<> &);

        void gsSolveKhfhSystem(gsSpaceTimeAssembler<real_t> &,
                               gsMatrix<>&,
                               gsMatrix<double>&,
                               int,
                               gsVector<real_t>&);

        void gsSolveKhwhfhSystem(gsSpaceTimeAssembler<real_t> &,
                                 gsMatrix<> &,
                                 gsMatrix<> &,
                                 int,
                                 gsVector<real_t>& );

        void gsSolveMajorantOptimalSystem(gsSparseMatrix<real_t> &, gsMatrix<> &, gsMatrix<> &,
                                          gsMatrix<double> &, int, int, gsVector<index_t> &,
                                          gsVector<real_t> & );

        void gsRecontructMajorantBasedOnOptimalFlux(int, gsMultiBasis<real_t> &, int,
                                                    gsMatrix<real_t> &, gsMultiPatch<real_t> &, gsVector<index_t> &,
                                                    const gsMultiPatch<real_t> &, const gsField<> &,
                                                    const gsMultiPatch<real_t> &, const gsField<> &,
                                                    gsVector<real_t> &stopcritVector,
                                                    const gsFunctionExpr<real_t> &, const gsFunctionExpr<real_t> &, real_t,
                                                    gsVector<real_t> &, gsVector<real_t> &, real_t,
                                                    gsVector<double> &, gsVector<double> &, gsVector<double> &, gsMatrix<double> &, gsVector<double> &, gsVector<double> &, gsVector<double> &,
                                                    gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &,
                                                    std::vector<real_t> &,
                                                    std::vector<real_t> &,
                                                    const gsVector<real_t>& e0,
                                                    int,
                                                    gsSpaceTimeAssembler<real_t> &,
                                                    std::vector<patchSide> &, std::vector<patchSide> &);

/*
        void gsRecontructMajorantBasedOnEquilibratedFlux(int, gsMultiBasis<real_t> &, int,
                                                         gsMatrix<real_t> &, gsMultiPatch<real_t> &, gsVector<index_t> &,
                                                         gsMultiPatch<real_t> &, gsField<> &, gsFunctionExpr<real_t> &,
                                                         gsVector<double> &, gsMatrix<double> &, gsVector<double> &,
                                                         gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &,
                                                         std::vector<real_t> &, std::vector<real_t> &,
                                                         int);
*/
        void gsRecontructV(int, gsSpaceTimeAssembler<real_t> & , gsBoundaryConditions<> &,
                                    gsMatrix<real_t> &,
                                    gsMultiPatch<real_t> &,
                                    gsField<> &,
                                    gsVector<index_t> &,
                                    gsVector<real_t> &,
                                    gsVector<double> &,
                                    gsMatrix<double> &);
        void gsRecontructW(int, gsSpaceTimeAssembler<real_t> &, gsBoundaryConditions<> &,
                           gsMatrix<real_t> &,
                           gsMultiPatch<real_t> &,
                           gsField<> &,
                           gsVector<index_t> &,
                           gsVector<real_t> &,
                           gsVector<double> &,
                           gsMatrix<double> &);

        /*
        void gsRecontructMajorantIIBasedOnOptimalW(int refCounter, gsBoundaryConditions& bcInfo,
                                                   gsMultiBasis<real_t> &basisY, int yDegree,
                                                                               gsMatrix<real_t> &yVector, gsMultiPatch<real_t> &mpY,
                                                                               gsVector<index_t> &yDOFs,
                                                                               real_t beta,
                                                                               gsMultiPatch<real_t> &mpV, const gsField<> &v,
                                                   gsMultiBasis<real_t> &basisW, int wDegree, gsVector<index_t> &wDOFs, gsMultiPatch<real_t> &mpW, gsField<> &w,
                                                                               gsVector<real_t> &stopcritVector,
                                                                               const gsFunctionExpr<real_t> &fFunc, real_t fL2NormSq,
                                                                               gsVector<real_t> &hmaxVector, gsVector<real_t> &hminVector, real_t theta,
                                                                               gsVector<double> &timeAsmblGradSpaceW,
                                                                               gsVector<double> &timeAsmblGradTW,
                                                                               gsVector<double> &timeAsmblW,
                                                                               gsMatrix<double> &timeSolvW,
                                                                               gsVector<double> &timeAsmblMajII,
                                                                               gsVector<real_t> &majIIVector,
                                                                               gsVector<real_t> &mdIIVector,
                                                                               gsVector<real_t> &meqIIVector,
                                                                               std::vector<real_t> &mIIdDistr,
                                                                               int elemNum);

            */
            /*
             *
             *
             * int refCounter,
                                                                gsSpaceTimeAssembler<real_t>& poissonAssemblerW,
                                                                gsBoundaryConditions<> &bcInfo,
                                                                gsMatrix<real_t> &wVector,
                                                                gsMultiPatch<real_t> &mpW,
                                                                gsField<>& w,
                                                                gsVector<index_t> &wDOFs,
                                                                gsVector<real_t> &stopcritVector,
                                                                gsVector<double> &timeAsmblW,
                                                                gsMatrix<double> &timeSolvW
             */

        void gsGenerateDivMatrixRhs(int refCounter,
                                    const gsDivPde<real_t> &divPde,
                                    const gsMultiBasis<real_t> &basisY, int yDegree,
                                    gsVector<index_t> &yDOFs,
                                    gsAssembler<> &divdivAssembler,
                                    gsVector<double> &timeAsmblDivDivY);
        void gsGenerateDivMatrixTwoRhs(int, const gsDivPde<real_t> &, const gsMultiBasis<real_t> &, int,
                                    gsVector<index_t> &yDOFs,
                                    gsAssembler<> &divdivAssembler,
                                    gsVector<double> &timeAsmblDivDivY);
        void gsGenerateMMMatrixRhs(int refCounter,
                                  const gsPoissonPde<> &dualPde,
                                  const gsMultiBasis<real_t> &basisY, int yDegree,
                                  gsVector<index_t> &yDOFs,
                                  gsAssembler<> &dualAssembler,
                                  gsVector<double> &timeAsmblVectMassY);

        void gsExecuteRefinement(gsSpaceTimeAssembler<real_t> &,
                                 gsMultiBasis<real_t> &, gsSpaceTimeAssembler<real_t> &,
                                 std::vector<gsMultiBasis<> > &, std::vector<gsMultiBasis<> > &,
                                 std::vector<real_t> &,
                                 MarkingStrategy,
                                 real_t,
                                 unsigned int, unsigned int, unsigned int);

        void gsCalculateDistribution(gsNorm<real_t> &, std::vector<real_t> &, int,
                                     gsVector<real_t> &, gsVector<double> &,
                                     int);
        void gsCalculateDistribution(gsEstimator<real_t> &, std::vector<real_t> &, int,
                                     gsVector<real_t> &, gsVector<double> &,
                                     int);

        void gsLogTestResults(int vDegree, int yDegree, int wDegree,
                              int m, int l,
                              int yBasisRefDelay, int wBasisRefDelay,
                              real_t markingParamTheta, int numInitUniformRef, int numTotalAdaptRef,
                              gsVector<index_t> &vDOFs, gsVector<index_t> &yDOFs, gsVector<index_t> &wDOFs,
                              gsVector<double> &timeAsmbV, gsVector<double> &timeAsmbDivDivY, gsVector<double> &timeAsmbMMY, gsVector<double> &timeAsmbY, gsVector<double> &timeAsmbW,
                              gsMatrix<double> &timeSolvV, gsMatrix<double> &timeSolvY, gsMatrix<double> &timeSolvW,
                              gsVector<double> &timeAsmbH1Error, gsVector<double> &timeAsmbSpaceTimeSolOperError, gsVector<double> &timeAsmbSpaceTimeDeltaxError, gsVector<double> &timeAsmbSpaceTimeDtError,
                              gsVector<double> &timeAsmbMajorant, gsVector<double> &timeAsmbDeltaHMajorant, gsVector<double> &timeAsmbMajII,
                              gsVector<double> &timeAsmbEtaIndicator, gsVector<double> &timeAsmbSpaceTimeErrorIdentity,
                              gsVector<real_t> &eL2Vector, gsVector<real_t> &eH1Vector,
                              gsVector<real_t> &eSpaceTimeVector,
                              gsVector<real_t> &eSpaceTimeSpaceGradVector, gsVector<real_t> &eFullSpaceTimeSpaceGradVector,
                              gsVector<real_t> &eSpaceTimeSolOperVector, gsVector<real_t> &eSpaceTimeDeltaxVector, gsVector<real_t> &eSpaceTimeDtVector,
                              gsVector<real_t> &relErrorVector, gsVector<real_t> &relError0Vector,
                              gsVector<real_t> &majVector, gsVector<real_t> &mdVector, gsVector<real_t> &meqVector, gsVector<real_t> &majhVector, gsVector<real_t> &majIIVector, gsVector<real_t> &majIIGapVector,
                              gsVector<real_t> &minVector, gsVector<real_t> &etaVector, gsVector<real_t> &eIdentVector);

        void gsLogAssemblingSolvingTime(gsVector<double> &, gsVector<double> &, int);

        void gsLogRefinementIterationErrorReport(const int refCount, const real_t theta, const gsVector<real_t> & hmaxVector, const gsVector<real_t> & hminVector,
                                               const gsVector<real_t> & eL2Vector, const gsVector<real_t> & eH1Vector,
                                               const gsVector<real_t> & eSpaceTimeSpaceGradVector, const gsVector<real_t> & eFullSpaceTimeSpaceGradVector,
                                               const gsVector<real_t> & eSpaceTimeVector, const gsVector<real_t> & eFullSpaceTimeVector,
                                               const gsVector<real_t> & eSpaceTimeSolOperVector, const gsVector<real_t> & eFullSpaceTimeSolOperVector,
                                               const gsVector<real_t> & eIdentVector, const gsVector<real_t> & eFullIdentVector,
                                               const gsVector<real_t> & gradxe0Vector, const gsVector<real_t> & gradxeTVector,
                                               const gsVector<real_t> & e0Vector, const gsVector<real_t> & eTVector, const gsVector<real_t> & ehSigmaTVector) const;

        void gsLogRefinementIterationInfo(int refCounter, gsVector<index_t> &, gsVector<index_t> &, gsVector<index_t> &,
                                          gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &,
                                          gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &, gsVector<real_t> &);

        void gsLogGeometryBasisData(gsGeometry<> &,
                                    int, int, int, int);

        void gsSaveToFileRefinementIterationInfo(bool save,
                                                 const gsField<> &v,
                                                 gsMultiBasis<> &basis,
                                                 std::vector<real_t> edDistr,
                                                 std::vector<real_t> mdDistr,
                                                 std::vector<real_t> eSolDistr,
                                                 std::vector<real_t> idDistr,
                                                 gsVector<real_t> & e0Vector, gsVector<real_t> & eTVector, gsVector<real_t> & gradxe0Vector, gsVector<real_t> & gradxeTVector, int refCount,
                                                 int refCounter, int refTotal,
                                                 const unsigned exampleNumber);

        void gsSaveToFileTestResults(bool save,
                                     gsVector<index_t> &vDOFs, gsVector<index_t> &yDOFs, gsVector<index_t> &wDOFs,
                                     gsVector<real_t> &eL2Vector, gsVector<real_t> &eH1Vector,
                                     gsVector<real_t> &eSpaceTimeSpaceGradVector, gsVector<real_t> &eSpaceTimeVector,
                                     gsVector<real_t> &eSpaceTimeSolOperVector,
                                     gsVector<real_t> &majVector, gsVector<real_t> &majhVector, gsVector<real_t> &majIIVector,  gsVector<real_t> &,
                                     gsVector<real_t> &eIdentVector,
                                     int refTotal,
                                     const unsigned exampleNumber);

        void gsSaveToFileDivDivMMMatrices(bool, gsSparseMatrix<real_t> &, gsSparseMatrix<real_t> &,
                                          int refCounter, const unsigned exampleNumber);

        void gsSaveToFileKMatrix(bool, gsSpaceTimeAssembler<real_t> &,
                                 int refCounter, const unsigned exampleNumber);

        void gsSetVVector(gsMatrix<> &vector){ this->vVector = vector; }
        void gsSetVRefVector(gsMatrix<> &vector){ this->vRefVector = vector; }
        gsMatrix<> getVRefVector(){ return this->vRefVector; }

        void gsSetWVector(gsMatrix<> &vector){ this->wVector = vector; }
        void gsSetWRefVector(gsMatrix<> &vector){ this->wRefVector = vector; }
        gsMatrix<> getWRefVector(){ return this->wRefVector; }

        void gsSetYVector(gsMatrix<> & vector){ this->yVector = vector; }
        void gsSetYRefVector(gsMatrix<> & vector){ this->yRefVector = vector; }
        gsMatrix<> gsGetYRefVector(){ return this->yRefVector; }

        void gsInitializeBasis(int, gsMultiBasis<> &, int);

        void interpolateToRefVectorWithDirichletBC(const gsSpaceTimeAssembler<real_t> &,
                                                   const gsMultiPatch<> &,
                                                   gsMatrix<> &,
                                                   gsMatrix<> &);
        void interpolateToRefVectorDDim(gsMultiPatch<> & mpY,
                                        gsMatrix<> &,
                                        gsMatrix<> &);

    public:
        gsMultiPatch<> patches;
        int dim;
        real_t cFriedrichs;

        std::string domainName; // domain \Omega
        std::string uExpr;      // solution
        std::string fExpr;      // right-had side
        std::string uDExpr;     // Dirichlet BC

        std::string resultFolder;     // Dirichlet BC

        // Define test-case parameters
        const unsigned exampleNumber;
        const bool isAdaptive;
        const bool withMajorant;
        const bool withMajorantOptimization;
        const bool withMajorantEquilibration;


    private:
        gsMatrix<> vVector;
        gsMatrix<> vRefVector;

        gsMatrix<> wVector;
        gsMatrix<> wRefVector;

        gsMatrix<> yVector;
        gsMatrix<> yRefVector;


    };


    template <typename T>
    std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
    {
        assert(a.size() == b.size());

        std::vector<T> result;
        result.reserve(a.size());

        std::transform(a.begin(), a.end(), b.begin(),
                       std::back_inserter(result), std::plus<T>());
        return result;
    }
    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::initGivenProblemData() {

        const double PI = 3.14159265359;
        const int NUM_PATCHES = 1;

        // d is d_total = d_x + 1
        switch (d) {
            case 2:
                // --------------------------------------------------------------------------------
                //  Unit-square domains examples
                // --------------------------------------------------------------------------------
                if (exampleNumber == 2 || exampleNumber == 3 || exampleNumber == 4 ||
                        exampleNumber == 5 || exampleNumber == 6) {
                    real_t unit_side(1.0);
                    patches = *gsNurbsCreator<>::BSplineSquareGrid(NUM_PATCHES, NUM_PATCHES, unit_side);
                    cFriedrichs = 1.0 / (math::sqrt((real_t) (dim - 1)) * PI);
                    domainName = "unit interval";
                } else if (exampleNumber == 7 || exampleNumber == 8) {
                    // --------------------------------------------------------------------------------
                    //  Rectangular domains examples
                    // --------------------------------------------------------------------------------
                    real_t lx(2.0), ly(1.0), x0(0.0), y0(0.0);
                    patches = *gsNurbsCreator<>::BSplineRectangle(x0, y0, lx, ly);
                    cFriedrichs = lx / PI;   // cFriedrichs <= 1 / (pi * sqrt(l_1^{-2} + ... + l_n^{-2}));
                    domainName = "rectangle $(0, 2) \\times (0, 1)$";

                } else if (exampleNumber == 9){
                    // --------------------------------------------------------------------------------
                    //  Rectangular domains examples
                    // --------------------------------------------------------------------------------
                    real_t lx(1.0), ly(2.0), x0(0.0), y0(0.0);
                    patches = *gsNurbsCreator<>::BSplineRectangle(x0, y0, lx, ly);
                    cFriedrichs = lx / PI;   // cFriedrichs <= 1 / (pi * sqrt(l_1^{-2} + ... + l_n^{-2}));
                    domainName = "rectangle $(0, 1) \\times (0, 2)$";

                }

                break;

            case 3:
                if (exampleNumber == 13 || exampleNumber == 17 || exampleNumber == 19 || exampleNumber == 20 ||
                        exampleNumber == 21 || exampleNumber == 25) {
                    patches = *gsNurbsCreator<>::BSplineCube(2);    // here 2 is the degree
                    cFriedrichs = 1.0 / (math::sqrt((real_t) (dim - 1)) * PI);
                    domainName = "unit cube";

                } else if (exampleNumber == 11) {
                    // --------------------------------------------------------------------------------
                    //  Circle domains examples
                    // --------------------------------------------------------------------------------
                    real_t radius(1.0);
                    patches = *gsNurbsCreator<>::BSplineFatCircle(radius, 0.0, 0.0);
                    cFriedrichs = 1.0 / (math::sqrt((real_t) 1 / (real_t) 2) * PI);
                    domainName = "cirle";

                } else if (exampleNumber == 27) {
                    // --------------------------------------------------------------------------------
                    //  Unit square moving in time interval [0, 2]
                    // --------------------------------------------------------------------------------
                    real_t lx(1.0), ly(1.0), x0(0.0), y0(0.0);
                    real_t tT(1.0);
                    //patches = *gsNurbsCreator<>::BSplineCuboid(x0, y0, t0, lx, ly, tT);
                    //patches = *gsNurbsCreator<>::BSplineCube(x0, y0, t0, lx, ly, tT);
                    gsTensorBSpline<2,real_t>::uPtr geo2D = gsNurbsCreator<>::BSplineRectangle(x0, y0, lx, ly);
                    patches = * gsNurbsCreator<>::lift3D(*geo2D, tT);
                    //patches.uniformRefine();
                    cFriedrichs = 1.0 / (math::sqrt((real_t) (dim - 1)) * PI);
                    domainName = "unit square x [0, 2]";
                }
                else if (exampleNumber == 28) {
                    // --------------------------------------------------------------------------------
                    //  Rectangular [0, 2] x [0, 1] moving in time interval [0, 1]
                    // --------------------------------------------------------------------------------
                    real_t lx(2.0), ly(1.0), x0(0.0), y0(0.0);
                    real_t tT(1.0);
                    gsTensorBSpline<2,real_t>::uPtr geo2D = gsNurbsCreator<>::BSplineRectangle(x0, y0, lx, ly);
                    patches = * gsNurbsCreator<>::lift3D(*geo2D, tT);
                    cFriedrichs = 1.0 / (math::sqrt((real_t) 5 / (real_t) 4) * PI);
                    domainName = "[0, 2] x [0, 1] x [0, 1]";
                }
                else if (exampleNumber == 29) {
                    // --------------------------------------------------------------------------------
                    //  Rectangular [0, 2] x [0, 1] moving in time interval [0, 1]
                    // --------------------------------------------------------------------------------
                    real_t lx(2.0), ly(3.0), x0(0.0), y0(0.0);
                    real_t tT(1.0);
                    gsTensorBSpline<2,real_t>::uPtr geo2D = gsNurbsCreator<>::BSplineRectangle(x0, y0, lx, ly);
                    patches = * gsNurbsCreator<>::lift3D(*geo2D, tT);
                    cFriedrichs = 1.0 / (math::sqrt((real_t) 1 / (real_t) 4 + (real_t) 1 / (real_t) 9) * PI);
                    domainName = "[0, 2] x [0, 3] x [0, 1]";
                }
                else if (exampleNumber == 14 || exampleNumber == 18) {

                    // --------------------------------------------------------------------------------
                    //  G-shape domain (inscribed in unit cube)
                    // --------------------------------------------------------------------------------
                    std::string fileSrc = gsFileManager::find( "volumes/GshapedVolume.xml" );
                    patches = static_cast<gsMultiPatch<> > (gsReadFile<real_t>(fileSrc));
                    cFriedrichs = 1.0 / (math::sqrt((real_t) dim) * PI);
                    domainName = "G-shape domain";

                } else if (exampleNumber == 23 || exampleNumber == 24 || exampleNumber == 26 || exampleNumber == 30) {
                    // --------------------------------------------------------------------------------
                    //  Annulus domains examples
                    // --------------------------------------------------------------------------------
                    gsTensorBSpline<2,real_t>::uPtr geo2D =  gsNurbsCreator<>::BSplineFatQuarterAnnulus();
                    patches = * gsNurbsCreator<>::lift3D(*geo2D, 1.0);

                    cFriedrichs = 1.0 / (math::sqrt((real_t) 1 / (real_t) 2) * PI);   // cFriedrichs <= 1 / (pi * sqrt(l_1^{-2} + ... + l_n^{-2}));
                    domainName = "quarter-annulus lifted in time";
                }
                else if (exampleNumber == 15)
                {
                    // --------------------------------------------------------------------------------
                    //  3d L-shaped domain
                    // --------------------------------------------------------------------------------
                    real_t z(1.0);
                    gsTensorBSpline<2,real_t>::uPtr geo2D = gsNurbsCreator<>::BSplineLShape_p2C1();
                    patches = * gsNurbsCreator<>::lift3D(*geo2D, z);
                    cFriedrichs = 1.0 / (math::sqrt(1.0 / 2.0) * PI);
                    domainName = "L-shape, 3d";
                }
                break;

            default:
                gsInfo << "WARNING: The geometry has not been prescribed.\n";
        }

        switch (exampleNumber) {
            case 2:
                uExpr = "(1 - x)*x*x*(1 - y)*y";
                fExpr = "(1 - x)*x*x*(1 - 2*y)-(2 - 6*x)*(1 - y)*y";
                uDExpr = uExpr;

                break;

            case 3:
                uExpr = "sin(pi*x)*sin(pi*y)";
                fExpr = "pi*sin(pi*x)*(cos(pi*y) + pi*sin(pi*y))";
                uDExpr = uExpr;

                break;

            case 4:
                uExpr = "sin(6.0*pi*x)*sin(3.0*pi*y)";
                fExpr = "3.0*pi*sin(6.0*pi*x)*(cos(3.0*pi*y)+12.0*pi*sin(3.0*pi*y))";
                uDExpr = uExpr;

                break;

            case 5:
                uExpr = "cos(x)*exp(y)";
                fExpr = "2.0*cos(x)*exp(y)";
                uDExpr = uExpr;

                break;

            case 6:
                uExpr = "(x^2 - x) * (y^2 - y) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))";
                fExpr = "-(y^2 - y) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2)) * "
                        "("
                        "(-100) * 2 * (x - 0.8) * ((2*x - 1) + (x^2 - x) * (-100) * 2 * (x - 0.8))"
                        "+ (2 + (-100) * 2 * ((2 * x - 1) * (x - 0.8) + (x^2 - x)))"
                        ")"
                        "+ "
                        "(2*y - 1) * (x^2 - x) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))"
                        "+ (x^2 - x) * (y^2 - y) * (-100) * 2 * (y - 0.05) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))";
                uDExpr = uExpr;

                break;

            case 7:
                uExpr = "(2 - x)*x*x*(1 - y)*y";
                fExpr = "(2 - x)*x*x*(1 - 2*y)-(4 - 6*x)*(1 - y)*y";
                uDExpr = uExpr;

                break;

            case 8:
                uExpr = "(x^2 - 2*x) * (y^2 - y) * exp(-100 * ((x - 1.4)^2 + (y - 0.95)^2))";
                fExpr = "-exp(-100 * ((x - 1.4)^2 + (y - 0.95)^2)) * "
                        "((-100) * 2 * (x - 1.4) * ((2*x - 2) * (y^2 - y) + (x^2 - 2*x) * (y^2 - y) * (-100) * 2*(x - 1.4))"
                        "+ ((2) * (y^2 - y) + (2*x - 2) * (y^2 - y) * (-100) * 2 *(x - 1.4) + (x^2 - 2*x) * (y^2 - y) * (-100) * 2))"
                        "+ "
                        "(x^2 - 2*x) * (2*y - 1) * exp(-100 * ((x - 1.4)^2 + (y - 0.95)^2)) "
                        "+ (x^2 - 2*x) * (y^2 - y) * (-100) * 2 * (y - 0.05) * exp(-100 * ((x - 0.8)^2 + (y - 0.05)^2))";
                uDExpr = uExpr;

                break;

            case 9:
            case 11:
            case 12:
            case 16:
                uExpr = "if( y > 0, (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) - pi)/3.0 ), "
                        "           (x^2 + y^2)^(1.0/3.0) * sin( (2.0*atan2(y,x) + 3.0*pi)/3.0 ) )";
                fExpr = "if( y > 1, "
                        "    1.0/3.0 * ((x - 1)^2 + (y - 1)^2)^(-2.0/3.0) * 2*(y - 1)* sin( (2.0*atan2(y - 1, x - 1) - pi)/3.0 ) + "
                        "    ((x - 1)^2 + (y - 1)^2)^(1.0/3.0) * 2.0/3.0*cos( (2.0*atan2(y - 1, x - 1) - pi)/3.0 )*(x - 1)/((x - 1)^2 + (y - 1)^2), "
                        "    1.0/3.0 * ((x - 1)^2 + (y - 1)^2)^(-2.0/3.0) * 2*(y - 1)* sin( (2.0*atan2(y - 1, x - 1) + 3.0*pi)/3.0 ) + "
                        "    ((x - 1)^2 + (y - 1)^2)^(1.0/3.0) * 2.0/3.0*cos( (2.0*atan2(y - 1, x - 1) + 3.0*pi)/3.0 )*(x - 1)/((x - 1)^2 + (y - 1)^2) "
                        ")";
                uDExpr = uExpr;

                break;

            case 10:
                uExpr = "cos(x)*exp(y)";
                fExpr = "2.0*cos(x)*exp(y)";
                uDExpr = uExpr;

                break;

            case 13:
                uExpr = "(1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2";
                fExpr = "(1 - x)*x^2*(1 - y)*y^2*(2*z - 3*z^2)-("
                        "(2 - 6*x)*(1 - y)*y^2*(1 - z)*z^2 "
                        " + (1 - x)*x^2*(2 - 6*y)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;

                break;

            case 14:
                // G-shape domain
                uExpr = "cos(x)*exp(y)*10*(z - 0.5)";
                fExpr = "1.0";
                uDExpr = uExpr;

                break;

            case 17:
                // unit cube domain
                uExpr = "cos(x)*exp(y)*(z + 1)";
                fExpr = "1.0";
                uDExpr = uExpr;

                break;
            case 18:
                // G-shape domain
                uExpr = "tanh(1 - 100*(x + 2*y + 4*z - 4))";
                fExpr = "100000 * tanh(1 - 100*(x + 2*y + 4*z - 4)) * (1 - (tanh(1 - 100*(x + 2*y + 4*z - 4)))^2)"
                        "-400 * (1 - (tanh(1 - 100*(x + 2*y + 4*z - 4)))^2)";
                uDExpr = uExpr;

                break;
            case 19:
                // unit cube domain
                uExpr = "tanh(1 - (x + 2*y + 4*z - 2))";
                fExpr = "10.0 * tanh(1 - (x + 2*y + 4*z - 2)) * (1 - (tanh(1 - (x + 2*y + 4*z - 2)))^2)"
                        "- 3.0 * (1 - (tanh(1 - (x + 2*y + 4*z - 2)))^2)";
                uDExpr = uExpr;

                break;

            case 20:
                uExpr = "sin(pi*x)*sin(pi*y)*sin(pi*z)";
                fExpr = "pi*sin(pi*x)*sin(pi*y)*cos(pi*z) + 2*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)";
                uDExpr = uExpr;

                break;
            case 21:
                uExpr = "cos(x)*exp(y)*sin(pi*z)";
                fExpr = "cos(x)*exp(y)*pi*cos(pi*z)";
                uDExpr = uExpr;

                break;
            case 22:
                uExpr = "sin(pi*x)*(1 - y)*exp(3*y)";
                fExpr = "sin(pi*x)*(2 - 3*y)*exp(3*y) + pi*pi*sin(pi*x)*(1-y)*exp(3*y)";
                uDExpr = uExpr;

                break;
            case 23:
                uExpr = "(1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2";
                fExpr = "(1 - x)*x^2*(1 - y)*y^2*(2*z - 3*z^2)"
                        "-("
                        "(2 - 6*x)*(1 - y)*y^2*(1 - z)*z^2 "
                        " + (1 - x)*x^2*(2 - 6*y)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;

            case 24:
                uExpr = "sin(pi*x)*sin(pi*y)*sin(pi*z)";
                fExpr = "pi*sin(pi*x)*sin(pi*y)*cos(pi*z) + 2*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)";
                uDExpr = uExpr;

                break;

            case 25:
                uExpr = "(1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2";
                fExpr = "(1 - x)*x^2*(1 - y)*y^2*(2*z - 3*z^2)-("
                        "(2 - 6*x)*(1 - y)*y^2*(1 - z)*z^2 "
                        " + (1 - x)*x^2*(2 - 6*y)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;

            case 26:
                uExpr = "(1 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2";
                fExpr = "(1 - x)*x^2*(1 - y)*y^2*(2*z - 3*z^2)-("
                        "(2 - 6*x)*(1 - y)*y^2*(1 - z)*z^2 "
                        " + (1 - x)*x^2*(2 - 6*y)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;

            case 27:
                uExpr = "(1 - x)*x^2*(1 - y)*y^2*(2 - z)*z^2";
                fExpr = "(1 - x)*x^2*(1 - y)*y^2*(4*z - 3*z^2)-("
                        "(2 - 6*x)*(1 - y)*y^2*(2 - z)*z^2 "
                        " + (1 - x)*x^2*(2 - 6*y)*(2 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;

            case 28:
                uExpr = "(2 - x)*x^2*(1 - y)*y^2*(1 - z)*z^2";
                fExpr = "(2 - x)*x^2*(1 - y)*y^2*(2*z - 3*z^2)-("
                        "(4 - 6*x)*(1 - y)*y^2*(1 - z)*z^2 "
                        " + (2 - x)*x^2*(2 - 6*y)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;

            case 29:
                uExpr = "(2 - x)*x^2*(3 - y)*y^2*(1 - z)*z^2";
                fExpr = "(2 - x)*x^2*(3 - y)*y^2*(2*z - 3*z^2)-("
                        "(4 - 6*x)*(3 - y)*y^2*(1 - z)*z^2 "
                        " + (2 - x)*x^2*(6 - 6*y)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;

            case 30:
                uExpr = "(x^2 + y^2 - 1)*(x^2 + y^2 - 4)*(1 - z)*z^2";
                fExpr = "(x^2 + y^2 - 1)*(x^2 + y^2 - 4)*(2*z - 3*z^2)"
                        "-("
                        "2*(-5 + 6*x^2 + 2*y^2)*(1 - z)*z^2 "
                        " + 2*(-5 + 6*y^2 + 2*x^2)*(1 - z)*z^2"
                        ")";
                uDExpr = uExpr;
                break;
            default :
                gsInfo << "WARNING: The data functions were prescribed.\n";
        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsCreateResultsFolder(bool save, int exampleNumber,
                                                  int vDegree, int yDegree, int wDegree,
                                                  int yBasisRefDelay, int wBasisRefDelay,
                                                  int numTotalAdaptRef, MarkingStrategy adaptRefCrit,
                                                  real_t markingParamTheta) {
        if (save) {
            std::string folder = "space-time-example-" + util::to_string(exampleNumber)
                                 + (this->isAdaptive ? ("-adapt-marker-" + util::to_string(adaptRefCrit)
                                                        + "-theta-" + util::to_string(markingParamTheta * 100)
                                                        + "-M-" + util::to_string(yBasisRefDelay)
                                                        + "-L-" + util::to_string(wBasisRefDelay))
                                                     : ("-uniform-M-" + util::to_string(yBasisRefDelay)))
                                 + "-v-" + util::to_string(vDegree)
                                 + "-y-" + util::to_string(yDegree)
                                 + "-w-" + util::to_string(yDegree)
                                 + "-total-ref-" + util::to_string(numTotalAdaptRef)
                                 + "-theta-0" //"-theta-hmin-squared" //
                                 + "-based-on-ed"
                                 ;
            //+ "-based-on-eta";
            struct stat st = {0};
            if (stat(folder.c_str(), &st) == -1) gsFileManager::mkdir(folder.c_str());

            resultFolder = folder;
        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsLogProblemData() {
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "Problem statement: \n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "u      : " << this->uExpr << "\n";
        gsInfo << "f      : " << this->fExpr << "\n";
        gsInfo << "domain : " << this->domainName << "\n";
        gsInfo << "dim    : " << this->dim << "\n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";
        if (this->isAdaptive)
            gsInfo << " Adaptive refinemnet strategy \n";
        else
            gsInfo << " Uniform refinemnet strategy \n";
        gsInfo << "%-------------------------------------------------------------------------------------------------%\n";

    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsLogRefinementBasisInfo(int refCounter, const int numPatches, int numTotalAdaptRef,
                                                     gsMultiBasis<> & thbMultBasisV,
                                                     gsMultiBasis<> & thbMultBasisY,
                                                     gsMultiBasis<> & thbMultBasisW) {
        gsInfo << "\n%-----------------------------------------------------------------\n";
        gsInfo << "Refinement step " << refCounter << " of " << numTotalAdaptRef << "\n";
        gsInfo << "%-----------------------------------------------------------------\n\n";

        for (int i = 0; i < numPatches; ++i) {
            std::cout << "Patch # " << i << "\n";
            gsInfo << "thbMultBasisV info: " << thbMultBasisV.basis(i) << "\n";
            gsInfo << "thbMultBasisY info: " << thbMultBasisY.basis(i) << "\n";
            gsInfo << "thbMultBasisW info: " << thbMultBasisW.basis(i) << "\n";
        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsInitializeProblemData(gsFunctionExpr<> &uDFunc, gsFunctionExpr<> &fFunc,
                                                    gsFunctionExpr<> &uFunc,
                                                    gsBoundaryConditions<> &bcInfo) {
        uDFunc = gsFunctionExpr<>(this->uDExpr, (int) d);
        fFunc = gsFunctionExpr<>(this->fExpr, (int) d);
        uFunc = gsFunctionExpr<>(this->uExpr, (int) d);

        //! [Set the Dirichlet Boundary Conditions]
        // constructor of gsFunctionExpr by an expression string and the domain dimension (real function)
        /*
        for (gsMultiPatch<>::const_biterator bit = this->patches.bBegin(); bit != this->patches.bEnd(); ++bit) {
            bcInfo.addCondition(*bit, condition_type::dirichlet, &uDFunc);
            gsInfo << "bit : \n" << *bit << "\n";
        }
        */

        if (d == 2){
            bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &uDFunc);
            bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &uDFunc);
            //bcInfo.addCondition( boundary::north, condition_type::dirichlet, &uDFunc); Sigma_T = top of the cylinder
            bcInfo.addCondition( boundary::south, condition_type::dirichlet, &uDFunc);
        }else if (d == 3){
            bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &uDFunc);
            bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &uDFunc);
            bcInfo.addCondition( boundary::north, condition_type::dirichlet, &uDFunc);
            bcInfo.addCondition( boundary::south, condition_type::dirichlet, &uDFunc);
            bcInfo.addCondition( boundary::front, condition_type::dirichlet, &uDFunc);
            // bcInfo.addCondition( boundary::back, condition_type::dirichlet, &uDFunc); Sigma_T = top of the cylinder

        }
        //! [Set the Dirichlet Boundary Conditions]
    }

    template<unsigned d>
    void
    gsTestSpaceTimeMajorant<d>::gsCalculateDistribution(gsNorm<real_t> &residual,
                                                        std::vector<real_t> &resDistr,
                                                        int elemNum,
                                                        gsVector<real_t> &resValsVector,
                                                        gsVector<double> &timeAsmb,
                                                        int refCounter) {
        gsCPUStopwatch clock;
        // Compute the residual error indicator
        clock.restart();
        resValsVector[refCounter] = residual.compute(true, elemNum);
        resDistr = residual.elementNorms();
        timeAsmb[refCounter] = clock.stop();
    }
/*
    template<unsigned d>
    void
    gsTestSpaceTimeMajorant<d>::gsCalculateDistribution(gsEstimator<real_t> &residual,
                                                        std::vector<real_t> &resDistr,
                                                        int elemNum,
                                                        gsVector<real_t> &resValsVector,
                                                        gsVector<double> &timeAsmb,
                                                        int refCounter) {
        gsCPUStopwatch clock;
        // Compute the residual error indicator
        clock.restart();
        resValsVector[refCounter] = residual.compute(false);
        resDistr = residual.elementNorms();
        timeAsmb[refCounter] = clock.stop();
    }
*/
    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsRecontructMajorantBasedOnOptimalFlux(int refCounter,
                                                                            gsMultiBasis<real_t> &basisY, int yDegree,
                                                                            gsMatrix<real_t> &yVector, gsMultiPatch<real_t> &mpY,
                                                                            gsVector<index_t> &yDOFs,
                                                                            const gsMultiPatch<real_t> &mpV, const gsField<> &v,
                                                                            const gsMultiPatch<real_t> &mpW, const gsField<> &w,
                                                                            gsVector<real_t> &stopcritVector,
                                                                            const gsFunctionExpr<real_t> &fFunc, const gsFunctionExpr<real_t> &uFunc, real_t fL2NormSq,
                                                                            gsVector<real_t> &hmaxVector, gsVector<real_t> &hminVector, real_t theta,
                                                                            gsVector<double> &timeAsmblDivDivY,
                                                                            gsVector<double> &timeAsmblVectMassY,
                                                                            gsVector<double> &timeAsmblY,
                                                                            gsMatrix<double> &timeSolvY,
                                                                            gsVector<double> &timeAsmblMaj, gsVector<double> &timeAsmblDeltaHMaj, gsVector<double> &timeAsmblMajII,
                                                                            gsVector<real_t> &majVector,
                                                                            gsVector<real_t> &mdVector,
                                                                            gsVector<real_t> &meqVector,
                                                                            gsVector<real_t> &majhVector,
                                                                            gsVector<real_t> &majIIVector, gsVector<real_t> &majIIGapVector,
                                                                            std::vector<real_t> &mdDistr, std::vector<real_t> &mIIdDistr,
                                                                            const gsVector<real_t>& e0Vector,
                                                                            int elemNum,
                                                                            gsSpaceTimeAssembler<real_t> &spaceTimeAssembler,
                                                                            std::vector<patchSide> &topSides,
                                                                            std::vector<patchSide> &bottomSides) {

        gsCPUStopwatch clock_, clock;

        gsBoundaryConditions<> freeBC;
        gsBoundaryConditions<> * pfreeBC = new gsBoundaryConditions<>(freeBC);  // get a deep copy of the free BC

        const gsPiecewiseFunction<real_t> vPiece = mpV.piece(0);
        //const gsPiecewiseFunction<real_t> vPiece = mpW.piece(0); // Experiment of improving mEq in majorant

        // Initialize the BVP for the generating of the optimal y for the majorant
        const gsDivPde<real_t> divPde(this->patches, freeBC, fFunc, vPiece);    // div(y) = f equation
        const gsPoissonPde<> dualPde(this->patches, *pfreeBC, vPiece);          // y = grad(v)

        // Initizlize assemblers
        gsAssembler<> divdivAssembler, dualAssembler;
        divdivAssembler.initialize(divPde, basisY); // div(y) = f with free BC
        dualAssembler.initialize(dualPde, basisY);  // y = grad(v) with free BC

        gsInfo << "---------------------------------------------------------------------------------\n";
        gsInfo << " Assembling optimal system for the majorant \n";
        gsInfo << "---------------------------------------------------------------------------------\n";

        clock.restart();
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                this->gsGenerateMMMatrixRhs(refCounter, dualPde, basisY, yDegree, yDOFs, dualAssembler, timeAsmblVectMassY);
            }
            #pragma omp section
            {
                this->gsGenerateDivMatrixRhs(refCounter, divPde, basisY, yDegree, yDOFs, divdivAssembler, timeAsmblDivDivY);
            }
        }
        timeAsmblY[refCounter] = clock.stop();

        gsSparseMatrix<real_t> divdivM = divdivAssembler.matrix();
        gsSparseMatrix<real_t> vectmassM = dualAssembler.matrix();
        gsMatrix<real_t> divdivRhs = divdivAssembler.rhs();
        gsMatrix<real_t> vectmassRhs = dualAssembler.rhs();

        this->gsSaveToFileDivDivMMMatrices(false, divdivM, vectmassM, refCounter, exampleNumber);

        // initialize components related to the maj and majh
        real_t mEq(0.0), mD(0.0), maj(0.0), divmD(0.0), majh(0.0), majUpwindSq(0.0), majSq(0.0);
        real_t beta(1.0), alpha(1.0);
        // initialize components related to the majII
        real_t mIIEq(0.0), mIID(0.0), mIIF(0.0), majII(0.0), mIIT(0.0), mII0(0.0);
        real_t majIIGap(1.0);
        real_t betaII(1.0);

        index_t currentSize = dualAssembler.matrix().innerSize();       // bigSize
        index_t newSize = currentSize * (this->dim - 1) / this->dim;    // smallSize
        yVector.setZero(currentSize, 1);                                // bigVector
        gsMatrix<real_t> yVector_(newSize, 1);                          // smallVector
        //divdivAssembler.constructSolution(yVector, mpY);
        dualAssembler.constructSolution(yVector, mpY);


        real_t ratioDualEq(10.0);                   // ratio for exiting criterion form the loop
        int iterMajOpt = 1;                         // number of loops to update beta and flux

        real_t h(hminVector[refCounter]);           // h_min = h_max in
        real_t delta_h = theta * h;

        gsInfo << "---------------------------------------------------------------------------------\n";
        gsInfo << " Solving optimal system for the majorant \n";
        gsInfo << "---------------------------------------------------------------------------------\n";


        if (this->withMajorant && this->withMajorantOptimization) {
            for (index_t i = 0; i < iterMajOpt; i++) {
                /*
                // old version
                // all d blocks (including d-th zero block are sent to the solver)
                // blocks are counter from 0 to d
                gsSparseMatrix<real_t> yM = dualAssembler.matrix() + math::pow(this->cFriedrichs, 2) / beta * divdivAssembler.matrix();
                gsMatrix<real_t>     yRhs = dualAssembler.rhs() - math::pow(this->cFriedrichs, 2) / beta * divdivAssembler.rhs();
                this->gsSolveMajorantOptimalSystem(yM, yRhs, yVector, timeSolvY, refCounter, i, yDOFs,
                                                   stopcritVector);
                */
                // new version
                // at this stage we can cut d-th block of yM and yRhs and just store the system responsible for 0, ..., d-1 blocks
                // solve it for 0, ..., d-1 blocks and add 0 block into yVector after gsSolveMajorantOptimalSystem procedure
                ///*

                gsSparseMatrix<real_t> yM_ = dualAssembler.matrix().block(0, 0, newSize, newSize)
                                             + math::pow(this->cFriedrichs, 2) / beta *
                                               divdivAssembler.matrix().block(0, 0, newSize, newSize);
                gsMatrix<real_t> yRhs_ = dualAssembler.rhs().block(0, 0, newSize, 1)
                                         - math::pow(this->cFriedrichs, 2) / beta *
                                           divdivAssembler.rhs().block(0, 0, newSize, 1);
                this->gsSolveMajorantOptimalSystem(yM_, yRhs_, yVector_, timeSolvY, refCounter, i, yDOFs,
                                                   stopcritVector);

                //gsInfo << "currentSize = " << currentSize << "\n";
                //gsInfo << "newSize = " << newSize << "\n";
                //gsInfo << "yDOFs[refCounter] = " << yDOFs[refCounter] << "\n";

                yVector.setZero(currentSize, 1);
                yVector.block(0, 0, newSize, 1) = yVector_;
                //*/
                /*
                index_t currentSize = dualAssembler.matrix().innerSize();
                index_t newSize = currentSize * (this->dim - 1) / this->dim; // oldSize
                gsInfo << "currentSize = " << currentSize << "\n";
                gsInfo << "newSize = " << newSize << "\n";
                gsInfo << "yDOFs[refCounter] = " << yDOFs[refCounter] << "\n";
                */
                yVector.resize(yDOFs[refCounter], this->dim);
                this->gsSetYRefVector(yVector);
                divdivAssembler.constructSolution(yVector, mpY);
                const gsField<real_t> y(divdivAssembler.patches(), mpY);

                /*
                gsMatrix<real_t> tmpMeq;
                clock_.restart();
                tmpMeq = yVector.transpose() * divdivM.toDense() * yVector + 2.0 * divdivRhs.transpose() * yVector;
                mEq_ = math::abs(tmpMeq(0, 0) + fL2NormSq);
                timeMeq += clock_.stop();

                clock_.restart();
                tmpMd = yVector.transpose() * vectmassM.toDense() * yVector - 2.0 * vectmassRhs.transpose() * yVector;
                mD_ = math::abs(tmpMd(0, 0) + vH1SeminormSq);
                timeMd += clock_.stop();

                // TODO: || grad v||^2 is diverging with respect to || grad u ||^2 if one compares these values
                gsInfo << "vH1SeminormSq = " << vH1SeminormSq << "\n";
                gsInfo << "uH1SeminormSq = " << uH1SeminormSq << "\n";
                gsInfo << "mD2_   = " << mD_ << "\n";
                gsInfo << "mEq2_  = " << mEq_ << "\n";

                gsInfo << "time for matrix evaluation of the mD: " << timeMd << " sec.\n";
                gsInfo << "time for matrix evaluation of the mEq: " << timeMeq << " sec.\n";

                */
                /*
                //divdivAssemblerMh.constructSolution(yVectorMh, mpYMh);
                //const gsField<real_t> yMh(divdivAssemblerMh.patches(), mpYMh);

                // Calculate the terms of the majorant separately
                //gsErrEstDivDualSpaceTimeMajorant<real_t> divmdIndicator(v, mpY);
                //gsErrEstDualSpaceTimeMajorantII<real_t> mIIdIndicator(v, y, w);
                //gsErrEstEquilSpaceTimeMajorant<real_t> meqIndicator(v, y, fFunc);
                //gsErrEstEquilSpaceTimeMajorant<real_t> mIIeqIndicator(w, y, fFunc);
                //gsErrEstFvwSpaceTimeMajorantII<real_t> mIIFIndicator(v, w, fFunc);
                */
                /*
                gsErrEstDivDualSpaceTimeMajorant<real_t> divmdMhIndicator(v, mpYMh);
                gsErrEstDualSpaceTimeMajorant<real_t> mdMhIndicator(v, mpYMh);
                gsErrEstEquilSpaceTimeMajorant<real_t> meqMhIndicator(v, yMh, fFunc);
                */

                gsFunctionExpr<> *fP1 = new gsFunctionExpr<>(fFunc);
                gsFunctionExpr<> *fP2 = new gsFunctionExpr<>(fFunc);
                gsCPUStopwatch clock_1, clock_2, clock_3, clock_4, clock_5, clock_6, clock_7;

                real_t timeAsmblMd(0.0), timeAsmblMeq(0.0), timeAsmblMIId(0.0), timeAsmblMIIeq(0.0);

#pragma omp parallel sections
                {
#pragma omp section
                    {
                        gsErrEstDualSpaceTimeMajorant<real_t> mdIndicator(v, mpY);
                        clock_2.restart();
                        mD = mdIndicator.compute(true, elemNum);
                        timeAsmblMd += clock_2.stop();
                        gsInfo << "time for element-wise evaluation of the mD: " << clock_2.stop() << " sec.\n";
                        mdDistr = mdIndicator.elementNorms();
                    }
#pragma omp section
                    {
                        if (theta != 0.0) {
                            gsErrEstDivDualSpaceTimeMajorant<real_t> divmdIndicator(v, mpY);
                            clock.restart();
                            divmD = divmdIndicator.compute(false, elemNum);
                            timeAsmblMaj[refCounter] += clock.stop();
                            gsInfo << "time for element-wise evaluation of the div(mD): " << clock.stop() << " sec.\n";

                            timeAsmblDeltaHMaj[refCounter] += clock.stop();
                        }
                    }

#pragma omp section
                    {
                        gsErrEstEquilSpaceTimeMajorant<real_t> meqIndicator(v, y, fFunc);
                        clock_1.restart();
                        mEq = meqIndicator.compute(true, elemNum);
                        timeAsmblMeq += clock_1.stop();
                        gsInfo << "time for element-wise evaluation of the mEq: " << timeAsmblMeq << " sec.\n";
                    }

#pragma omp section
                    {
                        gsErrEstResidualSigmaTSpaceTimeMajorantII<real_t> mIITIndicator(v, w,
                                                                                        spaceTimeAssembler.patches().interfaces(),
                                                                                        topSides);
                        clock_3.restart();
                        mIIT = mIITIndicator.compute(false);
                        gsInfo << "time for element-wise evaluation of the mIIT: " << clock_3.stop() << " sec.\n";
                    }
                    //*/
                    ///*
#pragma omp section
                    {
                        gsErrEstFvu0wSpaceTimeMajorantII<real_t> mII0Indicator(v, w, uFunc,
                                                                               spaceTimeAssembler.patches().interfaces(),
                                                                               bottomSides);
                        clock_4.restart();
                        mII0 = mII0Indicator.compute(false);
                        gsInfo << "time for element-wise evaluation of the mII0: " << clock_4.stop() << " sec.\n";
                    }
                    //*/
#pragma omp section
                    {
                        gsErrEstDualSpaceTimeMajorantII<real_t> mIIdIndicator(v, y, w);
                        clock_5.restart();
                        mIID = mIIdIndicator.compute(false, elemNum);
                        gsInfo << "time for element-wise evaluation of the mIID: " << clock_5.stop() << " sec.\n";
                        mIIdDistr = mIIdIndicator.elementNorms();
                    }
#pragma omp section
                    {
                        gsErrEstEquilSpaceTimeMajorant<real_t> mIIeqIndicator(w, y, *fP1);
                        clock_6.restart();
                        mIIEq = mIIeqIndicator.compute(false, elemNum);
                        gsInfo << "time for element-wise evaluation of the mIIEq: " << clock_6.stop() << " sec.\n";
                    }
#pragma omp section
                    {
                        gsErrEstFvwSpaceTimeMajorantII<real_t> mIIFIndicator(v, w, *fP2);
                        clock_7.restart();
                        mIIF = mIIFIndicator.compute(true, elemNum);
                        gsInfo << "time for element-wise evaluation of the mIIF: " << clock_7.stop() << " sec.\n";
                    }
                }
                timeAsmblMaj[refCounter] += (timeAsmblMd > timeAsmblMeq ? timeAsmblMd : timeAsmblMeq);
                gsInfo << "increment time for element-wise evaluation of the majorant: " << clock.stop() << " sec.\n";

                timeAsmblMajII[refCounter] += (timeAsmblMIId > timeAsmblMeq ? timeAsmblMIId : timeAsmblMIIeq);
                gsInfo << "time for element-wise evaluation of the majII: " << timeAsmblMajII[refCounter] << " sec.\n";

                // Update beta and alpha
                beta = this->cFriedrichs * mEq / mD;
                betaII = this->cFriedrichs * mIIEq / mIID;
                if (theta != 0.0) alpha = mEq / divmD;

                majSq = math::pow(e0Vector[refCounter], 2) + (1 + beta) * math::pow(mD, 2) +
                        (1 + 1 / beta) * math::pow(cFriedrichs * mEq, 2);
                majUpwindSq = delta_h * ((1 + alpha) * math::pow(divmD, 2) + (1 + 1 / alpha) * math::pow(mEq, 2));

                // Update error estimate majI
                maj = math::sqrt(majSq);

                // Update error estimate majh
                majh = math::sqrt(majSq + majUpwindSq);

                // Update error estimate majII
                majII = math::sqrt(math::pow(mIIT, 2) + 2 * mIIF + mII0
                                   + (1 + betaII) * math::pow(mIID, 2) +
                                   (1 + 1 / betaII) * math::pow(cFriedrichs * mIIEq, 2));
                majIIGap = sqrt(4 * (1 + betaII) - 2);
                std::cout << std::scientific;
                gsInfo << "iter " << i << ": \t"
                       << "maj   = " << maj << "\t"
                       << "mD  = " << mD << "\t"
                       << "mEq  = " << mEq << "\t"
                       << "maj0   = " << e0Vector[refCounter] << "\t";

                //gsInfo << "majSq = " << majSq << "\n";
                //gsInfo << "majUpwindSq = " << majUpwindSq << "\n";
                //gsInfo << "divmD = " << divmD << "\n";
                gsInfo << "beta  = " << beta << "\n";
                gsInfo << "iter " << i << ": \t" << "majh  = " << majh << "\t" << "divmD = " << divmD << "\n";
                gsInfo << "iter " << i << ": \t"
                       << "majII = " << majII << "\t"
                       << "mIID = " << mIID << "\t"
                       << "mIIEq = " << mIIEq << "\t"
                       << "mIIF = " << mIIF << "\t"
                       << "mII0 = " << mII0 << "\t"
                       << "mIIT = " << mIIT << "\t"
                       << "betaII = " << betaII << "\n";
                gsInfo << "comparison of maj = " << maj << " and majII = " << majII << "\n";
                gsInfo << "majIIGap = " << majIIGap << "\n";
                gsInfo << "comparison of maj = " << maj << " and majII (corrected with gap) = " << majII / majIIGap
                       << "\n\n";

                gsInfo << "h     = " << hmaxVector[refCounter] << "\n";
                gsInfo << "theta = " << theta << "\n";
                gsInfo << "beta  = " << beta << "\n\n";

                /*
                gsInfo << " -------------------------------------------------------------- \n";
                gsInfo << "iter " << i << ": \t" << "majMh  = " << majMh << "\t"
                       << "majDMh  = " << mDMh << "\t" << "majEqMh  = " << mEqMh << "\n";
                gsInfo << "majSqMh = " << majSqMh << "\n";
                gsInfo << "majUpwindSqMh = " << majUpwindSqMh << "\n";
                gsInfo << "divmDMh = " << divmDMh << "\n";
                gsInfo << " -------------------------------------------------------------- \n";
                */
                std::cout.unsetf(std::ios::scientific);

                // If the ratio between m_d and mEq is more then certain threshhold, do no continue to optimize
                if (mD / mEq >= ratioDualEq) iterMajOpt = i;
            }
            if (!this->withMajorantOptimization) {
                gsSparseMatrix<real_t> yM_ = dualAssembler.matrix().block(0, 0, newSize, newSize);
                gsMatrix<real_t> yRhs_ = dualAssembler.rhs().block(0, 0, newSize, 1);
                this->gsSolveMajorantOptimalSystem(yM_, yRhs_, yVector_, timeSolvY, refCounter, 0, yDOFs,
                                                   stopcritVector);

                yVector.setZero(currentSize, 1);
                yVector.block(0, 0, newSize, 1) = yVector_;
                yVector.resize(yDOFs[refCounter], this->dim);
                this->gsSetYRefVector(yVector);

                dualAssembler.constructSolution(yVector, mpY);
                const gsField<real_t> y(dualAssembler.patches(), mpY);

            }
            if (!this->withMajorantEquilibration) {
                gsSparseMatrix<real_t> yM_ = divdivAssembler.matrix().block(0, 0, newSize, newSize);
                gsMatrix<real_t> yRhs_ = - divdivAssembler.rhs().block(0, 0, newSize, 1);
                this->gsSolveMajorantOptimalSystem(yM_, yRhs_, yVector_, timeSolvY, refCounter, 0, yDOFs,
                                                   stopcritVector);

                yVector.setZero(currentSize, 1);
                yVector.block(0, 0, newSize, 1) = yVector_;

                yVector.resize(yDOFs[refCounter], this->dim);
                this->gsSetYRefVector(yVector);
                divdivAssembler.constructSolution(yVector, mpY);
                const gsField<real_t> y(divdivAssembler.patches(), mpY);
            }

        }

        majVector[refCounter] = maj;
        mdVector[refCounter] = mD;
        meqVector[refCounter] = mEq;
        majhVector[refCounter] = majh;
        majIIVector[refCounter] = majII; // = majII / majIIGap;
        majIIGapVector[refCounter] = majIIGap;

    }
    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsRecontructV(int refCount,
                                                    gsSpaceTimeAssembler<real_t>& spaceTimeAssemblerV,
                                                    gsBoundaryConditions<> &bcInfo,
                                                    gsMatrix<real_t> &vVector,
                                                    gsMultiPatch<real_t> &mpV,
                                                    gsField<>& v,
                                                    gsVector<index_t> &vDOFs,
                                                    gsVector<real_t> &stopcritVector,
                                                    gsVector<double> &timeAsmbV,
                                                    gsMatrix<double> &timeSolvV) {
        gsCPUStopwatch clock;
        spaceTimeAssemblerV.refresh();

        //! [Assemble System of Discretized Poisson Equation]
        // ---------------------------------------------------------------------------------------------------------- //
        clock.restart();
        spaceTimeAssemblerV.assemble();
        timeAsmbV[refCount] = clock.stop();
        //gsInfo << "time for assembling V-system: " << timeAsmbV[refCount] << " sec.\n";
        //! [Assemble System of Discretized Poisson Equation]


        //! [Dof mapper]
        // ---------------------------------------------------------------------------------------------------------- //
        gsDofMapper spaceTimeMapper; // Gets the indices mapped from Basis --> Matrix
        spaceTimeAssemblerV.multiBasis(0).getMapper((dirichlet::strategy) spaceTimeAssemblerV.options().getInt("DirichletStrategy"),
                                                    (iFace::strategy) spaceTimeAssemblerV.options().getInt("InterfaceStrategy"),
                                                    bcInfo, spaceTimeMapper, 0);
        vDOFs[refCount] = spaceTimeMapper.size();
        //! [Dof mapper]

        //! [Solve System]
        // ---------------------------------------------------------------------------------------------------------- //
        this->gsSolveKhfhSystem(spaceTimeAssemblerV, vVector, timeSolvV, refCount, stopcritVector);
        //! [SolvevSystem]

        //! [Recover the Approximation Field]
        // Construct the solution as gsMultiPatch and gsField
        spaceTimeAssemblerV.constructSolution(vVector, mpV);    // gsMultiPatch
        v = gsField<>(spaceTimeAssemblerV.patches(), mpV);      // gsField
        //! [Recover the Approximation Field]

    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsRecontructW(int refCounter,
                                                   gsSpaceTimeAssembler<real_t> &spaceTimeAssemblerW,
                                                   gsBoundaryConditions<> &bcInfo,
                                                   gsMatrix<real_t> &wVector,
                                                   gsMultiPatch<real_t> &mpW,
                                                   gsField<> &w,
                                                   gsVector<index_t> &wDOFs,
                                                   gsVector<real_t> &stopcritVector,
                                                   gsVector<double> &timeAsmblW,
                                                   gsMatrix<double> &timeSolvW) {
        gsCPUStopwatch clock;

        spaceTimeAssemblerW.refresh();

        //! [Assemble System of Discretized Poisson Equation]
        // ---------------------------------------------------------------------------------------------------------- //
        clock.restart();
        spaceTimeAssemblerW.assemble();
        timeAsmblW[refCounter] = clock.stop();
        //! [Assemble System of Discretized Poisson Equation]

        //! [Dof mapper]
        // ---------------------------------------------------------------------------------------------------------- //
        gsDofMapper spaceTimeMapper; // Gets the indices mapped from Basis --> Matrix
        spaceTimeAssemblerW.multiBasis(0).getMapper(
                (dirichlet::strategy) spaceTimeAssemblerW.options().getInt("DirichletStrategy"),
                (iFace::strategy) spaceTimeAssemblerW.options().getInt("InterfaceStrategy"),
                bcInfo, spaceTimeMapper, 0);
        wDOFs[refCounter] = spaceTimeMapper.size();
        //! [Dof mapper]

        //! [Solve System]
        // ---------------------------------------------------------------------------------------------------------- //
        this->gsSolveKhwhfhSystem(spaceTimeAssemblerW, wVector, timeSolvW, refCounter, stopcritVector);
        //! [SolvevSystem]

        //! [Recover the Approximation Field]
        // Construct the solution as gsMultiPatch and gsField
        spaceTimeAssemblerW.constructSolution(wVector, mpW);
        w = gsField<>(spaceTimeAssemblerW.patches(), mpW);
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsGenerateDivMatrixRhs(int refCounter,
                                                            const gsDivPde<real_t> &divPde,
                                                            const gsMultiBasis<real_t> &basisY, int yDegree,
                                                            gsVector<index_t> &yDOFs,
                                                            gsAssembler<> &divdivAssembler,
                                                            gsVector<double> &timeAsmblDivDivY) {
        int numMappers = 1;
        gsCPUStopwatch clock;

        gsDofMapper divdivMapper; // Gets the indices mapped from Basis --> Matrix
        divdivAssembler.options().setInt("InterfaceStrategy", iFace::conforming);
        basisY.getMapper((iFace::strategy) divdivAssembler.options().getInt("InterfaceStrategy"),
                         divdivMapper, 0);
        divdivMapper.finalize();
        yDOFs[refCounter] = divdivMapper.size();
        if(this->withMajorant) {

            std::vector<gsDofMapper> divdivMappers(numMappers, divdivMapper);
            // Generate the sparse matrix for div(y) * div(w) = f * div(w) system
            clock.restart();
            gsSparseSystem<> divdivSys(divdivMappers, this->dim, this->dim);

            //divdivSys.reserve(basisY.at(0), divdivAssembler.options(), divPde.numRhs(), this->dim); // new function added into gsSparseSystem.h for the vector-valued systems
            divdivSys.reserve(basisY.at(0), divdivAssembler.options(), divPde.numRhs());          // old alternative that reserve only the space based on the scalar variable
            // TODO: figure out the optimal way to reserve
            // TODO: numToReserve is sufficient for THB splines ?
            //real_t numToReserve = math::ipow(this->dim * (this->dim * yDegree + 1), this->dim);
            //divdivSys.reserve(numToReserve, divPde.numRhs());                                     // reserve based on the above calculated number
            divdivAssembler.setSparseSystem(divdivSys);             // set the spare matrix for the system
            divdivAssembler.push<gsVisitorDivDivSpaceTime<real_t> >();       // push the assembling procedure
            divdivAssembler.finalize();
            timeAsmblDivDivY[refCounter] = clock.stop();
        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsGenerateMMMatrixRhs(int refCounter,
                                                           const gsPoissonPde<> &dualPde,
                                                           const gsMultiBasis<real_t> &basisY, int yDegree,
                                                           gsVector<index_t> &yDOFs,
                                                           gsAssembler<> &dualAssembler,
                                                           gsVector<double> &timeAsmblVectMassY) {
        //gsBoundaryConditions<> freeBC;
        int numMappers = 1;
        gsCPUStopwatch clock;

        //if(this->withMajorant) {

            // Initialtize Dof poissonMapper for y = grad(v) equation with free BC
            gsDofMapper dualMapper;
            basisY.getMapper((iFace::strategy) dualAssembler.options().getInt("InterfaceStrategy"), dualMapper, 0);
            dualMapper.finalize();
            std::vector<gsDofMapper> dualMappers(numMappers, dualMapper);

            clock.restart();
            gsSparseSystem<> vectMassSys(dualMappers, this->dim, this->dim);
            //real_t numToReserve = math::ipow(this->dim * (this->dim * yDegree + 1), this->dim);
            //gsInfo << "numToReserve for vectMass = " << numToReserve << "\n";
            vectMassSys.reserve(basisY.at(0), dualAssembler.options(), dualPde.numRhs());
            //vectMassSys.reserve(basisY.at(0), dualAssembler.options(), dualPde.numRhs(), this->dim);
            //vectMassSys.reserve(numToReserve, dualPde.numRhs());
            dualAssembler.setSparseSystem(vectMassSys);
            dualAssembler.push<gsVisitorDualPoissonSpaceTime<real_t> >();
            dualAssembler.finalize();
            timeAsmblVectMassY[refCounter] = clock.stop();
        //}

    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsLogTestResults(int vDegree, int yDegree, int wDegree,
                                             int m, int l,
                                             int yBasisRefDelay, int wBasisRefDelay,
                                             real_t markingParamTheta, int numInitUniformRef, int numTotalAdaptRef,
                                             gsVector<index_t> &vDOFs, gsVector<index_t> &yDOFs, gsVector<index_t> &wDOFs,
                                             gsVector<double> &timeAsmbV, gsVector<double> &timeAsmbDivDivY, gsVector<double> &timeAsmbMMY, gsVector<double> &timeAsmbY, gsVector<double> &timeAsmbW,
                                             gsMatrix<double> &timeSolvV, gsMatrix<double> &timeSolvY, gsMatrix<double> &timeSolvW,
                                             gsVector<double> &timeAsmbH1Error,
                                             gsVector<double> &timeAsmbSpaceTimeSolOperError, gsVector<double> &timeAsmbSpaceTimeDeltaxError, gsVector<double> &timeAsmbSpaceTimeDtError,
                                             gsVector<double> &timeAsmbMajorant, gsVector<double> &timeAsmbDeltaHMajorant, gsVector<double> &timeAsmbMajII,
                                             gsVector<double> &timeAsmbEtaIndicator, gsVector<double> &timeAsmbSpaceTimeErrorIdentity,
                                             gsVector<real_t> &eL2Vector, gsVector<real_t> &eH1Vector,
                                             gsVector<real_t> &eSpaceTimeVector,
                                             gsVector<real_t> &eSpaceTimeSpaceGradVector, gsVector<real_t> &eFullSpaceTimeSpaceGradVector,
                                             gsVector<real_t> &eSpaceTimeSolOperVector, gsVector<real_t> &eSpaceTimeDeltaxVector, gsVector<real_t> &eSpaceTimeDtVector,
                                             gsVector<real_t> &relErrorVector, gsVector<real_t> &relError0Vector,
                                             gsVector<real_t> &majVector, gsVector<real_t> &mdVector, gsVector<real_t> &meqVector, gsVector<real_t> &majhVector, gsVector<real_t> &majIIVector, gsVector<real_t> &majIIGapVector,
                                             gsVector<real_t> &minVector, gsVector<real_t> &etaVector, gsVector<real_t> &eIdentVector) {

        gsInfo
                << "\n%-------------------------------------------------------------------------------------------------%\n";

        if (this->isAdaptive) gsInfo << "Adaptive refinement results: \n";
        else
            gsInfo << "Uniform refinement results: \n";
        gsInfo
                << "%-------------------------------------------------------------------------------------------------%\n";

        gsInfo << "v in S^{" << vDegree << ", " << vDegree << "}\n";
        gsInfo << "y in S^{" << yDegree << ", " << yDegree << "} + S^{" << yDegree << ", " << yDegree << "}\n";
        gsInfo << "w in S^{" << wDegree << "}\n\n";


        gsInfo << "Inital refiniments:     REF_0 = " << numInitUniformRef << "\n";
        gsInfo << "Total refiniments:      REF   = " << numTotalAdaptRef << "\n";
        gsInfo << "Difference in the degrees of v and y: m = " << m << "\n";
        gsInfo << "Difference in the degrees of v and w: l = " << l << "\n";


        if (this->isAdaptive) {
            gsInfo << "Bulk parameter:         theta = " << markingParamTheta * 100 << " % \n";
            gsInfo << "Delay in refining of y: M     = " << yBasisRefDelay << "\n";
            gsInfo << "Delay in refining of w: L     = " << wBasisRefDelay << "\n";
        } else {
            gsInfo << "Coarsening parameter: M = " << yBasisRefDelay << "\n";
            gsInfo << "Coarsening parameter: L = " << wBasisRefDelay << "\n";

        }
        index_t lastRefIndx = numTotalAdaptRef - 1;
        /*
        * Time comparison
        */
        gsInfo << "\n\n REF &    DOFs(v) &    DOFs(y) &    DOFs(w) & "
               << " t_{as}(v) &  t_{as}(y) &  t_{as}(w) &"
               << " t_{sol,  dir}(v) & t_{sol,  dir}(y) & t_{sol,  dir}(w) \\\\\n";
        for (int refCounter = 0; refCounter < numTotalAdaptRef; ++refCounter) {
            std::printf("%4u & %10u & %10u & %10u & %10.2e & %10.2e & %10.2e & %16.2e & %16.2e & %16.2e \\\\\n",
                        refCounter + 1,
                        vDOFs[refCounter],
                        yDOFs[refCounter] * this->dim,
                        wDOFs[refCounter],
                        timeAsmbV[refCounter],
                        timeAsmbY[refCounter], //timeAsmbDivDivY[refCounter] + timeAsmbDivDivY[refCounter],
                        timeAsmbW[refCounter],
                        timeSolvV(refCounter, 0),
                        timeSolvY(refCounter, 0),
                        timeSolvW(refCounter, 0));
        }
        gsInfo << "---------------------------------------------------------------------------------------\\\\\n";
        std::printf(" \t & \t & \t & \t & %10.2f & %10.2f & %10.2f & %16.2f & %16.2f & %16.2f \\\\\n",
                    timeAsmbV[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbY[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbW[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeSolvV(lastRefIndx, 0) / timeSolvW(lastRefIndx, 0),
                    timeSolvY(lastRefIndx, 0) / timeSolvW(lastRefIndx, 0),
                    timeSolvW(lastRefIndx, 0) /  timeSolvW(lastRefIndx, 0));


        gsInfo << "\n\n REF &    DOFs(v) &    DOFs(y) &    DOFs(w) & "
               << " t_{as}(v) &  t_{as}(y) &  t_{as}(w) &"
               << " t_{sol,  dir}(v) & t_{sol,  dir}(y) & t_{sol,  dir}(w) & ratio \\\\\n";
        for (int refCounter = 0; refCounter < numTotalAdaptRef; ++refCounter) {
            std::printf("%4u & %10u & %10u & %10u & %10.2e & %10.2e & %10.2e & %16.2e & %16.2e & %16.2e & %16.2f \\\\\n",
                        refCounter + 1,
                        vDOFs[refCounter],
                        yDOFs[refCounter] * this->dim,
                        wDOFs[refCounter],
                        timeAsmbV[refCounter],
                        timeAsmbY[refCounter], //timeAsmbDivDivY[refCounter] + timeAsmbDivDivY[refCounter],
                        timeAsmbW[refCounter],
                        timeSolvV(refCounter, 0),
                        timeSolvY(refCounter, 0),
                        timeSolvW(refCounter, 0),
                        (timeAsmbV[refCounter] + timeSolvV(refCounter, 0)) /
                                (timeAsmbY[refCounter] + timeAsmbW[refCounter] + timeSolvY(refCounter, 0) + timeSolvW(refCounter, 0)));
        }
        gsInfo << "---------------------------------------------------------------------------------------\\\\\n";
        std::printf(" \t & \t & \t & \t & %10.2f & %10.2f & %10.2f & %16.2f & %16.2f & %16.2f & %16.2f \\\\\n",
                    timeAsmbV[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbY[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbW[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeSolvV(lastRefIndx, 0) / timeSolvW(lastRefIndx, 0),
                    timeSolvY(lastRefIndx, 0) / timeSolvW(lastRefIndx, 0),
                    timeSolvW(lastRefIndx, 0) /  timeSolvW(lastRefIndx, 0),
                    (timeAsmbV[lastRefIndx] + timeSolvV(lastRefIndx, 0)) /
                    (timeAsmbY[lastRefIndx] + timeAsmbW[lastRefIndx] + timeSolvY(lastRefIndx, 0) + timeSolvW(lastRefIndx, 0)));

        gsInfo << "\n\n REF &    DOFs(v) &    DOFs(y) &    DOFs(w) & "
               << " t_{as}(v) &  t_{as}(y) &  t_{as}(w) &"
               << " t_{sol,  it}(v) & t_{sol,  it}(y) & t_{sol,  it}(w) \\\\\n";
        for (int refCounter = 0; refCounter < numTotalAdaptRef; ++refCounter) {
            std::printf("%4u & %10u & %10u & %10u & %10.2e & %10.2e & %10.2e & %16.2e & %16.2e & %16.2e \\\\\n",
                        refCounter + 1,
                        vDOFs[refCounter],
                        yDOFs[refCounter] * this->dim,
                        wDOFs[refCounter],
                        timeAsmbV[refCounter],
                        timeAsmbY[refCounter], //timeAsmbDivDivY[refCounter] + timeAsmbDivDivY[refCounter],
                        timeAsmbW[refCounter],
                        timeSolvV(refCounter, 1),
                        timeSolvY(refCounter, 1),
                        timeSolvW(refCounter, 1));
        }
        gsInfo << "---------------------------------------------------------------------------------------\\\\\n";
        std::printf(" \t & \t & \t & \t & %10.2f & %10.2f & %10.2f & %16.2f & %16.2f & %16.2f \\\\\n",
                    timeAsmbV[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbY[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbW[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeSolvV(lastRefIndx, 1) / timeSolvW(lastRefIndx, 1),
                    timeSolvY(lastRefIndx, 1) / timeSolvW(lastRefIndx, 1),
                    timeSolvW(lastRefIndx, 1) /  timeSolvW(lastRefIndx, 1));

        gsInfo << "\n\n REF &    DOFs(v) &    DOFs(y) &    DOFs(w) & "
               << " t_{as}(v) &  t_{as}(y) &  t_{as}(w) &"
               << " t_{sol,  it}(v) & t_{sol,  it}(y) & t_{sol,  it}(w) & ratio \\\\\n";
        for (int refCounter = 0; refCounter < numTotalAdaptRef; ++refCounter) {
            std::printf("%4u & %10u & %10u & %10u & %10.2e & %10.2e & %10.2e & %16.2e & %16.2e & %16.2e & %16.2f \\\\\n",
                        refCounter + 1,
                        vDOFs[refCounter],
                        yDOFs[refCounter] * this->dim,
                        wDOFs[refCounter],
                        timeAsmbV[refCounter],
                        timeAsmbY[refCounter], //timeAsmbDivDivY[refCounter] + timeAsmbDivDivY[refCounter],
                        timeAsmbW[refCounter],
                        timeSolvV(refCounter, 1),
                        timeSolvY(refCounter, 1),
                        timeSolvW(refCounter, 1),
                        (timeAsmbV[refCounter] + timeSolvV(refCounter, 1)) /
                        (timeAsmbY[refCounter] + timeAsmbW[refCounter] + timeSolvY(refCounter, 1) + timeSolvW(refCounter, 1)));
        }
        gsInfo << "---------------------------------------------------------------------------------------\\\\\n";
        std::printf(" \t & \t & \t & \t & %10.2f & %10.2f & %10.2f & %16.2f & %16.2f & %16.2f & %16.2f \\\\\n",
                    timeAsmbV[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbY[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeAsmbW[lastRefIndx] / timeAsmbW[lastRefIndx],
                    timeSolvV(lastRefIndx, 1) / timeSolvW(lastRefIndx, 1),
                    timeSolvY(lastRefIndx, 1) / timeSolvW(lastRefIndx, 1),
                    timeSolvW(lastRefIndx, 1) /  timeSolvW(lastRefIndx, 1),
                    (timeAsmbV[lastRefIndx] + timeSolvV(lastRefIndx, 1)) /
                    (timeAsmbY[lastRefIndx] + timeAsmbW[lastRefIndx] + timeSolvY(lastRefIndx, 1) + timeSolvW(lastRefIndx, 1)));

        gsInfo << "\n\n REF &    DOFs(v) &    DOFs(y) &    DOFs(w) & "
               << " t_{as}(|||e|||_L) &  t_{as}(|| Deltax u||) &  t_{as}(|| Dt u||) & t_{as}(ID) &"
               << " t_{as}(||e||_H1) & t_{as}(maj I) \\\\\n";
        for (int refCounter = 0; refCounter < numTotalAdaptRef; ++refCounter) {
            std::printf("%4u & %10u & %10u & %10u & %10.2e & %10.2e & %10.2e & %16.2e & %16.2e & %16.2e \\\\\n",
                        refCounter + 1,
                        vDOFs[refCounter],
                        yDOFs[refCounter] * this->dim,
                        wDOFs[refCounter],
                        timeAsmbSpaceTimeSolOperError[refCounter],
                        timeAsmbSpaceTimeDeltaxError[refCounter], //timeAsmbDivDivY[refCounter] + timeAsmbDivDivY[refCounter],
                        timeAsmbSpaceTimeDtError[refCounter],
                        timeAsmbSpaceTimeErrorIdentity[refCounter],
                        timeAsmbH1Error[refCounter],
                        timeAsmbMajorant[refCounter]);
        }

        gsInfo << "\n REF & || grad_x e || &   \t maj & I_{eff}(maj) & \t   majII & I_{eff}(majII) &    ||| e |||_h & \t  majh & I_{eff}(majh) &   "
               << "   ||| e |||_L &     idt &  I_eff(idt) &  order(||| e |||_h) &  order(||| e |||_L) \\\\\n";
        for (int refCounter = 0; refCounter < lastRefIndx; ++refCounter) {
            std::printf("%4u & %14.4e & %12.4e & %12.2f & %12.4e & %12.2f & %14.4e & %12.4e & %12.2f & %14.4e & %12.4e & %12.2f & %8.2f & %8.2f \\\\\n",
                        refCounter + 2,
                        cast<real_t, double>(eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1] / eSpaceTimeVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(math::log(eSpaceTimeVector[refCounter + 1] / eSpaceTimeVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))),
                        cast<real_t, double>(math::log(eSpaceTimeSolOperVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter]) /
                                         std::log(std::pow((double) vDOFs[refCounter + 1],
                                                           ((double) (-1) / (double) d)) /
                                                  std::pow((double) vDOFs[refCounter],
                                                           ((double) (-1) / (double) d)))));
        }


        gsInfo << "\n REF & || grad_x e || &   \t maj & I_{eff}(maj) & majII/gap & I(majII / gap) &   ||| e |||_h & \t majh & I_{eff}(majh) &   "
               << " \t ||| e |||_L &      idt & I_eff(idt) &  order(||| e |||_h) &  order(||| e |||_L) \\\\\n";
        for (int refCounter = 0; refCounter < lastRefIndx; ++refCounter) {
            std::printf("%4u & %14.4e & %12.4e & %12.2f & %12.4e & %12.2f & %14.4e & %12.4e & %12.2f & %14.4e & %12.4e & %12.2f & %8.2f & %8.2f \\\\\n",
                        refCounter + 2,
                        cast<real_t, double>(eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1] / majIIGapVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1] / majIIGapVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1] / eSpaceTimeVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(math::log(eSpaceTimeVector[refCounter + 1] / eSpaceTimeVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))),
                        cast<real_t, double>(math::log(eSpaceTimeSolOperVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))));
        }

        gsInfo << "\n REF & ||| e |||_L &  || Delta_x e ||_Q & || Dt e ||_Q  & \t    idt & I_eff(idt) &  order(||| e |||_L) \\\\\n";
        for (int refCounter = 0; refCounter < lastRefIndx; ++refCounter) {
            std::printf("%4u & %14.4e & %14.4e & %14.4e & %12.4e & %12.2f & %8.2f \\\\\n",
                        refCounter + 2,
                        cast<real_t, double>(eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeDeltaxVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeDtVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(math::log(eSpaceTimeSolOperVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))));
        }


        gsInfo << "\n REF & || grad_x e || & I_{eff}(maj) & I(majII / gap) &   ||| e |||_h & I_{eff}(majh) &   "
               << " \t ||| e |||_L & I_eff(idt) &  order(||| e |||_h) &  order(||| e |||_L) \\\\\n";
        for (int refCounter = 0; refCounter < lastRefIndx; ++refCounter) {
            std::printf("%4u & %14.4e & %12.2f & %12.2f & %14.4e & %12.2f & %14.4e & %12.2f & %8.2f & %8.2f \\\\\n",
                        refCounter + 2,
                        cast<real_t, double>(eSpaceTimeSpaceGradVector[refCounter + 1]),
                        //cast<real_t, double>(majVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        //cast<real_t, double>(majIIVector[refCounter + 1] / majIIGapVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1] / majIIGapVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeVector[refCounter + 1]),
                        //cast<real_t, double>(majhVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1] / eSpaceTimeVector[refCounter + 1]),
                        cast<real_t, double>(eSpaceTimeSolOperVector[refCounter + 1]),
                        //cast<real_t, double>(eIdentVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(math::log(eSpaceTimeVector[refCounter + 1] / eSpaceTimeVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))),
                        cast<real_t, double>(math::log(eSpaceTimeSolOperVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))));
        }

        gsInfo << "\n REF & \t maj & I_{eff}(maj) & majII/gap & I(majII / gap) & \t majh & I_{eff}(majh) &   "
               << " \t idt & I_eff(idt) &  order(||| e |||_h) &  order(||| e |||_L) \\\\\n";
        for (int refCounter = 0; refCounter < lastRefIndx; ++refCounter) {
            std::printf("%4u & %12.4e & %12.2f & %12.4e & %12.2f & %12.4e & %12.2f & %12.4e & %12.2f & %8.2f & %8.2f \\\\\n",
                        refCounter + 2,
                        //cast<real_t, double>(eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1]),
                        cast<real_t, double>(majVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1] / majIIGapVector[refCounter + 1]),
                        cast<real_t, double>(majIIVector[refCounter + 1] / majIIGapVector[refCounter + 1] / eSpaceTimeSpaceGradVector[refCounter + 1]),
                        //cast<real_t, double>(eSpaceTimeVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1]),
                        cast<real_t, double>(majhVector[refCounter + 1] / eSpaceTimeVector[refCounter + 1]),
                        //cast<real_t, double>(eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1]),
                        cast<real_t, double>(eIdentVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter + 1]),
                        cast<real_t, double>(math::log(eSpaceTimeVector[refCounter + 1] / eSpaceTimeVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))),
                        cast<real_t, double>(math::log(eSpaceTimeSolOperVector[refCounter + 1] / eSpaceTimeSolOperVector[refCounter]) /
                                             std::log(std::pow((double) vDOFs[refCounter + 1],
                                                               ((double) (-1) / (double) d)) /
                                                      std::pow((double) vDOFs[refCounter],
                                                               ((double) (-1) / (double) d)))));
        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsLogRefinementIterationErrorReport(const int refCount, const real_t theta, const gsVector<real_t> & hmaxVector, const gsVector<real_t> & hminVector,
                                                                       const gsVector<real_t> & eL2Vector, const gsVector<real_t> & eH1Vector,
                                                                       const gsVector<real_t> & eSpaceTimeSpaceGradVector, const gsVector<real_t> & eFullSpaceTimeSpaceGradVector,
                                                                       const gsVector<real_t> & ehQVector, const gsVector<real_t> & eFullSpaceTimeVector,
                                                                       const gsVector<real_t> & eSpaceTimeSolOperVector, const gsVector<real_t> & eFullSpaceTimeSolOperVector,
                                                                       const gsVector<real_t> & eIdentVector, const gsVector<real_t> & eFullIdentVector,
                                                                       const gsVector<real_t> & gradxe0Vector, const gsVector<real_t> & gradxeTVector,
                                                                       const gsVector<real_t> & e0Vector, const gsVector<real_t> & eTVector,
                                                                         const gsVector<real_t> & ehSigmaTVector) const {

        gsInfo << "\n|| e ||_Q        = " << eL2Vector[refCount] << "\n";
        gsInfo << "|| grad e ||_Q   = " << eH1Vector[refCount] << "\n\n";

        gsInfo << "|| grad_x e ||_Q                       = " << eSpaceTimeSpaceGradVector[refCount] << "\n";
        gsInfo << "(|| grad_x e ||^2_Q + || e ||^2_T)^1/2 = " << eFullSpaceTimeSpaceGradVector[refCount] << "\n\n";

        gsInfo << "|| e ||_{s, h, Q}   = " << ehQVector[refCount] << "\n";
        gsInfo << "|| e ||_{s, h, T}   = " << ehSigmaTVector[refCount] << "\n";
        gsInfo << "(|| e ||^2_T + theta * h * || grad_x e ||^2_T)^1/2   = " << math::sqrt(math::pow(eTVector[refCount], 2) + theta * hminVector[refCount] * math::pow(gradxeTVector[refCount], 2)) << "\n";
        gsInfo << "(|| e ||^2_{s, h, Q} + || e ||^2_T + theta * h * || grad_x e ||^2_T)^1/2 = " << eFullSpaceTimeVector[refCount] << "\n";
        gsInfo << "(|| e ||^2_{s, h, Q} + || e ||^2_{s, h, T})^1/2 = " << math::sqrt(math::pow(ehQVector[refCount], 2) + math::pow(ehSigmaTVector[refCount], 2)) << "\n\n";

        gsInfo << "(|| delta_x e ||^2_Q + || e_t ||^2_Q)^1/2 = " << eSpaceTimeSolOperVector[refCount] << "\n";
        gsInfo << "|| delta_x v + f - v_t ||                 = " << eIdentVector[refCount] << "\n\n";

        gsInfo << "(|| delta_x e ||^2_Q + || e_t ||^2_Q + || grad_x e ||^2_T)^1/2 = " << eFullSpaceTimeSolOperVector[refCount] << "\n";
        gsInfo << "(|| delta_x v + f - v_t ||^2 + || grad_x e ||^2_0)^1/2         = " << eFullIdentVector[refCount] << "\n\n";

        gsInfo << "|| e ||_0 = " << e0Vector[refCount] << "\n";
        gsInfo << "|| e ||_T = " << eTVector[refCount] << "\n";
        gsInfo << "|| grad_x e ||_0 = " << gradxe0Vector[refCount] << "\n";
        gsInfo << "|| grad_x e ||_T = " << gradxeTVector[refCount] << "\n\n";

        gsInfo << "theta = " << theta << "\n";
        gsInfo << "hmax = " << hmaxVector[refCount] << "\n";
        gsInfo << "hmin = " << hminVector[refCount] << "\n\n";

        //gsInfo << "|| e_w ||_Q        = " << eWL2Vector[refCount] << "\n";
        //gsInfo << "|| grad e_w ||_Q   = " << eWH1Vector[refCount] << "\n\n";
    }


    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsLogRefinementIterationInfo(int refCounter,
                                                         gsVector<index_t> &vDOFs, gsVector<index_t> &yDOFs, gsVector<index_t> &wDOFs,
                                                         gsVector<real_t> &eVector, gsVector<real_t> &eL2Vector,
                                                         gsVector<real_t> &eSpaceTimeVector,
                                                         gsVector<real_t> &eSpaceTimeSeminormVector, gsVector<real_t> &eFullSpaceTimeSeminormVector,
                                                         gsVector<real_t> &eSpaceTimeSolOperVector,
                                                         gsVector<real_t> &relErrorVector, gsVector<real_t> &relError0Vector, gsVector<real_t> &thetaVector, gsVector<real_t> &stopcritVector,
                                                         gsVector<real_t> &majVector, gsVector<real_t> &majhVector, gsVector<real_t> &majIIVector, gsVector<real_t> &majIIGapVector, gsVector<real_t> &minVector,
                                                         gsVector<real_t> &etaVector,  gsVector<real_t> &eIdentVector) {

        gsInfo
                << "\n%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << " Resulting values of the majorant and the residual error indicator: \n";
        gsInfo
                << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "DOFs(v)   = " << vDOFs[refCounter] << "\n"
               << "DOFs(y)   = " << yDOFs[refCounter] * this->dim << "\n"
               << "DOFs(w)   = " << wDOFs[refCounter] << "\n\n";
        gsInfo //<< "min   = " << minVector[refCounter] << "\t"
               //<< "iEff(min)   = " << minVector[refCounter] / eVector[refCounter] << "\n"
               << "||  e ||_L2  = " << eL2Vector[refCounter] << "\n\n"
               << "||| e |||_Q  = " << eFullSpaceTimeSeminormVector[refCounter] << "\n"
                << "maj         = " << majVector[refCounter] << "\t"
                << "iEff(maj) = " << majVector[refCounter] / eFullSpaceTimeSeminormVector[refCounter] << "\n"
                << "majII       = " << majIIVector[refCounter] << "\t"
                << "iEff(majII) = " << majIIVector[refCounter] / eSpaceTimeSeminormVector[refCounter] << "\n\n"
                << "majII / CgapII = " << majIIVector[refCounter] / majIIGapVector[refCounter]  << "\t"
                << "iEff(majII / CgapII) = " << majIIVector[refCounter] / majIIGapVector[refCounter] / eSpaceTimeSeminormVector[refCounter] << "\n\n"
                //
                << "||| e |||_h = " << eSpaceTimeVector[refCounter] << "\n"
                << "majh        = " << majhVector[refCounter] << "\t"
                << "iEff(majh) = " << majhVector[refCounter] / eSpaceTimeVector[refCounter] << "\n\n"
                << "|| e ||_L   = " << eSpaceTimeSolOperVector[refCounter] << "\n"
                << "idt         = " << eIdentVector[refCounter] << "\t"
                << "iEff(idt) = " << eIdentVector[refCounter] / eSpaceTimeSolOperVector[refCounter] << "\n\n";
        //gsInfo << "eta_K = " << etaVector[refCounter] << "\t"
        //       << "iEff(eta_K) = " << etaVector[refCounter] / eVector[refCounter] << "\n\n";

        if (refCounter > 0) {
            gsInfo << "Theta_i     = [e]_i / [e]_{i-1}   = " << thetaVector[refCounter-1] << "\n";
            gsInfo << "eps0_i      = || u_i - u_0 ||     = " << relErrorVector[refCounter-1] << "\n";
            gsInfo << "delta_i     = || u_i - u_{i-1} || = " << relError0Vector[refCounter-1] << "\n";
            gsInfo << "stop_crit_i = (1 - Theta_i) * TOL * rho * eps0_i = " << stopcritVector[refCounter-1] << "\n";

            gsInfo << "|| u - u_{i-1} || <= || u_i - u_{i-1} || / (1 - Theta): "
                   << relErrorVector[refCounter-1] << " <= "
                   << relErrorVector[refCounter-1] / (1 - eVector[refCounter] / eVector[refCounter - 1])<< "\n\n";
        };
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsInitializeBasis(int degree,
                                              gsMultiBasis<> &multBasis,
                                              int numInitUniformRef) {
        if (this->isAdaptive) {
            // Copy basis from the multi-patch geometry (one per patch) for v and flux
            gsTensorBSpline<d, real_t> *geo = dynamic_cast< gsTensorBSpline<d, real_t> * >( &(this->patches.patch(0)));
            gsTensorBSplineBasis<d, real_t> tensorBSplineBasis = geo->basis();

            if (patches.basis(0).degree(0) < patches.basis(0).degree(1) && d == 2) {
                tensorBSplineBasis.degreeElevate(1, 0);

            } else if (patches.basis(0).degree(0) > patches.basis(0).degree(1) && d == 2) {
                tensorBSplineBasis.degreeElevate(1, 1);
            }
            // Set the degree of the basises for the approximation
            tensorBSplineBasis.setDegree(degree);

            gsTHBSplineBasis<d, real_t> thbBasis(tensorBSplineBasis);
            multBasis = gsMultiBasis<>(thbBasis);

        } else {
            // Copy basis from the multi-patch geometry (one per patch) for v and flux
            multBasis = gsMultiBasis<>(patches);

            if (d == 2) {
                if (patches.basis(0).degree(0) < patches.basis(0).degree(1)) {
                    multBasis.degreeElevate(1, 0);
                } else if (patches.basis(0).degree(0) > patches.basis(0).degree(1)) {
                    multBasis.degreeElevate(1, 1);
                }
            }

            if (d == 3){
                if (patches.basis(0).degree(0) < patches.basis(0).degree(1)) {
                    multBasis.degreeElevate(1, 0);
                } else if (patches.basis(0).degree(0) > patches.basis(0).degree(1)) {
                    multBasis.degreeElevate(1, 1);
                }
                if (patches.basis(0).degree(2) < patches.basis(0).degree(1)){
                    multBasis.degreeElevate(1, 2);
                } else if (patches.basis(0).degree(2) > patches.basis(0).degree(1)){
                    multBasis.degreeElevate(1, 1);
                }
            }


            multBasis.setDegree(degree);
        }

        //! [Initial Refinement]
        for (int i = 0; i < numInitUniformRef; ++i)    multBasis.uniformRefine();
        //! [Initial Refinement]
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsGetInitialBasis(int vDegree, int yDegree, int wDegree,
                                              gsMultiBasis<> &multBasisV,
                                              gsMultiBasis<> &multBasisY,
                                              gsMultiBasis<> &multBasisW,
                                              int numInitUniformRefV,
                                              int numInitUniformRefY,
                                              int numInitUniformRefW) {

        gsInitializeBasis(vDegree, multBasisV, numInitUniformRefV);
        gsInitializeBasis(yDegree, multBasisY, numInitUniformRefY);
        gsInitializeBasis(wDegree, multBasisW, numInitUniformRefW);
        this->gsLogGeometryBasisData(this->patches.patch(0),
                                     vDegree, yDegree, wDegree, 1);

        gsInfo << "VBasis:      \n" << multBasisV.basis(0) << "\n";
        gsInfo << "YBasis:      \n" << multBasisY.basis(0) << "\n";
        gsInfo << "WBasis:      \n" << multBasisW.basis(0) << "\n";
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsLogGeometryBasisData(gsGeometry<> &pGeom,
                                                   int vDegree, int yDegree, int wDegree, int numPatches) {
        gsInfo
                << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "Geometry and basis info: \n";
        gsInfo
                << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "Geom. coeff: \n" << pGeom.coefs() << "\n"; // [N x dim] matrix of the [N control points] x [dim]
        gsInfo << "Geom. basis: \n" << pGeom.basis() << "\n";
        gsInfo
                << "%-------------------------------------------------------------------------------------------------%\n";
        gsInfo << "$v \\in S^{" << vDegree << ", " << vDegree << "}$\n";
        gsInfo << "$y \\in S^{" << yDegree << ", " << yDegree << "} + S^{" << yDegree << ", " << yDegree << "}$\n";
        gsInfo << "$w \\in S^{" << wDegree << ", " << wDegree << "}$\n";
        gsInfo
                << "%-------------------------------------------------------------------------------------------------%\n";
    }

    /*
     * poissonAssembler, basisY, basisW, basisYVector, basisWVector,
     * mdDistr,
                                             adaptRefCrit,
                                             markingParamTheta,
                                             refCounter, yBasisRefDelay,
                                             wBasisRefDelay
     */
    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsExecuteRefinement(gsSpaceTimeAssembler<real_t> &spaceTimeAssemblerV,
                                                gsMultiBasis<real_t> &thbMultBasisY, gsSpaceTimeAssembler<real_t> &spaceTimeAssemblerW,
                                                std::vector<gsMultiBasis<> > &thbMultBasisYVector, std::vector<gsMultiBasis<> > &thbMultBasisWVector,
                                                std::vector<real_t> &mdDistr, MarkingStrategy adaptRefCrit,
                                                real_t markingParamTheta,
                                                unsigned int refCounter,
                                                unsigned int yBasisRefDelay, unsigned int wBasisRefDelay) {

        gsCPUStopwatch clock;

        if (this->isAdaptive) {
            // ------------------------------------------------------------------------------------------------------ //
            // Refining spaceTimeAssemblerV mesh (basis for V)
            // ------------------------------------------------------------------------------------------------------ //
            //! [Marking procedure]
            // Get element-wise energy error indicator generated by majorant contributions
            std::vector<bool> elMarked(mdDistr.size());
            clock.restart();
            // Mark elements for refinement, based on the computed local indicators, the refinement-criterion and -parameter.
            gsMarkElementsForRef(mdDistr, adaptRefCrit, markingParamTheta, elMarked);
            gsInfo << "time for marking : " << clock.stop() << " sec.\n";
            //! [Marking procedure]

            //! [Refining procedure]
            // Refine the marked elements with a 1-ring of cells around marked elements
            clock.restart();
            gsRefineMarkedElements(spaceTimeAssemblerV.multiBasis(), elMarked);
            spaceTimeAssemblerV.multiBasis().repairInterfaces(this->patches.interfaces());
            gsInfo << "time for refining the basis for v : " << clock.stop() << " sec.\n";

            // ------------------------------------------------------------------------------------------------------ //
            // Saving basises for Y and W
            // ------------------------------------------------------------------------------------------------------ //
            if (refCounter == 0) {
                // Refinement with optimization of the effort
                for (unsigned int j = 0; j < yBasisRefDelay; j++) {
                    // First save the currect basis
                    thbMultBasisYVector.push_back(thbMultBasisY);
                }
                // Refinement with optimization of the effort
                for (unsigned int j = 0; j < wBasisRefDelay; j++) {
                    gsMultiBasis<>& basis = spaceTimeAssemblerW.multiBasis();
                    thbMultBasisWVector.push_back(basis);
                }

            }
            // ------------------------------------------------------------------------------------------------------ //
            // Refining basis for Y
            // ------------------------------------------------------------------------------------------------------ //

            // Make the refinement and safe the new basis for y
            clock.restart();
            gsRefineMarkedElements(thbMultBasisY, elMarked);
            thbMultBasisY.repairInterfaces(this->patches.interfaces());
            gsInfo << "time for refining the basis for y : " << clock.stop() << " sec.\n";

            real_t sizeBasisYVector = thbMultBasisYVector.size();
            gsInfo << "sizeBasisYVector = " << sizeBasisYVector << "\n";
            real_t lastStoredBasisYDOFs = (sizeBasisYVector <= 1) ? 0 : thbMultBasisYVector.at(sizeBasisYVector - 1).piece(0).anchors().size();
            gsInfo << "lastStoredBasisYDOFs = " << lastStoredBasisYDOFs << "\n";
            real_t nextBasisYDOFs = thbMultBasisY.piece(0).anchors().size();
            gsInfo << "nextBasisYDOFs       = " << nextBasisYDOFs << "\n";


            // Strategies of updating basis for Y and basis for W can be different
            // below is just one of them
            if (lastStoredBasisYDOFs < nextBasisYDOFs) {
                //if ( refCounter >= yBasisRefDelay - 1 && sizeBasisYVector >=2 &&
                //     thbMultBasisYVector.at(sizeBasisYVector - 1).piece(0).anchors().size() == thbMultBasisYVector.at(sizeBasisYVector - 2).piece(0).anchors().size())
                //    thbMultBasisYVector.pop_back();
                thbMultBasisYVector.push_back(thbMultBasisY);
            }
                /*
            else { // if the new basis is different from the original size
                gsMultiBasis<>& basis = thbMultBasisYVector.at(sizeBasisYVector - 1);
                thbMultBasisYVector.push_back(basis);
            }
                 */
            if (refCounter < thbMultBasisYVector.size())
                thbMultBasisY = thbMultBasisYVector[refCounter];
            else
                thbMultBasisY = thbMultBasisYVector[thbMultBasisYVector.size() - 1];


            // Make the refinement and safe the new basis for w
            clock.restart();
            gsRefineMarkedElements(spaceTimeAssemblerW.multiBasis(), elMarked);
            spaceTimeAssemblerW.multiBasis().repairInterfaces(this->patches.interfaces());
            gsInfo << "time for refining the basis for w : " << clock.stop() << " sec.\n";


            real_t sizeBasisWVector = thbMultBasisWVector.size();
            gsInfo << "sizeBasisWVector = " << sizeBasisWVector << "\n";
            real_t lastStoredBasisWDOFs = (sizeBasisWVector <= 1) ? 0 : thbMultBasisWVector.at(sizeBasisWVector - 1).piece(0).anchors().size() / this->dim;
            gsInfo << "lastStoredBasisWDOFs = " << lastStoredBasisWDOFs << "\n";
            real_t nextBasisWDOFs = spaceTimeAssemblerW.multiBasis().piece(0).anchors().size() / this->dim;
            gsInfo << "nextBasisWDOFs       = " << nextBasisWDOFs << "\n";

            if (lastStoredBasisWDOFs < nextBasisWDOFs) {
                //if (refCounter >= wBasisRefDelay - 1 && sizeBasisWVector >= 2 &&
                //    thbMultBasisWVector.at(sizeBasisWVector - 1).piece(0).anchors().size() ==
                //    thbMultBasisWVector.at(sizeBasisWVector - 2).piece(0).anchors().size())
                //    thbMultBasisWVector.pop_back();
                gsMultiBasis<> &basis = spaceTimeAssemblerW.multiBasis();
                thbMultBasisWVector.push_back(basis);
            }
                /*
            else {
                gsMultiBasis<>& basis = thbMultBasisWVector.at(sizeBasisWVector - 1);
                thbMultBasisWVector.push_back(basis);
            }
                 */
            if (refCounter < thbMultBasisWVector.size())
                spaceTimeAssemblerW.basisUpdate(thbMultBasisWVector[refCounter]);
            else
                spaceTimeAssemblerW.basisUpdate(thbMultBasisWVector.at(thbMultBasisWVector.size() - 1));

        } else {
            spaceTimeAssemblerV.multiBasis().uniformRefine();  // Uniform refinement of the basis for V
            if (refCounter >= yBasisRefDelay - 1)
                thbMultBasisY.uniformRefine();              // uniform refinement of basis of Y
            if (refCounter >= wBasisRefDelay - 1)
                spaceTimeAssemblerW.multiBasis().uniformRefine();       // uniform refinement of basis of W

        }
        spaceTimeAssemblerV.refresh();
        spaceTimeAssemblerW.refresh();
    }


    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsSaveToFileRefinementIterationInfo(bool save,
                                                                const gsField<> &v,
                                                                gsMultiBasis<> &basis,
                                                                std::vector<real_t> edDistr, std::vector<real_t> mdDistr, std::vector<real_t> eSolDistr, std::vector<real_t> idDistr,
                                                                gsVector<real_t> & e0Vector, gsVector<real_t> & eTVector, gsVector<real_t> & gradxe0Vector, gsVector<real_t> & gradxeTVector, int refCount,
                                                                int refCounter, int refTotal,
                                                                const unsigned exampleNumber) {
        if (save) {
            // Saving Paraview files
            std::string refTag = util::to_string(exampleNumber) + "-refnum-"
                                 + util::to_string(refCounter);
            if (this->isAdaptive){
                gsMesh<> mesh;
                gsTensorBSpline<d, real_t> *geo = dynamic_cast< gsTensorBSpline<d, real_t> * >( &(this->patches.patch(0)));
                makeMesh(basis.basis(0), mesh);
                std::string meshParam("param-domain-mesh-" + refTag),
                        meshPhys("phys-domain-mesh-" + refTag),
                        vPhys("v-" + refTag);

                gsWriteParaview(mesh, this->resultFolder + "/" + meshParam, true);
                geo->evaluateMesh(mesh);
                gsWriteParaview(mesh, this->resultFolder + "/" + meshPhys, true);
                gsWriteParaview(v, this->resultFolder + "/" + vPhys, 5001, true);
            }
            // Saving txt files
            std::string resultsFile = "results-" + refTag + ".txt";

            std::ofstream file((this->resultFolder + "/" + resultsFile).c_str());
            if (!file.is_open())
                std::cout << "Problem opening the result file!" << resultsFile << std::endl;
            else {
                int size = edDistr.size();
                file << "e_distr = [";
                for (int elemCount = 0; elemCount < size; ++elemCount) file << math::pow(edDistr[elemCount], 2) << ", ";
                file << "]; \n";

                file << "m_distr = [";
                for (int elemCount = 0; elemCount < size; ++elemCount) file << math::pow(mdDistr[elemCount], 2) << ", ";
                file << "]; \n";

                file << "eSolDistr = [";
                for (int elemCount = 0; elemCount < size; ++elemCount) file << math::pow(gradxeTVector[refCount], 2) + math::pow(eSolDistr[elemCount], 2) << ", ";
                file << "]; \n";

                file << "idDistr = [";
                for (int elemCount = 0; elemCount < size; ++elemCount) file << math::pow(gradxe0Vector[refCount], 2) + math::pow(idDistr[elemCount], 2) << ", ";
                file << "]; \n";
                file.close();
            }
        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsSaveToFileTestResults(bool save,
                                                             gsVector<index_t> &vDOFs, gsVector<index_t> &yDOFs, gsVector<index_t> &wDOFs,
                                                             gsVector<real_t> &eL2Vector, gsVector<real_t> &eH1Vector,
                                                             gsVector<real_t> &eSpaceTimeSpaceGradVector, gsVector<real_t> &eSpaceTimeVector, gsVector<real_t> &eSpaceTimeSolOperVector,
                                                             gsVector<real_t> &majVector, gsVector<real_t> &majhVector, gsVector<real_t> &majIIVector, gsVector<real_t> &majIIGapVector,
                                                             gsVector<real_t> &eIdentVector,
                                                             int refTotal, const unsigned exampleNumber) {

        if (save) {
            std::string refTag = util::to_string(exampleNumber);
            std::string resultsFile = "results-" + refTag + ".txt";

            std::ofstream file((this->resultFolder + "/" + resultsFile).c_str());
            if (!file.is_open())
                std::cout << "Problem opening the result file!" << resultsFile << std::endl;
            else {
                file << "vDOFs = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << vDOFs[refCounter] << ", ";
                file << "]; \n";

                file << "yDOFs = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << yDOFs[refCounter] * this->dim << ", ";
                file << "]; \n";

                file << "wDOFs = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << wDOFs[refCounter] << ", ";
                file << "]; \n";

                file << "ed = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << eSpaceTimeSpaceGradVector[refCounter] << ", ";
                file << "]; \n";

                file << "eh = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << eSpaceTimeVector[refCounter] << ", ";
                file << "]; \n";

                file << "maj = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << majVector[refCounter] << ", ";
                file << "]; \n";

                file << "majII = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << majIIVector[refCounter] << ", ";
                file << "]; \n";

                file << "majh = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << majhVector[refCounter] << ", ";
                file << "]; \n";

                file << "esol = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << eSpaceTimeSolOperVector[refCounter] << ", ";
                file << "]; \n";

                file << "id = [";
                for (int refCounter = 0; refCounter < refTotal; ++refCounter)   file << eIdentVector[refCounter] << ", ";
                file << "]; \n";

                file.close();

            }

        }
    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsSaveToFileDivDivMMMatrices(bool save,
                                                         gsSparseMatrix<real_t> &DivDiv, gsSparseMatrix<real_t> &MM,
                                                         int refCounter,
                                                         const unsigned exampleNumber) {
        {
            //if (save && refCounter <= 3 && ! this->isAdaptive) {
            if (save && refCounter <= 3) {
                // Saving Paraview files
                std::string refTag = "-refnum-" \
                                     + util::to_string(refCounter) \
                                     + "-example-"
                                     + util::to_string(exampleNumber);
                // Saving txt files
                std::string divdivFile = "DivDiv" + refTag + ".txt";
                std::string mmFile = "MM" + refTag + ".txt";

#pragma omp parallel sections
                {
#pragma omp section
                    {
                        std::ofstream streamDivDiv((this->resultFolder + "/" + divdivFile).c_str());
                        if (!streamDivDiv.is_open())
                            std::cout << "Problem opening the result file!" << divdivFile << std::endl;
                        else {
                            std::cout << std::scientific;
                            streamDivDiv << DivDiv.toDense() << " ";
                            std::cout.unsetf(std::ios::scientific);
                            streamDivDiv.close();

                        }
                    }
#pragma omp section
                    {

                        std::ofstream streamMM((this->resultFolder + "/" + mmFile).c_str());
                        if (!streamMM.is_open())
                            std::cout << "Problem opening the result file!" << mmFile << std::endl;
                        else {
                            std::cout << std::scientific;
                            streamMM << MM.toDense();
                            std::cout.unsetf(std::ios::scientific);
                            streamMM.close();
                        }
                    }
                }
            }
        }

    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsSaveToFileKMatrix(bool save,
                                                gsSpaceTimeAssembler<real_t> &assembler, int refCounter,
                                                const unsigned exampleNumber) {
        {
            if (save) {
                std::string refTag = "-refnum-" + util::to_string(refCounter) + "-example-"
                                     + util::to_string(exampleNumber);
                // Saving txt files
                std::string kFile = "K" + refTag + ".txt";
                std::ofstream streamDivDiv((this->resultFolder + "/" + kFile).c_str());
                if (!streamDivDiv.is_open())
                    std::cout << "Problem opening the result file!" << kFile << std::endl;
                else {
                    streamDivDiv << assembler.matrix().toDense() << " ";
                    streamDivDiv.close();
                }
            }
        }

    }

    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsSolveKhfhSystem(gsSpaceTimeAssembler<real_t>& assembler,
                                              gsMatrix<>& vVector,
                                              gsMatrix<double>& timeSolvV,
                                              int refCounter,
                                              gsVector<real_t>& stopcritVector) {

        real_t timeLinearSolver(0.0), timeDirect(0.0);
        gsMatrix<> vVectorIter0, vVectorIter, vVectorDir;

        int N = assembler.rhs().size();
        real_t hPowD = std::pow((double) N, (double)(- 1));
        index_t maxItLinSolver = 4 * N;
        real_t tolLinSolver = (refCounter > 1) ?
                              stopcritVector[refCounter-2] :
                              hPowD * std::pow(10.0, -2);

        if (assembler.numDofs()) {

            if (refCounter == 0) {
                gsSparseSolver<>::LU solverLU;
                solverLU.compute(assembler.matrix());
                vVectorIter0 = solverLU.solve(assembler.rhs());
                this->gsSetVVector(vVectorIter0);

            } else {

                vVectorIter0 = this->getVRefVector();
                this->gsSetVVector(vVectorIter0);
            }
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    gsCPUStopwatch clock_1;
                    clock_1.restart();
                    //gsInfo << "vVectorIter0 : \n" << vVectorIter0 << "\n";

                    gsSparseSolver<>::BiCGSTABIdentity solverIter;
                    solverIter.setMaxIterations(maxItLinSolver);
                    solverIter.setTolerance(tolLinSolver);
                    solverIter.compute(assembler.matrix());

                    real_t error0 = (assembler.matrix()*vVectorIter0 - assembler.rhs()).norm()/assembler.rhs().norm();
                    gsInfo << "error of vVectorIter0 = " << error0 << "\n";

                    vVectorIter = solverIter.solveWithGuess(assembler.rhs(), vVectorIter0);
                    timeLinearSolver = clock_1.stop();
                    ///*
                    real_t error = (assembler.matrix()*vVectorIter - assembler.rhs()).norm()/assembler.rhs().norm();
                    gsInfo << solverIter.detail() << "\n";
                    gsInfo << "error of vVectorIter = " << error << "\n";
                    gsInfo << "solverIter.error()   = " << solverIter.error() << "\n";

                    if ( solverIter.error() <= solverIter.tolerance() && error <= solverIter.tolerance() ) {
                        gsInfo <<" Test for v passed!\n";
                        vVectorDir = vVectorIter;
                    }
                    else {
                        gsInfo << " Test for v failed!\n";
                        clock_1.restart();

                        // DLDT is for SPD matrices - doesn't work
                        // QR is for any, rectangular matrices - works slower
                        // LU is for square matrices
                        gsSparseSolver<>::LU solverDirect;
                        solverDirect.compute(assembler.matrix());
                        vVectorDir = solverDirect.solve(assembler.rhs());
                        //real_t error = (assembler.matrix()*vVectorDir - assembler.rhs()).norm()/assembler.rhs().norm();
                        //gsInfo << "error of vVectorDir = " << error << "\n";
                        timeDirect = clock_1.stop();
                    }

                }

                #pragma omp section
                {
                    //if (this->dim == 2) {
                        gsCPUStopwatch clock_2;
                        clock_2.restart();

                        // DLDT is for SPD matrices - doesn't work
                        // QR is for any, rectangular matrices - works slower
                        // LU is for square matrices
                        gsSparseSolver<>::LU solverDirect;
                        solverDirect.compute(assembler.matrix());
                        vVectorDir = solverDirect.solve(assembler.rhs());
                        //real_t error = (assembler.matrix()*vVectorDir - assembler.rhs()).norm()/assembler.rhs().norm();
                        //gsInfo << "error of vVectorDir = " << error << "\n";
                        timeDirect = clock_2.stop();
                    //}
                }

            }
            timeSolvV(refCounter, 0) = timeDirect;
            timeSolvV(refCounter, 1) = timeLinearSolver;

            vVector = vVectorDir;  // return as the direct solver solution as a results

        }

    }


    template<unsigned d>
    void gsTestSpaceTimeMajorant<d>::gsSolveKhwhfhSystem(gsSpaceTimeAssembler<real_t>& assembler,
                                                gsMatrix<>& wVector, gsMatrix<> & timeSolvW,
                                                int refCounter, gsVector<real_t>& stopcritVector) {

        real_t timeLinearSolver(0.0), timeDirect(0.0);
        gsMatrix<> wVectorIter0, wVectorIter, wVectorDir;

        int N = assembler.rhs().size();
        real_t hPowD = std::pow((double) N, (double)(- 1));
        index_t maxItLinSolver = 4 * N;
        real_t tolLinSolver = (refCounter > 1) ?
                              stopcritVector[refCounter-2] * std::pow(10.0, -4) :
                              hPowD * std::pow(10.0, -4);

        if (assembler.numDofs())
        {
            if (refCounter == 0) {
                gsSparseSolver<>::LU solverLU;
                solverLU.compute(assembler.matrix());
                wVectorIter0 = solverLU.solve(assembler.rhs());
                //this->gsSetWVector(wVectorIter0);

            } else {
                wVectorIter0 = this->getWRefVector();
                //this->gsSetWVector(wVectorIter0);
            }
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    gsCPUStopwatch clock_1;
                    clock_1.restart();
                    gsSparseSolver<>::BiCGSTABIdentity solverIter;
                    solverIter.setMaxIterations(maxItLinSolver);
                    solverIter.setTolerance(tolLinSolver);
                    solverIter.compute(assembler.matrix());

                    real_t error0 = (assembler.matrix()*wVectorIter0 - assembler.rhs()).norm()/assembler.rhs().norm();
                    gsInfo << "error of wVectorIter0 = " << error0 << "\n";

                    wVectorIter = solverIter.solveWithGuess(assembler.rhs(), wVectorIter0);
                    timeLinearSolver = clock_1.stop();

                    real_t error = (assembler.matrix()*wVectorIter - assembler.rhs()).norm()/assembler.rhs().norm();
                    gsInfo << solverIter.detail() << "\n";
                    gsInfo << "error of wVectorIter = " << error << "\n";
                    gsInfo << "solverIter.error()   = " << solverIter.error() << "\n";

                    if ( solverIter.error() <= solverIter.tolerance() && error <= solverIter.tolerance() ) {
                        gsInfo << " Test for w passed!\n";
                        wVectorDir = wVectorIter;
                    }
                    else   {
                        gsInfo <<" Test for w failed : \n"
                               << solverIter.error() << " <= (?) " << solverIter.tolerance() << "\n"
                               <<  error << " <= (?) " << solverIter.tolerance() << "\n";
                        clock_1.restart();
                        // DLDT is for SPD matrices - doesn't work
                        // QR is for any, rectangular matrices - works slower
                        // LU is for square matrices
                        gsSparseSolver<>::LU solverDirect;
                        solverDirect.compute(assembler.matrix());
                        wVectorDir = solverDirect.solve(assembler.rhs());
                        timeDirect = clock_1.stop();
                    }
                }
                /*
                #pragma omp section
                {
                    gsCPUStopwatch clock_2;
                    clock_2.restart();

                    gsMatrix<> w0(wVectorIter0.rows(), wVectorIter0.cols());
                    w0.setZero();

                    gsSparseSolver<>::BiCGSTABIdentity solverIterWithZeroInitGuess;
                    solverIterWithZeroInitGuess.setMaxIterations(maxItLinSolver);
                    solverIterWithZeroInitGuess.setTolerance(tolLinSolver * std::pow(10.0, -2));
                    solverIterWithZeroInitGuess.compute(assembler.matrix());

                    real_t error0 = (assembler.matrix()*w0 - assembler.rhs()).norm()/assembler.rhs().norm();
                    gsInfo << "error of wVectorIter0 = " << error0 << "\n";

                    wVectorIter = solverIterWithZeroInitGuess.solveWithGuess(assembler.rhs(), w0);

                    real_t error = (assembler.matrix()*wVectorIter - assembler.rhs()).norm()/assembler.rhs().norm();
                    //gsInfo << solverIterWithZeroInitGuess.detail() << "\n";
                    gsInfo << "error of wVectorIter = " << error << "\n";
                    gsInfo << "solverIter.error()   = " << solverIterWithZeroInitGuess.error() << "\n";

                    if ( solverIterWithZeroInitGuess.error() <= solverIterWithZeroInitGuess.tolerance()
                         && error <= solverIterWithZeroInitGuess.tolerance() ) gsInfo <<" Test for v passed!\n";
                    else    gsInfo <<" Test for w failed!\n";

                    wVectorDir = wVectorIter;
                    timeDirect = clock_2.stop();
                }
                */

                #pragma omp section
                {
                    //if (this->dim == 2) {
                        gsCPUStopwatch clock_2;
                        clock_2.restart();
                        // DLDT is for SPD matrices - doesn't work
                        // QR is for any, rectangular matrices - works slower
                        // LU is for square matrices
                        gsSparseSolver<>::LU solverDirect;
                        solverDirect.compute(assembler.matrix());
                        wVectorDir = solverDirect.solve(assembler.rhs());
                        timeDirect = clock_2.stop();
                    //}
                }

            }
            timeSolvW(refCounter, 0) = timeDirect;
            timeSolvW(refCounter, 1) = timeLinearSolver;
            /*
            if (solverLinear.iterations() < maxItLinSolver) vVector = vVectorIter;
            else vVector = vVectorDir;
             */
            wVector = wVectorDir;  // return as the direct solver solution as a results
        }
    }

    template <unsigned d>
    void gsTestSpaceTimeMajorant<d>:: gsSolveMajorantOptimalSystem(gsSparseMatrix<real_t> & yM, gsMatrix<real_t> & yRhs, gsMatrix<real_t> & yVector,
                                                          gsMatrix<double> & timeSolvY, int refCounter, int iterMaj, gsVector<index_t> &yDOFs,
                                                          gsVector<real_t> & stopcritVector)
    {

        real_t timeCGSolver(0.0), timeDirectSolver(0.0);
        gsMatrix<> yVectorIter, yVectorDir;

        int N = yRhs.size();    // new size corresponding to d-1 blocks
        real_t hPowD = std::pow((double) N, -1.0);
        //real_t tolCGY = hPowD * std::pow(10.0, -6);
        real_t tolCGY = (refCounter > 1) ?
                        stopcritVector[refCounter-2] * std::pow(10.0, -2) :
                        hPowD * std::pow(10.0, -6);
        index_t maxIters = 3 * N;

        //gsInfo << "N = " << N << "\n";
        //gsInfo << "yDOFs[refCounter] = " << yDOFs[refCounter] << "\n";

        gsOptionList opt = gsIterativeSolver<real_t>::defaultOptions();
        opt.setInt ("MaxIterations", maxIters);
        opt.setReal("Tolerance"    , tolCGY);
        gsLinearOperator<>::Ptr preConMat = gsIdentityOp<>::make(N);

        if (refCounter == 0 && iterMaj == 0) {
            gsSparseSolver<>::LU solverLU;
            solverLU.compute(yM);
            yVectorIter = solverLU.solve(yRhs);

        } else{
            // we are getting here the vector of the size [(N_big / d) x d]
            // N = (d - 1) / d * N_big
            // but we need only d-1 columns of it since the last dimension is zero
            yVectorIter = (this->gsGetYRefVector()).block(0, 0, yDOFs[refCounter], this->dim-1);
            // we resize it to work with [N x 1] vector
            yVectorIter.resize(yDOFs[refCounter] * (this->dim - 1), 1);
        }

        if (this->dim == 2 || this->dim == 3) {
            #pragma omp parallel sections
            {

                #pragma omp section
                {
                //if (iterMaj == 0){
                    gsCPUStopwatch clock_1;
                    clock_1.restart();
                    gsSparseSolver<>::SimplicialLDLT solverSLDLT;
                    solverSLDLT.compute(yM);
                    yVectorDir = solverSLDLT.solve(yRhs);
                    timeDirectSolver = clock_1.stop();
                }

                #pragma omp section
                {
                    gsCPUStopwatch clock_2;
                    clock_2.restart();
                    gsConjugateGradient<> solverCG(yM, preConMat);
                    solverCG.setOptions(opt);
                    solverCG.solve(yRhs, yVectorIter);
                    timeCGSolver = clock_2.stop();
                    gsInfo << "gsConjugateGradient for Y: " << solverCG.detail() << "";

                    real_t error = (yM*yVectorIter - yRhs).norm() / yRhs.norm();
                    gsInfo << "error of yVectorIter = " << error << "\n";
                    gsInfo << "solverIter.error()   = " << solverCG.error() << "\n";
                    if ( solverCG.error() <= solverCG.tolerance() && error <= solverCG.tolerance() ) gsInfo <<" Test for y passed!\n";
                    else gsInfo <<" Test for y failed!\n";
                }
            }
        }
        else if (this->dim == 4)    {
            gsCPUStopwatch clock_2;
            clock_2.restart();
            gsConjugateGradient<> solverCG(yM, preConMat);
            solverCG.setOptions(opt);
            solverCG.solve(yRhs, yVectorIter);
            yVector = yVectorIter;
            timeCGSolver = clock_2.stop();
            gsInfo << "solver for Y: " << solverCG.detail() << "";

            real_t error = (yM*yVectorIter - yRhs).norm() / yRhs.norm();
            gsInfo << "error of yVectorIter = " << error << "\n";
            gsInfo << "solverIter.error()   = " << solverCG.error() << "\n";

            if ( solverCG.error() <= solverCG.tolerance() && error <= solverCG.tolerance() ) gsInfo <<" Test for y passed!\n";
            else {
                gsInfo <<" Test for y failed!\n";
                clock_2.restart();
                gsSparseSolver<>::SimplicialLDLT solverSLDLT;
                solverSLDLT.compute(yM);
                yVector = solverSLDLT.solve(yRhs);
                timeDirectSolver = clock_2.stop();
            }

            yVectorDir = yVector;
        }
        timeSolvY(refCounter, 0) += timeDirectSolver;
        timeSolvY(refCounter, 1) += timeCGSolver;

        //yVector = yVectorIter;
        yVector = yVectorDir;
    }

    template <unsigned d>
    void gsTestSpaceTimeMajorant<d>:: interpolateToRefVectorWithDirichletBC(const gsSpaceTimeAssembler<real_t> & assembler,
                                                                           const gsMultiPatch<> &mpV,
                                                                           gsMatrix<> &interpPoints,
                                                                           gsMatrix<> &freeInterpVector) {

        gsMatrix<real_t> interpVector = mpV.patch(0).eval(interpPoints);
        interpVector.transposeInPlace();
        //gsInfo << "full interpolated vector: \n" << interpVector << "\n";

        // Reconstruct solution coefficients on patch p
        int unk = 0;

        const gsDofMapper &mapper = assembler.system().colMapper(unk);
        freeInterpVector.resize(mapper.freeSize(), interpVector.cols());

        index_t j = 0;
        for (index_t i = 0; i < interpVector.size(); ++i) {
            if (mapper.is_free(i, 0)) freeInterpVector.row(j++) = interpVector.row(i);
        }
        //gsInfo << "free coeffs of interpolated vector: \n" << freeInterpVector << "\n";

    }
}

