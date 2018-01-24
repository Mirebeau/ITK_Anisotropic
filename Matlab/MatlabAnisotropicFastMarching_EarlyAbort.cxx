#include <iostream>

#include "MexInterface.h"
#include "AnisotropicFastMarching_EarlyAbort.h"

#include "Riemannian2DNorm.h"
#include "Finsler2DNorm.h"
#include "AsymmetricQuadratic2DNorm.h"
#include "Riemannian3DNorm.h"
#include "ExtendedNorm.h"


EXTERN_C
int __mexFunction__(int nlhs, mxArray *plhs[],
                    int nrhs, const mxArray *prhs[] ){
    MexMsg() << "AnisotropicFastMarching_EarlyAbort.\n";
    
    if(nrhs!=1) {
        MexWarnMsg() << "Exactly one input parameter needed";
        MexMsg() << "Input fields description \n"
        << "[NormType] (mandatory, string) among\n"
        << "{Riemannian2DNorm, Finsler2DNorm,  Riemannan3DNorm, Riemannian2DNormExt, Finsler2DNormExt, Riemannian3DNormExt}\n"
        
        << "\n[Metric] (mandatory) \n"
        << "Riemannian metric format : array of size  d*(d+1)/2,  [M11;M12; ... ; M1d; M22; ...; M2d; M33; ... ; Mdd].\n\n"
        << "Special metric types : \n"
        << "Riemannian metric made of a dxd and a 1x1 diagonal block : append corresponding coefficient at the end\n"
        << "2d asymmetric metric, sqrt(u'.M.u) - omega'.u : coefficients of M followed by the coefficients of omega\n"
        
        << "\nGeodesic data : \n"
        << "[Seeds] (mandatory) with format [x1; ... ; xd; value] \n"
        << "[Tips]  (optional)  with format [x1; ... ; xd] \n"
        
        << "\nGeometric data :\n"
        << "[Origin] (optional) default is [0,...,0]"
        << "[Spacing] (optional) default is [1,...,1]"
        << "[TransposeFirstTwoImageCoordinates] (optional) default is 1, following Matlab's convention"
        
        << "\nStopping criteria\n"
        << "[StopWhenTipsAreReached] (optional) \n"
        << "[StopAtDistance] (optional) refers to the geodesic distance\n"
        
        << "\nGeodesic euclidean length"
        << "[StopAtEuclideanDistance] (optional) refers to the euclidean length of the computed geodesic \n"
        << "[EuclideanSpacing] (optional) default is [Spacing]"
        
        << "\n[CheckData]"
        << "\n";
        // To do : export gradient ?
        return EXIT_SUCCESS;
    }
    
    if(nlhs!=1) {
        MexWarnMsg() << "Exactly one output parameter needed";
        MexMsg() << "Output fields description\n"
        << "[Distance] the geodesic distance\n"
        << "\nIf any tips specified\n"
        << "[Paths] the geodesic paths corresponding to the tips\n"
        << "\nIf StopAtEuclideanDistance specified\n"
        << "[EuclideanStoppingIndex]"
        << "[EuclideanDistance]"
        << "[GeodesicOfEuclideanStoppingIndex]"
        << "\n";
        return EXIT_SUCCESS;
    }

    mxIO IO(prhs[0],plhs);
    
    typedef itk::Riemannian2DNorm<double>   Riemannian2DNorm;
    typedef itk::Finsler2DNorm<double>      Finsler2DNorm;
    typedef itk::Riemannian3DNorm<double>   Riemannian3DNorm;
    typedef itk::AsymmetricQuadratic2DNorm<double> AsymQuad2DNorm;
    
    typedef itk::ExtendedNorm<Riemannian2DNorm> Riemannian2DNormExt;
    typedef itk::ExtendedNorm<Finsler2DNorm>    Finsler2DNormExt;
    typedef itk::ExtendedNorm<Riemannian3DNorm> Riemannian3DNormExt;
    
    try {
        const char * normName = IO.GetStringField("NormType");

        const int result =
        !strcmp(normName,"Riemannian2DNorm")    ? AnisotropicFastMarching_EarlyAbort<Riemannian2DNorm>(IO) :
        !strcmp(normName,"Finsler2DNorm")       ? AnisotropicFastMarching_EarlyAbort<Finsler2DNorm>(IO) :
        !strcmp(normName,"AsymmetricQuadratic2DNorm") ? AnisotropicFastMarching_EarlyAbort<AsymQuad2DNorm>(IO) :
        !strcmp(normName,"Riemannian3DNorm")    ? AnisotropicFastMarching_EarlyAbort<Riemannian3DNorm>(IO) :
        !strcmp(normName,"Riemannian2DNormExt") ? AnisotropicFastMarching_EarlyAbort<Riemannian2DNormExt>(IO) :
        !strcmp(normName,"Finsler2DNormExt")    ? AnisotropicFastMarching_EarlyAbort<Finsler2DNormExt>(IO) :
        !strcmp(normName,"Riemannian3DNormExt") ? AnisotropicFastMarching_EarlyAbort<Riemannian3DNormExt>(IO) :
        (MexWarnMsg() << "Unsupported norm type", EXIT_FAILURE);
        
        if(result!=EXIT_SUCCESS) MexWarnMsg() << "Failed (AnisotropicFastMarching_EarlyAbort)";
    } catch (std::exception & e) {
        MexWarnMsg() << "Exception caught !\n" << e.what() << "\n";
    }

    return EXIT_SUCCESS;
    }
