#include <iostream>

#define NDBUG

#include "MexInterface.h"
#include "AnisotropicDiffusion.h"

#include "Riemannian2DNorm.h"
#include "Riemannian3DNorm.h"

EXTERN_C
int __mexFunction__(int nlhs, mxArray *plhs[],
                    int nrhs, const mxArray *prhs[] ){
    
    if(nrhs!=1 || nlhs!=1){
        /*
        << "Riemannian metric format : Metric(:,i,j) = [xx; xy; yy] are the symmetric matrix coefficients at pixel (i,j)\n"
        << "Two and three dimensional metrics are supported\n"
        << "Multi-valued image support : Image(:,i,j) is the multivalued pixel of index (i,j)\n\n"
        ;*/

        MexWarnMsg() << "Anisotropic Diffusion : exactly one input and one output (structures).";
        MexMsg()
        << "Fields usage : "
        << "ImageDimension (defaults to 2)\n"
        << "PixelDimension (0 indicates scalar images, defaults to 0)\n";
        return EXIT_SUCCESS;
    }
    
    mxIO IO(prhs[0],plhs);

    const int imageDimension = IO.GetObject<double>("ImageDimension", 2);
    
    /* // For now, we forget about pixel dimension. Just call several times the linear filter if needed.
    const int pixelDimension = IO.GetObject<double>("PixelDimension",0);
    if(pixelDimension<0){
        MexWarnMsg() << "Field PixelDimension should be >=1";
    }*/
    
    typedef itk::Riemannian2DNorm<double> Riemannian2DNorm;
    typedef itk::Riemannian3DNorm<double> Riemannian3DNorm;
    
//    typedef itk::SymmetricSecondRankTensor<double,2> SSRT2; //itk::SSRT2<double>
//    typedef itk::SymmetricSecondRankTensor<double,3> SSRT3;
    
    try {
        switch (imageDimension) {
            case 2: return AnisotropicDiffusion<Riemannian2DNorm>(IO);
            case 3: return AnisotropicDiffusion<Riemannian3DNorm>(IO);
                                
            default:
                MexWarnMsg() << "Field ImageDimension should indicate 2 or 3\n";
                return EXIT_FAILURE;
        }
    } catch (std::exception & e) {
            MexWarnMsg() << "Exception caught !\n" << e.what() << "\n";
    }

/*} catch(itk::ExceptionObject & err) {
        MexWarnMsg() << "ITK exception ! \n" << err;
    } catch(...) {
        MexWarnMsg() << "Unknown exception caught.";
    }*/
    
    return EXIT_FAILURE;
}