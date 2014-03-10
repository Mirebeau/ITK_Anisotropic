#include "itkImage.h" 
#include <iostream>
#include <ctime>

using std::cout;
using std::endl;

#include "Tutorial.h"
#include "Tutorial3D.h"
//#include "TubularSegmentation2DTutorial.h" //Needs Fethallah's code.
#include "RescaleAndExport_Function.h"

#include "CurveNeighborhoodSource.h"
#include "Topographic_Test.h"
#include "Segmentation_Test.h"
#include "Tubular3D_Test.h"
#include "GenerateTestImages.h"

#include "Segmentation_Image.h"
#include "Topographic_Image.h"
#include "TubeCenterline2D_Image.h"
//#include "TubularBand_Test.h"
//#include "Topographic_Test.h"

int main(int argc, char* argv[])
{
    std::cout << "Hello world" << std::endl;
    
//    AnisotropicFastMarching_Segmentation_Test::GenerateImage();
//    AnisotropicFastMarching_Segmentation_Test::RunTest();
    
//    AnisotropicFastMarching_Topographic_Test::GenerateImage();
//    AnisotropicFastMarching_Topographic_Test::RunTest();
    
    AnisotropicFastMarching_TubeCenterline2D_Test::GenerateImage_Tree();
    AnisotropicFastMarching_TubeCenterline2D_Test::RunTest_Tree();
    
//    AnisotropicFastMarching_TubeCenterline2D_Test::RunTest_Multiscale();
    
    AnisotropicFastMarching_TubeCenterline2D_Test::GenerateImage();
//    AnisotropicFastMarching_TubeCenterline2D_Test::RunTest_ConstantRadius();
    
//    auto source = itk::CurveNeighborhoodSource<itk::Image<float,2> >::New();
//    GenerateSegmentationImage();
    
    
    //Tests::RescaleAndExport_Test();
//    auto imageSource = itk::TubularBandSource<float,2>::New();
    //TubularBand_3D_Test();
//    Topographic_Test();
    
    /*
     double F[6]={0.000232955113122171,0.000518318009049774,0.000194200676670343,0.00177825104072398,-5.91045537692347e-05,0.000553855015384615};
     

    
    typedef itk::Riemannian3DNorm<double> NormType;
    NormType M;
    for(int i=0; i<6; ++i) M.SetNthComponent(i, F[i]);
    
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            cout << M(i,j) << "\n";
    
    std::vector<NormType::CompressedOffsetType> l;
    M.GetCompressedStencil(l);
    
    
    typedef NormType::VectorType VectorType;
    cout << "\nl\n";
    for(int i=0; i<3; ++i) cout << VectorType(l[i]) << "\n";
    
    NormType::Stencil stencil(&l[0],l.size());
    cout << "\nstencil\n";
    for(int i=0; i<14; ++i) cout << VectorType(stencil[i]) << "\n";
    
    return EXIT_SUCCESS;
    */
    
    
//    return TutorialRiemannian3D();
//    return Finsler2DTutorial();
//    return TubularSegmentation2DTutorial(argc, argv);
//    return AFMIFTest2(argc, argv);
}


