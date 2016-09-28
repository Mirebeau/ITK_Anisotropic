#include <iostream>
#include <ctime> 

#ifndef DEBUG
#define NDEBUG
#endif

#include "LinearAnisotropicDiffusionCommandLine.h"
#include "CoherenceEnhancingDiffusionCommandLine.h"

#include "itkImage.h"
#include "AnisotropicDiffusionLBRImageFilter.h"
#include "CoherenceEnhancingDiffusionFilter.h"
#include "UnaryFunctorWithIndexImageFilter.h"
#include "SyntheticImageSource.h"
#include "DiffusionTest.h"
#include "DiffusionBasicTest.h"

using std::cout;
using std::endl;

/*
#include "DiffusionOperator.h"
#include "DiffusionTest.h"

#include "RescaleAndExportBMP.h"

#include "EdgeEnhancingDiffusionTensor.h"
 */

int main(int argc, char ** argv)
{
#ifdef DEBUG
    std::cout << "DEBUG" << std::endl;
#endif
  
    
//    Testing::SyntheticOscillations();
//    Testing::SyntheticTriangle();

//    Testing::CircleVectorField();
//    Testing::SyntheticTriangle();
    Testing::EdgeEnhancingDiffusion3D();
    
    Testing::Assemble3DSkull();
    //    Testing::NonUniformSpacing();
   
    
//    Testing::CoherenceEnhancingDiffusion2D();
//    Testing::EdgeEnhancingDiffusion2D();
    
//    DiffusionBasicTest::LinearDiffusion();
    //    Testing::Color();
    
    /*
    const int argc_ = 5+1;
    char * argv_[10];
//    ./itkDiffusion ../Images/VectorField_CircleOpposites.hdf5 VectorField_Circle_cEED.hdf5 10 0.05 cEED
    argv_[1] = "../Images/VectorField_CircleOpposites.hdf5";
    argv_[2] = "VectorField_Circle_cEED.hdf5";
    argv_[3] = "10";
    argv_[4] = "0.05";
    argv_[5] = "cEED";
    return CoherenceEnhancingDiffusionCommandLine::Execute(argc_, argv_);
    */
    
    try {
        //LinearAnisotropicDiffusionCommandLine::Execute(argc, argv);
        CoherenceEnhancingDiffusionCommandLine::Execute(argc, argv);
    } catch (std::exception & e) {
        std::cerr << "Exception caught !\n" << e.what() << "\n";
    }
    return EXIT_SUCCESS;
    
    
}

