#include "itkImage.h" 
#include <iostream>
#include <ctime>

using std::cout;
using std::endl;

//#include "RunAFMBench.h"


#include "CommandLineCall2.h"
#include "NonRegressionTests.h"


int main(int argc, char* argv[])
{

    return CommandLineMain(argc,argv);
    
//    ITKFM_NonRegressionTests::CompareTimings();
//    ITKFM_NonRegressionTests::CompareTiming< itk::SSRT2<float> >(100);
    
//    itk::Test<10>();
//    itk::Test<11>();
//    itk::Test<12>();
//    itk::Test<15>();
    
    //GenerateTestData::GenerateIdentityTensors();
    //NonRegressionTests::IdentityTest();
    return EXIT_SUCCESS;
    
    const int argc_ = 5;
    char * argv_[argc_] = {
        argv[0],
        (char*)"data/IdentityTensor.nii",
        (char*)"0",
        (char*)"output/Geodesic.mathematica",
        (char*)"output/DistanceImage.nii"
    };
    
    argc = argc_;
    argv = argv_;
    
    cout << argc << endl;
    for(int i=0; i<argc; ++i) cout << argv[i] << endl;
    cout << endl;
    
    return CommandLineMain(argc,argv);
    
}
