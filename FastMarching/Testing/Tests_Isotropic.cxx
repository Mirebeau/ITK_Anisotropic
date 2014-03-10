#include "itkImage.h" 
#include <iostream>
#include <ctime>

#include "AFMIFTest.cxx"

int main(int argc, char* argv[])
{
    AFMIFTest(argc, argv);
    
    auto Im = itk::Image<double,3>::New();
    
    std::cout << "Hello world" << std::endl;
    return EXIT_SUCCESS;
    
}
