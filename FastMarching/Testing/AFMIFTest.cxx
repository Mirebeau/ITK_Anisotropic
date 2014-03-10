#include "itkImage.h" 
#include <iostream>
#include <ctime>



int AFMIFTest(int argc, char* argv[])
{
    auto Im = itk::Image<double,3>::New();
    
    std::cout << "Hello world" << std::endl;
    return EXIT_SUCCESS;
    
}

int main(int argc, char* argv[]){
    return AFMIFTest(argc,argv);
}