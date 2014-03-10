
 
#define DEBUG 

#include "itkImage.h"
#include <iostream>
#include <ctime> 
 
using std::cout;
using std::endl;

#include "DiffusionTest.h"

int main()
{
    
#ifdef DEBUG
    std::cout << "DEBUG" << std::endl;
#endif
    Testing::CoherenceEnhancingDiffusion2D();
    Testing::EdgeEnhancingDiffusion2D();
    Testing::CoherenceEnhancingDiffusion3D();
  return 0;
}
