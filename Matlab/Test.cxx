/*
 *  HelloLib.cp
 *  HelloLib
 *
 *  Created by Jean-Marie Mirebeau on 02/04/13.
 *  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
 *
 */

#include <iostream>

//#include "MatlabITKFastMarching.h"
#include <sstream>

//#include "IsotropicFastMarching.h"
//#include "AnisotropicFastMarching.h"
//#include "MexMessageWrapper.h"

#include "Finsler2DNorm.h"
#include "Riemannian3DNorm.h"
#include "UnrollingConstructors.h"

int main(int argc, char* argv[])
{
    
    
    
    
    
    return EXIT_SUCCESS;
    
    using namespace UnrollingConstructors;
    using std::cout;
    using std::endl;
    
//    double V[3];
//    Set(V,1,2,3);
    
    itk::Vector<double,3> V;
    Set(V,1,2,3);
    
    cout << Make<itk::Vector<double> > (1,2,3) << endl;
    
    cout << V[0] << V[1] << V[2];
    
//    UnrollingConstructors::printf("Hello %i",2);
    
//    auto MyVec = Make<itk::Vector<double,3> > (1,2,3);
//    itk::Vector<double,3> MyVec = {1,2,3};
    
//    cout << MyVec << endl;
    
    
    
    
    return EXIT_SUCCESS;
    
    using namespace std;

    typedef itk::Finsler2DNorm<double> NormType;
    
    NormType N;
    N.SetNthComponent(0,1);
    N.SetNthComponent(1,0);
    N.SetNthComponent(2,1);
    N.SetNthComponent(3,0.5);
    N.SetNthComponent(4,0.5);
    
    cout << N.GetM() << " " << N.GetOmega() << endl;
    typedef itk::Vector<double,NormType::Dimension> VectorType;
    
    
    std::vector<NormType::CompressedOffsetType> l;
    N.GetCompressedStencil(l);
    
    for(auto it=l.begin(); it!=l.end(); ++it)
        cout << VectorType(*it) << endl;
    
    return EXIT_SUCCESS;
}