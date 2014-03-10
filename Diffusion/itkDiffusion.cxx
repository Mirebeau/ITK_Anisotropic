#define DEBUG

#include <iostream>
#include <ctime> 
#include "itkImage.h"
#include "AnisotropicDiffusionLBRImageFilter.h"
#include "CoherenceEnhancingDiffusionFilter.h"
#include "UnaryFunctorWithIndexImageFilter.h"
#include "SyntheticImageSource.h"
#include "DiffusionTest2.h"

using std::cout;
using std::endl;

/*
#include "DiffusionOperator.h"
#include "DiffusionTest.h"

#include "RescaleAndExportBMP.h"

#include "EdgeEnhancingDiffusionTensor.h"
 */


template<int repeats>
bool TestMultiThreadAccess(size_t s);


int main()
{
#ifdef DEBUG
    std::cout << "DEBUG" << std::endl;
#endif
    /*
    itk::Riemannian2DNorm<double, int> m;
    m(0,0)=0.6704657109251299;
    m(0,1)=0.3716795467893859;
    m(1,1)=0.2191640538129669;
    
    itk::Riemannian2DNorm<double, int>::VectorType superBase[3];
    for(int i=0; i<2; ++i)
        for(int j=0; j<2; ++j)
            superBase[i][j]=(i==j);
    m.ObtuseSuperBase(superBase);
    std::cout << superBase[0] << superBase[1] << std::endl;
    */
//    return 0;
//    Testing::StructureTensor();
//    Testing::LinearDiffusion();
//    Testing::FingerPrint();
    Testing::EdgeEnhancingDiffusion2D();
//    TestMultiThreadAccess<10>(10);
    return 0;
    /*
    itk::Test<1>();
//    itk::TestDiffusion3D();
//    itk::Test<3>();
    return 0;
    
    //itk::TestDiffusion2D();
    //itk::TestEED2D();
    //itk::Test<1>();
    //itk::TestEigenVals();
    
    //itk::TestDiffusion3D();
    for(int i=0; i<100; ++i)
        itk::Test<2>();
  return 0;
     */
}

// Just checking that (parallel) update of an array does not work.


#include <random>       // std::default_random_engine
#include <vector>
#include "itkImageRegionIterator.h"
#include "itkUnaryFunctorImageFilter.h"

template<int repeats>
struct AssignFunctorType {
    typedef itk::Vector<int,repeats> IndicesVectorType;
    int * vals;
    int k;
    IndicesVectorType operator()(IndicesVectorType v) const {
        IndicesVectorType output;
        for(int i=0; i<repeats; ++i){
            output[i] = v[i]*k;
//            std::cout << "("<<v[i]<<")" << std::endl;
            int & val = vals[v[i]];
            ++val;
//            std::cout << val << std::endl;
//            ++(vals[v[i]]);
        }
        return v;
    }
};

template<int repeats>
bool TestMultiThreadAccess(size_t s){
    
    static const int Dimension=2;
    
    size_t s2=s*s;
    std::vector<int> values;
    std::vector<int> indices;
    values.resize(s2,0);

    indices.reserve(s2*repeats);
    for(int i=0; i<s2; ++i)
        for(int j=0; j<repeats; ++j)
            indices.push_back(i);
    
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    auto r = std::default_random_engine(seed);
    std::shuffle(indices.begin(), indices.end(), r);

/*    for(auto it=indices.begin(); it!=indices.end(); ++it) {
        std::cout << (*it) << std::endl;
    }*/
    
    
    typedef itk::Vector<int,repeats> IndicesVectorType;
    typedef itk::Image<IndicesVectorType> IndicesImageType;
    
    itk::Size<Dimension> size;
    size.Fill(s);
    itk::Index<Dimension> index;
    index.Fill(0);
    
    auto indicesImage = IndicesImageType::New();
    indicesImage->SetRegions(itk::ImageRegion<Dimension>(index,size));
    indicesImage->Allocate();
    
    itk::ImageRegionIterator<IndicesImageType> it(indicesImage, indicesImage->GetBufferedRegion());
    auto v=indices.begin();
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        for(int j=0; j<repeats; ++j){
            it.Value()[j] = *(v++);
//            std::cout << it.Value()[j] << std::endl;
        }

    auto filter = itk::UnaryFunctorImageFilter<IndicesImageType, IndicesImageType, AssignFunctorType<repeats> >::New();
    
    filter->SetInput(indicesImage);
    filter->GetFunctor().vals = &values[0];
    filter->GetFunctor().k=2;
    filter->Update();
    
    for(auto it = values.begin(); it!=values.end(); ++it)
        std::cout << (*it) << std::endl;
    std::cout << "Repeats : " << repeats << ", observed" << *std::min_element(values.begin(), values.end()) << "\n";
    
    return true;
}