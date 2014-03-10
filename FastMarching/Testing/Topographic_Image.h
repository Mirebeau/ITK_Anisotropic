//
//  Topographic_Image.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 18/09/13.
//
//

#ifndef ITKFM_Topographic_Image_h
#define ITKFM_Topographic_Image_h

#include "Topographic_Test.h"
#include "itkBinomialBlurImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"

namespace AnisotropicFastMarching_Topographic_Test {
    void GenerateImage(){
        auto imageSource = itk::CurveNeighborhoodSource<ScalarImageType>::New();
        imageSource->SetOutputSpacing(path.spacing);
        imageSource->SetOutputOrigin(path.origin);
        
        std::vector<PointType> physicalPath;
        for(int i=0; i<=path.sampling; ++i)
            physicalPath.push_back(path(i/path.sampling));
        
        imageSource->AddPhysicalPath(physicalPath);
        imageSource->SetOutputRegionContainingPaths();
        imageSource->Update();
        
        auto curvilinearCoordinates = imageSource->ComputeCurvilinearCoordinates();
        
        auto calculator = itk::MinimumMaximumImageCalculator<ScalarImageType>::New();
        calculator->SetImage(curvilinearCoordinates);
        calculator->SetRegion(curvilinearCoordinates->GetBufferedRegion());
        calculator->ComputeMaximum();
        
        struct HeightOnPathFunctorType {
            ScalarType smax;
            ScalarType operator()(ScalarType s) const {
                return 1-fabsf(2*s/smax-1);
            }
        };
        
        auto heightOnPathFilter = itk::UnaryFunctorImageFilter<ScalarImageType, ScalarImageType, HeightOnPathFunctorType>::New();
        heightOnPathFilter->GetFunctor().smax = calculator->GetMaximum();
        heightOnPathFilter->SetInput(curvilinearCoordinates);
        
//        auto heightOnPath = heightOnPathFilter->GetOutput();
        auto pathNeighborhood = imageSource->TubeSmoothIndicatorFunction(path.width, path.width);
        
        auto blurFilter = itk::BinomialBlurImageFilter<ScalarImageType, ScalarImageType>::New();
        blurFilter->SetInput(heightOnPathFilter->GetOutput());
        
        struct BarycenterFunctorType {
            ScalarType operator()(ScalarType a, ScalarType b, ScalarType t){
                return (1-t)*a+t*b;
            }
        };
        
        auto barycenterFilter = itk::TernaryFunctorImageFilter<ScalarImageType, ScalarImageType, ScalarImageType, ScalarImageType, BarycenterFunctorType>::New();
        barycenterFilter->SetInput1(blurFilter->GetOutput());
        barycenterFilter->SetInput2(heightOnPathFilter->GetOutput());
        barycenterFilter->SetInput3(pathNeighborhood);
        
        barycenterFilter->Update(); 
        
        blurFilter->SetInput(barycenterFilter->GetOutput());
        for(int i=0; i<200; ++i)
            barycenterFilter->Update(); // Hopefully no infinite loop here.
        
        struct MultFunctorType {
            ScalarType alpha;
            ScalarType operator()(ScalarType s) const {return alpha*s;}
        };
        auto multFilter = itk::UnaryFunctorImageFilter<ScalarImageType, ScalarImageType, MultFunctorType>::New();
        multFilter->GetFunctor().alpha = 30;
        multFilter->SetInput(barycenterFilter->GetOutput());
        
        {
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(multFilter->GetOutput());
            writer->SetFileName(testPrefix+testImageFilename);
            writer->Update();
        }
        {
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(pathNeighborhood);
            writer->SetFileName(testPrefix+"TubeNeighborhood.hdf5");
            writer->Update();
        }
    }
}

#endif
