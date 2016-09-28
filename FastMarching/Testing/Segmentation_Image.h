//
//  Segmentation_Image.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 17/09/13.
//
//

#ifndef ITKFM_Segmentation_Image_h
#define ITKFM_Segmentation_Image_h

#include "Segmentation_Test.h"
#include "itkFileTools.h"

namespace AnisotropicFastMarching_Segmentation_Test {

    void GenerateImage(){
        // Creating image source
        auto imageSource = itk::CurveNeighborhoodSource<ScalarImageType>::New();
        imageSource->SetOutputSpacing(path.spacing);
        imageSource->SetOutputOrigin(path.origin);
        
        std::vector<PointType> physicalPath;
        for(int i=0; i<=path.sampling; ++i)
            physicalPath.push_back(path(i/path.sampling));
        
        imageSource->AddPhysicalPath(physicalPath);
        imageSource->SetOutputRegionContainingPaths();
        imageSource->Update();
        
        struct MultFunctorType {
            ScalarType alpha;
            ScalarType operator()(ScalarType s){return alpha*s;}
        };
        auto multFilter = itk::UnaryFunctorImageFilter<ScalarImageType, ScalarImageType, MultFunctorType>::New();
        multFilter->GetFunctor().alpha = itk::Math::pi;
        multFilter->SetInput(imageSource->InteriorSmoothIndicatorFunction(path.width));
                
        auto writer = itk::ImageFileWriter<ScalarImageType>::New();
        writer->SetFileName(testPrefix+testImageFilename);
        writer->SetInput(multFilter->GetOutput());
        
        itk::FileTools::CreateDirectory("Hello");
        
        writer->Update();
    }
    
} // end of namespace

#endif
