//
//  TubeCenterline2D_Image.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 19/09/13.
//
//

#ifndef ITKFM_TubeCenterline2D_Image_h
#define ITKFM_TubeCenterline2D_Image_h

#include "TubeCenterline2D_Test.h"

namespace AnisotropicFastMarching_TubeCenterline2D_Test {
    
    void GenerateImage(){
        itk::ImageRegion<Dimension> region;
        {// SingleScale
            auto imageSource = itk::CurveNeighborhoodSource<ScalarImageType>::New();
            imageSource->SetOutputOrigin(path.origin);
            imageSource->SetOutputSpacing(path.spacing);
            
            std::vector<PointType> pathPoints;
            for(int i=0; i<path.nPoints; ++i)
                pathPoints.push_back(path.GetPoint(i));
            imageSource->AddPhysicalPath(pathPoints);
            
            imageSource->SetOutputRegionContainingPaths();
            imageSource->Update();
            region = imageSource->GetOutput()->GetLargestPossibleRegion();
            
            struct FunctorType {
                ScalarType operator()(ScalarType s) const {return 1-s;}
            };
            const ScalarType radius = path.ConstantRadius;
            auto functorFilter = itk::UnaryFunctorImageFilter<ScalarImageType, ScalarImageType, FunctorType>::New();
            functorFilter->SetInput(imageSource->TubeSmoothIndicatorFunction(radius, radius/2));
            
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(functorFilter->GetOutput());
            writer->SetFileName(testPrefix+"TestImageConstantRadius.hdf5");
            writer->Update();
        }
       
        {// Multiscale
            
            // First half of path
            auto imageSource0 = itk::CurveNeighborhoodSource<ScalarImageType>::New();
            auto imageSource1 = itk::CurveNeighborhoodSource<ScalarImageType>::New();
            
            struct MaxFunctor {
                ScalarType operator()(ScalarType s, ScalarType t){return 1-std::max(s,t);}
            };
            auto maxFilter = itk::BinaryFunctorImageFilter<ScalarImageType, ScalarImageType, ScalarImageType, MaxFunctor>::New();
            
            for(int k=0; k<2; ++k){
                auto imageSource = (k==0 ? imageSource0 : imageSource1);
                imageSource->SetOutputOrigin(path.origin);
                imageSource->SetOutputSpacing(path.spacing);
                
                std::vector<PointType> pathPoints;
                std::vector<ScalarType> radiiWidths, boundaryWidths;
                for(int i=0; i<path.nPoints; ++i)
                    if( (i-path.nPoints/2)*(2*k-1) >=0 ){
                        pathPoints.push_back(path.GetPoint(i));
                        radiiWidths.push_back(path.Radius(i));
                        boundaryWidths.push_back(path.Radius(i)/2);
                    }
                imageSource->AddPhysicalPath(pathPoints);
                imageSource->SetOutputRegion(region);
                imageSource->Update();
                if(k==0) maxFilter->SetInput1(imageSource->TubeSmoothIndicatorFunction(radiiWidths, boundaryWidths));
                else maxFilter->SetInput2(imageSource->TubeSmoothIndicatorFunction(radiiWidths, boundaryWidths));
            }
            
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(maxFilter->GetOutput());
            writer->SetFileName(testPrefix+"TestImage.hdf5");
            writer->Update();
        }
        
    }
    
    
    const ScalarType PathType::points[nPoints][2] = {{-2.,0.},{-1.77415,0.117637},{-1.57564,0.230699},{-1.40298,0.339408},{-1.2547,0.443989},{-1.12933,0.544667},{-1.02541,0.641664},{-0.941451,0.735205},{-0.875989,0.825515},{-0.827552,0.912816},{-0.794667,0.997333},{-0.775861,1.07929},{-0.769664,1.15891},{-0.774603,1.23642},{-0.789205,1.31204},{-0.812,1.386},{-0.841515,1.45852},{-0.876277,1.52982},{-0.914816,1.60013},{-0.955659,1.66967},{-0.997333,1.73867},{-1.03837,1.80734},{-1.07729,1.87593},{-1.11263,1.94463},{-1.14291,2.0137},{-1.16667,2.08333},{-1.18272,2.15369},{-1.19109,2.22455},{-1.19211,2.29564},{-1.18608,2.36666},{-1.17333,2.43733},{-1.15419,2.50737},{-1.12896,2.57647},{-1.09797,2.64436},{-1.06155,2.71074},{-1.02,2.77533},{-0.973653,2.83785},{-0.922827,2.89799},{-0.86784,2.95548},{-0.809013,3.01002},{-0.746667,3.06133},{-0.68112,3.10913},{-0.612693,3.15311},{-0.541707,3.193},{-0.46848,3.2285},{-0.393333,3.25933},{-0.316587,3.28521},{-0.23856,3.30583},{-0.159573,3.32092},{-0.0799467,3.33018},{0.,3.33333},{0.0799467,3.33018},{0.159573,3.32092},{0.23856,3.30583},{0.316587,3.28521},{0.393333,3.25933},{0.46848,3.2285},{0.541707,3.193},{0.612693,3.15311},{0.68112,3.10913},{0.746667,3.06133},{0.809013,3.01002},{0.86784,2.95548},{0.922827,2.89799},{0.973653,2.83785},{1.02,2.77533},{1.06155,2.71074},{1.09797,2.64436},{1.12896,2.57647},{1.15419,2.50737},{1.17333,2.43733},{1.18608,2.36666},{1.19211,2.29564},{1.19109,2.22455},{1.18272,2.15369},{1.16667,2.08333},{1.14291,2.0137},{1.11263,1.94463},{1.07729,1.87593},{1.03837,1.80734},{0.997333,1.73867},{0.955659,1.66967},{0.914816,1.60013},{0.876277,1.52982},{0.841515,1.45852},{0.812,1.386},{0.789205,1.31204},{0.774603,1.23642},{0.769664,1.15891},{0.775861,1.07929},{0.794667,0.997333},{0.827552,0.912816},{0.875989,0.825515},{0.941451,0.735205},{1.02541,0.641664},{1.12933,0.544667},{1.2547,0.443989},{1.40298,0.339408},{1.57564,0.230699},{1.77415,0.117637},{2.,0.}};
        
    
    
    
    // Old
    
    void GenerateImage_Tree(){
        typedef unsigned char RGBComponentType;
        typedef itk::RGBPixel<RGBComponentType> RGBPixelType;
        typedef itk::Image<RGBPixelType,Dimension> RGBImageType;
        auto reader = itk::ImageFileReader<RGBImageType>::New();
        reader->SetFileName(testPrefix+"TreeTest.png");
        std::cout << reader->GetFileName() << std::endl;
        
        reader->UpdateOutputInformation();
        {
            const auto io = reader->GetImageIO();
            assert(io->GetNumberOfDimensions() == Dimension);
            assert(io->GetNumberOfComponents() == 3);
            assert(io->GetComponentType() == itk::ImageIOBase::UCHAR);
        }
        
        struct RGBToBWFunctor {
            ScalarType operator()(const RGBPixelType & rgb) const {
                return (ScalarType(rgb[0])+ScalarType(rgb[0])+ScalarType(rgb[0]))/(3*itk::NumericTraits<RGBComponentType>::max());
            }
        };
        
        auto functorFilter = itk::UnaryFunctorImageFilter<RGBImageType, ScalarImageType, RGBToBWFunctor>::New();
        functorFilter->SetInput(reader->GetOutput());
        
        auto writer = itk::ImageFileWriter<ScalarImageType>::New();
        writer->SetFileName(testPrefix+"TreeTest.hdf5");
        writer->SetInput(functorFilter->GetOutput());
        
        writer->Update();
    }
    
}

#endif
