//
//  DiffusionBasicTest.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 10/03/2014.
//
//

#ifndef itkDiffusion_DiffusionBasicTest_h
#define itkDiffusion_DiffusionBasicTest_h

#include "AnisotropicDiffusionLBRImageFilter.h"

namespace DiffusionBasicTest {
    using namespace itk;

    
    void StructureTensor(int resolution = 10){
        typedef float ScalarType;
        static const int Dimension = 2;
        typedef Image<ScalarType,Dimension> ImageType;
        typedef Point<double,Dimension> PointType;
        
        struct FunctorType {
            ScalarType operator()(const PointType & x){
                return x[0]+x[1];//sin(x[0]+0.5*x[1]);
            }
        } functor;
        
        typedef SyntheticImageSource<ImageType, FunctorType> ImageSourceType;
        ImageSourceType::Pointer imageSource = ImageSourceType::New();
        imageSource->SetFunctor(functor);
        
        {
            ImageType::IndexType index;
            index.Fill(0);
            ImageType::SizeType size;
            size.Fill(resolution);
            ImageType::SpacingType spacing;
            spacing.Fill(1./resolution);
            
            auto physicalData = imageSource->GetPhysicalData();
            physicalData->SetRegions(ImageType::RegionType(index, size));
            physicalData->SetSpacing(spacing);
        }
        
        {
            auto exporter = RescaleAndExportBMP<ImageType>::New();
            exporter->SetInput(imageSource->GetOutput());
            exporter->SetFileName("TestImage.bmp");
            exporter->Update();
        }
        
        auto structureTensorFilter =  StructureTensorImageFilter<ImageType>::New();
        structureTensorFilter->SetInput(imageSource->GetOutput());
        structureTensorFilter->Update();
        
        auto S = *structureTensorFilter->GetOutput()->GetBufferPointer();
        std::cout << S << std::endl;
    }
    
    
    
    void LinearDiffusion(){
        const std::string testName = "LinearDiffusionTest";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        typedef double ValueType;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        typedef SymmetricSecondRankTensor<ValueType,Dimension> TensorType;
        typedef Image<TensorType,Dimension> TensorImageType;
        
        typedef LinearAnisotropicDiffusionLBRImageFilter<ImageType, TensorImageType> DiffusionOperatorType;
        DiffusionOperatorType::Pointer DiffusionFilter;
        
        DiffusionFilter = DiffusionOperatorType::New();
        
        // Choose region
        typedef ImageRegion<Dimension> RegionType;
        RegionType BufferedRegion;
        {
            const long n=10;
            Size<Dimension> size; size.Fill(n);
            Index<Dimension> index; index.Fill(-n/2);
            BufferedRegion = RegionType(index,size);
        }
        
        // Setup the metric
        
        TensorImageType::Pointer TensorImage = TensorImageType::New();
        {
            TensorImage->SetRegions(BufferedRegion);
            TensorImage->Allocate();
            
            TensorType identity;
            identity.SetIdentity();
            TensorImage->FillBuffer(identity);
            
            TensorImageType::SpacingType spacing;
            spacing.Fill(1);
            
            TensorImage->SetSpacing(spacing);
        }
        
        // Setup the image
        ImageType::Pointer InputImage;
        {
            InputImage=ImageType::New();
            InputImage->SetRegions(BufferedRegion);
            InputImage->Allocate();
            
            typedef Index<Dimension> IndexType;
            IndexType centerIndex;
            centerIndex.Fill(0);
            InputImage->GetPixel(centerIndex) = 1;
        }
        
        
        // Setup the diffusion filter
        {
            DiffusionFilter->SetInputTensor(TensorImage);
            
            DiffusionFilter->SetInputImage(InputImage);
            DiffusionFilter->SetMaxDiffusionTime(1);
            DiffusionFilter->SetMaxNumberOfTimeSteps(5);
            DiffusionFilter->SetRatioToMaxStableTimeStep(0.5);
        }
        
        // Rescale and export
        {
            typedef RescaleAndExportBMP<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput(DiffusionFilter->GetOutput());
            
            //            writer->SetInput(InputImage);
            writer->SetFileName(testName+".bmp");
            writer->Update();
            
        }
        cout << "Finishing " << testName << endl;
    }

}
#endif
