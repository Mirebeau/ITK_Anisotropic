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
#include "RescaleAndExportBMP.h"
namespace DiffusionBasicTest {
    using namespace itk;
    using std::cout;
    using std::endl;
    
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
            
            ImageSourceType::PhysicalDataType::Pointer physicalData = imageSource->GetPhysicalData();
            physicalData->SetRegions(ImageType::RegionType(index, size));
            physicalData->SetSpacing(spacing);
        }
        
        {
            typedef RescaleAndExportBMP<ImageType> ExporterType;
            ExporterType::Pointer exporter = ExporterType::New();
            exporter->SetInput(imageSource->GetOutput());
            exporter->SetFileName("TestImage.bmp");
            exporter->Update();
        }
        
        typedef StructureTensorImageFilter<ImageType> StructureTensorImageFilterType;
        StructureTensorImageFilterType::Pointer structureTensorFilter =  StructureTensorImageFilterType::New();
        structureTensorFilter->SetInput(imageSource->GetOutput());
        structureTensorFilter->Update();
        
        typedef StructureTensorImageFilterType::TensorType TensorPixelType;
        TensorPixelType S = *structureTensorFilter->GetOutput()->GetBufferPointer();
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
        
        typedef LinearAnisotropicDiffusionLBRImageFilter<ImageType, ValueType> DiffusionOperatorType;
        DiffusionOperatorType::Pointer DiffusionFilter;
        
        DiffusionFilter = DiffusionOperatorType::New();
        
        // Choose region
        typedef ImageRegion<Dimension> RegionType;
        RegionType BufferedRegion;
        {
            const long n=100;
            Size<Dimension> size; size.Fill(n);
            Index<Dimension> index; index.Fill(-n/2);
            BufferedRegion = RegionType(index,size);
        }
        
        // Setup the metric
        
        TensorImageType::Pointer TensorImage = TensorImageType::New();
        {
            TensorImage->SetRegions(BufferedRegion);
            TensorImage->Allocate();
            
            TensorType tensor;
//            tensor.SetIdentity();
            tensor(0,0)=0.5;
            tensor(0,1)=-0.8;
            tensor(1,1)=2.;
            TensorImage->FillBuffer(tensor);
            
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
//            DiffusionFilter->SetMaxDiffusionTime(1);
//            DiffusionFilter->SetMaxNumberOfTimeSteps(5);
            DiffusionFilter->SetMaxDiffusionTime(100);
            DiffusionFilter->SetMaxNumberOfTimeSteps(500);
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
    
    void NumTraitsTest(){
        
        using namespace itk;
        
        typedef Vector<int,3> PixelType;
        typedef NumericTraits<PixelType>::RealType RealType;
        RealType a;
        NumericTraits<PixelType>::ScalarRealType b;
        std::cout << a << ", " << b << std::endl;
        
        
        typedef Image<PixelType,2> ImageType;
        typedef Image<RealType,2> RealImageType;
        typedef CastImageFilter<ImageType,RealImageType> CasterType;
        CasterType::Pointer caster = CasterType::New();
    }

    void RGBConvertTest(){
        {
            // Usage with uchar pixel type, cast to double.
            typedef Image<unsigned char,3> ImageType;
            typedef Image<double,3> FPImageType;
            typedef CastImageFilter<ImageType, FPImageType> CastFilterType;
            typedef CoherenceEnhancingDiffusionFilter<FPImageType> CEDFilterType;
            
            CastFilterType::Pointer castFilter = CastFilterType::New();
            CEDFilterType::Pointer cedFilter = CEDFilterType::New();
            cedFilter->SetInput(castFilter->GetOutput());
        }

        {
            // Usage with RGB pixel type, cast to Vector
            typedef Vector<float,3> VectorType;
            typedef RGBPixel<unsigned char> RGBType;
            struct CastFunctorType {
                VectorType operator()(const RGBType & x){
                    VectorType result; for(int i=0; i<3; ++i) result[i]=x[i]; return result;}
            };
            typedef Image<RGBType, 2> ImageType;
            typedef Image<VectorType, 2> FPImageType;
            typedef UnaryFunctorImageFilter<ImageType, FPImageType, CastFunctorType> CastFilterType;
            typedef CoherenceEnhancingDiffusionFilter<FPImageType, float> CEDFilterType;
            
            CastFilterType::Pointer castFilter = CastFilterType::New();
            CEDFilterType::Pointer cedFilter = CEDFilterType::New();
            cedFilter->SetInput(castFilter->GetOutput());
        }

    
    }
    
}
#endif
