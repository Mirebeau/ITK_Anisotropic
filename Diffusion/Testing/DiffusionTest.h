//
//  DiffusionTest2.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 07/03/2014.
//
//

#ifndef itkDiffusion_DiffusionTest2_h
#define itkDiffusion_DiffusionTest2_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRandomImageSource.h"
#include "itkGaussianDistribution.h"

#include "SyntheticImageSource.h"
#include "CoherenceEnhancingDiffusionFilter.h"


namespace Testing {
    using namespace itk;
    using std::cout;
    using std::endl;
    const std::string imageSourceDirectory = "../Images/";
    
    void EdgeEnhancingDiffusion2D(){ // test on a finger print image
        
        const std::string testName = "EdgeEnhancingDiffusion2D";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        typedef float ValueType;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        // Import Image, construct diffusion tensor
        
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+testName+"_TestImage.png");
        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> AnisotropicDiffusionFilterType;
        AnisotropicDiffusionFilterType::Pointer anisotropicDiffusionFilter = AnisotropicDiffusionFilterType::New();
        anisotropicDiffusionFilter->SetInput(reader->GetOutput());
        
        anisotropicDiffusionFilter->SetEnhancement(AnisotropicDiffusionFilterType::Edge);
        anisotropicDiffusionFilter->SetDiffusionTime(10);
        anisotropicDiffusionFilter->SetExponent(2);
        anisotropicDiffusionFilter->SetNoiseScale(0.5);
        anisotropicDiffusionFilter->SetFeatureScale(4);
        anisotropicDiffusionFilter->SetLambda(4);
        anisotropicDiffusionFilter->SetAlpha(0.001);
        
        auto writer = ImageFileWriter<ImageType>::New();
        writer->SetFileName(testName+".hdf5");
        writer->SetInput(anisotropicDiffusionFilter->GetOutput());
        writer->Update();
    }
    
    
    void CoherenceEnhancingDiffusion2D(){
        
        const std::string testName = "CoherenceEnhancingDiffusion2D";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        typedef float ValueType;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        // Import Image, construct diffusion tensor
        
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+testName+"_TestImage.png");

        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> AnisotropicDiffusionFilterType;
        AnisotropicDiffusionFilterType::Pointer anisotropicDiffusionFilter = AnisotropicDiffusionFilterType::New();
        anisotropicDiffusionFilter->SetInput(reader->GetOutput());
        
        anisotropicDiffusionFilter->SetDiffusionTime(10);
        anisotropicDiffusionFilter->SetExponent(2);
        anisotropicDiffusionFilter->SetNoiseScale(4);
        anisotropicDiffusionFilter->SetFeatureScale(4);
        anisotropicDiffusionFilter->SetLambda(4);
        anisotropicDiffusionFilter->SetAlpha(0.01);
        
        auto writer = ImageFileWriter<ImageType>::New();
        writer->SetFileName(testName+".hdf5");
        writer->SetInput(anisotropicDiffusionFilter->GetOutput());
        writer->Update();
        
        writer->Update();
    }
    
    void CoherenceEnhancingDiffusion3D(){
        const unsigned int Dimension = 3;
        typedef float ValueType;
        
        const std::string testName = "CoherenceEnhancingDiffusion3D";
        cout << "Beginning " << testName << endl;
        
        const long n=100;  // Number of pixels in each dimension
        const ValueType h = 1./n;
        const ValueType Variance = 0.5;
                
        typedef Image<ValueType,Dimension> ImageType;
        
        // Creating the region
        typedef ImageRegion<Dimension> RegionType;
        
        auto randomSource = RandomImageSource<ImageType>::New();
        randomSource->SetMax(1);
        randomSource->SetMin(0);
        randomSource->SetNumberOfThreads(1); // for reproducibility
        Size<Dimension> size;           size.Fill(n);       randomSource->SetSize(size);
        ImageType::SpacingType spacing; spacing.Fill(1./n); randomSource->SetSpacing(spacing);
        // index is automatically set to [0,0,0]
        
        struct FunctorType {
            ImageType * physicalData;
            Statistics::GaussianDistribution * gaussian;
            ValueType operator()(ValueType rand, const ImageType::IndexType & index) const {
                ImageType::PointType x;
                physicalData->TransformIndexToPhysicalPoint(index, x);
                ValueType r = x.GetVectorFromOrigin().GetNorm();
                
                return cos( vnl_math_cube(r/R) ) + gaussian->EvaluateInverseCDF(rand);
            }
            ValueType R=0.5;
        };
        
        auto functorFilter =  UnaryFunctorWithIndexImageFilter<ImageType, ImageType, FunctorType>::New();
        
        auto gaussian = Statistics::GaussianDistribution::New();
        gaussian->SetMean(0.);
        gaussian->SetVariance(Variance);
        functorFilter->GetFunctor().gaussian = gaussian;
        
        randomSource->Update();
        functorFilter->GetFunctor().physicalData = randomSource->GetOutput();
        functorFilter->SetInput(randomSource->GetOutput());
        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> DiffusionFilterType;
        DiffusionFilterType::Pointer diffusionFilter = CoherenceEnhancingDiffusionFilter<ImageType>::New();
        diffusionFilter->SetInput(functorFilter->GetOutput()); //functorFilter->GetOutput()
        
        // Do not forget to take the image scale h into account
        diffusionFilter->SetDiffusionTime(10*h*h);
        diffusionFilter->SetNoiseScale(5*h);
        diffusionFilter->SetFeatureScale(20*h);
        diffusionFilter->SetLambda(1);
        diffusionFilter->SetAlpha(0.01);
        

        auto writer = ImageFileWriter<ImageType>::New();
        writer->SetInput(diffusionFilter->GetOutput());
        writer->SetFileName(testName+".hdf5");
        writer->Update();
        
        
/*        {
            auto writer = ImageFileWriter<ImageType>::New();
            writer->SetInput(functorFilter->GetOutput());
            writer->SetFileName(testName+"_Noisy.hdf5");
            writer->Update();
        }*/
        cout << "finishing " << testName << endl;
    }
    
    
}
#endif
