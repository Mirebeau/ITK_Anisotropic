//
//  AnisotropicFastMarching.h
//  MatlabITKFastMarching
//
//  Created by Jean-Marie Mirebeau on 04/04/13.
//
//

#ifndef MatlabITKFastMarching_AnisotropicFastMarching_h
#define MatlabITKFastMarching_AnisotropicFastMarching_h

#include "itkImageRegionIterator.h"

#include "LinearAnisotropicDiffusionLBRImageFilter.h"
#include "CoherenceEnhancingDiffusionFilter.h"

#include "mexMessageWrapper.h"
#include "mxIO.h"

template<typename NormType>
int MatlabLinearAnisotropicDiffusion(mxIO & IO){
    
    MexMsg() << "Linear Anisotroic diffusion, with prescribed tensor field.\n";
    const int Dimension = NormType::Dimension;
    typedef double ValueType;
    typedef ValueType PixelType;
    
    typedef itk::Image<PixelType,Dimension> ImageType;
    typedef itk::Image<NormType,Dimension> MetricImageType;
    typedef itk::LinearAnisotropicDiffusionLBRImageFilter<ImageType, MetricImageType> DiffType;
    
    
    auto diffusionFilter = DiffType::New();
    diffusionFilter->SetMaxDiffusionTime(IO.GetObject<ValueType>("MaxDiffusionTime",1));
    diffusionFilter->SetMaxNumberOfTimeSteps(IO.GetObject<ValueType>("MaxNumberOfTimeSteps", 10));
    diffusionFilter->SetRatioToMaxStableTimeStep(IO.GetObject<ValueType>("RatioToMaxStableTimeStep",0.5));
    
    IO.transposeFirstTwoCoordinates = IO.GetObject<double>("TransposeFirstTwoImageCoordinates",1.);
    auto tensorImage = IO.GetImage<MetricImageType>("DiffusionTensorField");
    auto image = IO.GetImage<ImageType>("Image");
    
    
    typedef typename MetricImageType::SpacingType SpacingType;
    const SpacingType spacing = IO.GetObject<SpacingType>("Spacing",SpacingType(1.));
    tensorImage->SetSpacing(spacing);
    image->SetSpacing(spacing);
    
    diffusionFilter->SetInputTensor(tensorImage);
    diffusionFilter->SetInputImage(image);
    diffusionFilter->Update();
    
    /*
     MexMsg() << tensorImage->GetBufferedRegion() << "\n";
     MexMsg() << image->GetBufferedRegion() << "\n";
     MexMsg() << diffusionFilter->GetOutput()->GetBufferedRegion() << "\n";
     MexMsg() << "value : " << image->GetBufferPointer()[0] << ", tensor : " << tensorImage->GetBufferPointer()[0] << "\n";
     */
    
    IO.SetField("SmoothedImage", IO.mxImage<ImageType>(diffusionFilter->GetOutput()));
    IO.SetField("EffectiveDiffusionTime", IO.mxObject<double>(diffusionFilter->GetEffectiveDiffusionTime()));
    IO.SetField("EffectiveNumberOfTimeSteps", IO.mxObject<double>(diffusionFilter->GetEffectiveNumberOfTimeSteps()));
    
    return EXIT_SUCCESS;
}


template<typename NormType>
int MatlabCoherenceEnhancingDiffusion(mxIO & IO){
    
    MexMsg() << "Non-linear anisotropic diffusion, with structure tensor based tensor field.\n";
    const int Dimension = NormType::Dimension;
    typedef double ValueType;
    typedef ValueType PixelType;
    
    typedef itk::Image<PixelType,Dimension> ImageType;
    typedef itk::CoherenceEnhancingDiffusionFilter<ImageType> DiffusionFilterType;
    
    typedef typename ImageType::SpacingType SpacingType;
    const ValueType h = IO.GetObject<ValueType>("SpacingScale",1);
    const SpacingType spacing = IO.GetObject<SpacingType>("Spacing",SpacingType(h));
    
    
    auto diffusionFilter = DiffusionFilterType::New();
    diffusionFilter->SetDiffusionTime(IO.GetObject<ValueType>("DiffusionTime",1*h*h));
    diffusionFilter->SetMaxTimeStepsBetweenTensorUpdates(IO.GetObject<ValueType>("MaxTimeStepsBetweenTensorUpdates",10));
    diffusionFilter->SetRatioToMaxStableTimeStep(IO.GetObject<ValueType>("RatioToMaxStableTimeStep",0.5));
    diffusionFilter->SetNoiseScale(IO.GetObject<ValueType>("NoiseScale",1*h));
    diffusionFilter->SetFeatureScale(IO.GetObject<ValueType>("FeatureScale",5*h));
    
    if(IO.GetObject<ValueType>("EdgeEnhancement",0))
            diffusionFilter->SetEnhancement(DiffusionFilterType::Edge);
    else    diffusionFilter->SetEnhancement(DiffusionFilterType::Coherence);
    diffusionFilter->SetExponent(IO.GetObject<ValueType>("Exponent",2));
    diffusionFilter->SetLambda(IO.GetObject<ValueType>("Lambda", 4));
    diffusionFilter->SetAlpha(IO.GetObject<ValueType>("Alpha",0.01));
    
    IO.transposeFirstTwoCoordinates = IO.GetObject<double>("TransposeFirstTwoImageCoordinates",1.);
    auto image = IO.GetImage<ImageType>("Image");
    image->SetSpacing(spacing);
    
//    std::cout << image->GetBufferedRegion() << "\n" << image->GetRequestedRegion() << "\n";
//    return EXIT_SUCCESS;
    
    diffusionFilter->SetInput(image);
    diffusionFilter->Update();
        
    IO.SetField("SmoothedImage", IO.mxImage<ImageType>(diffusionFilter->GetOutput()));
    
    std::vector<std::pair<ValueType, ValueType> > linearFilterEffectiveTimesAndIterations;
    for(auto ti : diffusionFilter->GetLinearFilterEffectiveTimesAndIterations())
        linearFilterEffectiveTimesAndIterations.push_back({ti.first, ValueType(ti.second)});
    IO.SetField("LinearFilterEffectiveTimesAndIterations", IO.mxVector(linearFilterEffectiveTimesAndIterations));
    return EXIT_SUCCESS;
}


template<typename NormType>
int AnisotropicDiffusion(mxIO & IO){
    
    if(IO.HasField("DiffusionTensorField")) return MatlabLinearAnisotropicDiffusion<NormType>(IO);
    else                                    return MatlabCoherenceEnhancingDiffusion<NormType>(IO);
}





#endif
