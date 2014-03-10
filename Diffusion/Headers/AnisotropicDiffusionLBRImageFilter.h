//
//  AnisotropicDiffusionLBRImageFilter.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 28/02/2014.
//
//

#ifndef itkDiffusion_AnisotropicDiffusionLBRImageFilter_h
#define itkDiffusion_AnisotropicDiffusionLBRImageFilter_h

#include "LinearAnisotropicDiffusionLBRImageFilter.h"
#include "StructureTensorImageFilter.h"

namespace itk
{
    /** Implementation of Non-linear Anisotropic Diffusion.
     This class repeatedly calls the LinearAnisotropicDiffusionLBRImageFilter, with non-linear diffusion tensors built on the fly. These tensors are obtained by computing the image structure tensors, and appropriately modifying their eigenvalues with the method EigenValuesTransform. The latter method is not implemented, and needs to be provided in a subclass, such as CoherenceEnhancingDiffusionFilter.h
    */
    template<typename TImage>
    class AnisotropicDiffusionLBRImageFilter : public ImageToImageFilter< TImage, TImage>
    {
    public:
        typedef  AnisotropicDiffusionLBRImageFilter Self;
        typedef ImageToImageFilter< TImage, TImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        
        /// Method for creation through the object factory.
        itkNewMacro(Self);
        /// Run-time type information (and related methods).
        itkTypeMacro(AnisotropicDiffusionLBRImageFilter, Superclass);
        
        typedef TImage ImageType;
        static const unsigned int Dimension = ImageType::ImageDimension;
        
        // The float type in which computations are done.
        typedef typename NumericTraits<typename ImageType::PixelType>::RealType ScalarType;
        typedef Image<ScalarType, Dimension> ScalarImageType;
        
        typedef SymmetricSecondRankTensor<ScalarType,Dimension> TensorType;
        typedef Image<TensorType,Dimension> TensorImageType;
        
        typedef StructureTensorImageFilter<ScalarImageType, TensorImageType> StructureTensorFilterType;
        typedef LinearAnisotropicDiffusionLBRImageFilter<ScalarImageType, TensorImageType> LinearDiffusionFilterType;
        
        virtual void SetNoiseScale(ScalarType scale){   structureTensorFilter->SetNoiseScale(scale);}
        virtual void SetFeatureScale(ScalarType scale){ structureTensorFilter->SetFeatureScale(scale);}
        virtual void SetMaxTimeStepsBetweenTensorUpdates(int n){       linearDiffusionFilter->SetMaxNumberOfTimeSteps(n);}
        virtual void SetRatioToMaxStableTimeStep(ScalarType ratio){ linearDiffusionFilter->SetRatioToMaxStableTimeStep(ratio);}

        /*
        /// Access to SetNoiseScale(ScalarType), SetFeatureScale(ScalarType) of the structure tensor
        itkGetMacro(structureTensorFilter, StructureTensorFilterType*);
        /// Access to SetRatioToMaxStableTimeStep(ValueType ratio), SetMaxNumberOfTimeSteps(int n);
        itkGetMacro(linearDiffusionFilter, LinearDiffusionFilterType*);
        */
        
        itkSetMacro(DiffusionTime, ScalarType);
        itkGetConstMacro(DiffusionTime, ScalarType);
        
        typedef typename TensorType::EigenValuesArrayType EigenValuesArrayType;
        /// Should be specialized. Transformation of the Structure tensor eigenvalues into the diffusion tensor eigenvalues. (Structure tensor eigenvalues are sorted increasingly for convenience).
        virtual EigenValuesArrayType EigenValuesTransform(const EigenValuesArrayType &) const
        {itkExceptionMacro("Undefined tensor eigenvalues");};
        
        virtual typename TensorImageType::Pointer GetLastTensorImage(){return tensorImage;}
        typedef std::vector< std::pair<ScalarType, int> > EffectiveTimesAndIterationsType;
        itkGetConstReferenceMacro(LinearFilterEffectiveTimesAndIterations, EffectiveTimesAndIterationsType);
    protected:
        typename StructureTensorFilterType::Pointer structureTensorFilter;
        typename LinearDiffusionFilterType::Pointer linearDiffusionFilter;
        
        AnisotropicDiffusionLBRImageFilter(){
            structureTensorFilter = StructureTensorFilterType::New();
            linearDiffusionFilter = LinearDiffusionFilterType::New();
        };
        ~AnisotropicDiffusionLBRImageFilter(){};
        
        typename ScalarImageType::Pointer image;
        typename TensorImageType::Pointer tensorImage;
        virtual void ComputeDiffusionTensors();

        
        ScalarType m_DiffusionTime=0;
        virtual void GenerateData();
        
        EffectiveTimesAndIterationsType m_LinearFilterEffectiveTimesAndIterations;
    };
}

#include "AnisotropicDiffusionLBRImageFilter.hxx"
#endif
