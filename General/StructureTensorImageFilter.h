//
//  StructureTensorImageFilter.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 05/03/2014.
//
//

#ifndef itkDiffusion_StructureTensorImageFilter_h
#define itkDiffusion_StructureTensorImageFilter_h

#include "itkCastImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

namespace itk {
    /**
     Implementation of the structure tensor, defined by 
     \f[K_\rho (\nabla u_\sigma \otimes \nabla u_\sigma),\f] where \f$K_\rho\f$ denotes the gaussian kernel of standard deviation \f$\rho\f$, and \f$u_\sigma := K_\sigma * u\f$.
     */
    template<
    typename TImage,
    typename TTensorImage = Image<SymmetricSecondRankTensor<
            typename NumericTraits<typename TImage::PixelType>::RealType ,
            TImage::ImageDimension >,
        TImage::ImageDimension>
    >
    class StructureTensorImageFilter : public ImageToImageFilter<TImage, TTensorImage> {
    public:
        typedef StructureTensorImageFilter Self;
        typedef ImageToImageFilter< TImage, TImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        
        /// Method for creation through the object factory.
        itkNewMacro(Self);
        /// Run-time type information (and related methods).
        itkTypeMacro(StructureTensorImageFilter, Superclass);
        
        typedef TImage ImageType;
        static const unsigned int Dimension = ImageType::ImageDimension;
        typedef TTensorImage TensorImageType;
        typedef typename TensorImageType::PixelType TensorType;
        typedef typename TensorType::ComponentType ScalarType;
        typedef Image<ScalarType, Dimension> ScalarImageType;
        
        ///Parameter \f$\sigma\f$ of the structure tensor definition.
        itkSetMacro(NoiseScale, ScalarType);
        ///Parameter \f$\rho\f$ of the structure tensor definition.
        itkSetMacro(FeatureScale, ScalarType);

        itkGetConstMacro(NoiseScale, ScalarType);
        itkGetConstMacro(FeatureScale, ScalarType);
        
    protected:
        virtual void GenerateData();
        
        ScalarType m_FeatureScale=5, m_NoiseScale=1;
    };
    
    // No hxx file for such trivial composite filter.
    template<typename TI, typename TTI>
    void StructureTensorImageFilter<TI,TTI>::GenerateData(){
        typedef CastImageFilter<ImageType, ScalarImageType> CasterType;
        typename CasterType::Pointer caster = CasterType::New();
        caster->SetInput(this->GetInput());
        
        typedef GradientRecursiveGaussianImageFilter<ScalarImageType> GradientFilterType;
        auto gradientFilter = GradientFilterType::New();
        gradientFilter->SetInput(caster->GetOutput());
        gradientFilter->SetSigma(m_NoiseScale);
        
        typedef typename GradientFilterType::OutputImageType CovariantVectorImageType;
        typedef typename GradientFilterType::CovariantVectorType CovariantVectorType;
        
        struct RankOneTensorFunctor {
            TensorType operator()(const CovariantVectorType & u) const {
                TensorType m;
                for(int i=0; i<Dimension; ++i)
                    for(int j=0; j<Dimension; ++j)
                        m(i,j) = u[i]*u[j];
                return m;
            }
        };
        
        typedef UnaryFunctorImageFilter<CovariantVectorImageType, TensorImageType, RankOneTensorFunctor> ImageFunctorType;
        auto imageFunctor = ImageFunctorType::New();
        imageFunctor->SetInput(gradientFilter->GetOutput());
        
        typedef RecursiveGaussianImageFilter<TensorImageType> GaussianFilterType;
        auto gaussianFilter = GaussianFilterType::New();
        gaussianFilter->SetInput(imageFunctor->GetOutput());
        gaussianFilter->SetSigma(m_FeatureScale);

        gaussianFilter->Update();
        this->GraftOutput(gaussianFilter->GetOutput());        
    }
    
    
}

#endif
