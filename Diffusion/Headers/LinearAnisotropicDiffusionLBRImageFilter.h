//
//  LinearAnisotropicDiffusionLBRImageFilter.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 28/02/2014.
//
//

#ifndef itkDiffusion_LinearAnisotropicDiffusionLBRImageFilter_h
#define itkDiffusion_LinearAnisotropicDiffusionLBRImageFilter_h


#include "itkImageToImageFilter.h"
#include "Riemannian2DNorm.h"
#include "Riemannian3DNorm.h"

namespace itk
{
    /**
     Implementation of Anisotropic Diffusion
     \f[\partial_t u = {\rm div} (D \nabla u),\f]
     with Neumann boundary conditions. The numerical scheme is stable and satisfies the maximum principle, even for strongly anisotropic tensors, thanks to an adaptive discretization using arithmetic techniques (Lattice Basis Reduction, LBR).
     */
    template< typename TImage, typename TTensorImage>
    class LinearAnisotropicDiffusionLBRImageFilter : public ImageToImageFilter< TImage, TImage >
    {
    public:
        /** Standard class typedefs. */
        typedef LinearAnisotropicDiffusionLBRImageFilter    Self;
        typedef ImageToImageFilter< TImage, TImage >        Superclass;
        typedef SmartPointer< Self >                        Pointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(LinearAnisotropicDiffusionLBRImageFilter, ImageToImageFilter);
        
        /** The image to be inpainted in regions where the mask is white.*/
        void SetInputImage(const TImage* image);
        
        /** The mask to be inpainted. White pixels will be inpainted, black pixels will be passed through to the output.*/
        void SetInputTensor(const TTensorImage* tensorImage);
        
        typedef TImage ImageType;
        typedef TTensorImage TensorImageType;
        static const int Dimension = TensorImageType::ImageDimension;
        typedef typename TensorImageType::PixelType TensorType; // Symmetric second rank tensor
        typedef typename TensorType::ComponentType ValueType;
        typedef ImageRegion<Dimension> RegionType;
        
        void SetMaxDiffusionTime(ValueType time);
        void SetMaxNumberOfTimeSteps(int n);
        void SetRatioToMaxStableTimeStep(ValueType ratio);
        
        itkGetConstMacro(EffectiveDiffusionTime, ValueType);
        itkGetConstMacro(EffectiveNumberOfTimeSteps, int);
        
        itkGetConstMacro(SparseMatrixAssemblyTimeCost, ValueType);
        itkGetConstMacro(IterationsTimeCost, ValueType);
    protected:
        LinearAnisotropicDiffusionLBRImageFilter();
        ~LinearAnisotropicDiffusionLBRImageFilter(){}
        
        typename TImage::ConstPointer GetInputImage();
        typename TTensorImage::ConstPointer GetInputTensor();
        
        typedef Index<Dimension> IndexType;

        // ******* Containers for the stencils used in the discretization
        static const unsigned int HalfStencilSize = (Dimension == 2) ? 3 : 6;
        static const unsigned int StencilSize = 2*HalfStencilSize;
        
        typedef Vector<ValueType,HalfStencilSize> StencilCoefficientsType;
        typedef Offset<Dimension> OffsetType;
        typedef Vector<OffsetType, HalfStencilSize> StencilOffsetsType;
        
        typedef int InternalSizeT;
        typedef Vector<InternalSizeT,StencilSize> StencilBufferIndicesType;
        
        // ****************** Stencil support ******************
        
        template<unsigned int NDimension=Dimension, typename Dummy=void> struct SSRT_Traits;
        
        template<typename Dummy> struct SSRT_Traits<2,Dummy>{
            typedef Riemannian2DNorm<ValueType,int> SSRT;
            static void GetDiffusionStencil(const SSRT &D, StencilOffsetsType &Offsets, StencilCoefficientsType &Coefficients);
        };
        
        template<typename Dummy> struct SSRT_Traits<3,Dummy>{
            typedef Riemannian3DNorm<ValueType,int> SSRT;
            static void GetDiffusionStencil(const SSRT &D, StencilOffsetsType &Offsets, StencilCoefficientsType &Coefficients);
        };
                
        // *************** Computation *****************
        virtual void GenerateData();
        // These methods are called by generate data
        virtual void GenerateStencils();
        virtual void ImageUpdateLoop();
        
        typedef std::pair<StencilBufferIndicesType, StencilCoefficientsType> StencilType;
        typedef Image<StencilType,Dimension> StencilImageType;
        typename StencilImageType::Pointer stencilImage;
        
        typedef Image<ValueType,Dimension> ScalarImageType;
        typename ScalarImageType::Pointer diagonalCoefficients;
        
        virtual ValueType MaxStableTimeStep();
        
        ValueType diffusionTime = 1;
        ValueType ratioToMaxStableTimeStep = 1;
        int maximumNumberOfTimeSteps = 10;
        
        ValueType m_EffectiveDiffusionTime=0;
        int m_EffectiveNumberOfTimeSteps=0;

        virtual void ImageUpdate(ValueType delta);
        typename ScalarImageType::Pointer previousImage, nextImage;
        
        virtual RegionType GetRequestedRegion(){return GetInputImage()->GetRequestedRegion();}
        
        InternalSizeT OutsideBufferIndex() const {return NumericTraits<InternalSizeT>::max();}
        
        ValueType m_SparseMatrixAssemblyTimeCost, m_IterationsTimeCost;
        
    private:
        LinearAnisotropicDiffusionLBRImageFilter(const Self &); //purposely not implemented
        void operator=(const Self &);  //purposely not implemented
        
    };
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "LinearAnisotropicDiffusionLBRImageFilter.hxx"
#endif

#endif
