//
//  LinearAnisotropicDiffusionLBRImageFilter.hxx
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 28/02/2014.
//
//

#ifndef itkDiffusion_LinearAnisotropicDiffusionLBRImageFilter_hxx
#define itkDiffusion_LinearAnisotropicDiffusionLBRImageFilter_hxx

#include "itkUnaryFunctorImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "sys/time.h"


#include "UnaryFunctorWithIndexImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"

namespace itk
{
    
    template< typename TImage, typename TTensorImage>
    LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::LinearAnisotropicDiffusionLBRImageFilter()
    {
        this->SetNumberOfRequiredInputs(2);
    }
    
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::SetInputImage(const TImage* image)
    {
        this->SetNthInput(0, const_cast<TImage*>(image));
    }
    
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::SetInputTensor(const TTensorImage* tensorImage)
    {
        this->SetNthInput(1, const_cast<TTensorImage*>(tensorImage));
    }
    
    template< typename TImage, typename TTensorImage>
    typename TImage::ConstPointer LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::GetInputImage()
    {
        return static_cast< const TImage * >
        ( this->ProcessObject::GetInput(0) );
    }
    
    template< typename TImage, typename TTensorImage>
    typename TTensorImage::ConstPointer LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::GetInputTensor()
    {
        return static_cast< const TTensorImage * >
        ( this->ProcessObject::GetInput(1) );
    }
    
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::GenerateData()
    {
        timeval start, end;
        gettimeofday(&start, nullptr);
        
        GenerateStencils();
        
        gettimeofday(&end, nullptr);
        m_SparseMatrixAssemblyTimeCost = (end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/1000000.;
        start=end;
        
        ImageUpdateLoop();
        
        gettimeofday(&end, nullptr);
        m_IterationsTimeCost = (end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/1000000.;
    }
    
// **************************** Computation ***********************
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::GenerateStencils(){
                
        // Stencil type is a pair type because itk::UnaryFunctorImage filter
        // only produces one output
        typedef typename TensorImageType::SpacingType SpacingType;
        const RegionType region = GetRequestedRegion();


        struct StencilFunctor {
            void Initialize(RegionType region_, SpacingType spacing){
                region = region_;
                prod[0]=1;
                for(int i=1; i<Dimension; ++i)
                    prod[i]=prod[i-1]*region.GetSize()[i-1];
                for(int i=0; i<Dimension; ++i)
                   invSpacing[i] = ValueType(1)/spacing[i];
            }

            InternalSizeT BufferIndex(const IndexType & x) const {
                IndexValueType ans=0;
                for(int i=0; i<Dimension; ++i)
                    ans+=this->prod[i]*(x[i]-this->region.GetIndex()[i]);
                return ans;
            }
            
            StencilType operator()(const TensorType & tensor, const IndexType & x) const {
                StencilType stencil;
                StencilOffsetsType offsets;
                
                // The constructor SSRT(S,s) was meant for a norm S, and a spacing s.
                // Diffusion tensors are homogeneous to the inverse of norms, and are thus rescaled with an inverse spacing.
                typename SSRT_Traits<>::SSRT D(tensor, this->invSpacing);
                SSRT_Traits<>::GetDiffusionStencil(D, offsets, stencil.second);
                
                InternalSizeT * yIndex = &stencil.first[0];
                
                //Compute buffer offsets from geometrical offsets
                for(int i=0; i<HalfStencilSize; ++i){
                    for(int orientation = 0; orientation<2; ++orientation, ++yIndex){
                        const IndexType y = orientation ? x-offsets[i] : x+offsets[i];
                        if(this->region.IsInside(y)){
                            *yIndex = this->BufferIndex(y);
                        } else {
                            // Neumann boundary conditions.
                            *yIndex = this->OutsideBufferIndex();
                        } // if y
                    } // for eps
                } // for i
                return stencil;
            }
        protected:
            RegionType region;
            IndexType prod;
            SpacingType invSpacing;
            InternalSizeT OutsideBufferIndex() const {return NumericTraits<InternalSizeT>::max();}
        };
        
        auto filter = UnaryFunctorWithIndexImageFilter<TensorImageType, StencilImageType, StencilFunctor >::New();
        filter->SetInput(GetInputTensor());
        filter->GetFunctor().Initialize(region, GetInputTensor()->GetSpacing());
        filter->Update();
        stencilImage = filter->GetOutput();
        
        
        //setup diagonal coefficients. Cannot be parallelized due to non-local modifications of diagBuffer.
        
        diagonalCoefficients = ScalarImageType::New();
        diagonalCoefficients->CopyInformation(GetInputTensor());
        diagonalCoefficients->SetRegions(GetRequestedRegion());
        diagonalCoefficients->Allocate();
        diagonalCoefficients->FillBuffer(ValueType(0));

        ImageRegionConstIterator<StencilImageType> stencilIt(stencilImage,region);
        ImageRegionIterator<ScalarImageType> diagIt(diagonalCoefficients, region);
        ValueType * diagBuffer = diagonalCoefficients->GetBufferPointer();

        for(stencilIt.GoToBegin(), diagIt.GoToBegin();
            !stencilIt.IsAtEnd();
            ++stencilIt, ++diagIt)
            for(int i=0; i<StencilSize; ++i) {
                const InternalSizeT yIndex = stencilIt.Value().first[i];
                if(yIndex!=OutsideBufferIndex()){
                    const ValueType coefficient = stencilIt.Value().second[i/2];
                    diagIt.Value()      += coefficient;
                    diagBuffer[yIndex]  += coefficient;
                } // if y
            } // for i
        // for stencilIt, diagIt
    }
    
    template< typename TImage, typename TTensorImage>
    typename LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::ValueType
    LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::MaxStableTimeStep(){
        typedef MinimumMaximumImageCalculator<ScalarImageType> MaxCalculatorType;
        typename MaxCalculatorType::Pointer maximumCalculator = MaxCalculatorType::New();
        maximumCalculator->SetImage(diagonalCoefficients);
        maximumCalculator->SetRegion(GetRequestedRegion());
        maximumCalculator->ComputeMaximum();
        return 1./maximumCalculator->GetMaximum();
    }
    
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::SetMaxDiffusionTime(ValueType time){
        if(time<0)
            itkExceptionMacro("diffusion time must be finite and positive");
        diffusionTime = time;
    }
    
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::SetRatioToMaxStableTimeStep(ValueType ratio){
        if(ratio<=0 || ratio>1)
            itkExceptionMacro("Ratio to max time step " << ratio << "should be within ]0,1]");
        ratioToMaxStableTimeStep=ratio;
    }

    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::SetMaxNumberOfTimeSteps(int n){
        if(n<=0)
            itkExceptionMacro("Max number of time steps must be positive");
        maximumNumberOfTimeSteps=n;
    }

    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::ImageUpdateLoop(){
        ValueType delta = MaxStableTimeStep() * ratioToMaxStableTimeStep;
        int n = ceil(diffusionTime / delta);
        if(n>maximumNumberOfTimeSteps) {
            n=maximumNumberOfTimeSteps;
            m_EffectiveDiffusionTime = n*delta;
        } else {
            delta = diffusionTime/n;
            m_EffectiveDiffusionTime = diffusionTime;
        }
        m_EffectiveNumberOfTimeSteps=n;
        
        // Extraction of the region of interest is required for image buffer access.
        typedef ExtractImageFilter<ImageType, ScalarImageType> InputCasterType;
        typename InputCasterType::Pointer inputCaster = InputCasterType::New();
        inputCaster->SetInput(GetInputImage());
        inputCaster->SetExtractionRegion(GetRequestedRegion());
        inputCaster->SetDirectionCollapseToIdentity();
        inputCaster->Update();
        previousImage = inputCaster->GetOutput();
        
        nextImage = ScalarImageType::New();
        nextImage->CopyInformation(GetInputTensor());
        nextImage->SetRegions(GetRequestedRegion());
        nextImage->Allocate();
        
        for(int k=0; k<n; ++k){
            ImageUpdate(delta);
            std::swap(previousImage, nextImage);
        }
        
        typedef CastImageFilter<ScalarImageType, ImageType> OutputCasterType;
        typename OutputCasterType::Pointer outputCaster = OutputCasterType::New();
        outputCaster->SetInput(nextImage);
        
        outputCaster->Update();
        this->GraftOutput(outputCaster->GetOutput());
    }
    
    template< typename TImage, typename TTensorImage>
    void LinearAnisotropicDiffusionLBRImageFilter<TImage, TTensorImage>::ImageUpdate(ValueType delta){
    
    //Setting up iterators
    
    ImageRegion<Dimension> region = GetRequestedRegion();
    
    ImageRegionConstIterator<ScalarImageType>   inputIt(previousImage,region);
    ImageRegionIterator<ScalarImageType>        outputIt(nextImage,region);
    
    const ValueType * inputBuffer =  previousImage->GetBufferPointer();
    ValueType       * outputBuffer = nextImage->GetBufferPointer();
        
    ImageRegionConstIterator<ScalarImageType>           diagIt(diagonalCoefficients,    region);
    ImageRegionConstIterator<StencilImageType>          stencilIt(stencilImage,         region);
    
    // Doing the product.
    nextImage->FillBuffer(0.);
    
    // Taking care of Off-Diagonal matrix elements. Cannot be parallelized due to non-local modifications of outputBuffer
    for(inputIt.GoToBegin(), outputIt.GoToBegin(), stencilIt.GoToBegin();
        !inputIt.IsAtEnd();
        ++inputIt, ++outputIt, ++stencilIt){
        for(int i=0; i<StencilSize; ++i){
        const InternalSizeT   yIndex = stencilIt.Value().first[i];
            if(yIndex!=OutsideBufferIndex()){
                const ValueType coefficient = stencilIt.Value().second[i/2];
                outputIt.Value()        +=  coefficient * (inputBuffer[yIndex]);
                outputBuffer[yIndex]   += coefficient * inputIt.Value();
            }
        }
    }
  
    // parallelization of this specific step does not help much
    struct FunctorType {
        ValueType delta;
        ValueType operator()(ValueType output, ValueType input, ValueType diag){
            return this->delta*output + (ValueType(1)-this->delta*diag)*input;
        }
    };
    
    typedef TernaryFunctorImageFilter<ScalarImageType, ScalarImageType, ScalarImageType, ScalarImageType, FunctorType> ImageFunctorType;
    typename ImageFunctorType::Pointer imageFunctor = ImageFunctorType::New();
    imageFunctor->SetInput1(nextImage);
    imageFunctor->SetInput2(previousImage);
    imageFunctor->SetInput3(diagonalCoefficients);
    imageFunctor->GetFunctor().delta = delta;
    
    assert(imageFunctor->CanRunInPlace());
    imageFunctor->InPlaceOn();
    imageFunctor->Update();
    nextImage=imageFunctor->GetOutput();
      
    /*
    // Old Serial version for diagonal elements
    for(inputIt.GoToBegin(), outputIt.GoToBegin(), diagIt.GoToBegin();
        !inputIt.IsAtEnd();
        ++inputIt, ++outputIt, ++diagIt)
        outputIt.Value() = delta*outputIt.Value() + (1-delta*diagIt.Value())*inputIt.Value();
    */
    }
    
// **************************** subclass SSRT_Traits call method **************************

// 2D stencil construction.
template<typename TI,typename TTI> template<typename Dummy>
void LinearAnisotropicDiffusionLBRImageFilter<TI,TTI>::SSRT_Traits<2,Dummy>::
GetDiffusionStencil(const SSRT &D, StencilOffsetsType &Offsets, StencilCoefficientsType &Coefficients){
    
    typedef typename SSRT::ShortOffsetType ShortOffsetType;
    typedef Vector<ValueType,Dimension> VectorType;
    
    VectorType Basis[Dimension+1];
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<Dimension; ++j)
            Basis[i][j] = (i==j);
    
    D.ObtuseSuperBase(Basis);
    Basis[2] = -(Basis[0]+Basis[1]);
    
    for(int i=0; i<3; ++i){
        Coefficients[i] = (-0.5)*D.ScalarProduct(Basis[(i+1)%3],Basis[(i+2)%3]);
        assert(Coefficients[i]>=0);
        
        Offsets[i][0] = static_cast<OffsetValueType>(-Basis[i][1]);
        Offsets[i][1] = static_cast<OffsetValueType>( Basis[i][0]);
    }
    
} // 2D stencil construction


// 3D stencil construction

template<typename TI,typename TTI> template<typename Dummy>
void LinearAnisotropicDiffusionLBRImageFilter<TI,TTI>::SSRT_Traits<3,Dummy>::
GetDiffusionStencil(const SSRT &D, StencilOffsetsType &Offsets, StencilCoefficientsType &Coefficients){
    
    typedef typename SSRT::ShortOffsetType ShortOffsetType;
    
    using std::make_pair;
    using std::max;
    using std::min;
    
    
    typedef Vector<ValueType,Dimension> VectorType;
    VectorType Basis[Dimension+1];
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<Dimension; ++j)
            Basis[i][j] = (i==j);
    
    // Construction of a D-obtuse superbase.
    D.ObtuseSuperBase(Basis);
    Basis[3] = -(Basis[0]+Basis[1]+Basis[2]);
    
    // Computation of the weights
    itk::SymmetricSecondRankTensor<ValueType,Dimension+1> Weights;
    for(int i=1; i<Dimension+1; ++i)
        for(int j=0; j<i; ++j)
            Weights(i,j) = (-0.5)*D.ScalarProduct(Basis[i],Basis[j]);
    
    // Now that the obtuse superbasis has been created, generate the stencil.
    // First get the dual basis. Obtained by computing the comatrix of Basis[1..Dimension].
    
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<Dimension; ++j)
            Offsets[i][j] = Basis[(i+1)%Dimension][(j+1)%Dimension]*Basis[(i+2)%Dimension][(j+2)%Dimension]
            - Basis[(i+2)%Dimension][(j+1)%Dimension]*Basis[(i+1)%Dimension][(j+2)%Dimension];
    
    Offsets[Dimension  ] = Offsets[0]-Offsets[1];
    Offsets[Dimension+1] = Offsets[0]-Offsets[2];
    Offsets[Dimension+2] = Offsets[1]-Offsets[2];
    
    // The corresponding coefficients are given by the scalar products.
    for(int i=0; i<Dimension; ++i)
        Coefficients[i]=Weights(i,3);
    
    Coefficients[Dimension]   = Weights(0,1);
    Coefficients[Dimension+1] = Weights(0,2);
    Coefficients[Dimension+2] = Weights(1,2);
} // 3D stencil construction

}// end namespace


#endif
