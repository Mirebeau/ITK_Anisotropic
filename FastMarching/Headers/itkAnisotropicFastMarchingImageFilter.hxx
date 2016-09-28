//
//  itkAnisotropicFastMarchingImageFilter.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 08/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_itkAnisotropicFastMarchingImageFilter_hxx
#define ITKFM_itkAnisotropicFastMarchingImageFilter_hxx
 
namespace itk {
    
    template<typename TN, typename ist>
    void AnisotropicFastMarchingImageFilter<TN,ist>
    ::Initialize(LevelSetImageType *output)
    {
        Superclass::Initialize(output);
        const NormImageType *input = this->GetInput();
        
        Stencils.Init(input);
        if(GetGenerateUpwindGradient())
            Gradient.Init(this->m_BufferedRegion);
    }

    
    template<typename TN, typename ist>
    void AnisotropicFastMarchingImageFilter<TN,ist>
    ::UpdateNeighbors(const IndexType & index, const NormImageType *input, LevelSetImageType *output)
    {
        YoungestAliveIndex=index;
        typedef typename StencilData::ReverseStencilType ReverseStencilType;
        ReverseStencilType rStencil = Stencils.Reverse(index);
        
        for(int pos=0; pos<(int)rStencil.size(); ++pos){ 
            const IndexType neighIndex = 
            index-NormType::ShortOffset_to_Offset(rStencil[pos]);
            
            YoungestAliveIndex_i = rStencil._i(pos);
            
            const unsigned char label = this->GetLabelImage()->GetPixel(neighIndex);
            
            if ( ( label != this->AlivePoint ) &&
                ( label != this->InitialTrialPoint ) &&
                ( label != this->OutsidePoint ) )
                UpdateValue(neighIndex, input, output);
        }
    }
    
    
    
    template<typename TN, typename ist>
    double AnisotropicFastMarchingImageFilter<TN,ist>
    ::UpdateValue(const IndexType &x, const NormImageType * Norms, LevelSetImageType * Values)
    {
        typedef typename NormType::Stencil DirectStencilType;
        DirectStencilType xStencil = Stencils.Direct(x);
         
        //Take spacing and normalization into account.
        const NormType Norm(Norms->GetPixel(x),Norms->GetSpacing()*this->GetNormalizationFactor());
        
        //Accepted, or Alive, values are those that can be used in the computaion of new values.
        struct AcceptedImageType : public Self {
            const typename Superclass::OutputRegionType & GetBufferedRegion() const {
                return this->m_BufferedRegion;
            }
            bool GetPixel(const IndexType &z) const {
                return this->GetLabelImage()->GetPixel(z) == this->AlivePoint;
            }
        };
        
        AcceptedImageType *Accepted = static_cast<AcceptedImageType*>(this);
        
        //Container to store the upwind gradient.
        typename GradientRawData::RawType Grad;
        

        const ValueType EstimatedValue = Norm.Hopf_Lax(Values,x,Accepted,YoungestAliveIndex_i,
                                          xStencil,
                                          Grad.first, Grad.second);
        
        // Update distance
        const ValueType OldValue = Values->GetPixel(x);
        
        if ( EstimatedValue < OldValue ) 
        {
            // write solution to m_OutputLevelSet
            Values->SetPixel(x, EstimatedValue);
            
            // insert point into trial heap
            this->GetLabelImage()->SetPixel(x, Superclass::TrialPoint);
            typename Superclass::AxisNodeType node;
            node.SetValue(EstimatedValue);
            node.SetIndex(x );
            this->m_TrialHeap.push(node);
            
            // Save raw gradient.
            if(GetGenerateUpwindGradient())
                Gradient.Raw->GetPixel(x)=Grad;

        }
        
        return EstimatedValue;
    };
    
    
    
    // ************************ Gradient ***********************
    
    template <typename TN, typename TSS>
    void
    AnisotropicFastMarchingImageFilter<TN,TSS>::GradientRawData::Init(ImageRegionType BufferedRegion)
    {
        Raw=RawImageType::New();
        Raw->SetRegions(BufferedRegion);
        Raw->Allocate();
        
        RawType NullGradient;
        NullGradient.second.SetNullWeights();
        Raw->FillBuffer(NullGradient);
    }
    
    template <typename TN, typename TSS> 
    typename AnisotropicFastMarchingImageFilter<TN,TSS>::GradientImageType::Pointer 
    AnisotropicFastMarchingImageFilter<TN,TSS>::ComputeGradientImage()
    {
        if(!GetGenerateUpwindGradient()){
            itkExceptionMacro("AnisotropicFastMarchingImageFilter::ComputeGradientImage error : flag GenerateUpwindGradient was not set !");
        }
        
        const ImageRegionType & BufferedRegion=Stencils.BufferedRegion;
        
        typename GradientImageType::Pointer GradientImage;
        GradientImage=GradientImageType::New();
        GradientImage->SetRegions(BufferedRegion);
        GradientImage->Allocate();
        
        const NormImageType * Norms = this->GetInput();
        typedef typename StencilData::StencilStartType StencilStartType;
        typedef typename StencilData::CompressedOffsetType CompressedOffsetType;
        
        ImageRegionIterator<GradientImageType> GradIt(GradientImage,BufferedRegion);
        ImageRegionConstIterator<typename GradientRawData::RawImageType> GradRawIt(Gradient.Raw,BufferedRegion);
        ImageRegionConstIterator<typename StencilData::StartImageType> StartIt(Stencils.CompressedStart,BufferedRegion);
        ImageRegionConstIterator<NormImageType> NormIt(Norms,BufferedRegion);
        
        for(GradIt.GoToBegin(), GradRawIt.GoToBegin(), StartIt.GoToBegin(), NormIt.GoToBegin();
            !GradIt.IsAtEnd();
            ++GradIt, ++GradRawIt, ++NormIt){ //StartIt incremented inside
            const StencilStartType begin = StartIt.Value();
            ++StartIt;
            const StencilStartType end = StartIt.IsAtEnd() ? StencilStartType(Stencils.Compressed.size()) : StartIt.Value();
            
            const CompressedOffsetType * pCompressedStencil = &(Stencils.Compressed[begin]);
            const unsigned int compressedStencilSize = end-begin;
            
            const typename GradientRawData::IndicesType indices=GradRawIt.Value().first;
            const typename GradientRawData::WeightsType weights=GradRawIt.Value().second;
            
            VectorType grad; 
            grad.Fill(0);
            if(weights.IsNullWeights()) continue;
            
            DirectStencilType stencil(pCompressedStencil, compressedStencilSize);
            
            for(int i=0; i<Dimension; ++i){
                const ValueType w=weights(i);
                if(w==0) continue;
                const ShortOffsetType u = stencil[indices[i]];
                grad += w*VectorType(u);
            }
            
            const ValueType gradNorm = NormType(NormIt.Value(),Norms->GetSpacing()).Norm(grad);
            if(gradNorm>0) grad /= -gradNorm;
            
            for(int i=0; i<Dimension; ++i) GradIt.Value()[i]=grad[i];
        }
        return GradientImage;
    }
    
} //namespace itk



#endif
