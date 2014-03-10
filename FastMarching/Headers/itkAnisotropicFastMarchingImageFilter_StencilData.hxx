//
//  itkAnisotropicFastMarchingImageFilter_StencilData.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 14/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_itkAnisotropicFastMarchingImageFilter_StencilData_hxx
#define ITKFM_itkAnisotropicFastMarchingImageFilter_StencilData_hxx
 
namespace itk {
    
    template<typename TN,typename TSS>
    typename AnisotropicFastMarchingImageFilter<TN,TSS>::DirectStencilType
    AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::Direct
    (const IndexType & P) const
    {
        CompressedStartIt.SetIndex(P);
        const StencilStartType begin = CompressedStartIt.Value();
        ++CompressedStartIt;
        const StencilStartType end = CompressedStartIt.IsAtEnd() ? 
            StencilStartType(Compressed.size()) : CompressedStartIt.Value();
        
        const StencilStartType size=end-begin;
        return DirectStencilType(&Compressed[begin],size);
    }
    
    
    template<typename TN,typename TSS>
    const typename AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::ReverseStencilType
    AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::Reverse
    (const IndexType & P)
    {
        ReversedStartIt.SetIndex(P);
        const StencilStartType begin=ReversedStartIt.Value();
        ++ReversedStartIt;
        const StencilStartType end= ReversedStartIt.IsAtEnd() ? 
            StencilStartType(Reversed.size()) : ReversedStartIt.Value();
        
        ReverseStencilType rStencil(end-begin, &Reversed[begin], &Reversed_i[begin]);
        return rStencil;
    }
    
    
    // ************************ Initialization **************************
    
    template<typename TN, typename TSS>
    void
    AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::Init
    (const NormImageType* Norms)
    {
        BufferedRegion=Norms->GetBufferedRegion();
        Allocate();
        SetCompressed(Norms);
        SetReversed();
    }
    
    template<typename TN, typename TSS>
    void
    AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::Allocate()
    {
        CompressedStart=StartImageType::New();
        CompressedStart->SetRegions(BufferedRegion);
        CompressedStart->Allocate();
        
        CompressedStartIt=ImageRegionIterator<StartImageType>(CompressedStart,BufferedRegion);
        
        ReversedStart=StartImageType::New();
        ReversedStart->SetRegions(BufferedRegion);
        ReversedStart->Allocate();
        
        ReversedStartIt=ImageRegionIterator<StartImageType>(ReversedStart,BufferedRegion);
        
    }
    
    template <typename TN, typename TSS>
    void
    AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::SetCompressed
    (const NormImageType* Norms)
    {
        ImageRegionConstIterator<NormImageType> NormIt(Norms,BufferedRegion);
        
        for(NormIt.GoToBegin(), CompressedStartIt.GoToBegin();
            !NormIt.IsAtEnd();
            ++NormIt, ++CompressedStartIt){
            CompressedStartIt.Set(Compressed.size());
            NormType(NormIt.Value(),Norms->GetSpacing()).GetCompressedStencil(Compressed);
        }
    }
    
    template <typename TN, typename TSS>
    void
    AnisotropicFastMarchingImageFilter<TN,TSS>::StencilData::SetReversed()
    {
        ReversedStart->FillBuffer(0);
        ImageRegionConstIteratorWithIndex<StartImageType> CompIt(CompressedStart,BufferedRegion);
        
        for(CompIt.GoToBegin(); !CompIt.IsAtEnd(); ++CompIt){
            const IndexType & P = CompIt.GetIndex();
            
            const DirectStencilType stencil = Direct(P);
            unsigned int size = stencil.size();
            
            for(int i=0; i<size; ++i){
                const ShortOffsetType ShortOffset = stencil[i];
                const IndexType Q = P+NormType::ShortOffset_to_Offset(ShortOffset);
                if(!BufferedRegion.IsInside(Q)) continue;
                
                //                if(Q[0]==0 && Q[1]==0) cout << "Q " << Q << "; P " << P << "; Offset " << ShortOffset << endl;
                
                ReversedStart->GetPixel(Q)++;
            }
        }
        
        StencilStartType sum=0;
        for(ReversedStartIt.GoToBegin(); !ReversedStartIt.IsAtEnd(); ++ReversedStartIt){
            const StencilStartType size = ReversedStartIt.Value();
            ReversedStartIt.Set(sum);
            sum+=size;
        }
        
        Reversed.resize(sum); Reversed_i.resize(sum);
        
        for(CompIt.GoToBegin(); !CompIt.IsAtEnd(); ++CompIt){
            const IndexType & P = CompIt.GetIndex();
            
            const DirectStencilType stencil = Direct(P);
            unsigned int size = stencil.size();
            
            for(int i=0; i<size; ++i){
                const ShortOffsetType ShortOffset = stencil[i];
                const IndexType Q = P+NormType::ShortOffset_to_Offset(ShortOffset);
                if(!BufferedRegion.IsInside(Q)) continue;
                StencilStartType & RevStart = ReversedStart->GetPixel(Q);
                Reversed[RevStart]=ShortOffset;
                Reversed_i[RevStart]=i;
                RevStart++;
            }
        }
        
        for(ReversedStartIt.GoToEnd(), --ReversedStartIt; !ReversedStartIt.IsAtBegin(); --ReversedStartIt){
            --ReversedStartIt;
            const StencilStartType RevStart=ReversedStartIt.Value();
            ++ReversedStartIt;
            ReversedStartIt.Set(RevStart);
        }
        ReversedStartIt.Set(0);
    }
    
} //namespace itk

#endif
