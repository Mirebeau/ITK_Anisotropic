//
//  AdaptiveStencilRefinement2DNormBase.h
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 05/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_AdaptiveStencilRefinement2DNormBase_h
#define ITKFM_AdaptiveStencilRefinement2DNormBase_h 

#include "itkImage.h"

namespace itk {
    
    template<bool CentroSymmetric, typename TShortOffsetValue> class AdaptiveStencilRefinement2DNormBase
    {
    public:
        static const unsigned int Dimension = 2; 
        typedef TShortOffsetValue ShortOffsetValueType;
        
        typedef Vector<ShortOffsetValueType,Dimension> ShortOffsetType;
        typedef ShortOffsetType CompressedOffsetType;
        typedef ShortOffsetType IndicesType;
        typedef Offset<Dimension> OffsetType;
        
        static const ShortOffsetValueType IndexToIgnore(){return NumericTraits<ShortOffsetValueType>::max();}
        
        static OffsetType ShortOffset_to_Offset(const ShortOffsetType &u){OffsetType U; for(int i=0; i<(int)Dimension; ++i) U[i]=u[i]; return U;}
                
        
        struct Stencil {
            const ShortOffsetType * const pCompressedStencil;
            const int compressedStencilSize;
            inline int size() const
            {
                if(CentroSymmetric) return 2*compressedStencilSize;
                else return compressedStencilSize;
            }
            inline ShortOffsetType operator[](int i) const
            {
                assert(0 <= i && i < size());
                if(CentroSymmetric){
                    if(i<compressedStencilSize) return pCompressedStencil[i];
                    else return -pCompressedStencil[i-compressedStencilSize];
                } else
                    return pCompressedStencil[i];
            }
            Stencil(const ShortOffsetType * pCompressedStencil_, int compressedStencilSize_):
            pCompressedStencil(pCompressedStencil_),compressedStencilSize(compressedStencilSize_){};
        };
        
    protected:
        template<typename TNorm> void CompressedStencil(const TNorm &N, std::vector<CompressedOffsetType> &l) const;        
    };
}


#ifndef ITK_MANUAL_INSTANTIATION
#include "AdaptiveStencilRefinement2DNormBase.hxx"
#endif


#endif
