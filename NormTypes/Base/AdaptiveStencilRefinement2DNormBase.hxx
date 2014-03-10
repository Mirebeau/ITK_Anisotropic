//
//  AdaptiveStencilRefinement.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 05/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_AdaptiveStencilRefinement2DNormBase_hxx
#define ITKFM_AdaptiveStencilRefinement2DNormBase_hxx

namespace itk {
    
    template<bool CentroSymmetric, typename TSO> template<typename TNorm>
    void
    AdaptiveStencilRefinement2DNormBase<CentroSymmetric,TSO>::CompressedStencil
    (const TNorm &N, std::vector<CompressedOffsetType> &l) const
    {
        std::vector<CompressedOffsetType> m;
        {
            ShortOffsetType u;
            if(!CentroSymmetric){ 
                u[0]= 1; u[1]= 0; m.push_back(u);
                u[0]= 0; u[1]=-1; m.push_back(u);
            }
            u[0]=-1; u[1]= 0; m.push_back(u);
            u[0]= 0; u[1]= 1; m.push_back(u);
            
            u[0]= 1; u[1]=0; l.push_back(u);
        }
        
        while(!m.empty()){
            const CompressedOffsetType &u=l.back();
            const CompressedOffsetType &v=m.back();
            if(N.IsAcute(u,v)){
                l.push_back(v);
                m.pop_back();
            } else {
                m.push_back(u+v);
            }
        }
        l.pop_back();
    }
}

#endif
