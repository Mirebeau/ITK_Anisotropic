//
//  Finsler2DNorm.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 07/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_Finsler2DNorm_hxx
#define ITKFM_Finsler2DNorm_hxx

namespace itk {
    
    template<typename TC,typename TSO> 
    typename Finsler2DNorm<TC,TSO>::CovariantVectorType
    Finsler2DNorm<TC,TSO>::Gradient(const VectorType &u) const
    {
        const ValueType n=GetM().Norm(u); 
        assert(n>0);
        VectorType grad = GetM()*VectorType(u)/n-GetOmega();         
        return CovariantVectorType(grad.GetDataPointer());
    }
    
    
    template<typename TC,typename TSO> template<typename AcceptedImageType>
    TC Finsler2DNorm<TC,TSO>::Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                                       const AcceptedImageType * b, unsigned int y_i,
                                       const Stencil &stencil,
                                       IndicesType &gradient_indices, WeightsType &gradient_weights) const
    {
        const ShortOffsetType Iy=stencil[y_i];
        const VectorType Vy(Iy);
        const IndexType y=x+this->ShortOffset_to_Offset(Iy);
        const ValueType Dy=d->GetPixel(y);
        
        ValueType D = Dy+Norm(Vy);
        gradient_indices[0]=y_i;
        gradient_indices[1]=this->IndexToIgnore();
        
        gradient_weights[0]=1;
        //gradient_weights[1] defined implicitly
        
        const unsigned int stencilSize = stencil.size();
        for(int e=-1; e<=1; e+=2){
            int z_i = int(y_i)+e;
            z_i=(z_i+stencilSize)%stencilSize;
            const ShortOffsetType Iz=stencil[z_i];
            const VectorType Vz(Iz);
            const IndexType z=x+this->ShortOffset_to_Offset(Iz);
            
            if(!b->GetBufferedRegion().IsInside(z)) continue;
            if(!b->GetPixel(z)) continue;  //Point z must be accepted. 
            const ValueType Dz=d->GetPixel(z); //Since z is accepted, Dz < infinity
            
            //quadratic system to solve...
            
            const Riemannian2DNorm M(GetM().SquaredNorm(Vy),GetM().ScalarProduct(Vy,Vz),GetM().SquaredNorm(Vz));
            const Riemannian2DNorm N=M.GetInverse();
            
            VectorType ones; ones.Fill(1);
            VectorType dist; dist[0]=Dy - GetOmega()*Vy; dist[1]=Dz - GetOmega()*Vz;
            
            const ValueType ones2 = N.SquaredNorm(ones), dist2 = N.SquaredNorm(dist), ones_dist = N.ScalarProduct(ones, dist);
            const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
            if(delta<0) continue;
            const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
            
            if(DTest>=D) continue;
            
            const VectorType w=N*(DTest*ones-dist);
            if(w[0]<=0 || w[1]<=0) continue;
            
            D=DTest;
            gradient_weights[0]=w[0]/(w[0]+w[1]);
            //gradient_weights[1] defined implicitly

            //gradient_indices[0]=y_i; //already done
            gradient_indices[1]=z_i;
        }
        return D;
        
    }
    
}//namespace itk
#endif
