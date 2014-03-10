//
//  Riemannian2DNorm.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 05/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_Riemannian2DNorm_hxx
#define ITKFM_Riemannian2DNorm_hxx


namespace itk
{
    /*
     note that 
     typename Riemannian2DNorm<TC,TSO>::ValueType 
     is precisely
     TC
     */ 
    
    //Hopf-Lax update operator used in the fast marching algorithm
    template<typename TC, typename TSO>
    template<typename AcceptedImageType>
    TC
    Riemannian2DNorm<TC,TSO>::Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                                       const AcceptedImageType * b, unsigned int y_i,
                                       const Stencil &stencil,
                                       IndicesType &gradient_indices, WeightsType &gradient_weights) const
    {
        const ShortOffsetType Iy=stencil[y_i];
        const IndexType y=x+this->ShortOffset_to_Offset(Iy);
        const VectorType Vy(Iy);
        const ValueType Dy=d->GetPixel(y);
        
        ValueType D = Dy+this->Norm(Vy);
        gradient_indices[0]=y_i;
        gradient_indices[1]=this->IndexToIgnore();
        
        gradient_weights[0]=1;
        // gradient_weights[1] defined implicitly
        
        const unsigned int stencilSize = stencil.size();
        for(int e=-1; e<=1; e+=2){
            int z_i = int(y_i)+e;
            z_i = (z_i+stencilSize)%stencilSize;
            const ShortOffsetType Iz=stencil[z_i];
            const VectorType Vz(Iz);
            const IndexType z=x+this->ShortOffset_to_Offset(Iz);
            
            if(!b->GetBufferedRegion().IsInside(z)) continue;
            if(!b->GetPixel(z)) continue;  //Point z must be accepted. 
            const ValueType Dz=d->GetPixel(z); //Since z is accepted, Dz < infinity
            
            //quadratic system to solve...
            
            const Riemannian2DNorm M(this->SquaredNorm(Vy),this->ScalarProduct(Vy,Vz),this->SquaredNorm(Vz));
            const Riemannian2DNorm N=M.GetInverse();
            
            VectorType ones; ones.Fill(1);
            VectorType dist; dist[0]=Dy; dist[1]=Dz;
            
            const ValueType ones2 = N.SquaredNorm(ones), dist2 = N.SquaredNorm(dist), ones_dist = N.ScalarProduct(ones, dist);
            const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
            if(delta<0) continue;
            const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
            
            if(DTest>=D) continue;
            
            const VectorType w=N*(DTest*ones-dist);
            if(w[0]<=0 || w[1]<=0) continue;
            
            D=DTest;
            gradient_weights[0]= w[0]/(w[0]+w[1]);
            //gradient_weights[1] defined implicitly
            
            //gradient_indices[0]=y_i; already done
            gradient_indices[1]=z_i;
        }
        return D;
        
    }
    
    template<typename TC, typename TSO>
    void
    Riemannian2DNorm<TC,TSO>::ReducedBasis(VectorType Basis[Dimension]) const
    {
        VectorType &b0=Basis[0], &b1=Basis[1];
        
        assert( abs( b0[0]*b1[1]-b0[1]*b1[0] ) ==1 );
        using std::swap;
        
        // Squared norms
        ValueType n0=this->SquaredNorm(b0), n1=this->SquaredNorm(b1);
        
        if(n1<n0){swap(b0,b1);swap(n0,n1);}  //Sorting by norm
        
        while(true){
            b1 -=  ValueType(int(this->ScalarProduct(b0,b1)/n0))  *b0;
            n1= this->SquaredNorm(b1);
            if(n0<=n1) break;
            swap(b0,b1); swap(n0,n1);
        }
    }
    
    template<typename TC, typename TSO>
    void
    Riemannian2DNorm<TC,TSO>::ObtuseSuperBase(VectorType SuperBase[Dimension]) const
    {
        ReducedBasis(SuperBase);
        if(this->ScalarProduct(SuperBase[0],SuperBase[1]) > 0)
            SuperBase[1]*=-1;
    }
    
} //namespace itk

#endif
