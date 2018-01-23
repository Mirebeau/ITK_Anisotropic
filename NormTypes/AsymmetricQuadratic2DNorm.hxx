//
//  AsymmetricQuadratic2DNorm.hxx
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 02/10/2017.
//
//

#ifndef AsymmetricQuadratic2DNorm_hxx
#define AsymmetricQuadratic2DNorm_hxx

namespace itk {
    
    template<typename TC,typename TSO>
    typename AsymmetricQuadratic2DNorm<TC,TSO>::CovariantVectorType
    AsymmetricQuadratic2DNorm<TC,TSO>::Gradient(const VectorType &u) const
    {
        VectorType v = GetM()*u
        +std::max(GetOmega()*u,ScalarType(0.))*GetOmega();
        v /= sqrt(u*v);
        return CovariantVectorType(v.GetDataPointer());
    }
    
    template<typename TC,typename TSO>
    bool
    AsymmetricQuadratic2DNorm<TC,TSO>::IsAcute(const VectorType & u,
                                               const VectorType & v) const {
        
//        Equivalent to : return Gradient(v)*u>=0 && Gradient(u)*v>=0;

        const ScalarType
        scalUV = u*(GetM()*v),
        scalUO = u*GetOmega(),
        scalVO = v*GetOmega();

        return scalUV+std::min(scalUO*std::max(ScalarType(0.),scalVO),
                               scalVO*std::max(ScalarType(0.),scalUO))
        >=ScalarType(0.);
    }
    
    template<typename TC,typename TSO> template<typename AcceptedImageType>
    TC AsymmetricQuadratic2DNorm<TC,TSO>::Hopf_Lax(const DistanceImageType * d, const IndexType &x,
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
        
        const VectorType & om = GetOmega();
        const ScalarType y_om = Vy*om;
        
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
            
            const ScalarType z_om = Vz*om;
            
            //quadratic systems to solve...
            Riemannian2DNorm M(GetM().SquaredNorm(Vy),
                               GetM().ScalarProduct(Vy,Vz),
                               GetM().SquaredNorm(Vz));
            VectorType ones; ones.Fill(1);
            VectorType dist; dist[0]=Dy; dist[1]=Dz;
            
            for(int om_pos=0; om_pos<=1; ++om_pos){
                // This loop tries to solve the Riemannian Hopf-Lax operator with SSRTs
                // GetM() and GetM()+GetOmega() x GetOmega(). The suitable update is selected.
                
                // Wether the scalar product with GetOmega is positive.
                if(om_pos==0 && y_om>0 && z_om>0) continue;
                if(om_pos==1 && y_om<0 && z_om<0) continue;
                
                if(om_pos){
                    M(0,0)+=y_om*y_om;
                    M(0,1)+=y_om*z_om;
                    M(1,1)+=z_om*z_om;
                }
                
                const Riemannian2DNorm N=M.GetInverse();
                
                const ValueType
                ones2 = N.SquaredNorm(ones),
                dist2 = N.SquaredNorm(dist),
                ones_dist = N.ScalarProduct(ones, dist);
                
                const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                if(delta<0) continue;
                const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                
                if(DTest>=D) continue;
                
                const VectorType w=N*(DTest*ones-dist);
                if(w[0]<=0 || w[1]<=0) continue;
                
                const ScalarType t_om = w[0]*y_om+w[1]*z_om;
                if((om_pos==0 && t_om>0) || (om_pos==1 && t_om<0)) continue;
                
                D=DTest;
                gradient_weights[0]=w[0]/(w[0]+w[1]);
                //gradient_weights[1] defined implicitly
                
                //gradient_indices[0]=y_i; //already done
                gradient_indices[1]=z_i;
            }
        }
        return D;
    }
    
    
}

#endif /* AsymmetricQuadratic2DNorm_hxx */
