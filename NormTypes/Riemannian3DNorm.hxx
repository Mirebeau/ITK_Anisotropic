//
//  Riemannian3DNorm.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 07/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_Riemannian3DNorm_hxx
#define ITKFM_Riemannian3DNorm_hxx

namespace itk  {
    
#ifdef SuperBase
    
    template<typename TC, typename TSO>
    void Riemannian3DNorm<TC,TSO>::GetCompressedStencil(std::vector<CompressedOffsetType> &l) const
    {
        assert(IsDefinite());
        if(this->IsDiagonal()) return;
                
        VectorType Basis[Dimension];
        for(int i=0; i<(int)Dimension; ++i)
            for(int j=0; j<(int)Dimension; ++j)
                Basis[i][j]= (i==j);
        /*
         // Use of a static variable to help reduce the number of iterations in lattice basis reduction.
         // Check for safety that it is indeed a basis ? (Calling Determinant)
         // In practice, impact on CPU time is invisible -> removed.
         
         static bool Basis_is_initialized = false;
         if(!Basis_is_initialized){
         for(int i=0; i<Dimension; ++i)
         for(int j=0; j<Dimension; ++j)
         VectorBasis[i][j]= (i==j);
         Basis_is_initialized=true;
         }
         */
        
        ObtuseSuperBase(Basis);
        
        const int size = l.size();
        l.resize(size+Dimension);
        for(int i=0; i<(int)Dimension; ++i){
            l[size+i] = CompressedOffsetType(Basis[i]);
//            CompressedOffsetType &u = l[size+i];
//            const VectorType &U = Basis[i];
//            u=CompressedOffsetType(U);
//            for(int j=0; j<Dimension; ++j)
//                u[j] = static_cast<typename CompressedOffsetType::ValueType>(U[j]);
        }
    }
    
    
#else
    
    template<typename TC, typename TSO>
    void Riemannian3DNorm<TC,TSO>::GetCompressedStencil(std::vector<CompressedOffsetType> &l) const
    {
        if(this->IsDiagonal()) return;
        
        using std::swap;
         
        VectorType VectorBasis[Dimension];
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<Dimension; ++j)
                VectorBasis[i][j]= (i==j);
        /*
         // Use of a static variable to help reduce the number of iterations in lattice basis reduction.
         // Check for safety that it is indeed a basis ? (Calling Determinant)
         // In practice, impact on CPU time is invisible -> removed.
         
        static bool Basis_is_initialized = false;
        if(!Basis_is_initialized){
            for(int i=0; i<Dimension; ++i)
                for(int j=0; j<Dimension; ++j)
                    VectorBasis[i][j]= (i==j);
            Basis_is_initialized=true;
        }
         */
        
        ReducedBasis(VectorBasis);
        
        l.resize(l.size()+Dimension);
        CompressedOffsetType *Basis = &l[l.size()-Dimension];
        for(int i=0; i<Dimension; ++i) Basis[i] = CompressedOffsetType(VectorBasis[i]);
        
        ValueType 
        s01=this->ScalarProduct(Basis[0],Basis[1]), 
        s02=this->ScalarProduct(Basis[0],Basis[2]), 
        s12=this->ScalarProduct(Basis[1],Basis[2]);
        
        if(s01<0) {Basis[1]*=-1; s01*=-1; s12*=-1;}
        if(s02<0) {Basis[2]*=-1; s02*=-1; s12*=-1;}
        if(s12<0) {
            //case of two non-negative scalar products, and one negative. 
            //encoded by a negative determinant of the compressed stencil (i.e. the reduced basis)
            if(Stencil::Determinant(Basis)!=-1)
                for(int i=0; i<Dimension; ++i) Basis[i]*=-1;
            return;
        }
        
        //case of three non-negative scalar products. order them correctly
        if(s12<s01) {swap(Basis[0],Basis[2]); swap(s01,s12);}
        if(s01<s02) {swap(Basis[1],Basis[2]); swap(s01,s02);}
        if(s12<s01) {swap(Basis[0],Basis[2]); swap(s01,s12);}
        
        //encoded by a positive determinant of the compressed stencil (i.e. the reduced basis)
        if(Stencil::Determinant(Basis)!=1) 
            for(int i=0; i<Dimension; ++i) Basis[i]*=-1;
        
        assert(s12 >= s01 && s01 >= s02);
    }
#endif
    
    
    // *************** Stencil sub-structure *************
    
#ifdef SuperBase
    template<typename TC, typename TSO>
    Riemannian3DNorm<TC,TSO>::Stencil::Stencil(const ShortOffsetType * pCompressedStencil_, const unsigned int compressedStencilSize_):
    pCompressedStencil(pCompressedStencil_),compressedStencilSize(compressedStencilSize_)
    {}
#else
    template<typename TC, typename TSO>
    Riemannian3DNorm<TC,TSO>::Stencil::Stencil(const ShortOffsetType * pCompressedStencil_, const unsigned int compressedStencilSize_):
    pCompressedStencil(pCompressedStencil_),compressedStencilSize(compressedStencilSize_),
    topologyCase(compressedStencilSize_ == 0 ? Six :
      Determinant(pCompressedStencil_) == -1 ? Negative :
      Positive)
    {}
    
    template<typename TC, typename TSO>
    int Riemannian3DNorm<TC,TSO>::Stencil::Determinant(const ShortOffsetType Basis[3]){
        int det=0;
        for(int i=0; i<3; ++i) det+= int(Basis[0][i])*
            (int(Basis[1][(i+1)%3])*int(Basis[2][(i+2)%3]) - int(Basis[1][(i+2)%3])*int(Basis[2][(i+1)%3]));
        return det;
    }
#endif

    template<typename TC, typename TSO>
    unsigned int Riemannian3DNorm<TC,TSO>::Stencil::size() const
    {
        return compressedStencilSize == 0 ? 6 : 14;
    }
    
    
    template<typename TC,typename TSO>
    typename Riemannian3DNorm<TC,TSO>::ShortOffsetType
    Riemannian3DNorm<TC,TSO>::Stencil::operator[](int i) const
    {
        if(compressedStencilSize==0){ //topologyCase==Six
            ShortOffsetType x;
            x.Fill(0);
            switch (i) {
                case 0: x[0]= 1; break;
                case 1: x[0]=-1; break;
                case 2: x[1]= 1; break;
                case 3: x[1]=-1; break;
                case 4: x[2]= 1; break;
                case 5: x[2]=-1; break;
                default: assert(false);
                    
            }
            return x;
        }
        const ShortOffsetType * const Basis = pCompressedStencil;
        const ShortOffsetType &u=Basis[0], &v=Basis[1], &w=Basis[2];

#ifdef SuperBase
        const ShortOffsetType x=-(u+v+w);
        if(i<3)  return Basis[i];
        if(i==3) return x;
        i-=4;
        
        if(i<3)  return -Basis[i];
        if(i==3) return -x;
        i-=4;
        
        if(i<1) return Basis[1]+Basis[i];
        i-=1;
        
        if(i<2) return Basis[2]+Basis[i];
        i-=2;

        if(i<3) return x+Basis[i];
        i-=3;
        
#else
        switch (i) {
            case 0: return  u; 
            case 1: return -u; 
            case 2: return  v; 
            case 3: return -v; 
            case 4: return  w; 
            case 5: return -w; 
                
            case 6: return v-u; 
            case 7: return u-v; 
            case 8: return w-u; 
            case 9: return u-w;
        }
        
        if(topologyCase==Negative){ 
            switch (i) {
                case 10: return w+v;
                case 11: return -w-v;
                case 12: return w+v-u;
                case 13: return u-w-v;
            }
        } else {
            switch (i) {
                case 10: return w-v+u;
                case 11: return -w+v-u;
                case 12: return w-v;
                case 13: return v-w;
            }
        }
#endif
        assert(false);
        itkGenericExceptionMacro("Riemannian3DNorm<TC,TSO>::Stencil::operator[] Error : Index " << i << " exceeds stencil size 14.\n")
    }
    
    // ************** end of stencil substructure ************
    
    template<typename TC,typename TSO>
    typename Riemannian3DNorm<TC,TSO>::OffsetType 
    Riemannian3DNorm<TC,TSO>::ShortOffset_to_Offset(const ShortOffsetType &u)
    {
        OffsetType U; 
        for(int i=0; i<(int)Dimension; ++i) U[i]=u[i];
        return U;
    }
    
    
    template<typename TC,typename TSO> template<typename AcceptedImageType>
    TC Riemannian3DNorm<TC,TSO>::Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                                          const AcceptedImageType * b, unsigned int y_i,
                                          const Stencil &stencil,
                                          IndicesType &gradient_indices, WeightsType &gradient_weights) const
    {
        
        const ShortOffsetType Iy=stencil[y_i];
        const VectorType Vy(Iy);
        const IndexType y=x+ShortOffset_to_Offset(Iy);
        const ValueType Dy=d->GetPixel(y);
        
        ValueType D = Dy+this->Norm(Vy);
        gradient_indices[0]=y_i;
        gradient_indices[1]=this->IndexToIgnore();
        gradient_indices[2]=this->IndexToIgnore();
        
        gradient_weights[0]=1;
        gradient_weights[1]=0;
        //gradient_weights[2] defined implicitly
        
#ifdef SuperBase
        const int s = (stencil.compressedStencilSize == 0) ? 0 : 1;
#else
        //what is the type of stencil ? s=0 : classic orthogonal. s=1 negative scalar product. s=2 non-negative.
        const typename Stencil::TopologyCase s=stencil.topologyCase;
#endif
        const unsigned int nN=nNeigh[s][y_i];
        const unsigned int * N = Neigh[s][y_i];
        
        static const int max_nN = 6;
        bool B[max_nN];
//        ShortOffsetType I[max_nN];
        VectorType V[max_nN];
        ValueType d0[max_nN]; 
        for(int i=0; i<(int)nN; ++i){
            const unsigned int z_i = N[i];
            const ShortOffsetType Iz=stencil[z_i];
            V[i]=VectorType(Iz);
            const IndexType z=x+ShortOffset_to_Offset(Iz);
            if(!b->GetBufferedRegion().IsInside(z)) {B[i]=false; continue;}
            if(!b->GetPixel(z)) {B[i]=false; continue;}  //Point z must be accepted. 
            B[i]=true;
            d0[i]=d->GetPixel(z); //Since z is accepted, Dz < infinity
        }
        
        // Minimizations associated to edges
        for(int i=0; i<(int)nN; ++i){
            if(!B[i]) continue;
            const VectorType Vz=V[i];
            const ValueType Dz=d0[i];
            const int z_i = N[i]; //position in the stencil
            //quadratic system to solve...
            
            typedef Riemannian2DNorm<ValueType,ShortOffsetValueType> Riemannian2DNorm;

            Riemannian2DNorm M(this->SquaredNorm(Vy),this->ScalarProduct(Vy,Vz),this->SquaredNorm(Vz));
            M=M.GetInverse();
            
            typedef typename Riemannian2DNorm::VectorType VectorType2;
            VectorType2 ones; ones.Fill(1);
            VectorType2 dist; dist[0]=Dy; dist[1]=Dz;
            
            const ValueType ones2 = M.SquaredNorm(ones), dist2 = M.SquaredNorm(dist), ones_dist = M.ScalarProduct(ones, dist);
            const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
            if(delta<0) continue;
            const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
            
            if(DTest>=D) continue;
            
            const VectorType2 w=M*(DTest*ones-dist);
            if(w[0]<=0 || w[1]<=0) continue;
            
            D=DTest;
            gradient_weights[0] = w[0]/(w[0]+w[1]);
            gradient_weights[1]=1-gradient_weights[0];
            //gradient_weights[2] implicitly defined
            

//            gradient_indices[0]=y_i; //already done
            gradient_indices[1]=z_i;
            gradient_indices[2]=this->IndexToIgnore();
        }
        
        // Minimizations associated to faces
        for(int i=0; i<(int)nN; ++i){
            const int j=(i+1)%nN;
            if(!B[i] || !B[j]) continue;
            
            const VectorType &Vz = V[i];
            const ValueType Dz = d0[i];
            const int z_i = N[i];
            
            const VectorType &Vt = V[j];
            const ValueType Dt = d0[j];
            const int t_i = N[j];
            
            Riemannian3DNorm M; 
            M(0,0) = this->SquaredNorm(Vy);
            M(1,0) = this->ScalarProduct(Vy,Vz);   M(1,1) = this->SquaredNorm(Vz);
            M(2,0) = this->ScalarProduct(Vy,Vt);   M(2,1) = this->ScalarProduct(Vz,Vt);   M(2,2)=this->SquaredNorm(Vt); 
            
            M=M.GetInverse();
            VectorType ones; ones.Fill(1);
            VectorType dist; dist[0]=Dy; dist[1]=Dz; dist[2]=Dt;
            
            //same equation
            const ValueType ones2 = M.SquaredNorm(ones), dist2 = M.SquaredNorm(dist), ones_dist = M.ScalarProduct(ones, dist);
            const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
            if(delta<0) continue;
            const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
            
            if(DTest>=D) continue;
            
            const VectorType w=M*(DTest*ones-dist);
            if(w[0]<=0 || w[1]<=0 || w[2]<=0) continue;
            
            D=DTest;
            const ValueType wSum = w[0]+w[1]+w[2];
            gradient_weights[0] = w[0]/wSum;
            gradient_weights[1]=w[1]/wSum;
            // gradient_weights[2] implicitly defined
            
            //gradient_indices[0]=y_i; already done
            gradient_indices[1]=z_i;
            gradient_indices[2]=t_i;
        }
        return D;
    }
    
    
    // protected member functions
    
    template<typename TC,typename TSO>
    void Riemannian3DNorm<TC,TSO>::ReducedBasis(VectorType Basis[Dimension]) const
    {
        typedef Riemannian2DNorm<ValueType,ShortOffsetValueType> Riemannian2DNorm;
        typedef typename Riemannian2DNorm::VectorType R2;
        typedef typename Riemannian2DNorm::ShortOffsetType Z2;
        using std::swap;
        using std::min;
        
        VectorType &b0=Basis[0], &b1=Basis[1], &b2=Basis[2];
        ValueType n0=this->SquaredNorm(b0), n1=this->SquaredNorm(b1), n2=this->SquaredNorm(b2); // Squared norms
        
        if(n2<n1){swap(b1,b2);swap(n1,n2);} //Sorting the canonical basis by norm
        if(n1<n0){swap(b0,b1);swap(n0,n1);}
        if(n2<n1){swap(b1,b2);swap(n1,n2);}
        //        if(isDiagonal()) return;
        
        while(true){
            while(true){
                b1 -= ValueType(int(this->ScalarProduct(b0, b1)/n0))*b0;
                n1=this->SquaredNorm(b1);
                //        print_array(coutMath, Basis, Basis+3); cout << endl;        
                if(n0<=n1) break;
                swap(b0,b1); swap(n0,n1);
            }
            Riemannian2DNorm Gram(n0,this->ScalarProduct(b0,b1), n1);
            Gram = Gram.GetInverse();
            R2 scals; scals[0]=this->ScalarProduct(b0,b2); scals[1]=this->ScalarProduct(b1,b2); 
            R2 P = Gram*scals;
            Z2 Q; Q[0]=int(P[0]); Q[1]=int(P[1]);
            P-=R2(Q);
            b2 -= ValueType(Q[0])*b0+ValueType(Q[1])*b1;
            n2=this->SquaredNorm(b2);
            const VectorType b20 = b2-ValueType(sign(P[0]))*b0; const ValueType n20 = this->SquaredNorm(b20);
            const VectorType b21 = b2-ValueType(sign(P[1]))*b1; const ValueType n21 = this->SquaredNorm(b21);
            const double n2p=min(n2,min(n20,n21));
            if(n20==n2p) b2=b20;
            if(n21==n2p) b2=b21;
            n2=n2p;
            //        print_array(coutMath, Basis, Basis+3); cout << endl;        
            if(n1<=n2) break;
            swap(b1,b2); swap(n1,n2);
            if(n0<=n1) continue;
            swap(b0,b1); swap(n0,n1);
        }
    } // reduced basis
    
    template<typename TC,typename TSO>
    void Riemannian3DNorm<TC,TSO>::ObtuseSuperBase(VectorType Basis[Dimension]) const
    {
        ReducedBasis(Basis);
        
        ValueType
        s01 = this->ScalarProduct(Basis[0],Basis[1]),
        s02 = this->ScalarProduct(Basis[0],Basis[2]),
        s12 = this->ScalarProduct(Basis[1],Basis[2]);
        
        using std::swap;
        using std::min;
        const ValueType smin = min(fabs(s01),min(fabs(s02),fabs(s12)));
        if(smin == fabs(s02)){
            swap(s01,s02);
            swap(Basis[1],Basis[2]);
        } else if (smin == fabs(s12)) {
            swap(s01,s12);
            swap(Basis[0],Basis[2]);
        }
        assert( smin == fabs(s01) );
        
        if(s02>0){s02 *= -1; s12*=-1; Basis[2]*=-1;}
        if(s12>0){s12 *= -1; s01*=-1; Basis[1]*=-1;}
        assert( fabs(s01) <= min(-s02,-s12) );
        
        if(s01>0){
            Basis[0] *= -1;
            Basis[2] -= Basis[0];
        }
    } // obtuse superbase
    
    // non function members
    
#ifdef SuperBase
    template<typename TC, typename TSO> const unsigned int Riemannian3DNorm<TC,TSO>::nNeigh[2][14] = {
        {4, 4, 4, 4, 4, 4},
        {6,6,6,6,6,6,6,6,4,4,4,4,4,4}
    };
    
    template<typename TC, typename TSO> const unsigned int Riemannian3DNorm<TC,TSO>::Neigh[2][14][6] = {
        {
            {4, 2, 5, 3}, {4, 2, 5, 3}, {4, 0, 5, 1}, {4, 1, 5, 0},
            {0, 2, 1, 3}, {2, 0, 3, 1}
        }, {
            { 6,11, 5, 9, 7, 8},{ 6,12, 4,10, 7, 8},{ 5,13, 4,10, 7, 9},{ 5,13, 4,12, 6,11},
            {10, 2,13, 3,12, 1},{ 9, 2,13, 3,11, 0},{ 8, 1,12, 3,11, 0},{ 8, 1,10, 2, 9, 0},
            {6,1,7,0},{5,2,7,0},{4,2,7,1},
            {5,3,6,0},{4,3,6,1},{4,3,5,2}
        }
    };
#else
    template<typename TC, typename TSO> const unsigned int Riemannian3DNorm<TC,TSO>::nNeigh[3][14] = {
        {4, 4, 4, 4, 4, 4}, 
        {6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 6, 6}, 
        {6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 6, 6}
    };
    
    template<typename TC, typename TSO> const unsigned int Riemannian3DNorm<TC,TSO>::Neigh[3][14][6] = {
        {
            {4, 2, 5, 3}, {4, 2, 5, 3}, {4, 0, 5, 1}, {4, 1, 5, 0}, 
            {0, 2, 1, 3}, {2, 0, 3, 1}
        }, {
            {4, 10, 2, 9, 13, 7}, {5, 11, 3, 8, 12, 6}, {5, 9, 0, 10, 12, 6}, {4, 8, 1, 11, 13, 7}, 
            {0, 10, 12, 8, 3, 7}, {11, 1, 6, 2, 9, 13}, {5, 2, 12, 1}, {4, 3, 13, 0}, 
            {4, 12, 1, 3}, {5, 13, 0, 2}, {4, 0, 2, 12}, {5, 1, 3, 13}, 
            {4, 10, 2, 6, 1, 8}, {5, 11, 3, 7, 0, 9}
        }, {
            {4, 2, 13, 9, 7, 10}, {5, 3, 12, 8, 6, 11}, {4, 0, 13, 6, 8}, {5, 1, 12, 7, 9}, 
            {0, 2, 8, 12, 10}, {3, 1, 11, 13, 9}, {2, 8, 1, 11, 13}, {3, 9, 0, 10, 12}, 
            {4, 2, 6, 1, 12}, {5, 3, 7, 0, 13}, {4, 12, 7, 0}, {5, 13, 6, 1}, 
            {4, 8, 1, 3, 7, 10}, {5, 9, 0, 2, 6, 11}
        }
    };
#endif
} //namespace itk

#endif
