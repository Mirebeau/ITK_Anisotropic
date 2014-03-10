//
//  RiemannianNormBase.hxx
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 05/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_RiemannianNormBase_hxx
#define ITKFM_RiemannianNormBase_hxx

namespace itk {
    
    //Matrix Vector product.
    template<typename TC, unsigned int ND>  
    typename RiemannianNormBase<TC,ND>::VectorType 
    RiemannianNormBase<TC,ND>::Product(const VectorType &u) const
    {
        VectorType V;
        for(int i=0; i<Dimension; ++i){
            V[i]=0;
            for(int j=0; j<Dimension; ++j)
                V[i]+=coef(i,j)*u[j];
        }
        return V;
    }
    
    
    //Determinant, in dimension 2 and 3.
    template<typename TC, unsigned int ND> 
    TC RiemannianNormBase<TC,ND>::GetDeterminant() const
    {
        ValueType d;
        switch (Dimension) { 
            case 2:
                d=coef(0,0)*coef(1,1)-coef(1,0)*coef(1,0);
                break;
                
            case 3:
                d=0;
                for(int i=0; i<3; ++i) d+=coef(i,0)*(coef((i+1)%3,1)*coef((i+2)%3,2)-coef((i+2)%3,1)*coef((i+1)%3,2)); 
                return d;
                break;
                
            default:
                assert(false);
        }
        return d;
    }


    //Inverse, in dimension 2 and 3.
    template<typename TC, unsigned int ND> 
    RiemannianNormBase<TC,ND> RiemannianNormBase<TC,ND>::GetInverse() const
    {
        const ValueType d = GetDeterminant();
        assert(d!=ValueType(0));
        RiemannianNormBase inv;
        switch (Dimension) {
            case 2:
                inv(0,0)=coef(1,1)/d; 
                inv(1,0) = -coef(0,1)/d;
                inv(1,1) = coef(0,0)/d;
                break;
                
            case 3:
                for(int i=0; i<Dimension; ++i)
                    for(int j=i; j<Dimension; ++j)
                        inv(i,j) = (coef((i+1)%3,(j+1)%3)*coef((i+2)%3,(j+2)%3)-coef((i+1)%3,(j+2)%3)*coef((i+2)%3,(j+1)%3))/d;
                break;
                
            default:
                assert(false);
        }
        return inv;
    }

    template<typename TC, unsigned int ND>
    TC RiemannianNormBase<TC,ND>::WeightsType::operator()(int i) const
    {
        assert(0<=i && i<ND); 
        assert(!IsNullWeights());
        if(i<ND-1) return this->operator[](i);
        ValueType sum=this->operator[](0);
        for(int i=1; i<ND-1; ++i) sum+=this->operator[](i);
        return 1-sum;
    }
    
    template<typename TC, unsigned int ND>
    void RiemannianNormBase<TC,ND>::SetCanonicalBasis(VectorType Basis[Dimension])
    {
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<Dimension; ++j)
                Basis[i][j] = (i==j);
    }
    
} //namespace itk
#endif
