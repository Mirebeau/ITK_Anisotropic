//
//  RiemannianNormBase.h
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 05/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_RiemannianNormBase_h
#define ITKFM_RiemannianNormBase_h 


/** 
\class RiemannianNormBase
\brief A Symmetric Second Rank Tensor, with added functionality for geometry.
 
This class inherits from SymmetricSecondRankTensor, and has the same memory mapping. 
It implements some additional functionality, commonly used in geometry. 
 
\ingroup LevelSetSegmentation
\ingroup ITKFastMarching
 */

#include "itkSymmetricSecondRankTensor.h" 


namespace itk {
    
    template<typename TComponent, unsigned int NDimension = 3>
    class RiemannianNormBase :
    public SymmetricSecondRankTensor<TComponent, NDimension>
    {
    public:
        typedef TComponent ValueType;
        static const unsigned int Dimension=NDimension;
        typedef SymmetricSecondRankTensor<TComponent, NDimension> SuperClass;
        typedef Vector<ValueType,Dimension> VectorType;
        typedef Vector<double,Dimension> SpacingType; 
        
        RiemannianNormBase(){};
        RiemannianNormBase(const SuperClass &M):SuperClass(M){};
        RiemannianNormBase(const SuperClass &M, const SpacingType &s){for(int i=0; i<Dimension; ++i) for(int j=0; j<=i; ++j) coef(i,j)=M(i,j)*s[i]*s[j]; }
        
        ///Determinant is only implemented in dimensions 2 and 3.
        ValueType GetDeterminant() const;
        
        ///Inverse is only implemented in dimensions 2 and 3.
        RiemannianNormBase GetInverse() const;
        
        ///Is this a diagonal matrix ?
        bool IsDiagonal() const {for(int i=0; i<Dimension; ++i) for(int j=0; j<i; ++j) if(coef(i,j)!=0) return false; return true;}
        
        ///Matrix Vector Product
        VectorType operator*(const VectorType &u) const {return Product(u);}
        
        ///Scalar product associated to the symmetric matrix M = *this. Returns u*M*v.
        ValueType ScalarProduct(const VectorType &u, const VectorType &v) const {return u*Product(v);}
        
        
        ///Two vectors form an acute angle iff their scalar product is non-negative. Returns the boolean ( this->ScalarProduct(u,v) >= 0 ).
        ValueType IsAcute(const VectorType &u, const VectorType &v) const {return ScalarProduct(u,v)>=0;}
        
        ///Squared Norm of a vector, measured by the symmetric matrix M = *this (assumed here to be positive definite). Returns this->ScalarProduct(u,u).
        ValueType SquaredNorm(const VectorType &u) const {return ScalarProduct(u,u);}
        
        ///Norm of a vector, measured by the symmetric matrix M = *this (assumed here to be positive definite). Returns sqrt(this->SquaredNorm(u)).
        ValueType Norm(const VectorType &u) const {return sqrt(SquaredNorm(u));}
                
        ///Barycentric coordinates
        struct WeightsType : public Vector<ValueType,Dimension-1>
        {
            void SetNullWeights(){this->operator[](0)=-1;}
            bool IsNullWeights() const {return this->operator[](0)==-1;}
            ValueType operator()(int i) const;            
        };
        
        static void SetCanonicalBasis(VectorType[Dimension]);
        
        /// Constructs the rang one matrix "v*Transpose(v)"
        static RiemannianNormBase RankOneTensor(const VectorType & v)
        {
            RiemannianNormBase m;
            for(int i=0; i<Dimension; ++i)
                for(int j=0; j<Dimension; ++j)
                    m(i,j)=v[i]*v[j];
            return m;
        }
        
    protected:
        ValueType &         coef(unsigned int row, unsigned int col)        {return this->operator()(row,col);}
        const ValueType & 	coef(unsigned int row, unsigned int col) const  {return this->operator()(row,col);}
        
        ///Matrix Vector Product, with a vector of the same dimension
        VectorType Product(const VectorType &u) const;
    };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "RiemannianNormBase.hxx"
#endif

#endif
