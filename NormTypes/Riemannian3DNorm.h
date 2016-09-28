//
//  Riemannian3DNorm.h
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 19/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_Riemannian3DNorm_h
#define ITKFM_Riemannian3DNorm_h

/** \class Riemannian3DNorm
 * \brief A three dimensionnal Symmetric Second Rank Tensor, with added functionnality for Fast Marching.
 *
 * This class inherits from SymmetricSecondRankTensor.
 * Riemannian3DNorm<TComponent> has the same memory mapping as SymmetricSecondRankTensor<TComponent,3>.
 * It implements some additional methods for anisotropic fast marching, using the FM_ASR class.
 * An element M of type Riemannian3DNorm defines the norm x -> sqrt(x*M*x)
 *
 * This class is templated over the component type, like SymmetricSecondRankTensor. Dimension is fixed to 3.
 * An additional template parameter is defaulted to char. Change it to (short int) if you plan to use matrices with condition number >10 000. 
 (condition number is the ratio of the large and small eigenvalues)
 
 This class implements a construction of anisotropic stencils, for use with the FM_ASR class, as described in 
 * Anisotropic Fast Marching using Lattice Basis Reduction, J.-M. Mirebeau, 2012
 *
 * \ingroup LevelSetSegmentation
 * \ingroup ITKFastMarching
 */

#include "Riemannian2DNorm.h"  
#define SuperBase

namespace itk {

    template<typename TComponent=float, typename TShortOffsetValue=char> class Riemannian3DNorm :
    public RiemannianNormBase<TComponent,3> {
    public:
        
        static const unsigned int Dimension=3;
        typedef TComponent ValueType; 
        typedef RiemannianNormBase<ValueType,Dimension> Superclass;
        
        typedef Vector<ValueType,Dimension> VectorType; 
        typedef typename Superclass::SpacingType SpacingType;
        
        typedef TShortOffsetValue ShortOffsetValueType; 
        typedef Vector<ShortOffsetValueType,Dimension> ShortOffsetType;
        typedef ShortOffsetType CompressedOffsetType;
        typedef ShortOffsetType IndicesType;
        typedef typename Superclass::WeightsType WeightsType;
        
        static const ShortOffsetValueType IndexToIgnore(){return NumericTraits<ShortOffsetValueType>::max();}
        
        Riemannian3DNorm(){};
        Riemannian3DNorm(const Superclass &S):Superclass(S){};
        Riemannian3DNorm(const Superclass &M, const SpacingType &s)
        {
            for(int i=0; i<(int)Dimension; ++i) for(int j=i; j<(int)Dimension; ++j) this->coef(i,j)=M(i,j)*s[i]*s[j];
        }
        
        
        void GetCompressedStencil(std::vector<CompressedOffsetType> &l) const;

#ifdef SuperBase
        struct Stencil {
            const ShortOffsetType * const pCompressedStencil;
            const unsigned int compressedStencilSize;
            
            Stencil(const ShortOffsetType * pCompressedStencil_, const unsigned int compressedStencilSize);
            unsigned int size() const;
            ShortOffsetType operator[](int i) const;            
        };
#else
    struct Stencil {
        const ShortOffsetType * const pCompressedStencil;
        const unsigned int compressedStencilSize;
        enum TopologyCase {Six=0,Negative=1,Positive=2};
        const TopologyCase topologyCase;
        
        Stencil(const ShortOffsetType * pCompressedStencil_, const unsigned int compressedStencilSize);
        unsigned int size() const;
        ShortOffsetType operator[](int i) const;
        
        static int Determinant(const ShortOffsetType Basis[3]);
    };
#endif
    
    
        //Local minimisations : the Hopf-Lax update operator
        typedef Image<ValueType,Dimension> DistanceImageType;
//        typedef Image<bool,Dimension> FlagImageType;
        typedef Index<Dimension> IndexType;
        typedef Offset<Dimension> OffsetType;
        
        static OffsetType ShortOffset_to_Offset(const ShortOffsetType &u);
        
        //takes as input : d,x, b,y_i, v,s 
        
        template<typename AcceptedImageType>
        ValueType Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                           const AcceptedImageType * b, unsigned int y_i,
                           const Stencil &stencil,
                           IndicesType &gradient_indices, WeightsType &gradient_weights) const;
        
        
        void ReducedBasis(VectorType Basis[Dimension]) const;
        void ObtuseSuperBase(VectorType SuperBase[Dimension]) const;
        
        bool IsDefinite() const
        {
            ValueType Tr=0,Inv=0;
            for(int i=0; i<3; ++i){
                const int j=(i+1)%3;
                Tr+=this->coef(i,i);
                Inv+=this->coef(i,i)*this->coef(j,j)-this->coef(i,j)*this->coef(j,i);
            }
            return Tr>0 && Inv>0 && this->GetDeterminant()>0;
        }
        
    protected:
        static ShortOffsetValueType sign(ValueType x){return x<0 ? -1 : 1;}        
        
        ///description of the topology of the sub-stencils. Should probably be in the Stencil subclass.
        
#ifdef SuperBase
        static const unsigned int nNeigh[2][14];
        static const unsigned int Neigh[2][14][6];
#else
        static const unsigned int nNeigh[3][14];
        static const unsigned int Neigh[3][14][6];
#endif        
    }; //end of class Riemannian3DNorm
    
} //namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "Riemannian3DNorm.hxx"
#endif

#endif
