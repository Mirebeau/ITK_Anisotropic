//
//  Riemannian2DNorm.h
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 17/09/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_Riemannian2DNorm_h
#define ITKFM_Riemannian2DNorm_h 

/** \class Riemannian2DNorm
 * \brief A two dimensionnal Symmetric Second Rank Tensor, with added functionnality for Fast Marching.
 *
 * This class inherits from SymmetricSecondRankTensor.
 * Riemannian2DNorm<TComponent> has the same memory mapping as SymmetricSecondRankTensor<TComponent,2>.
 * It implements some additional methods for anisotropic fast marching, using the AnisotropicFastMarchingImageFilter class.
 * An element M of type Riemannian2DNorm defines the norm x -> sqrt(x*M*x)
 *
 * This class is templated over the ComponentType, like SymmetricSecondRankTensor. Dimension is fixed to 2.
 * An additional template parameter is defaulted to char. Change it to (short int) if you plan to use matrices with condition number >10 000. 
 (condition number is the ratio of the large and small eigenvalues)
 
 This class implements a construction of anisotropic stencils, for use with the AnisotropicFastMarchingImageFilter class, as described in
 * Efficient Fast Marching with Finsler Metrics, J.-M. Mirebeau, 2012
 *
 * \ingroup LevelSetSegmentation
 * \ingroup ITKFastMarching
 */ 

#include "RiemannianNormBase.h"  
#include "AdaptiveStencilRefinement2DNormBase.h" 

namespace itk
{

    template<typename TComponent=float, typename TShortOffsetValue=char> class Riemannian2DNorm :
    public RiemannianNormBase<TComponent,2>, 
    public AdaptiveStencilRefinement2DNormBase<true, TShortOffsetValue>
    {
    public:
        static const unsigned int Dimension=2;
        typedef TComponent ValueType; 
        typedef RiemannianNormBase<ValueType,Dimension> Superclass1;
        typedef AdaptiveStencilRefinement2DNormBase<true, TShortOffsetValue> Superclass2;
        
        typedef Vector<ValueType,Dimension> VectorType; 
        typedef typename Superclass1::SpacingType SpacingType; 
        
                
        Riemannian2DNorm(){};
        Riemannian2DNorm(ValueType xx, ValueType xy, ValueType yy){this->coef(0,0)=xx; this->coef(1,0)=xy; this->coef(1,1)=yy;}

        ///Pass through constructors
        Riemannian2DNorm(const Superclass1 &S, const SpacingType &s): Superclass1(S,s){};
        Riemannian2DNorm(const Superclass1 & S): Superclass1(S){};
        
        
        //Local minimisations : the Hopf-Lax update operator
        typedef Image<ValueType,Dimension> DistanceImageType;
//        typedef Image<bool,Dimension> FlagImageType;
        typedef Index<Dimension> IndexType;
        
        typedef typename Superclass2::ShortOffsetType ShortOffsetType;
        typedef typename Superclass2::CompressedOffsetType CompressedOffsetType;
        typedef typename Superclass2::IndicesType IndicesType;
        typedef typename Superclass2::Stencil Stencil;

        typedef typename Superclass1::WeightsType WeightsType;
        
        /**
         The Hopf-Lax update operator is used to propagate the distance estimate, from the seeds, 
         as the front advances in the fast marching algorithm.
         
         Input  
         d : image of the distances, from the seeds, that have been computed.
         x : point at which a new distance estimate is required.
         b : boolean image, indicating at which points distance estimates are final.
         y_i : index, in the direct stencil of x, of the front point which was just accepted 
               (by the fast marching algorithm) thus trigerring the evaluation of this Hopf-Lax update.
         pCompressedStencil, compressedStencilSize : data required to access the direct stencil of x.
         
         gradient_indices, gradient_weights : returns by reference, the virtual point, on a face of the 
           stencil of x, where the minimum in the Hopf-Lax update operator was reached. This is used to compute 
         the gradient of the distance image, and the geodesics.
         
         If the return value (the proposed distance estimate at x from the stencil position y_i) is smaller
         than the previous estimate at x, then a user executing the fast marching algorithm should update
         the image d at x wiTCth this value.
         */
        
        template<typename AcceptedImageType>
        ValueType Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                           const AcceptedImageType * b, unsigned int y_i,
                           const Stencil & stencil,
                           IndicesType &gradient_indices, WeightsType &gradient_weights) const;

        void GetCompressedStencil(std::vector<CompressedOffsetType> &l) const {
            assert(IsDefinite());
            return Superclass2::CompressedStencil(*this,l);
        }
        
        void ReducedBasis(VectorType Basis[Dimension]) const;
        void ObtuseSuperBase(VectorType SuperBase[Dimension]) const;
        
        bool IsDefinite() const
        {
            const ValueType Tr = this->coef(0,0) + this->coef(1,1);
            return Tr>0 && this->GetDeterminant()>0;
        }
    };

} //namespace itk 

#ifndef ITK_MANUAL_INSTANTIATION
#include "Riemannian2DNorm.hxx"
#endif

#endif
