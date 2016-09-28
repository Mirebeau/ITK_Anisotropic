//
//  ExtendedNorm.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 07/01/13.
//
//

#ifndef ITKFM_ExtendedNorm_h
#define ITKFM_ExtendedNorm_h

#include "Riemannian2DNorm.h"

namespace itk {
    /**
    Implementation of a \f$d+1\f$-dimensional norms, defined in terms of an \f$d\f$-dimensional norm and scalar. For use with AnisotropicFastMarchingImageFilter.
     */
    
    template<typename TPrimaryNorm>
    class ExtendedNorm :
    public std::pair<
    TPrimaryNorm,
    typename TPrimaryNorm::ValueType
    >
    {
    public:
        typedef TPrimaryNorm PrimaryNorm;
        typedef typename PrimaryNorm::ValueType ValueType;
        typedef typename std::pair<PrimaryNorm,ValueType> Superclass;
        typedef ExtendedNorm Self;
        
        static const unsigned int PrimaryDimension = PrimaryNorm::Dimension;
        typedef Vector<double,PrimaryDimension> PrimarySpacingType;
        
        
        static const unsigned int Dimension = PrimaryDimension + 1;
        typedef RiemannianNormBase<ValueType,Dimension> Traits1;
        typedef typename Traits1::VectorType VectorType;
        typedef typename Traits1::SpacingType SpacingType;
        typedef typename Traits1::WeightsType WeightsType;
        
        typedef typename PrimaryNorm::CompressedOffsetType CompressedOffsetType;
        typedef typename PrimaryNorm::ShortOffsetType PrimaryShortOffsetType;
        typedef typename PrimaryNorm::ShortOffsetValueType ShortOffsetValueType;
        typedef typename PrimaryNorm::VectorType PrimaryVectorType;
        typedef Vector<ShortOffsetValueType,Dimension> ShortOffsetType;
        typedef ShortOffsetType IndicesType;
        
        static const ShortOffsetValueType IndexToIgnore(){return PrimaryNorm::IndexToIgnore();}


        typedef Image<ValueType,Dimension> DistanceImageType;
        typedef Index<Dimension> IndexType;
        typedef Offset<Dimension> OffsetType;
        static OffsetType ShortOffset_to_Offset(const ShortOffsetType &u){
            OffsetType U; for(int i=0; i<(int)Dimension; ++i) U[i]=u[i]; return U;}
        static OffsetType PrimaryShortOffset_to_Offset(const PrimaryShortOffsetType &u){
            OffsetType U; for(int i=0; i<(int)PrimaryDimension; ++i) U[i]=u[i]; U[PrimaryDimension]=0; return U;}
        
        PrimaryNorm & GetPrimaryNorm()              {return this->first;}
        const PrimaryNorm & GetPrimaryNorm() const  {return this->first;}
        ValueType & GetScalar()                     {return this->second;}
        const ValueType & GetScalar() const         {return this->second;}
        
        
        ExtendedNorm(){};
        ExtendedNorm(const Superclass &S):Superclass(S){};
        ExtendedNorm(const Superclass &S, const SpacingType &s):
        Superclass(std::make_pair(PrimaryNorm(S.first, PrimarySpacingType(s.GetDataPointer() ) ), S.second*vnl_math_sqr(s[PrimaryDimension])) )
        {};
        
        struct Stencil {
            const typename PrimaryNorm::Stencil stencil;
            typedef typename PrimaryNorm::ShortOffsetType PrimaryShortOffsetType;
            typedef Self Parent;
            
            unsigned int size() const {return stencil.size()+2;}
            
            ShortOffsetType operator[](int i) const {
                ShortOffsetType U; U.Fill(0);
                const int j = i-stencil.size();
                if(j<0){
                    PrimaryShortOffsetType u = stencil[i];
                    for(int i=0; i<(int)PrimaryDimension; ++i) U[i]=u[i];
                } else if(j==0){
                    U[PrimaryDimension] = 1;
                } else {
                    U[PrimaryDimension]=-1;
                }
                return U;
            }
            
            Stencil(const CompressedOffsetType* pCompressedStencil, unsigned int compressedStencilSize):stencil(pCompressedStencil, compressedStencilSize){};
        };
        
        
        void GetCompressedStencil(std::vector<CompressedOffsetType> &l) const {return GetPrimaryNorm().GetCompressedStencil(l);}

        
        template<typename AcceptedImageType>
        ValueType Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                           const AcceptedImageType * b, unsigned int y_i,
                           const Stencil & stencil,
                           IndicesType &gradient_indices, WeightsType &gradient_weights) const;
        
        typedef ValueType ComponentType;
        static unsigned int GetNumberOfComponents(){return PrimaryNorm::GetNumberOfComponents()+1;}
        
        void SetNthComponent(int c, const ComponentType &v)
        {
            assert(0<=c && c<(int)GetNumberOfComponents());
            if(c<(int)PrimaryNorm::GetNumberOfComponents()) GetPrimaryNorm().SetNthComponent(c,v);
            else GetScalar() = v;
        }
        
        ComponentType GetNthComponent(int c)
        {
            assert(0<=c && c<(int)GetNumberOfComponents());
            if(c<(int)PrimaryNorm::GetNumberOfComponents()) return GetPrimaryNorm().GetNthComponent(c);
            else return GetScalar();
        }
        
        void SetIdentity(){GetPrimaryNorm().SetIdentity(); GetScalar() = 1;}
        
        bool IsDefinite() const {return GetPrimaryNorm().IsDefinite() && GetScalar()>0;}
        bool IsInfinite() const {return GetPrimaryNorm().IsInfinite();}
        
#pragma remark("Remove when old template mess is eliminated.")
/*
    protected:
        // I need to specialize here, depending on TPrimaryNorm -> Gets rather ugly.
        typedef Riemannian2DNorm<ValueType,ShortOffsetValueType>    Riemannian2DNorm;
        typedef Finsler2DNorm<ValueType,ShortOffsetValueType>       Finsler2DNorm;
        typedef Riemannian3DNorm<ValueType,ShortOffsetValueType>    Riemannian3DNorm;

        struct HL;
 */
        
    }; // End of class ExtendedNorm

    
    template<typename TPrimaryNorm>
    std::ostream &
    operator<<(std::ostream & os, const  ExtendedNorm<TPrimaryNorm> & s){
        os << "{" << s.GetPrimaryNorm() << "," << s.GetScalar() << "}";
        return os;
        }
    
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ExtendedNorm.hxx"
#endif

#endif
