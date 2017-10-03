//
//  AsymmetricQuadratic2DNorm.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 02/10/2017.
//
//

#ifndef AsymmetricQuadratic2DNorm_h
#define AsymmetricQuadratic2DNorm_h

#include "Riemannian2DNorm.h"

namespace itk
{
    
    /** \class AsymmetricQuadratic2DNorm
     * \brief A type of two dimensional asymmetric norm, also referred to as asymmetric quadratic norm, for use in Fast Marching.
     *
     * This class defines an asymmetric norm of the form \f[F(x) := sqrt{|x|_M^2+<l,x>^2 }.\f]
     For this norm to be well defined, M must be a positive definite symmetric second rank tensor.
     
     
     * This class is templated over the component type. Dimension is fixed to 2.
     * An additional template parameter is defaulted to char.
     Change it to (short int) if you plan to use norms with anisotropy ratio >100.
     
     The class AsymmetricQuadratic2DNorm<TComponent> inherits and has the same memory mapping as
     std::pair< SymmetricSecondRankTensor<TComponent,2>, Vector<TComponent,2> >
     
     This class implements a construction of anisotropic stencils, for use with the AnisotropicFastMarchingImageFilter class, as described in
     * Efficient Fast Marching with Finsler Metrics, J.-M. Mirebeau, 2012
     *
     * \ingroup LevelSetSegmentation
     * \ingroup ITKFastMarching
     */
    
    template<typename TComponent=float, typename TShortOffsetValue=char> class AsymmetricQuadratic2DNorm :
    public std::pair<Riemannian2DNorm<TComponent,TShortOffsetValue>, Vector<TComponent,2> >,
    public AdaptiveStencilRefinement2DNormBase<false, TShortOffsetValue> {
    public:
        static const unsigned int Dimension=2;
        typedef TComponent ValueType;
        typedef TShortOffsetValue ShortOffsetValueType;
        
        typedef Vector<ValueType,Dimension> VectorType;
        typedef CovariantVector<ValueType,Dimension> CovariantVectorType;
        typedef Riemannian2DNorm<ValueType,ShortOffsetValueType> Riemannian2DNorm;
        typedef std::pair<Riemannian2DNorm, VectorType> Superclass1;
        typedef AdaptiveStencilRefinement2DNormBase<false, ShortOffsetValueType> Superclass2;
        
        typedef typename Riemannian2DNorm::SpacingType SpacingType;
        
        
        Riemannian2DNorm & GetM() {return this->first;}
        const Riemannian2DNorm & GetM() const {return this->first;}
        VectorType & GetOmega() {return this->second;}
        const VectorType & GetOmega() const {return this->second;}
        
        ValueType Norm(const VectorType &u) const
        {
            return GetM().Norm(u)-GetOmega()*u;
        }
        
        CovariantVectorType Gradient(const VectorType &u) const;
        
        bool IsAcute(const VectorType &u, const VectorType &v) const
        {
            return Gradient(v)*u>=0 && Gradient(u)*v>=0;
        }
        
        AsymmetricQuadratic2DNorm(){};
        AsymmetricQuadratic2DNorm(const Superclass1 &S):Superclass1(S){};
        AsymmetricQuadratic2DNorm(const AsymmetricQuadratic2DNorm &S, const SpacingType &s)
        {
            GetM() = Riemannian2DNorm(S.GetM(),s);
            for(int i=0; i<(int)Dimension; ++i) GetOmega()[i] = S.GetOmega()[i]*s[i];
        }
        
        
        //Local minimisations : the Hopf-Lax update operator
        typedef Image<ValueType,Dimension> DistanceImageType;
        typedef Index<Dimension> IndexType;
        
        typedef typename Superclass2::ShortOffsetType ShortOffsetType;
        typedef typename Superclass2::CompressedOffsetType CompressedOffsetType;
        typedef typename Superclass2::IndicesType IndicesType;
        typedef typename Superclass2::Stencil Stencil;
        
        typedef typename Riemannian2DNorm::WeightsType WeightsType;
        
        template<typename AcceptedImageType>
        ValueType Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                           const AcceptedImageType * b, unsigned int y_i,
                           const Stencil &stencil,
                           IndicesType &gradient_indices, WeightsType &gradient_weights) const;
        
        void GetCompressedStencil(std::vector<CompressedOffsetType> &l) const
        {
            if(IsInfinite()) return;
            assert(IsDefinite());
            return Superclass2::CompressedStencil(*this,l);
        }
        
        typedef TComponent ComponentType;
        static unsigned int GetNumberOfComponents()
        {
            return Riemannian2DNorm::GetNumberOfComponents() + VectorType::GetNumberOfComponents();
        }
        
        void SetNthComponent(int c, const ComponentType &v)
        {
            assert(0<=c && c <(int)GetNumberOfComponents());
            if(c<(int)Riemannian2DNorm::GetNumberOfComponents()) GetM().SetNthComponent(c,v);
            else GetOmega().SetNthComponent(c-Riemannian2DNorm::GetNumberOfComponents(),v);
        }
        
        ComponentType GetNthComponent(int c)
        {
            assert(0<=c && c<(int)GetNumberOfComponents());
            if(c<(int)Riemannian2DNorm::GetNumberOfComponents()) return GetM().GetNthComponent(c);
            else return GetOmega()[c-Riemannian2DNorm::GetNumberOfComponents()];
        }
        
        void SetIdentity()
        {
            GetM().SetIdentity(); GetOmega().Fill(ValueType(0.));
        }
        
        bool IsInfinite() const {return GetM().IsInfinite();}
        
        bool IsDefinite() const
        {
            return GetM().IsDefinite() && GetM().GetInverse().SquaredNorm(this->second)<1;
        }
        
        AsymmetricQuadratic2DNorm DualNorm() const
        {
            assert(false);
        }
        static AsymmetricQuadratic2DNorm HalfDisk(ScalarType speed, const VectorType & v, ScalarType eps)
        {
            assert(false);
        }
    };
    
    
    
    /**
     * Print content to an ostream
     */
    template<typename TComponent, typename TShortOffsetValue>
    std::ostream &
    operator<<(std::ostream & os, const  AsymmetricQuadratic2DNorm<TComponent,TShortOffsetValue> & s) {
        os << "{" << s.GetM() << "," << s.GetOmega() << "}";
        return os;
    }
    
    
} //namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "AsymmetricQuadratic2DNorm.hxx"
#endif



#endif /* AsymmetricQuadratic2DNorm_h */
