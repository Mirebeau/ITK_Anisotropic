//
//  CoherenceEnhancingDiffusionFilter.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 06/03/2014.
//
//

#ifndef itkDiffusion_CoherenceEnhancingDiffusionFilter_h
#define itkDiffusion_CoherenceEnhancingDiffusionFilter_h

#include "AnisotropicDiffusionLBRImageFilter.h"


namespace itk {
/**
 Implementation of Coherence Enhancing Diffusion (CED), and Edge Enhancing Diffusion (EED), as described by Weickert. CED heuristically smoothes everywhere except accross image contours, while EED smoothes nowhere but tangentially to image contours.
 
 The non-linear diffusion tensor is defined in terms of the structure tensor.
 Denote by \f$\mu_i\f$ the structure tensor eigenvalues, at a given point \f$x\f$, with \f$0\leq i < d\f$. Let also \f$\mu_{\rm min}\f$ and \f$\mu_{\rm max}\f$, be the smallest and largest eigenvalues respectively. The diffusion tensor is defined by the same eigenvectors, but with modified with eigenvalues \f$\lambda_i\f$.
 
 Coherence Enhancing Diffusion : \f$\lambda_i := g(\mu_i - \mu_{\rm min})\f$, where \f$g(s) = 1 - (1-\alpha)*exp(-(\lambda/s)^m)\f$. Note the limit values \f$g(0) = 1\f$, \f$g(\infty) = \alpha\f$.\br

 Edge enhancing diffusion \f$\lambda_i := g(\mu_{\rm max} - \mu_i)\f$, where  \f$g(s) = \alpha + (1-\alpha)*exp(-(\lambda/s)^m)\f$. Note the limit values \f$g(0) = \alpha\f$, \f$g(\infty) = 1\f$.
 */
    template<typename TImage>
    class CoherenceEnhancingDiffusionFilter : public AnisotropicDiffusionLBRImageFilter<TImage> {
    public:
        typedef CoherenceEnhancingDiffusionFilter Self;
        typedef AnisotropicDiffusionLBRImageFilter<TImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        
        /// Method for creation through the object factory.
        itkNewMacro(Self);
        /// Run-time type information (and related methods).
        itkTypeMacro(CoherenceEnhancingDiffusionFilter, Superclass);
        
        typedef typename Superclass::EigenValuesArrayType EigenValuesArrayType;
        virtual EigenValuesArrayType EigenValuesTransform(const EigenValuesArrayType &) const;
        
        typedef typename Superclass::ScalarType ScalarType;
        /// Exponent m involved in the function g defining eigenvalues.
        itkSetMacro(Exponent, ScalarType);
        itkSetMacro(Lambda, ScalarType);
        itkSetMacro(Alpha, ScalarType);
        
        itkGetMacro(Exponent, ScalarType);
        itkGetMacro(Lambda, ScalarType);
        itkGetMacro(Alpha, ScalarType);
        
        enum EnhancementType {Coherence, Edge};
        /// Switch between CED and EED.
        itkSetEnumMacro(Enhancement, EnhancementType);
        itkGetEnumMacro(Enhancement, EnhancementType);
        
    protected:
        ScalarType m_Lambda = 4;
        ScalarType m_Exponent = 2;
        ScalarType m_Alpha = 0.01;
        EnhancementType m_Enhancement = Coherence;
        
        ScalarType g_CED(ScalarType s) const {return s<=0 ? 1 : 1 - (1-m_Alpha)*exp(-1/pow(s/m_Lambda,m_Exponent));}
        ScalarType g_EED(ScalarType s) const {return s<=0 ? m_Alpha : m_Alpha + (1-m_Alpha)*exp(-1/pow(s/m_Lambda,m_Exponent));}
    };
    
    template<typename TI>
    typename CoherenceEnhancingDiffusionFilter<TI>::EigenValuesArrayType
    CoherenceEnhancingDiffusionFilter<TI>::EigenValuesTransform(const EigenValuesArrayType & ev0) const {
        static const int Dimension = Superclass::Dimension;
        ScalarType evMin = ev0[0], evMax = ev0[Dimension-1];
        EigenValuesArrayType ev;
        switch(m_Enhancement){
            case Coherence:
                for(int i=0; i<Dimension; ++i)
                    ev[i] = g_CED(ev0[i]-evMin);
                break;
                
            case Edge:
                for(int i=0; i<Dimension; ++i)
                    ev[i] = g_EED(evMax-ev0[i]);
                break;
            default:
                itkExceptionMacro("Unsupported diffusion type");
        }
        return ev;
    }

    
}
#endif
