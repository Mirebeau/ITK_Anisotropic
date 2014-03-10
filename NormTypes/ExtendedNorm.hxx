//
//  ExtendedNorm.hxx
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 07/01/13.
//
//

#ifndef ITKFM_ExtendedNorm_hxx
#define ITKFM_ExtendedNorm_hxx

#include "Riemannian2DNorm.h"
#include "Finsler2DNorm.h"
#include "Riemannian3DNorm.h"

// The Hopf-Lax update operator unfortunately must be specialized.
// These specializations are put in an anonymous namespace.


namespace itk {
    
    namespace {
        template<typename TC, typename TSO, typename AcceptedImageType>
        typename ExtendedNorm<Finsler2DNorm<TC,TSO> >::ValueType
        Specialized_Hopf_Lax(
                             const ExtendedNorm<Finsler2DNorm<TC,TSO> > *pNorm,
                             const typename ExtendedNorm<Finsler2DNorm<TC,TSO> >::DistanceImageType * d,
                             const typename ExtendedNorm<Finsler2DNorm<TC,TSO> >::IndexType & x,
                             const AcceptedImageType * b,
                             unsigned int y_i,
                             const typename ExtendedNorm<Finsler2DNorm<TC,TSO> >::Stencil & stencil,
                             typename ExtendedNorm<Finsler2DNorm<TC,TSO> >::IndicesType & gradient_indices,
                             typename ExtendedNorm<Finsler2DNorm<TC,TSO> >::WeightsType & gradient_weights
                             ){
            typedef Finsler2DNorm<TC,TSO> Finsler2DNorm;
            typedef Riemannian2DNorm<TC,TSO> Riemannian2DNorm;
            typedef ExtendedNorm<Finsler2DNorm> ExtendedNorm;
            
            typedef typename ExtendedNorm::PrimaryVectorType PrimaryVectorType;
            typedef typename ExtendedNorm::ValueType ValueType;
            typedef typename ExtendedNorm::IndexType IndexType;
            typedef typename ExtendedNorm::OffsetType OffsetType;
            typedef typename ExtendedNorm::PrimaryNorm PrimaryNorm;
            typedef typename ExtendedNorm::IndicesType IndicesType;
            typedef typename ExtendedNorm::ShortOffsetType ShortOffsetType;
            
            const typename Finsler2DNorm::Stencil & primaryStencil = stencil.stencil;
            const unsigned int primarySize = primaryStencil.size();
            typedef typename Finsler2DNorm::ShortOffsetType PrimaryShortOffsetType;
            
            const Finsler2DNorm & primaryNorm = pNorm->GetPrimaryNorm();
            const Riemannian2DNorm & normM = primaryNorm.GetM();
            const PrimaryVectorType & normOmega = primaryNorm.GetOmega();
            const ValueType H = pNorm->GetScalar();
            
            ValueType D; // contains the estimated distance
            
            if(y_i<primarySize){
                
                const PrimaryShortOffsetType BSOy=primaryStencil[y_i];
                const OffsetType Oy = ExtendedNorm::PrimaryShortOffset_to_Offset(BSOy);
                const PrimaryVectorType Vy(BSOy);
                
                const IndexType y=x+Oy;
                const ValueType Dy=d->GetPixel(y);
                
                // Estimate from the point y alone
                D = Dy+primaryNorm.Norm(Vy);
                gradient_indices[0]=y_i;
                gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                gradient_indices[2]=PrimaryNorm::IndexToIgnore();
                
                gradient_weights[0]=1;
                gradient_weights[1]=0;
                //gradient_weights[Dimension-1] implicitly defined
                
                
                // Finding out who is the closest point p in parameter space
                // Storing the associated distance and index in Dp and IndexP
                ValueType Dp = NumericTraits<ValueType>::max();
                typename IndicesType::ValueType IndexP;
                
                for(int e=primarySize; e<stencil.size(); ++e){
                    const ShortOffsetType BSOp = stencil[e];
                    OffsetType Op = pNorm->ShortOffset_to_Offset(BSOp);
                    
                    const IndexType p=x+Op;
                    if(!b->GetBufferedRegion().IsInside(p)) continue;
                    if(!b->GetPixel(p)) continue;  //Point p must be accepted.
                    
                    const ValueType DpTest=d->GetPixel(p);
                    if(DpTest<Dp){
                        Dp=DpTest;
                        IndexP=e;
                    }
                } // for e in last two stencil elements
                
                
                
                // Estimate from the segment [y,p] joining y to p
                if(Dp<NumericTraits<ValueType>::max())
                    do {
                        const ValueType
                        N2y = normM.SquaredNorm(Vy),
                        N2p = H;
                        
                        const ValueType Dy_ = Dy - normOmega*Vy;
                        
                        //solving quadratic system a D^2 - 2 b D + c = 0;
                        const ValueType
                        a = 1./N2y+1./N2p,
                        b = Dy_/N2y + Dp/N2p,
                        c = vnl_math_sqr(Dy_)/N2y + vnl_math_sqr(Dp)/N2p - 1.;
                        
                        const ValueType delta = vnl_math_sqr(b) - a*c;
                        if(delta < 0) continue;
                        const ValueType DTest = ( b+sqrt(delta) ) / a;
                        if(DTest >= D) continue;
                        
                        PrimaryVectorType w;
                        w[0]=(DTest-Dy_)/N2y;
                        w[1]=(DTest-Dp)/N2p;
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[2]=IndexP;
                        
                        gradient_weights[0]=w[0]/(w[0]+w[1]);
                        gradient_weights[1]=0;
                        // gradient_weights[2] // automatically defined by 1- gradient_weights[0] -gradient_weights[0].
                        
                    } while (false);
                
                
                // Estimate from the segments [y,z] joining y to a planar neighbor z, and from the corresponding face [y,z,p] including p
                for(int e=-1; e<=1; e+=2){
                    int z_i = int(y_i)+e;
                    z_i = (z_i+primarySize)%primarySize;
                    
                    const PrimaryShortOffsetType BSOz=primaryStencil[z_i];
                    OffsetType Oz = pNorm->PrimaryShortOffset_to_Offset(BSOz);
                    PrimaryVectorType Vz(BSOz);
                    
                    const IndexType z=x+Oz;
                    if(!b->GetBufferedRegion().IsInside(z)) continue;
                    if(!b->GetPixel(z)) continue;  //Point z must be accepted.
                    const ValueType Dz=d->GetPixel(z); //Since z is accepted, Dz < infinity
                    
                    
                    
                    // estimate from the segment [y,z]
                    //quadratic system to solve...
                    
                    const Riemannian2DNorm M(normM.SquaredNorm(Vy),normM.ScalarProduct(Vy,Vz),normM.SquaredNorm(Vz));
                    const Riemannian2DNorm N=M.GetInverse();
                    
                    PrimaryVectorType ones; ones.Fill(1);
                    PrimaryVectorType dist; dist[0]=Dy - normOmega * Vy; dist[1]=Dz - normOmega * Vz;
                    
                    ValueType
                    ones2 = N.SquaredNorm(ones),
                    dist2 = N.SquaredNorm(dist),
                    ones_dist = N.ScalarProduct(ones, dist);
                    ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    do {
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const PrimaryVectorType w=N*(DTest*ones-dist);
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        gradient_weights[0]= w[0]/(w[0]+w[1]);
                        gradient_weights[1] = 1-gradient_weights[0];
                        //gradient_weights[Dimension-1] is implicitly defined
                        
                        //gradient_indices[0]=y_i; //aready set
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=PrimaryNorm::IndexToIgnore();
                    } while(false);
                    
                    
                    // estimate from the face [y,z,p]
                    // the quadratic system must be slightly updated to take p into account.
                    
                    if(Dp<NumericTraits<ValueType>::max()){
                        const ValueType N2p=H;
                        ones2+=1./N2p;
                        dist2+=vnl_math_sqr(Dp)/N2p;
                        ones_dist+=Dp/N2p;
                        
                        delta = ones_dist*ones_dist - ones2*(dist2-1);
                        
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const PrimaryVectorType w=N*(DTest*ones-dist);
                        const ValueType w2 = (DTest-Dp)/N2p;
                        
                        if(w[0]<=0 || w[1]<=0 || w2<=0) continue;
                        
                        D=DTest;
                        const ValueType wSum = w[0]+w[1]+w2;
                        gradient_weights[0]= w[0]/wSum;
                        gradient_weights[1] = w[1]/wSum;
                        //gradient_weights[Dimension-1] is implicitly defined
                        
                        //gradient_indiced[0]=y_i; //already set
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=IndexP;
                    }
                }
                
            } else { // y_i>=primarySize
                
                // reference point, of index y_i, is in a different parameter plane.
                // we call it p to be consistent with the first case.
                
                const unsigned int p_i = y_i;
                ShortOffsetType SOp = stencil[p_i];
                OffsetType Op = pNorm->ShortOffset_to_Offset(SOp);
                const IndexType p = x+Op;
                const ValueType Dp = d->GetPixel(p);
                const ValueType N2p = H;
                
                //direct test from p
                
                D = Dp+sqrt(N2p);
                gradient_indices[0]=p_i; gradient_indices[1]=0; gradient_indices[2]=0; //Arbitrary value at position 1 and 2, within GetStencilSize.
                gradient_weights[0]=1; gradient_weights[1]=0;
                
                // circle around points and segments.
                for(int y_i=0; y_i<primarySize; ++y_i){
                    //hides the previous instance of y_i on purpose (
                    PrimaryShortOffsetType BSOy = primaryStencil[y_i];
                    OffsetType Oy = pNorm->PrimaryShortOffset_to_Offset(BSOy);
                    IndexType y = x+Oy;
                    PrimaryVectorType Vy(BSOy);
                    
                    if(!b->GetBufferedRegion().IsInside(y)) continue;
                    if(!b->GetPixel(y)) continue;  //Point y must be accepted.
                    const ValueType Dy=d->GetPixel(y);
                    
                    const ValueType N2y = normM.SquaredNorm(Vy);
                    
                    
                    // Estimate from the segment [y,p] joining y to p
                    do {
                        const ValueType Dy_ = Dy - normOmega * Vy;
                        
                        //solving quadratic system a D^2 - 2 b D + c = 0;
                        const ValueType
                        a = 1./N2y+1./N2p,
                        b = Dy_/N2y + Dp/N2p,
                        c = vnl_math_sqr(Dy_)/N2y + vnl_math_sqr(Dp)/N2p - 1.;
                        
                        const ValueType delta = vnl_math_sqr(b) - a*c;
                        if(delta < 0) continue;
                        const ValueType DTest = ( b+sqrt(delta) ) / a;
                        if(DTest >= D) continue;
                        
                        PrimaryVectorType w;
                        w[0]=(DTest-Dy_)/N2y;
                        w[1]=(DTest-Dp)/N2p;
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        gradient_indices[0]=y_i;
                        gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[2]=p_i;
                        
                        gradient_weights[0]=w[0]/(w[0]+w[1]);
                        gradient_weights[1]=0;
                        // gradient_weights[2] // automatically defined by 1- gradient_weights[0] -gradient_weights[1].
                        
                    } while (false);
                    
                    // Estimate from the face [y,z,p]
                    
                    const int z_i= (y_i+1)%primarySize;
                    PrimaryShortOffsetType BSOz = primaryStencil[z_i];
                    OffsetType Oz = pNorm->PrimaryShortOffset_to_Offset(BSOz);
                    IndexType z = x+Oz;
                    PrimaryVectorType Vz(BSOz);
                    
                    if(!b->GetBufferedRegion().IsInside(z)) continue;
                    if(!b->GetPixel(z)) continue;  //Point z must be accepted.
                    const ValueType Dz=d->GetPixel(z);
                    
                    const Riemannian2DNorm M(normM.SquaredNorm(Vy),normM.ScalarProduct(Vy,Vz),normM.SquaredNorm(Vz));
                    const Riemannian2DNorm N=M.GetInverse();
                    
                    PrimaryVectorType ones; ones.Fill(1);
                    PrimaryVectorType dist; dist[0]=Dy - normOmega * Vy; dist[1]=Dz - normOmega * Vz;
                    
                    ValueType
                    ones2 = N.SquaredNorm(ones),
                    dist2 = N.SquaredNorm(dist),
                    ones_dist = N.ScalarProduct(ones, dist);
                    
                    // estimate from the face [y,z,p]
                    // the quadratic system must be slightly updated to take p into account.
                    
                    ones2+=1./N2p;
                    dist2+=vnl_math_sqr(Dp)/N2p;
                    ones_dist+=Dp/N2p;
                    
                    const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    
                    if(delta<0) continue;
                    const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                    
                    if(DTest>=D) continue;
                    
                    const PrimaryVectorType w=N*(DTest*ones-dist);
                    const ValueType w2 = (DTest-Dp)/N2p;
                    
                    if(w[0]<=0 || w[1]<=0 || w2<=0) continue;
                    
                    D=DTest;
                    const ValueType wSum = w[0]+w[1]+w2;
                    gradient_weights[0]= w[0]/wSum;
                    gradient_weights[1] = w[1]/wSum;
                    //gradient_weights[Dimension-1] is implicitly defined
                    
                    gradient_indices[0] = y_i;
                    gradient_indices[1] = z_i;
                    gradient_indices[2] = p_i;
                    
                } // for y_i (renamed..)
                
            } // if y_i
            
            return D;
        }
        
        // ************************************** Riemannian2DNorm *******************************************
        
        template<typename TC, typename TSO, typename AcceptedImageType>
        typename ExtendedNorm<Riemannian2DNorm<TC,TSO> >::ValueType
        Specialized_Hopf_Lax(
                             const ExtendedNorm<Riemannian2DNorm<TC,TSO> > *pNorm,
                             const typename ExtendedNorm<Riemannian2DNorm<TC,TSO> >::DistanceImageType * d,
                             const typename ExtendedNorm<Riemannian2DNorm<TC,TSO> >::IndexType & x,
                             const AcceptedImageType * b,
                             unsigned int y_i,
                             const typename ExtendedNorm<Riemannian2DNorm<TC,TSO> >::Stencil & stencil,
                             typename ExtendedNorm<Riemannian2DNorm<TC,TSO> >::IndicesType & gradient_indices,
                             typename ExtendedNorm<Riemannian2DNorm<TC,TSO> >::WeightsType & gradient_weights
                             ){
            typedef Riemannian2DNorm<TC,TSO> Riemannian2DNorm;
            typedef ExtendedNorm<Riemannian2DNorm> ExtendedNorm;
            
            typedef typename ExtendedNorm::PrimaryVectorType PrimaryVectorType;
            typedef typename ExtendedNorm::ValueType ValueType;
            typedef typename ExtendedNorm::IndexType IndexType;
            typedef typename ExtendedNorm::OffsetType OffsetType;
            typedef typename ExtendedNorm::PrimaryNorm PrimaryNorm;
            typedef typename ExtendedNorm::IndicesType IndicesType;
            typedef typename ExtendedNorm::ShortOffsetType ShortOffsetType;
            
            const unsigned int Dimension = ExtendedNorm::Dimension;
            
            const typename Riemannian2DNorm::Stencil & primaryStencil = stencil.stencil;
            const unsigned int primarySize = primaryStencil.size();
            typedef typename Riemannian2DNorm::ShortOffsetType PrimaryShortOffsetType;
            
            
            const Riemannian2DNorm & primaryNorm = pNorm->GetPrimaryNorm();
            const ValueType H = pNorm->GetScalar();
            
            ValueType D; // contains the estimated distance
            
            if(y_i<primarySize){
                
                const PrimaryShortOffsetType BSOy=primaryStencil[y_i];
                const OffsetType Oy = pNorm->PrimaryShortOffset_to_Offset(BSOy);
                const PrimaryVectorType Vy(BSOy);
                
                const IndexType y=x+Oy;
                const ValueType Dy=d->GetPixel(y);
                
                // Estimate from the point y alone
                D = Dy+primaryNorm.Norm(Vy);
                gradient_indices[0]=y_i;
                for(int i=1; i<Dimension; ++i) gradient_indices[i] = PrimaryNorm::IndexToIgnore();
                
                gradient_weights[0]=1;
                for(int i=1; i<Dimension-1; ++i) gradient_weights[i]=0;
                
                
                // Finding out who is the closest point p in parameter space
                // Storing the associated distance and index in Dp and IndexP
                ValueType Dp = NumericTraits<ValueType>::max();
                typename IndicesType::ValueType IndexP;
                
                for(int e=primarySize; e<stencil.size(); ++e){
                    const ShortOffsetType BSOp = stencil[e];
                    OffsetType Op = pNorm->ShortOffset_to_Offset(BSOp);
                    
                    const IndexType p=x+Op;
                    if(!b->GetBufferedRegion().IsInside(p)) continue;
                    if(!b->GetPixel(p)) continue;  //Point p must be accepted.
                    
                    const ValueType DpTest=d->GetPixel(p);
                    if(DpTest<Dp){
                        Dp=DpTest;
                        IndexP=e;
                    }
                } // for e in last two stencil elements
                
                
                
                // Estimate from the segment [y,p] joining y to p
                if(Dp<NumericTraits<ValueType>::max())
                    do {
                        const ValueType
                        N2y = primaryNorm.SquaredNorm(Vy),
                        N2p = H;
                        
                        //solving quadratic system a D^2 - 2 b D + c = 0;
                        const ValueType
                        a = 1./N2y+1./N2p,
                        b = Dy/N2y + Dp/N2p,
                        c = vnl_math_sqr(Dy)/N2y + vnl_math_sqr(Dp)/N2p - 1.;
                        
                        const ValueType delta = vnl_math_sqr(b) - a*c;
                        if(delta < 0) continue;
                        const ValueType DTest = ( b+sqrt(delta) ) / a;
                        if(DTest >= D) continue;
                        
                        PrimaryVectorType w;
                        w[0]=(DTest-Dy)/N2y;
                        w[1]=(DTest-Dp)/N2p;
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[2]=IndexP;
                        
                        gradient_weights[0]=w[0]/(w[0]+w[1]);
                        gradient_weights[1]=0;
                        // gradient_weights[2] // automatically defined by 1- gradient_weights[0] -gradient_weights[1].
                        
                    } while (false);
                
                
                // Estimate from the segments [y,z] joining y to a planar neighbor z, and from the corresponding face [y,z,p] including p
                for(int e=-1; e<=1; e+=2){
                    int z_i = int(y_i)+e;
                    z_i = (z_i+primarySize)%primarySize;
                    
                    const PrimaryShortOffsetType BSOz=primaryStencil[z_i];
                    OffsetType Oz = pNorm->PrimaryShortOffset_to_Offset(BSOz);
                    PrimaryVectorType Vz(BSOz);
                    
                    const IndexType z=x+Oz;
                    if(!b->GetBufferedRegion().IsInside(z)) continue;
                    if(!b->GetPixel(z)) continue;  //Point z must be accepted.
                    const ValueType Dz=d->GetPixel(z); //Since z is accepted, Dz < infinity
                    
                    
                    
                    // estimate from the segment [y,z]
                    //quadratic system to solve...
                    
                    const Riemannian2DNorm M(primaryNorm.SquaredNorm(Vy),primaryNorm.ScalarProduct(Vy,Vz),primaryNorm.SquaredNorm(Vz));
                    const Riemannian2DNorm N=M.GetInverse();
                    
                    PrimaryVectorType ones; ones.Fill(1);
                    PrimaryVectorType dist; dist[0]=Dy; dist[1]=Dz;
                    
                    ValueType
                    ones2 = N.SquaredNorm(ones),
                    dist2 = N.SquaredNorm(dist),
                    ones_dist = N.ScalarProduct(ones, dist);
                    ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    do {
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const PrimaryVectorType w=N*(DTest*ones-dist);
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        gradient_weights[0]= w[0]/(w[0]+w[1]);
                        gradient_weights[1] = 1-gradient_weights[0];
                        //gradient_weights[2] implicitly defined
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=PrimaryNorm::IndexToIgnore();
                    } while(false);
                    
                    
                    // estimate from the face [y,z,p]
                    // the quadratic system must be slightly updated to take p into account.
                    
                    if(Dp<NumericTraits<ValueType>::max()){
                        const ValueType N2p=H;
                        ones2+=1./N2p;
                        dist2+=vnl_math_sqr(Dp)/N2p;
                        ones_dist+=Dp/N2p;
                        
                        delta = ones_dist*ones_dist - ones2*(dist2-1);
                        
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const PrimaryVectorType w=N*(DTest*ones-dist);
                        const ValueType w2 = (DTest-Dp)/N2p;
                        
                        if(w[0]<=0 || w[1]<=0 || w2<=0) continue;
                        
                        D=DTest;
                        const ValueType wSum = w[0]+w[1]+w2;
                        gradient_weights[0]= w[0]/wSum;
                        gradient_weights[1] = w[1]/wSum;
                        //gradient_weights[2] implicitly defined
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=IndexP;
                    } // if p
                } // for edges
                
            } else { // y_i>=primarySize
                
                // reference point, of index y_i, is in a different parameter plane.
                // we call it p to be consistent with the first case.
                
                const unsigned int p_i = y_i;
                const ShortOffsetType SOp = stencil[p_i];
                const OffsetType Op = pNorm->ShortOffset_to_Offset(SOp);
                const IndexType p = x+Op;
                const ValueType Dp = d->GetPixel(p);
                const ValueType N2p = H;
                
                //direct test from p
                
                D = Dp+sqrt(N2p);
                gradient_indices[0]=p_i;
                gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                gradient_indices[2]=PrimaryNorm::IndexToIgnore(); //Arbitrary value at position 1 and 2, within GetStencilSize.
                
                gradient_weights[0]=1;
                gradient_weights[1]=0;
                
                // circle around points and segments.
                for(int y_i=0; y_i<primarySize; ++y_i){
                    //hides the previous instance of y_i on purpose (
                    const PrimaryShortOffsetType BSOy = primaryStencil[y_i];
                    const OffsetType Oy = pNorm->PrimaryShortOffset_to_Offset(BSOy);
                    const IndexType y = x+Oy;
                    const PrimaryVectorType Vy(BSOy);
                    
                    if(!b->GetBufferedRegion().IsInside(y)) continue;
                    if(!b->GetPixel(y)) continue;  //Point y must be accepted.
                    const ValueType Dy=d->GetPixel(y);
                    
                    const ValueType N2y = primaryNorm.SquaredNorm(Vy);
                    
                    
                    // Estimate from the segment [y,p] joining y to p
                    do {
                        //solving quadratic system a D^2 - 2 b D + c = 0;
                        const ValueType
                        a = 1./N2y+1./N2p,
                        b = Dy/N2y + Dp/N2p,
                        c = vnl_math_sqr(Dy)/N2y + vnl_math_sqr(Dp)/N2p - 1.;
                        
                        const ValueType delta = vnl_math_sqr(b) - a*c;
                        if(delta < 0) continue;
                        const ValueType DTest = ( b+sqrt(delta) ) / a;
                        if(DTest >= D) continue;
                        
                        PrimaryVectorType w;
                        w[0]=(DTest-Dy)/N2y;
                        w[1]=(DTest-Dp)/N2p;
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        gradient_indices[0]=y_i;
                        gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[2]=p_i;
                        
                        gradient_weights[0]=w[0]/(w[0]+w[1]);
                        gradient_weights[1]=0;
                        // gradient_weights[2] // automatically defined by 1- gradient_weights[0] -gradient_weights[1].
                        
                    } while (false);
                    
                    // Estimate from the face [y,z,p]
                    
                    const int z_i= (y_i+1)%primarySize;
                    const PrimaryShortOffsetType BSOz = primaryStencil[z_i];
                    const OffsetType Oz = pNorm->PrimaryShortOffset_to_Offset(BSOz);
                    const IndexType z = x+Oz;
                    const PrimaryVectorType Vz(BSOz);
                    
                    if(!b->GetBufferedRegion().IsInside(z)) continue;
                    if(!b->GetPixel(z)) continue;  //Point z must be accepted.
                    const ValueType Dz=d->GetPixel(z);
                    
                    const Riemannian2DNorm M(primaryNorm.SquaredNorm(Vy),primaryNorm.ScalarProduct(Vy,Vz),primaryNorm.SquaredNorm(Vz));
                    const Riemannian2DNorm N=M.GetInverse();
                    
                    PrimaryVectorType ones; ones.Fill(1);
                    PrimaryVectorType dist; dist[0]=Dy; dist[1]=Dz;
                    
                    ValueType
                    ones2 = N.SquaredNorm(ones),
                    dist2 = N.SquaredNorm(dist),
                    ones_dist = N.ScalarProduct(ones, dist);
                    
                    // the above coefficients must be modified to take p into account.
                    
                    ones2+=1./N2p;
                    dist2+=vnl_math_sqr(Dp)/N2p;
                    ones_dist+=Dp/N2p;
                    
                    const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    
                    if(delta<0) continue;
                    const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                    
                    if(DTest>=D) continue;
                    
                    const PrimaryVectorType w=N*(DTest*ones-dist);
                    const ValueType w2 = (DTest-Dp)/N2p;
                    
                    if(w[0]<=0 || w[1]<=0 || w2<=0) continue;
                    
                    D=DTest;
                    const ValueType wSum = w[0]+w[1]+w2;
                    gradient_weights[0]= w[0]/wSum;
                    gradient_weights[1] = w[1]/wSum;
                    //gradient_weights[2] implicitly defined
                    
                    gradient_indices[0] = y_i;
                    gradient_indices[1] = z_i;
                    gradient_indices[2] = p_i;
                    
                } // for y_i (renamed..)
                
            } // if y_i
            
            return D;
        }
        
        // ************************************** Riemannian3DNorm *******************************************

        struct Riemannian3DNorm_AccessToStencils : Riemannian3DNorm<> {
            typedef Riemannian3DNorm<> SuperClass;
            static const unsigned int nNeigh(int s, int y_i){return SuperClass::nNeigh[s][y_i];}
            static const unsigned int * Neigh(int s, int y_i){return SuperClass::Neigh[s][y_i];}
        };
        
        template<typename TC, typename TSO, typename AcceptedImageType>
        typename ExtendedNorm<Riemannian3DNorm<TC,TSO> >::ValueType
        Specialized_Hopf_Lax(
                             const ExtendedNorm<Riemannian3DNorm<TC,TSO> > *pNorm,
                             const typename ExtendedNorm<Riemannian3DNorm<TC,TSO> >::DistanceImageType * d,
                             const typename ExtendedNorm<Riemannian3DNorm<TC,TSO> >::IndexType & x,
                             const AcceptedImageType * b,
                             unsigned int y_i,
                             const typename ExtendedNorm<Riemannian3DNorm<TC,TSO> >::Stencil & stencil,
                             typename ExtendedNorm<Riemannian3DNorm<TC,TSO> >::IndicesType & gradient_indices,
                             typename ExtendedNorm<Riemannian3DNorm<TC,TSO> >::WeightsType & gradient_weights
                             ){
            
            typedef Riemannian2DNorm<TC,TSO> Riemannian2DNorm;
            typedef Riemannian3DNorm<TC,TSO> Riemannian3DNorm;
            typedef ExtendedNorm<Riemannian3DNorm> ExtendedNorm;
            
            typedef typename ExtendedNorm::PrimaryVectorType PrimaryVectorType;
            typedef typename ExtendedNorm::ValueType ValueType;
            typedef typename ExtendedNorm::IndexType IndexType;
            typedef typename ExtendedNorm::OffsetType OffsetType;
            typedef typename ExtendedNorm::PrimaryNorm PrimaryNorm;
            typedef typename ExtendedNorm::IndicesType IndicesType;
            typedef typename ExtendedNorm::ShortOffsetType ShortOffsetType;
            
            const unsigned int Dimension = ExtendedNorm::Dimension;
            
            const typename Riemannian3DNorm::Stencil & primaryStencil = stencil.stencil;
            const unsigned int primarySize = primaryStencil.size();
            typedef typename Riemannian3DNorm::ShortOffsetType PrimaryShortOffsetType;
            
            
            const Riemannian3DNorm & primaryNorm = pNorm->GetPrimaryNorm();
            const ValueType H = pNorm->GetScalar();
            
            
            ValueType D; // contains the estimated distance
            
            if(y_i<primarySize){
                
                const PrimaryShortOffsetType BSOy=primaryStencil[y_i];
                const OffsetType Oy = pNorm->PrimaryShortOffset_to_Offset(BSOy);
                const PrimaryVectorType Vy(BSOy);
                
                const IndexType y=x+Oy;
                const ValueType Dy=d->GetPixel(y);
                
                // Estimate from the point y alone
                D = Dy+primaryNorm.Norm(Vy);
                gradient_indices[0]=y_i;
                for(int i=1; i<Dimension; ++i) gradient_indices[i] = PrimaryNorm::IndexToIgnore();
                
                gradient_weights[0]=1;
                for(int i=1; i<Dimension-1; ++i) gradient_weights[i]=0;
                
                
                // Finding out who is the closest point p in parameter space
                // Storing the associated distance and index in Dp and IndexP
                ValueType Dp = NumericTraits<ValueType>::max();
                typename IndicesType::ValueType IndexP;
                
                for(int e=primarySize; e<stencil.size(); ++e){
                    const ShortOffsetType BSOp = stencil[e];
                    OffsetType Op = pNorm->ShortOffset_to_Offset(BSOp);
                    
                    const IndexType p=x+Op;
                    if(!b->GetBufferedRegion().IsInside(p)) continue;
                    if(!b->GetPixel(p)) continue;  //Point p must be accepted.
                    
                    const ValueType DpTest=d->GetPixel(p);
                    if(DpTest<Dp){
                        Dp=DpTest;
                        IndexP=e;
                    }
                } // for e in last two stencil elements
                
                
                
                // Estimate from the segment [y,p] joining y to p
                if(Dp<NumericTraits<ValueType>::max())
                    do {
                        const ValueType
                        N2y = primaryNorm.SquaredNorm(Vy),
                        N2p = H;
                        
                        //solving quadratic system a D^2 - 2 b D + c = 0;
                        const ValueType
                        a = 1./N2y+1./N2p,
                        b = Dy/N2y + Dp/N2p,
                        c = vnl_math_sqr(Dy)/N2y + vnl_math_sqr(Dp)/N2p - 1.;
                        
                        const ValueType delta = vnl_math_sqr(b) - a*c;
                        if(delta < 0) continue;
                        const ValueType DTest = ( b+sqrt(delta) ) / a;
                        if(DTest >= D) continue;
                        
                        PrimaryVectorType w;
                        w[0]=(DTest-Dy)/N2y;
                        w[1]=(DTest-Dp)/N2p;
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[2]=IndexP;
                        
                        gradient_weights[0]=w[0]/(w[0]+w[1]);
                        gradient_weights[1]=0;
                        // gradient_weights[2] // automatically defined by 1- gradient_weights[0] -gradient_weights[0].
                        
                    } while (false);
                
                // Estimates from the segments [y,z] joining y to a physical space neighbor z.
                
                
                
#ifdef SuperBase
                const int s= (primaryStencil.compressedStencilSize==0) ? 0 : 1;
#else
                //what is the type of stencil ? s=0 : classic orthogonal. s=1 negative scalar product. s=2 non-negative.
                const typename Riemannian3DNorm::Stencil::TopologyCase s=primaryStencil.topologyCase;
#endif
                const unsigned int nN=Riemannian3DNorm_AccessToStencils::nNeigh(s,y_i);
                const unsigned int * N = Riemannian3DNorm_AccessToStencils::Neigh(s,y_i);
                
                
                //caching some data
                static const int max_nN = 6;
                bool B[max_nN];
                PrimaryVectorType V[max_nN];
                ValueType d0[max_nN];
                for(int i=0; i<nN; ++i){
                    const unsigned int z_i = N[i];
                    const PrimaryShortOffsetType & BSOz=primaryStencil[z_i];
                    V[i] = PrimaryVectorType(BSOz);
                    const IndexType z=x+pNorm->PrimaryShortOffset_to_Offset(BSOz);
                    if(!b->GetBufferedRegion().IsInside(z)) {B[i]=false; continue;}
                    if(!b->GetPixel(z)) {B[i]=false; continue;}  //Point z must be accepted.
                    B[i]=true;
                    d0[i]=d->GetPixel(z); //Since z is accepted, Dz < infinity
                }
                
                // Minimizations associated to edges [y,z]
                for(int i=0; i<nN; ++i){
                    if(!B[i]) continue;
                    const PrimaryVectorType & Vz = V[i];
                    const ValueType Dz=d0[i];
                    const int z_i = N[i];
                    //quadratic system to solve...
                    
                    Riemannian2DNorm M(primaryNorm.SquaredNorm(Vy),primaryNorm.ScalarProduct(Vy,Vz),primaryNorm.SquaredNorm(Vz));
                    M = M.GetInverse();
                    
                    typedef typename Riemannian2DNorm::VectorType VectorType2;
                    VectorType2 ones; ones.Fill(1);
                    VectorType2 dist; dist[0]=Dy; dist[1]=Dz;
                    
                    ValueType
                    ones2 = M.SquaredNorm(ones),
                    dist2 = M.SquaredNorm(dist),
                    ones_dist = M.ScalarProduct(ones, dist);
                    ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    do {
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const VectorType2 w=M*(DTest*ones-dist);
                        if(w[0]<=0 || w[1]<=0) continue;
                        
                        D=DTest;
                        gradient_weights[0]= w[0]/(w[0]+w[1]);
                        gradient_weights[1]= 1-gradient_weights[0];
                        gradient_weights[2]= 0;
                        //gradient_weights[3] is implicitly defined
                        
                        //gradient_indices[0]=y_i; // already done
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[3]=PrimaryNorm::IndexToIgnore();
                        
                    } while (false);
                    
                    // estimate from the face [y,z,p]
                    // the quadratic system must be slightly updated to take p into account.
                    
                    if(Dp<NumericTraits<ValueType>::max()){
                        const ValueType N2p=H;
                        ones2+=1./N2p;
                        dist2+=vnl_math_sqr(Dp)/N2p;
                        ones_dist+=Dp/N2p;
                        
                        delta = ones_dist*ones_dist - ones2*(dist2-1);
                        
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const VectorType2 w=M*(DTest*ones-dist);
                        const ValueType w2 = (DTest-Dp)/N2p;
                        
                        if(w[0]<=0 || w[1]<=0 || w2<=0) continue;
                        
                        D=DTest;
                        const ValueType wSum = w[0]+w[1]+w2;
                        gradient_weights[0]= w[0]/wSum;
                        gradient_weights[1]= w[1]/wSum;
                        gradient_weights[2]= 0;
                        //gradient_weights[3] implicitly defined
                        
                        //gradient_indiced[0]=y_i; //already set
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=PrimaryNorm::IndexToIgnore();
                        gradient_indices[3]=IndexP;
                    } // if p
                    
                    
                } // for edges
                
                // Minimizations associated to faces [y,z,t]
                for(int i=0; i<nN; ++i){
                    const unsigned int j=(i+1)%nN;
                    if(!B[i] || !B[j]) continue;
                    
                    const int z_i = N[i];
                    const ValueType Dz = d0[i];
                    const PrimaryVectorType & Vz = V[i];
                    
                    const int t_i = N[j];
                    const ValueType Dt = d0[j];
                    const PrimaryVectorType & Vt = V[j];
                    
                    Riemannian3DNorm M;
                    M(0,0) = primaryNorm.SquaredNorm(Vy);
                    M(1,0) = primaryNorm.ScalarProduct(Vy,Vz);   M(1,1) = primaryNorm.SquaredNorm(Vz);
                    M(2,0) = primaryNorm.ScalarProduct(Vy,Vt);   M(2,1) = primaryNorm.ScalarProduct(Vz,Vt);   M(2,2)=primaryNorm.SquaredNorm(Vt);
                    
                    M=M.GetInverse();
                    PrimaryVectorType ones; ones.Fill(1);
                    PrimaryVectorType dist; dist[0]=Dy; dist[1]=Dz; dist[2]=Dt;
                    
                    //same equation
                    ValueType
                    ones2 = M.SquaredNorm(ones),
                    dist2 = M.SquaredNorm(dist),
                    ones_dist = M.ScalarProduct(ones, dist);
                    
                    ValueType
                    delta = ones_dist*ones_dist - ones2*(dist2-1);
                    do {
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const PrimaryVectorType w=M*(DTest*ones-dist);
                        if(w[0]<=0 || w[1]<=0 || w[2]<=0) continue;
                        
                        D=DTest;
                        const ValueType wSum = w[0]+w[1]+w[2];
                        gradient_weights[0]= w[0]/wSum;
                        gradient_weights[1]= w[1]/wSum;
                        gradient_weights[2]= w[2]/wSum;
                        // implicitly defined gradient_weights[3] equals 0.
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]= z_i;
                        gradient_indices[2]= t_i;
                        gradient_indices[3]=PrimaryNorm::IndexToIgnore();
                    } while(false);
                    
                    // Estimate from the 3D facet [y,z,t,p]
                    
                    if(Dp<NumericTraits<ValueType>::max()){
                        const ValueType N2p=H;
                        ones2+=1./N2p;
                        dist2+=vnl_math_sqr(Dp)/N2p;
                        ones_dist+=Dp/N2p;
                        
                        delta = ones_dist*ones_dist - ones2*(dist2-1);
                        
                        if(delta<0) continue;
                        const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                        
                        if(DTest>=D) continue;
                        
                        const PrimaryVectorType w=M*(DTest*ones-dist);
                        const ValueType w3 = (DTest-Dp)/N2p;
                        
                        if(w[0]<=0 || w[1]<=0 || w[2]<= 0 || w3<=0) continue;
                        
                        D=DTest;
                        const ValueType wSum = w[0]+w[1]+w[2]+w3;
                        gradient_weights[0]= w[0]/wSum;
                        gradient_weights[1]= w[1]/wSum;
                        gradient_weights[2]= w[2]/wSum;
                        //gradient_weights[3] implicitly defined
                        
                        //gradient_indices[0]=y_i; //already set
                        gradient_indices[1]=z_i;
                        gradient_indices[2]=t_i;
                        gradient_indices[3]=IndexP;
                    } // if p
                }
                
            } else { //y_i >= primarySize
                // reference point, of index y_i, is in a different parameter plane.
                // we call it p to be consistent with the first case.
                
                const unsigned int p_i = y_i;
                const ShortOffsetType SOp = stencil[p_i];
                const OffsetType Op = pNorm->ShortOffset_to_Offset(SOp);
                const IndexType p = x+Op;
                const ValueType Dp = d->GetPixel(p);
                const ValueType N2p = H;
                
                //direct test from p
                
                D = Dp+sqrt(N2p);
                
                for(int i=0; i<Dimension-1; ++i) gradient_indices[i]=PrimaryNorm::IndexToIgnore();
                gradient_indices[Dimension-1]=p_i;
                
                for(int i=0; i<Dimension-1; ++i) gradient_weights[i]=0;
                //gradient weights[Dimension-1] implicitly defined to 1
                
                //Cache values
                static const int max_N = 14; // max stencil size
                bool B[max_N];
                PrimaryVectorType V[max_N];
                ValueType d0[max_N];
                
                for(int i=0; i<primarySize; ++i){
                    const PrimaryShortOffsetType & BSOz=primaryStencil[i];
                    V[i] = PrimaryVectorType(BSOz);
                    const IndexType z=x+pNorm->PrimaryShortOffset_to_Offset(BSOz);
                    if(!b->GetBufferedRegion().IsInside(z)) {B[i]=false; continue;}
                    if(!b->GetPixel(z)) {B[i]=false; continue;}  //Point z must be accepted.
                    B[i]=true;
                    d0[i]=d->GetPixel(z); //Since z is accepted, Dz < infinity
                }
                
                
                // enumerate all stencil points in physical space.
                for(int y_i=0; y_i<primarySize; ++y_i){
                    //hides the previous instance of y_i on purpose
                    
                    if(!B[y_i]) continue;
                    const ValueType Dy = d0[y_i];
                    const PrimaryVectorType & Vy = V[y_i];
                    const ValueType N2y = primaryNorm.SquaredNorm(Vy);
                    
                    // Estimate from the segment [y,p] joining y to p
                    
                    //solving quadratic system a D^2 - 2 b D + c = 0;
                    const ValueType
                    a = 1./N2y+1./N2p,
                    b = Dy/N2y + Dp/N2p,
                    c = vnl_math_sqr(Dy)/N2y + vnl_math_sqr(Dp)/N2p - 1.;
                    
                    const ValueType delta = vnl_math_sqr(b) - a*c;
                    if(delta < 0) continue;
                    const ValueType DTest = ( b+sqrt(delta) ) / a;
                    if(DTest >= D) continue;
                    
                    PrimaryVectorType w;
                    w[0]=(DTest-Dy)/N2y;
                    w[1]=(DTest-Dp)/N2p;
                    if(w[0]<=0 || w[1]<=0) continue;
                    
                    D=DTest;
                    gradient_indices[0]=y_i;
                    gradient_indices[1]=PrimaryNorm::IndexToIgnore();
                    gradient_indices[2]=PrimaryNorm::IndexToIgnore();
                    // gradient_indices[3]=p_i; // already done
                    
                    gradient_weights[0]=w[0]/(w[0]+w[1]);
                    gradient_weights[1]=0.;
                    gradient_weights[2]=0.;
                    // gradient_weights[3] is implicitly defined
                    
                } // for y_i (renamed..)
                
                
                // number of segments in physical space, by topology type
#ifdef SuperBase
                static const int NumberOfSegments_[2] = {12,36};
                const int s= (primaryStencil.compressedStencilSize==0) ? 0 : 1;
#else
                static const int NumberOfSegments_[3] = {12,36,36};
                const typename Riemannian3DNorm::Stencil::TopologyCase s=primaryStencil.topologyCase;
#endif
                //                        const unsigned int nN=Riemannian3DNorm::nNeigh[s][y_i];
                //                        const unsigned int * N = Riemannian3DNorm::Neigh[s][y_i];
                
                const int NumberOfSegments = NumberOfSegments_[s];
                
                // All stencil segments in physical space
                for(int seg=0; seg<NumberOfSegments; ++seg){
                    
                    // list of all stencil segments in physical space
#ifdef SuperBase
                    static const int ListOfSegments[2][36][2] = {
                        {
                            {0, 2}, {0, 3}, {0, 4}, {0, 5},
                            {1, 2}, {1, 3}, {1, 4}, {1, 5},
                            {2, 4}, {2, 5}, {3, 4}, {3, 5}
                        }, {
                            { 8, 6},{ 9, 5},{ 8, 7},{11, 5},
                            { 9, 7},{11, 6},{10, 4},{12, 4},
                            {10, 7},{12, 6},{13, 4},{13, 5},
                            { 1,10},{ 2,10},{ 1,12},{ 3,12},
                            { 2,13},{ 3,13},{ 0, 9},{ 2, 9},
                            { 0,11},{ 3,11},{ 0, 8},{ 1, 8},
                            { 0, 6},{ 1, 6},{ 0, 7},{ 1, 7},
                            { 0, 5},{ 2, 5},{ 2, 7},{ 1, 4},
                            { 2, 4},{ 3, 5},{ 3, 6},{ 3, 4}
                        }
                    };
#else
                    static const int ListOfSegments[3][36][2] = {
                        {
                            {0, 2}, {0, 3}, {0, 4}, {0, 5},
                            {1, 2}, {1, 3}, {1, 4}, {1, 5},
                            {2, 4}, {2, 5}, {3, 4}, {3, 5}
                        }, {
                            {0, 2}, {0, 4}, {0, 7}, {0, 9},
                            {0,10}, {0,13}, {1, 3}, {1, 5},
                            {1, 6}, {1, 8}, {1,11}, {1,12},
                            {2, 5}, {2, 6}, {2, 9}, {2,10},
                            {2,12}, {3, 4}, {3, 7}, {3, 8},
                            {3,11}, {3,13}, {4, 7}, {4, 8},
                            {4,10}, {4,12}, {5, 6}, {5, 9},
                            {5,11}, {5,13}, {6,12}, {7,13},
                            {8,12}, {9,13},{10,12},{11,13}
                        }, {
                            {0, 2}, {0, 4}, {0, 7}, {0, 9},
                            {0,10}, {0,13}, {1, 3}, {1, 5},
                            {1, 6}, {1, 8}, {1,11}, {1,12},
                            {2, 4}, {2, 6}, {2, 8}, {2,13},
                            {3, 5}, {3, 7}, {3, 9}, {3,12},
                            {4, 8}, {4,10}, {4,12}, {5, 9},
                            {5,11}, {5,13}, {6, 8}, {6,11},
                            {6,13}, {7, 9}, {7,10}, {7,12},
                            {8,12}, {9,13},{10,12},{11,13}}
                    };
#endif
                    
                    const int y_i = ListOfSegments[s][seg][0];
                    const int z_i = ListOfSegments[s][seg][1];
                    
                    if(!B[y_i] || !B[z_i]) continue;
                    
                    const ValueType Dy = d0[y_i];
                    const ValueType Dz = d0[z_i];
                    
                    const PrimaryVectorType & Vy = V[y_i];
                    const PrimaryVectorType & Vz = V[z_i];
                    
                    // Estimate from the face [y,z,p]
                    
                    const Riemannian2DNorm M(primaryNorm.SquaredNorm(Vy),primaryNorm.ScalarProduct(Vy,Vz),primaryNorm.SquaredNorm(Vz));
                    const Riemannian2DNorm N=M.GetInverse();
                    
                    typedef typename Riemannian2DNorm::VectorType VectorType2;
                    
                    VectorType2 ones; ones.Fill(1);
                    VectorType2 dist; dist[0]=Dy; dist[1]=Dz;
                    
                    ValueType
                    ones2 = N.SquaredNorm(ones),
                    dist2 = N.SquaredNorm(dist),
                    ones_dist = N.ScalarProduct(ones, dist);
                    
                    // the above coefficients must be modified to take p into account.
                    
                    ones2+=1./N2p;
                    dist2+=vnl_math_sqr(Dp)/N2p;
                    ones_dist+=Dp/N2p;
                    
                    const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    
                    if(delta<0) continue;
                    const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                    
                    if(DTest>=D) continue;
                    
                    typedef typename Riemannian2DNorm::VectorType VectorType2;
                    const VectorType2 w=N*(DTest*ones-dist);
                    const ValueType w2 = (DTest-Dp)/N2p;
                    
                    if(w[0]<=0 || w[1]<=0 || w2<=0) continue;
                    
                    D=DTest;
                    const ValueType wSum = w[0]+w[1]+w2;
                    gradient_weights[0] = w[0]/wSum;
                    gradient_weights[1] = w[1]/wSum;
                    gradient_weights[2] = 0;
                    //                            gradient_weights[3] implicitly defined
                    
                    gradient_indices[0] = y_i;
                    gradient_indices[1] = z_i;
                    gradient_indices[2] = PrimaryNorm::IndexToIgnore();
                    //                            gradient_indices[3] = p_i; //already done
                }
                
                // All stencil faces in physical space
#ifdef SuperBase
                static const int NumberOfFaces_[2] = {8,24};
#else
                static const int NumberOfFaces_[3] = {8,24,24};
#endif
                
                const int NumberOfFaces = NumberOfFaces_[s];
                
                // All stencil faces in physical space
                for(int face=0; face<NumberOfFaces; ++face){
#ifdef SuperBase
                    static const int ListOfFaces[3][24][3]={
                        {
                            {0, 2, 4}, {0, 2, 5}, {0, 3, 4}, {0, 3, 5},
                            {1, 2, 4}, {1, 2, 5}, {1, 3, 4}, {1, 3, 5}
                        }, {
                            { 0, 8, 6},{ 0, 9, 5},{ 1, 8, 6},{ 1,10, 4},
                            { 2, 9, 5},{ 2,10, 4},{ 0, 8, 7},{ 0,11, 5},
                            { 1, 8, 7},{ 1,12, 4},{ 3,11, 5},{ 3,12, 4},
                            { 0, 9, 7},{ 0,11, 6},{ 2, 9, 7},{ 2,13, 4},
                            { 3,11, 6},{ 3,13, 4},{ 1,10, 7},{ 1,12, 6},
                            { 2,10, 7},{ 2,13, 5},{ 3,12, 6},{ 3,13, 5}
                        }
                    };
#else
                    static const int ListOfFaces[3][24][3]={
                        {
                            {0, 2, 4}, {0, 2, 5}, {0, 3, 4}, {0, 3, 5},
                            {1, 2, 4}, {1, 2, 5}, {1, 3, 4}, {1, 3, 5}
                        }, {
                            {0, 2, 9}, {0, 2,10}, {0, 4, 7}, {0, 4,10},
                            {0, 7,13}, {0, 9,13}, {1, 3, 8}, {1, 3,11},
                            {1, 5, 6}, {1, 5,11}, {1, 6,12}, {1, 8,12},
                            {2, 5, 6}, {2, 5, 9}, {2, 6,12}, {2,10,12},
                            {3, 4, 7}, {3, 4, 8}, {3, 7,13}, {3,11,13},
                            {4, 8,12}, {4,10,12}, {5, 9,13}, {5,11,13}
                        }, {
                            {0, 2, 4}, {0, 2,13}, {0, 4,10}, {0, 7, 9},
                            {0, 7,10}, {0, 9,13}, {1, 3, 5}, {1, 3,12},
                            {1, 5,11}, {1, 6, 8}, {1, 6,11}, {1, 8,12},
                            {2, 4, 8}, {2, 6, 8}, {2, 6,13}, {3, 5, 9},
                            {3, 7, 9}, {3, 7,12}, {4, 8,12}, {4,10,12},
                            {5, 9,13}, {5,11,13}, {6,11,13}, {7,10,12}
                        }
                    };
#endif
                    
                    const int y_i = ListOfFaces[s][face][0];
                    const int z_i = ListOfFaces[s][face][1];
                    const int t_i = ListOfFaces[s][face][2];
                    
                    if(!B[y_i] || !B[z_i] || !B[t_i]) continue;
                    
                    const ValueType Dy = d0[y_i];
                    const ValueType Dz = d0[z_i];
                    const ValueType Dt = d0[t_i];
                    
                    
                    const PrimaryVectorType & Vy = V[y_i];
                    const PrimaryVectorType & Vz = V[z_i];
                    const PrimaryVectorType & Vt = V[t_i];
                    
                    // Estimate from the face [y,z,t,p]
                    
                    Riemannian3DNorm M;
                    M(0,0) = primaryNorm.SquaredNorm(Vy);
                    M(1,0) = primaryNorm.ScalarProduct(Vy,Vz);   M(1,1) = primaryNorm.SquaredNorm(Vz);
                    M(2,0) = primaryNorm.ScalarProduct(Vy,Vt);   M(2,1) = primaryNorm.ScalarProduct(Vz,Vt);   M(2,2)=primaryNorm.SquaredNorm(Vt);
                    
                    const Riemannian3DNorm N=M.GetInverse();
                    
                    PrimaryVectorType ones; ones.Fill(1);
                    PrimaryVectorType dist; dist[0]=Dy; dist[1]=Dz; dist[2]=Dt;
                    
                    ValueType
                    ones2 = N.SquaredNorm(ones),
                    dist2 = N.SquaredNorm(dist),
                    ones_dist = N.ScalarProduct(ones, dist);
                    
                    // the above coefficients must be modified to take p into account.
                    
                    ones2+=1./N2p;
                    dist2+=vnl_math_sqr(Dp)/N2p;
                    ones_dist+=Dp/N2p;
                    
                    const ValueType delta = ones_dist*ones_dist - ones2*(dist2-1);
                    
                    if(delta<0) continue;
                    const ValueType DTest = (ones_dist+sqrt(delta))/ones2; //estimated distance
                    
                    if(DTest>=D) continue;
                    
                    const PrimaryVectorType w=N*(DTest*ones-dist);
                    const ValueType w3 = (DTest-Dp)/N2p;
                    
                    if(w[0]<=0 || w[1]<=0 || w[2]<=0 || w3<=0) continue;
                    
                    D=DTest;
                    const ValueType wSum = w[0]+w[1]+w[2]+w3;
                    gradient_weights[0] = w[0]/wSum;
                    gradient_weights[1] = w[1]/wSum;
                    gradient_weights[2] = w[2]/wSum;
                    //   gradient_weights[3] implicitly defined
                    
                    gradient_indices[0] = y_i;
                    gradient_indices[1] = z_i;
                    gradient_indices[2] = t_i;
                    //  gradient_indices[3] = p_i; // already done.
                }
                
                
            } // if y_i
            return D;
        } // end of Riemannian3DNorm case
        
    } // end of anonymous namespace
    
    
    template<typename TNorm>
    template<typename AcceptedImageType>
    typename ExtendedNorm< TNorm >::ValueType
    ExtendedNorm< TNorm >::Hopf_Lax(const DistanceImageType * d, const IndexType &x,
                                             const AcceptedImageType * b, unsigned int y_i,
                                             const Stencil & stencil,
                                             IndicesType &gradient_indices, WeightsType &gradient_weights) const
    {
        return Specialized_Hopf_Lax<ValueType, typename TNorm::ShortOffsetValueType, AcceptedImageType>
        (this,d,x,b,y_i,stencil,gradient_indices,gradient_weights);
    }
    
} // end of namespace ITK
#endif
