//
//  itkAnisotropicFastMarchingImageFilter_Geodesic.hxx
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 10/06/13.
//
//

#ifndef ITKFM_itkAnisotropicFastMarchingImageFilter_Geodesic_hxx
#define ITKFM_itkAnisotropicFastMarchingImageFilter_Geodesic_hxx

namespace itk {

    template<typename TN,typename TSS>
    typename AnisotropicFastMarchingImageFilter<TN,TSS>::ContinuousIndexType
    AnisotropicFastMarchingImageFilter<TN,TSS>::GeodesicContinuousIndex::GetContinuousIndex() const
    {
        ContinuousIndexType I;
        for(int k=0; k<Dimension; ++k) I[k]=this->first[k]+this->second[k];
        return I;
    }
    
    template<typename TN,typename TSS> template<typename TComponent>
    AnisotropicFastMarchingImageFilter<TN,TSS>::GeodesicContinuousIndex::
    GeodesicContinuousIndex
    (const ContinuousIndex<TComponent,Dimension> & I)
    {
        this->first.CopyWithRound(I);
        for(int k=0; k<Dimension; ++k) this->second[k] = I[k]-this->first[k];        
    }
    
    
    // ******************* Geodesic ********************
    
    //P : grid point. V : correction.
    template <typename TN, typename TSS>
    void AnisotropicFastMarchingImageFilter<TN,TSS>::Geodesic
    (std::vector<GeodesicContinuousIndex> & PV) const
    {
        if(!GetGenerateUpwindGradient())
            itkExceptionMacro("AnisotropicFastMarchingImageFilter::Geodesic error : GenerateUpwindGradient flag was not set !");
        
        if(PV.size()!=1)
            itkExceptionMacro("AnisotropicFastMarchingImageFilter::Geodesic error : exactly one starting point is required, "
                              << PV.size() << " were given.");
        
        
        IndexType p = PV[0].first;
        VectorType v = PV[0].second;
        
        if(!Stencils.BufferedRegion.IsInside(p))
            itkExceptionMacro("AnisotropicFastMarchingImageFilter::Geodesic error : the starting point " << p
                              << "is not inside the buffered region " << Stencils.BufferedRegion)
            
            while(true){
                
                const typename GradientRawData::RawType & G = Gradient.Raw->GetPixel(p);
                const typename GradientRawData::IndicesType & indices = G.first;
                const typename GradientRawData::WeightsType & weights = G.second;
                
                if(weights.IsNullWeights()) break;
                
                DirectStencilType stencil = Stencils.Direct(p);
                
                VectorType grad;
                grad.Fill(0);
                
                for(int i=0; i<Dimension; ++i){
                    if(indices[i]==NormType::IndexToIgnore()) continue;
                    const ShortOffsetType u = stencil[indices[i]];
                    grad+=weights(i)*VectorType(u);
                }
                
                //                cout << endl << "p : " << p << "; grad " << grad << "; v " << v << "; weights" << weights << endl;
                //            cout << "Distance value at p : " << this->GetOutput()->GetPixel(p) << endl;
                
                //            for(int i=0; i<Dimension; ++i){
                //                if(indices[i]==NormType::IndexToIgnore()) cout << "ToIgnore" << endl;
                //                else cout << VectorType(stencil[indices[i]]) << endl;
                //            }
                
                ShortOffsetType best_u;
                ValueType best_error=NumericTraits<ValueType>::max();
                VectorType best_correction;
                
                for(int i=0; i<Dimension; ++i){
                    if(indices[i]==NormType::IndexToIgnore()) continue;
                    const ShortOffsetType u = stencil[indices[i]];
                    const VectorType q = v - VectorType(u);
                    
                    const ValueType scalar = std::max(0., -(q*grad)/grad.GetSquaredNorm() );
                    const VectorType tested_correction = q+grad*scalar;
                    const ValueType tested_error = tested_correction.GetSquaredNorm();
                    
                    //                            cout << "v : " << v << "; q" << q << "; scalar " << scalar << "; tested_correction " << tested_correction << endl;
                    
                    if(best_error <= tested_error) continue;
                    best_u=u;
                    best_error=tested_error;
                    best_correction=tested_correction;
                }
                
                for(int i=0; i<Dimension; ++i) p[i]+=best_u[i];
                v=best_correction;
                PV.push_back(GeodesicContinuousIndex(p,v));
            }
    } //Geodesic
    
    /*
    template<typename TN, typename TSS>
    typename itk::AnisotropicFastMarchingImageFilter<TN,TSS>::PathType::Pointer
    itk::AnisotropicFastMarchingImageFilter<TN,TSS>::Geodesic
    (const PathContinuousIndexType & start) const
    {
        PathRawType pathRaw;
        {
            ContinuousIndexType p;
            p.CastFrom(start);
            pathRaw.push_back(GeodesicContinuousIndex(p));
        }
        Geodesic(pathRaw);
        typename PathType::Pointer path = PathType::New();
        
        for(typename PathRawType::const_iterator it=pathRaw.begin(); it!=pathRaw.end(); ++it)
        {
            PathContinuousIndexType p;
            p.CastFrom(it->GetContinuousIndex());
            path->AddVertex(p);
        }
        return path;
    }
    */
    
    template <typename TN, typename TSS>
    bool AnisotropicFastMarchingImageFilter<TN,TSS>::Geodesic
    (std::vector<ContinuousIndexType> & indices) const
    {
        if(indices.size()!=1)
            itkExceptionMacro("AnisotropicFastMarchingImageFilter::Geodesic error : exactly one starting point is required, "
                              << indices.size() << " were given.");
        
        const ContinuousIndexType & c = indices.back();
        if( !Stencils.BufferedRegion.IsInside(c)) return false;

        std::vector<GeodesicContinuousIndex> PV;
        PV.push_back(GeodesicContinuousIndex(c));
        Geodesic(PV);
        for(size_t i=1; i<PV.size(); ++i)
            indices.push_back(PV[i].GetContinuousIndex());

        return true;
    }
    
    template <typename TN, typename TSS>
    bool AnisotropicFastMarchingImageFilter<TN,TSS>::Geodesic
    (std::vector<PointType> & points) const
    {
        if(points.size()!=1)
            itkExceptionMacro("AnisotropicFastMarchingImageFilter::Geodesic error : exactly one starting point is required, "
                              << points.size() << " were given.");
        
        itk::ContinuousIndex<ValueType,Dimension> c;
        const NormImageType * Norms = this->GetInput();
        if( ! Norms->TransformPhysicalPointToContinuousIndex(points.back(),c) )
            return false;
        
        std::vector<ContinuousIndexType> indices;
        indices.push_back(c);
        Geodesic(indices);
        
        for(size_t i=1; i<indices.size(); ++i){
            PointType q;
            Norms->TransformContinuousIndexToPhysicalPoint(indices[i],q);
            points.push_back(q);
        }
        return true;
    }

}

#endif
