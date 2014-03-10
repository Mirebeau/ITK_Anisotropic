//
//  IsotropicGeodesic.h
//  MatlabITKFastMarching
//
//  Created by Jean-Marie Mirebeau on 05/04/13.
//
//

#ifndef MatlabITKFastMarching_IsotropicGeodesic_h
#define MatlabITKFastMarching_IsotropicGeodesic_h


// Vector stands here for (unavailable) continuous offset type.
template<typename ValueType, int Dimension>
int  IsoFM_Geodesic(std::vector<std::pair<itk::Index<Dimension>, itk::Vector<ValueType,Dimension> > > & Geodesic,
const typename itk::FastMarchingUpwindGradientImageFilter<itk::Image<ValueType,Dimension> >::GradientImageType * Gradient)
{
    if(Geodesic.size()!=1) {
        itkGenericExceptionMacro("Iso_FM_Geodesic error : exactly one starting point is required, "
                                 << Geodesic.size() << " were given.");
        return EXIT_FAILURE;
    }
    
    typedef itk::Index<Dimension> IndexType;
    typedef itk::Vector<ValueType,Dimension> VectorType;
    typedef itk::CovariantVector<ValueType,Dimension> CovariantVectorType;
    typedef itk::Point<ValueType,Dimension> PointType;
    typedef itk::ContinuousIndex<ValueType,Dimension> ContinuousIndexType;
    
    IndexType p = Geodesic[0].first;
    VectorType v = Geodesic[0].second;
    
    if( ! Gradient->GetBufferedRegion().IsInside(p) )
        itkGenericExceptionMacro("IsoFM_Geodesic error : starting point is out of domain bounds.");
    
    
    itk::FixedArray<VectorType, 2*Dimension> Stencil;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<2; ++j){
            Stencil[2*i+j].Fill(0);
            Stencil[2*i+j][i] = 1-2*j;
        }
    
    const itk::IndexValueType IndexToIgnore = itk::NumericTraits<itk::IndexValueType>::max();

#pragma message("To do: convert components to ValueType")
    typedef itk::Matrix<double,Dimension,Dimension> MatrixType;
    const MatrixType Direction = Gradient->GetDirection();
    const MatrixType DualDirection = MatrixType(MatrixType(Direction.GetTranspose()).GetInverse());
    
    int counter=0;
    while(true){
        if(counter++==10000) itkGenericExceptionMacro("IsoFM_Geodesic early abort : Stopping Geodesic computation after 10000 steps");
        
        IndexType indices;
        VectorType grad;
        {
            if( ! Gradient->GetBufferedRegion().IsInside(p) )
                itkGenericExceptionMacro("IsoFM_Geodesic error : geodesic led out of domain bounds.");
            
            const CovariantVectorType PhysicalGradient = Gradient->GetPixel(p);
            // Note the minus sign
            const CovariantVectorType LocalOppositeGradient = - (DualDirection*PhysicalGradient);
            
            for(int i=0; i<Dimension; ++i){
                ValueType &g=grad[i];
                g=LocalOppositeGradient[i];
                indices[i] = (g == 0) ? IndexToIgnore : (g>0 ? 2*i : 2*i+1);
            }
            
            if(grad.GetSquaredNorm() == 0) {
                Geodesic.push_back(std::pair<IndexType, VectorType>(p,VectorType(0.)));
                break;
            }
        }
                
        VectorType best_u;
        ValueType best_error=itk::NumericTraits<ValueType>::max();
        VectorType best_correction;
        
        for(int i=0; i<Dimension; ++i){
            if(indices[i]==IndexToIgnore) continue;
            const VectorType u = Stencil[indices[i]];
            const VectorType q = v - u;
            
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
        Geodesic.push_back(std::pair<IndexType, VectorType>(p,v));
    }
    
    return EXIT_SUCCESS;
}



template<typename ValueType, int Dimension>
int IsoFM_Geodesic(std::vector<itk::ContinuousIndex<ValueType, Dimension>  > & Geodesic,
                   const typename itk::FastMarchingUpwindGradientImageFilter<itk::Image<ValueType,Dimension> >::GradientImageType * Gradient)
{
    if(Geodesic.size()!=1) {
        itkGenericExceptionMacro("Iso_FM_Geodesic error : exactly one starting point is required, "
                                 << Geodesic.size() << " were given.");
        return EXIT_FAILURE;
    }
    
    typedef itk::Index<Dimension> IndexType;
    typedef itk::Vector<ValueType,Dimension> VectorType;
    typedef itk::ContinuousIndex<ValueType,Dimension> ContinuousIndexType;
    typedef std::pair<IndexType,VectorType> IndexVectorPairType;
    
    std::vector<IndexVectorPairType> rawGeodesic;
    rawGeodesic.resize(1);
    IndexType & p = rawGeodesic.back().first;
    VectorType & v = rawGeodesic.back().second;
    const ContinuousIndexType & InitialIndex = Geodesic.back();
    
    for(int i=0; i<Dimension; ++i){
        p[i] = itk::IndexValueType( InitialIndex[i] );
        v[i] = InitialIndex[i]-p[i];
    }

    IsoFM_Geodesic<ValueType,Dimension>(rawGeodesic,Gradient);
    
    Geodesic.resize(rawGeodesic.size());
    for(int i=1; i<Geodesic.size(); ++i){
        for(int j=0; j<Dimension; ++j){
            Geodesic[i][j] = rawGeodesic[i].first[j]+rawGeodesic[i].second[j];
        }
    }
    
    return EXIT_SUCCESS;
}

template<typename ValueType, int Dimension>
int IsoFM_Geodesic(std::vector<itk::Point<ValueType,Dimension> > & Geodesic,
                   const typename itk::FastMarchingUpwindGradientImageFilter<itk::Image<ValueType,Dimension> >::GradientImageType * Gradient)
{
    if(Geodesic.size()!=1) {
        itkGenericExceptionMacro("Iso_FM_Geodesic error : exactly one starting point is required, "
                                 << Geodesic.size() << " were given.");
        return EXIT_FAILURE;
    }
    
    typedef itk::ContinuousIndex<ValueType,Dimension> ContinuousIndexType;
    std::vector<ContinuousIndexType> rawGeodesic;
    
    rawGeodesic.resize(1);
    ContinuousIndexType & InitialIndex = rawGeodesic.back();
    if( !Gradient->TransformPhysicalPointToContinuousIndex(Geodesic[0],InitialIndex) ){
        itkGenericExceptionMacro("Iso_FM_Geodesic error : Initial point " << Geodesic[0]
                                 << " is out of image region.\n" << Gradient->GetBufferedRegion() );
    }
    
    IsoFM_Geodesic<ValueType,Dimension>(rawGeodesic, Gradient);
    
    Geodesic.resize(rawGeodesic.size());
    for(int i=0; i<Geodesic.size(); ++i)
        Gradient->TransformContinuousIndexToPhysicalPoint(rawGeodesic[i],Geodesic[i]);
    
    return EXIT_SUCCESS;
}

/*
template<typename ValueType, int Dimension>
int IsoFM_Geodesic(std::vector<itk::Point<ValueType,Dimension> > & Geodesic,
                   typename itk::FastMarchingUpwindGradientImageFilter<itk::Image<ValueType,Dimension> >::GradientImagePointer Gradient)
{    
    if(Geodesic.size()!=1) {
        itkGenericExceptionMacro("Iso_FM_Geodesic error : exactly one starting point is required, "
                                 << Geodesic.size() << " were given.");
        return EXIT_FAILURE;
    }
    
    typedef itk::Index<Dimension> IndexType;
    typedef itk::Vector<ValueType,Dimension> VectorType;
    typedef itk::CovariantVector<ValueType,Dimension> CovariantVectorType;
    typedef itk::Point<ValueType,Dimension> PointType;
    typedef itk::ContinuousIndex<ValueType,Dimension> ContinuousIndexType;
    
    IndexType p;
    VectorType v;
    {
        ContinuousIndexType InitialIndex;
        if( !Gradient->TransformPhysicalPointToContinuousIndex(Geodesic[0],InitialIndex) ){
            itkGenericExceptionMacro("Iso_FM_Geodesic error : Initial point " << Geodesic[0]
                                     << " is out of image region.\n" << Gradient->GetBufferedRegion() );
        }
        for(int i=0; i<Dimension; ++i){
            p[i] = itk::IndexValueType( InitialIndex[i] );
            v[i] = InitialIndex[i]-p[i];
        }
    }
    
    itk::FixedArray<VectorType, 2*Dimension> Stencil;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<2; ++j){
            Stencil[2*i+j].Fill(0);
            Stencil[2*i+j][i] = 1-2*j;
        }
    
    //MexMsg() << Stencil << "\n";
    
    const itk::IndexValueType IndexToIgnore = itk::NumericTraits<itk::IndexValueType>::max();
    
    //typedef typename itk::Image<itk::CovariantVector<ValueType,Dimension> >::DirectionType
    typedef itk::Matrix<ValueType,Dimension,Dimension> MatrixType;
    const MatrixType Direction = Gradient->GetDirection();
    const MatrixType DualDirection = MatrixType(MatrixType(Direction.GetTranspose()).GetInverse());
    
    
    int counter=0;
    while(true){
        if(counter++==10000) itkGenericExceptionMacro("IsoFM_Geodesic early abort : Stopping Geodesic computation after 10000 steps");
        
        IndexType indices;
        VectorType grad;
        {
            if( ! Gradient->GetBufferedRegion().IsInside(p) )
                itkGenericExceptionMacro("IsoFM_Geodesic error : geodesic led out of domain bounds");
            
            const CovariantVectorType PhysicalGradient = Gradient->GetPixel(p);
            // Note the minus sign
            const CovariantVectorType LocalOppositeGradient = - (DualDirection*PhysicalGradient);
            
            for(int i=0; i<Dimension; ++i){
                ValueType &g=grad[i];
                g=LocalOppositeGradient[i];
                indices[i] = (g == 0) ? IndexToIgnore : (g>0 ? 2*i : 2*i+1);
            }
            
            if(grad.GetSquaredNorm() == 0) {
                PointType FinalPoint;
                Gradient->TransformIndexToPhysicalPoint(p,FinalPoint);
                Geodesic.push_back(FinalPoint);
                break;
            }
        }
        
        //        MexMsg() << "indices : " << indices << ", grad : " << grad << "\n";
        //                cout << endl << "p : " << p << "; grad " << grad << "; v " << v << "; weights" << weights << endl;
        //            cout << "Distance value at p : " << this->GetOutput()->GetPixel(p) << endl;
        
        //            for(int i=0; i<Dimension; ++i){
        //                if(indices[i]==NormType::IndexToIgnore()) cout << "ToIgnore" << endl;
        //                else cout << VectorType(stencil[indices[i]]) << endl;
        //            }
        
        VectorType best_u;
        ValueType best_error=itk::NumericTraits<ValueType>::max();
        VectorType best_correction;
        
        for(int i=0; i<Dimension; ++i){
            if(indices[i]==IndexToIgnore) continue;
            const VectorType u = Stencil[indices[i]];
            const VectorType q = v - u;
            
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
        ContinuousIndexType IntermediateIndex;
        for(int i=0; i<Dimension; ++i) IntermediateIndex[i] = p[i]+v[i];
        PointType IntermediatePoint;
        Gradient->TransformContinuousIndexToPhysicalPoint(IntermediateIndex,IntermediatePoint);
        Geodesic.push_back(IntermediatePoint);
    }
    
    return EXIT_SUCCESS;
} //Geodesic
*/
#endif
