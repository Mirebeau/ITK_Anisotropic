//
//  TubularBand.hxx
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 01/09/13.
//
//

#ifndef ITKFM_TubularBand_hxx
#define ITKFM_TubularBand_hxx

#include <itkNumericTraits.h>
#include <itkMath.h>
//#include "IsoFM_Geodesic.h"

#include "itkUnaryFunctorImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkTernaryFunctorImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "RescaleAndExport_Function.h"

namespace itk {
    
    // *********** Core functionality : fast marching from a curve *********

    template<typename TLS>
    void CurveNeighborhoodSource<TLS>::SetOutputRegionContainingPaths(ComponentType marginInPixels){
        ContinuousIndexType bottomLeft, topRight;
        bottomLeft.Fill(NumericTraits<ComponentType>::max());
        topRight.Fill(NumericTraits<ComponentType>::min());
        for(const auto & path : paths){
            for(const auto & p : path){
                for(int i=0; i<Dimension; ++i){
                    bottomLeft[i]=std::min(bottomLeft[i],p[i]);
                    topRight[i] = std::max(topRight[i],p[i]);
                }
            }
        }
        if(bottomLeft[0]==NumericTraits<ComponentType>::max())
            itkExceptionMacro("SetOutputRegionContainingPaths : Error no path point");
        
        for(int i=0; i<Dimension; ++i){
            bottomLeft[i]-=marginInPixels;
            topRight[i]+=marginInPixels;
        }
        
        IndexType index;
        Size<Dimension> size;
        for(int i=0; i<Dimension; ++i){
            index[i] = itk::Math::Floor<IndexValueType>(bottomLeft[i]);
            size[i] = itk::Math::Ceil<IndexValueType>(topRight[i]) - index[i];
        }
        this->SetOutputRegion(ImageRegion<Dimension>(index,size));
    }

    template<typename TLS>
    const typename CurveNeighborhoodSource<TLS>::PathType &
    CurveNeighborhoodSource<TLS>::GetPath(size_t i) const {
        if(i>=paths.size())
            itkExceptionMacro("GetPath error: index " << i << " exceeds number of paths " << paths.size() << ".");
        return paths[i];
    }
    
    template<typename TLS>
    size_t
    CurveNeighborhoodSource<TLS>::AddPhysicalPath(const PhysicalPathType & physicalPath){
        this->GenerateOutputInformation();
        const LevelSetImageType * output = this->GetOutput();
        paths.resize(paths.size()+1);
        PathType & path = paths.back();
        
        for(const auto & point : physicalPath){
            ContinuousIndexType continuousIndex;
            output->TransformPhysicalPointToContinuousIndex(point,continuousIndex);
            path.push_back(continuousIndex);
        }
        return paths.size()-1;
    }
    

    template<typename TLS>
    void
    CurveNeighborhoodSource<TLS>::GenerateData(){
        std::cout << "Fast Marching" << std::endl;
        this->GenerateOutputInformation();
        const LevelSetImageType * output = this->GetOutput();
        
        // Prepare for the seeds.
        typename SuperClass::NodeContainerPointer seeds = SuperClass::NodeContainer::New();
        seeds->Initialize();
        size_t counter=0;
        
        // Set the seeds.
        IndexSetType indexSet;
        for(const auto & path : paths){
            for(int i=0; (i+1)<(int)path.size(); ++i){
                indexSet.clear();
                SegmentPoints(path[i], path[i+1], indexSet);
                PointType p,q;
                output->TransformContinuousIndexToPhysicalPoint(path[i],p);
                output->TransformContinuousIndexToPhysicalPoint(path[i+1],q);
                VectorType v = q-p;
                const ScalarType vNorm = v.GetNorm();
                if(vNorm==0) continue;
                else v/=vNorm;
                
                for(const IndexType & index : indexSet){
                    PointType r;
                    output->TransformIndexToPhysicalPoint(index,r);
                    const VectorType w = r-p;
                    
                    typename SuperClass::NodeType node;
                    node.SetIndex(index);
                    const ScalarType squaredDistanceToSegment = w.GetSquaredNorm() - vnl_math_sqr(w*v);
                    node.SetValue( sqrt(std::max(ScalarType(0.),squaredDistanceToSegment)));
                    seeds->InsertElement(counter++,node);
                } // for index around segment
            } // for segment
        } // for path
        
        this->SetTrialPoints(seeds);
        this->SetGenerateGradientImage(true);
        orderedIndices.clear();
        SuperClass::GenerateData();
    }
    
    // ************ Compute various images ************
    
    // At right predicate and signed distance, in 2D.
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::ComputeAtRightFuzzyPredicate() const {
        static_assert(Dimension==2,"ComputeAtRightFuzzyPredicate : The right of a curve is only defined in dimension 2.");
        const LevelSetImageType * output = this->GetOutput();
        
        ScalarImagePointerType atRightFuzzyPredicate = ScalarImageType::New();
        atRightFuzzyPredicate->CopyInformation(output); // Graft changed in version 4.4, doesn't work here.
        atRightFuzzyPredicate->SetRegions(GetRegion());
        atRightFuzzyPredicate->Allocate();
        atRightFuzzyPredicate->FillBuffer(ScalarType());
                
        IndexSetType indexSet;
        for(const auto & path : paths){
            for(int i=0; (i+1)<(int)path.size(); ++i){
                indexSet.clear();
                SegmentPoints(path[i], path[i+1], indexSet);
                PointType p,q;
                output->TransformContinuousIndexToPhysicalPoint(path[i],p);
                output->TransformContinuousIndexToPhysicalPoint(path[i+1],q);
                VectorType v = q-p;
                
                for(const IndexType & index : indexSet){
                    PointType r;
                    output->TransformIndexToPhysicalPoint(index,r);
                    const VectorType w = r-p;
                    
                    const ScalarType determinant2D = w[0]*v[1]-w[1]*v[0];
                    const ScalarType atRight = determinant2D >= 0 ? 1 : -1;
                    
                    atRightFuzzyPredicate->GetPixel(index) = atRight;
                } // for index around segment
            } // for segment
        } // for path
        
        ExtendDataToNeighborhood(atRightFuzzyPredicate.GetPointer());
        return atRightFuzzyPredicate;
    }
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::BooleanImageType::Pointer
    CurveNeighborhoodSource<TLS>::ComputeAtRightPredicate() const {
        typedef Functor::BinaryThreshold<ScalarType, bool> FunctorType;
        typedef UnaryFunctorImageFilter<ScalarImageType, BooleanImageType, FunctorType> ImageFunctorType;
        typename ImageFunctorType::Pointer imageFunctor = ImageFunctorType::New();
        FunctorType & functor = imageFunctor->GetFunctor();
        functor.SetInsideValue(true);
        functor.SetOutsideValue(false);
        functor.SetLowerThreshold(0);
        functor.SetUpperThreshold(NumericTraits<ScalarType>::max());
        
        imageFunctor->SetInput(ComputeAtRightFuzzyPredicate());
        imageFunctor->Update();
        return imageFunctor->GetOutput();
    }

    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::ComputeSignedDistance() const {
        
        struct FunctorType {
            bool operator!=(const FunctorType &) const {return false;}
            bool operator==(const FunctorType & other) const {return !( *this != other );}
            
            ScalarType operator()(ScalarType distance, ScalarType atRight){
                return atRight>=0 ? distance : -distance;
            }
        };
        
        typedef BinaryFunctorImageFilter<ScalarImageType, ScalarImageType, ScalarImageType, FunctorType> BinFunctorType;
        typename BinFunctorType::Pointer functor = BinFunctorType::New();
        functor->SetInput2(ComputeAtRightFuzzyPredicate());
        functor->SetInput1(this->GetOutput());
        functor->Update();
        return functor->GetOutput();
    }
    
    // Curvilinear coordinates and tangent vectors
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::ComputeCurvilinearCoordinates() const {
        const LevelSetImageType * output = this->GetOutput();
        
        std::vector<std::vector<ScalarType> > pathLengths;
        pathLengths.reserve(paths.size());
        ScalarType length = 0;
        for(const auto & path : paths){
            pathLengths.resize(pathLengths.size()+1);
            if(path.size()==0) continue;
            pathLengths.back().push_back(length);
            for(int i=0; (i+1)<(int)path.size(); ++i){
                PointType p,q;
                output->TransformContinuousIndexToPhysicalPoint(path[i],p);
                output->TransformContinuousIndexToPhysicalPoint(path[i+1],q);
                length += (p-q).GetNorm();
                pathLengths.back().push_back(length);
            } // for i
        } // for path
        return ExtendDataToNeighborhood<ScalarImageType>(pathLengths);
    }
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::VectorImageType::Pointer
    CurveNeighborhoodSource<TLS>::ComputeTangentVectors() const {
        LevelSetImageType * output = this->GetOutput();
        
        struct NormalizerType {
            bool operator!=(const NormalizerType &) const {return false;}
            bool operator==(const NormalizerType & other) const {return !( *this != other );}
            
            VectorType operator()(const VectorType & v){
                const ScalarType vNorm = v.GetNorm();
                if(vNorm>0) return v/vNorm;
                else return v;
            }

        } normalizer;
        
        std::vector<std::vector<VectorType> > pathTangents;
        pathTangents.reserve(paths.size());
        for(const auto & path : paths){
            pathTangents.resize(pathTangents.size()+1);
            if(path.size()==0) continue;
            if(path.size()==1)
                pathTangents.push_back(VectorType(0.));
            
            for(int i=0; i<path.size(); ++i){
                PointType p,q,r;
                if(i>0)             output->TransformContinuousIndexToPhysicalPoint(path[i-1],p);
                output->TransformContinuousIndexToPhysicalPoint(path[i],p);
                if(i+1<path.size()) output->TransformContinuousIndexToPhysicalPoint(path[i+1],q);
                
                VectorType tangents[3];
                if(i>0)             tangents[0] = normalizer(q-p);
                if(i+1<path.size()) tangents[1] = normalizer(r-q);
                tangents[2] = normalizer( (tangents[0]+tangents[1])/2. );
                
                
                if(i==0)
                    pathTangents.push_back(tangents[1]);
                else if (i+1==path.size())
                    pathTangents.push_back(tangents[0]);
                else pathTangents.push_Back(tangents[2]);
            } // for i
        } // for path
        
        
        typedef UnaryFunctorImageFilter<VectorImageType, VectorImageType, NormalizerType> FunctorType;
        typename FunctorType::Pointer normalizerFilter = FunctorType::New();
        normalizerFilter->SetInput(ExtendDataToNeighborhood(pathTangents));
        normalizerFilter->Update();
        return normalizerFilter->GetOutput();
    }
    
    // **************** Indicator functions **************
    
    // Curve interior
    
    namespace {
        template<typename ScalarType>
        struct SmoothTranslatedHeavisideFunctorType {
            ScalarType tubeWidth, boundaryWidth;
            inline ScalarType operator()(ScalarType s) const {return 0.5 + atan2(tubeWidth-s, boundaryWidth)/itk::Math::pi;}
        };
        
        template<typename ScalarType>
        struct SmoothTranslatedHeavisideTernaryFunctorType {
            inline ScalarType operator()(ScalarType s, ScalarType tubeWidth, ScalarType boundaryWidth) const {return 0.5 + atan2(tubeWidth-s, boundaryWidth)/itk::Math::pi;}
        };        
    }
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::InteriorSmoothIndicatorFunction(ScalarType width) const {
        typedef UnaryFunctorImageFilter<ScalarImageType, ScalarImageType, SmoothTranslatedHeavisideFunctorType<ScalarType> > FunctorType;
        typename FunctorType::Pointer functor = FunctorType::New();
        functor->GetFunctor().tubeWidth=0;
        functor->GetFunctor().boundaryWidth=width;
        functor->SetInput(ComputeSignedDistance());
        functor->Update();
        return functor->GetOutput();
    }
    
    // Curve neighborhood
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::TubeSmoothIndicatorFunction(ScalarType tubeWidth, ScalarType boundaryWidth) const {
        typedef UnaryFunctorImageFilter<ScalarImageType, ScalarImageType, SmoothTranslatedHeavisideFunctorType<ScalarType> > FunctorType;
        typename FunctorType::Pointer functor = FunctorType::New();
        functor->GetFunctor().tubeWidth=tubeWidth;
        functor->GetFunctor().boundaryWidth=boundaryWidth;
        functor->SetInput(this->GetOutput());
        functor->Update();
        return functor->GetOutput();
    }
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::TubeSmoothIndicatorFunction(const std::vector<ScalarType> & tubeWidths, const std::vector<ScalarType> & boundaryWidths) const {
        std::vector<std::vector<ScalarType> > tubeWidths_, boundaryWidths_;
        tubeWidths_.push_back(tubeWidths);
        boundaryWidths_.push_back(boundaryWidths);
        return TubeSmoothIndicatorFunction(tubeWidths_,boundaryWidths_);
    }
    
    template<typename TLS>
    typename CurveNeighborhoodSource<TLS>::ScalarImagePointerType
    CurveNeighborhoodSource<TLS>::TubeSmoothIndicatorFunction(const std::vector<std::vector<ScalarType> > & tubeWidths, const std::vector<std::vector<ScalarType> > & boundaryWidths) const {
        if(tubeWidths.size() != paths.size())
            itkExceptionMacro("TubeSmoothIndicatorFunction error : number of paths " << paths.size() << ", does not coincide with tube width data " << tubeWidths.size() << ".");
        if(boundaryWidths.size() != paths.size())
            itkExceptionMacro("TubeSmoothIndicatorFunction error : number of paths " << paths.size() << ", does not coincide with boundary width data " << boundaryWidths.size() << ".");
                              
        for(size_t i=0; i<paths.size(); ++i){
            if(paths[i].size() != tubeWidths[i].size())
                itkExceptionMacro("TubeSmoothIndicatorFunction error : length of path" << i << ", namely " << paths[i].size() << "does not coincide with tube widths data length for this path " << tubeWidths[i].size());
            if(paths[i].size() != boundaryWidths[i].size())
                itkExceptionMacro("TubeSmoothIndicatorFunction error : length of path" << i << ", namely " << paths[i].size() << "does not coincide with tube widths data length for this path " << boundaryWidths[i].size());
        }
        
        auto functorFilter = itk::TernaryFunctorImageFilter<ScalarImageType, ScalarImageType, ScalarImageType, ScalarImageType, SmoothTranslatedHeavisideTernaryFunctorType<ScalarType> >::New();
        functorFilter->SetInput1(this->GetOutput());
        functorFilter->SetInput2(ExtendDataToNeighborhood<ScalarImageType>(tubeWidths));
        functorFilter->SetInput3(ExtendDataToNeighborhood<ScalarImageType>(boundaryWidths));
        functorFilter->Update();
        return functorFilter->GetOutput();
    }

    
    // ************ Extend data to neighborhood ***********
    
    template<typename TLS>
    template<typename DataImageType>
    typename DataImageType::Pointer
    CurveNeighborhoodSource<TLS>::ExtendDataToNeighborhood(const std::vector<std::vector<typename DataImageType::PixelType> > & datas) const {
        if(datas.size() != paths.size())
            itkExceptionMacro("ExtendDataToNeighborhood error : size of paths and of datas are inconsistent.");
        typedef typename DataImageType::PixelType DataType;
        
        const LevelSetImageType * output = this->GetOutput();

        typedef typename DataImageType::Pointer DataImagePointerType;
        DataImagePointerType dataImage = DataImageType::New();
        dataImage->CopyInformation(output);
        dataImage->SetRegions(GetRegion());
        dataImage->Allocate();
        dataImage->FillBuffer(DataType());
        
        IndexSetType indexSet;
        for(size_t counter = 0; counter<paths.size(); ++counter){
            const PathType & path = paths[counter];
            const std::vector<DataType> & data = datas[counter];
            if(path.size() != data.size())
                itkExceptionMacro("ExtendDataToNeighborhood error : size of path and of data is inconsistent, at position " << counter << ".");
            
            for(int i=0; (i+1)<(int)path.size(); ++i){
                indexSet.clear();
                SegmentPoints(path[i], path[i+1], indexSet);
                PointType p,q;
                output->TransformContinuousIndexToPhysicalPoint(path[i],p);
                output->TransformContinuousIndexToPhysicalPoint(path[i+1],q);
                VectorType v = q-p;
                const ScalarType vSquaredNorm = v.GetSquaredNorm(); // Note that norm is squared here
                if(vSquaredNorm==0) continue;
                else v/=vSquaredNorm;
                
                for(const IndexType & index : indexSet){
                    PointType r;
                    output->TransformIndexToPhysicalPoint(index,r);
                    const VectorType w = r-p;
                    ScalarType t = v*w;
                    t = std::max(ScalarType(0.),std::min(ScalarType(1.),t));
                    
                    dataImage->GetPixel(index) = (1-t)*data[i]+t*data[i+1];
                } // for index around segment
            } // for segment
        } // for path
        
        typedef typename DataImageType::PixelType DataType;
        ExtendDataToNeighborhood<DataType>(dataImage); // To do : cast ITK Pointer to standard pointer
        return dataImage;
    }
    
    template<typename TLS>
    template<typename DataType>
    void
    CurveNeighborhoodSource<TLS>::ExtendDataToNeighborhood(Image<DataType,Dimension> * dataImage) const {
        typename SuperClass::GradientImagePointer gradientImage = this->GetGradientImage();
        typedef itk::Matrix<double,Dimension,Dimension> MatrixTypeD;
        const MatrixTypeD DirectionD = gradientImage->GetDirection();
        const MatrixTypeD DualDirectionD = MatrixTypeD(MatrixTypeD(DirectionD.GetTranspose()).GetInverse());
        typedef itk::Matrix<ScalarType, Dimension, Dimension> MatrixType;
        
        MatrixType DualDirection;
        for(int i=0; i<(int)Dimension; ++i)
            for(int j=0; j<(int)Dimension; ++j)
                DualDirection(i,j) = DualDirectionD(i,j);
                
        for(const IndexType & index : orderedIndices){
            const auto gradient = gradientImage->GetPixel(index);            
            const auto signedWeights = - (DualDirection*gradient); // To do : cast direction matric to avoid ambiguity
            
            ScalarType l1Norm = 0;
            for(int i=0; i<Dimension; ++i)
                l1Norm += fabs(signedWeights[i]);
            
            if(l1Norm==0) continue; // Case of a source point
            
            DataType & data = dataImage->GetPixel(index);
            if(data != 0) continue;
            data=0;
            
            for(int i=0; i<Dimension; ++i){
                IndexType indexNeigh = index;
                const ScalarType s = signedWeights[i]/l1Norm;
                if(s==0) continue;
                
                indexNeigh[i] += (s>0 ? 1 : -1);
                data += fabs(s)*dataImage->GetPixel(indexNeigh);
            }
        }
    }
    
    // *********** Segment points *********
    
    template<typename TLS>
    void
    CurveNeighborhoodSource<TLS>::SegmentPoints(const ContinuousIndexType & p, const ContinuousIndexType & q, IndexSetType & indexSet) const{
        
        typedef VectorType ContinuousOffsetType;
        ContinuousOffsetType v;
        for(int i=0; i<Dimension; ++i)
            v[i] = q[i]-p[i];
        
        for(int i=0; i<Dimension; ++i){
            const IndexValueType
            RiMin = itk::Math::Ceil<ScalarType>(std::min(p[i],q[i])),
            RiMax = itk::Math::Ceil<ScalarType>(std::max(p[i],q[i]));
            
            for(IndexValueType Ri = RiMin; Ri<RiMax; ++Ri){
                const double t = (ComponentType(Ri) - p[i]) / (q[i]-p[i]);
                ContinuousIndexType r;
                r.SetToBarycentricCombination(q,p,t);
                
                IndexType RMin, RMax;
                RMin[i]=Ri;
                RMax[i]=Ri;
                for(int j=0; j<Dimension; ++j){
                    if(j==i) continue;
                    RMin[j]=itk::Math::Ceil<IndexValueType>(r[j])-1;
                    RMax[j]=itk::Math::Floor<IndexValueType>(r[j])+1;
                }
                
                IndexType R = RMin;
                unsigned int j=0;
                while(j!=Dimension){
                    if(GetRegion().IsInside(R))
                        indexSet.insert(R);
                    for(j=0; j<Dimension; ++j){
                        if(R[j] < RMax[j]){
                            ++R[j];
                            break;
                        } else
                            R[j]=RMin[j];
                    } // for j
                } // while j
                
            } // for Ri
        } // for i dimension
    }
    
} // namespace itk
#endif
