//
//  TubularBand.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 01/09/13.
//
//

#ifndef ITKFM_TubularBand_h
#define ITKFM_TubularBand_h

#include <vector>
#include <set>
#include "itkFastMarchingUpwindGradientImageFilter.h"

/* 
 CurveNeighborhood source.
 Computes, using isotropic fast marching, the distance to a curve, and extends to the neighborhood quantities of interest, such as tangent vectors, left and right orientation, curve parametrization or any values on the curve.
 
 LevelSet image type must have floating point pixel type. This is used as component type for vectors, etc.
 
 Transformation from Image coordinates to Physical points and vectors, and conversely, is computed using using user specified Origin and Spacing. (Be sure to set SetOutputSpacing() and SetOutputOrigin() before using these, otherwise  default 1 to 1 mapping is used)
*/

namespace itk {
    
template<typename TLevelSet>
class CurveNeighborhoodSource :
public FastMarchingUpwindGradientImageFilter<TLevelSet>
{
public:
    using Self = CurveNeighborhoodSource;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;
    typedef FastMarchingUpwindGradientImageFilter<TLevelSet> SuperClass;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(Self, ImageSource);

    typedef TLevelSet LevelSetImageType;
    typedef LevelSetImageType OutputImageType;
    static const int Dimension = LevelSetImageType::ImageDimension;

    typedef LevelSetImageType ScalarImageType;
    typedef typename ScalarImageType::Pointer ScalarImagePointerType;
    typedef typename ScalarImageType::PixelType ScalarType;
    
    typedef ScalarType ComponentType;
    typedef ContinuousIndex<ComponentType,Dimension> ContinuousIndexType;
    typedef Vector<ComponentType,Dimension> VectorType;
    typedef Point<ComponentType,Dimension> PointType;
    
    typedef Index<Dimension> IndexType;
    typedef Image<bool,Dimension> BooleanImageType;
    
    // ********** Additional possible outputs **********
    
    // ComputeSignedDistance, ComputeAtRightPredicate and ComputeAtRightFuzzyPredicate obviously only work in 2D.
    
    ScalarImagePointerType
    ComputeSignedDistance() const;
    
    typename BooleanImageType::Pointer
    ComputeAtRightPredicate() const;
    
    ScalarImagePointerType
    ComputeAtRightFuzzyPredicate() const;
    
    // Data extension, and immediate specializations
    template<typename DataImageType>
    typename DataImageType::Pointer
    ExtendDataToNeighborhood(const std::vector<std::vector<typename DataImageType::PixelType> > &) const;

    ScalarImagePointerType
    ComputeCurvilinearCoordinates() const;
    
    typedef Image<VectorType,Dimension> VectorImageType;
    typename VectorImageType::Pointer
    ComputeTangentVectors() const;
    
    // ********* Smoothed indicator functions *********
    
    /// Smoothed indicator function of the interior (left of curve, in 2D). Parameter is the transition physical width.
    ScalarImagePointerType
    InteriorSmoothIndicatorFunction(ScalarType) const;
    
//    ScalarImagePointerType
//    InteriorSmoothIndicatorFunction(const std::vector<std::vector<ScalarType> > &) const;
    
    /// Smoothed indicator function of the neighborhood. Parameters are the tube radius and the boundary transition width.
    ScalarImagePointerType
    TubeSmoothIndicatorFunction(ScalarType, ScalarType) const;

    // Radius may vary along the curve(s)
    ScalarImagePointerType
    TubeSmoothIndicatorFunction(const std::vector<ScalarType> &, const std::vector<ScalarType> &) const;
    ScalarImagePointerType
    TubeSmoothIndicatorFunction(const std::vector<std::vector<ScalarType> > &, const std::vector<std::vector<ScalarType> > &) const;
    
    
    // ********** Path data ***********
    
    typedef std::vector<ContinuousIndexType> PathType;
    
    const PathType & GetPath(size_t i = 0) const;
    size_t AddPath(const PathType & path){paths.push_back(path); return paths.size()-1;}
    void ClearPaths(){paths.clear();}

    // Do not forget SetOutputOrigin and SetOutputSpacing before any use of physical coordinates.
    typedef std::vector<PointType> PhysicalPathType;
    size_t AddPhysicalPath(const PhysicalPathType &);
    
    void SetOutputRegionContainingPaths(ComponentType marginInPixels=10);
    
    void SetStoppingRadius(ScalarType r){m_stoppingRadius=r;}
    ScalarType GetStoppingRadius() const {return m_stoppingRadius;}
    
protected:
    std::vector<PathType> paths;
    
//    virtual void GenerateOutputInformation();
//    virtual void Initialize(){};
    virtual void GenerateData() override;
    const ImageRegion<Dimension> & GetRegion() const {return this->GetOutput()->GetRequestedRegion();}
    
    typedef std::set<IndexType, typename IndexType::LexicographicCompare> IndexSetType;
    void SegmentPoints(const ContinuousIndexType &, const ContinuousIndexType &, IndexSetType &) const;
    
    template<typename DataType>
    void
    ExtendDataToNeighborhood(Image<DataType,Dimension> *) const;
    
    
    // Early abort support, and index ordering
    std::vector<IndexType> orderedIndices;
    ScalarType m_stoppingRadius = std::numeric_limits<ScalarType>::max();
    virtual void UpdateNeighbors(const IndexType & index,
                                 const typename SuperClass::SpeedImageType * speed,
                                 LevelSetImageType * output) override {
        if(output->GetPixel(index) > this->GetStoppingRadius()) return;
        orderedIndices.push_back(index);
        SuperClass::UpdateNeighbors(index,speed,output);
    }
};

    
}

#include "CurveNeighborhoodSource.hxx"
    


#endif
