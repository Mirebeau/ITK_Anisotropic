//
//  GenerateTestImages.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 16/09/13.
//
//

#ifndef ITKFM_GenerateTestImages_h
#define ITKFM_GenerateTestImages_h

#include "CurveNeighborhoodSource.h"
#include "RescaleAndExport_Function.h"

void CurveNeighborhoodSource_Test(){
    using std::cout;
    using std::endl;
    
    // Testing tubular source.
    
    typedef float ScalarType;
    const unsigned int Dimension=2;
    typedef itk::Image<ScalarType,Dimension> ImageType;
    typedef itk::Point<ScalarType,Dimension> PointType;


    // Creating image source
    typedef itk::CurveNeighborhoodSource<ImageType> ImageSourceType;
    auto imageSource = ImageSourceType::New();
    
    const ScalarType h = 1/10.;
    ImageType::SpacingType spacing;
    spacing.Fill(h);
    imageSource->SetOutputSpacing(spacing);
    
    PointType origin;
    origin.Fill(0);
    imageSource->SetOutputOrigin(origin);
    
    
    std::vector<PointType> physicalPath;
    for(int i=0; i<2*itk::Math::pi/h; ++i){
        ScalarType t = h*i;
        PointType p;
        p[0] = cos(t);
        p[1] = sin(t);
        physicalPath.push_back(p);
    }
    physicalPath.push_back(physicalPath[0]);
    
    imageSource->AddPhysicalPath(physicalPath);
    imageSource->SetOutputRegionContainingPaths();
    imageSource->Update();
    
    std::cout << "Hi there" << imageSource->GetOutput() << std::endl;
    RescaleAndExport(imageSource->GetOutput(),"DistanceToCircle.bmp");
    

    
    std::cout << "Stopping radius : " << imageSource->GetStoppingRadius() << std::endl;
    auto atRightFuzzyPredicate = imageSource->ComputeAtRightFuzzyPredicate();
    RescaleAndExport(atRightFuzzyPredicate.GetPointer(), "AtRightFuzzyPredicate.bmp");

    /*
    struct Testing : public ImageSourceType {
        typedef ImageSourceType::IndexSetType IndexSetType;
        void SegmentPoints(ContinuousIndexType p, ContinuousIndexType q, IndexSetType & indexSet){
            ImageSourceType::SegmentPoints(p,q,indexSet);
        }
    } tester;
    
    const auto & path = imageSource->GetPath(0);
    for(int i=0; i<path.size()-1; ++i){
    
        Testing::IndexSetType indexSet;
        tester.SegmentPoints(path[i],path[i+1],indexSet);
    
        cout << "Extremities : " << path[i] << "," << path[i+1] << endl;
        for(auto index : indexSet)
            cout << index << endl;
    }
    */

    auto signedDistance = imageSource->ComputeSignedDistance();
    RescaleAndExport(signedDistance.GetPointer(), "SignedDistanceToCircle.bmp");
    
    auto atRightPredicate = imageSource->ComputeAtRightPredicate();
    RescaleAndExport(atRightPredicate.GetPointer(), "AtRightPredicate.bmp");

    auto curvilinearCoordinates = imageSource->ComputeCurvilinearCoordinates();
    RescaleAndExport(curvilinearCoordinates.GetPointer(), "CurvilinearCoordinates.bmp");
}




#endif
