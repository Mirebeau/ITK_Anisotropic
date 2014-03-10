//
//  AnisotropicFastMarching.h
//  MatlabITKFastMarching
//
//  Created by Jean-Marie Mirebeau on 04/04/13.
//
//

#ifndef MatlabITKFastMarching_AnisotropicFastMarching_h
#define MatlabITKFastMarching_AnisotropicFastMarching_h

#include "itkAnisotropicFastMarchingImageFilter.h"
#include "AnisotropicFastMarching_EarlyAbortClass.h"

#include "MexMessageWrapper.h"
#include "mxIO.h"

template<typename NormType>
int AnisotropicFastMarching_EarlyAbort(mxIO & IO){

    const int Dimension = NormType::Dimension;
    typedef double ValueType;

    typedef itk::Image<ValueType,Dimension> LevelSetType;
    typedef itk::Image<NormType,Dimension> MetricImageType;
    typedef itk::AnisotropicFastMarchingImageFilter_EarlyAbort<NormType> FMType;

    typedef itk::Point<ValueType,Dimension> PointType;
    typedef typename MetricImageType::SpacingType SpacingType;
    typedef std::pair<PointType,ValueType> SeedType;
    typedef itk::Index<Dimension> IndexType;
    
    // Get Metric, Seeds, Tips, and some options
    
    IO.transposeFirstTwoCoordinates = IO.GetObject<double>("TransposeFirstTwoImageCoordinates",1.);
    
    auto Metric = IO.GetImage<MetricImageType>("Metric");
    Metric->SetOrigin(IO.GetObject<PointType>("Origin", PointType(0.)));
    Metric->SetSpacing(IO.GetObject<SpacingType>("Spacing",SpacingType(1.)));

    const bool checkData = IO.GetObject<double>("CheckData", 1.);
    if(checkData){
        itk::ImageRegionConstIteratorWithIndex<MetricImageType> it(Metric, Metric->GetBufferedRegion());
        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
            if(!(it.Value().IsDefinite())){
                MexWarnMsg() << "Invalid data !" << it.Value() << " at position " << it.GetIndex() << ".";
                return EXIT_FAILURE;
            } // if not definite
    }
    
    // Create the fast marching
    auto FM = FMType::New();
    FM->SetInput(Metric);
    FM->SetOutputRegion(Metric->GetBufferedRegion());
    FM->SetGenerateUpwindGradient(true);
    

    // Vijaya's options
    
    if(IO.HasField("StopAtDistance"))
        FM->StopAtDistance(IO.GetObject<double>("StopAtDistance"));
    
    if(IO.GetObject<double>("StopWhenTipsAreReached",0.))
        for(auto Tip : IO.GetVector<PointType>("Tips")){
            itk::Index<Dimension> TipIndex;
            Metric->TransformPhysicalPointToIndex(Tip,TipIndex);
            FM->StopWhenLastIsAccepted(TipIndex);
        } // for it tips
    
    if(IO.HasField("StopWhenFirstIsAccepted")){
        auto stopWhenFirstIsAcceptedVector = IO.GetVector<PointType>("StopWhenFirstIsAccepted");
        for(const PointType & stopIfAcceptedPoint : stopWhenFirstIsAcceptedVector){
            IndexType stopIfAcceptedIndex;
            if(!Metric->TransformPhysicalPointToIndex(stopIfAcceptedPoint, stopIfAcceptedIndex)) continue;
            FM->StopWhenFirstIsAccepted(stopIfAcceptedIndex);
        }
    }
    
    // Chenda's options
    if(IO.HasField("StopAtEuclideanDistance")) {
        FM->StopAtEuclideanDistance(IO.GetObject<double>("StopAtEuclideanDistance"));
        FM->SetEuclideanSpacing(IO.GetObject<SpacingType>("EuclideanSpacing",Metric->GetSpacing() ) );
    }
    
    // Set the seeds
    {
        auto SeedsVector = IO.GetVector<SeedType>("Seeds");
        auto SeedsContainer = FMType::NodeContainer::New();
        SeedsContainer->Initialize();
        
        int counter=0;
        for(auto Seed: SeedsVector){
            typename FMType::NodeType Seed_;
            Seed_.SetValue(Seed.second);
            if( ! Metric->TransformPhysicalPointToIndex(Seed.first, Seed_.GetIndex()) )
                MexWarnMsg() << "Seed " << Seed.first << ", of index " << counter << " is out of bounds. Hence ignored.";
            else SeedsContainer->InsertElement(counter++, Seed_);
        }
        FM->SetTrialPoints(SeedsContainer);
    }
    
    // Run fast marching
    
    clock_t timing=-clock();
    FM->Update();
    timing+=clock();
    MexMsg() << "Anisotropic Fast Marching took : " << timing/double(CLOCKS_PER_SEC) << " seconds.\n";
    
    // Export Distance map
    IO.SetField("Distance", IO.mxImage<LevelSetType>(FM->GetOutput()));
    
    // Export euclidean distance
    if(IO.HasField("StopAtEuclideanDistance"))
        IO.SetField("EuclideanPathLengths", IO.mxImage<LevelSetType>(FM->GetEuclideanPathLengths()));
    
    std::string activeStoppingCriteria;
    switch (FM->GetActiveStoppingCriteria()) {
        case FMType::None:                  activeStoppingCriteria="None"; break;
        case FMType::Distance:              activeStoppingCriteria="StopAtDistance"; break;
        case FMType::EuclideanDistance:     activeStoppingCriteria="StopAtEuclideanDistance"; break;
        case FMType::LastAccepted:          activeStoppingCriteria="StopWhenLastIsAccepted"; break;
        case FMType::FirstAccepted:         activeStoppingCriteria="StopWhenFirstIsAccepted"; break;
        default: itkGenericExceptionMacro("Unhandled stopping criteria")
    }
    IO.SetField("ActiveStoppingCriteria", mxCreateString(activeStoppingCriteria.c_str()));
    
    {
        PointType stoppingPoint;
        Metric->TransformIndexToPhysicalPoint(FM->GetStoppingIndex(),stoppingPoint);
//        IO.SetField("StoppingPoint", IO.mxObject(stoppingPoint));
        std::vector<PointType> geodesic;
        geodesic.push_back(stoppingPoint);
        FM->Geodesic(geodesic);
        IO.SetField("GeodesicFromStoppingPoint", IO.mxVector(geodesic));
    }
    
    // Compute the paths
    if(IO.HasField("Tips")){
        auto Tips = IO.GetVector<PointType>("Tips");
        const mwSize nTips = Tips.size();
        mxArray * MatlabGeodesics = mxCreateCellArray(1,&nTips);
    
        int counter=0;
        for(auto Tip : Tips){
            std::vector<PointType> Geodesic;
            Geodesic.push_back(Tip);
            FM->Geodesic(Geodesic);
            mxSetCell(MatlabGeodesics, counter++, IO.mxVector<PointType>(Geodesic));
        }
        IO.SetField("Geodesics", MatlabGeodesics);
    }
    return EXIT_SUCCESS;
}

#endif
