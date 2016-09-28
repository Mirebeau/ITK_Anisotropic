//
//  IllustrationCode.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 21/05/2015.
//
//

#ifndef ITKFM_IllustrationCode_h
#define ITKFM_IllustrationCode_h

void HelloWorld(){
    
    // 1. Declare a norm type, and construct a metric.
    typedef itk::Riemannian2DNorm<float> NormType; // Alternative types
    typedef itk::Image<NormType, NormType::Dimension> MetricType;
    MetricType::Pointer metric = MetricType::New();
    /* TO DO : Resize, Allocate and Fill the metric as desired */
    
    // 2. Declare and initialize the Fast Marching object
    typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
    
    // Front origin is defined in a NodeContainer, inherited from itk::FastMarchingImageFilter
    FMType::NodeType node;
    FMType::IndexType index; index.Fill(0); const double value=0;
    node.SetValue(value); // value and index are user defined.
    node.SetIndex(index);
    
    FMType::NodeContainer::Pointer seeds = FMType::NodeContainer::New();
    seeds->Initialize();
    seeds->InsertElement(0,node); // Several seeds may be inserted
    
    // Fast marching filter
    FMType::Pointer fm = FMType::New();
    fm->SetGenerateUpwindGradient(true); // Required to extract geodesics
    fm->SetInput(metric);
    fm->SetTrialPoints(seeds);
    fm->Update();
    // fm->GetOutput() contains the front arrival times.

    // 3. Extract desired geodesics
    FMType::PointType start; /* TO DO : Choose geodesic start point coordinates */
    std::vector<FMType::PointType> geodesic;
    geodesic.push_back(start);
    fm->Geodesic(geodesic);
    
}

#endif
