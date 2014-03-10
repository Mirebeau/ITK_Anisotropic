//
//  TubularBand_Test.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 03/09/13.
//
//

#ifndef ITKFM_TubularBand_Test_h
#define ITKFM_TubularBand_Test_h

#include <ostream>
#include <algorithm>

#include "itkImageRegionIterator.h"

#include "itkAnisotropicFastMarchingImageFilter.h"
#include "Riemannian2DNorm.h"
#include "Riemannian3DNorm.h"
#include "Finsler2DNorm.h"

#include "TubularBandSource.h"

#include "ExportToMathematica.h"

template<typename ComponentType>
struct VectorToMetricCaster {
    void SetKappa(ComponentType kappa){
        if(kappa<=0) itkGenericExceptionMacro("VectorToMetricCaster error : kappa should be positive");
        delta = ComponentType(1)/kappa;
    }
    
    template<typename VectorType, typename NormType>
    void operator()(const VectorType & v, NormType & m) const {
        assert(v.GetNorm()<1.001);
        Convert(v,m);
        //assert(m.IsDefinite());
    }
    
    template<typename MetricType, typename VectorFieldType, typename Caster=VectorToMetricCaster>
    typename MetricType::Pointer
    MakeMetric(const VectorFieldType * vectorField, Caster * caster = nullptr){
        if(caster==nullptr) caster = this;
        
        auto region = vectorField->GetBufferedRegion();
        
        auto metric = MetricType::New();
        metric->SetRegions(region);
        metric->SetSpacing(vectorField->GetSpacing());
        metric->SetOrigin(vectorField->GetOrigin());
        metric->Allocate();
        
        itk::ImageRegionIterator<MetricType> metricIt(metric,region);
        itk::ImageRegionConstIterator<VectorFieldType> vectorFieldIt(vectorField,region);
        
        for(metricIt.GoToBegin(), vectorFieldIt.GoToBegin(); ! metricIt.IsAtEnd(); ++metricIt, ++vectorFieldIt)
            caster->operator()(vectorFieldIt.Value(), metricIt.Value());
        
        return metric;
    }
    
protected:
    ComponentType delta = 0.1;
    
    template<unsigned int Dimension>
    void Convert(const itk::Vector<ComponentType,Dimension> & v, ComponentType & m) const {
        if(v.GetNorm()==0) m = 1;
        else m=1/delta;
    }
    
    void Convert(const itk::Vector<ComponentType,2> & v, itk::Riemannian2DNorm<ComponentType> & m) const {
        const int Dimension = 2;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j) = (i==j) + (delta*delta - 1)*v[i]*v[j];
    }
    
    void Convert(const itk::Vector<ComponentType,3> & v, itk::Riemannian3DNorm<ComponentType> & m) const {
        const int Dimension = 3;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j) = (i==j) + (delta*delta - 1)*v[i]*v[j];
    }
    
    void Convert(const itk::Vector<ComponentType,2> & v, itk::Finsler2DNorm<ComponentType> & f) const {
        const int Dimension=2;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                f.GetM()(i,j) = (i==j);
        f.GetOmega() = (1-delta)*v;
    }
};


int
TubularBand_2D_Test(){
    const std::string testPrefix = "TubularBand_2D_Test/";
    
    const unsigned int Dimension = 2;
    typedef float ComponentType;
    typedef itk::TubularBandSource<ComponentType,Dimension> TubularBandSource;
    typedef TubularBandSource::PointType PointType;
    typedef TubularBandSource::SpacingType SpacingType;
    typedef itk::Image<ComponentType, Dimension> ScalarImageType;
    
    // Define the curve
    /*
    struct PathType {
        int sampling = 100;
        SpacingType spacing = SpacingType(0.01);
        ComponentType radius=0.2, margin=0.4;
        PointType operator()(ComponentType t){
            t*=2*M_PI;
            PointType p;
            p[0]=cos(t); p[1]=sin(t);
            return p;
        }
        PointType geodesicTip, geodesicSeed;
        PathType(){geodesicTip = this->operator()(0); geodesicSeed = this->operator()(0.5);}
    } path;
    */
    
    struct PathType {
        int sampling = 100;
        SpacingType spacing = SpacingType(0.01);
        ComponentType radius=0.4, margin=0.4;
        PointType operator()(ComponentType t){
            t*=2*M_PI;
            PointType p;
            p[0]=2*cos(t); p[1]=sin(t)+0.5*sin(3*t);
            return p;
        }
        PointType geodesicTip, geodesicSeed;
        PathType(){geodesicTip = this->operator()(1./12); geodesicSeed = this->operator()(7./12);}
    } path;
    
    // TubularBandSource
    
    auto tubularBandSource = TubularBandSource::New();
    tubularBandSource->radius = path.radius;
    tubularBandSource->margin = path.margin;
    tubularBandSource->spacing = path.spacing;
    tubularBandSource->paths.resize(1);
    for(int i=0; i<path.sampling; ++i)
        tubularBandSource->paths.back().push_back(path((i+0.)/path.sampling));
    
    tubularBandSource->Update();
    auto vectorField = tubularBandSource->GetOutput();
    auto region = vectorField->GetBufferedRegion();
    RescaleAndExport(tubularBandSource->GetDistance(), testPrefix+"Euclidean_distance_to_path.bmp");
    
    auto vectorFieldSupport = ScalarImageType::New();
    vectorFieldSupport->SetRegions(region);
    vectorFieldSupport->Allocate();
    itk::ImageRegionConstIterator<TubularBandSource::OutputImageType> vectorFieldIt(vectorField,region);
    itk::ImageRegionIterator<ScalarImageType> vectorFieldSupportIt(vectorFieldSupport,region);
    for(vectorFieldIt.GoToBegin(), vectorFieldSupportIt.GoToBegin(); !vectorFieldIt.IsAtEnd(); ++vectorFieldIt, ++vectorFieldSupportIt)
        vectorFieldSupportIt.Value() = vectorFieldIt.Value().GetNorm();
    RescaleAndExport(vectorFieldSupport.GetPointer(), testPrefix+"Vector_Field_Support.bmp");
    
    // Create metrics
    
    VectorToMetricCaster<ComponentType> converter;
    converter.SetKappa(5);
    
    typedef ScalarImageType IsotropicMetricType;
    typedef itk::Image<itk::Riemannian2DNorm<ComponentType>, Dimension> RiemannianMetricType;
    typedef itk::Image<itk::Finsler2DNorm<ComponentType>, Dimension> FinslerMetricType;
    
    typedef itk::FastMarchingImageFilter<IsotropicMetricType,ScalarImageType>::NodeContainer NodeContainer;
    NodeContainer::Element node;
    vectorField->TransformPhysicalPointToIndex(path.geodesicSeed, node.GetIndex());
    node.SetValue(0);
    
    auto seeds = NodeContainer::New();
    seeds->Initialize();
    seeds->InsertElement(0,node);
    
    std::vector<PointType> geodesic;
    geodesic.push_back(path.geodesicTip);
    
    std::ofstream data_out;
    MathematicaExporterType exporter(data_out);
    
    
    {// Isotropic test
        const std::string prefix = testPrefix+"Isotropic/";

        auto metric = converter.MakeMetric<IsotropicMetricType>(vectorField);
        auto fastMarching = itk::FastMarchingUpwindGradientImageFilter<IsotropicMetricType, ScalarImageType>::New();
        fastMarching->SetInput(metric);
        fastMarching->SetTrialPoints(seeds);
        fastMarching->SetGenerateGradientImage(true);
        
        fastMarching->Update();
        RescaleAndExport(fastMarching->GetOutput(), prefix+"Distance.bmp");
        
        auto gradient = fastMarching->GetGradientImage();
        IsoFM_Geodesic<ComponentType,Dimension>(geodesic, gradient);
        exporter.Open(prefix+"Geodesic.txt");
        for_each(geodesic.begin(), geodesic.end(), exporter);
    }
    
    
    {// Riemannian test.
        const std::string prefix = testPrefix+"Riemannian/";
        
        auto metric = converter.MakeMetric<RiemannianMetricType>(vectorField);
        auto fastMarching = itk::AnisotropicFastMarchingImageFilter<RiemannianMetricType::PixelType>::New();
        fastMarching->SetInput(metric);
        fastMarching->SetTrialPoints(seeds);
        fastMarching->SetGenerateUpwindGradient(true);
        
        fastMarching->Update();
        RescaleAndExport(fastMarching->GetOutput(), prefix+"Distance.bmp");

        geodesic.resize(1);
        fastMarching->Geodesic(geodesic);
        
        exporter.Open(prefix+"Geodesic.txt");
        for_each(geodesic.begin(), geodesic.end(), exporter);
    }

    {// Finsler test.
        const std::string prefix = testPrefix+"Finsler/";
        
        auto metric = converter.MakeMetric<FinslerMetricType>(vectorField);
        auto fastMarching = itk::AnisotropicFastMarchingImageFilter<FinslerMetricType::PixelType>::New();
        fastMarching->SetInput(metric);
        fastMarching->SetTrialPoints(seeds);
        fastMarching->SetGenerateUpwindGradient(true);
        
        fastMarching->Update();
        RescaleAndExport(fastMarching->GetOutput(), prefix+"Distance.bmp");
        
        geodesic.resize(1);
        fastMarching->Geodesic(geodesic);
        
        exporter.Open(prefix+"Geodesic.txt");
        for_each(geodesic.begin(), geodesic.end(), exporter);
        
#pragma message("To do: compute the path from tip to seed, which should be different")
    }
        
    return EXIT_SUCCESS;
}

int TubularBand_3D_Test(){
    const std::string testPrefix = "TubularBand_3D_Test/";
    
    const unsigned int Dimension = 3;
    typedef float ComponentType;
    typedef itk::TubularBandSource<ComponentType,Dimension> TubularBandSource;
    typedef TubularBandSource::PointType PointType;
    typedef TubularBandSource::SpacingType SpacingType;
    typedef itk::Image<ComponentType, Dimension> ScalarImageType;
    
    // Define the curve
    
    struct PathType {
        int sampling = 100;
        SpacingType spacing = SpacingType(0.1);
        ComponentType radius = 0.7, margin=0.8;
        PointType operator()(ComponentType t){
            t=(2*t-1)*M_PI;
            PointType p;
            p[0]=t; p[1]=cos(t); p[2]=sin(t);
            return p;
        }
        PointType geodesicTip, geodesicSeed;
        PathType(){geodesicTip = this->operator()(0); geodesicSeed = this->operator()(1);}
    } path;
    
    // TubularBandSource
    
    auto tubularBandSource = TubularBandSource::New();
    tubularBandSource->radius = path.radius;
    tubularBandSource->margin = path.margin;
    tubularBandSource->spacing = path.spacing;
    tubularBandSource->paths.resize(2);
    for(int i=0; i<path.sampling; ++i)
        tubularBandSource->paths[0].push_back(path((i+0.)/path.sampling));
    for(auto x : tubularBandSource->paths[0]){
        PointType oppositeX;
        oppositeX[0]=x[0];
        for(int i=1; i<Dimension; ++i) oppositeX[i]=-x[i];
        tubularBandSource->paths[1].push_back(oppositeX);
    }
    
    tubularBandSource->Update();
    auto vectorField = tubularBandSource->GetOutput();
    auto region = vectorField->GetBufferedRegion();
    
    std::cout << region << endl;
    
    //Which export format ?
    std::ofstream mathematicaImageOut;
    {
        const std::string filename = testPrefix+"Euclidean_distance_to_path.txt";
        mathematicaImageOut.open(filename.c_str());
        itk::ImageRegionConstIterator<ScalarImageType> distanceIt(tubularBandSource->GetDistance(),region);
        ExportTensorToTXT(mathematicaImageOut, distanceIt, region.GetSize());
        mathematicaImageOut.close();
        //RescaleAndExport(tubularBandSource->GetDistance(), "Euclidean_distance_to_path.bmp");
    }
    
    auto vectorFieldSupport = ScalarImageType::New();
    vectorFieldSupport->SetRegions(region);
    vectorFieldSupport->Allocate();
    itk::ImageRegionConstIterator<TubularBandSource::OutputImageType> vectorFieldIt(vectorField,region);
    itk::ImageRegionIterator<ScalarImageType> vectorFieldSupportIt(vectorFieldSupport,region);
    for(vectorFieldIt.GoToBegin(), vectorFieldSupportIt.GoToBegin(); !vectorFieldIt.IsAtEnd(); ++vectorFieldIt, ++vectorFieldSupportIt)
        vectorFieldSupportIt.Value() = vectorFieldIt.Value().GetNorm();
    {
        const std::string filename = testPrefix+"Vector_Field_Support.txt";
        mathematicaImageOut.open(filename.c_str());
        ExportTensorToTXT(mathematicaImageOut, vectorFieldSupportIt, region.GetSize());
        mathematicaImageOut.close();
//        RescaleAndExport(vectorFieldSupport.GetPointer(), "Vector_Field_Support.bmp");
    }
    // Create metrics
    
    VectorToMetricCaster<ComponentType> converter;
    converter.SetKappa(10);
    
    typedef ScalarImageType IsotropicMetricType;
    typedef itk::Image<itk::Riemannian3DNorm<ComponentType>, Dimension> RiemannianMetricType;
    
    typedef itk::FastMarchingImageFilter<IsotropicMetricType,ScalarImageType>::NodeContainer NodeContainer;
    NodeContainer::Element node;
    vectorField->TransformPhysicalPointToIndex(path.geodesicSeed, node.GetIndex());
    node.SetValue(0);
    
    auto seeds = NodeContainer::New();
    seeds->Initialize();
    seeds->InsertElement(0,node);
    
    std::vector<PointType> geodesic;
    geodesic.push_back(path.geodesicTip);
    
    std::ofstream data_out;
    MathematicaExporterType exporter(data_out);
    
    {// Isotropic test
        const std::string prefix = testPrefix+"Isotropic/";
        
        auto metric = converter.MakeMetric<IsotropicMetricType>(vectorField);
        auto fastMarching = itk::FastMarchingUpwindGradientImageFilter<IsotropicMetricType, ScalarImageType>::New();
        fastMarching->SetInput(metric);
        fastMarching->SetTrialPoints(seeds);
        fastMarching->SetGenerateGradientImage(true);
        
        fastMarching->Update();
//        RescaleAndExport(fastMarching->GetOutput(), prefix+"Distance.bmp");
        
        auto gradient = fastMarching->GetGradientImage();
        IsoFM_Geodesic<ComponentType,Dimension>(geodesic, gradient);
        exporter.Open(prefix+"Geodesic.txt");
        for_each(geodesic.begin(), geodesic.end(), exporter);
    }
    
    
    {// Riemannian test.
        const std::string prefix = testPrefix+"Riemannian/";
        
        auto metric = converter.MakeMetric<RiemannianMetricType>(vectorField);
        auto fastMarching = itk::AnisotropicFastMarchingImageFilter<RiemannianMetricType::PixelType>::New();
        fastMarching->SetInput(metric);
        fastMarching->SetTrialPoints(seeds);
        fastMarching->SetGenerateUpwindGradient(true);
        
        fastMarching->Update();
//        RescaleAndExport(fastMarching->GetOutput(), prefix+"Distance.bmp");
        
        geodesic.resize(1);
        fastMarching->Geodesic(geodesic);
        
        exporter.Open(prefix+"Geodesic.txt");
        for_each(geodesic.begin(), geodesic.end(), exporter);
    }
    
    return EXIT_SUCCESS;
}


#endif
