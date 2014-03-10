//
//  Tutorial.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 10/06/13.
//
//

#ifndef ITKFM_Tutorial_h
#define ITKFM_Tutorial_h

#include "itkAnisotropicFastMarchingImageFilter.h"
#include "Finsler2DNorm.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "math.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkPathToImageFilter.h"

#ifdef Export_SWC_Path
#pragma message("This include should be in itkPolyLineParametricTubularPath")
#include <fstream> 
#include "itkPolyLineParametricTubularPath.h"
#endif

int Finsler2DTutorial(){
    const int NumberOfSidePixels = 200;
    const float anisotropyCoefficient = 0.9; // should be 0 <= ... < 1.
    
    typedef itk::Finsler2DNorm<float> NormType;
    typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
    
    typedef FMType::NormImageType MetricType;
    
    // set region
    typedef MetricType::RegionType RegionType;
    RegionType region;
    {
        RegionType::SizeType size;
        size.Fill(NumberOfSidePixels);
        RegionType::IndexType index;
        index.Fill(-NumberOfSidePixels/2);
        region = RegionType(index,size);
        
        cout << region << endl;
    }
    
    typedef FMType::PointType   PointType;
    typedef FMType::VectorType  VectorType;
    typedef FMType::ValueType   ValueType;
    
    //Set the norms
    MetricType::Pointer metric = MetricType::New();
    {
        metric->SetRegions(region);
        metric->Allocate();
        metric->SetSpacing(1./NumberOfSidePixels);

        itk::ImageRegionIteratorWithIndex<MetricType> it(metric,region);
        for(it.GoToBegin(); !it.IsAtEnd(); ++it){
            PointType p;
            metric->TransformIndexToPhysicalPoint(it.GetIndex(), p);
            const ValueType scalar = sin(4*M_PI*p[0])*sin(4*M_PI*p[1]);

            VectorType v = p.GetVectorFromOrigin();
            if(v.GetSquaredNorm()>0) v.Normalize();
            
            it.Value().GetM().SetIdentity();
            it.Value().GetOmega() = anisotropyCoefficient*scalar*v;
        }
    }
    
    // Set the seeds
    typedef FMType::NodeContainer NodeContainer;
    NodeContainer::Pointer seeds = NodeContainer::New();
    {
        FMType::NodeType node;
        node.GetIndex().Fill(0);
        node.SetValue(0.);
        
        seeds->Initialize();
        seeds->InsertElement(0, node);
    }
    
    // Setup fast marching
    FMType::Pointer fm = FMType::New();
    {
        fm->SetInput(metric);
        fm->SetTrialPoints(seeds);
        fm->SetGenerateUpwindGradient(true);
        fm->SetOutputRegion(region);
    }
        
    // Set up the rescaling filter.
    typedef FMType::LevelSetImageType LevelSetImageType;
    typedef itk::Image<unsigned char, NormType::Dimension> OutputImageType;
    typedef itk::RescaleIntensityImageFilter<LevelSetImageType, OutputImageType >   CastFilterType;
    CastFilterType::Pointer caster = CastFilterType::New();
    {
        caster->SetInput(fm->GetOutput());
        caster->SetOutputMinimum(0);
        caster->SetOutputMaximum(255);
    }
    
    //Setting up the writer
    typedef  itk::ImageFileWriter< OutputImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    {
        writer->SetFileName("Tutorial_Finsler2DNorm.jpg");
        writer->SetInput(caster->GetOutput());
    }
    
    // Compute, and export image
    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cout << "Exception caught : " << e << std::endl;
    }

    /*
    // turning to geodesic computation
    typedef FMType::PathType PathType;
    PathType::Pointer path;
    {
        typedef FMType::PathContinuousIndexType PathContinuousIndexType;
        PathContinuousIndexType start;
        start.Fill(NumberOfSidePixels/3);
        path = fm->Geodesic(start);
    }
    
#ifdef EXPORT_SWC_Path
    const int Dimension = NormType::Dimension;
    typedef itk::PolyLineParametricTubularPath<Dimension> TubularPathType;
    TubularPathType::Pointer pathWithRadius = TubularPathType::New();

    {
        const PathType::VertexListType * vertexList = path->GetVertexList();
        
        cout << vertexList->Size() << endl;
        for(int i=0; i<vertexList->Size(); ++i) cout << vertexList->GetElement(i) << endl;
        
        pathWithRadius->Initialize();
        for(int i=0; i<vertexList->Size(); ++i) {
            itk::ContinuousIndex<double,Dimension> ci;
            for(int k=0; k<Dimension; ++k) ci[k] = vertexList->GetElement(i)[k] - fm->GetOutput()->GetLargestPossibleRegion().GetIndex()[k];
            pathWithRadius->AddVertex(ci);
        };
    }

    const std::string outputSWCFilename = "Tutorial_Finsler2DNorm.swc";
    {
        pathWithRadius->WriteSwcFile(outputSWCFilename, fm->GetOutput(), false);
    }
#endif
     */
    /*
    // Path drawing with ITK. Bof.
     
     
    typedef itk::PathToImageFilter<PathType, OutputImageType> PathToImageFilterType;
    PathToImageFilterType::Pointer pathToImage = PathToImageFilterType::New();
    {
        
        pathToImage->SetSize(region.GetSize());
        pathToImage->SetPathValue(255);
        pathToImage->SetBackgroundValue(0);
        pathToImage->SetInput(path);
        pathToImage->GetOutput()->SetRegions(region);

        writer->SetInput(pathToImage->GetOutput());
        writer->SetFileName("Tutorial_Finsler2DNorm_Path.png");
    }
     
    
    try {
        writer->Update();
    } catch (itk::ExceptionObject &e) {
        std::cout << "Exception caught : "<< e << std::endl;
    }
    */
    return EXIT_SUCCESS;
}

#endif
