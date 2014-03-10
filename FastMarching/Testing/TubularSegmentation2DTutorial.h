//
//  SegmentationExample.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 10/06/13.
//
//

#ifndef ITKFM_TubularSegmentation2DTutorial_h
#define ITKFM_TubularSegmentation2DTutorial_h

#include "itkImageFileReader.h"
#include "itkOrientedFluxMatrixImageFilter.h"
#include "itkPolyLineParametricTubularPath.h"
#include "itkImageRegionIterator.h"

int TubularSegmentation2DTutorial(int argc, char* argv[]){
    int argn=1;
    
    const unsigned int Dimension = 2;
    typedef unsigned char ExternalPixelType;
    typedef itk::Image<ExternalPixelType,Dimension> ExternalImageType;
    
    // Setup the reader
    typedef itk::ImageFileReader<ExternalImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    {
        if(argc<=argn) itkGenericExceptionMacro("Not enough arguments");
        reader->SetFileName(argv[argn++]);
        reader->UpdateOutputInformation();
        
        if(reader->GetImageIO()->GetNumberOfDimensions() != 2) itkGenericExceptionMacro("Image dimension should be 2");
        if(reader->GetImageIO()->GetPixelType() != itk::ImageIOBase::UCHAR) itkGenericExceptionMacro("Image pixels should be unsigned chars");
        
    //    reader->Update();
    }
    
    // Construct the tubularity measure
    typedef itk::OrientedFluxMatrixImageFilter<ExternalImageType> OrientedFluxFilterType;
    OrientedFluxFilterType::Pointer orientedFluxFilter = OrientedFluxFilterType::New();
    {
        double radius;
        if(argc<=argn) radius = atof(argv[argn++]);
        else radius = 5;
        cout << "Radius set to " << radius << endl;
        orientedFluxFilter->SetRadius(radius);
        
//        orientedFluxFilter->
        
        
        orientedFluxFilter->SetInput(reader->GetOutput());
        orientedFluxFilter->Update();
    }
    
    // Get the tensors
    typedef OrientedFluxFilterType::OutputImageType TensorImageType;
    TensorImageType * tensors = orientedFluxFilter->GetOutput();
    
    // Get the region
    typedef TensorImageType::RegionType RegionType;
    RegionType region = tensors->GetLargestPossibleRegion();
    
    
    typedef TensorImageType::PixelType::ComponentType ValueType;
    
    // Construct the metric
    typedef itk::Riemannian2DNorm<ValueType> NormType;
    typedef itk::Image<NormType,Dimension> MetricType;
    MetricType::Pointer metric = MetricType::New();
    {
        metric->SetRegions(region);
        metric->SetSpacing(tensors->GetSpacing());
        metric->Allocate();
    }
    
    // Set norms
    {
        itk::ImageRegionIterator<MetricType> metricIt(metric,region);
        itk::ImageRegionConstIterator<TensorImageType> tensorIt(tensors,region);
        
        for(metricIt.GoToBegin(), tensorIt.GoToBegin();
            !metricIt.IsAtEnd();
            ++metricIt,++tensorIt)
        {
            // Eigen analysis of the tensors.
            
            // Riemannian Norm with same eigenvectors.
        }
    }
    
    // Smooth out ? Anisotropically ?
    
    // Get starting point. Compute fast marching.
    
    
    // Get regions of large divergence. Get connected component of the origin.
    
    // Extract tree / paths to the origin.
    
    
    
}

#endif
