//
//  Tutorial3D.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 08/07/13.
//
//

#ifndef ITKFM_Tutorial3D_h
#define ITKFM_Tutorial3D_h

#include <fstream>
#include "Riemannian3DNorm.h"
#include "ExportToMathematica.h"

int TutorialRiemannian3D(){
    typedef itk::Riemannian3DNorm<double> NormType;
    const unsigned int Dimension=NormType::Dimension;
    
    typedef itk::Image<NormType,Dimension> MetricType;
    auto metric = MetricType::New();
    
    
    const unsigned int n=50;
    {
        itk::Index<Dimension> index; index.Fill(0);
        
        itk::Size<Dimension> size;
        size[0]=n; size[1]=n; size[2]=2.5*n;
        itk::ImageRegion<Dimension> region(index,size);
        
        metric->SetRegions(region);
        
        MetricType::SpacingType spacing;
        spacing[0]=2.2/n; spacing[1]=spacing[0]; spacing[2]=double(3)/size[2];
        metric->SetSpacing(spacing);
        
        MetricType::PointType origin;
        origin[0] = -1.1; origin[1]=origin[0]; origin[2]=0;
        metric->SetOrigin(origin);
    }

    metric->Allocate();
    
    {
        const double k=0.2; // rotation parameter
        const double delta=0.2; // tunnel width
        const double lambda=0.02*0.02;
        
        itk::ImageRegionIteratorWithIndex<MetricType> metricIt(metric, metric->GetBufferedRegion());
        
        for(metricIt.GoToBegin(); !metricIt.IsAtEnd(); ++metricIt){
            metricIt.Value().SetIdentity();
            
            MetricType::PointType p;
            metric->TransformIndexToPhysicalPoint(metricIt.GetIndex(), p);
            
            const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
            if( fabs(r-1) < delta
               && fabs( cos(p[2]/k) - p[0] ) < delta
               && fabs( sin(p[2]/k) - p[1] ) < delta
               ) {
                
                itk::Vector<double,Dimension> v; // tangent vector to the test curve
                v[0]=-p[1]/r; v[1]=p[0]/r; v[2]=k;
                v.Normalize();
                
                for(int i=0; i<3; ++i)
                    for(int j=i; j<3; ++j)
                        metricIt.Value()(i, j) -= (1-lambda)*v[i]*v[j];
            } // if on spiral
        }// for metric It
        
        std::ofstream of;
        of.open("Metric_Riemannian3DNorm.txt");
        ExportTensorToTXT(of,metricIt,metric->GetBufferedRegion().GetSize());
        
    }
    
    typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
    auto FM = FMType::New();
    FM->SetInput(metric);
    
    {
        FMType::NodeType node;
        node.GetIndex()[0]=n/2;
        node.GetIndex()[1]=node.GetIndex()[0];
        node.GetIndex()[2]=0;
        node.GetValue()=0;
        
        auto nodeContainer=FMType::NodeContainer::New();
        nodeContainer->Initialize();
        nodeContainer->InsertElement(0, node);
        
        FM->SetTrialPoints(nodeContainer);

    }
    clock_t t = -clock();
    FM->Update();
    t+=clock();
    cout << "FM run took " << t/double(CLOCKS_PER_SEC) << " s\n";
    {
        std::ofstream of;
        of.open("Distance_Riemannian3DNorm.txt");
        auto it = itk::ImageRegionIterator<FMType::LevelSetImageType>(FM->GetOutput(), FM->GetOutput()->GetBufferedRegion() );
        ExportTensorToTXT(of,it,FM->GetOutput()->GetBufferedRegion().GetSize());
    }
    
    return EXIT_SUCCESS;
}

#endif
