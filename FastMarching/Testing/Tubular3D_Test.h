//
//  Tubular3D_Test.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 13/09/13.
//
//

#ifndef ITKFM_Tubular3D_Test_h
#define ITKFM_Tubular3D_Test_h

#include <itkRecursiveGaussianImageFilter.h>
#include <itkSubtractImageFilter.h>
#include "itkBinaryFunctorImageFilter.h"
#include <itkHessianRecursiveGaussianImageFilter.h>
#include "ExtendedNorm.h"


void Tubular3D_Test(){
    
    // Three different models for motion planning
    const std::string testPrefix = "Tubular_Test/";
#pragma message("Extension is not fine")
    const std::string imageFilename = "Segmentation_Test_Image.bmp";
    
    const unsigned int Dimension = 3;
    typedef float ComponentType;
    //    typedef itk::Point<ComponentType,Dimension> PointType;
    typedef itk::ContinuousIndex<ComponentType,Dimension> ContinuousIndexType;
    typedef itk::Image<ComponentType, Dimension> ScalarImageType;
    typedef itk::Image<unsigned char, Dimension> SourceImageType;
    
    // Get the image gradient
    
    auto reader = itk::ImageFileReader<SourceImageType>::New();
    reader->SetFileName(testPrefix+imageFilename);
	reader->UpdateOutputInformation();
    {
        const auto io = reader->GetImageIO();
        assert(io->GetNumberOfDimensions() == Dimension);
        assert(io->GetNumberOfDimensions() == 1);
        assert(io->GetComponentType() == itk::ImageIOBase::UCHAR);
    }
    
    // Alternatively,
    // itk::HessianRecursiveGaussianImageFilter< TInputImage, TOutputImage >,
    // itk::Hessian3DToVesselnessMeasureImageFilter< TPixel > Class Template,
    // #include <itkMultiScaleHessianBasedMeasureImageFilter.h>
    // Anisotropic fast marching is not an extension of a specific methods,
    // but a framework to propose new ideas.
    
    /*
    typedef itk::RecursiveGaussianImageFilter<SourceImageType,ScalarImageType> GaussianFilterType;
    auto gaussianFilter = GaussianFilterType::New();
    auto gaussianHalfFilter = GaussianFilterType::New();
    
    typedef itk::Functor::Sub2<ComponentType> SubstractFunctorType;
    typedef itk::BinaryFunctorImageFilter<ScalarImageType, ScalarImageType, ScalarImageType, SubstractFunctorType> SubstractFilterType;
    auto substractFilter = SubstractFilterType::New();
    
    substractFilter->SetInput1(gaussianFilter->GetOutput());
    substractFilter->SetInput2(gaussianHalfFilter->GetOutput());
    */
    
    
    typedef itk::HessianRecursiveGaussianImageFilter<SourceImageType, ScalarImageType> HessianFilterType;
    auto hessianFilter = HessianFilterType::New();
    hessianFilter->SetInput(reader->GetOutput());
    typedef itk::Riemannian3DNorm<float> TensorType;
    typedef itk::Image<TensorType,Dimension> TensorImageType;
    
    struct HessianToMetricType {
        bool operator!=(const HessianToMetricType &) const {return false;}
        bool operator==(const HessianToMetricType & other) const {return !( *this != other );}
        
        ComponentType alpha, beta, gamma, anisotropy_ratio;
        
        
        TensorType operator()(const TensorType & hessian){
                        
            struct EigenValuesComp {
                TensorType::EigenValuesArrayType eigenVals;
                bool operator()(int i, int j){return fabs(eigenVals[i]) <= fabs(eigenVals[j]);}
            } comp;
            
            auto & eigenVals = comp.eigenVals;
            TensorType::EigenVectorsMatrixType eigenVecs;
            hessian.ComputeEigenAnalysis(eigenVals, eigenVecs);

            int order[Dimension] = {0,1,2};
            std::sort(order,order+Dimension,comp);
            
            ComponentType vesselness;
            const ComponentType lambda1 = eigenVals[order[0]], lambda2 = eigenVals[order[1]], lambda3 = eigenVals[order[2]];
            
            if(lambda2>=0 || lambda3>=0) vesselness=0;
            else {
                const ComponentType rA = fabs(lambda2/lambda3), rB=fabs(lambda1/(lambda2*lambda3)), s=sqrt(lambda1*lambda1+lambda2*lambda2+lambda3*lambda3);
                vesselness = (1.-exp(-0.5*vnl_math_sqr(rA/alpha)))*exp(-0.5*vnl_math_sqr(rB/beta))*(1-exp(-0.5*vnl_math_sqr(s/gamma)));
            }
            const ComponentType inverseSquaredSpeed = 1-vesselness;
            
            TensorType output;
            for(int i=0; i<Dimension; ++i){
                for(int j=0; j<i; ++j)
                    output(i,j)=0;
                if(i==order[0])
                    output(i,i) = inverseSquaredSpeed;
                else
                    output(i,i) = std::min(ComponentType(1),vnl_math_sqr(anisotropy_ratio)*inverseSquaredSpeed);
            }
            return static_cast<TensorType>( output.Rotate(eigenVecs.GetTranspose()) );
        }
    };
    
    auto hessianToMetricFilter = itk::UnaryFunctorImageFilter<TensorImageType, TensorImageType, HessianToMetricType>::New();
    auto & hessianToMetricFunctor = hessianToMetricFilter->GetFunctor();
    // Set parameters of functor
    
    std::ofstream data_out;
    itk::MathematicaExporterType exporter(data_out);

    {// Single scale
        // Set scale of hessian
        
        typedef TensorType NormType;
        typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
        auto fm = FMType::New();
        fm->SetInput(hessianToMetricFilter->GetOutput());
        
        
        FMType::NodeType node;
        // set node index
        node.SetValue(0);
        
        
        auto seeds = FMType::NodeContainer::New();
        seeds->Initialize();
        seeds->InsertElement(0, node);
        fm->SetTrialPoints(seeds);
        
        fm->Update();
        // Export distance
        
        typedef itk::ContinuousIndex<ComponentType,Dimension> ContinuousIndexType;
        std::vector<ContinuousIndexType> geodesic;
        // set geodesic start
        
        fm->Geodesic(geodesic);
        exporter.open(testPrefix+"SingleScaleGeodesic.txt");
        for_each(geodesic.begin(), geodesic.end(), exporter);
    }
    
    { // Multi-scale
        const size_t numberOfScales = 5;
        
        typedef itk::ExtendedNorm<TensorType> NormType;
        typedef itk::Image<NormType, NormType::Dimension> MetricType;
        auto metric = MetricType::New();
        
        // Set the metric
        const auto baseRegion = reader->GetOutput()->GetBufferedRegion();
        
        MetricType::IndexType index;
        MetricType::SizeType size;
        for(int i=0; i<Dimension; ++i){
            index[i] = baseRegion.GetIndex()[i];
            size[i] = baseRegion.GetSize()[i];
        }
        index[Dimension]=0;
        size[Dimension] =numberOfScales;
        
        MetricType::RegionType region(index,size);
        metric->SetRegions(region);
        metric->Allocate();
        
        // Set the metric
        for(size_t iScale = 0; iScale<numberOfScales; ++iScale){
            // Set scale
            
            hessianToMetricFilter->Update();
            
        }
    
        
        typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
        auto fm = FMType::New();

    }
    
}

#endif
