//
//  Segmentation_Test.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 13/09/13.
//
//

#ifndef ITKFM_Segmentation_Test_h
#define ITKFM_Segmentation_Test_h

#include "IsoFM_Geodesic.h"

namespace AnisotropicFastMarching_Segmentation_Test {

    const std::string testPrefix = "Segmentation_Test/";
    const std::string testImageFilename = "TestImage.hdf5";

    typedef float ScalarType;
    const unsigned int Dimension = 2;
    typedef itk::Image<ScalarType,Dimension> ScalarImageType;
    typedef ScalarImageType::SpacingType SpacingType;
    typedef ScalarType ComponentType;
    typedef itk::Point<ComponentType,Dimension> PointType;
    
    struct PathType {
        ScalarType sampling = 100, width=0.4;
        const SpacingType spacing = SpacingType(1/8.);
        PointType operator()(ScalarType t){
            t*=2*M_PI;
            PointType p;
            p[0]=2*cos(t); p[1]=sin(t)+0.6*sin(3*t);
            return p;
        }
        PointType geodesicTip, geodesicSeed, origin;
        PathType(){geodesicTip = this->operator()(1./12); geodesicSeed = this->operator()(7./12); origin.Fill(0);}
    } path;

    inline ScalarType TangentSpeed(ScalarType n){
        const ScalarType alpha = 2.;
        return ScalarType(1.+pow(n,alpha)) ;
    }
    
    inline ScalarType OrthogonalSpeed(ScalarType tangentSpeed){
        const ScalarType maxAnisotropy = 10;
        return std::max(ScalarType(1), tangentSpeed/maxAnisotropy);
    }
    
    void RunTest(){
        // Get the image gradient
        auto reader = itk::ImageFileReader<ScalarImageType>::New();
        reader->SetFileName(testPrefix+testImageFilename);
        reader->UpdateOutputInformation();
        {
            const auto io = reader->GetImageIO();
            assert(io->GetNumberOfDimensions() == Dimension);
            assert(io->GetNumberOfComponents() == 1);
            assert(io->GetComponentType() == itk::ImageIOBase::FLOAT);
        }
        const ScalarImageType * physicalImage = reader->GetOutput();
        const auto region = physicalImage->GetLargestPossibleRegion();
        
        typedef itk::GradientImageFilter<ScalarImageType,float,float> GradientImageFilterType;
        typedef GradientImageFilterType::OutputImageType GradientImageType;
        typedef GradientImageType::PixelType CovariantVectorType;
        
        auto gradientFilter = GradientImageFilterType::New();
        gradientFilter->SetInput(reader->GetOutput());
        
        // Prepare for geodesic extraction
        
        typedef itk::FastMarchingImageFilter<ScalarImageType,ScalarImageType>::NodeContainer NodeContainer;
        NodeContainer::Element node;
        node.SetValue(0);
        physicalImage->TransformPhysicalPointToIndex(path.geodesicSeed, node.GetIndex());
        
        auto seeds = NodeContainer::New();
        seeds->Initialize();
        seeds->InsertElement(0,node);
        
        typedef itk::ContinuousIndex<ComponentType,Dimension> ContinuousIndexType;
        std::vector<ContinuousIndexType> geodesic;
        ContinuousIndexType geodesicTipCIndex;
        physicalImage->TransformPhysicalPointToContinuousIndex(path.geodesicTip, geodesicTipCIndex);
        geodesic.push_back(geodesicTipCIndex);
        assert(region.IsInside(geodesic.back()));
        
        std::ofstream data_out;
        itk::MathematicaExporterType exporter(data_out);
        
        {// Isotropic metric            
            struct FunctorType {
                ScalarType operator()(CovariantVectorType g) const {
                    return TangentSpeed(g.GetNorm());
                }
            };
            
            auto imageFunctor = itk::UnaryFunctorImageFilter<GradientImageType, ScalarImageType, FunctorType>::New();
            imageFunctor->SetInput(gradientFilter->GetOutput());
            
            auto fm = itk::FastMarchingUpwindGradientImageFilter<ScalarImageType>::New();
            fm->SetInput(imageFunctor->GetOutput());
            fm->SetGenerateGradientImage(true);
            fm->SetTrialPoints(seeds);
            
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(fm->GetOutput());
            writer->SetFileName(testPrefix+"IsotropicDistance.hdf5");
            
            writer->Update();
            
            //RescaleAndExport<ScalarImageType,ScalarImageType>(fm->GetOutput(), testPrefix+"RiemannianDistance.hdf5");
            
            geodesic.resize(1);
            IsoFM_Geodesic<ComponentType,Dimension>(geodesic, fm->GetGradientImage());
            exporter.open(testPrefix+"IsotropicGeodesic.txt");
            std::for_each(geodesic.begin(),geodesic.end(),exporter);
            exporter.close();
            
        }
        
        {// Riemannian metric
            typedef itk::Riemannian2DNorm<float> NormType;
            
            struct FunctorType {
                NormType operator()(CovariantVectorType g) const {
                    NormType m;
                    const ScalarType tangentSpeed = TangentSpeed(g.GetNorm());
                    const ScalarType orthogonalSpeed = OrthogonalSpeed(tangentSpeed);
                    
                    const ScalarType tangentSpeedInv2 = 1/(tangentSpeed*tangentSpeed);
                    const ScalarType orthogonalSpeedInv2 = 1/(orthogonalSpeed*orthogonalSpeed);
                    if(g.GetNorm()!=0) g.Normalize();
                    
                    for(int i=0; i<Dimension; ++i)
                        for(int j=0; j<=i; ++j)
                            m(i,j) = tangentSpeedInv2*(i==j) + (orthogonalSpeedInv2-tangentSpeedInv2)*g[i]*g[j];
                    assert(m.IsDefinite());
                    return m;
                }
            };
            
            typedef itk::Image<NormType,Dimension> MetricType;
            auto imageFunctor = itk::UnaryFunctorImageFilter<GradientImageType, MetricType, FunctorType>::New();
            imageFunctor->SetInput(gradientFilter->GetOutput());
            
            auto fm = itk::AnisotropicFastMarchingImageFilter<NormType>::New();
            fm->SetInput(imageFunctor->GetOutput());
            fm->SetGenerateUpwindGradient(true);
            fm->SetTrialPoints(seeds);
            fm->Update();
            
            RescaleAndExport<ScalarImageType,ScalarImageType>(fm->GetOutput(), testPrefix+"RiemannianDistance.hdf5");
            
            geodesic.resize(1);
            fm->Geodesic(geodesic);
            exporter.open(testPrefix+"RiemannianGeodesic.txt");
            std::for_each(geodesic.begin(), geodesic.end(),exporter);
            exporter.close();            
        }

        for(int reversed=0; reversed<2; ++reversed){// Finsler metric
            if(reversed){
                physicalImage->TransformPhysicalPointToContinuousIndex(path.geodesicSeed, geodesic[0]);
                physicalImage->TransformPhysicalPointToIndex(path.geodesicTip, seeds->operator[](0).GetIndex());
            }
            const std::string metricName = reversed ? "FinslerRev" : "Finsler";
            
            typedef itk::Finsler2DNorm<float> NormType;
            
            struct FunctorType {
                NormType operator()(CovariantVectorType g){
                    NormType m;
                    const ScalarType tangentSpeed = TangentSpeed(g.GetNorm());
                    const ScalarType orthogonalSpeed = OrthogonalSpeed(tangentSpeed);

                    if(g.GetNorm()!=0) g.Normalize();
                    std::swap(g[0],g[1]); // rotate g
                    g[0]*=-1;
                    
                    m.SetIdentity();
                    m.GetM() /= orthogonalSpeed*orthogonalSpeed;
                    for(int i=0; i<Dimension; ++i)
                        m.GetOmega()[i] = (1./orthogonalSpeed-1./tangentSpeed)*g[i];
                    return m;
                }
            };

            typedef itk::Image<NormType,Dimension> MetricType;
            auto imageFunctor = itk::UnaryFunctorImageFilter<GradientImageType, MetricType, FunctorType>::New();
            imageFunctor->SetInput(gradientFilter->GetOutput());
                        
            auto fm = itk::AnisotropicFastMarchingImageFilter<NormType>::New();
            fm->SetInput(imageFunctor->GetOutput());
            fm->SetGenerateUpwindGradient(true);
            fm->SetTrialPoints(seeds);
            fm->Update();
            
            RescaleAndExport<ScalarImageType,ScalarImageType>(fm->GetOutput(), testPrefix+metricName+"Distance.hdf5");
            
            geodesic.resize(1);
            fm->Geodesic(geodesic);
            exporter.open(testPrefix+metricName+"Geodesic.txt");
            std::for_each(geodesic.begin(), geodesic.end(),exporter);
            exporter.close();
        }
    }

}


#endif
