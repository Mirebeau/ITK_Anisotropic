//
//  Topographic_Test.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 04/09/13.
//
//

#ifndef ITKFM_Topographic_Test_h
#define ITKFM_Topographic_Test_h

#include "CurveNeighborhoodSource.h"
#include "itkAnisotropicFastMarchingImageFilter.h"
#include "Finsler2DNorm.h"
//#include "TubularBand_Test.h"
#include "itkImageFileReader.h"
#include "itkGradientImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "ExportToMathematica.h"
#include "RescaleAndExport_Function.h"

namespace AnisotropicFastMarching_Topographic_Test {

    const std::string testPrefix = "Topographic_Test/";
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
            t=(2*t-1)*15;
            PointType p;
            p[0]=t; p[1]= (t>0?1:-1)*pow(fabs(t),2./3.)*cos(t);
            return p;
        }
        PointType geodesicTip, geodesicSeed, origin;
        PathType(){geodesicTip = this->operator()(0); geodesicSeed = this->operator()(1); origin.Fill(0);}
    } path;
    
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
        
        {// What is the max norm ? (Hence anisotropy of Riemannian metric)
            auto magnitudeFilter = itk::UnaryFunctorImageFilter<GradientImageType, ScalarImageType, itk::Functor::VectorMagnitude<CovariantVectorType, ScalarType> >::New();
            magnitudeFilter->SetInput(gradientFilter->GetOutput());
            magnitudeFilter->Update();
            
            auto calculator = itk::MinimumMaximumImageCalculator<ScalarImageType>::New();
            calculator->SetImage(magnitudeFilter->GetOutput());
            calculator->SetRegion(region);
            calculator->ComputeMaximum();
            
            const ScalarType maxNorm = calculator->GetMaximum();
            const ScalarType maxKappa = sqrt(1+maxNorm*maxNorm);
            std::cout << "Max magnitude : " << maxNorm << ", max anisotropy of Riemannian : " << maxKappa << std::endl;
        }
        
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
        
        {// Define the Riemannian metric
            typedef itk::Riemannian2DNorm<float> NormType;
            
            struct FunctorType {
                NormType operator()(CovariantVectorType g){
                    NormType m;
                    for(int i=0; i<(int)Dimension; ++i)
                        for(int j=0; j<=i; ++j)
                            m(i,j) = (i==j)+0.3*g[i]*g[j];
                    return m;
                }
            };
            
            typedef itk::Image<NormType,Dimension> MetricType;
            auto functorFilter = itk::UnaryFunctorImageFilter<GradientImageType, MetricType, FunctorType>::New();
            functorFilter->SetInput(gradientFilter->GetOutput());
            
            typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
            auto fm = FMType::New();
            fm->SetInput(functorFilter->GetOutput());
            fm->SetGenerateUpwindGradient(true);
            fm->SetTrialPoints(seeds);

            
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetFileName(testPrefix+"RiemannianDistance.hdf5");
            writer->SetInput(fm->GetOutput());
            writer->Update();
            
            geodesic.resize(1);
            fm->Geodesic(geodesic);
            exporter.open(testPrefix+"RiemannianGeodesic.txt");
            std::for_each(geodesic.begin(), geodesic.end(),exporter);
            exporter.close();
            
        }
        
        {// Define the Finsler metric            
            typedef itk::Finsler2DNorm<float> NormType;
            
            struct FunctorType {
                ScalarType alpha;
                NormType operator()(const CovariantVectorType & g) const {
                    itk::Vector<ScalarType,Dimension> v;
                    v[0]=g[0]; v[1]=g[1];
                    itk::Riemannian2DNorm<ScalarType> m;
                    m.SetIdentity();
                    
                    ScalarType s = 1./(ScalarType(1)+std::pow(ScalarType(v.GetNorm()),alpha) );
                    s=std::min(s,ScalarType(50));
                    if(v.GetNorm()!=0) v.Normalize();
                    auto result = NormType::TranslatedEllipse(m,(ScalarType(1)-s)*v);
                    assert(result.IsDefinite());
                    return result;
                }
            };
            
            typedef itk::Image<NormType,Dimension> MetricType;
            auto functorFilter = itk::UnaryFunctorImageFilter<GradientImageType, MetricType, FunctorType>::New();
            functorFilter->GetFunctor().alpha=2;
            functorFilter->SetInput(gradientFilter->GetOutput());
            
            auto fm = itk::AnisotropicFastMarchingImageFilter<NormType>::New();
            fm->SetInput(functorFilter->GetOutput());
            fm->SetGenerateUpwindGradient(true);
            fm->SetTrialPoints(seeds);
            
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetFileName(testPrefix+"FinslerDistance.hdf5");
            writer->SetInput(fm->GetOutput());
            writer->Update();
                        
            geodesic.resize(1);
            fm->Geodesic(geodesic);
            exporter.open(testPrefix+"FinslerGeodesic.txt");
            std::for_each(geodesic.begin(), geodesic.end(),exporter);
            exporter.close();
        }
    }
    
} // end of namespace

#endif
