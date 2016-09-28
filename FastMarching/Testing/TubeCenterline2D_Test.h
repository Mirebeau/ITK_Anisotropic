//
//  TubularCenterline2D_Test.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 19/09/13.
//
//

#ifndef ITKFM_TubeCenterline2D_Test_h
#define ITKFM_TubeCenterline2D_Test_h

#include "ExportToMathematica.h"

namespace AnisotropicFastMarching_TubeCenterline2D_Test {
    
    const std::string testPrefix = "TubeCenterline2D_Test/";
    const std::string testImageFilename = "TestImage.hdf5";
    
    typedef float ScalarType;
    const unsigned int Dimension = 2;
    typedef itk::Image<ScalarType,Dimension> ScalarImageType;
    typedef ScalarImageType::SpacingType SpacingType;
    typedef ScalarType ComponentType;
    typedef itk::Point<ComponentType,Dimension> PointType;
    
    struct PathType {
        PointType origin;
        SpacingType spacing;
        static const int nPoints = 101;
        ScalarType ConstantRadius = 0.8;
        
        PointType GetPoint(int i){
            assert(0<=i && i<nPoints);
            PointType p;
            p[0]=points[i][0];
            p[1]=points[i][1];
            return p;
        }
        ScalarType Radius(int i){
            const ScalarType s = 1.-i/100.;
            return 1/(1+2*s*s);
//            return 0.8/(1+s*s);
        }
        
        PathType(){origin.Fill(0); spacing.Fill(1./8);}
    protected:
        static const ScalarType points[nPoints][2];
    } path;

    void RunTest_ConstantRadius(){
        const ScalarType sigma = path.ConstantRadius/4;
        const ScalarType deltaMinFixed = 0.01;
        const ScalarType eVal1Threshold = 0.02/(sigma*sigma);
        const ScalarType kappa = 4;
        typedef itk::Riemannian2DNorm<ScalarType> NormType;
        
        auto reader = itk::ImageFileReader<ScalarImageType>::New();
        reader->SetFileName(testPrefix+"TestImageConstantRadius.hdf5");
        reader->UpdateOutputInformation();
        const auto * physicalImage = reader->GetOutput();
        const auto region = reader->GetOutput()->GetLargestPossibleRegion();
        
        typedef itk::HessianRecursiveGaussianImageFilter<ScalarImageType> HessianFilterType;
        typedef HessianFilterType::OutputImageType TensorImageType;
        typedef TensorImageType::PixelType TensorType;
        auto hessianFilter = HessianFilterType::New();
        hessianFilter->SetInput(reader->GetOutput());
        hessianFilter->SetSigma(sigma);
        
        struct HessianToMetricFunctor {
            ScalarType deltaMinFixed;
            ScalarType deltaMinEncountered;
            ScalarType eVal1Threshold;
            ScalarType kappa;
            bool Show;
            NormType operator()(const TensorType & H) {
                TensorType::EigenValuesArrayType eVal;
                TensorType::EigenVectorsMatrixType eVec;
                H.ComputeEigenAnalysis(eVal, eVec);
                eVec = eVec.GetTranspose();
                
                if(fabs(eVal[0]) > fabs(eVal[1])){
                    std::swap(eVal[0], eVal[1]);
                    std::swap(eVec(0,0),eVec(0,1));
                    std::swap(eVec(1,0),eVec(1,1));
                }
                
                const ScalarType tangentSpeed = exp( std::max(ScalarType(eVal[1]- 2*fabs(eVal[2]) ) / eVal1Threshold, ScalarType(0)) );
                const ScalarType normalSpeed = std::max(ScalarType(1),tangentSpeed/kappa);
                
                NormType norm;
                norm(0,0) = 1/vnl_math_sqr(tangentSpeed);
                norm(0,1) = 0;
                norm(1,1) = 1/vnl_math_sqr(normalSpeed);
                norm = NormType( norm.Rotate(eVec) );
                
                if(Show)
                    std::cout << "H : " << H << ", eVal : " << eVal << ", eVec : " << eVec << ", tangent Speed " << tangentSpeed << ", normal Speed : " << normalSpeed << ", norm : " << norm << endl;
                
                return norm;
            }
        };
        
        typedef itk::Image<NormType,Dimension> MetricType;
        auto hessianToMetricFilter = itk::UnaryFunctorImageFilter<TensorImageType, MetricType, HessianToMetricFunctor>::New();
        hessianToMetricFilter->SetInput(hessianFilter->GetOutput());
        auto & functor = hessianToMetricFilter->GetFunctor();
        functor.deltaMinEncountered = 1;
        functor.deltaMinFixed = deltaMinFixed;
        functor.eVal1Threshold = eVal1Threshold;
        functor.kappa = kappa;
        functor.Show = false;
        
        typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
        auto fm = FMType::New();
        fm->SetInput(hessianToMetricFilter->GetOutput());
        fm->SetGenerateUpwindGradient(true);
        
        auto seeds = FMType::NodeContainer::New();
        seeds->Initialize();
        FMType::NodeType node;
        node.SetValue(0);
        itk::Index<Dimension> seedIndex;
        physicalImage->TransformPhysicalPointToIndex(path.GetPoint(0), seedIndex);
        node.SetIndex(seedIndex);
        seeds->InsertElement(0,node);
        fm->SetTrialPoints(seeds);
        
        {
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(fm->GetOutput());
            writer->SetFileName(testPrefix+"ConstantRadiusDistance.hdf5");
            writer->Update();
        }
        
        std::cout << "deltaMinEncountered : " << hessianToMetricFilter->GetFunctor().deltaMinEncountered << std::endl;
        
        struct EVFunctor {
            int eVIndex;
            ScalarType operator()(const TensorType & H){
                TensorType::EigenValuesArrayType eVal;
                H.ComputeEigenValues(eVal);
                if(fabs(eVal[0]) > fabs(eVal[1])) std::swap(eVal[0],eVal[1]);
                return eVal[eVIndex];
            }
        };
        
        for(int eVIndex = 0; eVIndex<2; ++eVIndex){
            auto eVFilter = itk::UnaryFunctorImageFilter<TensorImageType, ScalarImageType, EVFunctor>::New();
            eVFilter->GetFunctor().eVIndex = eVIndex;
            eVFilter->SetInput(hessianFilter->GetOutput());
            
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(eVFilter->GetOutput());
            std::ostringstream  s;
            s << testPrefix << "eV" << eVIndex << "ConstantRadius.hdf5";
            writer->SetFileName(s.str().c_str());
            writer->Update();
        }
        
        typedef itk::ContinuousIndex<ScalarType,Dimension> ContinuousIndexType;
        ContinuousIndexType cindex;
        physicalImage->TransformPhysicalPointToContinuousIndex(path.GetPoint(path.nPoints-1), cindex);
        std::vector<ContinuousIndexType> geodesic;
        geodesic.push_back(cindex);
        fm->Geodesic(geodesic);
        
        {
            std::ofstream data_out;
            const std::string filename = testPrefix+"ConstantGeodesic.txt";
            data_out.open(filename.c_str());
            PrintRangeAsList(data_out,geodesic.begin(),geodesic.end());
        }
        
        {
            typedef itk::FixedArray<ScalarType,3> NormExportType;
            typedef itk::Image<NormExportType,Dimension> MetricExportType;
            auto caster = itk::CastImageFilter<MetricType, MetricExportType>::New();
            caster->SetInput(hessianToMetricFilter->GetOutput());
            
            
            auto writer = itk::ImageFileWriter<MetricExportType>::New();
            writer->SetInput(caster->GetOutput());
            writer->SetFileName(testPrefix+"ConstantRadiusMetric.hdf5");
            writer->Update();
        }
        
        functor.Show = true;
        functor(hessianFilter->GetOutput()->GetPixel({{12,14}}));
    }
    
    
    // Old
    
    void RunTest_Tree(){
        const ScalarType sigma = 3;
        const ScalarType deltaMinFixed = 0.01;
        const ScalarType eVal1Threshold = 0.1/(sigma*sigma);
        const ScalarType kappa = 10;
        typedef itk::Riemannian2DNorm<ScalarType> NormType;
        
        auto reader = itk::ImageFileReader<ScalarImageType>::New();
        reader->SetFileName(testPrefix+"TreeTest.hdf5");
        reader->UpdateOutputInformation();
        const auto region = reader->GetOutput()->GetLargestPossibleRegion();
        
        typedef itk::HessianRecursiveGaussianImageFilter<ScalarImageType> HessianFilterType;
        typedef HessianFilterType::OutputImageType TensorImageType;
        typedef TensorImageType::PixelType TensorType;
        auto hessianFilter = HessianFilterType::New();
        hessianFilter->SetInput(reader->GetOutput());
        hessianFilter->SetSigma(sigma);
        
        struct HessianToMetricFunctor {
            ScalarType deltaMinFixed;
            ScalarType deltaMinEncountered;
            ScalarType eVal1Threshold;
            ScalarType kappa;
            NormType operator()(const TensorType & H) {
                TensorType::EigenValuesArrayType eVal;
                TensorType::EigenVectorsMatrixType eVec;
                H.ComputeEigenAnalysis(eVal, eVec);
                eVec = eVec.GetTranspose();
                
                if(fabs(eVal[0]) > fabs(eVal[1])){
                    std::swap(eVal[0], eVal[1]);
                    std::swap(eVec(0,0),eVec(0,1));
                    std::swap(eVec(1,0),eVec(1,1));
                }
                
                NormType norm;
                if(eVal[1] <= eVal1Threshold){
                    norm.SetIdentity();
                    return norm;
                }
                
                const ScalarType ratio = fabsf(eVal[0])/eVal[1];
                
                ScalarType delta = vnl_math_sqr(ratio);
                deltaMinEncountered = std::min(delta, deltaMinEncountered);
                delta = std::max(delta,deltaMinFixed);
                const ScalarType gamma = std::min(ScalarType(1),kappa*delta);
                
                norm(0,0) = vnl_math_sqr(delta);
                norm(0,1) = 0;
                norm(1,1) = vnl_math_sqr(gamma);
                norm = NormType( norm.Rotate(eVec) );
                
                return norm;
            }
        };
        
        typedef itk::Image<NormType,Dimension> MetricType;
        auto hessianToMetricFilter = itk::UnaryFunctorImageFilter<TensorImageType, MetricType, HessianToMetricFunctor>::New();
        hessianToMetricFilter->SetInput(hessianFilter->GetOutput());
        auto & functor = hessianToMetricFilter->GetFunctor();
        functor.deltaMinEncountered = 1;
        functor.deltaMinFixed = deltaMinFixed;
        functor.eVal1Threshold = eVal1Threshold;
        functor.kappa = kappa;
        
        typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
        auto fm = FMType::New();
        fm->SetInput(hessianToMetricFilter->GetOutput());
        
        auto seeds = FMType::NodeContainer::New();
        seeds->Initialize();
        FMType::NodeType node;
        node.SetValue(0);
        itk::Index<Dimension> index;
        index[0]=86;
        index[1]=162;
        node.SetIndex(index);
        seeds->InsertElement(0,node);
        fm->SetTrialPoints(seeds);
        
        {
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(fm->GetOutput());
            writer->SetFileName(testPrefix+"TreeDistance.hdf5");
            writer->Update();
        }
        
        std::cout << "deltaMinEncountered : " << hessianToMetricFilter->GetFunctor().deltaMinEncountered << std::endl;

        struct EVFunctor {
            int eVIndex;
            ScalarType operator()(const TensorType & H){
                TensorType::EigenValuesArrayType eVal;
                H.ComputeEigenValues(eVal);
                if(fabs(eVal[0]) > fabs(eVal[1])) std::swap(eVal[0],eVal[1]);
                return eVal[eVIndex];
            }
        };
        
        auto eVFilter = itk::UnaryFunctorImageFilter<TensorImageType, ScalarImageType, EVFunctor>::New();
        eVFilter->GetFunctor().eVIndex = 1;
        eVFilter->SetInput(hessianFilter->GetOutput());
        
        {
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(eVFilter->GetOutput());
            writer->SetFileName(testPrefix+"eV1Tree.hdf5");
            writer->Update();
        }
    }

    void RunTest_Multiscale(){
        
        const ScalarType sigmaMin = 2;
        const ScalarType sigmaMax = 10;
        const int sigmaScales = 5;
        
        const ScalarType maxInverseSpeed = 0.01;
        const ScalarType eVal1ThresholdNormalized = 0.05;
        const ScalarType kappaMax = 10;
        const ScalarType verticalLogarithmicSpeed = 10;
        
        typedef itk::Index<Dimension+1> ExtendedIndexType;
        const ExtendedIndexType seedIndex = {{356,728,sigmaScales-1}};
        const int numberOfTips = 5;
        const ExtendedIndexType tips[numberOfTips] = { {{310,528,sigmaScales-1}} , {{428,49,0}}, {{91,227,0}}, {{625,300,0}}, {{639,425,0}} };

        
        const ScalarType sigmaRatio = pow(sigmaMax/sigmaMin, ScalarType(1)/std::max(sigmaScales-1,1));
        const ScalarType timeBetweenScales = verticalLogarithmicSpeed*log(sigmaRatio);
        itk::FixedArray<ScalarType,sigmaScales> sigmas;
        sigmas[0]=sigmaMin;
        for(int i=1; i<sigmaScales; ++i)
            sigmas[i] = sigmas[i-1]*sigmaRatio;
        std::cout << "Scales" << sigmas << std::endl;
        
        
        typedef unsigned char RGBComponentType;
        typedef itk::RGBPixel<RGBComponentType> RGBPixelType;
        typedef itk::Image<RGBPixelType> RGBImageType;
        auto reader = itk::ImageFileReader<RGBImageType>::New();
        reader->SetFileName(testPrefix+"Real-DSA-01.png");
        reader->UpdateOutputInformation();
        {
            const auto io = reader->GetImageIO();
            assert(io->GetNumberOfDimensions() == Dimension);
            assert(io->GetNumberOfComponents() == 3);
            assert(io->GetComponentType() == itk::ImageIOBase::UCHAR);
        }
        const RGBImageType * physicalImage = reader->GetOutput();
        const auto region = physicalImage->GetLargestPossibleRegion();
        std::cout << "Spacing : " << physicalImage->GetSpacing() << std::endl;

        struct RGBToBWFunctor {
            ScalarType operator()(const RGBPixelType & rgb){
                return (rgb[0]+rgb[1]+rgb[2])/ScalarType(3*itk::NumericTraits<RGBComponentType>::max());
            }
        };

        auto rgb2bwFilter = itk::UnaryFunctorImageFilter<RGBImageType, ScalarImageType, RGBToBWFunctor>::New();
        rgb2bwFilter->SetInput(reader->GetOutput());
        
        { // Export for visualization
            auto writer = itk::ImageFileWriter<ScalarImageType>::New();
            writer->SetInput(rgb2bwFilter->GetOutput());
            writer->SetFileName(testPrefix+"Real-DSA-01.hdf5");
            writer->Update();
        }
        
        typedef itk::Riemannian2DNorm<ScalarType> NormType;
        typedef itk::ExtendedNorm<NormType> ExtendedNormType;
        const unsigned int ExtendedDimension = Dimension+1;
        
        typedef itk::Image<NormType,Dimension> MetricType;
        typedef itk::Image<ExtendedNormType,ExtendedDimension> ExtendedMetricType;
        auto extendedMetric = ExtendedMetricType::New();
        
        itk::Index<ExtendedDimension> extendedIndex;
        itk::Size<ExtendedDimension> extendedSize;
        for(int i=0; i<(int)Dimension; ++i){
            extendedIndex[i] = region.GetIndex(i);
            extendedSize[i] = region.GetSize(i);
        }
        extendedIndex[Dimension]=0;
        extendedSize[Dimension]=sigmaScales;
        extendedMetric->SetRegions(itk::ImageRegion<ExtendedDimension>(extendedIndex,extendedSize));
        extendedMetric->Allocate();
        
        for(int scale=0; scale<sigmaScales; ++scale){
            std::cout << "Assembling metric at scale " << scale << std::endl;
            typedef itk::HessianRecursiveGaussianImageFilter<ScalarImageType> HessianFilterType;
            typedef HessianFilterType::OutputImageType TensorImageType;
            typedef TensorImageType::PixelType TensorType;
            auto hessianFilter = HessianFilterType::New();
            hessianFilter->SetInput(rgb2bwFilter->GetOutput());
            hessianFilter->SetSigma(sigmas[scale]);
            
            struct HessianToMetricFunctor {
                ScalarType deltaMinFixed;
                ScalarType eVal1Threshold;
                ScalarType kappa;
                ScalarType speedMultiplier;
                NormType operator()(const TensorType & H) {
                    TensorType::EigenValuesArrayType eVal;
                    TensorType::EigenVectorsMatrixType eVec;
                    H.ComputeEigenAnalysis(eVal, eVec);
                    eVec = eVec.GetTranspose();
                    
                    if(fabs(eVal[0]) > fabs(eVal[1])){
                        std::swap(eVal[0], eVal[1]);
                        std::swap(eVec(0,0),eVec(0,1));
                        std::swap(eVec(1,0),eVec(1,1));
                    }
                    
                    const ScalarType eValRatio = 1.4;
                    const ScalarType Bonus = 4;
                    NormType norm;
                    norm.SetIdentity();
                    if(eVal[1] <= eVal1Threshold) return norm;
                    if(fabs(eVal[0]) >= fabs(eVal[1])/eValRatio) return norm;
                    
                    const ScalarType tangentSpeed = Bonus*Bonus*vnl_math_sqr(eVal[1]/eVal1Threshold);
                    const ScalarType normalSpeed = tangentSpeed/Bonus;
                    norm(0,0)=1/vnl_math_sqr(tangentSpeed);
                    norm(0,1)=0;
                    norm(1,1)=1/vnl_math_sqr(normalSpeed);
                    norm = NormType( norm.Rotate(eVec) );
//                    norm /= vnl_math_sqr(Bonus*vnl_math_sqr(eVal[1]/eVal1Threshold)); //exp( eVal[1]/eVal1Threshold -ScalarType(1));
                    return norm;
                                        
                    /*
                    const ScalarType ratio = fabsf(eVal[0])/eVal[1];
                    ScalarType delta = vnl_math_sqr(ratio)/speedMultiplier;
                    delta = std::max(delta,deltaMinFixed);
                    const ScalarType gamma = std::min(ScalarType(1),kappa*delta);
                    
                    norm(0,0) = vnl_math_sqr(delta);
                    norm(0,1) = 0;
                    norm(1,1) = vnl_math_sqr(delta);
//                    norm(1,1) = vnl_math_sqr(gamma);
                     */
                    
                }
            };
            
            auto hessianToMetricFilter = itk::UnaryFunctorImageFilter<TensorImageType, MetricType, HessianToMetricFunctor>::New();
            hessianToMetricFilter->SetInput(hessianFilter->GetOutput());
            auto & functor = hessianToMetricFilter->GetFunctor();
            functor.deltaMinFixed = maxInverseSpeed;
            functor.eVal1Threshold = eVal1ThresholdNormalized/vnl_math_sqr(sigmas[scale]);
//            functor.speedMultiplier = exp(sigmas[scale]);
            functor.speedMultiplier = 1;
            functor.kappa = kappaMax;
            
            hessianToMetricFilter->Update();
            
            // Copy the result into the metric
            auto sliceIndex = extendedIndex;
            sliceIndex[Dimension] = scale;
            auto sliceSize = extendedSize;
            sliceSize[Dimension] = 1;
            const itk::ImageRegion<ExtendedDimension> sliceRegion(sliceIndex,sliceSize);
            
            itk::ImageRegionIterator<ExtendedMetricType> extendedIt(extendedMetric,sliceRegion);
            itk::ImageRegionConstIterator<MetricType> it(hessianToMetricFilter->GetOutput(),region);
            for(extendedIt.GoToBegin(),it.GoToBegin();
                !extendedIt.IsAtEnd();
                ++extendedIt,++it){
                extendedIt.Value().GetPrimaryNorm() = it.Value();
                extendedIt.Value().GetScalar() = 1/timeBetweenScales;
            }
            
            { // Export for visualization
                struct EVFunctor {
                    int eVIndex;
                    ScalarType operator()(const TensorType & H){
                        TensorType::EigenValuesArrayType eVal;
                        H.ComputeEigenValues(eVal);
                        if(fabs(eVal[0]) > fabs(eVal[1])) std::swap(eVal[0],eVal[1]);
                        return eVal[eVIndex];
                    }
                };
                
                auto eVFilter = itk::UnaryFunctorImageFilter<TensorImageType, ScalarImageType, EVFunctor>::New();
                eVFilter->GetFunctor().eVIndex = 1;
                eVFilter->SetInput(hessianFilter->GetOutput());
                
                auto writer = itk::ImageFileWriter<ScalarImageType>::New();
                writer->SetInput(eVFilter->GetOutput());
                std::ostringstream s;
                s <<testPrefix << "Multiscale2D_EV1_" << scale << ".hdf5";
                writer->SetFileName(s.str());
                writer->Update();
            }
            
            { // Compute geodesics at this specific scale
                std::cout << "Running local fast marching" << std::endl;
                typedef itk::Index<Dimension> IndexType;
                typedef itk::AnisotropicFastMarchingImageFilter<NormType> FMType;
                auto fm = FMType::New();
                fm->SetInput(hessianToMetricFilter->GetOutput());
                
                auto seeds = FMType::NodeContainer::New();
                seeds->Initialize();
                FMType::NodeType node;
                node.SetValue(0);
                IndexType index;
                index[0]=seedIndex[0];
                index[1]=seedIndex[1];
                node.SetIndex(index);
                seeds->InsertElement(0, node);
                fm->SetTrialPoints(seeds);
                fm->SetGenerateUpwindGradient(true);
                fm->Update();
                
                for(int i=0; i<numberOfTips; ++i){
                    typedef itk::ContinuousIndex<ScalarType,Dimension> ContinuousIndexType;
                    std::vector<ContinuousIndexType> geodesic;
                    ContinuousIndexType tip;
                    tip[0] = tips[i][0];
                    tip[1] = tips[i][1];
                    geodesic.push_back(tip);
                    std::cout << "Extracting geodesic" << std::endl;
                    fm->Geodesic(geodesic);
                    
                    std::ofstream data_out;
                    std::ostringstream s;
                    s << testPrefix << "MultiscaleGeodesic_scale" << scale << "_Tip" << i << ".txt";
                    data_out.open(s.str().c_str());
                    PrintRangeAsList(data_out, geodesic.begin(), geodesic.end());
                    data_out.close();
                } // for i number of tips
            }
            
            {
                // Export metric for visualization
                struct MetricToScalarFunctor {
                    ScalarType operator()(const NormType & m){
                        return m(0,0);
                    }
                };
                auto metricToScalarFilter = itk::UnaryFunctorImageFilter<MetricType, ScalarImageType, MetricToScalarFunctor>::New();
                metricToScalarFilter->SetInput(hessianToMetricFilter->GetOutput());
                
                auto writer = itk::ImageFileWriter<ScalarImageType>::New();
                writer->SetInput(metricToScalarFilter->GetOutput());
                std::ostringstream s;
                s << testPrefix << "MetricAtScale" << scale << "EV1.hdf5";
                writer->SetFileName(s.str());
                writer->Update();
            }
        }
        typedef itk::AnisotropicFastMarchingImageFilter<ExtendedNormType> FMType;
        auto fm = FMType::New();
        fm->SetInput(extendedMetric);
        
        auto seeds = FMType::NodeContainer::New();
        seeds->Initialize();
        FMType::NodeType node;
        node.SetValue(0);
        itk::Index<ExtendedDimension> nodeIndex = seedIndex;
        node.SetIndex(nodeIndex);
        seeds->InsertElement(0,node);
        fm->SetGenerateUpwindGradient(true);
        fm->SetTrialPoints(seeds);
        
        {
            auto writer = itk::ImageFileWriter<itk::Image<ScalarType,ExtendedDimension> >::New();
            writer->SetInput(fm->GetOutput());
            writer->SetFileName(testPrefix+"Multiscale2DDistance.hdf5");
            std::cout << "Running fast marching" << endl;
            writer->Update();
        }
        for(int i=0; i<numberOfTips; ++i){
            typedef itk::ContinuousIndex<ScalarType,ExtendedDimension> ExtendedContinuousIndex;
            std::vector<ExtendedContinuousIndex> geodesic;
            geodesic.push_back(ExtendedContinuousIndex(tips[i]) );
            std::cout << "Extracting geodesic" << std::endl;
            fm->Geodesic(geodesic);
        
            std::ofstream data_out;
            std::ostringstream s;
            s << testPrefix << "MultiscaleGeodesic" << i << ".txt";
            data_out.open(s.str().c_str());
            PrintRangeAsList(data_out, geodesic.begin(), geodesic.end());
            data_out.close();
        } // for i number of tips
    }
}


#endif
