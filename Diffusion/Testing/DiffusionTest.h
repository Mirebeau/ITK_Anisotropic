//
//  DiffusionTest2.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 07/03/2014.
//
//

#ifndef itkDiffusion_DiffusionTest2_h
#define itkDiffusion_DiffusionTest2_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRandomImageSource.h"
#include "itkGaussianDistribution.h"
#include "itkPermuteAxesImageFilter.h"


#include "SyntheticImageSource.h"
#include "CoherenceEnhancingDiffusionFilter.h"


namespace Testing {
    using namespace itk;
    using std::cout;
    using std::endl;
    const std::string imageSourceDirectory = "../Images/";
    const std::string FPExtension = ".vtk"; // or hdf5
    
    typedef double ValueType;
    
    void CoherenceEnhancingDiffusion2D(){ // test on a fingerprint image
        
        const std::string testName = "CoherenceEnhancingDiffusion2D";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        // Import Image, construct diffusion tensor
        
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+testName+"_TestImage.png");
        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> AnisotropicDiffusionFilterType;
        AnisotropicDiffusionFilterType::Pointer anisotropicDiffusionFilter = AnisotropicDiffusionFilterType::New();
        anisotropicDiffusionFilter->SetInput(reader->GetOutput());
        
        anisotropicDiffusionFilter->SetEnhancement(AnisotropicDiffusionFilterType::CED);
        anisotropicDiffusionFilter->SetDiffusionTime(10);
        anisotropicDiffusionFilter->SetExponent(2);
        anisotropicDiffusionFilter->SetNoiseScale(0.5);
        anisotropicDiffusionFilter->SetFeatureScale(4);
        anisotropicDiffusionFilter->SetLambda(4);
        anisotropicDiffusionFilter->SetAlpha(0.01);
        anisotropicDiffusionFilter->SetAdimensionize(false);
        
        {
            typedef ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(testName+FPExtension);
            writer->SetInput(anisotropicDiffusionFilter->GetOutput());
            writer->Update();
        }
        cout << "Finishing " << testName << endl;
    }
    
    
    void EdgeEnhancingDiffusion2D(){
        
        const std::string testName = "EdgeEnhancingDiffusion2D";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        // Import Image, construct diffusion tensor
        
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+testName+"_TestImage.png");

        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> AnisotropicDiffusionFilterType;
        AnisotropicDiffusionFilterType::Pointer anisotropicDiffusionFilter = AnisotropicDiffusionFilterType::New();
        anisotropicDiffusionFilter->SetInput(reader->GetOutput());
        
        anisotropicDiffusionFilter->SetEnhancement(AnisotropicDiffusionFilterType::cEED);
        anisotropicDiffusionFilter->SetDiffusionTime(10);
        anisotropicDiffusionFilter->SetExponent(2);
        anisotropicDiffusionFilter->SetNoiseScale(4);
        anisotropicDiffusionFilter->SetFeatureScale(4);
        anisotropicDiffusionFilter->SetLambda(4);
        anisotropicDiffusionFilter->SetAlpha(0.01);
        anisotropicDiffusionFilter->SetAdimensionize(false);

        
        typedef ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(testName+FPExtension);
        writer->SetInput(anisotropicDiffusionFilter->GetOutput());
        writer->Update();
        
        cout << "Finishing " << testName << endl;
    }
    
    
    struct EED3D_FunctorType {
        typedef Image<ValueType,3> ImageType;

        ImageType * physicalData;
        Statistics::GaussianDistribution * gaussian;
        ValueType operator()(ValueType rand, const ImageType::IndexType & index) const {
            ImageType::PointType x;
            physicalData->TransformIndexToPhysicalPoint(index, x);
            ValueType r = x.GetVectorFromOrigin().GetNorm();
            
            return cos( vnl_math_cube(r/R) ) + gaussian->EvaluateInverseCDF(rand);
        }
        ValueType R;
        EED3D_FunctorType(){R=0.5;}
    };
    
    struct CastFunctorType {
        unsigned char operator()(ValueType r) const {
            return (unsigned char) 255*std::min(std::max((1.+r)/2.,0.),1.);
        }
    };
    
    void EdgeEnhancingDiffusion3D(){
        const unsigned int Dimension = 3;
        
        const std::string testName = "EdgeEnhancingDiffusion3D";
        cout << "Beginning " << testName << endl;
        
        const long n=100;  // Number of pixels in each dimension
        const ValueType h = 1./n;
        const ValueType Variance = 0.5;
                
        typedef Image<ValueType,Dimension> ImageType;
        
        // Creating the region
        typedef ImageRegion<Dimension> RegionType;
        
        typedef RandomImageSource<ImageType> RandomImageSourceType;
        RandomImageSourceType::Pointer randomSource = RandomImageSourceType::New();
        randomSource->SetMax(1);
        randomSource->SetMin(0);
        randomSource->SetNumberOfThreads(1); // for reproducibility
        Size<Dimension> size;           size.Fill(n);       randomSource->SetSize(size);
        ImageType::SpacingType spacing; spacing.Fill(1./n); randomSource->SetSpacing(spacing);
        // index is automatically set to [0,0,0]
        
        typedef UnaryFunctorWithIndexImageFilter<ImageType, ImageType, EED3D_FunctorType> FunctorFilterType;
        FunctorFilterType::Pointer functorFilter =  FunctorFilterType::New();
        
        typedef Statistics::GaussianDistribution DistributionType;
        DistributionType::Pointer gaussian = DistributionType::New();
        gaussian->SetMean(0.);
        gaussian->SetVariance(Variance);
        functorFilter->GetFunctor().gaussian = gaussian;
        
        randomSource->Update();
        functorFilter->GetFunctor().physicalData = randomSource->GetOutput();
        functorFilter->SetInput(randomSource->GetOutput());
        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> DiffusionFilterType;
        DiffusionFilterType::Pointer diffusionFilter = CoherenceEnhancingDiffusionFilter<ImageType>::New();
        diffusionFilter->SetInput(functorFilter->GetOutput()); //functorFilter->GetOutput()
        
        // Do not forget to take the image scale h into account
        diffusionFilter->SetEnhancement(DiffusionFilterType::cEED);
        diffusionFilter->SetDiffusionTime(10*h*h);
        diffusionFilter->SetNoiseScale(5*h);
        diffusionFilter->SetFeatureScale(20*h);
        diffusionFilter->SetLambda(1);
        diffusionFilter->SetAlpha(0.01);
        diffusionFilter->SetAdimensionize(false);

        
        typedef ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput(diffusionFilter->GetOutput());
        writer->SetFileName("Cos3D_Smoothed"+FPExtension);
//        writer->Update();
        
        typedef Image<unsigned char,Dimension> OutputImageType;
        typedef UnaryFunctorImageFilter<ImageType, OutputImageType, CastFunctorType> CastFilterType;
        CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(functorFilter->GetOutput());
        
        typedef ImageFileWriter<OutputImageType> OutputWriterType;
        {
            OutputWriterType::Pointer outputWriter = OutputWriterType::New();
            outputWriter->SetInput(castFilter->GetOutput());
            outputWriter->SetFileName("Cos3D_Noisy"+FPExtension);
            outputWriter->Update();
        }
        
        {
            functorFilter->GetFunctor().gaussian->SetVariance(0.);
            functorFilter->Modified();
            OutputWriterType::Pointer outputWriter = OutputWriterType::New();
            outputWriter->SetFileName("Cos3D_Source"+FPExtension);
            outputWriter->SetInput(castFilter->GetOutput());
            outputWriter->Update();
        }
        
/*
        {
            auto writer = ImageFileWriter<ImageType>::New();
            writer->SetInput(functorFilter->GetOutput());
            writer->SetFileName(testName+"_Noisy.hdf5");
            writer->Update();
        }
*/
        // Export a cropped 8 bit version.
        
        
        
        cout << "finishing " << testName << endl;
    }
    
    // Tests by Jerome Fehrenbach.
    
    void CoherenceLena(){
        const std::string testName = "CoherenceLena";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        // Import Image, construct diffusion tensor
        
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+"LenaBW.png");
        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> AnisotropicDiffusionFilterType;
        AnisotropicDiffusionFilterType::Pointer anisotropicDiffusionFilter = AnisotropicDiffusionFilterType::New();
        anisotropicDiffusionFilter->SetInput(reader->GetOutput());
        
        anisotropicDiffusionFilter->SetEnhancement(AnisotropicDiffusionFilterType::CED);
        anisotropicDiffusionFilter->SetDiffusionTime(50);
        anisotropicDiffusionFilter->SetExponent(2);
        anisotropicDiffusionFilter->SetNoiseScale(0.5);
        anisotropicDiffusionFilter->SetFeatureScale(4);
        anisotropicDiffusionFilter->SetLambda(sqrt(1e-5));
        anisotropicDiffusionFilter->SetAlpha(0.01);
        anisotropicDiffusionFilter->SetAdimensionize(false);

        {
            typedef ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(testName+FPExtension);
            writer->SetInput(anisotropicDiffusionFilter->GetOutput());
            writer->Update();
        }
        cout << "Finishing " << testName << endl;
    }
    
    void EdgeSynthetic(){
        const std::string testName = "EdgeSynthetic";
        cout << "Beginning " << testName << endl;
        
        const unsigned int Dimension = 2;
        
        typedef Image<ValueType,Dimension> ImageType;
        
        // Import Image, construct diffusion tensor
        
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+"SyntheticImage.png");
        
        typedef CoherenceEnhancingDiffusionFilter<ImageType> AnisotropicDiffusionFilterType;
        AnisotropicDiffusionFilterType::Pointer anisotropicDiffusionFilter = AnisotropicDiffusionFilterType::New();
        anisotropicDiffusionFilter->SetInput(reader->GetOutput());
        
        // scaling : 140
        
        anisotropicDiffusionFilter->SetEnhancement(AnisotropicDiffusionFilterType::cEED);
        anisotropicDiffusionFilter->SetDiffusionTime(5);
        anisotropicDiffusionFilter->SetExponent(4);
        anisotropicDiffusionFilter->SetNoiseScale(0.5);
        anisotropicDiffusionFilter->SetFeatureScale(4);
        anisotropicDiffusionFilter->SetLambda(200.);
        anisotropicDiffusionFilter->SetAlpha(0.02);
        anisotropicDiffusionFilter->SetAdimensionize(false);

        {
            typedef ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(testName+FPExtension);
            writer->SetInput(anisotropicDiffusionFilter->GetOutput());
            writer->Update();
        }
        cout << "Finishing " << testName << endl;
    }
    
    
    void Color(){
        
        typedef float ScalarType;
        typedef Vector<float,3> PixelType;
        typedef Image<PixelType,2> ImageType;
        typedef SymmetricSecondRankTensor<float,2> TensorType;
        typedef Image<TensorType,2> TensorImageType;
        {
            typedef StructureTensorImageFilter<ImageType,TensorImageType> FilterType;
            FilterType::Pointer filter = FilterType::New();
        }

        {
            typedef LinearAnisotropicDiffusionLBRImageFilter<ImageType,ScalarType> FilterType;
            FilterType::Pointer filter = FilterType::New();
        }

        {
            typedef AnisotropicDiffusionLBRImageFilter<ImageType,ScalarType> FilterType;
            FilterType::Pointer filter = FilterType::New();

        }
        
        
//        typedef GradientRecursiveGaussianImageFilter<ImageType> GradientFilterType;
//        GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
        typedef CoherenceEnhancingDiffusionFilter<ImageType,ScalarType> DiffusionFilterType;
        typename DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
                
        typedef ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imageSourceDirectory+"Lena.png");
        reader->Update();
        std::cout << * reader->GetOutput()->GetBufferPointer() << "\n";

        diffusionFilter->SetDiffusionTime(20);
        diffusionFilter->SetInput(reader->GetOutput());
        diffusionFilter->Update();
        
        typedef Image<Vector<unsigned char, 3> > ExportImageType;
        typedef CastImageFilter<ImageType, ExportImageType> CasterType;
        CasterType::Pointer caster = CasterType::New();
        caster->SetInput(diffusionFilter->GetOutput());
        
        
        typedef ImageFileWriter<ExportImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName("Lena_smoothed.png");
        writer->SetInput(caster->GetOutput());
        writer->Update();
        
    }
    
    
    struct AffineFunctorType {
        static const int Dimension=2;
        Vector<double, Dimension> gradient;
        Point<double, Dimension> center;
        double operator()(const Point<double,Dimension> & p){
            return gradient * (p-center);
        }
/*        double radius;
        bool disc = false;
        
        double operator()(const Point<double,Dimension> & p){
            return disc ? (p-center).GetNorm()<radius : gradient * (p-center);
        }*/
        AffineFunctorType(){gradient.Fill(1.); center.Fill(0.);}
    };
    
    void NonUniformSpacing(){
        const int Dimension=2;
        const int n=21;
        typedef Image<ValueType,Dimension> ImageType;
        // Objective : check that the method remains reasonable with non-uniform spacing.
        typedef SyntheticImageSource<ImageType, AffineFunctorType> ImageSourceType;
        ImageSourceType::Pointer imageSource = ImageSourceType::New();
        
        ImageType::SpacingType spacing;
        spacing.Fill(1.);
        spacing[0]=2.;
        imageSource->GetPhysicalData()->SetSpacing(spacing);
        std::cout << spacing << std::endl;
        
        ImageType::RegionType region;
        Size<Dimension> size;       size[0]=n;  size[1]=2*n;      region.SetSize(size);
        Index<Dimension> index;     index.Fill(0);      region.SetIndex(index);
        imageSource->GetPhysicalData()->SetRegions(region);

        Index<Dimension> index0 = index, index1=index, indexCenter;
        index0[0]=1; index1[1]=1; indexCenter[0]=n/2; indexCenter[1]=n;

        imageSource->Update();
        ImageType::Pointer image = imageSource->GetOutput();
        std::cout
        << image->GetPixel(index)  << ", "
        << image->GetPixel(index0) << ", "
        << image->GetPixel(index1) << ".\n";
        
        // Checking structure tensors.
        
        typedef StructureTensorImageFilter<ImageType> StructureTensorFilterType;
        StructureTensorFilterType::Pointer structureTensorFilter = StructureTensorFilterType::New();
        structureTensorFilter->SetInput(imageSource->GetOutput());
        structureTensorFilter->Update();
        
        std::cout << structureTensorFilter->GetOutput()->GetPixel(indexCenter) << "\n";
        
        
        // Now checking that diffusion is well rescaled. Seems conclusive.
        
        image->FillBuffer(0);
        image->GetPixel(indexCenter)=1;
        typedef LinearAnisotropicDiffusionLBRImageFilter<ImageType> DiffusionFilterType;
        
        typedef DiffusionFilterType::TensorImageType TensorImageType;
        TensorImageType::Pointer tensorImage = TensorImageType::New();
        tensorImage->CopyInformation(image);
        tensorImage->SetRegions(image->GetRequestedRegion());
        tensorImage->Allocate();
        TensorImageType::PixelType tensor;
        tensor.SetIdentity();
        tensorImage->FillBuffer(tensor);
        
        DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
        diffusionFilter->SetInputImage(image);
        diffusionFilter->SetInputTensor(tensorImage);
        diffusionFilter->SetMaxNumberOfTimeSteps(100);
        diffusionFilter->SetMaxDiffusionTime(20);
        
        typedef ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput(diffusionFilter->GetOutput());
        writer->SetFileName("Non-UniformSpacing"+FPExtension);
        writer->Update();
    }
    
    
    struct SliceFunctorType {
        typedef unsigned char PixelType;
        typedef Image<PixelType, 2> SliceType;
        std::vector<SliceType::Pointer> slices;
        typedef Index<3> IndexType;
        typedef Index<2> SliceIndexType;
        Statistics::GaussianDistribution * gaussian;
        
        PixelType operator()(float rand, const IndexType & index) const {
            SliceIndexType sIndex;
            for(int i=0; i<2; ++i) sIndex[i]=index[i];
            int s=index[2];
            assert(0<=s && s<= (int) slices.size());
//            return slices[s]->GetPixel(sIndex);

            double value = slices[s]->GetPixel(sIndex) / double(255) + gaussian->EvaluateInverseCDF(rand);
            value = std::max(0., std::min(1., value));
            return (unsigned char) (255*value);
        }
    };
    
    void Assemble3DSkull(){
        const int nSlices = 99;
        const double Variance=0.01;
        typedef unsigned char PixelType;
        typedef Image<double,3> RandomImageType;
        typedef Image<PixelType,3> ImageType;
        typedef Image<PixelType,2> SliceType;
        typedef ImageFileReader<SliceType> ReaderType;
        

        
        typedef UnaryFunctorWithIndexImageFilter<RandomImageType, ImageType, SliceFunctorType> SliceFilterType;
        SliceFilterType::Pointer sliceFilter = SliceFilterType::New();
        SliceFunctorType & sliceFunctor = sliceFilter->GetFunctor();
        sliceFunctor.slices.resize(nSlices);
        
        typedef Statistics::GaussianDistribution DistributionType;
        DistributionType::Pointer gaussian = DistributionType::New();
        gaussian->SetMean(0.);
        gaussian->SetVariance(Variance);
        sliceFunctor.gaussian = gaussian;

        
        for(int i=1; i<=nSlices; ++i){
            std::string filename = (i<10 ? "mrbrain00" : "mrbrain0") + std::to_string(i) + ".tif";
            
            ReaderType::Pointer reader = ReaderType::New();
            reader->SetFileName(imageSourceDirectory+"resultatedge/"+filename);
            reader->UpdateOutputInformation();
            reader->Update();
            
            sliceFunctor.slices[i-1] = reader->GetOutput();
        }
        
        
        SliceType::Pointer slice = sliceFunctor.slices[0];
        SliceType::RegionType sliceRegion = slice->GetRequestedRegion();
        
        Size<3> size;
//        ImageType::RegionType imageRegion;
        for(int i=0; i<2; ++i){
            size[i]=sliceRegion.GetSize(i);
            assert(sliceRegion.GetIndex(i)==0);
//            imageRegion.SetSize(i, sliceRegion.GetSize(i));
//            imageRegion.SetIndex(i, sliceRegion.GetIndex(i));
        }
        size[2]=nSlices;
//        imageRegion.SetSize(2, nSlices);
//        imageRegion.SetIndex(2, 0);
        
        typedef RandomImageSource<RandomImageType> RandomImageSourceType;
        RandomImageSourceType::Pointer randomSource = RandomImageSourceType::New();
        randomSource->SetMax(1);
        randomSource->SetMin(0);
        randomSource->SetSize(size);
        ImageType::SpacingType spacing; spacing.Fill(1.); randomSource->SetSpacing(spacing);

        randomSource->Update();
        
//        ImageType::Pointer image = ImageType::New();
//        image->SetRegions(imageRegion);
//        image->Allocate();
        
        sliceFilter->SetInput(randomSource->GetOutput());
        sliceFilter->Update();
        
        {
            typedef ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput(sliceFilter->GetOutput());
            writer->SetFileName("mrbrain_noisy"+FPExtension);
            writer->Update();
        }
        
        {
            typedef ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput(sliceFilter->GetOutput());
            writer->SetFileName("mrbrain"+FPExtension);
            sliceFunctor.gaussian->SetVariance(0.);
            sliceFilter->Modified();
            writer->Update();
        }

        
        /*
        const ImageIOBase * io = reader->GetImageIO();
//        const int ImageDimension = io->GetNumberOfDimensions();
        const itk::ImageIOBase::IOComponentType componentType = io->GetComponentType();
//        const int nComponents = io->GetNumberOfComponents();

        std::cout << componentType << ", " << (componentType==itk::ImageIOBase::UCHAR) << "\n";*/
        
    }
    
    struct PolygonFunctorType {
        ValueType s;
        static const int Dimension=2;
        typedef Point<ValueType,Dimension> PointType;
        typedef Vector<ValueType,Dimension> VectorType;
        
        std::vector<PointType> points;
        Statistics::GaussianDistribution * gaussian;
        
        bool Inside(const PointType & x) const {
            for(int i=0; i<(int)points.size(); ++i){
                const PointType p = points[i], q = points[(i+1)%points.size()];
                const VectorType v = q-p;
                VectorType w; w[0]=-v[1]; w[1]=v[0];
                if(w*(x-p) < 0) return false;
            }
            return true;
        }
        unsigned char operator()(ValueType rand, const Index<Dimension> & index) const {
            PointType x;
            for(int i=0; i<Dimension; ++i)
                x[i] = s*index[i];
            ValueType val = (Inside(x) ? 0.2 : 0.8) + gaussian->EvaluateInverseCDF(rand);
            return (unsigned char) (255*std::max(ValueType(0), std::min(ValueType(1), val)));
        }
    };
    
    void SyntheticTriangle(){
        const int Dimension=2;
        const ValueType xyRatio = 3./5.;
        const int n=40;
        const double Variance = 0.02;
        Size<Dimension> size;
        size[0]=n;
        size[1]=(int)ceil(n*xyRatio);
        
        typedef Image<ValueType,2> ImageType;
        typedef RandomImageSource<ImageType> SourceType;
        SourceType::Pointer randomSource = SourceType::New();
        randomSource->SetSize(size);
        randomSource->SetMax(1.);
        randomSource->SetMin(0.);
        
        typedef Image<unsigned char, 2> OutputImageType;
        typedef UnaryFunctorWithIndexImageFilter<ImageType,OutputImageType,PolygonFunctorType> FunctorFilterType;
        FunctorFilterType::Pointer functorFilter = FunctorFilterType::New();
        functorFilter->SetInput(randomSource->GetOutput());
        PolygonFunctorType & functor = functorFilter->GetFunctor();
        
        functor.s = 5./n;
        const int nPts = 3;
        double pts[nPts][2] = {{1.24457, 0.337341}, {4.11057, 1.2239}, {0.579645, 2.48685}};
        functor.points.resize(nPts);
        for(int i=0; i<nPts; ++i)
            for(int j=0; j<Dimension; ++j)
                functor.points[i][j] = pts[i][j];
        
        double barycenter_[2] = {1.97826, 1.34936};
        Point<ValueType,Dimension> barycenter(barycenter_);
        std::cout << "Bary inside: " << functor.Inside(barycenter) << std::endl;
        
        typedef Statistics::GaussianDistribution DistributionType;
        DistributionType::Pointer gaussian = DistributionType::New();
        gaussian->SetMean(0.);
        gaussian->SetVariance(Variance);
        functor.gaussian = gaussian;

        typedef PermuteAxesImageFilter<OutputImageType> PermuteFilterType;
        PermuteFilterType::Pointer permuteFilter = PermuteFilterType::New();
        FixedArray<unsigned int, 2> order;
        order[0] = 1; order[1] = 0;
        permuteFilter->SetOrder(order);
        permuteFilter->SetInput(functorFilter->GetOutput());
        
        typedef ImageFileWriter<OutputImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput(permuteFilter->GetOutput());
        writer->SetFileName("Triangle.png");
        writer->Update();
        
    }
    
    
    void SyntheticOscillations(){
        static const int Dimension=2;
        const int n=50;
        typedef double ScalarType;
        typedef unsigned char OutputPixelType;
        typedef Point<ScalarType,Dimension> PointType;
        typedef Vector<ScalarType, Dimension> VectorType;
        
        typedef Statistics::GaussianDistribution DistributionType;
        DistributionType::Pointer gaussian = DistributionType::New();
        gaussian->SetMean(0.);

        struct OscillationsFunctorType {
            ScalarType k_p, k_m;
            VectorType v;
            DistributionType * gaussian;
            OutputPixelType operator()(const PointType & p) const {
                VectorType w; w[0]=-v[1]; w[1]=v[0];
                const VectorType z = p.GetVectorFromOrigin();
                const ScalarType sv=v*z, sw = w*z;
                const ScalarType k = sv>0 ? k_p : k_m;
                
                ScalarType result =
                (1+cos(k*sv+sw))/2. + gaussian->EvaluateInverseCDF(random()/double(RAND_MAX));
                result = std::max(0.,std::min(1.,result));
                return 255*result;
            }
        };

        typedef Image<OutputPixelType,Dimension> ImageType;
        typedef SyntheticImageSource<ImageType, OscillationsFunctorType> ImageSourceType;
        ImageSourceType::Pointer imageSource = ImageSourceType::New();
        OscillationsFunctorType & functor = imageSource->GetFunctor();
        functor.gaussian = gaussian;
        functor.v[0]=1.; functor.v[1]=-0.3;
        functor.k_m=-1.2;
        functor.k_p=2.1;
        
        
        ImageType::SpacingType spacing;
        spacing.Fill(25./n);
        imageSource->GetPhysicalData()->SetSpacing(spacing);
        
        Size<Dimension> size; size[0]=(2*n)/3; size[1]=n;
        Index<Dimension> index;
        for(int i=0; i<Dimension; ++i)
            index[i]=-(int)size[i]/2;
        ImageType::RegionType region;
        region.SetSize(size);
        region.SetIndex(index);
        imageSource->GetPhysicalData()->SetRegions(region);
        
        for(int i=1; i<=2; ++i){
            const ScalarType variance = 0.1*i;
            gaussian->SetVariance(variance);
            imageSource->Modified();
            
            typedef ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput(imageSource->GetOutput());
            writer->SetFileName("Oscillations_Noisy"+std::to_string(i)+".png");
            writer->Update();
        }
    }
    
    void CircleVectorField(){
        const int Dimension=2;
        const int n=50;
        typedef float ScalarType;
        const ScalarType variance = 1.4;
        
        typedef Statistics::GaussianDistribution DistributionType;
        DistributionType::Pointer gaussian = DistributionType::New();
        gaussian->SetMean(0.);
        gaussian->SetVariance(variance);
        
        typedef Point<ScalarType,Dimension> PointType;
        typedef Vector<ScalarType, Dimension> VectorType;
        struct NoisyTangentVector {
            DistributionType * gaussian;
            VectorType operator()(const PointType & p) const {
                const ScalarType norm = p.GetVectorFromOrigin().GetNorm();
                VectorType result;
                result[0]= norm==0 ? 0 : -p[1]/norm;
                result[1]= norm==0 ? 0 : p[0]/norm;
                for(int i=0; i<Dimension; ++i)
                    result[i]+=gaussian->EvaluateInverseCDF(random()/double(RAND_MAX));
                return norm>1 ? result : -result;
            }
        };
        
        typedef Image<VectorType,Dimension> VectorImageType;
        typedef SyntheticImageSource<VectorImageType, NoisyTangentVector> VectorSourceType;
        VectorSourceType::Pointer vectorSource = VectorSourceType::New();
        vectorSource->GetFunctor().gaussian = gaussian;
        
        typedef SymmetricSecondRankTensor<ScalarType,Dimension> TensorType;
        struct CircleTangentTensors {
            enum {Constant, Isotropic, Anisotropic} choice;
            ScalarType alpha=0.01;
            ScalarType lambda=0.25;
            ScalarType spacing;
            TensorType operator()(const PointType & p) const {
                const ScalarType norm = p.GetVectorFromOrigin().GetNorm();
                const ScalarType ev = std::min(std::max(ScalarType(fabs(norm-ScalarType(1)))/lambda,alpha),ScalarType(1));
                TensorType result;
                result.SetIdentity();
                switch (choice) {
                    case Constant: break;
                    case Isotropic: result*=ev; break;
                    case Anisotropic:{
                        const VectorType q = p.GetVectorFromOrigin() / (norm==0 ? 1 : norm);
                        for(int i=0; i<Dimension; ++i)
                            for(int j=i; j<Dimension; ++j)
                                result(i,j) = (i==j) - (1.-ev)*q[i]*q[j];
                    }
                }
                return result*spacing*spacing;
            }
        };
        typedef Image<TensorType,Dimension> TensorImageType;
        typedef SyntheticImageSource<TensorImageType, CircleTangentTensors> TensorSourceType;
        TensorSourceType::Pointer tensorSource = TensorSourceType::New();
        
        VectorImageType::SpacingType spacing;
        spacing.Fill(2.6/n);
        vectorSource->GetPhysicalData()->SetSpacing(spacing);
        tensorSource->GetPhysicalData()->SetSpacing(spacing);
        
        Size<Dimension> size; size.Fill(n);
        Index<Dimension> index; index.Fill(-n/2);
        VectorImageType::RegionType region;
        region.SetSize(size);
        region.SetIndex(index);
        vectorSource->GetPhysicalData()->SetRegions(region);
        tensorSource->GetPhysicalData()->SetRegions(region);
        
        typedef ImageFileWriter<VectorImageType> VectorWriterType;
        VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
        vectorWriter->SetInput(vectorSource->GetOutput());
        vectorWriter->SetFileName("VectorField_CircleOpposites"+FPExtension);
        vectorWriter->Update();
        
        
        for(int i=0; i<3; ++i){
            typedef ImageFileWriter<TensorImageType> TensorWriterType;
            TensorWriterType::Pointer tensorWriter = TensorWriterType::New();
            tensorWriter->SetInput(tensorSource->GetOutput());
            tensorSource->GetFunctor().choice =
            i==0 ? CircleTangentTensors::Constant :
            i==1 ? CircleTangentTensors::Isotropic :
            CircleTangentTensors::Anisotropic;
            const std::string filename = (std::string)"TensorField_Circle_" +
            (i==0 ? "Constant" : i==1 ? "Isotropic" : "Anisotropic" ) +
            FPExtension;
            tensorSource->GetFunctor().spacing=spacing[0];
            tensorSource->Modified();
            tensorWriter->SetFileName(filename);
            tensorWriter->Update();
        }
        
        
        
        
        

        
        
        
    }
    
    
}
#endif
