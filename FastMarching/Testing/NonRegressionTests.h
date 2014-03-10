//
//  NRRD_Examples.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 01/02/13.
//
//

#ifndef ITKFM_NonRegressionTests_h
#define ITKFM_NonRegressionTests_h

#include "itkNrrdImageIO.h"
#include "itkMetaImageIO.h"
#include "CommandLineCall2.h"
#include "itkFastMarchingImageFilter.h"
#include <time.h>

namespace ITKFM_NonRegressionTests
{
    void GenerateIdentityTensors();
    template <typename NormType> void ExportIdentityTensor(std::string, int, bool=true);
    int IdentityTest();
    
    // ************************ Test ***************************
    
    int IdentityTest()
    {
        const std::string dirData = "data/Testing";
        const std::string dirOutput = "output/Testing";
        const std::string extensionTensor = "nii";
        const std::string extensionGeodesic = "mathematica";
//        const std::string extensionDistances = "nii";
        const std::string extensionDistances = "mathematica";
        
        const int nIdentifiers = 12;
        const std::string identifiers[nIdentifiers] = {
            "2D_double.",
            "2D_float.",
            "2DAsym_double.",
            "2DAsym_float.",
            "3D_double.",
            "3D_float.",
            "2DTubular_double.",
            "2DTubular_float.",
            "2DAsymTubular_double.",
            "2DAsymTubular_float.",
            "3DTubular_double.",
            "3DTubular_float."
        };

        const int argc = 5;
        const int max_arg_len = 255;
        char argv_[argc][max_arg_len];
        char *argv[argc]; for(int i=0; i<argc; ++i) argv[i] = &argv_[i][0];
        
        for(int i=0; i<nIdentifiers; ++i) //nIdentifiers
//        int i=4;
        {
            strcpy(argv[0], "");
            strcpy(argv[1], (dirData+"/Identity"+identifiers[i]+extensionTensor).c_str());
            strcpy(argv[2], "0");
            strcpy(argv[3], (dirOutput+"/Geodesic"+identifiers[i]+extensionGeodesic).c_str());
            strcpy(argv[4], (dirOutput+"/Distances"+identifiers[i]+extensionDistances).c_str());
            
            CommandLineMain(argc, argv);
        }
        return EXIT_SUCCESS;
    }
    
    // ************************* Generate tensors ***************************
    
    void GenerateIdentityTensors()
    {
        const std::string dirname = "data/Testing";
        const std::string extension = "nii";
        
        ExportIdentityTensor<itk::SymmetricSecondRankTensor<double,2> >(dirname+"/Identity2D_double."+extension, 151, false);
        ExportIdentityTensor<itk::SymmetricSecondRankTensor<float ,2> >(dirname+"/Identity2D_float." +extension, 151, false);
        
        ExportIdentityTensor<itk::Finsler2DNorm<double> >(dirname+"/Identity2DAsym_double."+extension, 151);
        ExportIdentityTensor<itk::Finsler2DNorm<float > >(dirname+"/Identity2DAsym_float." +extension, 151);
        
        ExportIdentityTensor<itk::SymmetricSecondRankTensor<double,3> >(dirname+"/Identity3D_double."+extension, 29, false);
        ExportIdentityTensor<itk::SymmetricSecondRankTensor<float ,3> >(dirname+"/Identity3D_float." +extension, 29, false);
        
        ExportIdentityTensor<itk::ExtendedNorm<itk::Riemannian2DNorm<double> > >(dirname+"/Identity2DTubular_double."+extension, 29);
        ExportIdentityTensor<itk::ExtendedNorm<itk::Riemannian2DNorm<float > > >(dirname+"/Identity2DTubular_float." +extension, 29);
        
        ExportIdentityTensor<itk::ExtendedNorm<itk::Finsler2DNorm<double> > >(dirname+"/Identity2DAsymTubular_double."+extension, 29);
        ExportIdentityTensor<itk::ExtendedNorm<itk::Finsler2DNorm<float > > >(dirname+"/Identity2DAsymTubular_float." +extension, 29);
        
        ExportIdentityTensor<itk::ExtendedNorm<itk::Riemannian3DNorm<double> > >(dirname+"/Identity3DTubular_double."+extension, 21);
        ExportIdentityTensor<itk::ExtendedNorm<itk::Riemannian3DNorm<float > > >(dirname+"/Identity3DTubular_float." +extension, 21);
    } // Generate Identity tensors
    
    template <typename NormType> void ExportIdentityTensor(std::string filename, int side_size, bool convert_to_FixedArray)
    {
        const int Dimension = NormType::Dimension;
        itk::Size<Dimension> size;      size.Fill(side_size);
        itk::Index<Dimension> index;    index.Fill(0);
        itk::ImageRegion<Dimension> region(index,size);
        
        
        NormType Identity; Identity.SetIdentity();
        
        if(!convert_to_FixedArray){
            typedef itk::Image<NormType,Dimension> ImageType;
            auto Image = ImageType::New();
            Image->SetRegions(region);
            Image->Allocate();
            Image->FillBuffer(Identity);
            
            auto writer = itk::ImageFileWriter<ImageType>::New();
            writer->SetFileName(filename);
            writer->SetInput(Image);
            try {writer->Update();}
            catch(itk::ExceptionObject &e){std::cout << e << endl;}
            
            
        } else {
            
            typedef typename NormType::ValueType ValueType;
            const int ArraySize = sizeof(NormType)/sizeof(ValueType);
            typedef itk::Vector<ValueType,ArraySize> ArrayType;
            ArrayType ArrayIdentity;
            
            for(int i=0; i<ArraySize; ++i) ArrayIdentity.SetElement(i,Identity.GetNthComponent(i));
            
            typedef itk::Image<ArrayType,Dimension> ImageType;
            auto Image = ImageType::New();
            Image->SetRegions(region);
            Image->Allocate();
            
            Image->FillBuffer(ArrayIdentity);
            
            auto writer = itk::ImageFileWriter<ImageType>::New();
            writer->SetFileName(filename);
            writer->SetInput(Image);
            try {writer->Update();}
            catch(itk::ExceptionObject &e){std::cout << e << endl;}
        }
        
    } // Export Identity tensor
    
    template<unsigned int Dimension>
    void GeneratePointPairs()
    {
        std::stringstream identifier_; identifier_ << Dimension << "D";
        const std::string identifier = identifier_.str();
        
        typedef itk::Point<double,Dimension> PointType;
        std::vector<PointType> Points;
        
        PointType Origin; Origin.Fill(0.);
        PointType Ten; Ten.Fill(10.);
        Points.push_back(Origin);
        Points.push_back(Ten);
        
        CommandLineMain_Subroutines::ExportPointListToFile("data/Testing/Point_Origin_Ten_"+identifier+".txt",Points);
//        CommandLineMain_Subroutines::ExportPointListToFile("data/Testing/Point_Origin_Ten_"+identifier+".mathematica",Points);
        CommandLineMain_Subroutines::ExportPointListToFile("data/Testing/Point_Origin_Ten_"+identifier+".nii",Points);
        
        std::vector<PointType> ImportedPoints;
        CommandLineMain_Subroutines::ImportPointListFromFile("data/Testing/Point_Origin_Ten_"+identifier+".txt", ImportedPoints);
        CommandLineMain_Subroutines::ImportPointListFromFile("data/Testing/Point_Origin_Ten_"+identifier+".nii", ImportedPoints);
        for(auto it=ImportedPoints.begin(); it!=ImportedPoints.end(); ++it) {cout << *it << endl;}
    }

    template<typename NormType>
    void CompareTiming(int side_size) // in the basic euclidean case
    { 
        const int Dimension = NormType::Dimension;
        typedef typename NormType::ValueType ValueType;
        
        cout << "Dimension : " << Dimension << ", side_size : " << side_size << endl;
        
        typedef itk::FastMarchingImageFilter<itk::Image<ValueType,Dimension> > IsoFMType;
        typedef itk::AnisotropicFastMarchingImageFilter<NormType> AnisoFMType;
        
        itk::Size<Dimension> size;
        size.Fill(side_size);
        itk::Index<Dimension> index;
        index.Fill(0);
        const itk::ImageRegion<Dimension> region(index,size);
        
        typename IsoFMType::NodeType node;
        node.SetValue(0.);
        node.SetIndex(index);
        auto seeds = IsoFMType::NodeContainer::New();
        seeds->Initialize();
        seeds->InsertElement(0,node);
        
        {
            NormType Identity;
            Identity.SetIdentity();
            
            auto Metric = itk::Image<NormType,Dimension>::New();
            Metric->SetRegions(region);
            Metric->Allocate();
            Metric->FillBuffer(Identity);
            
            auto FM = AnisoFMType::New();
            FM->SetInput(Metric);
            FM->SetTrialPoints(seeds);

            clock_t timing = -clock();
            FM->Update();
            timing+=clock();
            cout << "Anisotropic Fast Marching took : " << timing/double(CLOCKS_PER_SEC) << " seconds" << endl;
        }
        
        {
            auto FM = IsoFMType::New();
            FM->SetOutputRegion(region);
            FM->SetTrialPoints(seeds);
            
            clock_t timing=-clock();
            FM->Update();
            timing+=clock();
            cout << "Isotropic Fast Marching took : " << timing/double(CLOCKS_PER_SEC) << " seconds" << endl;
        }
    }

    void CompareTimings()
    {
        cout << "itk::Riemannian2DNorm<float>" << endl;
        CompareTiming< itk::Riemannian2DNorm<float> >(1000);
        
        cout << "itk::Finsler2DNorm<float>" << endl;
        CompareTiming< itk::Finsler2DNorm<float> >(1000);

        cout << "itk::ExtendedNorm<itk::Riemannian2DNorm<float> > " << endl;
        CompareTiming< itk::ExtendedNorm<itk::Riemannian2DNorm<float> > >(100);

        cout << "itk::Riemannian3DNorm<float>" << endl;
        CompareTiming< itk::Riemannian3DNorm<float> >(100);
        
        cout << "itk::ExtendedNorm<itk::Riemannian3DNorm<float> > " << endl;
        CompareTiming< itk::ExtendedNorm<itk::Riemannian3DNorm<float> > >(30);
    }
    
} // namespace GenerateTestData

#endif
