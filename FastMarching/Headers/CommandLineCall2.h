//
//  CommandLineCall2.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 31/01/13.
//
//

#ifndef ITKFM_CommandLineCall2_h
#define ITKFM_CommandLineCall2_h

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>

#include "Riemannian2DNorm.h"
#include "Finsler2DNorm.h"
#include "Riemannian3DNorm.h"
#include "ExtendedNorm.h"

#include "itkAnisotropicFastMarchingImageFilter.h"
#include "ExportToMathematica.h"

namespace CommandLineMain_Subroutines {
    void Usage();
    template <typename PixelType> int Execute(int,char**);
    
    template<typename PointType>
    int ImportPointListFromFile(std::string,std::vector<PointType> &);
    
    template<typename PointType>
    int ExportPointListToFile(std::string,const std::vector<PointType> &);
    
    std::string GetExtension(std::string fn){return fn.substr(fn.find_last_of(".") + 1);}
    using std::cerr;
}

int CommandLineMain(int argc, char * argv[])
{
    using std::cerr;
    using std::endl;
    using namespace itk;
    using namespace CommandLineMain_Subroutines;
    
    if(argc<3+1) {Usage(); return EXIT_FAILURE;}
    const int maxDimension = 4;
    const char *ImageFileName =argv[0+1];
    auto reader = ImageFileReader<Image<unsigned char,maxDimension> >::New();
    reader->SetFileName(ImageFileName);
    
	reader->UpdateOutputInformation(); //try, catch ?
	
    const auto io = reader->GetImageIO();
    
    const int NormInternalDimension = io->GetNumberOfComponents();
    const int ImageDimension = io->GetNumberOfDimensions();
    
    const bool is_a_SSRT = (io->GetPixelType() == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR);
    auto ComponentType = io->GetComponentType();
    enum {FLOAT = itk::ImageIOBase::FLOAT, DOUBLE = itk::ImageIOBase::DOUBLE};
    
    switch( NormInternalDimension ){
            
        case 3:
            if(ImageDimension!=2) itkGenericExceptionMacro("Unsupported dimension");
            switch (ComponentType) {
                case FLOAT:  return Execute<Riemannian2DNorm<float>  >(argc, argv);
                case DOUBLE: return Execute<Riemannian2DNorm<double> >(argc, argv);
                default: return EXIT_FAILURE;
            }
            
            
        case 4:
            if(ImageDimension!=3) itkGenericExceptionMacro("Unsupported dimension");
            switch (ComponentType) {
                case FLOAT:  return Execute<ExtendedNorm<Riemannian2DNorm<float>  > >(argc, argv);
                case DOUBLE: return Execute<ExtendedNorm<Riemannian2DNorm<double> > >(argc, argv);
                default: return EXIT_FAILURE;
            }
            
        case 5:
            if(ImageDimension!=2) itkGenericExceptionMacro("Unsupported dimension");
            switch (ComponentType) {
                case FLOAT:  return Execute<Finsler2DNorm<float>  >(argc, argv);
                case DOUBLE: return Execute<Finsler2DNorm<double> >(argc, argv);
                default: return EXIT_FAILURE;
            }

        case 6:
            if(ImageDimension!=3) itkGenericExceptionMacro("Unsupported dimension");
            if(is_a_SSRT){
                cout << "Riemannian3DNorm detected" << endl;
                switch (ComponentType) {
                    case FLOAT:  return Execute<Riemannian3DNorm<float>  >(argc, argv);
                    case DOUBLE: return Execute<Riemannian3DNorm<double> >(argc, argv);
                    default: return EXIT_FAILURE;
                }
            } else {
                switch (ComponentType) {
                    case FLOAT:  return Execute<ExtendedNorm<Finsler2DNorm<float>  > >(argc, argv);
                    case DOUBLE: return Execute<ExtendedNorm<Finsler2DNorm<double> > >(argc, argv);
                    default: return EXIT_FAILURE;
                }
            }
            
        case 7:
            if(ImageDimension!=4) itkGenericExceptionMacro("Unsupported dimension");
            switch (ComponentType) {
                case FLOAT:  return Execute<ExtendedNorm<Riemannian3DNorm<float>  > >(argc, argv);
                case DOUBLE: return Execute<ExtendedNorm<Riemannian3DNorm<double> > >(argc, argv);
                default: return EXIT_FAILURE;
            }
            
        default: return EXIT_FAILURE;
    } // switch NormInternalDimension
    
//    return EXIT_SUCCESS;
}

namespace CommandLineMain_Subroutines {
    
    void Usage()
    {
        using std::cerr;
        using std::endl;
        
        cerr << "First argument : metric image filename (read)" << endl;
        cerr << "Second argument : geodesic path endpoints filename (read)" << endl;
        cerr << "third argument : geodesic path filename (write)" << endl;
        cerr << "Fourth argument (optional) : ouput image filename (write)" << endl;

        cerr << endl;
        
        /*TODO : several seeds and tips*/
    }
    
    template<typename NormType>
    int Execute(int argc, char* argv[])
    {
        using std::cerr;
        using std::endl;
        
        const unsigned int Dimension = NormType::Dimension;
        typedef typename NormType::ValueType ValueType;
        typedef itk::Image<NormType,Dimension> NormImageType;
        typedef typename NormImageType::Pointer NormImagePointerType;
        typedef itk::Index<Dimension> IndexType;
        typedef itk::Point<ValueType,Dimension> PointType;
        
        // Import the image
        
        int argumentOffset = 1;
        const char * ImageFileName  = argv[argumentOffset++];
        const char * EndPointsFileName = argv[argumentOffset++];
        const char * GeodesicFileName = argv[argumentOffset++];
        
        typedef itk::ImageFileReader<NormImageType> NormImageReaderType;
        typename NormImageReaderType::Pointer reader = NormImageReaderType::New();
        reader->SetFileName(ImageFileName);
        reader->Update();
        
        NormImagePointerType NormImage = reader->GetOutput();
        const auto RequestedRegion = NormImage->GetRequestedRegion();
        cerr << "RequestedRegion : " << RequestedRegion << endl;
        
        
        std::vector< PointType > PathEndPoints;
        if( ! (std::string(EndPointsFileName) == "0") ){ // Try to import specified values
            
            const char * EndPointsFileName = argv[argumentOffset++];
            try { ImportPointListFromFile(EndPointsFileName, PathEndPoints); }
            catch(itk::ExceptionObject &e){
                cerr << "Invalid endpoints (in) file " << EndPointsFileName << " \n " << e << endl;
                return EXIT_FAILURE;
            }
        } else { // default : opposite corners
            IndexType Center;
            for(int i=0; i<(int)Dimension; ++i)
                Center[i] = (RequestedRegion.GetIndex()[i]+RequestedRegion.GetUpperIndex()[i])/2;
            
            //PathEndIndices.push_back(RequestedRegion.GetIndex());
            PathEndPoints.resize(2);
            NormImage->TransformIndexToPhysicalPoint(Center,PathEndPoints[0]);
            NormImage->TransformIndexToPhysicalPoint(RequestedRegion.GetUpperIndex(),PathEndPoints[1]);
        }
        
        /*
        for(auto it = PathEndIndices.begin(); it!=PathEndIndices.end(); ++it)
            if( !RequestedRegion.IsInside(*it))
                itkGenericExceptionMacro("Point " << *it << " is outside requested region " << RequestedRegion);*/
        
        if(PathEndPoints.size() < 2)
            itkGenericExceptionMacro("Not enough endpoints");
        
        cout << "Path endpoints " << PathEndPoints[0] << ", " << PathEndPoints[1] << endl;
        
        typedef itk::AnisotropicFastMarchingImageFilter<NormType> AFMType;
        
        auto AFM = AFMType::New();
        AFM->SetInput(NormImage);
        AFM->SetGenerateUpwindGradient(true);
        
        auto seeds = AFMType::NodeContainer::New();
        typename AFMType::NodeType node;
        
        node.SetValue( 0. );
        if(! NormImage->TransformPhysicalPointToIndex(PathEndPoints[0], node.GetIndex()) )
            itkGenericExceptionMacro("Start point not in image.");
        
        
        seeds->Initialize();
        seeds->InsertElement( 0, node );
        
        AFM->SetTrialPoints(  seeds  );
        
        // Running the fast marching algorithm
        
        cout << "Spacing : " << NormImage->GetSpacing() << endl;
        cout << "Running the Fast Marching algorithm" << endl;
        AFM->Update();
        
        if(argumentOffset < argc){
            const char * DistanceImageFilename = argv[argumentOffset++];
            
            if(GetExtension(DistanceImageFilename)=="mathematica"){
                auto Image = AFM->GetOutput();
                itk::ImageRegionConstIterator<itk::Image<ValueType,Dimension> > it(Image,Image->GetRequestedRegion());
                it.GoToBegin();
                
                std::ofstream f; f.open(DistanceImageFilename);
                ExportTensorToTXT(f,it,Image->GetRequestedRegion().GetSize());
            } else {
            
                auto writer = itk::ImageFileWriter<itk::Image<ValueType,Dimension> >::New();
                
                writer->SetFileName(DistanceImageFilename);
                writer->SetInput(AFM->GetOutput());
                try { writer->Update(); }
                catch(itk::ExceptionObject &e){
                    cerr << e << endl;
                }
            }
        }
        
        // Extracting the geodesic
        std::vector<PointType> geodesic;
        geodesic.push_back(PathEndPoints[1]);
        
        cout << "Extracting the geodesic" << endl;
        AFM->Geodesic(geodesic);
        cout << "done" << endl;
        // Write the path...
        // Possibilities : SWC (up to 3D), meta object (painful, uselessly general), 1D image of points (chosen here).
        
        ExportPointListToFile<PointType>(GeodesicFileName,geodesic);
        
        return EXIT_SUCCESS;
    }
        
    // ************************* Import Point list routine **************************
    
    template<typename PointType>
    int ImportPointListFromFile(std::string filename,
                                std::vector<PointType> & Points)
    {
        const unsigned int Dimension = PointType::PointDimension;
        const std::string extension = GetExtension(filename);
        
        if(extension=="meta")
        {
            auto spatialObjectPointArrayReader = itk::SpatialObjectReader<Dimension>::New();
            spatialObjectPointArrayReader->SetFileName( filename );
            try
            {
                spatialObjectPointArrayReader->Update();
            }
            catch( itk::ExceptionObject & excp )
            {
                std::cerr << "Can not read the point list file: ";
                std::cerr << filename << std::endl;
                std::cerr << excp << std::endl;
                return EXIT_FAILURE;
            }
            
            auto PointsSpatialObject = static_cast
            <itk::LandmarkSpatialObject<Dimension>*>(spatialObjectPointArrayReader->GetScene()->GetObjects()->front().GetPointer());
            auto spatialObjectPointArray = PointsSpatialObject->GetPoints();
            for(auto it = spatialObjectPointArray.begin();
                it != spatialObjectPointArray.end();
                it++)
                Points.push_back(it->GetPosition());
        }
        else if(extension == "txt") // One point per line, coordinates separated by spaces
        {
            std::ifstream f;
            f.open(filename.c_str());
            std::string line;
            while(std::getline(f,line)){
                std::stringstream linestream(line);
                itk::Point<double,Dimension> P;
                for(int i=0; i<(int)Dimension; ++i){
                    if(linestream >> P[i]) continue;
                    itkGenericExceptionMacro("ImportPointListFromFile error ! Requires one point per line, coordinates separated by spaces.");
                }
                Points.push_back(P);
            }
        }
        else //******** assuming an image format ************
        {
            auto reader = itk::ImageFileReader<itk::Image<itk::Vector<double,Dimension>, 1> >::New();
            reader->SetFileName(filename);
            reader->UpdateOutputInformation();
            
            const auto io = reader->GetImageIO();
            if( io->GetPixelType() != itk::ImageIOBase::VECTOR )
                itkGenericExceptionMacro("Invalid format of points in " << filename << " : should be Vector.");
            if( io->GetNumberOfComponents() != Dimension )
                itkGenericExceptionMacro("Invalid dimension of points in " << filename << " : should be " << Dimension << ".");
            if( io->GetComponentType() != itk::ImageIOBase::DOUBLE )
                itkGenericExceptionMacro("Invalid component type of points in " << filename << " : should be double.");
            if( io->GetNumberOfDimensions() != 1 )
                itkGenericExceptionMacro("Invalid dimension of " << filename << " : should be 1.");
            
            reader->Update();
            
            auto EndPointsImage = reader->GetOutput();
            itk::Point<double,Dimension> Origin(0.);
            const auto BufferStart =EndPointsImage->GetBufferPointer();
            const auto BufferEnd = BufferStart + EndPointsImage->GetRequestedRegion().GetSize()[0];
            
            for(auto it = BufferStart; it != BufferEnd; ++it)
                Points.push_back(Origin + *it);
            //    PathEndPoints.assign(BufferStart,BufferEnd);
            
        }
        
        return EXIT_SUCCESS;
    }
    
    
    
    // ************************* Export Point List routine ************************
    
    template<typename PointType>
    int ExportPointListToFile(std::string filename,
                              const std::vector<PointType> & Points)
    {
        const unsigned int Dimension = PointType::PointDimension;
        const std::string extension = GetExtension(filename);
        
        if(extension=="txt")
        {
            std::ofstream f;
            f.open(filename.c_str());
            if(!f) return EXIT_FAILURE;
            for(auto it=Points.begin(); it!=Points.end(); ++it){
                for(int i=0; i<(int)Dimension; ++i)
                    f << it->operator[](i) << " ";
                f << endl;
            }
        }
        else if(extension=="mathematica")
        {
            std::ofstream f;
            f.open(filename.c_str());
            if(!f) return EXIT_FAILURE;
            f << "{";
            for(auto it = Points.begin(); it!=Points.end(); ++it){
                if(it!=Points.begin()) f << ",";
                f << "{";
                for(int i=0; i<(int)Dimension; ++i){
                    if(i!=0) f << ",";
                    f << (*it)[i];
                }
                f << "}";
            }
            f << "}";
        }
        else // ****** Assuming image format ******
        {
            // Note : images of vectors are used due to an itk bug with images of points.
            typedef itk::Image< itk::Vector<double,Dimension> , 1> PathImageType;
            auto PathImage = PathImageType::New();
            itk::Size<1> size; size.Fill(Points.size());
            itk::Index<1> index; index.Fill(0);
            itk::ImageRegion<1> region(index,size);
            PathImage->SetRegions(region);
            PathImage->Allocate();
            
            auto it1 = Points.begin();
            itk::ImageRegionIterator<PathImageType> it2(PathImage,region);
            for(it2.GoToBegin(); it1!=Points.end(); ++it1,++it2)
                it2.Value() = it1->GetVectorFromOrigin();
            
            auto pathWriter = itk::ImageFileWriter<PathImageType>::New();
            pathWriter->SetFileName(filename);
            pathWriter->SetInput(PathImage);
            try{
                pathWriter->Update();
            } catch(itk::ExceptionObject &e) {
                cout << e << endl;
            }
        }
#pragma message("To do : export as meta file")

        return EXIT_SUCCESS;
    }

    
        
        
        
        
}
#endif
