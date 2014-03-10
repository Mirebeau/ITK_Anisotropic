//
//  mxIO.hxx
//  MatlabITK
//
//  Created by Jean-Marie Mirebeau on 17/10/13.
//
//

#ifndef MatlabITK_mxIO_hxx
#define MatlabITK_mxIO_hxx

#include "MexMessageWrapper.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPermuteAxesImageFilter.h"
#include <numeric>

// %%%%%%%%%%%%%% As data source %%%%%%%%%%%%%%%%%%
bool
mxIO::HasField(std::string fieldName) const {
    return mxGetFieldNumber(mxInput, fieldName.c_str()) >= 0;
//    return mxGetField(mxInput,0,fieldName.c_str());
}

std::string mxIO::GetNameOfClass() const {
    return "mxIO" + ((dataName.size()>0) ? "-"+dataName : "");
}

mxIO::mxIO(const mxArray * mxInput_, mxArray ** mxOutput_):mxInput(mxInput_), mxOutput(mxOutput_){
    if(!mxIsStruct(mxInput))
        itkExceptionMacro("Input not a struct");
    
    
    for(size_t i = 0; i<mxGetNumberOfDimensions(mxInput); ++i)
        if(mxGetDimensions(mxInput)[i] != 1)
            itkExceptionMacro("Input structure is not flat");
}

mxIO::~mxIO(){
    // copy all field names into a strong and static structure...
    const int totalLen = std::accumulate(output.begin(), output.end(), 0,
                                   [](int i, const OutputType & o){return i+o.first.length();} )
    +output.size();
    
    char * allStrings = new char [totalLen];
    std::fill(allStrings, allStrings+totalLen, 0);
    
    std::vector<char*> stringStarts;
    stringStarts.push_back(allStrings);
    
    for(const OutputType & o : output){
        strcpy(stringStarts.back(), o.first.c_str());
        stringStarts.push_back(stringStarts.back()+o.first.length()+1);
    }
    stringStarts.back()=nullptr;
    
    const mwSize dims[]={1,1};
    *mxOutput = mxCreateStructArray(2,dims,output.size(),const_cast<const char **>(&stringStarts[0]));
    for(int i=0; i<output.size(); ++i)
        mxSetField(*mxOutput, 0, stringStarts[i], output[i].second);
    
    delete allStrings;
    
    /*
     // For some reason, the old compact version stopped working ... ?
    std::vector<const char *> fieldNames;
    for(auto nameValue:output)
                fieldNames.push_back(nameValue.first.c_str());
    
    if(*mxOutput != NULL) { // Useless it seems
        MexMsg() << "Destroyed previous data\n";
        mxDestroyArray(*mxOutput);
    };
    
    const mwSize dims[]={1,1};
    *mxOutput = mxCreateStructArray(2,dims,output.size(),&fieldNames[0]);
    
    for(auto stringArray:output)
        mxSetField(*mxOutput, 0, stringArray.first.c_str(), stringArray.second);
     */
    
    /* // Testing
     const char * fieldnames[2]={"Hello", "Bye"};
     const mwSize dims[]={1,1};
     *mxOutput = mxCreateStructArray(2, dims, 2, fieldnames);
     mxSetField(*mxOutput, 0, fieldnames[0], output[0].second); //mxCreateString("World"));
     mxSetField(*mxOutput, 0, fieldnames[1], output[1].second); //mxCreateDoubleScalar(18.68));
     
     MexMsg() << output[0].first << "," << output[1].first << "\n";
     
     std::string s="helllo";
     */
}

template<typename FieldType>
FieldType
mxIO::GetObject(std::string fieldName) const {
    const auto vector = GetVector<FieldType>(fieldName);
    if(vector.size()!=1)
        itkExceptionMacro("(GetObject) Size of field " << fieldName << " (" << vector.size() << ") is incorrect.");
    return vector[0];
}

template<typename FieldType>
FieldType
mxIO::GetObject(std::string fieldName, FieldType defaultValue) const {
    if(!HasField(fieldName)) {
        if(verbose)
            MexMsg() << "Default value " << defaultValue << " used for " << fieldName << ".\n";
        return defaultValue;
    }
    return GetObject<FieldType>(fieldName);
}

template<typename FieldType>
std::vector<FieldType>
mxIO::GetVector(std::string fieldName) const {
    typedef itk::Image<FieldType,1> FieldImageType;
    auto image = GetImageWOTranspose<FieldImageType>(fieldName);
    
    std::vector<FieldType> vector;
    if(image.IsNull())
        itkExceptionMacro("GetVector " << fieldName << " error : Null image pointer.");

    itk::ImageRegionConstIterator<FieldImageType> it(image,image->GetBufferedRegion());
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        vector.push_back(it.Value());
    return vector;
}

template<typename ImageType>
typename ImageType::Pointer
mxIO::GetImage(std::string fieldName) const {
    return PermuteCoordinates<ImageType>(GetImageWOTranspose<ImageType>(fieldName));
}

const char *
mxIO::GetStringField(std::string fieldName) const {
    
    const mxArray * field = mxGetField(mxInput,0,fieldName.c_str());
    if(field==NULL)
        itkExceptionMacro("GetStringField error : field " << fieldName << " does not exist.\n");
    
    const char * ans = mxArrayToString(field);
    if(ans==NULL)
        itkExceptionMacro("GetStringField error : field " << fieldName << " is not a string.\n");

    return ans;
}

// %%%%%%%%%%%%%%%% As data sink %%%%%%%%%%%%%%%%%

template<typename FieldType>
mxArray *
mxIO::mxObject(FieldType value) const {
    std::vector<FieldType> values;
    values.push_back(value);
    return mxVector(values);
}

template<typename FieldType>
mxArray *
mxIO::mxVector(const std::vector<FieldType> & values) const {
    typedef itk::Image<FieldType,1> ImageType;
    auto image=ImageType::New();
    itk::ImageRegion<1> region;
    region.SetIndex({{0}});
    region.SetSize({{values.size()}});
    image->SetRegions(region);
    image->Allocate();
    itk::ImageRegionIteratorWithIndex<ImageType> it(image, image->GetBufferedRegion());
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        it.Value() = values[it.GetIndex()[0]];
    return mxImageWOTranspose<ImageType>(image);
}

template<typename ImageType>
mxArray *
mxIO::mxImage(typename ImageType::Pointer image) const {
    return mxImageWOTranspose<ImageType>(PermuteCoordinates<ImageType>(image));
}

void
mxIO::SetField(std::string fieldName, mxArray * value){
    output.push_back({fieldName,value});
}
// Protected methods

template<typename ImageType>
typename ImageType::Pointer
mxIO::PermuteCoordinates(typename ImageType::Pointer inputImage) const
{
    const int Dimension = ImageType::ImageDimension;
    if(Dimension<2)
        itkExceptionMacro("PermuteCoordinates error : image dimension < 2.\n" << inputImage);
    
    if(!transposeFirstTwoCoordinates) return inputImage;
    
//    return TransposeFirstTwoImageCoordinates<ImageType>(image);
    auto inputRegion = inputImage->GetRequestedRegion();
    
    auto outputSize = inputRegion.GetSize();
    std::swap(outputSize[0],outputSize[1]);
    
    auto outputIndex = inputRegion.GetIndex();
    std::swap(outputIndex[0],outputIndex[1]);
    
    // Physical offset
    
    auto outputSpacing = inputImage->GetSpacing();
    std::swap(outputSpacing[0],outputSpacing[1]);
    
    auto outputOrigin = inputImage->GetOrigin();
    std::swap(outputOrigin[0],outputOrigin[1]);
    
    itk::ImageRegion<Dimension> outputRegion(outputIndex,outputSize);
    auto outputImage = ImageType::New();
    outputImage->SetRegions(outputRegion);
    outputImage->SetSpacing(outputSpacing);
    outputImage->SetOrigin(outputOrigin);
    outputImage->Allocate();
    
    itk::ImageRegionConstIteratorWithIndex<ImageType> inputIt(inputImage, inputImage->GetLargestPossibleRegion());
    for(inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt){
        itk::Index<Dimension> outputIndex = inputIt.GetIndex();
        std::swap(outputIndex[0],outputIndex[1]);
        
        outputImage->GetPixel(outputIndex) = inputIt.Value();
    }
    return outputImage;
    
    /* Curiously does not work
    typedef itk::PermuteAxesImageFilter<ImageType> FilterType;
    typename FilterType::PermuteOrderArrayType permutation;
    permutation[0]=1;
    permutation[1]=0;
    for(int i=2; i<Dimension; ++i)
        permutation[i]=i;
    
    auto filter = FilterType::New();
    filter->SetOrder(permutation);
    filter->SetInput(image);
    
    filter->Update();
    return filter->GetOutput();
     */
}


template<typename ImageType>
typename ImageType::Pointer
mxIO::GetImageWOTranspose(std::string fieldName) const
{
    std::string errMsg = "GetImageWOTranspose(" + fieldName + ") error : ";
    
    if(!HasField(fieldName))
        itkExceptionMacro("Field " << fieldName << " does not exist.");
    
    mxArray * p = mxGetField(mxInput,0,fieldName.c_str());
    typedef double mxComponentType;
    
    if(!mxIsNumeric(p)) itkExceptionMacro(<< errMsg << "Not a numeric array.");
    
    if(mxIsComplex(p)) itkExceptionMacro(<< errMsg << "Complex arrays not supported.");
    if(mxGetClassID(p) != mxDOUBLE_CLASS) itkExceptionMacro(<< errMsg << "Array must be of double type.");
    
    typedef typename ImageType::PixelType PixelType;
    const bool PixelIsAtomic = (typeid(PixelType)==typeid(mxComponentType));
    
    const mwSize PixelAtoms = sizeof(PixelType)/sizeof(mxComponentType);
    if(sizeof(PixelType)%sizeof(mxComponentType) != 0)
        itkExceptionMacro(<< errMsg << "InvalidPixelType.");
    
    const mwSize * mxSize = mxGetDimensions(p);
    
    if( ! PixelIsAtomic )
        if(mxSize[0] != PixelAtoms)
            itkExceptionMacro(<< errMsg << "First mxNumericArray dimension, " << mxSize[0] << ", does not match pixel size, " << PixelAtoms << ".");
    
    const mwSize FirstDimension = (PixelIsAtomic ? 0 : 1);
    const mwSize Dimension = ImageType::ImageDimension;
    const mwSize ImageFullDimension = FirstDimension+Dimension;
    const mwSize ArrayDimension = mxGetNumberOfDimensions(p);
    
    if(ImageFullDimension > ArrayDimension)
        itkExceptionMacro(<< errMsg << "Dimensions of image and array do not match. "
                          << "Respectively : " << FirstDimension + Dimension << " and "
                          << ArrayDimension);
    if(ImageFullDimension < ArrayDimension){
        for(int i=ImageFullDimension; i<ArrayDimension; ++i)
            if(mxSize[i]!=1)
                itkExceptionMacro(<< errMsg << "Dimensions of image and array do not match. "
                                  << "Respectively : " << FirstDimension + Dimension << " and "
                                  << ArrayDimension << ". Additional dimension(s) are not unit." );
    }
    
    typename ImageType::SizeType Size;
    for(int i=0; i<Dimension; ++i)
        Size[i] = mxSize[FirstDimension+i];
    
    typename ImageType::IndexType Origin;
    Origin.Fill(0);
    
    itk::ImageRegion<Dimension> Region(Origin,Size);
    auto Image = ImageType::New();
    Image->SetRegions(Region);
    Image->Allocate();
    
    mxComponentType * ImageBuffer = reinterpret_cast<mxComponentType *>( Image->GetBufferPointer() );
    const mxComponentType * mxArrayBuffer = mxGetPr(p);
    const mxComponentType * const mxArrayBufferEnd = mxArrayBuffer+mxGetNumberOfElements(p);
    
    for(;mxArrayBuffer!=mxArrayBufferEnd; ++mxArrayBuffer, ++ImageBuffer)
        (*ImageBuffer) = (*mxArrayBuffer);
    
    return Image;
}


template<typename ImageType>
mxArray *
mxIO::mxImageWOTranspose(typename ImageType::Pointer image) const {
    typedef double mxComponentType;
    typedef typename ImageType::PixelType PixelType;
    
    const bool PixelIsAtomic = (typeid(PixelType)==typeid(mxComponentType));
    const mwSize PixelAtoms = sizeof(PixelType)/sizeof(mxComponentType);
    if(sizeof(PixelType)%sizeof(mxComponentType) != 0)
        itkExceptionMacro("itkImageWOTranspose error : InvalidPixelType.\n" << image);
    
    const mwSize FirstDimension = PixelIsAtomic ? 0 : 1;
    const mwSize ImageDimension = ImageType::ImageDimension;
    const mwSize ArrayDimension = FirstDimension + ImageDimension;
    
    std::vector<mwSize> mxSize_;
    mxSize_.resize(ArrayDimension);
    mwSize * mxSize = &mxSize_[0];
    
    if(!PixelIsAtomic) mxSize[0] = PixelAtoms;
    for(mwSize i=0; i<ImageDimension; ++i)
        mxSize[FirstDimension+i] = image->GetBufferedRegion().GetSize()[i];
    
    mxArray * p = mxCreateNumericArray(ArrayDimension,mxSize,mxDOUBLE_CLASS,mxREAL);
    
    mxComponentType * mxArrayBuffer = mxGetPr(p);
    const mxComponentType * const mxArrayBufferEnd = mxArrayBuffer + mxGetNumberOfElements(p);
    const mxComponentType * ImageBuffer = reinterpret_cast<const mxComponentType *>(image->GetBufferPointer());
    
    for(;mxArrayBuffer!=mxArrayBufferEnd; ++mxArrayBuffer, ++ImageBuffer)
        (*mxArrayBuffer) = (*ImageBuffer);
    
    return p;
}
#endif
