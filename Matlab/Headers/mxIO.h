//
//  mxIO.h
//  MatlabITK
//
//  Created by Jean-Marie Mirebeau on 17/10/13.
//
//

#ifndef MatlabITK_mxIO_h
#define MatlabITK_mxIO_h

#include "itkImage.h"
#include "mex.h"
#include <string>
#include <vector>

struct mxIO {
    mxIO(const mxArray *, mxArray **);
    mxIO(const mxIO &) = delete;
    ~mxIO();
    // mxStruct Input
    bool HasField(std::string) const;
    
    template<typename FieldType>
    FieldType
    GetObject(std::string) const;
    
    template<typename FieldType>
    FieldType
    GetObject(std::string,FieldType) const;

    template<typename FieldType>
    std::vector<FieldType>
    GetVector(std::string) const;
    
    template<typename ImageType>
    typename ImageType::Pointer
    GetImage(std::string) const;
    
    const char *
    GetStringField(std::string) const;
    
    // mxStruct Output
    template<typename FieldType>
    mxArray *
    mxObject(FieldType) const;

    template<typename FieldType>
    mxArray *
    mxVector(const std::vector<FieldType> &) const;
    
    template<typename ImageType>
    mxArray *
    mxImage(typename ImageType::Pointer) const;

    void SetField(std::string, mxArray *);
    
    std::string dataName;
    virtual std::string GetNameOfClass() const;
    bool verbose = true;
    bool transposeFirstTwoCoordinates = true;

protected:
    const mxArray * mxInput;
    mxArray ** mxOutput;
    typedef std::pair<std::string,mxArray*> OutputType;
    std::vector<OutputType> output;
//    template<typename FieldType>
//    bool IsScalar() const;
    
    template<typename ImageType>
    typename ImageType::Pointer
    GetImageWOTranspose(std::string) const;

    template<typename ImageType>
    mxArray *
    mxImageWOTranspose(typename ImageType::Pointer) const;
    
    template<typename ImageType>
    typename ImageType::Pointer
    PermuteCoordinates(typename ImageType::Pointer) const;
};

#include "mxIO.hxx"

#endif
