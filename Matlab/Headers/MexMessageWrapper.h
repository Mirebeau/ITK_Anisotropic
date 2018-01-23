//
//  MexMessageWrapper.h
//  MatlabITKFastMarching
//
//  Created by Jean-Marie Mirebeau on 04/04/13.
//
//

#ifndef MatlabITKFastMarching_MexMessageWrapper_h
#define MatlabITKFastMarching_MexMessageWrapper_h


#include <sstream>
#include "mex.h"

template<int MessageGrade=0>
struct MexMessageWrapper {
    enum MexMessageGrade {Normal,Warning,Error};
    
    template<typename DataType>
    MexMessageWrapper & operator << (DataType data){oss << data; return (*this);}
    
    ~MexMessageWrapper(){
        const char * msg = oss.str().c_str();
        switch (MessageGrade) {
            case Warning:
                mexWarnMsgIdAndTxt("ITKFM:Warning","%s",oss.str().c_str());
                break;
                
            case Error:
                mexWarnMsgIdAndTxt("ITKFM:Error","%s",oss.str().c_str());
                break;
                
            default:
                mexPrintf("%s",oss.str().c_str());
                break;
        };
    };
    
    std::ostringstream oss;
};
        
typedef MexMessageWrapper<MexMessageWrapper<>::Warning> MexWarnMsg;
typedef MexMessageWrapper<MexMessageWrapper<>::Normal>  MexMsg;
        
        
        

#endif
