%mex -output AnisotropicDiffusion  ../Interface/MexInterface.cpp  ../Release/libMatlabAnisotropicDiffusion.dylib
mex CXXFLAGS="-std=c++11" CXXLIBS="\$CXXLIBS -lc++" -output AnisotropicDiffusion  ../Interface/MexInterface.cpp  ../Release/libMatlabAnisotropicDiffusion.dylib

% -v for verbose