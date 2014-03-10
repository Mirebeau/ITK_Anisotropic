Usage of ITKFM
Image formats .nii or .hdr (or .nrrd) are suggested but not required. 

Example. Once the program is built, type in a terminal :
./ITKFM data/Identity2D_float.nii 0 geodesic.txt distances.nii

----- First argument : metric image (read) -----
This image can be 2D, 3D or 4D. 
The pixeltype of this image can be
 - 2D symmetric second rank tensor
 - 3D symmetric second rank tensor
 - 2D symmetric second rank tensor followed by a 2D vector (Asymmetric norm |u|_M - omega*u)
 - Any of the above followed by a scalar (ND+1 tubular metric)

Since pixel types cannot be arbitrary, an image of vectors of the appropriate dimension should be used in the last two cases.
The pixel component type can be float or double.


----- Second component : geodesic path endpoints filename (read) -----
This file should specify a list of points. The first point will be used as a seed for the fast marching algorithm, and a geodesic path from the second point to the first will be extracted. Further points will be ignored.

Supported file extensions : 
 - txt : one point per line. Coordinates separated by spaces.
 - meta
 - any suitable image extension : the required point list can be given as a one dimensional image of vectors (Note : itk does not support images of points so far.)

 If you enter 0 instead of a filename, the center and the upper corner will be used.


----- Third argument : geodesic path filename (write) -----
Supported file extensions are 
 - txt : same as above.
 - mathematica : path can be imported in Mathematica via ToExpression[Import[filename]]
 - any suitable image extension : see above. 

 ----- Fourth argument (optional) : ouput image filename (write) -----
 Supported file extensions are :
  - any suitable image extension : The pixels of this image are of the same type, float or double, as the pixel components of the metric image.
  - mathematica : see above.