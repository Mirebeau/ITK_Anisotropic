% Linear diffusion Demo
if false
    clear('input');
    n=11;

    input.ImageDimension=2;
    input.PixelDimension=0; %Scalar image
    
    input.Image = 0*ones(n,n);
    input.Image((n+1)/2,(n+1)/2)=1;
    
    input.DiffusionTensorField = ones(3,n,n);
    input.DiffusionTensorField(2,:,:)=0;
    
    output = AnisotropicDiffusion(input);
    imshow(output.SmoothedImage/max(output.SmoothedImage(:)));
end

%Non-linear anisotropic diffusion demo
if true
    clear('input');
    input.Image = imread('EdgeEnhancingDiffusion2D_TestImage.png');
    input.Image = double(input.Image);
%    input.Image = input.Image/max(input.Image(:));
    
    input.DiffusionTime = 8;
    input.EdgeEnhancement = 1;
    input.NoiseScale = 0.5;
    input.Alpha=0.01;
    input.Lambda = 0.02;
    input.FeatureScale=4;
    output = AnisotropicDiffusion(input);
    imshow(output.SmoothedImage/max(output.SmoothedImage(:))); 
end

