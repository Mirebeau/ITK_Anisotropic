test.input_info = false;
test.output_info = false;
test.missing_data = false;
test.Iso2D = false;
test.TipsAbort = false;
test.EuclideanDistanceAbort = false;
test.Aniso3D = false;
test.StopWhenFirstIsAccepted = false;
test.AsymmetricQuadratic2D = true;

if test.input_info 
    AnisotropicFastMarching_EarlyAbort(); 
end

if test.output_info
    AnisotropicFastMarching_EarlyAbort(0); 
end

if test.missing_data
    clear('input');
    input.NormType = 'Riemannian2DNorm'; % Missing data
    output = AnisotropicFastMarching_EarlyAbort(input);
end

if test.Iso2D
    disp('--- TEST --- 2D Isotropic fast marching.');
    clear('input'); %clear('output');
    input.NormType = 'Riemannian2DNorm';

    n = 200; 
    r = linspace(0,1,n);
    [x,y] = meshgrid(r,r); 
    input.Origin=[x(1,1);y(1,1)];
    input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)]; 
    input.TransposeFirstTwoImageCoordinates=1; % Now on by default
    
    Speed = ones(n,n) + (x>0.5);
    input.Metric = zeros([3,n,n]);
    input.Metric(1,:,:) = 1./(Speed.*Speed);
    input.Metric(3,:,:) = input.Metric(1,:,:);
    
    input.Seeds = [0.1,0.9,0.2;0.5,0.2,0.8; 0,0,0]; % x coords, y coords, values
    
    s=10;
    R=linspace(0,1,s);
    [Y,X] = meshgrid(R,R); %Standard coordinates for seeds, Tips, etc
    input.Tips = [X(:)';Y(:)'];
    
    output = AnisotropicFastMarching_EarlyAbort(input);
    imagesc(output.Distance);
    
    for i=1:size(output.Geodesics,1)
        rescaledGeodesic = RescaledCoords(output.Geodesics{i},input.Origin,input.Spacing);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;

    fprintf('\n\n');

    pause;
    clf;
end

if test.TipsAbort
    disp('--- DEMO --- 2D Asymmetric fast marching, aborted when tips are accepted.');
    clear('input'); clear('output');
    input.NormType = 'Finsler2DNorm';
    n=100; s=10;

    r = linspace(-2,2,n);
    [x,y]=meshgrid(r,r); 
    
    input.Origin=[x(1,1);y(1,1)];
    input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
    input.TransposeFirstTwoImageCoordinates=1;

    Speed = ones(5,n,n);
    Speed(2,:,:)=0;

    scals = sqrt(1+x.*x+y.*y).^-1;
    Speed(4,:,:)=y.*scals;
    Speed(5,:,:)=-x.*scals;
    input.Metric = Speed;
    
    input.Seeds = [0;0;0];

    input.Tips = [1.2,-0.5;1.2,-0.2];

    input.StopWhenTipsAreReached = 1;
    
    clf;
    output = AnisotropicFastMarching_EarlyAbort(input);
    
    % To get a better image, set to zero the points not reached by the Fast
    % Marching algorithm (by default DBL_MAX)
    Distance = output.Distance;
    NotReached = Distance == max(Distance(:));
    Distance(NotReached) = 0;
    imagesc(Distance);

    for i=1:size(output.Geodesics,1)
        rescaledGeodesic = RescaledCoords(output.Geodesics{i},input.Origin,input.Spacing);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;
    fprintf('\n\n');
    pause;
end;

if test.EuclideanDistanceAbort
    disp('--- DEMO --- 2D Anisotropic fast marching, aborted at some distance from the origin.');
    clear('input'); clear('output');
    input.NormType='Riemannian2DNorm';
    n=200;
    r = linspace(-1,1,n);
    [x,y] = meshgrid(r,r);
    clear('options');
    input.Origin=[x(1,1);y(1,1)];
    input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
    input.TransposeFirstTwoImageCoordinates = 1;
    input.Seeds = [0;0;0];
    
    Radius = sqrt(x.*x+y.*y) + 1e-8; % avoid division by zero
    Ux = x./Radius;
    Uy = y./Radius;
    
    Theta = atan2(y,x);
    
    % Construct a nice spiral
    k=0.05;
    delta = 0.03;
    Spiral = zeros(size(x));
    for i=0:5
        Thetai = Theta+2*i*pi;
        Spirali = (abs(k*Thetai - Radius)< delta);
        
        Theta(Spirali) = Thetai(Spirali);
        Spiral = (Spiral + Spirali)>0;
    end
    
    % Tangent vector to the spiral
    Vx = k*Ux - k*Theta.*Uy;
    Vy = k*Uy + k*Theta.*Ux;
    RV = sqrt(Vx.*Vx+Vy.*Vy);
    Vx = Vx./RV;
    Vy = Vy./RV;
    
    % Construct the riemannian metric
    
    IdentityMetric = ones(3,n,n);
    IdentityMetric(2,:,:) = 0;
    
    Lambda1 = (3e-2)^2;
    Lambda2 = 1;
    RotatingMetric = zeros(3,n,n);
    RotatingMetric(1,:,:) = Lambda1*Vx.*Vx + Lambda2*Vy.*Vy;
    RotatingMetric(2,:,:) = Lambda1*Vx.*Vy - Lambda2*Vx.*Vy;
    RotatingMetric(3,:,:) = Lambda1*Vy.*Vy + Lambda2*Vx.*Vx;
    
    input.Metric = zeros(3,n,n);
    Spiral = reshape(Spiral,1,n,n);
    for i=1:3
        input.Metric(i,:,:) = IdentityMetric(i,:,:).*(1-Spiral)+RotatingMetric(i,:,:).*Spiral;
    end
    Spiral = reshape(Spiral,n,n);

    input.StopAtEuclideanDistance = 5;
    output = AnisotropicFastMarching_EarlyAbort(input);
    
    for i=1:2
        if i==1 
            Distance = output.Distance;
        else
            Distance = output.EuclideanPathLengths;
        end
        
        NotReached = Distance == max(Distance(:));
        Distance(NotReached) = 0;
        
        imagesc(Distance);
        rescaledGeodesic = RescaledCoords(output.GeodesicFromStoppingPoint,input.Origin,input.Spacing);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
        fprintf('\n\n');
        pause;
    end
end


if test.Aniso3D
    disp('--- DEMO --- 3D Anisotropic fast marching.');
    clear('input'); clear('output');
    input.NormType = 'Riemannian3DNorm';
    n=50;
    r=linspace(-1.1,1.1,n);
    rz = linspace(0,3,2.5*n);
    [x,y,z] = meshgrid(r,r,rz);
    Size = size(x);
    input.Origin = [x(1,1,1);y(1,1,1);z(1,1,1)];
    input.Spacing = [x(2,2,2)-x(1,1,1); y(2,2,2)-y(1,1,1); z(2,2,2)-z(1,1,1)];
    input.TransposeFirstTwoImageCoordinates=1;
    
    radius = sqrt(x.*x+y.*y)+1e-8;
    
    k=0.2;
    TangentVector = zeros([3,Size]);
    TangentVector(1,:,:,:) = reshape(-y./radius,[1,Size]);
    TangentVector(2,:,:,:) = reshape( x./radius,[1,Size]);
    TangentVector(3,:,:,:) =  k;
    
    RV = sqrt(1+k^2);
    TangentVector = TangentVector/RV;
    
    delta = 0.2;
    Spiral = (abs(radius-1)<delta) & (abs(cos(z/k)-x)<delta) & (abs(sin(z/k)-y)<delta);
    
    % Metric coefficients : 11, 12, 13, 22, 23 ,33
    IdentityMetric = ones([6,Size]);
    IdentityMetric(2,:,:,:)=0;
    IdentityMetric(3,:,:,:)=0;
    IdentityMetric(5,:,:,:)=0;
    
    RotatingMetric = zeros([6,Size]);
    RotatingMetric(1,:,:,:) = TangentVector(1,:,:,:).*TangentVector(1,:,:,:);
    RotatingMetric(2,:,:,:) = TangentVector(1,:,:,:).*TangentVector(2,:,:,:);
    RotatingMetric(3,:,:,:) = TangentVector(1,:,:,:).*TangentVector(3,:,:,:);
    RotatingMetric(4,:,:,:) = TangentVector(2,:,:,:).*TangentVector(2,:,:,:);
    RotatingMetric(5,:,:,:) = TangentVector(2,:,:,:).*TangentVector(3,:,:,:);
    RotatingMetric(6,:,:,:) = TangentVector(3,:,:,:).*TangentVector(3,:,:,:);

    lambda = (0.02)^2;
    input.Metric = zeros([6,Size]);
    Spiral = reshape(Spiral,[1,Size]);
    for i=1:6
        input.Metric(i,:,:,:) = IdentityMetric(i,:,:,:)-(1-lambda)*RotatingMetric(i,:,:,:).*Spiral;
    end
    Spiral = reshape(Spiral,Size);
    
    input.Seeds = [0;0;0;0];
    input.Tips = [0;0;3];
    
    output = AnisotropicFastMarching_EarlyAbort(input);
    clf;
    isosurface(x,y,z,output.Distance,2);
    isosurface(x,y,z,output.Distance,1);
    axis equal;
    geodesic = output.Geodesics{1};
    line(geodesic(1,:),geodesic(2,:),geodesic(3,:));
    pause;
end

if(test.StopWhenFirstIsAccepted)
    disp('--- DEMO --- Stop when first is accepted.');
    clear('input'); clear('output');
    input.NormType='Riemannian2DNorm';
    n=100;
    r = linspace(-1,1,n);
    [x,y] = meshgrid(r,r);
    clear('options');
    input.Origin=[x(1,1);y(1,1)];
    input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
    input.TransposeFirstTwoImageCoordinates = 1;
    
    input.Metric = ones(3,n,n); % Identity metric
    input.Metric(2,:,:) = 0;
    
    
    input.Seeds = [0;0;0];
    input.Tips = [0.1,-0.3,0.7,-0.9;0.3,0.4,-0.5,-0.9];
    input.StopWhenFirstIsAccepted = [0.6,-0.5,0.9; 0.7,0.8,-0.5];

    output = AnisotropicFastMarching_EarlyAbort(input);
    
    
    NotReached = output.Distance == max(output.Distance(:));
    output.Distance(NotReached) = 0;
    
    imagesc(output.Distance);
    for i=1:size(output.Geodesics,1)
        rescaledGeodesic = RescaledCoords(output.Geodesics{i},input.Origin,input.Spacing);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;
    
    rescaledGeodesic = RescaledCoords(output.GeodesicFromStoppingPoint,input.Origin,input.Spacing);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));

    fprintf('\n\n');
    pause;
end

if(test.AsymmetricQuadratic2D)
    % This test demonstrates fast marching w.r.t. asymmetric quadratic norms, which take the form
    %   F_x(v) = sqrt( v.M(x).v + max(0,omega(x).v)^2)
    
    % An interesting special case is the so called half-disk model, defined by 
    % M(x) = Id/speed(x)^2
    % omega(x) = (-1/eps) omega0(x)/speed(x)
    % where 
    % - speed is the speed function. 
    % - eps is a relaxation parameter, e.g. eps = 1./10
    % - omega0 is a field of unit vectors.
    % In the limiting model, as eps->0, the length of a path is the ordinary euclidean length, weighted by 1/speed(x).
    % However, paths which tangent has, at any point, a negative scalar product with omegav0(x), are rejected.
    
    disp('---- DEMO ---- AsymmetricQuadratic2DNorm (Half disk model)');
        clear('input'); clear('output');
    input.NormType='AsymmetricQuadratic2DNorm';
    n=100;
    r = linspace(-1,1,n);
    [x,y] = meshgrid(r,r);
    clear('options');
    input.Origin=[x(1,1);y(1,1)];
    input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
    input.TransposeFirstTwoImageCoordinates = 1;
    
    input.Metric = zeros(5,n,n); 
    input.Metric(1,:,:) = 1; % Tensor field M(x). Here the identity tensor.
    input.Metric(3,:,:) = 1;
    
    input.Metric(4,:,:) = 5; % Vector fied omega(x). Here a constant vector field.
    input.Metric(5,:,:) = 5;
    
    input.Seeds = [0;0;0];
    input.Tips = [0.1,-0.3,0.7,-0.9;0.3,0.4,-0.5,-0.9];
    %input.StopWhenFirstIsAccepted = [0.6,-0.5,0.9; 0.7,0.8,-0.5];

    output = AnisotropicFastMarching_EarlyAbort(input);
    
    NotReached = output.Distance == max(output.Distance(:));
    output.Distance(NotReached) = 0;
    
    imagesc(output.Distance);
    for i=1:size(output.Geodesics,1)
        rescaledGeodesic = RescaledCoords(output.Geodesics{i},input.Origin,input.Spacing);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;
    
    rescaledGeodesic = RescaledCoords(output.GeodesicFromStoppingPoint,input.Origin,input.Spacing);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));

    fprintf('\n\n');
    pause;
end