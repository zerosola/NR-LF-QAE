function filter = Newthree1_gaborconvolve(height, width, depth, nscale, minWaveLength, ...
                                  sigmaOnf, thetaSigma, angtheta_set, angphi_set, norient)

if mod(width, 2)
    xrange = (-(width - 1) / 2:(width - 1) / 2) / (width - 1);
else
    xrange = (-width / 2:(width / 2 - 1)) / width;	
end

if mod(height, 2)
    yrange = (-(height - 1) / 2:(height - 1) / 2) / (height - 1);
else
    yrange = (-height / 2:(height / 2 - 1)) / height;	
end

if mod(depth, 2)
    zrange = (-(depth - 1) / 2:(depth - 1) / 2) / (depth - 1);
else
    zrange = (-depth / 2:(depth / 2 - 1)) / depth;	
end

[x,y,z] = meshgrid(xrange, yrange, zrange);
x(x==0)=eps;
y(y==0)=eps;
z(z==0)=eps;

radius = sqrt(x.^2 + y.^2 + z.^2);       
theta = atan2(y,x);
phi = asin(z./radius);

radius = ifftshift(radius);     
theta  = ifftshift(theta);
phi  = ifftshift(phi);

for s = 1:nscale                                          % For each scale.
    
    for o = 1:norient
        
        angphi =  angphi_set(s);
        angtheta = angtheta_set(o);

        tk =[cos(angphi)*cos(angtheta), cos(angphi)*sin(angtheta), sin(angphi)];
        dthetaphi =acos(cos(phi).*cos(theta)*tk(1)+cos(phi).*sin(theta)*tk(2)+sin(phi)*tk(3));
        Angread = exp( -(dthetaphi.^2) / (2 * sigmaOnf^2) );      % Calculate the angular filter component.
        
        
        fo = (2/minWaveLength)*3.^0;                         % Centre frequency of filter
        logGabor = exp( -(log(radius/fo).^2) / (2 * log(thetaSigma)^2) );
        
        filter{s,o} =logGabor.* Angread;  
        
        
    end 
    
end



