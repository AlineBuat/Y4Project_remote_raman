function fit_2d= fit_2d_gaussian_image(im)
% Fit a 2D Gaussian surface to an image input im

im=double(im);
x=1:size(im,1);
y=1:size(im,2);


% define fit parameters:
fit_eq='a*exp(- ((x-b1)/c1)^2/2 -((y-b2)/c2)^2/2)'; % fit equation: 2d Gaussian

[xo,yo,zo] = prepareSurfaceData(x,y,im);

% start points for regression
startPoints=[max(im,[],'all'),length(x)/2,length(y)/2,length(x)/30,length(y)/30];

% exclude saturated region
exclude_zone= zo>=255;
%exclude_zone=[];
fit_2d=fit([xo,yo],zo,fit_eq,'Start',startPoints,'Exclude',exclude_zone);

end