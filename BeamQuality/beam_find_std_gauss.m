%% Find standard deviation of a beam image

close all
clear;

%filename = "Beam-around-focus\beam_010.tiff";
%filename = "test_beam_spots\test_3.tif";
filename = "test_beam_spots2\test2_0.tif";


im=imread(filename);

% test if greyscale image (pixel intensity 
if size(im,3) >1
    im = (im(:,:,1)+im(:,:,2)+im(:,:,3))/3; % make image greyscale
end



figure;
subplot(1,2,1);
image(im); %colormap(gray(256));

subplot(1,2,2);
imagesc(im);

% elim background
M = size(im,1); N = size(im,2);
bckgnd = mean2(im(M-round(M/5):end, N-round(N/5):end ));
figure;
imagesc(im); title('bckgnd removed');

th_p=1/10;
p_padding=1.5;
[im_redsize,row,col]=crop_im_around_spot(im,th_p, p_padding);


figure;
imagesc(im_redsize); colormap(gray(256));

%% fit reduced size image
% surface fit to gaussian (ignores values higher than 255)

fit_2d=fit_2d_gaussian_image(im_redsize);
disp(fit_2d);

mx=fit_2d.b1; my=fit_2d.b2;
sx=abs(fit_2d.c1); sy=abs(fit_2d.c2);

% plot the fit
figure;  
plot(fit_2d); hold on;
xlabel('x'); ylabel('y');zlabel('Image Intensity');
x=1:size(im_redsize,1);
y=1:size(im_redsize,2);
[X,Y]=meshgrid(x,y);  
plot3(X,Y,im_redsize-mean2(im)); hold on;





%% display centre?
figure; 
imagesc(im_redsize); colormap(gray(256));
hold on; axis on;

[M,N]=size(im_redsize);


theta=linspace(0.01,2*pi,200);
% equation for an oval: x^2/a^2+y^2/b^2=1

% x_t=mx+sqrt(sx^2+sy^2).*cos(theta);
% y_t=my+sqrt(sx^2+sy^2).*sin(theta);
x_t=mx+sx.*cos(theta);
y_t=my+sy.*sin(theta);

plot(my,mx,'b+');
plot(y_t,x_t,'r');



