%% Finding beam waist and size

close all;
clear all;
im=imread("Beam-around-focus\beam_070.tiff");
figure; imagesc(im); colormap(gray(256));

% sum image along each direction to have 1D sum signal
im_s_x=sum(im, 1);
im_s_y=sum(im,2);

% smoothing mask for 1d signal
l_mask=30;
mask_d=ones(l_mask,1)/l_mask;
% smooth signal (smooth background) with square window
im_x_red=conv(im_s_x, mask_d,'valid'); % careful: size change
im_y_red=conv(im_s_y, mask_d,'valid');
% could add some padding with average background value to the sides to
% ensure 'same' signal size

% plot original sum and smoothed
figure;
hold on;
plot(im_s_x);
plot(im_x_red);

% minimal value of smoothed = background in non-signal regions
min_val=min(im_x_red);
im_s_x_adj=im_s_x-min_val;

min_val=min(im_y_red);
im_s_y_adj=im_s_y - min_val;

figure;
plot(im_s_x_adj);


% fit gaussian to the sum signal -- find center of spots
% attention: some spots have destructive interference at their centres
% after going through the focal point

x=[1:length(im_s_x_adj)];
y=im_s_x_adj;

f_x=fit(x.',y.','gauss1')

x=[1:length(im_s_y_adj)];
y=im_s_y_adj;

f_y=fit(x.',y,'gauss1')
figure; plot(f_y,x,y);



%%  plot centre of beam on image
figure(1); hold on;
[m_x,I_x]=max(im_s_x);
[m_y,I_y]=max(im_s_y);
%plot(ones(size(im,2)).*I_y)
%plot(ones(size(im,1)).*I_x,1:size(im,1));

mu_x=f_x.b1;
mu_y=f_y.b1;
hold on
plot(ones(size(im,2)).*mu_y);
plot(ones(size(im,1)).*mu_x , 1:size(im,1), '-r');


%% take a slice along the centre


%% Cropping the image

[N,M]=size(im);
p=1/20; % proportion to keep

% index bounds for image to keep, centred around spot centre
N_red_lb=ceil(mu_y-N*p)
N_red_hb=floor(mu_y+N*p)
M_red_lb=ceil(mu_x-M*p)
M_red_hb=floor(mu_x+M*p)


% crop image
im_red=im(N_red_lb:N_red_hb,M_red_lb:M_red_hb);
figure;
imagesc(im_red); colormap(gray(256));

hold on;


%%  fit 2D gaussian to image
% fit parameters: a, b,c,d,e,f
d_im_red=double(im_red);
x=1:size(d_im_red,1);
y=1:size(d_im_red,2);
[xo,yo,zo] = prepareSurfaceData(x,y,d_im_red);

fit_eq='a*exp(- ((x-b1)/c1)^2 -((y-b2)/c2)^2)+d';
%fit_gauss2d_im=fit([xo,yo],zo,'a*(exp(-((x-mux)/sx)^2) + exp(-((y-muy)/sy)^2 )) + b')
fit_gauss2d_im=fit([xo,yo],zo,fit_eq);

% assign fit parameters to variables


figure; % plot fit along x direction
%plot(a*exp(-((x-mux)/sx).^2)+b);
plot(fit_gauss2d_im);
hold on

[X,Y]=meshgrid(x,y);
plot3(X,Y,d_im_red)
