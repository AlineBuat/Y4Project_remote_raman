%% ----- Analyse a stack of spot images in a directory -------------
%% Find their standard deviation and plot std against distance (given in a list)
%% intermediary steps : crop the image around the spot to reduce unecesary computations

% Aline Buat



close all; clear;

lambda = 488E-9; % wavelength (in metre): lambda= 488 nm
pixel_size = 3.45E-6; % pixel size (in metre) 

% path for all image files
% dirname = "test_beam_spots\";
% extension_type = ".tif";
% dist = 0:10:110;

% dirname = "test_beam_spots2\";
% extension_type = ".tif";
% dist = 0:10:110;

% dirname = "20220211_KG\";
% extension_type = ".tif";
% dist = 0:5:110; % in mm


% dirname = "20220221_AB\";
% extension_type = ".tiff";
% dist = 0:5:110; % in mm

% dirname = "20220221_AB_2\";
% extension_type = ".tiff";
% dist = [0:5:110, 80:5:110]; % in mm

% dirname = "20220221_AB_3\";
% extension_type = ".tiff";
% dist = 0:5:110; % in mm

%dirname = "Beam-around-focus\";
%extension_type = ".tiff";
%dist=[0,10,50,60,70:80,80:90,90:100,110];


filenames = dirname+"*"+extension_type;
files = dir(filenames);

% tiling of image
a=length(files); m=floor(sqrt(a)); n=ceil(sqrt(a));
if m*n<a
    n=n+1;
end


i=0;


th_p=1/8;
p_padding=1.5;
w=700;

fitstack={};
moment_stack={};
for file=files'
    i=i+1;
    filename=dirname+file.name;
    im=imread(filename);

    % transform to greyscale image if it is not
    if size(im,3) >1
        im = (im(:,:,1)+im(:,:,2)+im(:,:,3))/3; % make image greyscale
    end

    [im_redsize,row,col]=crop_im_around_spot(im,th_p,p_padding);
    row=round(row); col=round(col);
    
    fit_2d=fit_2d_gaussian_image(im_redsize);
    fitstack{i}=fit_2d;
    
    moment_stack{i}=image_moments(im_redsize);

%     lb_x=max([1,row-w/2]); lb_y=max([1,col-w/2]);
%     hb_x=min(size(im,1),row+w/2); hb_y=min(size(im,2),col+w/2);
%     im_redsize=im(lb_x:hb_x , lb_y:hb_y);

%     figure;
%     subplot(1,2,1);imagesc(im);
%     subplot(1,2,2);imagesc(im_redsize);
    
    imstack{i}=im_redsize;
end

%% Show the cropped image stack
figure; montage(imstack,'DisplayRange',[0, 110])
title('Laser Spots montage');

%%
f_multiplot=figure;
for i=1:length(imstack)
    subplot(m,n,i);
    imagesc(imstack{i});
    axis image;
    title(files(i).name,'Interpreter','none');
end

%% analyse beam width changes (std)
std_1=zeros(length(fitstack),1); % stores fits from gauss
std_2=zeros(length(fitstack),1); % stores fits from moment calculations

% run over the entire fit stack
for i=1:length(fitstack)
    f=fitstack{i};
    std_1(i)=mean(abs([f.c1,f.c2]));

    moments=moment_stack{i};
    s_x=sqrt(moments(3,1));
    s_y=sqrt(moments(1,3));
    m_x = moments(1,2);
    m_y = moments(2,1);
    std_2(i)=mean([s_x,s_y]);

    fprintf("File: %s: std=%d, std2=%d\n",files(i).name , std_1(i),std_2(i) );
end


figure;
hold on
plot(dist,std_1,'.','LineWidth',3);
plot(dist,std_2,'.r','LineWidth',3);
legend('std-Gauss fit','std-im moments');
xlabel('Distance');
ylabel('std of spot');
%ylim([0, inf]);
title('Standard Deviation of beam spot vs Distance');


%% sort file by proper name (number)
M=[];
for i=1:length(files)
    M=[M ; string(files(i).name)];
end
M


num=round(linspace(0,250,12));
num_str=string(num);

%% plot std circle
figure;
%colormap("gray");
for i=1:length(imstack)

    subplot(m,n,i);
    
    imagesc(imstack{i}); % plot image
    
    hold on
    axis on

    m_y = moment_stack{i}(1,2);
    m_x = moment_stack{i}(2,1);
    s_x = std_1(i);
    s_y = std_1(i);

    % plot mean and std dev around spot
    plot(m_y,m_x,'+r');
    hold on;
    theta=linspace(0.01,2*pi,200);
    x_t=m_x+s_x.*cos(theta);
    y_t=m_y+s_y.*sin(theta);
    plot(y_t,x_t, 'r');

    title("z="+dist(i)+"mm");
end

% OBSERVATION: fitting to a Gaussian does not provide accurate results

%% Final step: fit measured data points to:
% W^2(z) = W0^2 + M^4*(lambda/pi/W0)^2*(z-z0)^2
% where W(z) = D4sigma/2 

% here: z is in mm (described by the vector "dist"

% defined at the top:
% pixel_size: 3.45 um
% wavelength: lambda=488 nm

%D4s = 4*std_1*pixel_size;   % pixel size is in metre/pixel, std in pixels

D4s = 4*std_2*pixel_size;   % pixel size is in metre/pixel, std in pixels

z = dist*1E-3;  %convert mm to m and use as variable

W = D4s/2;  % (m)

[W0, ind] = min(W);
z0 = z(ind); % dist is in mm

figure;
subplot(1,2,1);
hold on;
plot(dist, D4s,'.', 'DisplayName', 'D4\sigma [m]' );
title("D4\sigma from beam spots");
xlabel("z [mm]");
ylabel("D4\sigma [m]");

subplot(1,2,2);
plot(dist, W*10^3,'.', 'DisplayName', 'w [m]');
title("\omega(z) from beam spots");
xlabel("z [mm]");
ylabel("w [mm]");
ylim([0,inf]);




W2 = W.^2; % variable to fit is W^2(z)
W02 = W0^2      % (m^2)

%fitEq = 'a + b*(x - c)^2';
startPoints = [W02 , 1, z0];
fitEq = 'poly2';
fit_w2 = fit(z', W2, fitEq)


b = fit_w2.p1;
z0 = -fit_w2.p2/2/b
W02 = fit_w2.p3 - b*z0^2
W0 = sqrt(W02)

M4 = b*(pi*W0/lambda)^2;
M2 = sqrt(M4)   % M^2 MUST be >1

theta_divergence = M2*lambda/pi/W0  % beam divergence angle
zR = pi*W02/lambda/M2    % Rayleigh Length


figure;
plot(z, W2,'.','LineWidth',3);
hold on
plot(fit_w2); 

xlabel('z [m]');
ylabel('w^2(z) [m]');
title_text = "Fitting for M^2 factor measurement: M^2=" + num2str(M2);
title(title_text);




fprintf("Beam waist radius w0 = %5.5e m\n", W0);
fprintf("Beam Quality factor M^2 = %g\n", M2);
fprintf("Beam Divergence angle = %.4g rad\n", theta_divergence);
fprintf("Rayleigh length zR = %5.5e m\n", zR);




