close all;
clear;

% to do: substract background!

n=10;
%n_str=num2str(n,'%03g');
%filename="Beam-around-focus\beam_"+n_str+".tiff";
%filename = "test_beam_spots\test_2.tif";
filename = "test_beam_spots2\test2_01.tif";


im=imread(filename);
im=double(im);

% transform to greyscale if it is not
if size(im,3) >1
    tic
    im = (im(:,:,1)+im(:,:,2)+im(:,:,3))/3;
    %im = mean(im,3); % make image greyscale
    toc
end

figure;
imagesc(im); colorbar;
im=double(im);

M=size(im,1);
N=size(im,2);


% eliminate background noise by analysing corner of image
%bckgnd=mean(im(1:round(M/5), 1:round(N/5) ),'all');


%im=im-bckgnd; %remove background
%im(im<0)=0;

%% resize image
th_p=1/8;
p_padding=1.5;

[im, ~,~]=crop_im_around_spot(im,th_p,p_padding);
figure; imagesc(im);

%%
M00=sum(im,'all');
fprintf("\nfind centroids using function\n");
moments=image_moments(im)
m_x=moments(2,1);
m_y=moments(1,2);
s_x=sqrt(moments(3,1));
s_y=sqrt(moments(1,3));


fprintf("centroid x=%g;\t centroid y=%g; \n", m_x,m_y);
fprintf("std dev_x=%g;\t std dev_y=%g\n",s_x,s_y);


%% Display
figure;
image(im);
hold on
axis on

% plot mean and std dev
plot(m_y,m_x,'+r','LineWidth',1);
hold on
theta=linspace(0.01,2*pi,200);
x_t=m_x+s_x.*cos(theta);
y_t=m_y+s_y.*sin(theta);
plot(y_t,x_t, 'r','LineWidth',1);




