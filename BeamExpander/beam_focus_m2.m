close all; clear;


filename="Beam-around-focus/beam_000.tiff";
[im]=imread(filename);

imagesc(im); colormap(gray(256));

im_data_x=sum(im,1);
im_data_y=sum(im, 2);
figure; hold on;
plot(im_data_x);
plot(im_data_y);
max_x_i=find(im_data_x);

n=size(im,1); m=size(im, 2);
im_redsize=im(ceil(n/3):floor(n*(1-1/3)), ceil(m/3):floor(m*(1-1/3)));


imstack=[];
for i=1:31
    n=69+i;
    n_str=num2str(n,'%03g');
    filename="Beam-around-focus\beam_"+n_str+".tiff"
    [im, map]=imread(filename_text);
    imstack(:,:,i)=im;
    
    imagesc(im)
end
