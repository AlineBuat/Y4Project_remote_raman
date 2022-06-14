function [im_redsize,row_c,col_c]=crop_im_around_spot(im, th_p, p_padding)
% Ouputs thresholded version of image
% INPUTS: crop_im_around_spot(im,th_p,p_padding)
% number of pixels is reduced: real resolution (in number of pixels) is
% conserved
% only image size (NxM) is modified

% inputs:
% im -- image
% th_p -- proportion of the original image maximum used as a threshold
% p_padding -- padding around the detected thresholded images (as a
% proportion of the detected image)

th=max(im,[],'all')*th_p + max(im,[],'all')/20; % threshold (absolute)

im_B=(double(im)>th); % binary image: pixels above threshold

% debug: check cropped image
%figure; imagesc(im_B);


% find limits for cropping
sum_x=sum(im_B,1);    % sum along x for y-dir
sum_y=sum(im_B,2);    % sum along y for x-dir


% debug:
% figure; hold on;
% plot(sum_x); plot(sum_y);

% find the bounds of the thresholded zone
b=find(sum_x>1); % zones along the y direction
if isempty(b)
    error("No zone detected. Increase threshold to select zone.");
end
[m_y,I_y]=max(sum_x);
w_y=b(end)-b(1); % width of zone detected above threshold
bw= round(p_padding.*w_y);
lb_y=max([1,I_y(1)-bw]); % lower bound
hb_y=min([size(im,2),I_y(1)+bw]); % higher bound of region

col_c=round(mean(lb_y,hb_y)); %centre of selected region

b=find(sum_y>1);
if isempty(b)
    error("No zone detected. Increase threshold to select zone.");
end

[m_x,I_x]=max(sum_y);
w_x=b(end)-b(1); % width of zone detected above threshold
bw= round(p_padding.*w_x); % width of zone to select
lb_x=max([1,I_x(1)-bw]);
hb_x=min([size(im,1),I_x(1)+bw]);
row_c=round(mean([lb_x,hb_x]));

% select image region 
im_redsize=double(im(lb_x:hb_x,lb_y:hb_y));


end