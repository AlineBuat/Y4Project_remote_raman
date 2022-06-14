function moments = image_moments(im)
% returns a matrix of image moments.
% moments higher than second order are normalized by image weight

im=double(im); % make sure image is in double format

M00=sum(im,'all'); % find total weight of image
M=size(im,1);  % size in the m direction (how many rows)
N=size(im,2); % size in the y direction (how many columns)

m_vect=1:M; n_vect=1:N;

%moments:
% Mpq=sum_(n,m){m^p * n^q * im(m,n)}

m_mean=sum(m_vect*im)/M00; %M10
n_mean=sum(im*n_vect')/M00; %M01
M11=m_vect*im*n_vect'; %M11

mu20=sum( ((m_vect-m_mean).^2)*im)/M00;
mu02=sum(im*((n_vect'-n_mean).^2))/M00;




moments=[M00 ,n_mean, mu20; m_mean, M11/M00,0;mu02,0,0];




