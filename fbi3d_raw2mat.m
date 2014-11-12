% Tool to convert raw images into .mat file for MATLAB
% Joaquin L. Herraiz (herraiz@mit.edu)
% This tool is part of the FBI3D code

% -- IMAGE PARAMETERS 
RES = 75;
NZS = 61;

% -- INPUT IMAGES ----
filename='image_ref.raw';
fid = fopen(filename,'r');
IMAGE_REF = reshape(fread(fid,Inf,'float'),[RES RES NZS]);
fclose(fid);

filename='image_ref_senscor.raw';
fid = fopen(filename,'r');
IMAGE_REF_NORM = reshape(fread(fid,Inf,'float'),[RES RES NZS]);
fclose(fid);

filename='image.raw';
fid = fopen(filename,'r');
IMAGE = reshape(fread(fid,Inf,'int'),[RES RES NZS]);
fclose(fid);

filename='image_norm.raw';
fid = fopen(filename,'r');
IMAGE_NORM = reshape(fread(fid,Inf,'float'),[RES RES NZS]);
fclose(fid);

filename='image_norm_med.raw';
fid = fopen(filename,'r');
IMAGE_NORM_MED = reshape(fread(fid,Inf,'float'),[RES RES NZS]);
fclose(fid);

save('FBI3D_output.mat')
