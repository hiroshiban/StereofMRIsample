function field=nf_CreateWedgesField(fieldSize,radiusMin,radiusMax,height,nwedges,...
                                    wedgeangle,rot_angle,pix_per_deg,fine_coefficient)

% function field=nf_CreateWedgesField(fieldSize,radiusMin,radiusMax,height,nwedges,...
%                                     wedgeangle,rot_angle,pix_per_deg,fine_coefficient)
%
% Creates wedge-shaped (Baumkuchen) field
% modified from CreateWedgesField for speeding-up the processing
% 
% [input]
% fieldSize   : the whole image size in deg, [row,col]
% radiusMin   : the min size of the wedge in degrees, [val]
% radiusMax   : the max size of the wedge in degrees, [val]
% height      : field height, [val]
% nwedges     : the number of wedges in the field, [val]
% wedgeangle  : angle of each wedge in deg, [val]
% rot_angle   : wedge rotation angle along the center of the field in deg, [val]
% pix_per_deg : pixels per degree.
% fine_coefficient : (optional) if larger, the generated field becomes finer
%                    along x-axis but comsumes much CPU power. [val]
%                    (default=1, as is, no tuning)
% 
% [output]
% field       : generated wedge field
%               the drawing start from right horizontal meridian, counter-clockwise
%  
% !!! NOTICE !!!
% for speeding up image generation process,
% I will not put codes for nargin checks.
% Please be careful.
% 
% Created    : "2010-08-05 01:15:18 ban"
% Last Update: "2010-08-05 10:27:24 ban"

% convert from deg to pixels
fieldSize=round(fieldSize.*pix_per_deg);

radiusMin=round(radiusMin/2*pix_per_deg);
radiusMax=round(radiusMax/2*pix_per_deg);

% calculate distance & angles
step=1/fine_coefficient;
[x,y]=meshgrid(0:step:fieldSize(2),0:1:fieldSize(1));
x=x-fieldSize(2)/2; y=y-fieldSize(1)/2;
if mod(size(x,1),2), x=x(1:end-1,:); y=y(1:end-1,:); end
if mod(size(x,2),2), x=x(:,1:end-1); y=y(:,1:end-1); end
r=sqrt(x.*x+y.*y);
theta=mod(180*atan2(y,x)./pi+360+rot_angle,360); % adjust 0-360 deg

% calculate gap between adjacent wedges
gap=(360-nwedges*wedgeangle)/nwedges;

% generate wedge field
field=zeros(size(y));
for ii=1:1:nwedges % wedges
  field( (ii*gap + (ii-1)*wedgeangle <= theta) & ...
         (theta <= ii*gap + ii*wedgeangle) & ...
         (radiusMin < r) & (r < radiusMax) )=height;
end

return
