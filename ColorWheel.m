close all
clear all
clc

load('Vars.mat');

filename = 'chanlocs_64.csd';
content = fileread(filename);
data = textscan(content, '%s %f %f %f %f %f %f %s', ...
                     'CommentStyle','//');
                 
% Set channel number
ch = 38;

% Set parameters (these could be arguments to a function)
rInner = 50;     % inner radius of the colour ring
rOuter = 200;    % outer radius of the colour ring

% Get polar coordinates of each point in the domain
[x, y] = meshgrid(-rOuter:rOuter);
[theta, rho] = cart2pol(x, y);

huelin = phitFS1m(ch,:);
hue = theta;

for n = 1:size(huelin,2);
    % convert channel number to theta
    u = n*(2*pi/67)+(-pi);
    l = (n-1)*(2*pi/67)+(-pi);
    for i = 1:size(hue,1);
        for j = 1:size(hue,2);
            if (theta(i,j)>= l && theta(i,j) < u)
                hue(i,j) = huelin(1,n);
            end
        end
    end
end

% for i = 1:size(hue,1);
%     for j = 1:size(hue,2);
%         % convert theta to channel number
%         chi = 1;
%         hue(i,j) = huelin(1,chi);
%     end
% end
%hue = phitFS1m(38,:);
%hue = (theta + pi) / (2 * pi);     % hue into range (0, 1]
%hue = ceil(hue * ncols) / ncols;   % quantise hue 
saturation = ones(size(hue));      % full saturation
brightness = double(rho >= rInner & rho <= rOuter);  % black outside ring

% Convert to rgb space for display
rgb = hsv2rgb(cat(3, hue, saturation, brightness));

% Check result
figure(1);
imshow(rgb);
title(['PLV: Channel Number: ',data{1}{ch}]);

huelinchop = huelin(1,1:64);
figure(2);
pie(huelinchop,data{1});
title(['PLV: Channel Number: ',data{1}{ch}]);
