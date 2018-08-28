%% Plot EEG topography of connectivity

% Readapted by Matthew Ning
% Adapted by Jasmine Kwasa
% Originally written by Matthew Ning
% Auditory Neuroscience Lab
% Professor Barbara Shinn-Cunningham

close all
clear 
clc

filename = 'chanlocs_64.csd';
content = fileread(filename);
data = textscan(content, '%s %f %f %f %f %f %f %s', ...
                     'CommentStyle','//');

subjINFO = 'ASA067_XX_NETWORK_PLI_All.mat';
home = pwd;
cd ..
load(subjINFO)
cd(home)
titleName = subjINFO;
stuff = zeros(64,64,3);
stuff(:,:,1) = mean(PLI(1:64,1:64,FN),3);
stuff(:,:,2) = mean(PLI(1:64,1:64,BN),3);
stuff(:,:,3) = mean(PLI(1:64,1:64,FN),3)-mean(PLI(1:64,1:64,BN),3);

%%
for k = 1:3
    phittemp = squeeze(stuff(:,:,k));

    for i = 1: size(phittemp,1)
        for j = 1: size(phittemp,2)
            if(phittemp(i,j)==1)
                phittemp(i,j) = -Inf;
            end
        end
    end
    phittemporary = phittemp;


    % Plot Matt's way
    % numch = size(Phi,1);
    % figure
    % set(gcf,'Name',sprintf('%d-channel EEG Montage',numch));
    % m = 100;
    % t = [0:pi/100:2*pi]; 
    % r = m/2 + 0.5;
    % head = [sin(t)*r + m/2+1; cos(t)*r + m/2+1]' - m/2;
    % scrsz = get(0,'ScreenSize');
    % d = min(scrsz(3:4)) / 2;
    % set(gcf,'Position',[scrsz(3)/2 - d/2 scrsz(4)/2 - d/2 d d]); 
    % %whitebg('w');
    % axes('position',[0 0 1 1]);
    % set(gca,'Visible','off');
    % line(head(:,1),head(:,2),'Color','k','LineWidth',1); 
    % mark = '\bullet';
    % l = sqrt((Xcoor-0.5).^2 + (Ycoor-0.5).^2) * 2;
    % r = (r - 3.5) / (max(l) / max([max(Xcoor) max(Ycoor)]));

    % colormap jet
    % colmap = jet(256);
    % hold on
    % 
    % % grab the top 3 for each channel
    % [val, chind] = max(phittemp,[],2); 
    % colorval = ceil(val*256);
    % phittemp(chind(chind),chind) = -Inf;
    % [val2, chind2] = max(phittemp,[],2);
    % colorval2 = ceil(val2*256);
    % phittemp(chind2(chind2),chind2) = -Inf;
    % [val3, chind3] = max(phittemp,[],2);
    % colorval3 = ceil(val3*256);



    % for e = 1:numch
    %     plot([Xcoor(e)*2*r - r + 0.5 Xcoor(chind(e))*2*r - r + 0.5], [Ycoor(e)*2*r - r + 2.5 Ycoor(chind(e))*2*r - r + 2.5], 'color', colmap(colorval(e),:),'LineWidth', 2);
    %     plot([Xcoor(e)*2*r - r + 0.5 Xcoor(chind2(e))*2*r - r + 0.5], [Ycoor(e)*2*r - r + 2.5 Ycoor(chind2(e))*2*r - r + 2.5], 'color', colmap(colorval2(e),:),'LineWidth', 2);
    %     plot([Xcoor(e)*2*r - r + 0.5 Xcoor(chind3(e))*2*r - r + 0.5], [Ycoor(e)*2*r - r + 2.5 Ycoor(chind3(e))*2*r - r + 2.5], 'color', colmap(colorval3(e),:),'LineWidth', 2);
    % 
    %     text(Xcoor(e)*2*r - r + 0.5,Ycoor(e)*2*r - r + 2.5,mark);
    %     text(Xcoor(e)*2*r - r + 1, ...
    %          Ycoor(e)*2*r - r , ...
    %          Label(e), ...
    %           'FontSize',8, ...
    %           'FontWeight','bold', ...
    %           'VerticalAlignment','middle', ...
    %           'HorizontalAlignment','center');
    % 
    % end
    % 
    % hold off;
    % colorbar;


    % Jasmine's way: Plot the most coherent... 
    phittemporary(phittemporary==0) = NaN;
    phittemporary(isnan(phittemporary)) = -Inf;

    Label = data{1};
    Theta = data{2};
    Phi = data{3};
    Radius = data{4};
    Xcoor = data{5};
    Ycoor = data{6};
    phiT = 90 - Phi;                    % calculate phi from top of sphere
    theta2 = (2 * pi * Theta) / 360;    % convert degrees to radians
    phi2 = (2 * pi * phiT) / 360;
    [x,y] = pol2cart(theta2,phi2);      % get plane coordinates
    xy = [x y];
    xy = xy/max(max(xy));               % set maximum to unit length
    xy = xy/2 + 0.5;                    % adjust to range 0-1
    Xcoor = xy(:,1);
    Ycoor = xy(:,2);

    figure
    set(gcf,'Name',sprintf('%d coherence',k));
    m = 100;
    t = [0:pi/100:2*pi]; 
    r = m/2 + 0.5;
    head = [sin(t)*r + m/2+1; cos(t)*r + m/2+1]' - m/2;
    scrsz = get(0,'ScreenSize');
    d = min(scrsz(3:4)) / 2;
    set(gcf,'Position',[scrsz(3)/2 - d/2 scrsz(4)/2 - d/2 d d]); 
    %whitebg('w');
    axes('position',[0 0 1 1]);
    set(gca,'Visible','off');
    line(head(:,1),head(:,2),'Color','k','LineWidth',1); 
    mark = '\bullet';
    l = sqrt((Xcoor-0.5).^2 + (Ycoor-0.5).^2) * 2;
    r = (r - 3.5) / (max(l) / max([max(Xcoor) max(Ycoor)]));
    colormap jet
    colmap = jet(256);
    hold on

    %plot it
    hold on
    xcords = Xcoor*2*r - r + 0.5;
    ycords = Ycoor*2*r - r + 2.5;

    allPositiveValues = abs(phittemporary(:));
    
    [yy,ii]=sort(allPositiveValues,'ascend'); %search for max of full matrix
    
    %calculate the original position of those maxima in the matrix... 
    [row,column]=ind2sub([64 64],ii);

%     normYY = yy;
    normYY = yy/max(yy(~isinf(yy)));
    colorval = ceil(normYY*256);
    f = 1;

    while ~isinf(yy(f))
        hold on
        plot([xcords(row(f)) xcords(column(f))],[ycords(row(f)) ycords(column(f))],...
            'color', colmap(colorval(f),:),'LineWidth', 2)
        f = f+1;
    end
    text(Xcoor*2*r - r + 0.5,Ycoor*2*r - r + 2.5,mark);
    text(Xcoor*2*r - r + 1, Ycoor*2*r - r , Label, ...
          'FontSize',8, ...
          'FontWeight','bold', ...
          'VerticalAlignment','middle', ...
          'HorizontalAlignment','center');
    colorbar

    title(num2str(k))
    caxis([0 max(yy(~isinf(yy)))])

end