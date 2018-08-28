%% Matthew Ning

% This calculates PLI as defined by Stam (in time domain)
% Also please avoid using 3D matrix due to memory size limit!!!
% Actual error msg below:
% Requested 1153x4489x240 (9.3GB) array exceeds maximum array size preference.

function [PLI,pliTrial] = calcSync_Stam(input,varargin)
    %input = [ch x times x trials]
    %x = [time x ch x trials]
    
    timing = [640:1792]; %the "active" portion of each trial
    x = permute(input,[2 1 3]); % the function hilbert acts on each column, so we want column to represent time
    x = x(timing,:,:);
    [ts,chs,trs] = size(x);
    xht = zeros(ts,chs,trs);
    phi = xht;
    pliTrial = zeros(chs,chs,trs);
    
    % Compute Hilbert Transform
    for triali = 1:trs
        tempx = x(:,:,triali);
        tempxht = hilbert(tempx);
        xht(:,:,triali) = tempxht;
        % Phase must be wrapped from [-pi pi] as per Stam paper
        %phi(:,:,triali) = wrapToPi(atan(imag(xht(:,:,triali))./real(xht(:,:,triali))));
        phi(:,:,triali) = wrapToPi(atan2(imag(xht(:,:,triali)),real(xht(:,:,triali))));
    end
    
    % Number of pairs of channels
    PLI = zeros(chs,chs,nargin-1);

    % avg over time first
    for chi = 1:chs
        for chj = 1:chs
            for triali = 1:trs
                tempinstaPHI = phi(:,chi,triali)-phi(:,chj,triali);
                pliTrial(chi,chj,triali) = abs(mean(sign(tempinstaPHI)));
            end
        end
    end
    
    for task = 1:nargin-1
        pliTask = pliTrial(:,:,varargin{task});
        PLI(:,:,task) = mean(pliTask,3);
    end

end