%% Matthew Ning

% This calculates PLI as defined by Vinck (in frequency domain)

function [freq, PLV, plvTrial, PLI, pliTrial, WPLI,wpliTrial, Coherency] = calcSync_Vinck(input,varargin)
    % z is the analytic signal computed from the Hilbert Transform fxn: ch x times x trials
    
    Fs = 1000;
    
    timing = [640:1792]; %the "active" portion of each trial
    s = permute(input,[2 1 3]); % the function fft acts on each column, so we want column to represent time
    s = s(timing,:,:);
    [ts,chs,trs] = size(s);
    n = 2^nextpow2(ts);
    Z = zeros(n/2+1,chs,trs);
    Zconj = Z;
    freq = 0:Fs/n:Fs/2;
    cTrial = zeros(chs,chs,trs);
    plvTrial = cTrial;
    pliTrial = cTrial;
    wpliTrial = cTrial;
    
    % Compute F.T.
    for triali = 1:trs
        tempz = s(:,:,triali);
        for chn = 1:chs
            % Scaling from DFT to CFT
            tempzft = fft(tempz(:,chn),n)/Fs;
            % Single-side F.T.
            tempzft = tempzft(1:n/2+1);
            % Multiply by 2 and remove DC component
            tempzft(2:end-1) = 2*tempzft(2:end-1);
            tempzft(1) = 0;
            Z(:,chn,triali) = tempzft;
            Zconj(:,chn,triali) = conj(tempzft);
        end
    end
    
    % Should return vector of instantaneous phase at time t, channel ch
    % conj function is not really needed, can just compute phi1-phi2
    %phi=wrapToPi(atan(imag(Z)./real(Z)));
    %phiconj=wrapToPi(atan(imag(Zconj)./real(Zconj)));
    phi=wrapToPi(atan2(imag(Z),real(Z)));
    phiconj=wrapToPi(atan2(imag(Zconj),real(Zconj)));
    
    % Number of pairs of channels
    PLV = zeros(chs,chs,nargin-1);
    PLI = PLV;
    WPLI = PLV;
    Coherency = PLV;
    
    % avg over frequency first
    for chi = 1:chs
        for chj = 1:chs
            for triali = 1:trs
                tempX = Z(:,chi,triali).*Zconj(:,chj,triali);
                tempEMsq1 = mean(abs(Z(:,chi,triali)).^2,1);
                tempEMsq2 = mean(abs(Z(:,chj,triali)).^2,1);
                tempphidiff = phi(:,chi,triali)-phiconj(:,chj,triali);
                
                cTrial(chi,chj,triali) = mean(tempX,1)/(sqrt(tempEMsq1*tempEMsq2));
                plvTrial(chi,chj,triali) = abs(nanmean(exp(1i*tempphidiff),1));
                pliTrial(chi,chj,triali) = abs(mean(sign(imag(tempX)),1));
                %wpliTrial(chi,chj,triali) = abs(mean(imag(tempX),1))/mean(abs(imag(tempX)),1);
                wpliTrial(chi,chj,triali) = abs(mean(abs(imag(tempX)).*sign(imag(tempX)),1))/mean(abs(imag(tempX)),1);
            end
        end
    end
    
    % avg over trial
    for task = 1:nargin-1
        cTask = cTrial(:,:,varargin{task});
        Coherency(:,:,task) = mean(cTask,3);
        
        plvTask = plvTrial(:,:,varargin{task});
        PLV(:,:,task) = mean(plvTask,3);
        
        pliTask = pliTrial(:,:,varargin{task});
        PLI(:,:,task) = mean(pliTask,3);
        
        wpliTask = wpliTrial(:,:,varargin{task});
        WPLI(:,:,task) = mean(wpliTask,3);
    end
end