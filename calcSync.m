%% Matthew Ning
% Auditory Neuroscience Lab
% Professor Barbara Shinn-Cunningham
% Supervised by Jasmine Kwasa

function [phi, phidiff, phipair,PLIpair,WPLIpair,Cpair] = calcSync(trialiHT,CH)
    % Should return vector of instantaneous phase at time t, channel ch
    % conj function is not really needed, can just compute phi1-phi2
    trialiHTconj = conj(trialiHT);
    phi=atan(imag(trialiHT)./real(trialiHT));
    phiconj=atan(imag(trialiHTconj)./real(trialiHTconj));
    
    % Number of pairs of channels
    %pairs = CH*(CH-1)/2 + CH;
    pairs = CH^2;
    lent = size(trialiHT,2);
    phipair = zeros(CH);
    PLIpair = phipair;
    WPLIpair = phipair;
    Cpair = phipair;
    phidiff = zeros(pairs,lent);
    trialiPLV = zeros(pairs,1);
    trialiPLI = trialiPLV;
    trialiWPLI = trialiPLV;
    trialiC = trialiPLV;
    for chi = 1:CH
        for chj = 1:CH
        %for chj = chi:CH
            idex = sub2ind(size(phipair),chi,chj);
            x = trialiHT(chi,:).*trialiHTconj(chj,:);
            phidiff(idex,:) = phi(chi,:)+phiconj(chj,:);
            phidiff(idex,:) = 1i.*phidiff(idex,:);
            trialiPLV(idex) = abs(mean(exp(phidiff(idex,:)),2));
            trialiPLI(idex) = abs(mean(sign(imag(x)),2));
            trialiWPLI(idex) = abs(mean(imag(x),2))/mean(abs(imag(x)),2);
            trialiC(idex) = mean(x,2)/sqrt(mean(abs(trialiHT(chi,:)).^2,2)*(mean(abs(trialiHT(chj,:)).^2,2)));
        end
        meanphi = mean(phi(chi,:),2);
    end
    
    % Find best pair for analysis
    PLVsorted = sort(trialiPLV,'descend');
    PLIsorted = sort(trialiPLI,'descend');
    WPLIsorted = sort(trialiWPLI,'descend');
    Csorted = sort(trialiC,'descend');
    
    [CHI,CHJ] = ind2sub([CH CH],(1:pairs));
    for id = 1:pairs
        phipair(CHI(id),CHJ(id)) = trialiPLV(id);
        PLIpair(CHI(id),CHJ(id)) = trialiPLI(id);
        WPLIpair(CHI(id),CHJ(id)) = trialiWPLI(id);
        Cpair(CHI(id),CHJ(id)) = trialiC(id);
    end
end