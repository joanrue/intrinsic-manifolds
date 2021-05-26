function Xsurr = phaseRandomize(X)
% Returns a phase-randomized surrogate of the multivariate time-series X. 
%   INPUT:
%           * X: matrix(ROIs, samples) containing fMRI time series

    [N,L]=size(X);
    % Get spectrum
    Y=fft(X,[],2);
    
    % Add random phase shifts (negative for conjugates) and preserve DC
    % offset. Same for each ROI
    rnd_theta= -pi + (2*pi).*rand(1,L/2-1);
    rnd_theta= repmat(rnd_theta,N,1); 
    
    % Make negative frequencies' phases equal to their positive counterpart
    % to have a real surrogate signal. 
    Y(:,2:L/2)=Y(:,2:L/2).*exp(1i*rnd_theta);
    Y(:,L/2+2:L)=Y(:,L/2+2:L).*exp(-1i*flip(rnd_theta,2));
    
    % return phase-randomized data
    Xsurr =ifft(Y,[],2);
end