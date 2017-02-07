function S = taylorstats(r, f)
%TAYLORSTATS Calculate statistics needed for a Taylor diagram
%
% S = taylorstats(r, f)
%
% This function calculates the statistics needed to produce a Taylor
% diagram, as described in Taylor, 2001.  These diagrams can be used to
% visualize closeness of fit between multiple models and a target dataset.
%
% Input variables:
%
%   r:  n x 1 array, values of reference (observed) variable
%
%   f:  n x m array, values of modeled variable, for m different models
%
% Output variables:
%
%   S:  1 x 1 structure with the following fields; all are 1 x m+1 arrays,
%       where the first element corresponds to the reference data and 2:end
%       correspond to the models.
%
%       std:    standard deviation, normalized to n.
%
%       cor:    correlation coefficient between f and r
%
%       rmsd:   root mean squared difference between f and r
%
%       rms:    centered pattern RMS
%
%       bias:   bias (not needed for plot, but a useful stat

% Copyright 2013-2017 Kelly Kearney

[npt, nf] = size(f);
if ~isvector(r) || ~isequal(npt, length(r))
    error('Length of r must match size(f,1)');
end

% Standard deviation

S.std = [std(r,1) std(f,1)];

% Correlation coef

fbar = mean(f);
rbar = mean(r);

S.cor(1) = 1;
for ii = 1:nf
    S.cor(ii+1) = (sum((f(:,ii) - fbar(ii)) .* (r - rbar))./npt)./(S.std(1).*S.std(ii+1));
end

% RMSD

S.rmsd = sqrt(sum((bsxfun(@minus, f, r)).^2)/npt);

% Centered pattern RMS

S.rms(1) = 0;
for ii = 1:nf
    S.rms(ii+1) = sqrt(sum(((f(:,ii) - fbar(ii)) - (r - rbar)).^2)./npt);
end

% Bias

S.bias = [0 bsxfun(@minus, mean(f,1), mean(r))];
