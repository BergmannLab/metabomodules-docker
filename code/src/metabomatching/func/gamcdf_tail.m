function [p,plo,pup] = gamcdf_tail(x,a,b,pcov,alpha)
%GAMCDF_TAIL Tail of gamma cumulative distribution function.

if nargin < 2
    error('stats:gamcdf:TooFewInputs',...
          'Requires at least two input arguments.');
elseif nargin < 3
    b = 1;
end

% More checking if we need to compute confidence bounds.
if nargout > 1
    if nargin < 4
        error('stats:gamcdf:TooFewInputs',...
              'Must provide covariance matrix to compute confidence bounds.');
    end
    if ~isequal(size(pcov),[2 2])
        error('stats:gamcdf:BadCovariance',...
              'Covariance matrix must have 2 rows and columns.');
    end
    if nargin < 5
        alpha = 0.05;
    elseif ~isnumeric(alpha) || numel(alpha) ~= 1 || alpha <= 0 || alpha >= 1
        error('stats:gamcdf:BadAlpha',...
              'ALPHA must be a scalar between 0 and 1.');
    end
end

% Return NaN for out of range parameters.
a(a <= 0) = NaN;
b(b <= 0) = NaN;
x(x < 0) = 0;

try
    z = x ./ b;
    p = gammainc(z, a, 'upper');
catch
    error('stats:gamcdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end
p(z == Inf) = 1;

% Compute confidence bounds if requested.
if nargout >= 2
    % Approximate the variance of p on the logit scale
    logitp = log(p./(1-p));
    dp = 1 ./ (p.*(1-p)); % derivative of logit(p) w.r.t. p
    da = dgammainc(z,a) .* dp; % dlogitp/da = dp/da * dlogitp/dp
    db = -exp(a.*log(z)-z-gammaln(a)-log(b)) .* dp; % dlogitp/db = dp/db * dlogitp/dp
    varLogitp = pcov(1,1).*da.^2 + 2.*pcov(1,2).*da.*db + pcov(2,2).*db.^2;
    if any(varLogitp(:) < 0)
        error('stats:gamcdf:BadCovariance',...
              'PCOV must be a positive semi-definite matrix.');
    end
    
    % Use a normal approximation on the logit scale, then transform back to
    % the original CDF scale
    halfwidth = -norminv(alpha/2) * sqrt(varLogitp);
    explogitplo = exp(logitp - halfwidth);
    explogitpup = exp(logitp + halfwidth);
    plo = explogitplo ./ (1 + explogitplo);
    pup = explogitpup ./ (1 + explogitpup);
end
