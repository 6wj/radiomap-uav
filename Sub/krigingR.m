function [zi,s2zi] = krigingR(vstruct,x,z,xi,x0,chunksize)

% size of input arguments
sizest = size(xi(:,1));
numest = length(xi);
numobs = length(x);
if ~exist('chunksize', 'var')
    chunksize = 100;
end
global userPos0
userPos0 = x0;
% check if the latest version of variogramfit is used
if ~isfield(vstruct, 'func')
    error('please download the latest version of variogramfit from the FEX')
end
% variogram function definitions
switch lower(vstruct.model)    
    case {'whittle' 'matern'}
        error('whittle and matern are not supported yet');
    case 'stable'
        stablealpha = vstruct.stablealpha; %#ok<NASGU> % will be used in an anonymous function
end
% distance matrix of locations with known values
Dx = pdist2(x,x,@distpolar);
% if we have a bounded variogram model, it is convenient to set distances
% that are longer than the range to the range since from here on the
% variogram value remains the same and we don�t need composite functions.
switch vstruct.type
    case 'bounded'
        Dx = min(Dx,vstruct.range);
    otherwise
end
% now calculate the matrix with variogram values 
A = vstruct.func([vstruct.range vstruct.sill],Dx);
if ~isempty(vstruct.nugget)
    A = A+vstruct.nugget;
end
% the matrix must be expanded by one line and one row to account for
% condition, that all weights must sum to one (lagrange multiplier)
A = [[A ones(numobs,1)];ones(1,numobs) 0];
% A is often very badly conditioned. Hence we use the Pseudo-Inverse for
% solving the equations
A = pinv(A);
% we also need to expand z
z  = [z;0];
% allocate the output zi
zi = nan(numest,1);
if nargout == 2
    s2zi = nan(numest,1);
    krigvariance = true;
else
    krigvariance = false;
end
% parametrize engine
nrloops   = ceil(numest/chunksize);
% now loop 
for r = 1:nrloops
    % built chunks
    if r<nrloops
        IX = (r-1)*chunksize +1 : r*chunksize;
    else
        IX = (r-1)*chunksize +1 : numest;
        chunksize = numel(IX);
    end
    
    % build b
    b = pdist2(x,xi(IX,:),@distpolar);
    % again set maximum distances to the range
    switch vstruct.type
        case 'bounded'
            b = min(vstruct.range,b);
    end
    
    % expand b with ones
    b = [vstruct.func([vstruct.range vstruct.sill],b);ones(1,chunksize)];
    if ~isempty(vstruct.nugget)
        b = b+vstruct.nugget;
    end
    
    % solve system
    warning off;
    lambda = A*b;
    warning on;
    
    % estimate zi
    zi(IX)  = lambda'*z;
    
    % calculate kriging variance
    if krigvariance
        s2zi(IX) = sum(b.*lambda,1);
    end
    
end
% reshape zi
zi = reshape(zi,sizest);
if krigvariance
    s2zi = reshape(s2zi,sizest);
end
