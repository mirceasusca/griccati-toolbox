function [V,W,J,F,X] = kpys(sigma,X)
%KPYS Solves KPYS or KSPYS for the popov triplet "sigma".
%
% [V,W,J,F] = kpys(sigma)
% [V,W,J,F] = kpys(sigma,X), if the Riccati solution is already precomputed
%

if sigma.discr == true
    if nargin == 2
        [V,W,J,F,X] = kspys_slv(sigma,X);
    else
        [V,W,J,F,X] = kspys_slv(sigma);
    end
else
    if nargin == 2
        [V,W,J,F,X] = kpys_slv(sigma,X);
    else
        [V,W,J,F,X] = kpys_slv(sigma);
    end
end

end

