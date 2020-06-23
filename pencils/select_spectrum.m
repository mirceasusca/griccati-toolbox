function select = select_spectrum(ev,discr,smarg,stable,tol)
%SELECT_SPECTRUM Given the eigenvalues vector ev and a set of options,
% returns a vector of clusters for the Cg/Cb ordering of the matrix
% spectrum.
%
% select = select_spectrum(ev,discr,smarg,stable), where tola=0
% select = select_spectrum(ev,discr,smarg,stable,tola)
%

if nargin == 4
    tol = 0;
end

if stable

    % stable spectrum
    if discr
        assert(smarg <= 1);
        select = (abs(ev) <= smarg);
    else % 
        assert(smarg <= 0);
        select = (real(ev) <= smarg) | (abs(ev)*tol > 1);
    end

else
    
    % antistable spectrum
    if discr
        assert(smarg >= 1);
        select = (abs(ev) >= smarg);
    else
        assert(smarg >= 0);
        select = (real(ev) >= smarg)  | (abs(ev)*tol > 1);
    end

    
end

end

