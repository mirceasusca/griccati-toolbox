function [krd,krc] = extract_kronecker_structure(info,structure)
%EXTRACT_REGULAR_PART
% krd - kronecker row dimension
% krc - kronecker column dimension

if nargin < 2
    structure = 'right';
end

if strcmp(structure,'right')
    ridx = 1;
    cidx = 1;
    if ~isempty(info.mr)
        for i = 1:length(info.mr)
            ridx = ridx + (i-1)*(info.nr(i)-info.mr(i));
            cidx = cidx + i*(info.nr(i)-info.mr(i));
        end
    end

    krd = ridx-1;
    krc = cidx-1;
elseif strcmp(structure,'left')
    ridx = 1;
    cidx = 1;
    if ~isempty(info.ml)
        ll = length(info.ml);
        for i = ll:-1:1
            cidx = cidx + (ll-i)*(info.ml(i)-info.nl(i));
            ridx = ridx + (ll-(i-1))*(info.ml(i)-info.nl(i));
        end
    end

    krd = ridx-1;
    krc = cidx-1;
else
    error('Unknown requested structure.');
end

end
