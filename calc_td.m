function [td] = calc_td(tp,n,DC)

% case ** page 99 Labbook 1
if n<2
    td=0;
else
td=tp*(1/DC-1);
end;
% case * page 99 Labbook 1
% td=n/(n-1)*tp*(1/DC-1);
