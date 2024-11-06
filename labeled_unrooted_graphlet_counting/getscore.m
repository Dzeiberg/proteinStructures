function [s] = getscore (string, alphabet)

% This function gives a strong and an alphabet and it convers this string
% into a score using the base system of the size of alphabet. We need this
% function to tell us what is the automorphic structure to actually be
% counted.

s = find(alphabet == string(1));
for i = 2 : length(string)
    q = find(alphabet == string(i));
    s = (s - 1) * length(alphabet) + q;
end

%s = string(1);
%for i = 2 : length(string)
%    s = (s - 1) * length(alphabet) + string(i);
%end


return