function [hac, parameter] = read_hac(filename)
% Read Heat Flux AutoCorelation Function file and out put hac matrix
%   Input: filename (char)
%   Output: reshaped hac matrix (double) & parameter structure (double)
%           `parameter` includes: Tstart, Nc, Ns, M

[parameter.Tstart, parameter.Nc] = textread(filename, '%n%n', 1, 'headerlines', 3);
[index, TimeDelta, Ncount, c_1, c_2, c_3] = ...
    textread(filename, '%n%n%n%n%n%n', 'headerlines', parameter.Nc+4, 'emptyvalue', 0);

parameter.Ns = TimeDelta(3,1) - TimeDelta(2, 1);

hac = [index TimeDelta Ncount c_1 c_2 c_3];
clear index TimeDelta Ncount c_1 c_2 c_3;


nall = length(hac) / (parameter.Nc+1);
for j = 1 : nall
    hac((j-1)*parameter.Nc+1, :) = [];
end

parameter.M = size(hac, 1) / parameter.Nc;

hac = reshape(mean(hac(:, 6), 2), parameter.Nc, parameter.M);

end