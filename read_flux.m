function [flux] = read_flux(filename,runsteps,runtimes)
% Read Heat Flux data from LAMMPS `output` file and output flux matrix
%   Input: filename (char), run steps (double) & run times (double)
%   Output: reshaped flux matrix (double) 

% NOTE: `runsteps` and `runtimes` should not be incomed via manual input
%       instead, they should be incomed via file

try
    dump = fopen(filename,'r');
catch
    error('Dumpfile not found!');
end
fclose(dump);

% NOTE: it should use `textscan` to speed up, however it's hard to
%       identify the headerlines albeit some unnecessary parameters are
%       inputted

fileID = fopen(filename,'r');

format longG;
flux = zeros(runsteps,3,runtimes); % preallocation RAM

i = 1; % i = 1 : runtimes
while feof(fileID) == 0
    line_content = fgetl(fileID);
    switch line_content
        case 'Step Temp c_flux[3] ' % DONOT REMOVE THE SPACE AT THE END OF THE LINE
            for j = 1 : runsteps
                flux_content = fgetl(fileID);
                flux(j, :, i) = cell2mat(textscan(flux_content,'%n %n %n'));
            end
            i = i + 1;
    end
end

fclose(fileID);

end

