function [wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo(path, name)
% parseHdrInfo
% [wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo(path,filename)
% tint = integration time
% settings = band intensities for Acuity files

                %Enter a path (e.g. 'M:\Engineering\...') and .hdr filename
f = fullfile(path, name);       
fp = fopen(f, 'r');             %Opens file as read-only

spatial =0;                     %Initialize spatial and spectral
spectral=0;
waveFlag = 0;                   %Flag changes to 1 when wavelengths vector is collected
wavelengths = zeros(1,spectral);
settings = zeros(1,spectral);

line = fgetl(fp);               %Sets variable line to the first line (and then the next, the next) of the file
while(~feof(fp)) && (waveFlag == 0)   %While not end of file and spatial or waveflag is 0
    ind = strfind(line, 'samples');
    if numel(ind)
        spatial = str2double(line((strfind(line, '=')+1):size(line,2)));
    else
        ind = strfind(line, 'lines');
        if numel(ind) && ~exist('frames','var')
            frames = str2double(line((strfind(line,'=')+1):size(line,2)));
        else
            ind = strfind(line, 'bands');
            if (ind == 1)
                spectral = str2double(line((strfind(line, '=')+1):size(line, 2)));
            else
                ind = strfind(line, 'tint');
                ind2 = strfind(line, 'Exposure time');
                ind3 = strfind(line, 'Integration Time');
                if numel(ind) || numel(ind2) || numel(ind3)
                    tint = str2double(line((strfind(line, '=')+1):size(line, 2)));
                else
                    ind = strfind(line, 'Wavelength');
                    ind2 = strfind(line,'wavelength');     %If line has the word Wavelength, ind = character # of location of 'Wavelength' (1, b/c line starts w/"Wavelength")
                    if numel(ind) || numel(ind2)
                        i = 1;
                        while i<=spectral
                            line=fgetl(fp);
                            wlgths = str2num(line); %allows multiple wavelengths on one line
                            N = numel(wlgths);
                            wavelengths(1,i:i+N-1) = wlgths;
                            i = i+N;
                        end
                        
%                         line=fgetl(fp);
%                         for k = 1:spectral
%                             wavelengths(1,k) = str2double(line(1:size(line,2)));
%                             %Str2double converts argument text > vector
%                             %Use characters from beginning to end of line (2nd dim.)
%                             line=fgetl(fp);
%                         end
%                         waveFlag = 1;
                    else
                        ind = strfind(line, 'Band Intensities');
                        if numel(ind)
                            line=fgetl(fp);
                            for k = 1:spectral
                                settings(1,k) = str2double(line(1:size(line,2)));
                                %Str2double converts argument text > vector
                                %Use characters from beginning to end of line (2nd dim.)
                                line=fgetl(fp);
                            end
                            waveFlag = 1;
                        end
                    end
                end
            end
        end
    end
    line=fgetl(fp);
end
fclose(fp);

wavelengths = sort(wavelengths);
    
                       

                    
    
    
    
    