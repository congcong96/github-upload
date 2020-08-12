%THRESHOLD2MAT get threshold crossing for .raw files in given directory
%   [THRESH,OUTFILE] = THRESHOLD2MAT(DPATH,NSTD,STIMTYPE,FS) similar in
%   function to spikesort2mat, this function returns threshold crossing
%   times instead of spike times.  NSTD determines the threshold level in
%   number of recording standard deviations.
%
%   Written by Jonathan Shih 2010-11-10
function [thresh,outfile] = intan_raw2thresh(dpath, nstd, stimtype)

%Checking input parameters
if ( nargin ~= 3 )
    error('You need 3 input args.');
end


%Creating THRESH structs
str = struct(...
    'file',        [], ...
    'exp',         [], ...
    'site',        [], ...
    'amplifier',   [], ...
    'chan',        [], ...
    'depth',       [], ...
    'stim',        [], ...
    'atten',       [], ...
    'spiketimes',  [], ...
    'fs',          []);


%Removing backslash at end of directory path if necessary
if ( length(dpath) > 0 )
   if(dpath(end) == '\')
       dpath = dpath(1:(end-1));
   end
end


% Example file: 140910-site1-1606um-30db-rn1-fs20000-A-000.raw
dirfiles = gfn(fullfile(dpath, sprintf('*-%s-*min*.raw', stimtype)), 1);

if isempty(dirfiles)
    dirfiles = gfn(fullfile(dpath, sprintf('*-%s-*.raw', stimtype)), 1);
end

% Don't want the trigger channels, which are 'ADC' files
trigidx = ~contains(dirfiles, 'ADC');
datafiles = dirfiles(trigidx);


% for i = 1:length(datafiles)
%    fprintf('%s\n', datafiles(i).name);
% end % (for i)


% example_file = datafiles(1).name;
% exp = regexp(example_file, '^\d{6}_\d{6}','match','once');
% site = regexp(example_file, 'site\d{1,2}','match','once');
% depth = regexp(example_file, '\w+um','match','once');
% atten = regexp(example_file, '\d{1,4}db','match','once');
% stimlen = regexp(example_file, '\d{1,3}min','match','once');
% probe = regexp(example_file, '(?<=(min-)\S+(?=(-fs))','match','once');
% fs = regexp(example_file, 'fs\d{3,6}','match','once');
% amplifier = regexp(example_file, '\w{1}(?=(-\d{3}.raw))','match','once');
% chan = regexp(example_file, '\d{3}(?=(.raw))','match','once');



% 
% index = findstr(example_file, '-');
% 
% exp = example_file(1:index(1)-1);
% 
% index_site = findstr(example_file, 'site');
% site = str2num( example_file(index_site+4:index(2)-1) );
% 
% index_um = findstr(example_file, 'um');
% depth = str2num( example_file(index(2)+1:index_um-1) );
% 
% index_db = findstr(example_file, 'db');
% depth = str2num( example_file(index(3)+1:index_db-1) );
% 
% index_fs = findstr(example_file, 'fs');
% fs = str2num( example_file(index_fs+2:index(6)-1) );
% 
% amplifier = example_file(index(6)+1:index(7)-1);
% 
% index_dot = findstr(example_file, '.raw');
% chantemp = str2num(example_file(index(7)+1:index_dot-1) );


thresh = [];
%Going through each channel to get the threshold crossing times
ii = length(datafiles);
rawfile = datafiles{ii};
[~,b,c] = fileparts(rawfile);
basefile = [b c];
basename = regexp(basefile, '\S+(?=(-\w{1}-\d{3}.raw))','match','once');
outfile = fullfile(dpath, [basename '-thresh.mat']);
d = dir(outfile);
if ~isempty(d)
    thresh = [];
    return
end
for ii = 1:length(datafiles)
    %Current file being processed
    rawfile = datafiles{ii};
    [~,b,c] = fileparts(rawfile);
    basefile = [b c];
    
    details = regexp(basefile, ['(?<exp>\d{6}_\d{6})-site(?<site>\d{1,2})-'...
        '(?<depth>\d{3,})um-(?<atten>\w+db)-' sprintf('(?<stim>%s)-',stimtype)...
        '(?<stimlen>\d{1,2})min-(?<probe>\S+)-fs'...
        '(?<fs>\d{3,6})-(?<amplifier>[ABCD])-(?<chan>\d{3}(?=(.raw)))'],...
        'names');
    
    if isempty(details)
        
        details = regexp(basefile, ['(?<exp>\d{6}_\d{6})-site(?<site>\d{1,2})-'...
            '(?<depth>\d{3,})um-(?<atten>\w+db)-' sprintf('(?<stim>%s)-',stimtype)...
            '(?<probe>\S+)-fs'...
            '(?<fs>\d{3,6})-(?<amplifier>[ABCD])-(?<chan>\d{3}(?=(.raw)))'],...
            'names');
        
    end
    
    if isempty(details)
        
        details = regexp(basefile, ['(?<exp>\d{6}_\d{6})-site(?<site>\d{1,2})-'...
            '(?<depth>\d{3,})um-(?<atten>\w+db)-' sprintf('(?<stim>%s)-',stimtype)...
            'fs(?<fs>\d{3,6})-(?<amplifier>[ABCD])-(?<chan>\d{3}(?=(.raw)))'],...
            'names');
        
    end
    
    
%     exp = regexp(rawfile, '^\d{6}_\d{6}','match','once');
%     site = regexp(rawfile, 'site\d{1,2}','match','once');
%     depth = regexp(rawfile, '\w+um','match','once');
%     atten = regexp(rawfile, '\d{1,4}db','match','once');
%     stimlen = regexp(rawfile, '\d{1,3}min','match','once');
%     probe = regexp(rawfile, '(?<=(min-))\S+(?=(-fs))','match','once');
%     fs = str2double(regexp(rawfile, '(?<=(fs))\d{3,6}','match','once'));
%     amplifier = regexp(rawfile, '\w{1}(?=(-\d{3}.raw))','match','once');
%     chan = regexp(rawfile, '\d{3}(?=(.raw))','match','once');
%    [exp, site, depth, atten, stim, probe, fs, amplifier, chan] = ...
%       textscan(rawfile,'%ssite%d%dum%ddb%sfs%d%s%d.raw','Delimiter','-');
% 
%    % Convert from cell array if needed
%    if (iscell(exp)), exp = exp{1}; end;
%    if (iscell(site)), site = site{1}; end;
%    if (iscell(depth)), depth = depth{1}; end;
%    if (iscell(atten)), atten = atten{1}; end;
%    if (iscell(fs)), fs = fs{1}; end;
%    if (iscell(amplifier)), amplifier = amplifier{1}; end;
%    if (iscell(chan)), chan = chan{1}; end;
   
   fs = str2double(details.fs);
   fidin = fopen(rawfile, 'r');
   nsamples = 60 * fs;
   [samples, count] = fread(fidin, nsamples, 'int16');
   fclose(fidin);


    %Getting STD of .raw file
    stdsample = std(samples);

    
    %Initializing struct
    s = str;
    s.file = basefile;
    s.exp = details.exp;
    s.site = str2double(details.site);
    s.amplifier = details.amplifier;
    s.depth = str2double(details.depth);
    s.chan = str2double(details.chan); % assign channel number
    s.stim = stimtype;
    if isfield(details, 'stimlen')
        s.stimlen = str2double(details.stimlen);
    end
    s.atten = details.atten;
    if isfield(details, 'probe')
        s.probe = details.probe;
    else
        s.probe = 'a1x32-poly3';
    end
    s.spiketimes = get_threshold_crossings(rawfile,nstd*stdsample,fs);
    s.fs = fs; % save sampling frequency
    

    %Updating user
    fprintf('Finished %s\n', basefile);
    
    thresh = [thresh s];
    clear('s');
end




% outfile = sprintf('%s-site%d-%dum-%ddb-%s-fs%d-thresh.mat', ...
%    exp, site, depth, atten, stimtype, fs);

return;




%Helper Function: get_threshold_crossings
%----------------------------------------
%Gets the threshold crossings from a .raw file and converts them into msec
function crossings = get_threshold_crossings(rawpath,thresh,fs)

%Opening .raw file
fid = fopen(rawpath,'r');

%Number of data values to read in each iteration
blockSize = 10*1024;

init_value = 0;
nblocks = 0;
crossings = [];
while(~feof(fid))
    %Getting current block of data
    blockData = fread(fid,blockSize,'int16');
    blockData = -blockData; % invert waveform because large spike peaks are downward in Intan
                            % software; they are upward in neuralynx
                            % software
    
    %Appending last value of previous block to current block
    blockData = [init_value; blockData];
    
    %Finding threshold crossings in the current block
    ind_supra = find(blockData(2:end) > thresh);
    ind_cross = ind_supra(find(blockData(ind_supra) <= thresh));
    crosstimes = (nblocks*blockSize + ind_cross - 1)/fs*1000; %converting samples into msec
    crossings = [crossings; crosstimes];
    
    %Updating value to append to beginning of next block
    init_value = blockData(end);
    
    %Updating number of blocks
    nblocks = nblocks + 1;
    
%     plot(blockData), hold on;
%     plot(xlim,thresh*[1 1],'k--');
%     scatter(ind_cross,blockData(ind_cross+1),'rx'), hold off;
%     pause;
end

%Closing .raw file
fclose(fid);

return;



