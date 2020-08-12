function [trigger, triptrig] = intan_findtrig_longfile(trigpath, Thresh, M)
%findtrig_longfile  Find triggers from a .raw or .mich int16 file
%
%  [trigger, triptrig] = findtrig_longfile(trigpath, Thresh, M)
%  This function is used to find the triggers within a long trigger file,
%  where multiple stimuli were presented. This function differs from 
%  findtrig.m in that it returns the position of the triple triggers, and
%  does not remove the extra triggers for each triple trigger.
%
%  trigpath : filename 
%
%  Thresh : Threshhold : [-1 , 1]. Default: Thresh==.75
%
%  M : Block Size. Default: M==1024*256. Unless you really care, I don't
%     recommend supplying this input argument. The default works nicely.
%
%  trigger : Returned Trigger Time Vector (in sample number)
%  triptrig : vector of triple trigger times (in sample number)
%
% Craig Atencio
% 6-8-2012
%
% [trigger, triptrig] = findtrig_longfile(trigpath, Thresh, M);

error(nargchk(1,3,3));

%Checking Inputs
if nargin<2
	M = 1024*256;
	Thresh = .75;
elseif ( nargin < 3 )
	M = 1024*256;
end 


% fileinfo = dir('board-ADC-00.dat');
% num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
% fid = fopen('board-ADC-00.dat', 'r');
% v = fread(fid, num_samples, 'uint16');
% fclose(fid);
% v = v * 0.000050354; % convert to volts


%Checking if trigpath is a string
if ( ~ischar(trigpath) )
    error('TRIGPATH must be a string.');
else
    %Remove leading and trailing whitespaces
    trigpath = strtrim(trigpath);
end


% %If trigpath is a .ncs file, converts it to a .raw file
% if ( strcmpi( trigpath((end-2):end) , 'ncs') )
%     ncs2raw(trigpath,[trigpath(1:(end-3)) 'raw'],0);
%     trigpath = [trigpath(1:(end-3)) 'raw'];
% end

fileinfo = dir(trigpath);
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(trigpath, 'r');
X = fread(fid, num_samples, 'uint16');
fclose(fid);
%X = X * 0.000050354; % convert to volts
MeanX = mean(X);
MaxX = max(X);
X = (X-MeanX)/MaxX;
ThreshN = max(X)*Thresh; %max(X)*0.70; % add new threshold 


% % Opening trigger file
% fid = fopen(trigpath,'r');
% 
% % Finding the Aproximate Mean in the File. Opening file and finding the 
% % Mean by Searching 10 Segments
% MeanX=0;
% for j=1:10
%    %Reading Data
%    X = fread(fid,M,'int16')';
%    MeanX = MeanX + mean(X);
% end
% MeanX = MeanX / 10;
% 
% 
% %Finding the Max in the File. 
% frewind(fid);
% 
% MaxX = 0;
% while ~feof(fid)
%    %Reading Data
%    X = fread(fid,M,'int16')';
%    MaxX = max([MaxX abs(X-MeanX)]);
% end
% fclose(fid);


%Closing Files
fclose('all');


% Finding Triggers
trigger = [];
count = 0;

%Opening trigger file
fid = fopen(trigpath,'r');

%Loading data and finding Triggers
while ~feof(fid)

   % Reading Input File - normalize to between [-1, 1]
   X = ( fread(fid, M,'int16')'-MeanX ) / MaxX; % changed to this for snr 10Oct17 Natsumi
% % 	X = ( fread(fid, M,'int16')'-MeanX ) / MaxX;
% 
% 	X = fread(fid, M,'uint16');

%    % Correct for biphasic pulse triggers; zero out the phase of the pulse
%    % that won't be used to find the trigger time
%    if ( Thresh > 0 )
%       X(X<0) = 0;
%    else
%       X(X>0) = 0;
%       X = -X;
%    end



   % Set the triggers to a maximum/minimum value to eliminate the Gibbs
   % phenomenon effect
%    X(X>0.5*MaxX) = 1;
%    X( X >= Thresh ) = 1;
%    X( X < Thresh ) = 0;
   X = X(:)';

% plot(X)
% size(X)
% ylim([0 3.3]);
% pause


   %Adding First Element to avoid missing trigger 
   %if trigger falls exactly on edge
   if ( count == 0 )
      X0 = X(1);
      X = [X0 X];
   elseif ( length(X)>0 )
      X = [X0 X];	
   end


   %Finding Triggers
   if ( length(X) > 0 )

%       %Setting Anything < Tresh to zero
%       index = find( X < abs(Thresh) );
%       X(index) = zeros(1,length(index));

      X( X >= ThreshN ) = 1;
      X( X < ThreshN ) = 0;
%       X( X >= Thresh ) = 1;
%       X( X < Thresh ) = 0;


% clf;
% plot(X);
% ylim([0 3.3]);
% pause

      %Finding Edges in X
      D = diff(X);
      indexD = find( D >= abs(Thresh) );	%Finding Trigger Locations

      %Converting to Trigger Times
      trigger = [trigger M*count+indexD];

      %Finding the Last Element of X. This is placed as the first 
      % element in the new X Array
      X0 = X(length(X));

      %Incrementing Counter
      count = count + 1;
   end
end


%Re-calculating inter-trigger spacing
Dtrigger = diff(trigger);
maxD = max(Dtrigger);
maxD = median(Dtrigger); % just trying to find what is a normal trigger interval

%Finding triple trigger(s)
triple_ind = find(Dtrigger(1:end-1) < 0.34*maxD & Dtrigger(2:end) < 0.34*maxD );
%triple_ind = find(Dtrigger(1:end-2) < 0.34*maxD & Dtrigger(2:end-1) < 0.34*maxD & Dtrigger(3:end) < 0.34*maxD); % chenged natsumi 10Oct17
triptrig = trigger(triple_ind); % this should tell us where the triple triggers are


%Removing additional triggers if extra triple triggers are found
if(isempty(triple_ind))

    disp('Could not find initial triple trigger.');

else %Triple triggers found

   % Multiple triple triggers found. We will assume that the first and 
   % last triple trigger were for the same stimulus. This is the case
   % for the static ripple/moving ripple stimuli, but not for other stimuli
   if diff(triple_ind) > 1 %( length(triple_ind) == 2 )  % changed  this to avoid over trimming when triple_ind = [1,2] Natsumi
   %         trigger = trigger( triple_ind(1):(triple_ind(2)-1) );
     trigger = trigger( triple_ind(1):triple_ind(2) );

   %There was only one triple trigger
   else
     trigger = trigger(triple_ind(1):end);
   end
   
   if ( length(triple_ind) == 2 &&diff(triple_ind) > 1) % add this instead of above lines Natsumi
       trigger = [trigger(1) trigger(5:end)];
   else
       %Making triple trigger count as one trigger
       trigger = [trigger(1) trigger(4:end)];
   end

end

return;



