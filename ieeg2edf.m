function ieeg2edf(dataset, range, varargin)
  %IEEG2EDF  Downloads data from IEEG-Portal to EDF format.
  %
  %  IEEG2EDF(dataset,range) saves the data in the dataset IEEGDATASET
  %  object to edf. the range is a vector of 1x2 with the start and
  %  end-index that should be saved to EDF.
  %
  %  IEEG2EDF(dataset,range,'filePath') same as above but with supplied
  %  EDF path.
  %
  %     example:
  %         >> session = IEEGSession('Study 005', 'username', 'pwdFile');
  %         >> ieeg2edf(session.data, [1 100000]);   
  %
  %
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Copyright 2013 Trustees of the University of Pennsylvania
  % 
  % Licensed under the Apache License, Version 2.0 (the "License");
  % you may not use this file except in compliance with the License.
  % You may obtain a copy of the License at
  % 
  % http://www.apache.org/licenses/LICENSE-2.0
  % 
  % Unless required by applicable law or agreed to in writing, software
  % distributed under the License is distributed on an "AS IS" BASIS,
  % WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  % See the License for the specific language governing permissions and
  % limitations under the License.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  % Check Input
  if nargin > 2
    fileName = varargin{1};
  else
    [fileName, pathName, ~] = uiputfile('*.edf', ...
          'Select a location to save the EDF file','newFile.edf');
    fileName = fullfile(pathName,fileName);    
  end
  
  % pre-alloc array.
  nrChannels = length(dataset.channels);
  metaStruct = struct(...
    'channels',[], ...
    'samplingrate', dataset.channels(1).sampleRate,...
    'numchannels', nrChannels);
  metaStruct.channels = {dataset.channels.label};
  
  BLOCKSIZE = 5000; % Request 10000 samples per iteration (arbitrary)
  nblocks = max([1 ceil(diff(range)./BLOCKSIZE)]);
  
  % Populate array with data.
  firstCall = true;
  startIdx = range(1);
  try
    for iBlock = 1: nblocks
      display(sprintf('Saving block: %i',iBlock));
      idx = startIdx:(startIdx + BLOCKSIZE);
      data = dataset.getvalues(idx, 1:nrChannels);     
      writeEDF(fileName, data', metaStruct, firstCall);
      firstCall = false;
      startIdx = idx(end)+1;
    end
  catch ME
    rethrow(ME)
  end
  
end

function writeEDF(filename, data, header, firstCall)
  %  Parts of this method originates from the Mathworks File
  %  Exchange and is protected under a BSD licence, which is provided
  %  below:
  %
  %   Copyright (c) 2012, fhz
  %   All rights reserved.
  % 
  %   Redistribution and use in source and binary forms, with or without 
  %   modification, are permitted provided that the following conditions are 
  %   met:
  % 
  %     * Redistributions of source code must retain the above copyright 
  %       notice, this list of conditions and the following disclaimer.
  %     * Redistributions in binary form must reproduce the above copyright 
  %       notice, this list of conditions and the following disclaimer in 
  %       the documentation and/or other materials provided with the distribution
  %       
  %   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
  %   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
  %   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
  %   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
  %   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
  %   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
  %   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
  %   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
  %   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
  %   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
  %   POSSIBILITY OF SUCH DAMAGE.

  % The Streambuffer is used to buffer partial seconds between iterations
  % as EDF is stored in 1 second blocks.
  persistent streamBuffer
  
  samplesPerRecord = floor(header.samplingrate);
  
  [nChans,N] = size(data);

  if firstCall
    streamBuffer = []; %#ok<NASGU>
        
    % Create Labels
    labels = char(32*ones(nChans, 16));
    for n=1:nChans
      ln = length(header.channels{n});
      if ln > 16
        fprintf(1, 'Warning: truncating label %s to %s\n', header.channels{n}, header.channels{n}(1:16));
        ln = 16;
      end
      labels(n,1:ln) = header.channels{n}(1:ln);
    end

    % Scale and convert to in16 (data)
    maxV = max(data, [], 2);
    minV = min(data, [], 2);
    if max(maxV) > 32760 || min(minV) < -32760
        warning('Some data cannot be represented as signed 16-bit integers, extreme values are floored');
        data(data > 32760) = 32760; 
        data(data < -32760) = -32760;
    end
    
    % Scale returned by Webservices is 1.
    maxV = int16(repmat(32767,size(data,1),1));
    minV = int16(repmat(-32767,size(data,1),1));
    maxVc = repmat(32767,size(data,1),1);
    minVc = repmat(-32767,size(data,1),1);
    clearvars Scale

    if ~isa(data,'int16')
      data = int16(data);
    end
    
    % Find integer number of Records in data.
    nrRecords = floor(size(data,2)./samplesPerRecord);
    toBuffer = size(data,2)- (nrRecords * samplesPerRecord);
    streamBuffer = data(:, end-toBuffer+1:end);
    data = data(:,1 : (end-toBuffer));
    
    % Create events channel
    if isfield(header,'events')
        minV(end+1,1) = -32768;
        maxV(end+1,1) = 32767;
        minVc(end+1,1) = -32768;
        maxVc(end+1,1) = 32767;
        recordstmp = zeros(1,ceil(N/header.samplingrate));
        eventchannel =repmat(uint8(00),120,ceil(N/header.samplingrate));
        for i = 1:ceil(N/header.samplingrate)
            tmp = unicode2native(num2str(i-1));
            tmp = [43 tmp 20 20 00]; %#ok<AGROW>
            eventchannel(1:length(tmp),i) = tmp';
            recordstmp(1,i) = length(tmp);
        end
        clearvars tmp
        for i = 1:size(header.events.POS,2);
            tmpPOS = unicode2native(num2str(round(1000* header.events.POS(1,i) / header.samplingrate) / 1000));
            tmpDUR = unicode2native(num2str(round(1000* header.events.DUR(1,i) / header.samplingrate) / 1000));
            tmpTYP = unicode2native(header.events.TYP{1,i});
            eventchanneltmp = [43 tmpPOS 21 tmpDUR 20 tmpTYP 20 00];
            j = ceil(header.events.POS(1,i) / header.samplingrate);
            eventchannel(recordstmp(1,j)+1:recordstmp(1,j) + length(eventchanneltmp),j) = eventchanneltmp';
            recordstmp(1,j) =  recordstmp(1,j) + length(eventchanneltmp);
        end
        clearvars tmpPOS tmpDUR tmpTYP eventchanneltmp
        lengthevents = size(eventchannel,1) / 2;
        if ceil(size(eventchannel,1) / 2) > (size(eventchannel,1) / 2)
            eventchannel = cat(1,eventchannel,repmat(uint8(00),1,ceil(N/header.samplingrate)));
        end
        eventchannel = reshape(eventchannel,1,size(eventchannel,1)*size(eventchannel,2));
        eventchannel = typecast(eventchannel','int16')';
        eventchannel = reshape(eventchannel,lengthevents,length(eventchannel)/lengthevents);
        data = reshape(data,size(data,1),samplesPerRecord,nrRecords);
        data = permute(data,[2 1 3]);
        data = reshape(data,size(data,1) * size(data,2),size(data,3));
        data = cat(1,data,eventchannel);
    else
        data = reshape(data, size(data, 1), samplesPerRecord, nrRecords);
        data = permute(data,[2 1 3]);
        data = reshape(data,size(data,1) * size(data,2),size(data,3));
    end

    %Convert digMin digMax physMin physMax to char-array
    digMin = sprintf('%-8i', minVc);
    digMax = sprintf('%-8i', maxVc);
    physMin = sprintf('%-8i', minV);
    physMax = sprintf('%-8i', maxV);

    % Control EEG recording timestamp
    header.year   = 1970;
    header.month  = 1;
    header.day    = 1;
    header.hour   = 0;
    header.minute = 0;
    header.year   = 1970;
    header.month  = 1;
    header.day    = 1;

    [~,monthstr] = month([num2str(header.day,'%02i') '-' num2str(header.month,'%02i') '-' num2str(header.year)],'dd-mm-yyyy');
    clearvars tmp
    monthstr = upper(monthstr);
    if header.year > 2000
        yearshort = header.year - 2000;
    elseif header.year > 1900
        yearshort = header.year - 1900;
    elseif header.year > 100
        yearshort = header.year - (100*floor(header.year/100));
    else
        yearshort = header.year;
    end

    % Control Recording and Subject info
    header.ID = 'X';
    header.technician = 'X';
    header.equipment = 'X';
    header.subject = [];
    header.subject.ID = 'X';
    header.subject.sex = 'X';
    header.subject.name = 'X';
    header.subject.birthdate = 'X';
    header.subject.ID = regexprep(header.subject.ID,' ','_');
    header.subject.name = regexprep(header.subject.name,' ','_');
    header.technician = regexprep(header.technician,' ','_');
    header.equipment = regexprep(header.equipment,' ','_');

    % Create physdim-info (uV for all channels)
    physdim = char([117 86 32 32 32 32 32 32]);
    physdim = repmat(physdim, nChans, 1);

    % Write edf
    fid = fopen(filename, 'wb', 'ieee-le');
    fprintf(fid, '0       ');   % version
    fprintf(fid, '%-80s', ...
      [header.subject.ID ' ' ...
      header.subject.sex ' ' ...
      header.subject.birthdate ' ' ...
      header.subject.name]);

    fprintf(fid,'%-80s', ...
      ['Startdate ' num2str(header.day,'%02i') '-' monthstr '-' ...
      num2str(header.year,'%04i') ' ' header.ID ' ' ...
      header.technician ' ' header.equipment]);

    fprintf(fid, '%02i.%02i.%02i', header.day, header.month, yearshort); % date as dd.mm.yy
    fprintf(fid, '%02i.%02i.%02i', header.hour, header.minute, 0); % time as hh.mm.ss

    fprintf(fid, '%-8i', 256*(1+nChans));  % number of bytes in header
    if isfield(header,'events')
        fprintf(fid, '%-44s', 'EDF+C'); % reserved (44 spaces)
    else
        fprintf(fid, '%44s', ' '); % reserved (44 spaces)
    end
    fprintf(fid, '%-8i', ceil(N/header.samplingrate));  % number of EEG records
    fprintf(fid, '%8f', 1/header.samplingrate * samplesPerRecord);  % duration of EEG record (=Fs)
    fprintf(fid, '%-4i', nChans);  % number of signals = channels

    fwrite(fid, labels', 'char*1'); % labels
    fwrite(fid, 32*ones(80,nChans), 'uint8'); % transducer type (all spaces)
    fwrite(fid, physdim', 'char*1'); % phys dimension (all spaces)
    fwrite(fid, physMin', 'char*1'); % physical minimum
    fwrite(fid, physMax', 'char*1'); % physical maximum
    fwrite(fid, digMin', 'char*1'); % digital minimum
    fwrite(fid, digMax', 'char*1'); % digital maximum
    fwrite(fid, 32*ones(80,nChans), 'uint8'); % prefiltering (all spaces)
    if isfield(header,'events')
        for k=1:(nChans-1)
            fprintf(fid, '%-8i', samplesPerRecord); % samples per record 
        end
        fprintf(fid, '%-8i', lengthevents); % samples for annotations per record
    else
        for k=1:nChans
            fprintf(fid, '%-8i', samplesPerRecord); % samples per record 
        end
    end
    fwrite(fid, 32*ones(32,nChans), 'uint8'); % reserverd (32 spaces / channel)
    fwrite(fid, data, 'int16');
    fclose(fid);
  else
    
    % Reformat data
    if ~isa(data,'int16')
      data = int16(data);
    end
        
    % concatenate Buffer from previous iteration with new data
    data  = [streamBuffer data];
    
    % Find integer number of Records in data.
    nrRecords = floor(size(data,2)./samplesPerRecord);
    toBuffer = size(data,2)- (nrRecords * samplesPerRecord);
    streamBuffer = data(:, end-toBuffer+1:end);
    data = data(:,1 : (end-toBuffer));
    
    data = reshape(data,size(data,1),samplesPerRecord, nrRecords);
    data = permute(data,[2 1 3]);
    data = reshape(data,size(data,1) * size(data,2),size(data,3));
    
    % Not the first call, assert file exists and append data to file.
    fid = fopen(filename, 'r+', 'ieee-le');
    
    %Adding data to end of file    
    fseek(fid,0,'eof');
    fwrite(fid, data, 'int16');
    
    %Updating header with new number of records in file.
    fseek(fid,236,-1);
    curNrRecords = str2double(fread(fid, 8,'*char'));
    newNrRecords = curNrRecords + nrRecords;
    newNrRecordsStr = sprintf('%-8i',newNrRecords);
    fseek(fid,236,-1);
    fwrite(fid, newNrRecordsStr,'char*1');
    
    fclose(fid);
  end
  

end