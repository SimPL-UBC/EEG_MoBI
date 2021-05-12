function P300_sound_EEG_sync_loop(DataDir,ScriptDir,savedFolder)
%Read timestamp of P300 sounds (high freq and low freq)
%Save into %P300_sound variable
% DataDir = 'C:\Users\Cidnee\Documents\UBC Masters Program\Data\20200731_Cyton_EyeOpenCloseBlink_AudioStimuli_Walk\text files';
% ScriptDir = 'C:\Users\Cidnee\Documents\UBC Masters Program\Data\20200731_Cyton_EyeOpenCloseBlink_AudioStimuli_Walk\Timestamps';

% Specify and create a new folder where the plots will go, based on the savedFolder path.
newpath = [savedFolder '\synced_data'];
mkdir(newpath)

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(DataDir)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
%picking out all the .mat files
filePattern_data = fullfile(DataDir, '*.txt');
theFiles_data = dir(filePattern_data);

filePattern_stamp = fullfile(ScriptDir, '*.txt');
theFiles_stamp = dir(filePattern_stamp);

theFiles_data = natsortfiles({theFiles_Data.name});
theFiles_stamp = natsortfiles({theFiles_stamp.name});

%this loop reads in each of the files in the listed directory, and applies
%the import and cropping functions, for both the data and script directory
for k = 1 : length(theFiles_data)
  baseFileName_data = theFiles_data{k};
  fullFileName_data = fullfile(DataDir, baseFileName_data);
  fprintf(1, 'Now reading %s\n', fullFileName_data);
  [~,name_data,~] = fileparts(theFiles_data{k});
  
  baseFileName_stamp = theFiles_stamp{k};
  fullFileName_stamp = fullfile(ScriptDir, baseFileName_stamp);
  fprintf(1, 'Now reading %s\n', fullFileName_stamp);
  %[~,name_stamp,~] = fileparts(theFiles_stamp{k});
  
  % Calling Import_txt_data_2, to convert these files from text data to .mat files
  P300_sound_EEG_sync(fullFileName_data,fullFileName_stamp,newpath,name_data);
end

P300_timestamp_filename = 'Test 1 Trial 4 P300_sound.txt';
P300_timestamp_file = fullfile(ScriptDir, P300_timestamp_filename);
P300_sound = readtable(P300_timestamp_file);
[num_row_P300, num_col_P300] = size(P300_sound);

%Read EEG data file
EEG_filename = 'Trial4.txt';
EEG_file_dir = fullfile(DataDir, EEG_filename);
EEG_file = fopen(EEG_file_dir, 'rt');
EEG_data = readtable(EEG_file_dir);
[num_row_EEG, num_col_EEG] = size(EEG_data);

%Create a .txt file to save the synced EEG data
sync_file = fopen( 'EEG_P300_sound_synced.txt', 'wt' );
fprintf(sync_file, 'EEG_sync_P300_stimuli\n');

P300_row = 2;
EEG_row = -5;
while true
  stimuli_tag = ', 0';
  this_line = fgetl(EEG_file);
  if ~ischar(this_line); break; end
  if((EEG_row < 1)||(P300_row >= num_row_P300))
      this_line = append(this_line, stimuli_tag, '\n');
  else          
      if(P300_sound{P300_row,3} <= EEG_data{EEG_row, 13})
          if (P300_sound{P300_row,3}+seconds(1)>= EEG_data{EEG_row, 13})%(P300_sound{P300_row+1,3}>= EEG_data{EEG_row, 13})
            if(strcmp(P300_sound{P300_row,1}, 'High'))
                      stimuli_tag = ', 1';
            elseif (strcmp(P300_sound{P300_row,1}, 'Low'))
                      stimuli_tag = ', -1';
            end            
          else
                P300_row = P300_row + 2;
          end
      end
      
      this_line = append(this_line, stimuli_tag,'\n');
  end
  if(EEG_row < num_row_EEG)
        EEG_row = EEG_row+1;
  end

  fprintf(sync_file, this_line);
  
end
fclose(EEG_file);
fclose(sync_file);