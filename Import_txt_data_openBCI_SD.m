% PURPOSE: 
%   Import EEG and EMG data from OBCI_XX.TXT
% INPUT:
%   file_location: string pointing to text file from OpenBCI
% OUTPUT: 
%   (1) workspace variables EEGData_XX, IMUData_XX with all EEG, EMG, and accelerometer data as a
%       numeric matrix
%   (2) tEEG_XX and tIMU_XX, the time vectors that are the same length as EEGData and IMUData
% WARNING: must input required columns on lines 29-31 (approx.)


function [EEGData, IMUData, tEEG, tIMU] = Import_txt_data_openBCI_SD(file_location)

%read file in
fileID = fopen(file_location,'r');
A = textscan(fileID,'%s');
size_A = size(A{1,1});
C = cell(size_A(1),12);
for i = 4:size_A(1)
    D = strsplit(A{1,1}{i,1},',');
    D_size = size(D); 
    for j = 1:D_size(2)
        C{i-3,j} = D{1,j};
    end
end
fclose(fileID);

%Remove the last 100 lines because they are used for SD card documentation
C(end-100:end,:) = [];

%Add time to everything
fs = 500;
t=(0:1/fs:length(C)/fs);
tEEG=t(1:length(C));
tEEG=num2cell(tEEG);
C = [C tEEG'];

% Change this to the channels that contain data!!!
% 8-channel Cyton: 2-9 are eeg or emg, 10-12 are accelerometer, 13 is time
% 16-channel Cyton: 2-17 are eeg or emg, 18-20 are accelerometer, 21 is time
C_eeg = C(:,[2:17,21]);
C_imu = C(:,18:21);

% This section removes empty rows from C_accel, since the  accelerometer only 
% samples at 25Hz

C_imu(~cellfun(@ischar,C_imu(:,1)),:) = [];

% Converting eeg and accel data to decimal, then applying the scale factor
% that is on the website. Both are in volts.
% Note: assumes EMG is on last channel, which was always the case for
% Thomas' testing
C_eeg = [hextodec_matrix(C_eeg(:,1:end-1)) cell2mat(C_eeg(:,end))];
EEG_scale = 4.5 / 24 / (2^23 - 1);
EMG_scale = 4.5 / 1 / (2^23 - 1);
C_eeg(:,1:end-2) = C_eeg(:,1:end-2)*EEG_scale;
C_eeg(:,end-1) = C_eeg(:,end-1)*EMG_scale;

C_imu = [hextodec_matrix(C_imu(:,1:end-1)) cell2mat(C_imu(:,end))];
Accel_scale = 0.002 / 2^4;
C_imu(:,1:end-1) = C_imu(:,1:end-1)*Accel_scale;


% Converting from cell to numeric matrix format, now that only numbers are
% in matrix
EEGData = C_eeg(:,1:end-1);
IMUData = C_imu(:,1:end-1);
tEEG = C_eeg(:,end);
tIMU = C_imu(:,end);
end

function [dec] = hextodec_matrix( matrix )

a = size(matrix);
% matrix = matrix(1:(a(1)-1),:);
dec = zeros(a(1),a(2));
    for i = 1:a(2)
        dec(:,i) = hex2dec(matrix(:,i));
    end

end
