% PURPOSE: Contains a number of one-off functions that were useful, but not 
%   part of the standard data processing pipeline
       
%% For plotting the adafruitIMU
%change upper and lower limit for the graphing, or just use zoom functions
lowerLim=1;
upperLim=length(adaf_11.imu1);

subplot(6,1,1)
plot(adaf_11.imu1(:,2))
xlim([lowerLim upperLim])

subplot(6,1,2)
plot(adaf_11.imu1(:,3))
xlim([lowerLim upperLim])

subplot(6,1,3)
plot(adaf_11.imu1(:,4))
xlim([lowerLim upperLim])

subplot(6,1,4)
plot(adaf_11.imu1(:,5))
xlim([lowerLim upperLim])

subplot(6,1,5)
plot(adaf_11.imu1(:,6))
xlim([lowerLim upperLim])

subplot(6,1,6)
plot(adaf_11.imu1(:,7))
xlim([lowerLim upperLim])

% resultant = zeros(length(adafruitIMU.imu1),1);
% for i=1:length(adafruitIMU.imu1)
%     resultant(i,1) = sqrt(adafruitIMU.imu1(i,2)^2 + adafruitIMU.imu1(i,3)^2 + adafruitIMU.imu1(i,4)^2);
% end
% figure(2)
% title('Resultant Acceleration, sqrt(x^2 + y^2 + z^2')
% plot(resultant)
% xlim([lowerLim upperLim])
figure(2)
plot(adafruitIMU.imu1(:,3))
xlim([8*10^5 8.2*10^5])

%% Plotting the Cyton IMU data
lowerLim=1;
upperLim=length(IMUData_11);
%lowerLim=3*10^4;
%upperLim=3.02*10^4;

subplot(3,1,1)
plot(tIMU_11,IMUData_11(:,1))
%xlim([lowerLim upperLim])

subplot(3,1,2)
plot(tIMU_11,IMUData_11(:,2))
%xlim([lowerLim upperLim])

subplot(3,1,3)
plot(tIMU_11,IMUData_11(:,3))
%xlim([lowerLim upperLim])



%% Adding time to IMUData and EEGData
%Removing the first 2 lines of IMUData and 
IMUDataTime = zeros(length(IMUData),1);
for i=0:(length(IMUData)-1)
    IMUDataTime(i+1,1)=i/25;
end
EEGDataTime = zeros(length(EEGData),1);
for i=0:(length(EEGData)-1)
    EEGDataTime(i+1,1)=i/500;
end
plot(IMUDataTime,IMUData(:,2),EEGDataTime,EEGData(:,16)/10^6)
xlim([500 501])

%% Chopping a matrix into Trial1, Trial2, etc.
% % Make sure to use combined IMU, EEG matrix.
% 
% function [C_trials] = ChopTheData(C_eeg,C_imu,tEEG_Han,tIMU_Han)
% 
% %Getting the bounds of the trials. Minimum trial length = 170s
% 
%     MinLength = 170;
%     UpperLim_Trial1 = 627;
%     LowerLim_Trial1 = UpperLim_Trial1-MinLength;
%     UpperLim_Trial2 = 908;
%     LowerLim_Trial2 = UpperLim_Trial2-MinLength;
%     UpperLim_Trial3 = 1160;
%     LowerLim_Trial3 = UpperLim_Trial3-MinLength;
%     UpperLim_Trial4 = 1422;
%     LowerLim_Trial4 = UpperLim_Trial4-MinLength;
% 
%     ULim_T1_eeg=(UpperLim_Trial1*500)+1;
%     LLim_T1_eeg=(LowerLim_Trial1*500)+1;
%     ULim_T2_eeg=(UpperLim_Trial2*500)+1;
%     LLim_T2_eeg=(LowerLim_Trial2*500)+1;
%     ULim_T3_eeg=(UpperLim_Trial3*500)+1;
%     LLim_T3_eeg=(LowerLim_Trial3*500)+1;
%     ULim_T4_eeg=(UpperLim_Trial4*500)+1;
%     LLim_T4_eeg=(LowerLim_Trial4*500)+1;
% 
%     ULim_T1_imu=find(tIMU_Han==UpperLim_Trial1);
%     LLim_T1_imu=find(tIMU_Han==LowerLim_Trial1);
%     ULim_T2_imu=find(tIMU_Han==UpperLim_Trial2);
%     LLim_T2_imu=find(tIMU_Han==LowerLim_Trial2);
%     ULim_T3_imu=find(tIMU_Han==UpperLim_Trial3);
%     LLim_T3_imu=find(tIMU_Han==LowerLim_Trial3);
%     ULim_T4_imu=find(tIMU_Han==UpperLim_Trial4);
%     LLim_T4_imu=find(tIMU_Han==LowerLim_Trial4);
% 
%     C_trials.eeg_T1 = C_eeg(LLim_T1_eeg:ULim_T1_eeg,:);
%     C_trials.eeg_T2 = C_eeg(LLim_T2_eeg:ULim_T2_eeg,:);
%     C_trials.eeg_T2 = C_eeg(LLim_T3_eeg:ULim_T3_eeg,:);
%     C_trials.eeg_T2 = C_eeg(LLim_T4_eeg:ULim_T4_eeg,:);
% 
%     C_trials.imu_T2 = C_imu(LLim_T1_imu:ULim_T1_imu,:);
%     C_trials.imu_T2 = C_imu(LLim_T2_imu:ULim_T2_imu,:);
%     C_trials.imu_T2 = C_imu(LLim_T3_imu:ULim_T3_imu,:);
%     C_trials.imu_T2 = C_imu(LLim_T4_imu:ULim_T4_imu,:);
% end

%% plotting a trial of the cyton IMU data
lowerLim=1;
upperLim=length(IMUData_Han);
%lowerLim=3*10^4;
%upperLim=3.02*10^4;

subplot(3,1,1)
plot(tIMU_Han,IMUData_Han(:,1))
%xlim([lowerLim upperLim])

subplot(3,1,2)
plot(tIMU_Han,IMUData_Han(:,2))
%xlim([lowerLim upperLim])

subplot(3,1,3)
plot(tIMU_Han,IMUData_Han(:,3))
%xlim([lowerLim upperLim])

%% Finding local maxima
findpeaks(C_trials.imu_T1(:,2),C_trials.imu_T1(:,4),'MinPeakProminence',0.25,'Annotate','extents')

%% Plotting things
%EMG spike should happen BEFORE the heel strike (during swing), pick odd or
%even
plot(C_trials.eeg_T1(:,17),Trial1_EMG_processed);
xline(locs(1:2:end)); %odd ones
xline(locs(2:2:end),':'); %even ones
xlabel('Time (s)');ylabel('Volts (V)'); title('Processed EMG Plotted With Heel Strikes');

offset = C_trials.eeg_T1(1,17)-1/500; %this is the time that the EMG processed data starts at;
incr = 5;
for i=1:floor(length(locs)/2)
        x(i) = mean(Trial1_EMG_processed((locs(2*i)-offset)*500+1-(incr*2):(locs(2*i+1)-offset)*500+1));
end
mean(x)

for i=1:floor(length(locs)/2)
    x(i) = mean(Trial1_EMG_processed((locs(2*i-1)-offset)*500+1:(locs(2*i)-offset)*500+1));
end
mean(x)

%% printing average around heel strike to determine activation length
offset = C_trials.eeg_T1(1,17)-1/500; %this is the time that the EMG processed data starts at;
for j=1:474
    for i=1:floor(length(locs)/2)
%         test=j-312;
%         test2=locs((i*2)-1)-offset;
%         test3=(locs((i*2)-1)-offset)*500+1;
%         test4=(locs((i*2)-1)-offset)*500+1+(j-312);
         f(i) = Trial1_EMG_processed(round((locs((i*2)-1)-offset)*500+1+(j-312)));
         test(i)=round((locs((i*2)-1)-offset)*500+1+(j-312));
    end
    g(j)=mean(f);
end
domain=(-311:1:162)/500;
plot(domain,g)
title('Averaged EMG over each step')
xlim([-0.7 0.3])
xlabel('Offset from heel strike at 0 ms')
xline(0);

%% Making epochs!
offset = C_trials.eeg_T1(1,17)-1/500; %this is the time that the EMG processed data starts at;

EMG_epochs=zeros(floor(length(locs)/2),400);
EEG_epochs=zeros(floor(length(locs)/2),400,14); %14 is the number of EEG channels
for i=1:floor(length(locs)/2)
    EMG_epochs(i,:) = Trial1_EMG_processed((locs((i*2)-1)-offset)*500+1-299:(locs((i*2)-1)-offset)*500+1+100);
end

for j=1:14
    for i=1:floor(length(locs)/2)
        EEG_epochs(i,:,j) = Trial1_EEG_processed((locs((i*2)-1)-offset)*500+1-299:(locs((i*2)-1)-offset)*500+1+100,j);
        hold on;
        plot(EEG_epochs(i,:,j))
    end
end



