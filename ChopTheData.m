% PURPOSE:
%   To divide EEGData_XX and IMUData_XX into trials, and output them all
%   into struct C_trials_XX. 
% 
% INPUT:
%   EEGData_XX, the raw EEG Data from participant XX
%   IMUData_XX, the raw IMU data from participant XX
%   tEEG_XX, the time vector from participant XX
%   tIMU_XX, the time vector from participant XX
%
% OUTPUT:
%   C_trials_XX, a struct where: 
%       eeg_TX variables have 17 channels; 16 EEG and 1 time. 
%       imu_TX have 4 channels; 3 IMU and 1 time. 
% 
% Trial 0 is start of the first static contraction to the end of last
%   static contraction.
% Trial 1 and 3 are focused trials (no math task)
% Trial 2 and 4 are distracted trials (math task) 
%
% First plot the EEGData (function below), input trial limits, then run the 
% function in the command window.


function [C_trials] = ChopTheData(EEGData,IMUData,tEEG,tIMU)

%%% RUN FUNCTION BELOW:
% plot(tEEG_XX,EEGData_XX(:,16)); grid minor;

%%% INPUT LIMITS BELOW to the closest second for each trial. 
% Record as XXX.00. Go narrow rather than wide if in doubt (treadmill needs
% to slow down and speed up).

    LowerLim_Trial0 = 38.00;
    UpperLim_Trial0 = 299.00;
    LowerLim_Trial1 = 370.00;
    UpperLim_Trial1 = 552.00;
    LowerLim_Trial2 = 644.00;
    UpperLim_Trial2 = 826.00;
    LowerLim_Trial3 = 945.00;
    UpperLim_Trial3 = 1122.00;
    LowerLim_Trial4 = 1188.00;
    UpperLim_Trial4 = 1370.00;

    ULim_T0_eeg=find(abs(tEEG-UpperLim_Trial0)<0.001);
    LLim_T0_eeg=find(abs(tEEG-LowerLim_Trial0)<0.001);
    ULim_T1_eeg=find(abs(tEEG-UpperLim_Trial1)<0.001);
    LLim_T1_eeg=find(abs(tEEG-LowerLim_Trial1)<0.001);
    ULim_T2_eeg=find(abs(tEEG-UpperLim_Trial2)<0.001);
    LLim_T2_eeg=find(abs(tEEG-LowerLim_Trial2)<0.001);
    ULim_T3_eeg=find(abs(tEEG-UpperLim_Trial3)<0.001);
    LLim_T3_eeg=find(abs(tEEG-LowerLim_Trial3)<0.001);
    ULim_T4_eeg=find(abs(tEEG-UpperLim_Trial4)<0.001);
    LLim_T4_eeg=find(abs(tEEG-LowerLim_Trial4)<0.001);

    ULim_T0_imu=find(abs(tIMU-UpperLim_Trial0)<0.005);
    LLim_T0_imu=find(abs(tIMU-LowerLim_Trial0)<0.005);
    ULim_T1_imu=find(abs(tIMU-UpperLim_Trial1)<0.005);
    LLim_T1_imu=find(abs(tIMU-LowerLim_Trial1)<0.005);
    ULim_T2_imu=find(abs(tIMU-UpperLim_Trial2)<0.005);
    LLim_T2_imu=find(abs(tIMU-LowerLim_Trial2)<0.005);
    ULim_T3_imu=find(abs(tIMU-UpperLim_Trial3)<0.005);
    LLim_T3_imu=find(abs(tIMU-LowerLim_Trial3)<0.005);
    ULim_T4_imu=find(abs(tIMU-UpperLim_Trial4)<0.005);
    LLim_T4_imu=find(abs(tIMU-LowerLim_Trial4)<0.005);

    C_trials.eeg_T0 = cat(2,EEGData(LLim_T0_eeg:ULim_T0_eeg,:),tEEG(LLim_T0_eeg:ULim_T0_eeg,1));
    C_trials.eeg_T1 = cat(2,EEGData(LLim_T1_eeg:ULim_T1_eeg,:),tEEG(LLim_T1_eeg:ULim_T1_eeg,1));
    C_trials.eeg_T2 = cat(2,EEGData(LLim_T2_eeg:ULim_T2_eeg,:),tEEG(LLim_T2_eeg:ULim_T2_eeg,1));
    C_trials.eeg_T3 = cat(2,EEGData(LLim_T3_eeg:ULim_T3_eeg,:),tEEG(LLim_T3_eeg:ULim_T3_eeg,1));
    C_trials.eeg_T4 = cat(2,EEGData(LLim_T4_eeg:ULim_T4_eeg,:),tEEG(LLim_T4_eeg:ULim_T4_eeg,1));

    C_trials.imu_T0 = cat(2,IMUData(LLim_T0_imu:ULim_T0_imu,:),tIMU(LLim_T0_imu:ULim_T0_imu,1));
    C_trials.imu_T1 = cat(2,IMUData(LLim_T1_imu:ULim_T1_imu,:),tIMU(LLim_T1_imu:ULim_T1_imu,1));
    C_trials.imu_T2 = cat(2,IMUData(LLim_T2_imu:ULim_T2_imu,:),tIMU(LLim_T2_imu:ULim_T2_imu,1));
    C_trials.imu_T3 = cat(2,IMUData(LLim_T3_imu:ULim_T3_imu,:),tIMU(LLim_T3_imu:ULim_T3_imu,1));
    C_trials.imu_T4 = cat(2,IMUData(LLim_T4_imu:ULim_T4_imu,:),tIMU(LLim_T4_imu:ULim_T4_imu,1));
end

