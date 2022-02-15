% PURPOSE:
%   Calculates the cortico-muscular coherence (CMC) for every EEG electrode 
%   compared with the EMG data at specific intervals along the gait cycle. 
%   Plots this data.
%   Also calculates the significance threshold and the imaginary component
%   of coherence, but these metrics need more testing. 
% INPUT: 
%   EMGepochs_TX, A 2D matrix of EMG epochs for each dominant-leg heel 
%       strike during trial X.
%   EEGepochs_TX, A 3D matrix of EEG epochs for each dominant-leg heel
%       strike during trial X. 
% OUTPUT: 
%   C_trials_p, the processed EMG and EEG data
%   C_trials, the raw Cyton IMU data and the time vector
%
% Before running the code, you need to check areas marked with ALL CAPS and
% three '%', which may require manual editing. 
% Often you need to restart Matlab before running this code. EEGLAB seems
% to mess up the nexttile function

%%% INSERT CORRECT CHANNEL LOCS
%These are correct channel locs for participants 4-12.
chanlocs{1}='Fp1';  
chanlocs{2}='Fp2';
chanlocs{3}='F3';
chanlocs{4}='Fz';
chanlocs{5}='F4';
chanlocs{6}='C3';
chanlocs{7}='Cz';
chanlocs{8}='C4';
chanlocs{9}='P3';
chanlocs{10}='Pz';
chanlocs{11}='P4';
chanlocs{12}='Oz';
chanlocs{13}='T3';
chanlocs{14}='T4';
chanlocs{15}='C1';
chanlocs=string(chanlocs);


%These are the Undistracted trials, U
% EEGepochs_U=cat(1,EEGepochs_T1,EEGepochs_T3);
% EMGepochs_U=cat(1,EMGepochs_T1,EMGepochs_T3);
EEGepochs_U=EEGepochs_T3;
EMGepochs_U=EMGepochs_T3;


%These are the Distracted trials, D
EEGepochs_D=cat(1,EEGepochs_T2,EEGepochs_T4);
EMGepochs_D=cat(1,EMGepochs_T2,EMGepochs_T4);


plot_U = 1; %1 if you want undistracted coherences plotted
plot_D = 1; %1 if you want distracted coherences plotted


[xU,yU,zU]=size(EEGepochs_U); %x = epochs or gait cycles = ~150, y = samples per epoch = 1000, z = electrodes=15
[xD,yD,zD]=size(EEGepochs_D); 


fs=500;
window_samples=125; %Length of one window to calculate coherence over is 125 data points (or 250ms at 500Hz)
window = hamming(window_samples);   %Use Hamming window
% overlap = 100;           % 125 - 100 = 25 points @ 500Hz =~ 50ms increment
noverlap = 0;
frequencies = linspace(1,60,60);    %Calculate CMC at these discrete frequencies. 


%Across an epoch of 1000 points, we want to calculate coherence every 25
%points (every 50ms), with a window of 125 points (250ms). Note that
%windows overlap by 100 data points. Therefore...
numWindows = round((StrideLengthSamples-window_samples)/25) + 1;
%might need to remove +1... in this case, change domain below!


%Initializing coherence matrices
CoherenceMat_U = zeros(numWindows,length(frequencies),15); 
CoherenceMat_D = zeros(numWindows,length(frequencies),15); 
PhaseMat_U = zeros(numWindows,length(frequencies),15); 
PhaseMat_D = zeros(numWindows,length(frequencies),15); 


% %Calculating the average stride length in ms, so we can plot domain as % of stride
% AvgStrideLength_U_odd = mean(cat(1,diff(locs_T1(1:2:end)),diff(locs_T3(1:2:end))));
% AvgStrideLength_U_even = mean(cat(1,diff(locs_T1(2:2:end)),diff(locs_T3(2:2:end))));
AvgStrideLength_U_odd = mean(diff(locs_T3(1:2:end)));
AvgStrideLength_U_even = mean(diff(locs_T3(2:2:end)));

AvgStrideLength_D_odd = mean(cat(1,diff(locs_T2(1:2:end)),diff(locs_T4(1:2:end))));
AvgStrideLength_D_even = mean(cat(1,diff(locs_T2(2:2:end)),diff(locs_T4(2:2:end))));

AvgStridems = ((AvgStrideLength_D_odd + AvgStrideLength_D_even + AvgStrideLength_U_odd + AvgStrideLength_U_even)/4)*1000;


%Starting the tiled layout figures for plotting later
if (plot_U == 1) 
    figure(7)
    tlo1 = tiledlayout(5,3,'TileSpacing','compact','Padding','compact');
end
if (plot_D == 1)
    figure(8)
    tlo2 = tiledlayout(5,3,'TileSpacing','compact','Padding','compact');
end


clear i j k;
for k=1:15
    for  i=0:(numWindows-1) 
     
        %Make a long vector of every 125-point window from every epoch that is offset from
        %the heel strike by a specific amount
     
        temp_EEGwindows_U = reshape(EEGepochs_U(1:end,(1+i*25):(window_samples+i*25),k)',[xU*window_samples,1]);
        temp_EMGwindows_U = reshape(EMGepochs_U(1:end,(1+i*25):(window_samples+i*25))',[xU*window_samples,1]);

        temp_EEGwindows_D = reshape(EEGepochs_D(1:end,(1+i*25):(window_samples+i*25),k)',[xD*window_samples,1]);
        temp_EMGwindows_D = reshape(EMGepochs_D(1:end,(1+i*25):(window_samples+i*25))',[xD*window_samples,1]);

        args = {window,noverlap,frequencies,fs}; 


        %Calculating fxx, fyy, and fxy, the power and auto spectral
        %densities. x = EEG signals, y = EMG signals. 
        [fxx_U,f1_U,fxxc_U] = pwelch(temp_EEGwindows_U,args{:});
         fyy_U = pwelch(temp_EMGwindows_U,args{:});
         fxy_U = cpsd(temp_EEGwindows_U,temp_EMGwindows_U,args{:});

        [fxx_D,f1_D,fxxc_D] = pwelch(temp_EEGwindows_D,args{:});
         fyy_D = pwelch(temp_EMGwindows_D,args{:});
         fxy_D = cpsd(temp_EEGwindows_D,temp_EMGwindows_D,args{:});


         %Getting the coherence (same as mscohere) and the coherency (imaginary
         %and real parts)
         Rxy2_U = (abs(fxy_U).^2)./bsxfun(@times,fxx_U,fyy_U);          %Coherence (coherency squared)
         Rxy_U = fxy_U./(sqrt(bsxfun(@times,fxx_U,fyy_U)));             %Coherency (complex number)
        
         Rxy2_D = (abs(fxy_D).^2)./bsxfun(@times,fxx_D,fyy_D);          %Coherence (coherency squared)
         Rxy_D = fxy_D./(sqrt(bsxfun(@times,fxx_D,fyy_D)));             %Coherency (complex number)


         %Saving to CoherenceMat
         CoherenceMat_U(i+1,:,k) = Rxy2_U;
         PhaseMat_U(i+1,:,k) = angle(Rxy_U);

         CoherenceMat_D(i+1,:,k) = Rxy2_D;
         PhaseMat_D(i+1,:,k) = angle(Rxy_D);


%          %Using the function on line 54-55 of welch.m to get the number of
%          %windows, L
%           args = {window,noverlap,[],fs};
%           [temp,~,~,temp2,~,~,winName,winParam,noverlap,k1,L,options] = signal.internal.spectral.welchparse(temp_EEGwindows_U,'psd',args{:});
%           L = double(k1); % L should be 22, but this is a good just-in-case check
%           
%           %Getting the variance
%           Rxy_Z = atanh(abs(Rxy_U));  %Z-transformed magnitude of coherency
%           Variance_Rxy_Z = 1/(2*L);
% 
%            SignificanceLevelMsCohere = 1 - (0.05/5)^(1/(L-1));
% %           SignificanceLevelMsCohere = 1 - 0.05^(1/(220-1));
% %           SignificanceLevelZ = 1.96/sqrt(2*L);
% %           SignificanceLevelZ = mean(Rxy_Z) + 1.96/sqrt(2*L);
% %           SignificanceLevelIm;
    end
    
    %TODO: This (and all other code for creating figures) should be moved to
    %their own function, not inserted into the CMC.m function.

    if (plot_U == 1) 
        figure(7)
        nexttile(tlo1);
        %The first window spans from -499 to -375 data points offset from
        %the heel strike, so centered about ~-437 or -874ms. 
        %The last window spans from 375 to 500 data points offset from the
        %heel strike, so centered about ~438 or 876ms.
        %We divide by AvgStridems to translate domain to % of gait cycle.
        domain=(linspace(-874/AvgStridems*100,(876/AvgStridems*100),numWindows)); 
        colormap jet;
        contourf(domain,frequencies(10:50),CoherenceMat_U(:,10:50,k)',linspace(0,0.1,20));
        caxis([0 0.05]);
        xlabel(tlo1, '% Stride Length')
        title('' + chanlocs(k) + '')
        title(tlo1, 'Participant 3: Undistracted T3')
        xline(0)
        xlim([-50,50])
        ylabel(tlo1, 'Frequency (Hz)')
    end

    if (plot_D == 1) 
        figure(8)
        nexttile(tlo2);
        domain=(linspace(-874/AvgStridems*100,(876/AvgStridems*100),numWindows)); 
        colormap jet;
        contourf(domain,frequencies(10:50),CoherenceMat_D(:,10:50,k)',linspace(0,0.1,20));
        caxis([0 0.05]);
        xlabel(tlo2, '% Stride Length')
        title('' + chanlocs(k) + '')
        title(tlo2, 'Participant 3: Distracted T2 and T4')
        xline(0)
        xlim([-50,50])
        ylabel(tlo2, 'Frequency (Hz)')
    end


end

if (plot_U == 1) 
    figure(7);
    cb1=colorbar;
    cb1.Layout.Tile = 'east';
end
if (plot_D == 1)
    figure(8);
    cb2=colorbar;
    cb2.Layout.Tile = 'east';
end

%%% SAVE FIGURES AS (un)distracted.fig AND WORKSPACE AS CoherenceMat.mat
%%% PASTE FIGURES INTO FinalResults.pptx



% %This code was used for one long trial (no epochs based on heel strikes). Could
% %be used for T0. 
%
% for k=1:depth 
%         
%     %Calculating fxx, fyy, and fxy
%     [fxx,f1,fxxc] = welch(EEG_epochs(1,:,k),'pwelch',window,overlap,[],fs);
%      fyy = welch(EMG_epochs(1,:),'pwelch',window,overlap,[],fs);
%      fxy = welch({EEG_epochs(1,:,k),EMG_epochs(1,:)},'cpsd',window,overlap,[],fs);
%      
%      %Getting the coherence (same as mscohere) and the coherency (imaginary
%      %and real parts)
%      Rxy2 = (abs(fxy).^2)./bsxfun(@times,fxx,fyy);          %Coherence (coherency squared)
%      Rxy = fxy./(sqrt(bsxfun(@times,fxx,fyy)));             %Coherency
%      
%      %Using the function on line 54-55 of welch.m to get the number of
%      %windows, L
%       args = {window,overlap,frequencies,fs};
%       [temp,~,~,temp2,~,~,winName,winParam,noverlap,k1,L,options] = signal.internal.spectral.welchparse(EEG_epochs(1,:,k),'psd',args{:});
%       L = double(k1);
% 
%       %Getting the variance
%       Rxy_Z = atanh(abs(Rxy));  %Z-transformed magnitude of coherency
%       Variance_Rxy_Z = 1/(2*L);
%       
% %      SignificanceLevelMsCohere = 1 - 0.05^(1/(L-1));
%       SignificanceLevelMsCohere = 1 - 0.05^(1/(220-1));
% %      SignificanceLevelZ = 1.96/sqrt(2*L);
%       SignificanceLevelZ = mean(Rxy_Z) + 1.96/sqrt(2*L);
% %       SignificanceLevelIm;
%       
% %       figure(2);
% %       subplot(2,7,i)
% %       sgtitle('Magnitude squared coherence, abs(Rxy)^2')
% %       plot(f1,Rxy2); 
% %       yline(SignificanceLevelMsCohere);
% %       title(chanlocs{i});
% %       xlim([0 60]);
% % 
%       figure(3);
%       subplot(2,7,k)
%       sgtitle('Z-transformed coherence, arctanh(abs(Rxy))')
%       plot(f1,Rxy_Z); 
%       yline(SignificanceLevelZ);
%       title(chanlocs{k});
%       xlim([0 60]);
% %       
% %       figure(4);      
% %       subplot(2,7,i)
% %       sgtitle('Imaginary part of coherence')
% %       plot(f1,imag(Rxy)); 
% % %     yline(SignificanceLevelMsCohere);
% %       title(chanlocs{i});
% %       xlim([0 60]);
%       
%       
% end
