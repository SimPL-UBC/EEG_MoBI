This repo contains code for analyzing data from the **Mo**bile **B**rain and Body **I**maging (**MoBI**) system developed in the final four months of Thomas Lynch's bachelors degree. 

<p align="center">
  <img src="https://i.imgur.com/Y9gSyUQ.png" width="700">
</p>


We examined cortico-muscular coherence (CMC) during treadmill walking under varying attentional conditions. Specifically, we tracked how CMC differed when participants focused their attention on the active tibia muscle versus when their attention was diverted by a concurrent cognitive task.

You can view the project proposal for more detail here: (LINK)

If there are any questions, please contact Thomas Lynch at lyncht248@gmail.com.


# Data Processing Pipeline
	
	## Import_txt_data_openBCI_SD.m
		Imports EEG, EMG, and IMU data into Matlab from the OBCI_XX.TXT files that are
		recorded on the SD card by the OpenBCI board. 

	## ChopTheData.m
		Divides the EEG, EMG, and IMU data into trials with manually-defined
		boundaries. Trial 0 is the seated isometric contractions, Trial 1 and 3 are 
		the focused walking trials, and Trials 2 and 4 are the walking trials with 
		a math task. 
		Outputs all data into a struct, C_trials.

	## EEGLabScripting.mlx
		Live script with all the steps neccessary to process the EEG and EMG data in
		EEGLAB (band-ass filter, ICA, etc.). Each section requires some manual input
		before running; these spots are marked by '%%%' and ALL-CAPS. 
		Outputs struct C_trials_p, processed EEG and EMG data. 

	## EpochingByFootsteps.mlx
		Live script that uses IMU data to epoch the EEG and EMG data. Each epoch is 
		some number of data points either side of the heel strike from the dominant 
		leg which has the EMG electrodes attached.
		Outputs EMGepochs_TX, a 2D matrix, and EEGepochs_TX, a 3D matrix.

	## CMC.m
		Calculates the coherence matrix of every EEG channel compared to the EMG 
		channel. It combines epochs from Trials 1 and 3 (Undistracted trials) 
		and Trials 2 and 4 (Distracted trials) into two separate coherence 
		matrices, and plots them both. 
		Experimentially, this function also computs the imaginary portion of the
		coherence phasor, and the significance level. These functions haven't been
		fully tested.
		

# Other Files

	## Playground.m
		Contains a number of one-off functions that were (and may be) useful, but 
		aren't explicitly part of the data-processing pipeline. 

	## Dependencies
		This folder contains all the files that CMC.m requires to run. The top-most 
		files is welch.m, which needs to be in the Matlab path for CMC.m to work.

	## OpenBCI_32bit_Library_IMUXX_SRXXX
		These are firmware files used to change the OpenBCI board's sample rate.
		The IMU can sample at either 25Hz or 100Hz, and the electrodes can either
		sample at 250Hz or 500Hz. However, the OpenBCI GUI can only display data
		at 25Hz for the IMU and 250Hz for the electrodes.

© Thomas Lynch
