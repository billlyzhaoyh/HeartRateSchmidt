#implementation to get heart rate from PCG

# Get heatrate:
# From Schmidt:
# "The duration of the heart cycle is estimated as the time from lag zero
# to the highest peaks between 500 and 2000 ms in the resulting
# autocorrelation"
# This is performed after filtering and spike removal:

#function [heartRate systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)

def butterworth_low_pass_filter(original_signal,order,cutoff,sampling_frequency, figure=False):
	from scipy import signal
	b, a = signal.butter(order, 2*cutoff/sampling_frequency, 'low', analog=True)
	low_pass_filtered_signal = signal.filtfilt(b, a, original_signal)
	return low_pass_filtered_signal

def butterworth_high_pass_filter(original_signal,order,cutoff,sampling_frequency, figures=False):
	from scipy import signal
	b, a = signal.butter(order, 2*cutoff/sampling_frequency, 'high', analog=True)
	high_pass_filtered_signal = signal.filtfilt(b, a, original_signal)
	return high_pass_filtered_signal

def schmidt_spike_removal(original_signal, fs):
	from statistics import median

	#Find the window size
	windowsize = round(fs/2)
	#Find any samples outside of a integer number of windows
	trailingsamples = len(original_signal)%windowsize
	#Reshape the signal into a number of windows
	sampleframes=[original_signal[i:i + windowsize] for i in range((len(original_signal)-trailingsamples)/windowsize)]
	#Find the MAAs
	MAAs = [ max([abs(number) for number in individual_frame]) for individual_frame in superframes]



	#While there are still samples greater than 3* the median value of the
	#MAAs, then remove those spikes:
	while any(MAA >median(MAAs)*3 for MAA in MAAs):
		#Find the window with the max MAA
		

	while(~isempty(find((MAAs>median(MAAs)*3))))
    
    	%Find the window with the max MAA:
    	[val window_num] = max(MAAs);
    	if(numel(window_num)>1)
        	window_num = window_num(1);
    	
    
    	%Find the postion of the spike within that window:
    	[val spike_position] = max(abs(sampleframes(:,window_num)));
    
    	if(numel(spike_position)>1)
        	spike_position = spike_position(1);
   
    
    
     	#Finding zero crossings (where there may not be actual 0 values, just a change from positive to negative):
    zero_crossings = [abs(diff(sign(sampleframes(:,window_num))))>1; 0];
    
    	#Find the start of the spike, finding the last zero crossing before
    	#spike position. If that is empty, take the start of the window:
    spike_start = max([1 find(zero_crossings(1:spike_position),1,'last')]);
    
    	#Find the end of the spike, finding the first zero crossing after
    	#spike position. If that is empty, take the end of the window:
    zero_crossings(1:spike_position) = 0;
    spike_end = min([(find(zero_crossings,1,'first')) windowsize]);
    
    	#Set to Zero
    sampleframes(spike_start:spike_end,window_num) = 0.0001;

    	#Recaclulate MAAs
    MAAs = max(abs(sampleframes));



	return despiked_signal

def getHeartRate(audio_data,Fs,figures=False):
	# 25-400Hz 4th order Butterworth band pass
	audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs)
	audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs)
	# Spike removal from the original paper:
	audio_data = schmidt_spike_removal(audio_data,Fs)
	return heartRate, systolicTimeInterval
