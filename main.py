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
    import numpy as np
    b, a = signal.butter(order, 2*cutoff/sampling_frequency, 'low')
    low_pass_filtered_signal = signal.filtfilt(b, a, np.squeeze(original_signal),padlen=3*max(len(a), len(b))-3)
    return low_pass_filtered_signal

def butterworth_high_pass_filter(original_signal,order,cutoff,sampling_frequency, figures=False):
    from scipy import signal
    import numpy as np
    b, a = signal.butter(order, 2*cutoff/sampling_frequency, 'high')
    high_pass_filtered_signal = signal.filtfilt(b, a,np.squeeze(original_signal),padlen=3*max(len(a), len(b))-3)
    return high_pass_filtered_signal

#not tested
def schmidt_spike_removal(original_signal,fs,figures=False):
    from statistics import median
    import numpy as np 
    #Find the window size
    windowsize = round(fs/2)
    #Find any samples outside of a integer number of windows
    trailingsamples = len(original_signal)%windowsize
    #Reshape the signal into a number of windows

    sampleframes=[original_signal[i:(i+windowsize)] for i in range((len(original_signal)-(trailingsamples))//windowsize)]
    #Find the MAAs
    MAAs = [ max([abs(number) for number in individual_frame]) for individual_frame in sampleframes]

    #While there are still samples greater than 3* the median value of the
    #MAAs, then remove those spikes:
    while any(MAA >median(MAAs)*3 for MAA in MAAs):
        #Find the window with the max MAA
        val = max(MAAs)
        window_num=MAAs.index(val)
        if len(window_num)>1:
            window_num = window_num[0]
        #Find the postion of the spike within that window
        abssampleframe=[abs(number) for number in sampleframes[window_num]]
        spike_position=abssampleframe.index(max(abssampleframe))
        if len(spike_position)>1:
            spike_position = spike_position[0]
        #Finding zero crossings (where there may not be actual 0 values, just a change from positive to negative)
        selected_window=np.array(sampleframes[window_num])
        zero_crossings=np.abs(np.diff(np.sign(selected_window)))
        zero_crossings.append(0)
        zero_crossings=(zero_crossings>1).astype(int)
        #Find the start of the spike, finding the last zero crossing before
        #spike position. If that is empty, take the start of the window:
        zero_crossings_list=zero_crossings.tolist()
        start_nonzeros=[i for i, e in enumerate(zero_crossings_list[:spike_position]) if e != 0]
        if not start_nonzeros:
            spike_start=start_nonzeros[-1]
        else:
            spike_start=0

        #Find the end of the spike, finding the first zero crossing after
        #spike position. If that is empty, take the end of the window:
        zero_crossings_list[:spike_position] = 0
        after_nonzeros==[i for i, e in enumerate(zero_crossings_list) if e != 0]
        if not start_nonzeros:
            spike_end=start_nonzeros[-1]
        else:
            spike_end=windowsize

        #Set to Zero
        sampleframes[window_num][spike_start:spike_end] = 0.0001

        #Recaclulate MAAs
        MAAs = [ max([abs(number) for number in individual_frame]) for individual_frame in sampleframes]

    despiked_signal = np.array(sampleframes)

    # Add the trailing samples back to the signal:
    despiked_signal=np.append(despiked_signal,np.array(original_signal[despiked_signal.size:-1]))

    return despiked_signal

	
def Homomorphic_Envelope_with_Hilbert(input_signal, samplingFrequency,lpf_frequency=8,figures=False):
	#8Hz, 1st order, Butterworth LPF
	from scipy import signal
	import numpy as np
	B_low,A_low = signal.butter(1,2*lpf_frequency/samplingFrequency,'low')
	hilbert_coeff=np.array(signal.hilbert(input_signal))
	homomorphic_envelope = np.exp(signal.filtfilt(B_low,A_low,np.log(np.abs(hilbert_coeff))))
	#Remove spurious spikes in first sample:
	homomorphic_envelope[0] = homomorphic_envelope[1]

	return homomorphic_envelope


#not tested
def getHeartRate(audio_data,Fs,figures=False):
	# 25-400Hz 4th order Butterworth band pass
	audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs)
	audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs)
	# Spike removal from the original paper:
	audio_data = schmidt_spike_removal(audio_data,Fs)
	# Find the homomorphic envelope
	homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs)
	# Find the autocorrelation:
	import numpy as np
	h_e = np.array(homomorphic_envelope)
	y=h_e-np.mean(h_e)
	auto_coeff = np.correlate(y, y, mode='full')
	c = auto_coeff/(np.square(np.max(auto_coeff)))
	signal_autocorrelation = c[len(homomorphic_envelope):-1]
	min_index = int(0.5*Fs)
	max_index = int(2*Fs)

	index = np.argmax(signal_autocorrelation[min_index:max_index])
	true_index = index+min_index-1

	heartRate = 60/(true_index/Fs)

	# Find the systolic time interval:

	max_sys_duration = int(round(((60/heartRate)*Fs)/2))
	min_sys_duration = int(round(0.2*Fs))

	pos = np.argmax(signal_autocorrelation[min_sys_duration:max_sys_duration])
	systolicTimeInterval = (min_sys_duration+pos)/Fs;

	return heartRate, systolicTimeInterval


if __name__ == '__main__':
	from scipy.io import loadmat
	mat = loadmat('test.mat')
	audio_data=mat['audio_data']
	audio_Fs=1000
	heartRate, systolicTimeInterval = getHeartRate(audio_data,audio_Fs)
	print(heartRate,systolicTimeInterval)
	#test result here: 60.24096385542169 0.497
	#matlab result: 70.921985815602840 0.285000000000000


