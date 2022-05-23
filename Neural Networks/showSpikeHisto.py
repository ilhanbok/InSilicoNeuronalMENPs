'''
showSpikeHisto.py
Description: Gather and output spike burst data from a .pkl file
             output from NetPyNE
Usage: Fill out the indicated sections and run the script
Author(s): Ilhan Bok
Last Modified: Jan. 19, 2022
'''

import pandas as pd
import matplotlib.pyplot as plt

''' v FILL THIS OUT v '''
obj = pd.read_pickle(r'net_lfp_raster_old.pkl')
'''
(87, 171, 171) = E4
(107, 212, 150) = E3
(107, 171, 214) = E2
(230, 194, 0) = I2
(230, 82, 0) = I3
(230, 150, 0) = I4
'''
Es = [(87, 171, 171), (107, 212, 150), (107, 171, 214)]
Is = [(230, 194, 0), (230, 82, 0), (230, 150, 0)]

numExc = 495-400 + 588-400 + 570-400 # Number of excitatory cells in network
numInh = 495+400 + 588+400 + 570+400 # Number of inhibitory cells in network

spike_ranges = [(36,42),(68,73),(92,105),(121,139)] # Ad-hoc identification of spike ranges from histogram

l = 0.140 # length of recording time in seconds
''' ^ FILL THIS OUT ^ '''

plt.hist(obj['spkTimes'], bins=[x+0.5 for x in list(range(0,150))])
plt.show()

# Spike peaks (ms): 37, 69, 101, 129
colorset = set(map(tuple,obj['spkColors']))
for color in colorset:
    print(tuple(round(i * 255) for i in color))
    
# Filter out the first ten seconds as per our analysis
ostfilt = [i for i in obj['spkTimes'] if i >= 10]
colfilt = obj['spkColors'][-len(ostfilt):]

sum_sr = 0
for t in spike_ranges:
    sum_sr += t[1]-t[0]
sum_sr = sum_sr / 1000

ostfilt_spk = []
colfilt_spk = []

# Get all spike time points within the given spike time ranges
for ran in spike_ranges:
    appost = [v for v in obj['spkTimes'] if ran[0] <= v <= ran[1]]
    ostfilt_spk += appost
    appcol = [i for i,v in enumerate(obj['spkTimes']) if ran[0] <= v <= ran[1]]
    colfilt_spk += appcol

# Get average excitatory, inhibitory, and total frequency based on MANUAL color mappings
total_frequency = len(ostfilt)/l/len(obj['cellGids'])
print('Total Frequency:')
print(str(total_frequency))

tf_spk = len(ostfilt_spk)/sum_sr/len(obj['cellGids'])
print('Total Spike Frequency:')
print(str(tf_spk))

sanicols = []
colors = map(tuple,colfilt)
for color in colors:
    sanicols.append(tuple(round(i * 255) for i in color))
    
sanicols_spk = []
colors_spk = map(tuple,obj['spkColors'])
for color in colors_spk:
    sanicols_spk.append(tuple(round(i * 255) for i in color))
    
# Excitatory
esum = 0
for e in Es:
    esum += sanicols.count(e)

exc_frequency = esum/l/numExc
print('Exc Frequency:')
print(str(exc_frequency))

esum_spk = 0
for e in Es:
    esum_spk += [sanicols_spk[i] for i in colfilt_spk].count(e)

exc_frequency_spk = esum_spk/sum_sr/numExc
print('Exc Spike Frequency:')
print(str(exc_frequency_spk))

# Inhibitory
isum = 0
for i in Is:
    isum += sanicols.count(i)

inh_frequency = isum/l/numInh
print('Inh Frequency:')
print(str(inh_frequency))

isum_spk = 0
for i in Is:
    isum_spk += [sanicols_spk[i] for i in colfilt_spk].count(i)

inh_frequency_spk = isum_spk/sum_sr/numInh
print('Inh Spike Frequency:')
print(str(inh_frequency_spk))

print('All done.')