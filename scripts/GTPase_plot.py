#
# Authors: Cindy Tan and Dhananjay Bhaskar
# Last Modified: 29 Jun, 2017
#

from lxml import etree as et
import numpy as np
import sys
import matplotlib.pyplot as plt
import math

if len(sys.argv) == 1:
    print '\nError: no file path specified\n'
    exit()

tree = et.parse(sys.argv[1])
root = tree.getroot()
timeframe_elements = root.findall(".//time")

# Parsing data (rows = timeframes, columns = cells)
data = 0 
timeframe_num = len(timeframe_elements)
cell_num = len(timeframe_elements[0].findall(".//cell"))

for timeframe in timeframe_elements:
    cell_elements = timeframe.findall(".//cell")
    timeframe_data = np.empty(cell_num)

    for cell in cell_elements:
        timeframe_data[cell_elements.index(cell)] = cell.get('G')
    if timeframe_elements.index(timeframe) == 0:
        data = [timeframe_data]
    else:
        data = np.concatenate((data, [timeframe_data]))

# Plot data
data = np.transpose(data)

# Rows contain cells and columns are time frames
# Cell IDs correspond to the middle of the 20 x 20 monolayer
plt.figure(1)
plt.plot(data[189], label = "Cell ID 189")
plt.plot(data[190], label = "Cell ID 190")
plt.plot(data[209], label = "Cell ID 209")
plt.plot(data[210], label = "Cell ID 210")
plt.gca().legend(loc="upper right")
plt.xlabel('Time')
plt.ylabel('GTPase')
plt.title('Rho GTPase Activity')
plt.ylim([0.5,0.9])
plt.grid(True)
plt.savefig('GTPase_plot.png')

(num_cells, timelimit) = np.shape(data)

# Get rid of DC component
for cell in np.arange(0, num_cells, 1):
	data[cell] = data[cell] - np.mean(data[cell])

fft_res = np.fft.rfft(data[189])

x_range = len(fft_res);

# Frequency plot
plt.figure(2)
frequencies = np.empty(x_range)
for i in range(x_range):

	# to convert fft indices to frequency values:
	# divide by number of indices and multiply by # of samples per second (freq)
	# ODE solver params: 'stepSize':0.01
	# CHASTE Sampling Rate = 200
    # RHS Multiplier (time scale) = 0.25
	frequencies[i] = float(float(i)/timelimit)*(1/(200*0.01));
		
plt.plot(frequencies, np.absolute(fft_res))

# The following are equivalent: np.sqrt(np.multiply(fft_res, np.conj(fft_res))) and np.absolute(fft_res)

plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT of Rho GTPase Activity (Cell ID 189)')
plt.xlim([0.0,0.1])
plt.grid(True)
plt.savefig('FFT.png')

# Find index corresponding to dominant frequency (i.e. frequency with largest amplitude)
indx = np.argmax(np.absolute(fft_res));
print indx

# Print result
freq = float(float(indx)/timelimit)*(1/(200*0.01));
print freq
