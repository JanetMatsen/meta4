#!/gscratch/esci/dacb/anacondainstall/bin/python

import re
import os  # to make a dir for plots

inputFile = "extract_CIGAR.LakWasMet55_HOW8_2.sorted.sam.dat"
inputFile = "extract_CIGAR.LakWasMet55_HOW8_2.sorted.sam.dat.testing"

# counters
totalReads = 0
unmappedReads = 0

# initialize dictionaries
lengthFreqDict = { 'M' : {}, 'I' : {}, 'D' : {}, 'N' : {}, 'S' : {}, 'H' : {}, 'P' : {}, 'X' : {}, '=' : {} }
elementFreqDict = { 'M' : 0, 'I' : 0, 'D' : 0, 'N' : 0, 'S' : 0, 'H' : 0, 'P' : 0, 'X' : 0, '=' : 0 }

# open the file
with open(inputFile, "r") as inF:
    # extract each line from the file and process
    for line in inF:
        totalReads = totalReads + 1
        # progress bar
        if totalReads % 100000 == 0:
            print("parsed 100000 reads, {0:d} in total".format(totalReads))
        # remove trailing newline
        line = line.rstrip()
        # do something with the line
        # lines that contain a * are unmapped
        if line == "*":
            unmappedReads = unmappedReads + 1
        else:
            # find the matching elements
            elements = re.findall("([0-9]+[MIDNSHPX=])", line)
            # iterate over them
            for element in elements:
                # parse out the element type (e.g. MIDN) and length
                elementType = element[-1]
                elementLength = element[:-1]
                # add one to the length histogram if the type has been seen before or set it to 1 if not
                if elementLength in lengthFreqDict[elementType].keys():
                    lengthFreqDict[elementType][elementLength] += 1
                else:
                    lengthFreqDict[elementType][elementLength] = 1
                # update elementType histogram
                elementFreqDict[elementType] += 1

# output
print("Found {0:d} total reads of which {1:d} were unmapped".format(totalReads, unmappedReads))

# element types per read across all mapped reads
print("Frequency of each element type across all mapped reads:")
for elementType in elementFreqDict.keys():
	print("{0:s}\t{1:d}".format(elementType, elementFreqDict[elementType]))

# per element histogram
print("Frequency of lengths for each element type:")
for elementType in lengthFreqDict.keys():
	print("{0:s}:".format(elementType))
	# return the lengths in the dictionary after mapping each element of keys
	# to an integer and sorting the list
	for lenKey in sorted(map(int, lengthFreqDict[elementType].keys())):
		# print out the line, not the conversion from integer 
		#  which is what lenKey is to string which is what required for dict access
		print("\t{0:d}\t{1:d}".format(lenKey, lengthFreqDict[elementType][str(lenKey)]))

# Bad form to import low in script, so move to top once things work. 
# http://stackoverflow.com/questions/5926061/plot-histogram-in-python
import matplotlib
# Force matplotlib to not use any Xwindows backend.
        # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# make dir for all figures
plot_dir = "./parse_CIGAR_plots/"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir, 0755)

# generate plot for frequency of each CIGAR letter's appearance across reads. 
pos = np.arange(len(elementFreqDict.keys()))
width = 1.0     # gives histogram aspect to the bar diagram
ax = plt.axes()
ax.set_xticks(pos + (width / 2))
ax.set_xticklabels(elementFreqDict.keys())

print(elementFreqDict.keys())
print(elementFreqDict.values())
ax.set_xlabel('CIGAR letter')
ax.set_ylabel('frequency across reads')
ax.set_title('CIGAR fequencies: BWA')

plt.bar(pos, elementFreqDict.values(), width, color='g')
plt.savefig(plot_dir + "elementFreq.pdf", dpi=None, facecolor='w', edgecolor='w',
        transparent=True, bbox_inches=None, pad_inches=0.1,
        frameon=None)

## clear figure between plots
#plt.clf()
#
#pos = np.arange(len(lengthFreqDict['M'].keys()))
#print(lengthFreqDict['M'].keys())
#print(lengthFreqDict['M'].values())
#ax = plt.axes()
#ax.set_xlabel('length for CIGAR letter ZZZ')
#ax.set_ylabel('frequency across reads')
#ax.set_title('length of ZZZ sequence of CIGARS: BWA')
#plt.bar(pos, lengthFreqDict['M'].values(), width, color='b')
#plt.savefig(plot_dir + "lengthFreq.pdf", dpi=None, facecolor='w', edgecolor='w',
#        transparent=True, bbox_inches=None, pad_inches=0.1,
#        frameon=None)


print "use plot_freq_for_cigar_val(key)"
def plot_freq_for_cigar_val(key, plot_dir):
	print "lengthFreqDict[key].items()" 
	print lengthFreqDict[key].items()
	d = {int(k):int(v) for k,v in lengthFreqDict[key].items()}
	print "d:"
	print d
	# clear plotting history
	plt.clf() 
	pos = np.arange(len(d.keys()))
	print "pos:"
	print pos
	# define axes labels
	ax = plt.axes()
	#ax.set_xticks(pos + (width / 2))
	#ax.set_xticklabels(lengthFreqDict[key].keys())
	ax.set_xlabel('lengths for CIGAR letter ' + key)
	ax.set_ylabel('frequency across reads')
	ax.set_title('length for flag ' + key +  ' in CIGARS: BWA')
	# generate plot
	#plt.bar(pos, lengthFreqDict[key].values(), width, color='b')

	print "d.keys():"
	print d.keys()
	min_bin = np.min(d.keys()) # doesn't work if all values are 0? 
	max_bin = np.max(d.keys())
	bins = np.arange(min_bin, max_bin + 1)
	vals = np.zeros(max_bin - min_bin + 1)
	for k,v in d.items():
		vals[k - min_bin] = v
	print bins
	print vals
	plt.bar(bins, vals)
	# plt.hist(lengthFreqDict[key].keys(), weights = lengthFreqDict[key].values())  # do we need something like bins=range(50) from   http://stackoverflow.com/questions/19212508/plotting-a-histogram-from-pre-counted-data-in-matplotlib ? 
	plt.savefig(plot_dir + "lengthFreq-"+key+ ".pdf", dpi=None, facecolor='w', edgecolor='w',
		transparent=True, bbox_inches=None, pad_inches=0.1,
		frameon=None)

plot_freq_for_cigar_val(key="M", plot_dir=plot_dir)
# Didn't work for all zero ones.  Will fix the np.zeroes thing on Monday.
#for key in lengthFreqDict.keys():
#	print key
#	plot_freq_for_cigar_val(key=key, plot_dir=plot_dir)

# temporary fix: avoid keys that have all zero values.
for key in ["S", "M", "I", "H", "D"]:
	print key
	plot_freq_for_cigar_val(key=key, plot_dir=plot_dir)
