#!/gscratch/esci/dacb/anacondainstall/bin/python

import re

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
