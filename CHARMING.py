import sys
from os.path import isfile
from scipy.stats.mstats import gmean
from scipy.stats import pearsonr
import math
import matplotlib.pyplot as plt
import numpy as np
import random


def calculateMinMax(sequence, aaFreqDict, freqDict, mapDict, windowSize):
    """This function calculates %MinMax values for a given sequence, returned as a list of floats"""

    minMaxValues = []
    
    
    for i in range(int(windowSize/2)):
        minMaxValues.append(0)
    
    #Using the specified sliding window size (windowSize/2 - 1 on either side of the central codon), min/max is calculated
    for i in range(len(sequence)-windowSize+1):
        window = sequence[i:i+windowSize] #list of the codons in the current window

        Actual = 0.0     #average of the actual codon frequencies
        Max = 0.0        #average of the min codon frequencies
        Min = 0.0        #average of the max codon frequencies
        Avg = 0.0        #average of the averages of all the frequencies associated with each amino acid

        #Sum the frequencies
        for codon in window:
            frequencies = aaFreqDict[mapDict[codon]] #list of all frequencies associated with the amino acid this codon encodes
            Actual += freqDict[codon]
            Max += max(frequencies)
            Min += min(frequencies)
            Avg += sum(frequencies)/len(frequencies)

        #Divide by the window size to get the averages
        Actual = Actual/windowSize
        Max = Max/windowSize
        Min = Min/windowSize
        Avg = Avg/windowSize

        percentMax = ((Actual-Avg)/(Max-Avg))*100
        percentMin = ((Avg-Actual)/(Avg-Min))*100

        if(percentMax >= 0):
            minMaxValues.append(round(percentMax,2))
        else:
            minMaxValues.append(round(-percentMin,2))

    #fills in values for codons where window size makes min/max unable to be calculated
    if windowSize % 2 == 1:
        for i in range(int(windowSize/2)):
            minMaxValues.append(0)
    else:
        for i in range(int(windowSize/2)-1):
            minMaxValues.append(0)

    return minMaxValues

def generateCodonSeq(stringSeq):
    """Takes in an mRNA sequence as a string and returns it as a list with 
    each element in the list being a specific codon
    """
    codonSeq = [] #list of codons to be returned
    extras = ""
    for line in stringSeq:
        line = line.rstrip()
        string = str(extras) + str(line)
        i=0
        j=3
        while j<=len(string):
            codonSeq.append(string[i:j])
            i+=3
            j+=3
        extras = str(string[i:])
    return codonSeq


def calculateWindowedCAI(CAIFreqDict, codonSeqList, windowSize):
    """Calculate CAI using sliding windows akin to %MinMax"""
    caiValues = []
    for i in range(int(windowSize/2)):
        caiValues.append(0)
    for i in range(len(codonSeqList)):
        if i + windowSize <= len(codonSeqList):
            windowValues = []
            for j in range(i,i+windowSize):
                windowValues.append(CAIFreqDict[codonSeqList[j]])
            caiValues.append(gmean(windowValues))
    if windowSize % 2 == 1:
        for i in range(int(windowSize/2)):
            caiValues.append(0)
    else:
        for i in range(int(windowSize/2)-1):
            caiValues.append(0)
    return caiValues


def initializeHarmonizedSequence(codonSeq, fromFreqDict, toFreqDict, toAAFreqDict, toCAIFreqDict, windowSize, method):
    """do harmonization initialization (Rodriguez Harmonization) by choosing codons
    of the same rank in each organism
    """
    initSeq = []
    for codon in codonSeq:#chose codon of equal rank in other species
        numBelow = 0
        optionList = aaDict[mapDict[codon]]
        for i in optionList:
            if fromFreqDict[i] < fromFreqDict[codon]:
                numBelow += 1
        toFreqOptions = sorted(toAAFreqDict[mapDict[codon]])          
        desiredFreq = toFreqOptions[numBelow]
        flag = 0 #needed to add flag because sometimes different codons have the same frequency
        for i in optionList:
            if toFreqDict[i] == desiredFreq and flag == 0:
                initSeq.append(i)
                flag = 1
    if method == "MM":
        initialValues = calculateMinMax(initSeq,toAAFreqDict,toFreqDict,mapDict, windowSize)
    elif method == "GM":
        initialValues = calculateWindowedCAI(toCAIFreqDict, initSeq, windowSize)
    return initSeq, initialValues
    
def harmonizedSeqErrorCorrection(codonSeq, harmonizedValues, targetValues, toFreqDict, toAAFreqDict, toCAIFreqDict, windowSize, method):
    """The main harmonization function - finds areas where there is a lot of error between the target and the
    current model (harmonized) values, and attempts to make educated codon changes in these regions to decrease the total
    deviation between the harmonized values and the target values.
    """
    iteration = 0
    flag = 0
    while flag < 5:
        totalDist = 0
        for i in range(len(harmonizedValues)):
            totalDist += abs(harmonizedValues[i]-targetValues[i])
        start = totalDist #to check at end of iteration if sequence improved or not
        
        aboveBelow = [] #0 is below, 1 is above, 2 is equal, tracks for regions of 'continuous error'
        for i in range(len(codonSeq)):
            if targetValues[i]>harmonizedValues[i]:
                aboveBelow.append(0)
            elif targetValues[i]<harmonizedValues[i]:
                aboveBelow.append(1)
            else:
                aboveBelow.append(2)
                
        last = -1
        consecutive = 0
        for i in range(len(aboveBelow)):
            if aboveBelow[i] == last:
                consecutive += 1
            else:
                EC_size = 10
                if consecutive >= EC_size and last != 2:
                    codonSeq2 = list(codonSeq) #copy of codon sequence for editing
                    numReplace = int(consecutive/EC_size) #number of positions in continuous region to change
                    replaceStart = i - consecutive 
                    replacePos = []
                    for j in range(numReplace):
                        replacePos.append(replaceStart + int(consecutive/(numReplace+1))*(j+1) - 2 + (iteration%5))
                    for j in replacePos:
                        replaceCodon = codonSeq[j] #codon at position chosen for alteration
                        replaceAA = mapDict[replaceCodon] #aa where codon is being replaced
                        optionalCodons = aaDict[replaceAA] #synonymous codons for that AA
                        replaceFreq = toFreqDict[replaceCodon]
                        currentChoice = None
                        for k in optionalCodons:
                            if aboveBelow[i-1] == 0: #Then the last window was too low and needs to be raised
                                if toFreqDict[k]>replaceFreq: #if optional codon would raise it
                                    if currentChoice == None:
                                        currentChoice = k
                                    elif toFreqDict[k] < toFreqDict[currentChoice]:
                                        currentChoice = k #to change to the next highest synonymous codon relative to 'replaceCodon'
                            else: #if last window was too high, needs to be lowered
                                if toFreqDict[k] < replaceFreq:
                                    if currentChoice == None:
                                        currentChoice = k
                                    elif toFreqDict[k] > toFreqDict[currentChoice]:
                                        currentChoice = k #to change to the next lowest synonymous codon relative to 'replaceCodon'
                        if currentChoice != None:
                            codonSeq2[j] = currentChoice
                            
                        if method == "GM":
                            harmonizedValues2 = calculateWindowedCAI(toCAIFreqDict,codonSeq2,windowSize)

                        elif method == "MM":
                            harmonizedValues2 = calculateMinMax(codonSeq2, toAAFreqDict, toFreqDict, mapDict, windowSize)

                        else:
                            print("Error, method not specified")
                            return None

                        totalDist2 = 0
                        for j in range(len(harmonizedValues2)):
                            totalDist2 += abs(harmonizedValues2[j]-targetValues[j])
                        if totalDist2<totalDist: #If there is improvement, reinitialize everything for next iteration
                            codonSeq = list(codonSeq2)
                            harmonizedValues = list(harmonizedValues2)
                            totalDist = totalDist2
                            flag = 0
                            
                last = aboveBelow[i]#says sign of difference at last codon position for next conescutive error window
                consecutive = 1 #resets the count

        if start == totalDist:#then no improvement was found for this iteration
            flag += 1
        iteration += 1
    return codonSeq, harmonizedValues

def readCUBTable(path):
    """
    This function reads in a CUB table where the format for a line is:
    <Codon> <Frequency>
    repeated for each codon.
    """
    cubDict = dict()
    cubTable = open(path, "r")
    for line in cubTable:
        line = line.replace("U", "T")
        line = line.split()
        if len(line) == 2:
            if line[0] in mapDict:
                cubDict[line[0]] = float(line[1])
            
    if len(cubDict) != len(mapDict):
        print("The following codons are missing/in the incorrect format in your input CUB table:")
        for codon in mapDict:
            if codon not in cubDict:
                print(codon)
        print("Please fix these and rerun the script.")
        sys.exit()
    return cubDict

def generateAAFreqDict(freqDict, aaDict):
    """
    Takes in a dictionary mapping codons to frequencies and amino acids to codons
    and returns a dictionary mapping amino acids to a list of theirrespective 
    codon frequencies
    """
    aaFreqDict = dict()
    
    for aa in aaDict:
        aaFreqDict[aa] = []
        for codon in aaDict[aa]:
            aaFreqDict[aa].append(freqDict[codon])
            
    return aaFreqDict

def generateRandomSeqs(wtSequence, numHarmonize):
    """
    Generate the desired number of random starting seeds for harmonization
    """
    seeds = []

    for i in range(numHarmonize):
        newSeq = []
        for position in range(len(wtSequence)):
            aa = mapDict[wtSequence[position]]
            options = aaDict[aa]
            choice =  options[random.randint(0,len(options)-1)]
            newSeq.append(choice)
        seeds.append(newSeq)
    
    return seeds

def readInputSeq(pathToSeqFile):
    """
    Go to the file and read the input sequence. This function looks for seqeuences
    in FASTA format. If there are multiple sequences in the file then only the first
    sequence is read
    """
    file = list(open(pathToSeqFile))
    stringSeq = ""
    seqName = file[0][1:].strip()
    for lineNum in range(1,len(file)):
        line = file[lineNum]
        if line[0] == ">":
            break
        else:
            stringSeq += line.strip()
    
    stringSeq = stringSeq.upper()
    As = stringSeq.count("A")
    Ts = stringSeq.count("T")
    Cs = stringSeq.count("C")
    Gs = stringSeq.count("G")
    
    if As + Ts + Cs + Gs != len(stringSeq):
        print("ERROR: An unknown character was found in your sequence.")
        print("ERROR: Please rerun the script inputting a sequence only containing A, T, C, and G.")
        sys.exit()
    elif len(stringSeq)%3 != 0:
        print("ERROR: The input sequence is not of a length divisible by 3.")
        print("ERROR: Please rerun the script inputting a sequence with a length that is divisible by 3.")
        sys.exit()
    else:
        return stringSeq, seqName
    
def checkArgs(argList):
    """
    This function checks each input argument to determine whether it is conforms to the 
    necessary format for that argument. If not, the flag is set to 1 and the program exits.
    """
    flag = 0 
    if len(argList) == 7 or len(argList) == 8:
        if not isfile(argList[1]):
            print("ERROR: No file exists at the given path for the sequence file.")
            flag = 1
        try:
            numMutants = int(argList[2])
        except:
            print("ERROR: The argument given for the number of desired output mutants could not be converted to an integer.")
            flag = 1
        try:
            windowSize = int(argList[3])
            if windowSize < 1:
                print("ERROR: The window size must be a positive integer.")
                flag = 1
        except:
            print("ERROR: The argument given for CUB metric's window size could not be converted to an integer.")
            flag = 1
        if not (argList[4] == "MM" or argList[4] == "GM"):
            print("ERROR: The CUB metric argument was not recognized.")
            flag = 1
        if not isfile(argList[6]):
            print("ERROR: No file exists at the given path for the origin CUB table.")
            flag = 1
        if len(argList) == 8:
            if not isfile(argList[7]):
                print("ERROR: No file exists at the given path for the destination CUB table.")
                flag = 1
    else:
        print("ERROR: The code did not receive the correct number of input arguments.")
        flag = 1
        
    if flag == 1:
        print()
        print("The program ecountered the above error(s) when attempting to read the input arguments.")
        print("Please see the README file for correct input argument formats and rerun the script.")
        print()
        sys.exit()

        
            
    

#dictionaries for mapping between codons and amino acids

mapDict = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}

aaDict = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
          'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
          'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
          'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
          'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
          'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
          'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}

# grab nucleotide sequence and CUB table path arguments
args = sys.argv

print()


checkArgs(args)

# arg order script, seqFile, numMutants, windowSize, metric, CUB1, *CUB2
#if harmonizing sequence to native species
if len(args) == 7:
    
    nucSeq, seqName = readInputSeq(args[1])
    numHarmonized = int(args[2])
    windowSize = int(args[3])
    model = args[4]
    outName = args[5]
    freqDict = readCUBTable(args[6])
    aaFreqDict = generateAAFreqDict(freqDict, aaDict)
    codonSeq = generateCodonSeq(nucSeq)
    randomSeeds = generateRandomSeqs(codonSeq,numHarmonized*10)#generate 10 times as many, choose only best results
    outputs = []
    counter = 0
    
    if model == "MM":
        targetVals = calculateMinMax(codonSeq, aaFreqDict, freqDict, mapDict, windowSize)
        for sequence in randomSeeds:
            initVals = calculateMinMax(sequence, aaFreqDict, freqDict, mapDict, windowSize)
            finalSeq, finalVals = harmonizedSeqErrorCorrection(sequence, initVals, targetVals, freqDict, aaFreqDict, None, windowSize, model)
            
            finalDist = 0
            for val in range(len(finalVals)):
                finalDist += abs(finalVals[val] - targetVals[val])
                
            outputs.append((finalDist,finalSeq,finalVals))
            counter += 1
            if counter%10 == 0:
                print("Harmonized mutant #" + str(int(counter/10)) + " complete")
            
    elif model == "GM":
        targetVals = calculateWindowedCAI(freqDict, codonSeq, windowSize)
        for sequence in randomSeeds:
            initVals = calculateWindowedCAI(freqDict, sequence, windowSize)
            finalSeq, finalVals = harmonizedSeqErrorCorrection(sequence, initVals, targetVals, freqDict, None, freqDict, windowSize, model)
            
            finalDist = 0
            for val in range(len(finalVals)):
                finalDist += abs(finalVals[val] - targetVals[val])
                
            outputs.append((finalDist,finalSeq,finalVals))
            counter += 1
            if counter%10 == 0:
                print("Harmonized mutant #" + str(int(counter/10)) + " complete")
        
#when harmonizing from one species for production in another    
elif len(args) == 8:
    
    nucSeq, seqName = readInputSeq(args[1])
    numHarmonized = int(args[2])
    windowSize = int(args[3])
    model = args[4]
    outName = args[5]
    freqDictNative = readCUBTable(args[6])
    freqDictHost = readCUBTable(args[7])
    aaFreqDictNative = generateAAFreqDict(freqDictNative, aaDict)
    aaFreqDictHost = generateAAFreqDict(freqDictHost, aaDict)
    codonSeq = generateCodonSeq(nucSeq)
    randomSeeds = generateRandomSeqs(codonSeq,numHarmonized*10)#generate 10 times as many, choose only best results
    
    #swap first two sequences in randomSeeds with the WT sequence
    randomSeeds[0] = codonSeq
    if len(randomSeeds) > 1:
        randomSeeds[1] = codonSeq
    
    outputs = []
    counter = 0
    
    if model == "MM":
        targetVals = calculateMinMax(codonSeq, aaFreqDictNative, freqDictNative, mapDict, windowSize)
        unharmonizedVals = calculateMinMax(codonSeq, aaFreqDictHost, freqDictHost, mapDict, windowSize)
        for sequence in randomSeeds:
            
            if counter == 0: #use rodriguez harmonization to initialize first sequence
            
                initSeq, initVals = initializeHarmonizedSequence(sequence, freqDictNative, freqDictHost, aaFreqDictHost, None, windowSize, model)
                finalSeq, finalVals = harmonizedSeqErrorCorrection(initSeq, initVals, targetVals, freqDictHost, aaFreqDictHost, None, windowSize, model)
                
            else:
                initVals = calculateMinMax(sequence, aaFreqDictHost, freqDictHost, mapDict, windowSize)
                finalSeq, finalVals = harmonizedSeqErrorCorrection(sequence, initVals, targetVals, freqDictHost, aaFreqDictHost, None, windowSize, model)
            
            finalDist = 0
            for val in range(len(finalVals)):
                finalDist += abs(finalVals[val] - targetVals[val])
                
            outputs.append((finalDist,finalSeq,finalVals))
            counter += 1
            if counter%10 == 0:
                print("Harmonized mutant #" + str(int(counter/10)) + " complete")
        
    elif model == "GM":
        targetVals = calculateWindowedCAI(freqDictNative, codonSeq, windowSize)
        unharmonizedVals = calculateWindowedCAI(freqDictHost, codonSeq, windowSize)
        for sequence in randomSeeds:
            
            if counter == 0: #use Rodriguez harmonization to initialize first sequence
            
                initSeq, initVals = initializeHarmonizedSequence(sequence, freqDictNative, freqDictHost, aaFreqDictHost, freqDictHost, windowSize, model)
                finalSeq, finalVals = harmonizedSeqErrorCorrection(initSeq, initVals, targetVals, freqDictHost, aaFreqDictHost, freqDictHost, windowSize, model)
                
            else:
                initVals = calculateWindowedCAI(freqDictHost, sequence, windowSize)
                finalSeq, finalVals = harmonizedSeqErrorCorrection(sequence, initVals, targetVals, freqDictHost, aaFreqDictHost, freqDictHost, windowSize, model)
                
            finalDist = 0
            for val in range(len(finalVals)):
                finalDist += abs(finalVals[val] - targetVals[val])
                
            outputs.append((finalDist,finalSeq,finalVals))
            counter += 1
            if counter%10 == 0:
                print("Harmonized mutant #" + str(int(counter/10)) + " complete")
            
outputs.sort()
outFile = open(outName + "_" + model + "_harmonized_sequences.txt", "w")
outFile2 = open(outName + "_" + model + "_harmonized_values.txt", "w")

if len(args) == 8:
    outFile2.write("Unharmonized CUB metric values: " + str(unharmonizedVals) + "\n")
outFile2.write("Target CUB metric values: " + str(targetVals) + "\n")

for seqNum in range(numHarmonized):
    
    seqCorr, sig = pearsonr(outputs[seqNum][2], targetVals)
    outStr = ""
    for codon in outputs[seqNum][1]:
        outStr += codon
    gc = (outStr.count("G") + outStr.count("C"))/len(outStr)
    
    outFile.write(">Harmonized Sequence " + str(seqNum + 1) + ", %GC = " + str(round(gc*100,1)) + " Net distance from target values = " + str(round(outputs[seqNum][0],1)) + ", Pearson correlation " + str(round(seqCorr,4)) + " with target sequence.\n")
    outFile.write(outStr + "\n")
    
    outFile2.write("Harmonized Sequence " + str(seqNum + 1) + " CUB metric values: " + str(outputs[seqNum][2]) + "\n")
    
    # create and save a plot for each harmonized sequence
    plt.plot(range(len(targetVals)), targetVals, c='r', label = "WT")
    plt.plot(range(1,len(outputs[seqNum][2]) + 1), outputs[seqNum][2], c='b', label = "Harmonized")
    plt.xlabel("Codon Position")
    if model == "MM":
        plt.ylim(-100,100)
        plt.ylabel("%MinMax Values")
    else:
        plt.ylabel("Geometric mean values")
    plt.legend()
    plt.savefig(outName + "_" + model + "_harmonization_seq_" + str(seqNum + 1) + "_plot", dpi = 200)
    plt.close()
outFile.close()  
outFile2.close()
    


