#!/usr/bin/env python

import argparse
import numpy as np
import re
import csv
import os
import sys
import multiprocessing
import math
from math import remainder

def main():

	#text 'A program to create a matrices with mzs as rows and (All)Samples-Scans as columns. Desired use of output is to calculate fragment ion(mz) cooccurance accross samples. See samenamed pdf for full program docs and usage'
	parser=argparse.ArgumentParser()
	#parser.add_argument('echo', help='test arg 1')
	parser.add_argument('--inputTargetList','-i1', help='provide a file containg the paths to target files SampleScans files', required='True' )
	parser.add_argument('--inputSampleList','-i2', help='provide a file containg all the <biologicalSampleNumber>_<scanNumber> ids, one per line', required='True' )
	parser.add_argument('--sizeMzIndexSlice','-i3', help='provide the desired number output lines(rows) to be written in output files. default is 1200', default=2400, type=int )
	parser.add_argument('--batchSize','-i4', help='provide the desired number of jobs to be put in a batch. default is 10', default=10, type=int )
	parser.add_argument('--outputDir', '-o', help='provide a dir (full path) where the outputs will be written, dont add a / at the end', required='True')
	parser.add_argument('--MZLower', help='the lower bound (from here and including) of the mass filter', default=300 , type=int) #mz-range as per machine
	parser.add_argument('--MZUpper', help='the upper bound (up to and including) of the mass filter', default=1500, type=int)
	args = parser.parse_args()	
	#print (args.echo)
	targetList = args.inputTargetList
	sampleList = args.inputSampleList
	chunkSize = args.sizeMzIndexSlice
	batchSize = args.batchSize
	outDir = args.outputDir
	mzLow = args.MZLower
	mzUp = args.MZUpper

	errIn1 = ("%s does not exist, check filename and permissions?" % targetList)
	assert os.path.exists(targetList), errIn1
	errIn2 = ("%s does not exist, check filename and permissions?" % sampleList)
	assert os.path.exists(sampleList), errIn2
	
	errO = ("cannot make %s, check dir path and permissions?" % outDir)
	warn1 = ("%s does exist, samename contents will be overwritten" % outDir)
	if not os.path.exists(outDir):
		print ("%s does not exist, attemting to making it" % outDir)
		os.makedirs(outDir)
		assert os.path.exists(outDir), errO
	else:
		print (warn1)
		
	outDirBinaryVectors = outDir + '/' + 'BinaryVectors'
	warn5 = ("%s does exist, samename contents will be overwritten" % outDirBinaryVectors)
	if not os.path.exists(outDirBinaryVectors):
		print ("%s does not exist, attemting to making it" % outDirBinaryVectors)
		os.makedirs(outDirBinaryVectors)
		assert os.path.exists(outDirBinaryVectors), errO
	else:
		print (warn5)

	outDirEmpties = outDir + '/' + 'EmptyVecs'
	warn2 = ("%s does exist, samename contents will be overwritten" % outDirEmpties)
	if not os.path.exists(outDirEmpties):
		print ("%s does not exist, attemting to making it" % outDirEmpties)
		os.makedirs(outDirEmpties)
		assert os.path.exists(outDirEmpties), errO
	else:
		print (warn2)

	outDirIon = outDir + '/' + 'IonFrequencies'
	warn4 = ("%s does exist, samename contents will be overwritten" % outDirIon)
	if not os.path.exists(outDirIon):
		print ("%s does not exist, attemting to making it" % outDirIon)
		os.makedirs(outDirIon)
		assert os.path.exists(outDirIon), errO
	else:
		print (warn4)
	cwd = os.getcwd()
	#print (cwd) #prints with no fwd slash
	#print sys.argv[0]
	progName = sys.argv[0] 
	#print (progName) #prints full path and progname
	#print (os.path.basename(progName))
	progNBase = (os.path.basename(progName))
	progNBase = re.sub(r".py", "", progNBase)
	#print (progNBase)
	#logName = ('%s.log' % progNBase)
	logName = ('%s/%s.log' % (outDir, progNBase) )
	print (logName)
        #ErrFname = ('%s.Err' % progNBase)
	ErrFname = ('%s/%s.Err' % (outDir, progNBase))
	print (ErrFname)
	LOG = open (logName, "w+") 
	LOG.write('Starting program log for: %s\n' % progNBase) 
	fullcmd = sys.argv #capture the command that was enetered on the terminal
	#logName = ('%s.log' % progNBase)
	LOG.write('Ran the following command: %s\n' % fullcmd) 
	ERR = open (ErrFname, 'w+')
	ERR.write('Starting prog err log for %s\n' % progNBase)
	#lineoad a file containing the path to each sample file, one target file per line, 38,820 files total
	print ('Running program with ChunkSize: %i and BatchSize: %i \n' % (chunkSize, batchSize) )
	LOG.write('Running program with ChunkSize: %i and BatchSize: %i \n' % (chunkSize, batchSize))
	print ('Running program for mzRange: %i to %i \n' % (mzLow, mzUp) )
	LOG.write('Running program for mzRange: %i to %i \n' % (mzLow, mzUp) )

	def loadAllTargets(fileWithPaths):
		f = open (fileWithPaths, "r")
		fContent = f.read().splitlines()
		listContent = []
		for line in fContent:
			#print (li
			 listContent.append(line)	
		f.close()
		#print (len(listContent))
		numOfTargets = (len(listContent))
		LOG.write('Number of Target files read in %i \n' % numOfTargets)
		return listContent

	#load a target file (from targets provided by loadAlltargets)  
	def getListMzSamplePairs(target):
		#print ("begin read of target file: %s" %target)
		t = open (target, "r")
		content_SF = t.readlines()[1:]
		#content_SF = content_SF.splitlines()
		#content_SF = t.read()[1:].splitlines()
		listOfTuples = []
		for threeColInput in content_SF:
			threeColInput = threeColInput.rstrip()
			#print (3colInput)
			mzId, ssId, Intens  = threeColInput.split("\t")
			#print (mzId)
			#mzId = mzSampleIdPair[0]
			#sampleString = mzSampleIdPair[1]
			uCount = mzId.count('_')
			#if (re.search("\_", mzId)):
			if (uCount == 2):
				literalMz, intMass, decMass = mzId.split("_")
				#print (intMass)
				#print (intMass, decMass)
				#reform the float before filtering, maybe something like mzFloat =  float("%s.%s" % (intMass, decMass)) 
				intMass = int(intMass)
		##t.close()
			else:
		#print(mzId)
				literalMz, intMass = mzId.split("_")
				decMass = 0
				intMass = int(intMass)
                                
			if (intMass >= mzLow and intMass <= mzUp):
					#print (intMass, decMass)
					tup = (mzId, ssId)
					listOfTuples.append(tup)	
		t.close()
		print ("finished read of target file: %s" %target)
		#print (listOfTuples)
		#print (len(listOfTuples))
		return listOfTuples
	#loading all possible sample ids into a dictionary whith the sample ids as the keys, indexes as interger values
	def loadAllSampleIDsIntoDict(filepath):
		dictAllSample = {}
		index=0
		with open (filepath) as f:
			for line in f:
				#index = index + 1
				line = line.rstrip('\n')
				dictAllSample[line] = index
				index = index + 1	
		return dictAllSample
			
	#creates a 'long' vector filled with zeros, uses the information in short vector to change some of the zeros in long vector to ones 
	def createVector(shortVector, vecSize):
	#create list of size (size ) and fill with zeros
		longVec = [0] * vecSize
		modified = False
		#go over each element on the list, and save its index
		for value in shortVector:
			if (value >= 0 and value < vecSize): 
				longVec[value] = 1
				modified = True
			else:
				#print ("Error: Abort, index out of bounds\n" )
				ERR.write("createVector Error: Abort, index out of bounds\n" )
				sys.exit("createVector Error: Abort, index out of bounds\n")
				
		return longVec, modified

	#formats the output by introducing some placeholders, modding the id a bit etc... 
	def formatOutput (indexKey, valueArray ):
		indexKey = re.sub(r"mz_", "", indexKey)
		indexKey = re.sub(r"_", "-", indexKey)
		#print (indexKey)
		indexKey = indexKey.zfill(8)
		indexKey = ('_' + indexKey + '_') 
		#print (type(indexKey))
		#tmpArr = valueArray[0:5]
		arrString = ' '.join([str(v) for v in valueArray])
		#print (indexKey, arrString)
		tpedString = (str(indexKey) + ' ' + arrString)
		return (tpedString)

	def makeMzSampleMatrix(): #this function returns 
		allMzs = {}
		targetCounter = 0
		tupleCounter = 0
		for target in listAllTargets:
			#print (targetCounter, len(listAllTargets))
			##print ("begin making mzSampleMAtrix for target file in list
			MzSamplePairsList = getListMzSamplePairs(target) #rerturns a list of tuples. Each tup has two elements
			#print ('size of target:',len(targetTuples))
			#tupleCounter = 0
			for MzSamplePairTup in MzSamplePairsList: #each MzSamplePairTup has two elements and in this loop we separate those elements 
				#print (tmpTuple)		
				#print (tmpTuple[0])#tmpTuple[0] is the mz-id
				mzId = MzSamplePairTup[0]
				sampleScanIdString = MzSamplePairTup[1]
				
				sampleScanIndx = allSamplesNumericDict[sampleScanIdString] #look up the ss string in the numeric index you created earlier (where?) 
				#print (tmpIndex)
				if mzId in allMzs:
					#print (tmpTuple[0])
					allMzs[mzId].append(sampleScanIndx)
				else:
					allMzs[mzId] = []
					allMzs[mzId].append(sampleScanIndx)
					#print (tmpTuple[0])
				#tupleCounter += 1
			targetCounter += 1
		sizeDictAllMzs = len(allMzs)
		AllMzsKeysIndex = list(dict.keys(allMzs)) #where going to use the position of keys in the list as an index, and we look up the mz by listnumber when we go over range in line 229
		print ("Dictionary construction finished. Number of uniq Mzs in dataset ie size of AllMzsDictsizeDictAllMzs: %i \n" % sizeDictAllMzs )
		LOG.write("Dictionary construction finished. Number of uniq Mzs in dataset ie size of AllMzsDictsizeDictAllMzs: %i \n" % sizeDictAllMzs )
		return (allMzs, sizeDictAllMzs, AllMzsKeysIndex ) 

	def reFormatMzSampleMatrix (dictAllMzs, chunkID, leftIndex, rightIndex, head):
		#mzLow = AllMzsKeysIndex[leftIndex] #is a stringid, looks like 'mz_301_421' (no quotes)
		#if (rightIndex >= sizeDictAllMzs):
		#	mzHigh = AllMzsKeysIndex[(sizeDictAllMzs-1)]
		#else:
		#	mzHigh = AllMzsKeysIndex[rightIndex]
			#mzHigh = 'End'
			#is a stringid #fails on last job lookup, right index is outof bounds - put in try -else here
		#print ("reFormatMzSampleMatrix: what is my mz-low lookup: %s \n" % mzLow)
		#outFile = ("%s/%i-BinaryVectorsOutIndexRange%ito%i_MZRange%sto%s.txt" % (outDir, chunkID, leftIndex, rightIndex, mzLow, mzHigh) )
		#outFileIons = ("%s/%i-IonFrequenciesIndexRange%ito%i_MZRange%sto%s.txt" % (outDirIon, chunkID, leftIndex, rightIndex, mzLow, mzHigh) )
		#outfEmpt = ("%s/%i-EmptyVectorsOutIndexRange%ito%i_MZRange%sto%s.txt" % (outDirEmpties, chunkID, leftIndex, rightIndex, mzLow, mzHigh) )
		outFile = ("%s/%i-BinaryVectorsOutIndexRange%ito%i.txt" % (outDirBinaryVectors, chunkID, leftIndex, rightIndex) )
		outFileIons = ("%s/%i-IonFrequenciesIndexRange%ito%i.txt" % (outDirIon, chunkID, leftIndex, rightIndex) )
		outfEmpt = ("%s/%i-EmptyVectorsOutIndexRange%ito%i.txt" % (outDirEmpties, chunkID, leftIndex, rightIndex) )
		FREQ = open (outFileIons, 'w+')
		OUTVEC = open (outFile, 'w+')
		#OUTVEC_BIN = open (outFile, 'wb') #OUTVEC_BIN is binary encoded python outfile
		EMPTID = open (outfEmpt, 'w+')
		##pHolder1 = 'placeHolder1'
		mzIDcol = 'mzIDColumns'
		##pHolder2 = 'placeHolder2'
		##pHolder3 = 'placeHolder3'

		#header for each output file
		SampleWS = ' '.join([str(v) for v in allSamplesNumericDict])
		#print (pHolder1, indexKey, pHolder2, pHolder3, arrString)
		headerString = ( str(mzIDcol) + ' ' + str(SampleWS) )
		#OUTVEC.write('%s\n' % headerString) #uncomment this to write a header to each outfile
		if (head == True):
			headFile = ("%s/HeaderForBinaryVectorOutputs.txt" % outDir )
			HEAD = open (headFile, 'w+')
			HEAD.write('%s\n' % headerString)

		interval = rightIndex - leftIndex
		warn3 = "reFormatMzSampleMatrix: Warning - interval not equal to chunksize for range:%i to %i on jobNum:%i. This should only occurr if sizeDictAllMzs is not properly divisible by chunkSize - if latter is True, this should only happen when the loop reaches the last interval- ie a 'remainder' job. \n" % ( leftIndex, rightIndex, chunkID) 
		if ( interval != chunkSize):
			print (warn3)
		
		CountEmptyVecs = 0
		CountFullVecs = 0
			
		for start in range(leftIndex, rightIndex):
			#print("!!!, " ,start, leftIndex, rightIndex)				
			mz = AllMzsKeysIndex[start]
			#for mz, shortList in DictAllMzs.items():
			shortList = dictAllMzs[mz]
			lenghtShortList = len(shortList)
			FREQ.write('%s\t%i\t%i\t%i\t%i\n' % (mz,lenghtShortList, chunkID, leftIndex, rightIndex))
			tmpFullVec, modified = createVector(shortList,sizeAllSamVect)
			#print(mz, len(tmpFullVec))
			if (modified == True ):
				stringToPrint=formatOutput(mz, tmpFullVec)
				testString = 'I am test string'
				#print ('This %s would have gone to file: %s' %(mz, outFileNoP))
				#OUTVEC.write('%s\n' % (testString)) #uncomment to allow TEST vec print
				OUTVEC.write('%s\n' % (stringToPrint)) #uncomment to allow full vec print
				#OUTVEC_VEC.write('%s\n' % (stringToPrint)) #uncomment to allow full vec print
				CountFullVecs += 1	
			else:
				CountEmptyVecs += 1
				EMPTID.write('%s\t%i\t%i\t%i\n' % (mz, chunkID, leftIndex, rightIndex)) #grab the name of the empty vector here and write it to a file
		sumOfVecCounts = CountFullVecs + CountEmptyVecs
		if (interval != sumOfVecCounts):
			print ('reFormatMzSampleMatrix: VecSum Error. Number of Vectors not adding up to chunksize for chunkID %i, leftIndex %i, rightIndex %i . Was this a reaminder batch? did all the jobs complete?\n' % (chunkID, leftIndex, rightIndex))
			ERR.write('reFormatMzSampleMatrix: VecSum Error. Number of Vectors not adding up to chunksize for chunkID %i, leftIndex %i, rightIndex %i . Was this a reaminder batch? did all the jobs complete?\n' % (chunkID, leftIndex, rightIndex))
		else:
			LOG.write('reFormatMzSampleMatrix: VecSum Correct. Number of Vectors are adding up to interval for chunkID %i, leftIndex %i, rightIndex %i \n' % (chunkID, leftIndex, rightIndex))
			LOG.write('reFormatMzSampleMatrix: jobID %i, leftIndex %i, rightIndex %i completed sucessfully \n' % (chunkID, leftIndex, rightIndex))
		##EMPTID.write("reFormatMzSampleMatrix: Number of empty vectors for %s : %i \n" %(outFile, CountEmptyVecs) ) #this should theoretically be zero - always

	def MakeJobList (chunkSize):
		LOG.write('MakeJobList: Starting creation of job pool for multiprocessing \n')
		rawNumOfChunks = sizeDictAllMzs / chunkSize
		dictRemain = remainder(sizeDictAllMzs,chunkSize)
		ceilNumOfChunks = math.ceil( (sizeDictAllMzs/ chunkSize) )
		dictSizeIfAllChunksWereFull = chunkSize * ceilNumOfChunks
		max_number_processes = 100 #doesn't look like this is doing anything?
		print ('MakeJobList: sizeDictAllMzs / chunkSize = NumOfChunks and remainder(dash-seperated) : %i%i\n' % (rawNumOfChunks,dictRemain))
		print ("MakeJobList: Number of Chunks in Dictionary rounded to ceiling: %i \n" % ceilNumOfChunks)
		LOG.write('MakeJobList: Number of Chunks in dictionary and remainder: %i%i\n' % (rawNumOfChunks,dictRemain))
		LOG.write("MakeJobList: Number of Chunks in Dictionary aka number of jobs needed: %i \n" % ceilNumOfChunks)
		jobList=[] 
		for chunk in range (ceilNumOfChunks):
			if (chunk == 1):
				head = True
			else:
				head = False

			if (chunk != ceilNumOfChunks-1):
				leftIndex = (chunk * chunkSize)
				rightIndex = (chunk+1) * chunkSize
			else:
				leftIndex = dictSizeIfAllChunksWereFull - chunkSize
				rightIndex = int(sizeDictAllMzs)
			print ('MakeJobList: Sending this to jobsList array: ChunkID %i , leftIndex %i , rightIndex %i' % (chunk,leftIndex,rightIndex) )
			LOG.write ('MakeJobList: Sending this to jobsList array: ChunkID %i , leftIndex %i , rightIndex %i' % (chunk,leftIndex,rightIndex) )
			p = multiprocessing.Process(target=reFormatMzSampleMatrix, args=( allMzs, chunk, leftIndex, rightIndex, head))
			jobList.append(p)
		#print (jobList)
		lenJobList = len(jobList)
		print ("MakeJobList: Number of jobs in jobList:%i\n" % lenJobList)
		LOG.write("MakeJobList: Number of jobs in jobList: %i \n" % lenJobList)
			
		if (ceilNumOfChunks != lenJobList):
			ERR.write('MakeJobList: Number of chunks not equal to num of jobs in jobList. Fatal error , exiting program\n')
			sys.exit('MakeJobList: Number of chunks not equal to num of jobs in jobList. Fatal error , exiting program\n')
		return jobList, lenJobList

	def jobLauncher(jobs, lenJobList):
		#take a list of jobs, subselect on <batchSize>, send selection to a cpu core, do this for range <numberOfBatches> 
		#I want 500 jobs per processsor, so a batchSize=500 (2 files open per function),(1024 is limit for number of files open)
		LOG.write('jobLauncher: Starting creation of batches \n')
		print ('jobLauncher: Starting creation of batches \n')
		
		rawNumOfBatches = lenJobList/batchSize
		jobRemainder = int(remainder(lenJobList, batchSize))
		ceilNumOfBatches = math.ceil(lenJobList/batchSize)
		lenJobListTheor = batchSize * ceilNumOfBatches
		
		LOG.write('jobLauncher: batchSize %i \n' % batchSize)
		print ('Number of batches not rounded:%i \n' % rawNumOfBatches )
		LOG.write ('Number of batches not rounded:%i \n' % rawNumOfBatches )
		print ('Number of batches rounded to ceiling:%i \n' % ceilNumOfBatches )
		LOG.write ('Number of batches rounded to ceiling:%i \n' % ceilNumOfBatches )
		print ("Number of jobs to be in the last batch:%i" % jobRemainder)
		LOG.write ("Number of jobs to be in the last batch:%i" % jobRemainder)
 		
		#if (jobRemainder != 0):
		for batchID in range(ceilNumOfBatches): #ceilNumOfBatches is already the ceiling so we don't need to +1
			if (batchID != (ceilNumOfBatches-1)): #if not last value, else we reach the last value of this range (remember range is zero-based, so if we want acces to the 4th batch we stop at i=3) 
				lIndex  = (batchID * batchSize)
				rIndex = (batchID+1) * batchSize
				#print ('jobLauncher: ')
				for j in range(lIndex, rIndex):
					print("Launching job %i of batchID %i " % (j,batchID))
					print("jobindexes", lIndex, rIndex)
					jobs[j].start()
				for j in range(lIndex,rIndex):
					#print("Grouping job %i of batchID %i " % (j,batchID))
				#print("join ", leftIndex, rightIndex)
					jobs[j].join()
			else:
				lIndex  = int(lenJobListTheor-batchSize)
				#print (lIndex)
				rIndex = lenJobList
				#print (rIndex)
				for j in range(lIndex, rIndex):
					print("Launching job %i of batchID %i " % (j,batchID))
					print("jobindexes", lIndex, rIndex)
					jobs[j].start()
				for j in range(lIndex,rIndex):
					#print("Grouping job %i of batchID %i " % (j,batchID))
					#print("join ", leftIndex, rightIndex)
					jobs[j].join()


	###MAIN###

	#Global Variables

	####functions called:
	listAllTargets = loadAllTargets(targetList)
	allSamplesNumericDict = loadAllSampleIDsIntoDict(sampleList)
	sizeAllSamVect = len(allSamplesNumericDict)
	allMzs, sizeDictAllMzs, AllMzsKeysIndex = makeMzSampleMatrix()
	SampleListIdenties = (list(allSamplesNumericDict))
	lenSampleListIdenties = (len(allSamplesNumericDict))
	LOG.write('Number of SampleIDs as columns:%s' % lenSampleListIdenties)
	jobList, lenJobList = MakeJobList(chunkSize)
	jobLauncher(jobList, lenJobList)

	#a job constitutes a writing task - a function call to reFormatMzSampleMatrix with a chunck as the mz range being processed
	# There are multiple jobs - listed in joblist. Joblauncher looks at joblist   
if __name__ == '__main__': #__name__ is python internal, basically if the name of the program is given on the command line , the function main will be executed
	main()

