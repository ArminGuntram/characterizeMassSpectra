# Synopsis
A program to create a matrix(s) with mzs' as rows and SampleScans as columns. It can be used to filter mzs' by mass range and/or mz resolution (defined by numbers after the decimal, typically 3rd or 4th depending on the instrumentation). By definition all samples are considered. The desired output are matrices that can be used as inputs to CoMet workflows on Summit. In particular, the output matrices from this program are to be used as input for mz occurrence calculations on Summit. This program Records and writes out ion frequencies for entire dataset. This program forms part of a larger MS2 analysis works. 

# Quick run
```
buildMatrix-FragIonPresenceAbsence.py --inputTargetList <pathsTargetFilesSampleScans> --inputSampleList <FileSampleNameScanNumList> --outputDir <PathtoOutDir>
```
## Other options:
--sizeMzIndexSlice provide the desired number output lines(rows) to be written in output files. default is 1200, default=2400  
--batchSize provide the desired number of jobs to be put in a batch. default is 10, default=10  
--MZLower the lower bound (from here and including) of the mass filter, default=300  
--MZUpper the upper bound (up to and including) of the mass filter, default=1500  

# Description
A program to create a binary matrix(s) with mzs as rows and SampleScans as columns. Each value (1 or 0) representing the presence/absence of an mz for every biological SampleScan combination in the input data. The program can filter mzs' by mass range. The desired output are matrices that can be used as inputs to CoMet workflows on Summit. In particular, the output matrices from this program are to be used as input for mz occurrence calculations on Summit. 

This program also filters out empty (all 0) vectors from the desired outputs. The details of emty vectors (if any) are recorded and written to 'emptyVecs' output.  Ion frequencies for entire dataset as well as other summary stats are also part of this program's output. This program forms part of a larger MS2 analysis works.

This program takes the output from a shotgun-proteomics experiment as it's input. The desired input is reformatted information acquired from a peaklist. In this document this reformatted information is referred to a as a SampleScan file(s).
Each line of such an input contains an ion (given as an mz), by definition if it is in the file it is 'present' ie '1' in our future binary matrix. On each line, the SampleScan number is also given. This is the origin of the fragment ion. 

The program also requires a non-redundant list of SampleScan ids (for use as an index). The program has an err file and log. Some checkpoints are also written to std out. 

The intended use case had 1.8 million scans from 36 biological samples. Ion accuracy was taken at 3 decimals. (Determined at step prior to this program).
The 2.3 million SampleScans are the columns of the program output. Their order is constant. The number of rows (mzs) was determined at 2.3 million when not filtered by mz range. The detectable instrument range is (300mz-1500mz) although one may want to cnsidered ion below or above that. Note that there is a 1000x multiplier here since we are looking at 3 decimal accuracy. 

The output matrices were written in files of approximately 10GB. Each output matrix has the full compliment of SampleScan ids with same order. The number of rows is determined by 'chunkSize' . This user-defined parameter controls how many mz's (number of rows) will be written to every file. (See exemptions for empty vectors.)

# Inputs
Program has two types of input, namely a list of targets and a list of "SampleScan" names.

## --inputTargetList
I1, A text file containing the paths to SampleScans files, has header. 
## --inputSampleList
I2, file containing all the Sample ids. This can be produced on the cmd.

# Multiprocessing -quick overview
Optional command-line argument --batchSize (default=500) controls how the work will be distributed. Maximum 'batchSize' is 500 (default best setting for Merlot). A decrease in 'batchSize' will increase 'numberOfBatches' that need to be run. 'numberOfBatches' = Number of cpu's that will be used. Future add-on will take numCPUs and solve batchsize.

# Outputs
For details of output format or to change the output format see function formatOutput.

# Performance
On Andes the 3 decimal mzs ~1.5 million scans took 12 hours. The same 2 decimal 10 hour was enough although one can probably use spreadsheet to optimize.

# Language 
Python 3.7.2

# Dependencies
argparse; numpy as np; re; csv; os; multiprocessing

# Verbose Notes,Comments for important functions

## loadAllTargets
Takes a file with the locations of multiple target files. Returns list of targets.
## makeMzSampleMatrix
Is the function that performs the main aim of the program.
Function is called without any arguments, instead it uses/needs global variables, namely 'listAllTargets' and 'allSamplesNumericDict'. The main purposes of this function is to create the dictionary of all the mzs that are encountered in the dataset (namely all targets) and importantly in what samples (from SampleList) they were encountered. Because we have such a large amount of mzs and such a large amount of samplescans to evaluate if a given mz was in fact in that sample we need to perform this evaluation as we construct the dictionary (in other words, on-the-fly).

To facilitate this we need a 3 dimensional dictionary. Another way to phrase is a dictionary where the value of every key is actually a list.
As will become clear later we needed the list because retaining order was important(hindsight not sure if it was necessary here). Could we have used a tuple here? or would immutability be a problem? So, the key->value pair would be 'mz':'ListofSampleScanIDs' where the list is constantly being appended to.

We could leave the dictionary in this form, however to increase efficiency we created a numeric index value for every mz using an earlier function(see \ref{SampleDict} ). Thus we have our dictionary 'allSamplesNumericDict' available as global variable. So all that needs to happen in 'makeMzSampleMatrix' is a simple dictionary look-up of every SampleScanID that we encounter in 'target' given in line201: ``sampleScanIndx = allSamplesNumericDict[sampleScanIdString]``

Lines 203 through 207 evaluates and populates the dictionary (Note the simplistic syntax!).
mz entry if not already there, appends the corresponding sampleScanIndx to the list. 'allMzs'. Note that if mz is not yet defined , it and and empty list is added. Note that you definitely do need the 'append' statement in the 'else' portion of the loop. this was a major bug fix. If you don't have it here you miss the 1st instance of every mz!  

It returns three new global variables 'allMzs', 'sizeDictAllMzs'and 'AllMzsKeysIndex'. 'makeMzSampleMatrix' does call function 'getListMzSamplePairs'. It writes size of 'allMzs' (number of keys (MZs) in the dictionary) to LOG. This is a non-redundant count values of all the uniq mzs encountered in the dataset 

'targetCounter' and 'tupleCounter' are debugging count variables. 'targetCounter' keeps track of the number of target input files that have been read in. 'tupleCounter' keeps track of how many 'MzSamplePairTup' have been processed for each target file.

Clarification note: There are two different dictionaries used in this function. The most prominent is 'allMzs' since the main purpose if this function is to construct it. The second is 'allSamplesNumericDict' which we only need to look-up the corresponding index value of 'sampleScanIdString'

## loadAllSampleIDsIntoDict
This function creates an empty dictionary, opens a file and reads it line by line. for every line an index is created, and inserted as an entry into the dictionary. The key:value pairs are sampleScanID:numericIndex. The index is 1-based. The function returns this dictionary 'dictAllSample'. Function 'makeMzSampleMatrix' needs this dictionary and currently gets it as global variable 'allSamplesNumericDict'. Alternatively, 'loadAllSampleIDsIntoDict' may be called from within 'makeMzSampleMatrix' since this is the only place where 'allSamplesNumericDict' is used. For instance calling it at line 157 - before looping over 'listAllTargets'.

## getListMzSamplePairs

This function is called inside function 'makeMzSampleMatrix'. It receives variable named 'target' during a for loop for 'target in listAllTargets'. target is a filename for a the 3-col target file which this function parses. from the 1st col it gets the mz and sample id, then separates out the mz part splits that into the integer and decimal components. 'literalMZ' is just the literal 'mz prefix'. We keep the interger part and make sure it is correctly formatted as such because we are using it as part of our conditional mass filter.

## reFormatMzSampleMatrix
This function is a main program function and is the subject of the multiprocessing.
The 1st ~10 lines are mostly for creating uniq outfile names and making Filehandles. \texttt{mzLow = AllMzsKeysIndex[leftIndex]; mzHigh = AllMzsKeysIndex[rightIndex]} are lookups to find out what the actual mz ids are for the batch about to be processed by the the function 'reFormatMzSampleMatrix'.
'reFormatMzSampleMatrix' calls two other functions, namely 'createVector' and 'formatOutput'. 
'reFormatMzSampleMatrix' Requires/expects dict 'allMzs' as 1st arg. 2nd arg is a 'chunkID' and 3rd and 4th args are the bounds of a range that will be used to access the index of Dict 'AllMzs'. This range of numbers represents a slice of mzs in the dictionary - so literally this range corresponds to the number of mz vectors that form part of one job. 'chunkSize' literally controls how many mz's we are going to pass to 'createVector' then 'formatOutput' and consequently write to its own file.

Importantly, after 'createVector' has been called and returned it's outputs - namely 'tmpFullVec, modified' line (210) - we evaluate 'modified'. If 'True' we pass it to 'formatOutput' else we increment the empty vector count 'CountEmptyVecs' and write the mz-ID to a corresponding outputfile. Note that an empty vector should only exist if that mz is not present in ANY of the input files. This can only occur if we applied a Mass filter in 'getListMzSamplePairs'. an mz vector should only be empty if that mz originated from a sample scan(s) that had no ions in the selected range.

This function was kept as simple as possible but to understand it we must consider the context of the program. This function assumes that 'dictAllMZs' exists and is complete. we also have a variable called 'sizeDictAllMzs'. This is an integer value for the number of keys we have in 'dictAllMZs', ie a count starting from 1 to n. Thus n is literally the number of unique mz's we have. An important thing to remember here is that 'sizeDictAllMzs' is a 1-based count and the indices we are working with are all zero-based. Similarly, all of our 'range' based for loops are also zero-based and the range operator starts and includes 'start' co-ordinate and goes up-to but not including the 'end' co-ordinate.

### err prevention
'reFormatMzSampleMatrix' gets values from a dictionary for a interval provided to it by the two index args. It follows logically that for our intervals to be used, they must exist within the range defined by 'sizeDictAllMzs'. An important thing to remember here is that 'sizeDictAllMzs' is a 1-based count and the indices we are working with are all zero-based. Similarly, all of our 'range' based for loops are also zero-based and the range operator starts and includes 'start' co-ordinate and goes up-to but not including the 'end' co-ordinate. So if we wanted to loop over every index value of 'dictAllMZs' we can simply write: \texttt{for i in range(sizeDictAllMzs):} since we start at zero and don't include the last value we don't need to make any increment/decrement modifications. 
From this we know that the largest integer value rightIndex many be is 'sizeDictAllMzs'-1.
By definition, we don't have any data beyond this index so when the program tries to look-up a value \texttt{line 220 mz = AllMzsKeysIndex[start]} using something that is larger will result in an out-of-bounds error. 

This problem is an artifact that results from the fact that we must work in multiples of 'chunkSize' to establish our intervals.  If everything goes correctly, 'reFormatMzSampleMatrix' should never receive an out-of-bound index since 'makeJobList' is responsible for the calculation of the intervals. However, if there is a problem the condition in line 219 \texttt{if (rightIndex < sizeDictAllMzs):} will prevent the execution of 'createVector', print the details of the offending function call and args to std out and err log.

'warn3': At the beginning of the function we calculate 'interval', simply the right Index minus the left. If the arguments to 'reFormatMzSampleMatrix' are correct 'interval' should be the same as the user-defined 'chunkSize' - the variable that controls the size of slice from the dictionary. This should only be triggered in the event that there is a remainder left from the division of 'sizeDictAllMzs' by 'chunkSize'.

Checking if 'reFormatMzSampleMatrix' completed fully:
Line 215 and 216 defines 'count' variables for Empty and FullVecs and we increment these depending on the state of 'modified'. 'CountFullVecs' and 'CountEmptyVecs'
should always sum to the same value as 'interval' even for remainder jobs. If they are not the same it means that some of the index values in the 'chunk' did not complete. We can't know yet if it is 'CountFullVecs' and/or 'CountEmptyVecs' that is short but at least the condition will allow us to catch the interval and details if there is a problem here and write it to Err.

### Filehandles and output direction
'reFormatMzSampleMatrix' creates local filehandles 'OUTVEC' (line 205) and 'EMPTID' (line 206). Our primary output, namely a binary vector is written to OUTVEC. The actual filename is a concatenation of an integer 'chunkID'(aka 'jobID' <a jobnumber> ), leftIndex and RightIndex. Ideally we would have liked this interval to correspond to actual mz values instead of just the index values that represent them since this is more informative when looking at the outputs (assuming everything went according plan). We deal with the indexes since we are still in development. You definitely don't want your file opening to happen inside the for loop, that means you cant put it inside 'createVector' or 'formatOutput'. If you need this, just make a function that looks-up the mz values corresponding to the right and left indexes respectively. call this function before creating the file handles and have it return the two ids, then use them when you do the concatenation. This may not be necessary since the mz-ids are the 1st col of the output already. 

### Clarifications 
Note that this function takes place in the context of a slice of dict 'allMzs'. Not to be confused with batch and 'batchSize' - see job launcher \ref{jLaunch}. There is a similar adjustment made for the same type of scenario in \ref{jLaunch}, key differences are that instead of a dict it is 'jobList' and the indices are defined by 'batchSize'.

Under correct operation the value of 'ChunkID' and 'JobID' should be the same since these all equal the same numeric value. There is an exit condition and warning if this is not the case. 

## createVector
Function has two args: 'ShortVector' and 'vecSize'
General: creates a 'long' vector filled with zeros, uses the information in short vector to change some of the zeros in long vector to '1's.
Specific: The length of 'longVec' is equal to numeric arg 'vecSize'. This is literally the number of 'sampleScans' and is stored in 'sizeAllSamVect'. This is also equal to the number of columns in our desired final matrix. Each column is a sampleScanid -represented as an index number. Our relation or reference is literally 'allSamplesNumericDict'. the index number to SampleScanID relation is how we know our column order - it doesn't matter that this is a dictionary.

In 'createVector' we go over each element of 'shortVector'(line 111) and check if there is a match in longVec - if there is a match we modify the entry of longVec from '0' to '1' (line 113). This is the main operation of this function. Note that the 'conditional if' in line 112 has nothing to do with the eval and modification done in line 113. The check in line 112 is to prevent errors. In essence, if everything is correct then all the values that one can encounter in 'shortVec' must exist within the range of 'longVec'. If this con dition returns 'False' (even for just one instance of shortVec ) then we have a problem with how we are creating either short or longVec - hence a fatal stop to the program and error to ERR.  
In 'createVector' we create our the long vector that will eventually end up in the final matrix, however we don't want all-zero rows in our output, but at the same time we do need to know their identity and quantity. 'modified' (line 114) is 'False' by default. We change this to 'True' as soon as we have a hit from 'shortVector'.
If no problems, function returns a modified 'longVec' and status of 'modified'

### formatOutput
A simple function to assemble and properly format a string to be written as a single line in the output file. Function requires 'indexKey', 'valueArray' as args. Returns a single string variable (containing row elements separated by white space). Currently, the output is human-readable binary file. Some commented code is left here which can be uncommented if TPED output is desired.

### makeJobList
a function to create a list of jobs - analogous to creation of 'master.txt' or 'pool'. This function creates the argument sets for 'reFormatMzSampleMatrix' and add each job to a list of jobs. Each instance of 'reFormatMzSampleMatrix' that goes into the pool is discrete and has it's argument set. The arguments that every instance of 'reFormatMzSampleMatrix' requires are 'args=( allMzs, i, leftIndex, rightIndex)'. See \ref{reFormatMzSampleMatrix}   

Important information for understanding the correct operation of this function:
-'Chunks' refers to a groups/sections of dictAllMzs via index. a 'Chunk' consists of <ChunkSize> slice.  
-The usage of the numeric value stored in 'sizeDictAllmzs'
-How the python 'range' operator works
-Handling exceptions and what that means in this context. 
-Why we need the remainder, why do we calculate the ceiling to get our numOfChunks

Function takes one arg, namely 'chunkSize'. 'chunkSize' is a user-defined global variable, (must be an integer) and it is the the number of mzs that we will process in a single chunk. To create a pool of jobs, we 1st need to calculate how many chunks there will be.\begin{verbatim}line 243 rawNumOfChunks = (sizeDictAllMzs / chunkSize )\end{verbatim}.

numOfChunks is directly proportional to the 'sizeDictAllMzs' so for a constant 'chunkSize', the larger the dictionary - the more chunks we will have. Conversely 'numOfChunks' is inversely proportional to the 'chunkSize'. If 'sizeDictAllMzs' is constant, 'numOfChunks' will increase if we decrease 'chunkSize'. 

We take the ceiling of 'sizeDictAllMzs' divided by 'chunkSize' since if there is a remainder, we need to round 'numOfChunks' to the closest ceiling integer. The reason for this become clear when we define the bounds of out chunks in/for the dictionary. We calculate this remainder in \texttt{line 247 dictRemain = remainder(sizeDictAllMzs,chunkSize)} and write it to the log since it is useful for debugging and needed to check if all outputs were created successfully.     

### creating left -and right index numbers
For every chunk we must calculate the coordinates (as index numbers) for an interval of size 'chunkSize' - hence the 'for loop' with range operator in (line 256). The python 'range' operator, when given a single value, will start at zero and iterate up to but not including that value. The condition \texttt{if (i != ceilNumOfChunks-1):} separates the last element in the range via 'else'. 

The default operation - applied to all but the last entry in range - is Create 'leftIndex' and 'rightIndex' where 'leftIndex' is the lowerbound of the interval and 'rightIndex' is the upperbound. \begin{verbatim}line 258 leftIndex = (i * chunkSize)\end{verbatim}
 \begin{verbatim}line 259 rightIndex = (i+1) * chunkSize\end{verbatim}

If i=0 (ie 1st iteration of the loop) 'leftIndex' will remain zero and 'rightIndex' will be equal to the chunkSize. When \ref{reFormatMzSampleMatrix} receives these two arguments, it will retrieve dictionary entries starting and including 'leftIndex' and up to but not including 'rightIndex'. On the second loop iteration of this for loop, 'i' becomes equal to 2, we check if it is the last entry in the range - if not we again times 'i' by 'chunkSize' to get 'lefIndex' and increment by 1 for 'rightIndex'. This continues until 'i' reaches the last value of range at which point we must make an exception: \texttt{lines 264,265 leftIndex = dictSizeIfAllChunksWereFull - chunkSize; rightIndex = int(sizeDictAllMzs)}
'dictSizeIfAllChunksWereFull' is simply the predicted size of the dictionary if we used out ceiling numOfChunks and ChunkSize. It is an overestimate of sizeDictAllMzs and we use it only once here as an easy way to find our left index for the exception case. Our rightIndex will be the same as 'sizeDictAllMzs' since this is the end of our range. Importantly, we explicitly use 'sizeDictAllMzs' since its value is +1 of the last possible index value. The range operator in the receiving function will go up to but not include this value. 

### multiprocessing
\texttt{ line 261 p = multiprocessing.Process(target=reFormatMzSampleMatrix, \\ args=( allMzs, i, leftIndex, rightIndex))}.
This line is how one constructs a job for the multiprocessing pool. Generic form is: 
\texttt{multiprocessing.Process(target=<YourFunctionNameHere>, args= <Argsfor>, <YourFunction>, <Here>)}

"multiprocessing.Process" is a module function. It expects a target with arguments. The target is the name of the function you want to put in the pool. The arguments are what will be passed to your function when it is executed. Once we have our 'job' constructed, we append it to a list. We create this list before we enter the for loop, namely:\texttt{line 255 jobList=[]}. left and right index are created under the main and exception respectively. We only need to define (line 265 ) the jobline and then append (line 266) to the existing joblist for every instance of numOfChunks. The joblines are printed to log and std out.The next step is explained in \section{jLaunch}.

### error check}
\texttt{if (ceilNumOfChunks != lenJobList):} causes program to exit. If this condition not met, something is probably wrong with the way jobList is compiled. This can be directly related to chunkSize or inputSize (by inference dictSize) and mz-range.

## jobLauncher

Function takes two arguments, namely a list of 'jobs' and the lenght of that jobsList then subselects on <batchSize> and sends each selection to a processors and this is done for range of <number of batches>. 'numberOfBatches' is how many 'batchjobs' will be sent to the OS. 

### Making and launching batches
The logic is similar to what we did with 'MakeJobList'.
parallels:
User-defined:
batchSize || chunkSize

Data-Object:
JobList || DictAllMzs (indirect)

Range:
numberOfBatches || ceilNumOfChunks (numberOfBatches is already a ceiling)

incremented numeric identifier:
batchID || chunk

Conditions:
An 'if not equal to' conditions to separate out the last range instance. 

Index definitions:
BatchId * batchSize are used to form the left(lower bound)
(batchID+1) * batchSize for the right - same as 'makeJList'.

For the exception we again use the remainder to get the left index for the last batch. 

Note that here we have to use the paired for loops -each pair called for the main and then once for the exception.
