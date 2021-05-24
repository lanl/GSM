
########################################################
#
# This script tests for differences from 2 files
#
#
# Written by CMJ, XCP-3, July 2018
#
########################################################

# Set amount of output
verbose = 2
# Determine what markers to use when comparing
useL = True
# When output stored, print overall results to stderr
printStdErr = True
# Calculate difference in computation times
compTimeDiff = True

# Imports
import datetime
import sys
import os
import time
import math

# -----------------------------------------------------------------------------------
def clearScreen(numlines=100):
# Thanks to Steven D'Aprano, http://www.velocityreviews.com/forums
# Argument 'numlines' is optional and defaults to 100

  if os.name == "posix":
    # Unix/Linux/MacOS/BSD/etc
    os.system('clear')
  elif os.name in ("nt", "dos", "ce"):
    # DOS/Windows
    os.system('CLS')
  else:
    # Fallback for other operating systems.
    print('\n' * numlines)





# -----------------------------------------------------------------------------------
def getInputName():
    # Check for input argument
    if len(sys.argv) >= 2:
        # Argument passed in, Load input file
        input = sys.argv[1]
        inputName = input.strip() # Remove unnecessary spaces
        
        # Check for valid input type
        if not len(inputName) > 0:
            # User gave spaces for error message, print correction.
            inputError()
        else:
            return inputName


    else:
        # No input file given; correct user and exit script
        inputError()



# ============================================================
def exit():
    # Exit script
  if ( printStdErr ):
    stdErrPrint("\n")
  else:
    print("\n")
  sys.exit()


# -----------------------------------------------------------------------------------
def inputError():
    # Prints error message to user regarding how to use program.
    print("----------------------------------------------------------------------------------------------------")
    print("Run this script to determine what variables are stored in global storage for Fortran source code.\n")
    print("This script requires 1 argument, being an input file name. The input file should list all of the")
    print("files to be compared, where the format is as follows:")
    print("----------------")
    print("[Instructions]")
    print("[Instructions]")
    print("[Instructions]")
    print("[Instructions]")
    print("")
    print("FileA-1")
    print("FileA-2")
    print("OutputNameA.txt")
    print("")
    print("FileB-1")
    print("FileB-2")
    print("OutputNameB.txt")
    print("   .")
    print("   .")
    print("   .")
    print("----------------")
    print("Script is ran by typing 'python Comparer.py [InputFileName]'")
    print("\n")
    print("For questions regarding this script, contact Chase Juneau (chasemjun@lanl.gov) for more information.")
    print("----------------------------------------------------------------------------------------------------")
    exit()


# -----------------------------------------------------------------------------------
def findFile(fileName):
    try:
        file = open(fileName, 'r', encoding='ascii', errors='ignore')
        # Determine number of lines in file
        k = 0
        for k, l in enumerate(file):
            pass
        file.close()
        return (k+1) # Return length of file
    except IOError:
        print("Warning: file ", fileName.strip(), " does not exist in this directory.")
        exit()
    file.close()
        

# -----------------------------------------------------------------------------------
def printTabs(num):
    # Print tabs for easy visualization
    for i in range(0, num, 1):
        # Print tabs
        print('   ', end='')




# -----------------------------------------------------------------------------------
# Copied from "https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python"
def stdErrPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



# -----------------------------------------------------------------------------------
def parseTime(src):
  # Parses out how long the computation took (in seconds)
  text = src.lower().strip()
  text = text[ text.find("=")+len("=") : ].strip().lower() # Remove flag and equals size

  # Get hours
  hourT = 0
  if ( "hr" in text ):
    hourIndx = text.find("hr")
    hourT = float( text[ : hourIndx ].strip() ) # Obtain number of hours
    text = text[ hourIndx + len("hr,") + 1 : ].strip()   # Remove "XX.YY hour," from text

  minT = 0
  if ( "min" in text ):
    minIndx = text.find("min")
    minT = float( text[ : minIndx ].strip() ) # Obtain number of minutes
    text = text[ minIndx + len("min") : ].strip() # Remove "XX.YY min"
    if ( "," in text ):
      # Remove a "," from text
      text[ text.find(",")+1 : ].strip()


  text = text [ text.find("and") + len("and") : ].strip()
  secT = 0
  if ( "sec" in text ):
    secIndx = text.find("sec")
    secT = float( text[ : secIndx ].strip() ) # Obtain number of seconds


  compTime = secT + (60 * minT) + (3600 * hourT)
  return compTime




# -----------------------------------------------------------------------------------
def compare(nameA, nameB, outName):
    global compTimeList
    global numCompTimes

    # Compare files 'nameA' and 'nameB', output goes to 'outName'

    # Obtain length of files
    lengthA = findFile(nameA)
    lengthB = findFile(nameB)
    print("========================================================")
    print("\t\t1. %15s contains %5d line(s)." % (nameA, lengthA) )
    print("\t\t2. %15s contains %5d line(s)." % (nameB, lengthB) )


    # Open a new output file
    try:
        outFile = open(outName, 'w', encoding='ascii', errors='ignore')
    except IOError:
        print("Warning: file '%s' was unable to be created." % outName )
        exit()

    # Print information above to output file
    outFile.write("\n")
    outFile.write("This file contains all of the differences between the files '%s' and '%s'.\n" % (nameA, nameB) )
    outFile.write("\n")
    outFile.write("This script was ran on %s.\n" % datetime.datetime.now())
    outFile.write("\n")
    outFile.write("========================================================\n")
    outFile.write("   1. %15s contains %5d lines.\n" % (nameA, lengthA) )
    outFile.write("   2. %15s contains %5d lines.\n" % (nameB, lengthB) )


    # Open files (existence known by above call to 'findFile')
    fileA = open(nameA, 'r', encoding='ascii', errors='ignore')
    fileB = open(nameB, 'r', encoding='ascii', errors='ignore')


    # Obtain number of differences and file information
    lineA = []
    lineB = []
    lineDiffers = [] # Marks which lines are different
    lineMarkers = [] # Marks differences with a 'v', i.e. will look like "    vvvvv     vvvvv    "
    numDifferences = 0
    noOutput = "*-*-*-*-*- no output -*-*-*-*-*"
    # Set what marker to use for differences
    if ( useL ):
        marker = '|'
    else:
        marker = 'v'

    for i in range(0, max(lengthA, lengthB), 1):

        # Obtain line information
      if ( i < lengthA ):
        lineA.append ( fileA.readline().strip().lower() )
      else:
        lineA.append ( noOutput )
      if ( i < lengthB ):
        lineB.append ( fileB.readline().strip().lower() )
      else:
        lineB.append ( noOutput )


      # Obtain line differences
      if ( not lineA[i] == lineB[i] ):
        lineDiffers.append ( i )
        numDifferences += 1

        # See what the differences are and 'mark' them
        minSize = min( len(lineA[i]), len(lineB[i]) )
        maxSize = max( len(lineA[i]), len(lineB[i]) )
        tempStr = " " * maxSize # Create empty string as long as largest line
        tempA = lineA[i]
        tempB = lineB[i]
        numLDiffs = 0
        for j in range (0, maxSize, 1):
          if ( j > minSize ):
            numLDiffs += 1
            tempStr = tempStr[:j] + marker + tempStr[j+1:]
          else:
            if not ( tempA[j:j+1] == tempB[j:j+1] ):
              numLDiffs += 1
              tempStr = tempStr[:j] + marker + tempStr[j+1:]

        lineMarkers.append ( tempStr + "  (" + str(numLDiffs)  + " diffs)" )


    # Close the files once done reading (information stored in lineA/B lists)
    fileA.close()
    fileB.close()


    # Print out differences
    nameA = "'" + nameA + "'"
    nameB = "'" + nameB + "'"
    maxAllowedDiffs = 3
    # Print to Std Err if desired 
    if ( numDifferences > maxAllowedDiffs ):
        theSummary = "W:   There are %10d differences between %15s (A) and %15s (B)." % ( numDifferences, nameA, nameB )
    else:
        theSummary = "     Effectively NO differences exist between %15s (A) and %15s (B)." % (nameA, nameB)


    # Calculate difference in computation time
    if ( compTimeDiff ):
      flag = "elapsed cpu time (computation)"

      # Find line for comp. time in file
      searchSize = 10
      timeA = -1
      for i in range(max(lengthA - searchSize, 0), lengthA, 1):
        if ( flag in lineA[i] ):
          # Found flag, determine comp time
          timeA = parseTime(lineA[i]) # Get time in seconds
          break
      timeB = -1
      for i in range(max(lengthB - searchSize, 0), lengthB, 1):
        if ( flag in lineB[i] ):
          # Found flag, determine comp time
          timeB = parseTime(lineB[i]) # Get time in seconds
          break

      # Print time information if not the same
      if ( timeA < 0 or timeB < 0 ):
        # Unable to compare
        timeSummary = " [Couldn't calculate comp. times]"
      elif ( not timeA == timeB ):
        compTimeList.append ( 100 * ( (timeB/timeA) - 1) )
        numCompTimes += 1
        if ( timeA > timeB ):
          # B took less time
          type = 1
          timeDiff = 1 - timeB/timeA    # B took this fraction of timeA less
          timeSummary = "[-"
        else:
          # B took more time
          type = 2
          timeDiff = timeB/timeA - 1   # B took this fraction of timeA more
          timeSummary = "[+"

        timeDiff = 100*timeDiff    # Convert to percentage
        timeSummary += "%7.2f%% in comp. time for (B) relative to (A)]" % (timeDiff)
      else:
        # Times are equal
        timeSummary = " [No Time Difference]"

      theSummary  += " %s" % (timeSummary)


    # Print summary information
    if ( printStdErr ):
        stdErrPrint(theSummary)
    else:
        print(theSummary)

    # Print all information
    print("\n   There are %d differences between the two files." % (numDifferences), end='')
    if ( numDifferences <= maxAllowedDiffs ):
        print("\nThey are...")
        print("==================")
        for i in range(0, numDifferences, 1):
            print("\nDifference %d (line %d):" % (i+1, lineDiffers[i]) )
            print("%18s: %s"    % ("-" * 18, lineMarkers[i]) )
            print("   %15s: %s" % (nameA, lineA[ lineDiffers[i] ]) )
            print("   %15s: %s" % (nameB, lineB[ lineDiffers[i] ]) )
    else:
        print(" No more than %d differences are recommended." % (maxAllowedDiffs) )

    # Place some extra space at end
    print("\n\n")

    # Print to output file
    outFile.write("\n   There are %d differences between the two files." % (numDifferences))
    if ( not numDifferences <= maxAllowedDiffs ):
        outFile.write(" No more than %d differences are recommended.\n" % (maxAllowedDiffs) )
    else:
        outFile.write("\n")

    # Write summary information
    outFile.write( ("\n" + theSummary + "\n\n") )
    # Write the differences to the output file
    outFile.write("===========================\n")
    outFile.write("The differences (%d) are:\n" % numDifferences )
    outFile.write("===========================\n")
    for i in range(0, numDifferences, 1):
        outFile.write("\n\n")
        outFile.write("Difference %d (line %d):\n" % (i+1, lineDiffers[i]) )
        outFile.write("%18s: %s\n"    % ("-" * 18, lineMarkers[i]) )
        outFile.write("   %15s: %s\n" % (nameA, lineA[ lineDiffers[i] ]) )
        outFile.write("   %15s: %s\n" % (nameB, lineB[ lineDiffers[i] ]) )

    if ( numDifferences > maxAllowedDiffs ):
      return 1
    else:
      return 0

        

#####################################################################################
#
# Start of actual script
#
#####################################################################################
# Clear terminal
# clearScreen()

# Print some whitespace
scriptStart = "This script was ran on %s.\n" % datetime.datetime.now()
if ( printStdErr ):
  stdErrPrint("\n")
  stdErrPrint(scriptStart)
else:
  print("\n")
  print(scriptStart)

# Determine input name, print to user
inputName = getInputName()

# Determine if file exists
inputlength = findFile(inputName)
if ( verbose >= 1 ):
    print("Reading source file information from input file '", inputName.strip(), "'.\n")


# Read input file for files to read
inputfile = open(inputName, 'r', encoding='ascii', errors='ignore')

# Read through first 4 lines
for i in range(0, 4, 1):
    null = inputfile.readline()

# Now read in sets of 2
numDiffs = 0
numAttempts = 0
compTimeList = []
numCompTimes = 0
for i in range(0, inputlength, 1):

    # Find first file to compare
    searchingA = True
    while ( searchingA ):
        if ( i == inputlength ):
          i = inputlength + 1
          break
        newNameA = inputfile.readline().strip()
        i += 1
        if not ( newNameA == "" ):
            # Found a line with the file name on it
            if ( verbose >= 1 ):
                print("")
                print("======================================================================================")
                print("The following files will be compared:")
                print("======================================================================================")
                print("\t1. %s" % newNameA )
            searchingA = False

    if ( i >= inputlength ):
      break


    # Find second file to compare
    searchingB = True
    while ( searchingB ):
        if ( i >= inputlength ):
            print("Error: File name for differences was not provided on line %d. Please correct." % (i+1) )
            exit()
        newNameB = inputfile.readline().strip()
        i += 1
        if not ( newNameB == "" ):
            # Found a line with the file name on it
            if ( verbose >= 1 ):
                print("\t1. %s" % newNameB )
            searchingB = False


    # Find second file to compare
    searchingO = True
    while ( searchingO ):
        if ( i >= inputlength ):
            print("Error: File name for differences was not provided on line %d. Please correct." % (i+1) )
            exit()
        outName = inputfile.readline().strip()
        i += 1
        if not ( outName == "" ):
            # Found a line with the output file name on it; append '.txt.'
            outName = outName + ".txt"
            if ( verbose >= 1 ):
                print("--> Results will be printed to  %s" % outName )
            searchingO = False


    # Have 2 files, now compare the lines
    numAttempts += 1
    numDiffs += compare(newNameA, newNameB, outName)



thisText = "\n%d of %d comparison tests failed. See resulting output files for details.\n" % (numDiffs, numAttempts)
if ( compTimeDiff and numCompTimes > 0 ):

  # Obtain average
  compAvg = 0
  for i in range(0, numCompTimes, 1):
    compAvg += compTimeList[i]

  # Obtains average computation time (as a fraction), relative to first file
  compAvg = (compAvg / numCompTimes)

  # Compute Std. Dev.
  compStdDev = 0
  for i in range(0, numCompTimes, 1):
    compStdDev += pow( (compTimeList[i] - compAvg), 2 )

  compStdDev = compStdDev / ( numCompTimes - 1 )
  compStdDev = math.sqrt(compStdDev)

  timeSummary = "Average time difference is [%7.2f%%], with a Std. Dev. of [%7.2f%%].\n" % (compAvg, compStdDev)
  thisText += timeSummary


if ( printStdErr ):
    stdErrPrint(thisText)
else:
    print(thisText)


print("All comparisons have been made. Script exiting...")


scriptEnd = "\n-------------------\nThis script completed its task on %s.\n" % datetime.datetime.now()
if ( printStdErr ):
  stdErrPrint(scriptEnd)
else:
  print(scriptEnd)



exit()
