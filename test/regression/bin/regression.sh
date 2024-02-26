###########################################################
#
# This script performs all regression test commands
#
###########################################################
#
#
# NOTE: THIS SCRIPT IS TO BE RAN IN THE TOP-LEVEL GSM
#       DIRECTORY!
#
#
###########################################################
#
# TO USE, TYPE THE FOLLOWING IN A SHELL (LINUX) PROMPT
# FROM THE TOP-LEVEL GSM DIRECTORY (quotations omitted):
#         "sh -e regression.sh"
#
###########################################################
set -e   # If errors occur, exit cleanly and return the error message:
reset
echo ""
echo "     #============================================#"
echo "     ||                                          ||"
echo "     ||                                          ||"
echo "     ||        -----   PERFORMING   -----        ||"
echo "     ||              GSM Regression              ||"
echo "     ||                   Tests                  ||"
echo "     ||                                          ||"
echo "     ||                                          ||"
echo "     #============================================#"
echo ""



# Build library:
###################################
echo " "
echo " "
echo "================================== Building GSM ========="
# Remove the build directory completely and re-make everything:
if [ -d "./build" ]; then
    rm -rf ./build
fi
mkdir ./build
cd ./build
if [ -z "${1}" ]; then
  cmake -DCMAKE_BUILD_TYPE="DEBUG" ../
else
  echo "Specified regression compiler: ${1}"
  cmake -DCMAKE_BUILD_TYPE="DEBUG" -DCMAKE_Fortran_COMPILER="${1}" ../
fi
make
reset
make install
make my_gsm

# Go to simulation directory and change name of executible (allows more to be generated if desired)
cd ../my_gsm
cp xgsm1 xgsmReg1



# Run all simulations
############i#######################
function test_input()
{
  inp="${1}"
  echo "...${inp}..."
  ./xgsmReg1 -v=10 -i="../test/regression/inputs/${inp}.g.inp" > ${inp}.g.txt
}
echo " "
echo " "
echo "================================== Performing Simulations ========="
# They are:
# AlC.g.inp  bLi.g.inp  liAu.g.inp  nC.g.inp   pAg.g.inp	pBe.g.inp
# AuO.g.inp  dU.g.inp   nAl.g.inp   nSi.g.inp  pAu.g.inp	SiCu.g.inp
# gAm.g.inp  gBe.g.inp  gCu.g.inp
# piBi.g.inp piC.g.inp  piCr.g.inp
#
#
# Start w/ photon-induced
#################
echo " "
echo "---------------------------------"
echo "Testing photon-induced reactions"
echo "---------------------------------"
# Brems. photons
test_input "gBe"
test_input "gXe"
test_input "gU"
# Others
test_input "gAm"
test_input "gCu"
#
#
# Move to pion-induced
#################
echo " "
echo "---------------------------------"
echo "Testing pion-induced reactions"
echo "---------------------------------"
test_input "piBi"
test_input "piC"
test_input "piCr"
#
#
# Move to neutron-induced
#################
echo " "
echo "---------------------------------"
echo "Testing neutron-induced reactions"
echo "---------------------------------"
test_input "nC"
test_input "nSi"
test_input "nAl"
test_input "nLi"
test_input "nCu"
test_input "nBi"
#
#
# Move to p-induced
#################
echo " "
echo "---------------------------------"
echo "Testing proton-induced reactions"
echo "---------------------------------"
test_input "pAg"
test_input "pAu"
test_input "pBe"
test_input "pF"
test_input "pPb"
test_input "pU"
#
#
# Move to light-ion induced
#################
echo " "
echo "---------------------------------"
echo "Testing light-ion induced reactions"
echo "---------------------------------"
test_input "dU"
test_input "bLi"
test_input "liAu"
test_input "cTi"
test_input "oC"
test_input "dNi"
#
#
# Move to heavy-ion induced reactions
#################
echo " "
echo "---------------------------------"
echo "Testing heavy-ion induced reactions"
echo "---------------------------------"
test_input "SiCu"
test_input "AlC"
test_input "AuO"
test_input "ArCa"
test_input "NeTa"
test_input "UC"


# Move all files for testing
###################################
echo " "
echo " "
echo "================================== Moving Files ========="
if [ ! -d "../test/regression/resultsCurrent" ]; then
   mkdir ../test/regression/resultsCurrent
fi
mv ./*.aux ../test/regression/resultsCurrent/
mv ./*.out ../test/regression/resultsCurrent/
mv ./*.txt ../test/regression/resultsCurrent/
rm xgsmReg1
# Moving files to the directory in which to test:
###################################
if [ ! -d "../test/regression/compare" ]; then
   mkdir ../test/regression/compare
fi
cd ../test/regression/compare/
# Copy from ../bin/:
cp ../bin/comp.prm          ./
cp ../bin/Comparer.py       ./
# Copy results to current directory:
cp ../resultsCurrent/*.out  ./
if [ ! -d "../resultsPrevious" ]; then
  echo "ERROR: No previous results folder was found (../resultsPrevious)."
  echo "ERROR: Cannot perform regression checks on the output."
  echo "ERROR: Exiting..."
  exit 2
else
  cp ../resultsPrevious/*.out ./
fi



# Compare old and new output files:
###################################
echo " "
echo " "
echo "================================= Comparing Files ========="
python3 Comparer.py comp.prm > totalDifferences.txt
# Copy results to ../differences/ directory:
if [ ! -d "../differences" ]; then
   mkdir ../differences
fi
cp *.txt ../differences/
rm *.txt *.out
cd ../differences/


# Finish:
echo ""
echo ""
echo "     #============================================#"
echo "     ||                                          ||"
echo "     ||                  Done                    ||"
echo "     ||                                          ||"
echo "     #============================================#"
echo ""
