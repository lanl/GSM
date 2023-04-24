###########################################################
#
#
# This script is used to copy all files
# from ./resultsCurrent/ to ./resultsPrevious
#
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
#         "sh -e updateOutputs.sh"
#
###########################################################
reset
echo ""
echo "     #============================================#"
echo "     ||                                          ||"
echo "     ||                                          ||"
echo "     ||          -----  UPDATING  -----          ||"
echo "     ||              GSM Regression              ||"
echo "     ||                   Tests                  ||"
echo "     ||                                          ||"
echo "     ||                                          ||"
echo "     #============================================#"
echo ""


# Ensure all directories exist:
###########################################
echo " "
echo " "
echo "============================== Ensuring directories exist ========="
if [ ! -d "./test" ]; then
    echo "Generating test directory..."
    mkdir ./test/
fi
if [ ! -d "./test/regression" ]; then
    echo "Generating regression directory..."
    mkdir ./test/regression/
fi
if [ ! -d "./test/regression/resultsPrevious-BACKUP" ]; then
    echo "Generating backup directory..."
    mkdir ./test/regression/resultsPrevious-BACKUP/
fi
if [ ! -d "./test/regression/resultsPrevious" ]; then
    echo "Generating previous results directory..."
    mkdir ./test/regression/resultsPrevious/
fi
if [ ! -d "./test/regression/resultsCurrent" ]; then
    echo "Generating current results directory..."
    mkdir ./test/regression/resultsCurrent/
fi



# Copy all files from ./resultsCurrent/ to ./resultsPrevious:
###########################################
cd ./test/regression/
echo " "
echo " "
echo "============================== Copying generated output files ========="
### Copy current to previous results:
cp ./resultsCurrent/*.out ./resultsPrevious/
cp ./resultsCurrent/*.txt ./resultsPrevious/
### Copy current previousResults to a backup directory:
cp ./resultsPrevious/*.out ./resultsPrevious-BACKUP/
cp ./resultsPrevious/*.txt ./resultsPrevious-BACKUP/



# Move new result files to those of the old result files
echo " "
echo " "
echo "============================== Renaming files ========="

# Move all current results to previous results
cd ./resultsPrevious
### ( for AlC )
mv AlC.g.out AlC.go.out
mv AlC.g.txt AlC.go.txt
### ( for AuO )
mv AuO.g.out AuO.go.out
mv AuO.g.txt AuO.go.txt
### ( for bLi )
mv bLi.g.out bLi.go.out
mv bLi.g.txt bLi.go.txt
### ( for dU )
mv dU.g.out dU.go.out
mv dU.g.txt dU.go.txt
### ( for liAu )
mv liAu.g.out liAu.go.out
mv liAu.g.txt liAu.go.txt
### ( for nAl )
mv nAl.g.out nAl.go.out
mv nAl.g.txt nAl.go.txt
### ( for nC )
mv nC.g.out nC.go.out
mv nC.g.txt nC.go.txt
### ( for nSi )
mv nSi.g.out nSi.go.out
mv nSi.g.txt nSi.go.txt
### ( for pAg )
mv pAg.g.out pAg.go.out
mv pAg.g.txt pAg.go.txt
### ( for pAu )
mv pAu.g.out pAu.go.out
mv pAu.g.txt pAu.go.txt
### ( for pBe )
mv pBe.g.out pBe.go.out
mv pBe.g.txt pBe.go.txt
### ( for SiCu )
mv SiCu.g.out SiCu.go.out
mv SiCu.g.txt SiCu.go.txt
### ( for pF )
mv pF.g.out pF.go.out
mv pF.g.txt pF.go.txt
### ( for pPb )
mv pPb.g.out pPb.go.out
mv pPb.g.txt pPb.go.txt
### ( for pU )
mv pU.g.out pU.go.out
mv pU.g.txt pU.go.txt
### ( for nLi )
mv nLi.g.out nLi.go.out
mv nLi.g.txt nLi.go.txt
### ( for nCu )
mv nCu.g.out nCu.go.out
mv nCu.g.txt nCu.go.txt
### ( for nBi )
mv nBi.g.out nBi.go.out
mv nBi.g.txt nBi.go.txt
### ( for cTi )
mv cTi.g.out cTi.go.out
mv cTi.g.txt cTi.go.txt
### ( for oC )
mv oC.g.out oC.go.out
mv oC.g.txt oC.go.txt
### ( for dNi )
mv dNi.g.out dNi.go.out
mv dNi.g.txt dNi.go.txt
### ( for ArCa )
mv ArCa.g.out ArCa.go.out
mv ArCa.g.txt ArCa.go.txt
### ( for NeTa )
mv NeTa.g.out NeTa.go.out
mv NeTa.g.txt NeTa.go.txt
### ( for UC )
mv UC.g.out UC.go.out
mv UC.g.txt UC.go.txt
### ( for gAm )
mv gAm.g.out gAm.go.out
mv gAm.g.txt gAm.go.txt
### ( for gXe )
mv gXe.g.out gXe.go.out
mv gXe.g.txt gXe.go.txt
### ( for gU )
mv gU.g.out gU.go.out
mv gU.g.txt gU.go.txt
### ( for gBe )
mv gBe.g.out gBe.go.out
mv gBe.g.txt gBe.go.txt
### ( for gCu )
mv gCu.g.out gCu.go.out
mv gCu.g.txt gCu.go.txt
### ( for piBi )
mv piBi.g.out piBi.go.out
mv piBi.g.txt piBi.go.txt
### ( for piC )
mv piC.g.out piC.go.out
mv piC.g.txt piC.go.txt
### ( for piCr )
mv piCr.g.out piCr.go.out
mv piCr.g.txt piCr.go.txt
# --------------
# SAMPLE
# --------------
# # ( for XXX )
# mv XXX.g.out XXX.go.out
# mv XXX.g.txt XXX.go.txt
# --------------


# DONE:
echo ""
echo ""
echo "     #============================================#"
echo "     ||                                          ||"
echo "     ||                  Done                    ||"
echo "     ||                                          ||"
echo "     #============================================#"
echo ""
