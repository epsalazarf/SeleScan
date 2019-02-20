#Concatenates all PopGenome results to directory


echo "STARTED> Plot empirical data:"
cd PopGenome/
echo "Plotting Neutrality..."
Rscript NEUTplot.R
echo "Plotting SFS..."
Rscript SFSplot.R
echo "Plotting SFS..."
Rscript SFSplot.empxsim.R
echo "Plotting Fst..."
Rscript FSTplot.R
echo "Plotting paired Fst..."
Rscript pairFSTplot.R
cd ../rEHH/
echo "Plotting iHS..."
Rscript iHSplot.R
echo "Plotting XP-EHH..."
Rscript XPEHHplot.R
cd ..
echo "FINISHED> All plottings for empirical data"
