# Run components of Wenk19 inference
export PYTHONPATH="../"
python3 ../FGPGM/mainFiles/FitzHughNagumo/createExperiments.py
python3 ../FGPGM/mainFiles/FitzHughNagumo/getHyperparams.py
python3 ../FGPGM/mainFiles/FitzHughNagumo/doFGPGM.py

# Now can open R and run the "plot-Wenk-FN.R" to generate visualization
