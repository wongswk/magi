# Run components of Wenk19 inference
export PYTHONPATH="../"
python3 ../FGPGM/mainFiles/ProteinTransduction/createExperiments.py
python3 ../FGPGM/mainFiles/ProteinTransduction/getHyperparams.py
python3 ../FGPGM/mainFiles/ProteinTransduction/doFGPGM.py

# Now can open R and run the "plot-Wenk-PTrans.R" to generate visualization
