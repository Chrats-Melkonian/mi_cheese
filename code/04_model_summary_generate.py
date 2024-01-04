import sys
import glob
import cobra
import reframed
from reframed import load_cbmodel
from reframed import Environment
from reframed import FBA

# Loop to read in models in hardcoded folder with .xml extension
for file in glob.iglob(r'*.xml'):
	
	# Use try to handle exceptions with reading SMBL model
	try:

		# Load model
		model = load_cbmodel(file, flavor='fbc2')

		# open file
		sys.stdout = open(file+"_summary.txt", "w")

		# get summary
		summary = model.summary()


	# Catch cobra.io.sbml.CobraSBMLError and continue loop    
	except cobra.io.sbml.CobraSBMLError:

		# Print message with model ID for log
		print("The model",file,"raised an SBML error and could not be read in by COBRA")
		pass