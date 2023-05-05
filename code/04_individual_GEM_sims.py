import sys
import glob
import cobra
import reframed
from reframed import load_cbmodel
 
from reframed import FBA
from reframed import FVA

# Loop to read in models in hardcoded folder or working directory with .xml extension
for file in glob.iglob(r'*.xml'):
	
	# Use try to handle exceptions with reading SMBL model
	try:

		# Load model
		model = load_cbmodel(file, flavor='fbc2')
		# Close exchange fluxes
		Environment.empty(model,inplace=True)


		# SIM 1: Use full milk media (aerobic)
		env = Environment.from_compounds(["h2o","o2","co2","ca2","cl","cobalt2","cu2","fe2","fe3","h","k","mg2","mn2","mobd","na1","nh4","ni2","pi","so4","zn2","ala__L","asn__L","asp__L","glu__L","gln__L","gly","his__L","ile__L","leu__L","lys__L","orn","phe__L","peamn","pro__L","ser__L","thr__L","trp__L","tyr__L","val__L","lcts","glc__D","gal","gal_bD","cit","lac__D","lac__L","for","ac","oxa","pydx","cbl1","thm","pnto__R","fol","ribflv","nac","btn","but","caproic","octa","dca","ddca","ttdca","ptdca","hdca","ocdca","arach","ttdcea","hdcea","ocdcea","lnlc","arachd","ade","gua","ins","thymd","ura","xan"])
		
		# SIM 1: Solve FVA in full milk media (aerobic)
		sol = FVA(model, constraints=env)
		# SIM 1: Create text file with .fva extension
		text_file = open(file+"_milk_aer.fva", "w")
		text_file.write(str(sol))
		text_file.close()

		# SIM 1: Solve FBA in full milk media (aerobic)
		sol = FBA(model, constraints=env)
		# SIM 1: Create text file with .fba extension
		sys.stdout = open(file+"_milk_aer.fba", "w")
		fluxes = sol.show_values(sort=True)


		# SIM 2: Use Milk media env anaerobic
		env = Environment.from_compounds(["h2o","co2","ca2","cl","cobalt2","cu2","fe2","fe3","h","k","mg2","mn2","mobd","na1","nh4","ni2","pi","so4","zn2","ala__L","asn__L","asp__L","glu__L","gln__L","gly","his__L","ile__L","leu__L","lys__L","orn","phe__L","peamn","pro__L","ser__L","thr__L","trp__L","tyr__L","val__L","lcts","glc__D","gal","gal_bD","cit","lac__D","lac__L","for","ac","oxa","pydx","cbl1","thm","pnto__R","fol","ribflv","nac","btn","but","caproic","octa","dca","ddca","ttdca","ptdca","hdca","ocdca","arach","ttdcea","hdcea","ocdcea","lnlc","arachd","ade","gua","ins","thymd","ura","xan"])
		
		# SIM 2: Solve FVA in full milk (anaerobic)
		sol = FVA(model, constraints=env)
		# SIM 2: Create text file with .fva extension
		text_file = open(file+"_milk_ana.fva", "w")
		text_file.write(str(sol))
		text_file.close()

		# SIM 2: Solve FBA in full milk (anaerobic)
		sol = FBA(model, constraints=env)
		# SIM 2: Create text file with .fba extension
		sys.stdout = open(file+"_milk_ana.fba", "w")
		fluxes = sol.show_values(sort=True)


		# SIM 3: Use minimal milk media (aerobic)
		env = Environment.from_compounds(["btn","ca2","cbl1","cl","cobalt2","cu","cu2","fe2","fe3","fol","k","mg2","mn2","nac","pi","pnto__R","pydx","ribflv","so4","thm","zn2","lcts","ade","gua","thymd","ura","h2o","co2","h","na1","cit","o2"])
		
		# SIM 3: Solve FVA
		sol = FVA(model, constraints=env)
		# SIM 3: Create text file with .fva extension
		text_file = open(file+"_mm_aer.fva", "w")
		text_file.write(str(sol))
		text_file.close()

		# SIM 3: Solve FBA
		sol = FBA(model, constraints=env)
		# SIM 3: Create text file with .fba extension
		sys.stdout = open(file+"_mm_aer.fba", "w")
		fluxes = sol.show_values(sort=True)


		# SIM 4: Use minimal milk media (anaerobic)
		env = Environment.from_compounds(["btn","ca2","cbl1","cl","cobalt2","cu","cu2","fe2","fe3","fol","k","mg2","mn2","nac","pi","pnto__R","pydx","ribflv","so4","thm","zn2","lcts","ade","gua","thymd","ura","h2o","co2","h","na1","cit"])
		
		# SIM 4: Solve FVA
		sol = FVA(model, constraints=env)
		# SIM 4: Create text file with .fva extension
		text_file = open(file+"_mm_ana.fva", "w")
		text_file.write(str(sol))
		text_file.close()

		# SIM 4: Solve FBA
		sol = FBA(model, constraints=env)
		# SIM 4: Create text file with .fba extension
		sys.stdout = open(file+"_mm_ana.fba", "w")
		fluxes = sol.show_values(sort=True)

		
	# Catch cobra.io.sbml.CobraSBMLError and continue loop    
	except cobra.io.sbml.CobraSBMLError:

		# Print message with model ID for log
		print("The model",file,"raised an SBML error and could not be read in by COBRA")
		pass
