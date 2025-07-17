**** "trajectory analysis.Rproj" ****
R project file. USE THIS FILE TO OPEN R STUDIO and run the ".R" scripts

Output files, figures and example input files can be found in the [data publication](https://doi.org/10.17632/fdr7znz5dh.2)
associated with [Blake et al. 2025](https://doi.org/10.1242/jeb.250318).
		
		
#---- "R code" ----
R code used to perform single rep/run analyses and produce the figures related to
the trajectory analysis

	#---- "color functions.R" ----
	R source file with custom R functions for creating a list colors to associate with given  
	stimuli (IE 430_0.20). Uses the same naming scheme as used for stimuli in the name of
	Arduino sketches in the "arduino code" folder. Also includes lookup function to go
	from a color (IE all, gray, or 450) and intensity (IE 0.00, 0.50, 1.00) to a hex
	color to use in the graphs. This source file is called in the "wind tunnel functions.R"
	source file.
	
	#---- "wind tunnel functions.R" ----
	Custom functions for importing and graphing trajectory data
	
	#---- "relative intensities.csv" ----
	These intensities were determined from photographic measurements of wind tunnel stimuli 
	using DCRaw import than preserved sensor linearity. Values were drawn from calculations 
	performed in the two following spread sheets. These intensities are used in the analysis 
	of the color intensity ramps.
	"photography/color stimuli/2023-05-10/spectral sweep w projectors on.xlsx"
	"photography/gray stimuli/2023-05-11/gray-0.10 to gray-1.50 light spot photo measurment.xlsx"
	
	#---- "gray sigmoid coefficients.csv" ----
	Saved coefficients for a sigmoid line fit to the gray stimuli (see Fig. 6A). Used in
	Figs. 6,S6.

	#---- "main_single_exp_LED.R" ----
	Imports "input files" from a single run/rep and produces all the files within a
	"output files/[exp_series]/[YYYY][MM][DD]_[HH][MM][SS]" subfolder. 
	
	#---- "figure 1E.R" ----
	R code to import the data from one of the original data file (.h5) and plot individual
	trajectories. Used to create Fig. 1E.
	
	
	#---- "figure 2.R" ----	
	#---- "figure 3.R" ----
	#---- "figure 4.R" ----
	#---- "figure 5.R" ----
	#---- "figure 6.R" ----
	#---- "figure S2.R" ----
	#---- "figure S3.R" ----
	#---- "figure S4.R" ----
	#---- "figure S5.R" ----
	#---- "figure S6.R" ----
	R code necessary to import intermediate data files ("event data.csv" or 
	"tracks data (resp).csv"), perform statistical analyses, and generate figure visuals.
	The list of figures that each data folders is associated with can be found above. 
	"figure S3" and "figure S4" also contain activation analysis of the spectral sweep and
	intensity ramp experiments respectively. 
