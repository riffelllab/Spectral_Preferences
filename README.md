This repositories contains code related to [Spectral preferences of mosquitos
are altered by odors](https://doi.org/10.1242/jeb.250318). The article focuses on 
wind tunnel behavioral bioassay investigations of mosquito spectral preferences. The methods describing the generation of this data are described in the associated article with additional methodological details given in subfolder "read me.txt" files. The
repository includes Arduino and python code that was used to generate odor and visual
stimuli as well as R code to analyse the 3d mosquito trajectories recorded during the 
bioassays. Data files, final figures, and intermediate fules can be found in the [data
publication](https://doi.org/10.17632/fdr7znz5dh.2) associated with published article.


#---- "arduino code" ----
Arduino sketches uploaded to LED synths by the python scripts in the "python
code" folder to generate the visual cues.

#---- "python code" ----
Python scripts that were used to control the timing of and generate odor and
visual cues.

#---- "trajectory analysis" ----
R code, example input files, intermediate files, and output files related to the
windtunnel behavioral bioassays. Also has R code for the analysis and creation
of all figures with the exceptions of Figs. 1 & S1.



Most analyses were performed in Rstudio. The version information can be found
below. The version of packages used in the analysis can be found in "required R
packages.csv"

#--- R version ---#                 
platform       aarch64-apple-darwin20 arch           aarch64 os            
darwin20 system         aarch64, darwin20 status major          4 minor         
4.1 year           2024 month          06 day            14 svn rev        86737
language       R version.string R version 4.4.1 (2024-06-14) nickname       Race
for Your Life

#--- RStudio version ---# 
RStudio 2024.09.1+394 "Cranberry Hibiscus" Release
(a1fe401fc08c232d470278d1bc362d05d79753d9, 2024-11-02) for macOS Mozilla/5.0
(Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko)
RStudio/2024.09.1+394 Chrome/124.0.6367.243 Electron/30.4.0 Safari/537.36,
Quarto 1.5.57
