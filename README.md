This dataset is a supplement to the article "Spectral preferences of mosquitos
are altered by odors". The article focuses on wind tunnel behavioral bioassay
investigations of mosquito spectral preferences, along with supporting
photographic and spectrographic measurements. The methods describing the
generation of this data are described in the associated article with additional
methodological details given in the subfolder "read me.txt" files. The dataset
also includes Arduino and python code that was used to generate odor and visual
stimuli.

Article Abstract: "Vision underlies many important behaviors in insects
generally and in mosquitos specifically. Mosquito vision plays a role in
predator avoidance, mate finding, oviposition, locating vertebrate hosts, and
vectoring disease. Recent work has shown that when sensitized to CO2, the visual
responses of Aedes aegypti are wavelength-dependent, but little is known about
how other olfactory stimuli can modulate visual responses. The visual cues
associated with flowers, vertebrate hosts, or oviposition sites differs
substantially and it is possible that odors might prime the mosquito visual
system to respond to these different resources. To investigate the interplay of
olfactory and visual cues, we adapted previously used wind tunnel bioassays to
use quasi-monochromatic targets (390-740 nm) created with a novel LED synth. We
coupled these visual targets with CO2 and the odors representative of vertebrate
hosts, floral nectar or oviposition sites and assessed responses via 3D tracking
of female mosquitos. When CO2 alone is present, we observe a lower preference
for wavelengths in the green portion of the visible spectrum with a gradual
increase as wavelengths moved towards the violet and red ends of the spectrum.
However, when odors associated both with flowers and oviposition sites were
present, we observed significant increases in mosquito preference for green
(475-575 nm) stimuli. In contrast when vertebrate host odor was present, we saw
increased preference for stimuli across the entire visible spectrum. These odor
shifts in the mosquito spectral preferences suggest these preferences are not
fixed and shift depending on behavioral context."



This dataset contains (see subfolder "read me.txt" file for additional details):

#---- "arduino code" ----
Arduino sketches uploaded to LED synths by the python scripts in the "python
code" folder to generate the visual cues.

#---- "figures" ----
Folder containing final and intermediate version of figures. R scripts export
figure versions to this folder.

#---- "humidity measurements" ----
Relative humidity measurements of the air plume in the wind tunnel used in Fig.
S2a.

#---- "odor source diagrams" ----
SVG drawings of odor sources used in the figures

#---- "photography" ----
Photographic measurements used to calibrate the colored and gray visual stimuli.

#---- "python code" ----
Python scripts that were used to control the timing of and generate odor and
visual cues.

#---- "spectroscopy" ----
Spectrographic measurements of surfaces and light sources in the wind tunnel.
This data is used in Figs. 1, S1.

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
