These are the python scripts that were used to control the timing of and
generate odor and visual cues. The script also outputs a list of time stamps as
a csv file that are used as input for the trajectory analysis in combination
with the H5 file output by the Braid 3d tracking system. The Arduino stimuli
(ex 'a_470_3.00_b_gray_0.50') used in the for loops refer to Arduino sketches
in the "arduino code" folder. The stimuli were presented in two different
orders (IE BGR and RGB) to help account for and sequence effect between stimuli
and any diminished response in the later parts of a bioassay run.

#---- "windtunnel_color_ramps_vs_gray_0.50_CO2odor_(BGR|RGB|...).py" ----
These are the python scripts used in the color intensity ramps bioassays. The
three letters at the end of the file name indicate which colors were included. R
= red (625 nm), G = green (525 nm), B = blue (470 nm), Y = yellow (545 nm), V =
violet (430 nm), Gr = composite "gray"

#---- "windtunnel_spectral_sweep_vs_all_0.50_CO2odor_(up|down)sweep.py" ----
These are the python scripts used in the spectral sweep bioassays.

#---- "windtunnel_odor_series_vs_gray_0.50_air_CO2_odor_CO2odor.py" ----
#---- "windtunnel_odor_series_vs_gray_0.50_CO2odor_odor_CO2_air.py" ----
These python scripts were used in bioassays testing the effect of humidity and
foot odor on mosquito response and visual preference.

#---- "windtunnel_odor_series_vs_gray_0.50_air_CO2-01_CO2-05_CO2-10.py" ----
#---- "windtunnel_odor_series_vs_gray_0.50_CO2-10_CO2-05_CO2-01_air.py" ----
These python scripts were used in a bioassay testing the effect of CO2
concentration on mosquito response and visual preference.
