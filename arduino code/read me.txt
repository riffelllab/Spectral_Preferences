These are the Arduino sketches uploaded to LED synths by the python scripts in
the "python code" folder to generate the visual cues. Each time the visual cues
change a new sketch is uploaded with the LED pwm remaining constant in the
sketch after upload. The pwm values for each of the 17 colors channels (385-740)
on both side a and side b are specified at the top of each sketch. These are
then assigned to individual pins on 1 of 4 PCA9685 boards. In some cases, some
LED channels will have multiple LEDs in order to help account for intensity
differences between different colors.

ex. "a_430_0.20_b_gray_0.50.ino"
The naming convention gives the stimulus color (385-740 | gray | all) and intensity
(0.00-3.00) on sides a and b.


#---- "color, intensity ramps" ----
These are sketches used in the color intensity ramp bioassays.

#---- "color, vs gray_0.50" ----
There are sketches used in the spectral sweep bioassays. All other odor and CO2 
bioassays also used a subset of these sketches.