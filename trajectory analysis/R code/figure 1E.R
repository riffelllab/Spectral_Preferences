#---------- Loads required packages & functions ----------
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(stringr)
#sets dplyr version of select to be called with "select()"
select<-dplyr::select


#Loads custom R functions
source("R code/wind tunnel functions.R")

#Sets the output directory
output_dir<-"../figures/"


#---------- Experiment Variables ----------

#experiment subfolders
inputFolder<-"input files (examples)"
outputFolder<-"output files"
subFolder<-"color vs gray 0.50, CO2"

#display list of experiment files
exp_file_list<-
  file.path(inputFolder,subFolder) %>%
  list.files(pattern=".csv") %>%
  gsub(pattern=".csv","",.) %>%
  print()
length(exp_file_list)

# name of h5 file and accessory csv file (no file extension)
fileID<-exp_file_list[1]

#stimulus information
testOdor<-"CO2"
neutral_color<-"gray"
neutral_intensity<-"1.00"
control_color<-"gray"
control_intensity<-"0.50"
off_color<-"all"
off_intensity<-"0.00"

# mosquito information
mType<-"wt"
mSex<-"f"

# Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5

# Test section X,Y,Z limits (start_x= -lim$x and end_x= +lim$x --> Total x distance is lim$x*2)
# These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim<-data.frame(x=0.9144, y=0.3048, z=0.6096)

#Size of hypothetical cylinder where if entered the mosquito is 
#considered to have responded to the stimulus. For the odor source this is
#a sphere with same radius
resp_cyl_radius<-0.07 #radius of the cylinder
resp_cyl_height<-0.04 #height of the cylinder from the arena floor

# cue locations
# The locations for "rand" left as NA but defined immediately below
cuesSetup<-data.frame(
  row.names = c("odor","a","b", "rand"),
  x=c(-0.7644,-0.4344,-0.4344,NA),
  y=c(0.00,-0.09,0.09,NA),
  z=c(0.20,0.01,0.01,NA),
  size=c(NA,0.015,0.015,NA)
)

#Choose a random spot ("rand") in the wind tunnel to compare its number of visits
#with the two stimuli positions. The spot is constrained to be downwind of the
#test stimuli by at least twice the "resp_cyl_radius" and to have it's
#"resp_cyl_radius" be completely contained within the limits of the arena
cuesSetup["rand","x"] <- 
  #Generates between the min and the max
  runif(1, min=cuesSetup["a","x"]+2*resp_cyl_radius, max=lim$x-resp_cyl_radius)
cuesSetup["rand","y"] <-
  runif(1, min=-lim$y+resp_cyl_radius, max=lim$y-resp_cyl_radius)


#---------- Create Event Data ----------
#Process the csv file into a table ("event_data") of stimuli and time codes
#Events are segments of an experiment with distinct stimuli.
#IE 505 nm LED light at 100% intensity with C02

# Concatenate file path
file_path<-file.path(inputFolder,subFolder,paste0(fileID,".csv"))
# Read in accessory csv file with time stamps
timeCodes_raw<-read.csv(file_path)

#Assign final output to "event_data"
#"%>% pipes the output of one line to the next
event_data<-
  #Sets initial input to timeCodes_raw
  timeCodes_raw %>%
  #Filters to lines containing "start
  filter(str_detect(event, "start")) %>%
  #renames "timecode" column to "start"
  rename(start=timecode) %>%
  #creates "end" column using the "timecode" from  
  #the rows immediately below (+) the "start" rows 
  mutate(end=timeCodes_raw$timecode[
    grep(pattern="start", timeCodes_raw$event)+1]) %>%
  #creates columns with the genotype and sex of the mosquitos
  mutate(genotype=rep(mType, length(event)), .before=event) %>%
  mutate(sex=rep(mSex, length(event)), .before=event) %>%
  #adds the file name as the rep
  mutate(rep=rep(fileID, length(event)), .before=event) %>%
  #creates a column with the length of the event in minutes
  mutate(event_length=round((end-start)/60)) %>%
  #creates an "order" columns (1 to n) before the "event" column
  mutate(order=seq(1:length(event)), .before=event) %>%
  #remove "start" from text of "event" column
  mutate(event=gsub(pattern="start_", replacement="", event)) %>%
  #Add leds information (color, intensity) split from "event" column 
  #for both light box "a" and "b" (facing upwind a is on the left)
  mutate( a_color=str_split(event,"_",simplify = TRUE)[,3], .before=date.time) %>%
  mutate( a_intensity=str_split(event,"_",simplify = TRUE)[,4], .before=date.time) %>%
  mutate( b_color=str_split(event,"_",simplify = TRUE)[,6], .before=date.time) %>%
  mutate( b_intensity=str_split(event,"_",simplify = TRUE)[,7], .before=date.time) %>%
  #Edits "event" column to remove the LED info by splitting the original "event" column
  mutate( event=str_split(event,"_",simplify = TRUE)[,1]) %>%
  #Creates an odor stimulus column. Sets wait time to "none" and pre and post
  #CO2 time to "air". All other all other values set to "testOdor"
  mutate( odor=case_when(
    event=="WaitTime" ~ "none",
    event=="PreCO2Time" ~ "air",
    event=="resp check, initial no odor" ~"air",
    event=="PostCO2Time" ~ "air",
    TRUE ~ testOdor ), 
    .after=event ) %>%
  #Creates "trt_side" column by looking for the neutral stimulus
  #When neutral stimulus is on only one side, the treatment side is assigned
  #to the oposite side from the neutral stimulus. If the off stimulus or
  #neutral stimulus appear of both sides "trt_side" is set to NA
  mutate( trt_side=case_when(
    #checks for "control_intensity" on both sides during intensity ramp
    #Treatment side is arbitrary, odd "order" #s in upsweep are a, evens are b
    #downsweep is the reverse.
    a_color == control_color & a_intensity == control_intensity &
      b_color == control_color & b_intensity == control_intensity &
      (str_detect(event, "intensity ramp") | 
         str_detect(event, "spectral upsweep")) &
      #checks for even order number
      (order %% 2) == 0 ~ "b",
    a_color == control_color & a_intensity == control_intensity &
      b_color == control_color & b_intensity == control_intensity &
      (str_detect(event, "intensity ramp") | 
         str_detect(event, "spectral upsweep")) &
      #checks for even order number
      (order %% 2) == 1 ~ "a",
    a_color == control_color & a_intensity == control_intensity &
      b_color == control_color & b_intensity == control_intensity &
      (str_detect(event, "intensity ramp") | 
         str_detect(event, "spectral downsweep")) &
      #checks for even order number
      (order %% 2) == 0 ~ "a",
    a_color == control_color & a_intensity == control_intensity &
      b_color == control_color & b_intensity == control_intensity &
      (str_detect(event, "intensity ramp") | 
         str_detect(event, "spectral downsweep")) &
      #checks for even order number
      (order %% 2) == 1 ~ "b",
    #checks for "control_intensity" on both sides
    a_color == control_color & a_intensity == control_intensity &
      b_color == control_color & b_intensity == control_intensity
    ~ as.character(NA),
    #checks for "off_intensity" on both sides 
    a_color == off_color & a_intensity == off_intensity &
      b_color == off_color & b_intensity == off_intensity
    ~ as.character(NA),
    #checks for "neutral_intensity" on both sides 
    a_color == neutral_color & a_intensity == neutral_intensity &
      b_color == neutral_color & b_intensity == neutral_intensity
    ~ as.character(NA),
    #checks for black vs gray 1.00 response check
    a_color == "all" & a_intensity == "0.00" &
      b_color == neutral_color & b_intensity == neutral_intensity &
      event != "intensity ramp" ~ "a",
    b_color == "all" & b_intensity == "0.00" &
      a_color == neutral_color & a_intensity == neutral_intensity &
      event != "intensity ramp" ~ "b",
    #checks for "neutral_intensity" on side a and b
    a_color == control_color & a_intensity == control_intensity ~ "b",
    b_color == control_color & b_intensity == control_intensity ~ "a",
    #If none of these true set value to NA
    TRUE ~ as.character(NA)), 
    .before=a_color) %>%
  #Renames "event" "stim_series"
  rename(stim_series=event) %>%
  #Adds wavelength or intensity stimulus information to the event columns
  #depending on stimulus series. For spectral sweeps, the wavelength info
  #of the "trt_side" is added. For PLogIs (preference vs log intensity) the
  #intensity of the "trt_side" is added.
  mutate( event=case_when(
    #PLogI cases
    str_detect(stim_series, "PLogI") & trt_side=="a" ~
      paste0(stim_series,", intensity ", a_intensity),
    str_detect(stim_series, "PLogI") & trt_side=="b" ~
      paste0(stim_series,", intensity ", b_intensity),
    #intensity ramp cases
    str_detect(stim_series, "intensity ramp") & trt_side=="a" ~
      paste0(stim_series,", intensity ", a_intensity),
    str_detect(stim_series, "intensity ramp") & trt_side=="b" ~
      paste0(stim_series,", intensity ", b_intensity),
    #black vs black
    str_detect(a_color, "all") & str_detect(a_intensity, "0.00") &
      str_detect(b_color, "all") & str_detect(b_intensity, "0.00") ~
      paste0("black vs black"),
    #black vs 100% gray, resp checks
    str_detect(a_color, "all") & str_detect(a_intensity, "0.00") &
      !is.na(trt_side) ~
      stim_series,
    #Spectral sweep cases (both up and downsweeps)
    str_detect(stim_series, "spectral") & trt_side=="a" ~
      paste0(stim_series,", ", a_color," nm"),
    str_detect(stim_series, "spectral") & trt_side=="b" ~
      paste0(stim_series,", ", b_color," nm"),
    #For all other cases keep current stim_series label
    TRUE ~ stim_series),
    .before = stim_series
  )


#---------- Create Spot Data ----------
#These data are associated with a singular insect and singular time point

# Concatenate file path
file_path<-file.path(inputFolder,subFolder,paste0(fileID,".mainbrain.h5"))

#sets the time limits for the experiment
#experiment starts at the "start" of "PreCO2Time 
#and ends at the "end" of "PostCO2Time"
lim_start<-subset(event_data, event==event[2])$start
lim_end<-subset(event_data, event==tail(event,1))$end

spot_data<-
  #Calls H5 import function
  h5read(file_path) %>%
  #Filters out spots beyond the arena test section limits
  filter( x < lim$x & x > -lim$x) %>%
  filter( y < lim$y & y > -lim$y) %>%
  filter( z < lim$z & z > 0) %>%
  #Filters out spots outside the experiment time limits
  filter( timestamp > lim_start & timestamp < lim_end) %>%
  #Do the operations below on each obj_id separately
  group_by(obj_id) %>% 
  #calculate time differences between adjacent spots on the same track
  #The differences are all calculated between the current spot
  #and the previous spot
  mutate( time_diff = timestamp - lag(timestamp) ) %>%
  #calculate XYZ (3D) distance between current spot and previous spot
  mutate( spot_dist = sqrt( (x-lag(x))^2 + (y-lag(y))^2 + (z-lag(z))^2 ) ) %>%
  #Calculates speed from xyz distance and time difference
  mutate( spot_speed = spot_dist/time_diff ) %>%
  #clear the grouping from above
  ungroup() %>%
  #creates empty columns for the for loop below
  mutate( genotype=NA, sex=NA, rep=NA, 
          order=NA, event=NA, stim_series=NA, 
          trt_side=NA, odor=NA,
          a_color=NA, a_intensity=NA, 
          b_color=NA, b_intensity=NA)

#Adds "event_data" information to "spot_data"
#Will produce several "Unknown or uninitialised column" errors
for(i in 1:length(event_data$event)) {
  #Estimate the time stamps values matching each "event"
  event_index<- spot_data$timestamp >= event_data$start[i] & 
    spot_data$timestamp <= event_data$end[i]
  #For all the rows matching the "event_index" write in
  #the corresponding values from "event_data"
  #Spots in between stimuli will have NA values
  spot_data$genotype[event_index] <- event_data$genotype[i]
  spot_data$sex[event_index] <- event_data$sex[i]
  spot_data$rep[event_index] <- event_data$rep[i]
  spot_data$order[event_index] <- event_data$order[i]
  spot_data$event[event_index] <- event_data$event[i]
  spot_data$stim_series[event_index] <- event_data$stim_series[i]
  spot_data$trt_side[event_index] <- event_data$trt_side[i]
  spot_data$odor[event_index] = event_data$odor[i];
  spot_data$a_color[event_index] = event_data$a_color[i];
  spot_data$a_intensity[event_index] = event_data$a_intensity[i];
  spot_data$b_color[event_index] = event_data$b_color[i];
  spot_data$b_intensity[event_index] = event_data$b_intensity[i];
}

spot_data<- spot_data %>%
  #Moves "event_data" columns to the beginning of the table
  relocate(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
           a_color, a_intensity, b_color, b_intensity) %>%
  #Erase the entries recorded between stimulus event_data
  #These entries will have a NA value
  filter(!is.na(event)) %>%
  #Calculates the distances between stimuli locations and the CO2 source
  mutate(odor_dist = sqrt(  (x-cuesSetup["odor","x"])^2 + 
                              (y-cuesSetup["odor","y"])^2 + 
                              (z-cuesSetup["odor","z"])^2 ) ) %>%
  #Calculates the xy distances between stimuli locations and the stimuli
  mutate(a_dist = sqrt(   (x-cuesSetup["a","x"])^2 + 
                            (y-cuesSetup["a","y"])^2) ) %>%
  mutate(b_dist = sqrt(   (x-cuesSetup["b","x"])^2 + 
                            (y-cuesSetup["b","y"])^2 ) ) %>%
  mutate(rand_dist = sqrt(  (x-cuesSetup["rand","x"])^2 + 
                              (y-cuesSetup["rand","y"])^2 ) )


#---------- Create Track Data ----------
#These data are associated with a singular insect track

track_data<-
  #Begins with spot_data as input
  spot_data %>%
  #Groups the spot data by the variables below, the new variables below are
  #calculated for each unique combination of the grouping variables
  group_by(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
           a_color, a_intensity, b_color, b_intensity, obj_id) %>%
  #creates new summary variables (left), using the listed functions
  #from the variables on the right from the original data
  summarise(
    #Calculates the nFrames of the track
    nframes = max(frame, na.rm=TRUE)-min(frame, na.rm=TRUE),
    #Calculates the timestamp of the end of the track
    min_timestamp = min(timestamp, na.rm=TRUE),
    #Calculates the timestamp of the end of the track
    max_timestamp = max(timestamp, na.rm=TRUE),
    #calculates the "duration" of the track in seconds
    duration = max(timestamp, na.rm=TRUE)-min(timestamp, na.rm=TRUE),
    #calculates the total distance of the track
    track_dist = sum(spot_dist, na.rm=TRUE),
    #Removes the grouping specified above 
    .groups = "drop") %>%
  #calculates an average speed along the track from the distance and duration
  mutate( track_speed = track_dist / duration )%>%
  #removes tracks with duration less than "flightTimeLimit"
  filter(duration >= flightTimeLimit)

#removes "spot_data" with a "obj_id" not matching those in the filter 
#"track_data". Removes "spot_data" from tracks with a duration < "flightTimeLimit"
spot_data<- spot_data %>%
  filter(obj_id %in% track_data$obj_id)


#-----Calculates the time spent in the vicinity of each stimulus

#-- For the odor stimulus
track_data <- #outputs to "track_data"
  #input is spot_data
  spot_data %>%
  #filters to row where the distance to the odor source is 
  #less than "resp_cyl_radius"
  filter(odor_dist<=resp_cyl_radius) %>%
  #Groups the spot data by the variables below, the new variables below are
  #calculated for each unique combination of the grouping variables
  group_by(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
           a_color, a_intensity, b_color, b_intensity, obj_id) %>%
  #create a new variable, summing the "time_diff"s of each frame/observation
  #that was within the response sphere. If a track starts inside the 
  #response sphere and exits before the next frame, the value for "time_diff"
  #for the first frame will be NA, and with the "na.rm=TRUE" from above this
  #will result in a 0.
  summarise( odor_time = sum(time_diff, na.rm=TRUE),
             .groups = "drop") %>%
  #joins this data back to "track_data" with NAs being created for those 
  #"obj_id" where no time was spent in the response cylinder
  left_join(x=track_data, y=.) %>%
  #This line will set all the NAs created above to 0s
  mutate( odor_time = replace(odor_time, is.na(odor_time), 0)) 

#-- For the random point (used for comparison with stimuli, see above)
track_data <- #outputs to "track_data"
  #input is spot_data
  spot_data %>%
  #filters to rows where the xy distance to the b stimulus  is 
  #less than "resp_cyl_radius" and the height in the arena is below 4 cm
  filter(rand_dist<=resp_cyl_radius & z<=resp_cyl_height) %>%
  #Groups the spot data by the variables below, the new variables below are
  #calculated for each unique combination of the grouping variables
  group_by(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
           a_color, a_intensity, b_color, b_intensity, obj_id) %>%
  #create a new variable, summing the "time_diff"s of each frame/observation
  #that was within the response cylinder. If a track starts inside the 
  #response cylinder and exits before the next frame, the value for "time_diff"
  #for the first frame will be NA, and with the "na.rm=TRUE" from above this
  #will result in a 0.
  summarise( rand_time = sum(time_diff, na.rm=TRUE),
             .groups = "drop") %>%
  #joins this data back to "track_data" with NAs being created for those 
  #"obj_id" where no time was spent in the response cylinder
  left_join(x=track_data, y=.) %>%
  #This line will set all the NAs created above to 0s
  mutate( rand_time = replace(rand_time, is.na(rand_time), 0)) 

#-- For the a stimulus
track_data <- #outputs to "track_data"
  #input is spot_data
  spot_data %>%
  #filters to rows where the xy distance to the a stimulus  is 
  #less than "resp_cyl_radius" and the height in the arena is below 4 cm
  filter(a_dist<=resp_cyl_radius & z<=resp_cyl_height) %>%
  #Groups the spot data by the variables below, the new variables below are
  #calculated for each unique combination of the grouping variables
  group_by(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
           a_color, a_intensity, b_color, b_intensity, obj_id) %>%
  #create a new variable, summing the "time_diff"s of each frame/observation
  #that was within the response cylinder. If a track starts inside the 
  #response cylinder and exits before the next frame, the value for "time_diff"
  #for the first frame will be NA, and with the "na.rm=TRUE" from above this
  #will result in a 0.
  summarise( a_time = sum(time_diff, na.rm=TRUE),
             .groups = "drop") %>%
  #joins this data back to "track_data" with NAs being created for those 
  #"obj_id" where no time was spent in the response cylinder
  left_join(x=track_data, y=.) %>%
  #This line will set all the NAs created above to 0s
  mutate( a_time = replace(a_time, is.na(a_time), 0)) 

#-- For the b stimulus
track_data <- #outputs to "track_data"
  #input is spot_data
  spot_data %>%
  #filters to rows where the xy distance to the b stimulus  is 
  #less than "resp_cyl_radius" and the height in the arena is below 4 cm
  filter(b_dist<=resp_cyl_radius & z<=resp_cyl_height) %>%
  #Groups the spot data by the variables below, the new variables below are
  #calculated for each unique combination of the grouping variables
  group_by(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
           a_color, a_intensity, b_color, b_intensity, obj_id) %>%
  #create a new variable, summing the "time_diff"s of each frame/observation
  #that was within the response cylinder. If a track starts inside the 
  #response cylinder and exits before the next frame, the value for "time_diff"
  #for the first frame will be NA, and with the "na.rm=TRUE" from above this
  #will result in a 0.
  summarise( b_time = sum(time_diff, na.rm=TRUE),
             .groups = "drop") %>%
  #joins this data back to "track_data" with NAs being created for those 
  #"obj_id" where no time was spent in the response cylinder
  left_join(x=track_data, y=.) %>%
  #This line will set all the NAs created above to 0s
  mutate( b_time = replace(b_time, is.na(b_time), 0)) 

#Creates new track variables based on these new occupancy variables added above
track_data <- 
  #input is "track_data"
  track_data %>%
  #Using "rowwise" allow the mutate command to take the sum "a_time" and 
  #"b_time" on that row rather than the entire column. The "case_when"
  #is necessary as we want to use "na.rm=TRUE" but we want NAs for the
  #case when both "a_time" and "b_time" are NAs. Otherwise we would get
  #a zero here.
  rowwise() %>%
  #Calculates the total time each track spends in each response cylinder
  mutate(stim_time =  sum(a_time,b_time,na.rm=TRUE)) %>%
  #Due the following mutate, row by row
  #rowwise() %>%
  #Calculates the total time each track spends in the treatment cylinder
  #or in the a cylinder if "trt_side" is NA
  mutate(trt_time = case_when(
    is.na(trt_side) ~ a_time,
    trt_side=="a" ~ a_time,
    trt_side=="b" ~ b_time,
    #without "as.numeric" the datatype will not match
    TRUE ~ as.numeric(NA)),
    .before=stim_time) %>%
  #Calculates the total time each track spends in the treatment cylinder
  #or in the a cylinder if "trt_side" is NA
  mutate(ctr_time = case_when(
    is.na(trt_side) ~ b_time,
    trt_side=="a" ~ b_time,
    trt_side=="b" ~ a_time,
    #without "as.numeric" the datatype will not match
    TRUE ~ as.numeric(NA)),
    .before=stim_time) %>%
  #calculates preference index from treatment and control times
  mutate(pref_index = case_when(
    #without "as.numeric" the datatype will not match
    (trt_time+ctr_time)<=0 ~ as.numeric(NA),
    TRUE ~ (trt_time-ctr_time)/(trt_time+ctr_time)
  )) %>%
  #Filter out tracks that do no respond to either stimuli
  filter(stim_time>0) %>%
  #limits to these columns
  select(genotype, sex, rep, order, event, stim_series, trt_side, odor, 
         a_color, a_intensity, b_color, b_intensity, obj_id, duration, 
         trt_time, ctr_time, stim_time, pref_index) %>%
  #Determines treatment color
  mutate(trt_color = case_when(
    #When no trt_side is specified a side is marked as treatment
    is.na(trt_side) ~ a_color,
    trt_side=="a" ~ a_color,
    trt_side=="b" ~ b_color),
    .before = a_color) %>%
  #Determines treatment intensity
  mutate(trt_intensity = case_when(
    #When no trt_side is specified a side is marked as treatment
    is.na(trt_side) ~ a_intensity,
    trt_side=="a" ~ a_intensity,
    trt_side=="b" ~ b_intensity),
    .before = a_color) %>% 
  mutate(ctr_color = case_when(
    #When no trt_side is specified a side is marked as treatment
    is.na(trt_side) ~ b_color,
    trt_side=="a" ~ b_color,
    trt_side=="b" ~ a_color),
    .before = a_color) %>%
  #Determines treatment intensity
  mutate(ctr_intensity = case_when(
    #When no trt_side is specified a side is marked as treatment
    is.na(trt_side) ~ b_intensity,
    trt_side=="a" ~ b_intensity,
    trt_side=="b" ~ a_intensity),
    .before = a_color) %>%
  #Drops a & b stimulus info
  select(!c(a_color, a_intensity, b_color, b_intensity)) %>%
  #moves "obj_id" to first column
  relocate(obj_id)


#----- Creates a version of spot data where the "x" value for portions of 
#trajectories outside the response cylinder are "NA". This change makes
#it so that only this portion of the line will be plotted if this version
#of the data is used
# spot_data_resp_only<-
#   spot_data %>%
#   mutate(x = case_when(
#     !((a_dist<=resp_cyl_radius | b_dist<=resp_cyl_radius) & z<=resp_cyl_height) 
#     ~ as.numeric(NA),
#     TRUE ~ x))


#----- filter responding tracks to identify example tracks
black_test_resp<-
  track_data %>%
  #filter to specfic stimuli with treatment on side
  filter(  trt_color=="all" & 
           trt_intensity=="0.00" & 
           stim_series=="intensity ramp" & 
           trt_side=="b") %>%
  #sort tracks by amount of time spent near the test stimulus
  arrange(desc(pref_index), desc(trt_time))

green_mixed_resp<-
  track_data %>%
  #filter to specfic stimuli with treatment on side
  filter(  trt_color=="545" & 
             trt_intensity=="1.00" & 
             trt_side=="b") %>%
  #sort tracks by amount of time spent near the test stimulus
  arrange(pref_index, desc(ctr_time))

red_ctr_resp<-
  track_data %>%
  #filter to specfic stimuli with treatment on side
  filter(  trt_color=="660" & 
             trt_intensity=="1.00" & 
             trt_side=="b") %>%
  #sort tracks by amount of time spent near the test stimulus
  arrange(pref_index, desc(ctr_time))




#----- figure 1d track1 plot

#sets the track to use in the graph
track_id <- black_test_resp$obj_id[5]

#subsets "spot_data" to the specified track by "obj_id"
track_sub_line <- subset(spot_data, obj_id==track_id)

#subsets the "spot_data" to the specified track by "obj_id"
track_sub <- subset(spot_data, obj_id==track_id)

#Looks up "a" and "b" colors from "color.table" defined with functions
#The "[1]" index uses the first row for the look up. A 70% alpha is then added.
a.color <- 
  color_lookup(track_sub_line$a_color[1], track_sub_line$a_intensity[1]) %>%
  color.add.alpha(70)
b.color <- 
  color_lookup(track_sub_line$b_color[1], track_sub_line$b_intensity[1]) %>%
  color.add.alpha(70)


# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure 1E track1.pdf"),
    width=4, height=1.5, family="sans", pointsize=8)

#sets the margins to 0
par(mai=c(0,0,0,0), omi=c(0,0,0,0))

#creates a black plot
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=c(-0.9144,0.9144), ylim=c(-0.3048,0.3048), xaxs="i", yaxs="i")

#draw side a and b response cylinders
plotrix::draw.circle(x=cuesSetup["a","x"], y=cuesSetup["a","y"],
                     radius=resp_cyl_radius, 
                     border="black", col=NA, lwd=1)
plotrix::draw.circle(x=cuesSetup["b","x"], y=cuesSetup["b","y"],
                     radius=resp_cyl_radius,
                     border="black", col=NA, lwd=1)

#draw side a and b stimuli
plotrix::draw.circle(x=cuesSetup["a","x"], y=cuesSetup["a","y"],
                     radius=cuesSetup["a","size"],
                     border="black", col=a.color, lwd=1)
plotrix::draw.circle(x=cuesSetup["b","x"], y=cuesSetup["b","y"],
                     radius=cuesSetup["b","size"],
                     border="black", col=b.color, lwd=1)


#plots the track
lines(x=track_sub_line$x, y=track_sub_line$y,
      col=hsv(0,0,0.5,0.3), lwd=1.5)

#subsets the version of "spot_data" with the added "NA"s (see above)
#to the specified track by "obj_id"
# track_sub_line_resp <- subset(spot_data_resp_only, obj_id==track_id)
#plots the track but only the portions within the response cylinders
# lines(x=track_sub_line_resp$x, y=track_sub_line_resp$y,
#       col=hsv(0,0,0.2,0.5), lwd=2.5)

#plots a triangle to indicate track starrt
points(x=track_sub_line$x[1], y=track_sub_line$y[1], pch=17, 
       col=hsv(0,0,0.5,0.5), cex=0.8)

#adds 10 cm by 10 cm scale "L"
lines(x=c(-0.8,-0.8,-0.7), y=c(-0.1,-0.2,-0.2), lwd=1.5, col="black")

#Close Graph "Device"
dev.off()




#----- figure 1d track2 plot

#sets the track to use in the graph
track_id <- green_mixed_resp$obj_id[8]

#subsets "spot_data" to the specified track by "obj_id"
track_sub_line <- subset(spot_data, obj_id==track_id)

#subsets the "spot_data" to the specified track by "obj_id"
track_sub <- subset(spot_data, obj_id==track_id)

#Looks up "a" and "b" colors from "color.table" defined with functions
#The "[1]" index uses the first row for the look up. A 70% alpha is then added.
a.color <- 
  color_lookup(track_sub_line$a_color[1], track_sub_line$a_intensity[1]) %>%
  color.add.alpha(70)
b.color <- 
  color_lookup(track_sub_line$b_color[1], track_sub_line$b_intensity[1]) %>%
  color.add.alpha(70)


# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure 1E track2.pdf"),
    width=4, height=1.5, family="sans", pointsize=8)

#sets the margins to 0
par(mai=c(0,0,0,0), omi=c(0,0,0,0))

#creates a black plot
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=c(-0.9144,0.9144), ylim=c(-0.3048,0.3048), xaxs="i", yaxs="i")

#draw side a and b response cylinders
plotrix::draw.circle(x=cuesSetup["a","x"], y=cuesSetup["a","y"],
                     radius=resp_cyl_radius, 
                     border="black", col=NA, lwd=1)
plotrix::draw.circle(x=cuesSetup["b","x"], y=cuesSetup["b","y"],
                     radius=resp_cyl_radius,
                     border="black", col=NA, lwd=1)

#draw side a and b stimuli
plotrix::draw.circle(x=cuesSetup["a","x"], y=cuesSetup["a","y"],
                     radius=cuesSetup["a","size"],
                     border="black", col=a.color, lwd=1)
plotrix::draw.circle(x=cuesSetup["b","x"], y=cuesSetup["b","y"],
                     radius=cuesSetup["b","size"],
                     border="black", col=b.color, lwd=1)


#plots the track
lines(x=track_sub_line$x, y=track_sub_line$y,
      col=hsv(0,0,0.5,0.3), lwd=1.5)

#subsets the version of "spot_data" with the added "NA"s (see above)
#to the specified track by "obj_id"
# track_sub_line_resp <- subset(spot_data_resp_only, obj_id==track_id)
#plots the track but only the portions within the response cylinders
# lines(x=track_sub_line_resp$x, y=track_sub_line_resp$y,
#       col=hsv(0,0,0.2,0.5), lwd=2.5)

#plots a triangle to indicate track starrt
points(x=track_sub_line$x[1], y=track_sub_line$y[1], pch=17, 
       col=hsv(0,0,0.5,0.5), cex=0.8)

#adds 10 cm by 10 cm scale "L"
lines(x=c(-0.8,-0.8,-0.7), y=c(-0.1,-0.2,-0.2), lwd=1.5, col="black")

#Close Graph "Device"
dev.off()




#----- figure 1E track3 plot

#sets the track to use in the graph
track_id <- red_ctr_resp$obj_id[3]

#subsets "spot_data" to the specified track by "obj_id"
track_sub_line <- subset(spot_data, obj_id==track_id)

#subsets the "spot_data" to the specified track by "obj_id"
track_sub <- subset(spot_data, obj_id==track_id)

#Looks up "a" and "b" colors from "color.table" defined with functions
#The "[1]" index uses the first row for the look up. A 70% alpha is then added.
a.color <- 
  color_lookup(track_sub_line$a_color[1], track_sub_line$a_intensity[1]) %>%
  color.add.alpha(70)
b.color <- 
  color_lookup(track_sub_line$b_color[1], track_sub_line$b_intensity[1]) %>%
  color.add.alpha(70)


# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure 1E track3.pdf"),
    width=4, height=1.5, family="sans", pointsize=8)

#sets the margins to 0
par(mai=c(0,0,0,0), omi=c(0,0,0,0))

#creates a black plot
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=c(-0.9144,0.9144), ylim=c(-0.3048,0.3048), xaxs="i", yaxs="i")

#draw side a and b response cylinders
plotrix::draw.circle(x=cuesSetup["a","x"], y=cuesSetup["a","y"],
                     radius=resp_cyl_radius, 
                     border="black", col=NA, lwd=1)
plotrix::draw.circle(x=cuesSetup["b","x"], y=cuesSetup["b","y"],
                     radius=resp_cyl_radius,
                     border="black", col=NA, lwd=1)

#draw side a and b stimuli
plotrix::draw.circle(x=cuesSetup["a","x"], y=cuesSetup["a","y"],
                     radius=cuesSetup["a","size"],
                     border="black", col=a.color, lwd=1)
plotrix::draw.circle(x=cuesSetup["b","x"], y=cuesSetup["b","y"],
                     radius=cuesSetup["b","size"],
                     border="black", col=b.color, lwd=1)


#plots the track
lines(x=track_sub_line$x, y=track_sub_line$y,
      col=hsv(0,0,0.5,0.3), lwd=1.5)

#subsets the version of "spot_data" with the added "NA"s (see above)
#to the specified track by "obj_id"
# track_sub_line_resp <- subset(spot_data_resp_only, obj_id==track_id)
#plots the track but only the portions within the response cylinders
# lines(x=track_sub_line_resp$x, y=track_sub_line_resp$y,
#       col=hsv(0,0,0.2,0.5), lwd=2.5)

#plots a triangle to indicate track starrt
points(x=track_sub_line$x[1], y=track_sub_line$y[1], pch=17, 
       col=hsv(0,0,0.5,0.5), cex=0.8)

#adds 10 cm by 10 cm scale "L"
lines(x=c(-0.8,-0.8,-0.7), y=c(-0.1,-0.2,-0.2), lwd=1.5, col="black")

#Close Graph "Device"
dev.off()
