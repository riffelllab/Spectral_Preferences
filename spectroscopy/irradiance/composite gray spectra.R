#Loads required libraries
library(dplyr)


#Function for extracting, averaging and converting spectra
led.spectra.extract.full<-function(spec.id,file.name) {
  #Loads spectra
  read.table(paste0("irradiance/LED synth channels/full int/", file.name)) %>%
    #discard extra wavelength columns
    select(V1,V2,V4,V6,V10) %>%
    #renames first column
    rename(wl=V1) %>%
    #takes the mean of all five spectra
    mutate(spec.id=rowMeans(.[,-1])) %>%
    select(wl,spec.id) %>%
    #Converts spectra from µwatts (default unit in spectrasuite) to photons 
    #(see Johnsen 2011, The Optics of Life, Chpt 2)
    #The "1000000" here is to convert µwatts to watts
    mutate(spec.id=spec.id*wl*(5050000000000000/1000000)) %>%
    select(last_col())
}

#normalization functions
norm_by_var<-function(x, y){ x / mean(y, na.rm=TRUE) }


#----------------------------------------------------------------------------
#Load final gray stimuli pwm values determined from photographic measurements
#calibrating led channels both between side a and b and among each other
PWM<-readxl::read_excel("../photography/gray stimuli/2023-04-23/grey-1.00 light spot photo measurment v4.xlsx", 
                range = "K55:M62") %>%
  #rename columns
  rename(wl=1, a_pwm=2, b_pwm=3)


#---------------  Imports full intensity spectra (pwm=4096) -----------------
#Generates a list of spectra files
LED_synth_channel_list.full<-
  data.frame(file.name=list.files("irradiance/LED synth channels/full int/", pattern="*.txt")) %>%
  mutate(spec.id=gsub(".txt", "", file.name)) %>%
  mutate(spec.id=gsub(" ", ".", spec.id))

#Creates empty data frame with wl column, used in for loop below
LED_synth_channels.full<-
  data.frame(wl= read.table( 
    paste0("irradiance/LED synth channels/full int/", 
           LED_synth_channel_list.full$file.name[1]) )[,1]
  )

#For loop to extract and concatenate all the spectra
for(i in seq(1,length(LED_synth_channel_list.full[,1]),1)) {
  LED_synth_channels.full[,i+1] <- led.spectra.extract.full(LED_synth_channel_list.full$spec.id[i], 
                                                        LED_synth_channel_list.full$file.name[i])
}

#renames the columns to match the list
names(LED_synth_channels.full)<-c("wl",LED_synth_channel_list.full$spec.id)

#narrows the spectra to the extent of the led emmission to avoid incorporating noise
#side a
LED_synth_channels.full$a.385[!(LED_synth_channels.full$wl>=370 & LED_synth_channels.full$wl<=435)]<-NA
LED_synth_channels.full$a.395[!(LED_synth_channels.full$wl>=375 & LED_synth_channels.full$wl<=435)]<-NA
LED_synth_channels.full$a.415[!(LED_synth_channels.full$wl>=385 & LED_synth_channels.full$wl<=465)]<-NA
LED_synth_channels.full$a.430[!(LED_synth_channels.full$wl>=395 & LED_synth_channels.full$wl<=475)]<-NA
LED_synth_channels.full$a.450[!(LED_synth_channels.full$wl>=415 & LED_synth_channels.full$wl<=500)]<-NA
LED_synth_channels.full$a.470[!(LED_synth_channels.full$wl>=440 & LED_synth_channels.full$wl<=535)]<-NA
LED_synth_channels.full$a.505[!(LED_synth_channels.full$wl>=466 & LED_synth_channels.full$wl<=560)]<-NA
LED_synth_channels.full$a.525[!(LED_synth_channels.full$wl>=485 & LED_synth_channels.full$wl<=585)]<-NA
LED_synth_channels.full$a.545[!(LED_synth_channels.full$wl>=505 & LED_synth_channels.full$wl<=615)]<-NA
LED_synth_channels.full$a.570[!(LED_synth_channels.full$wl>=540 & LED_synth_channels.full$wl<=600)]<-NA
LED_synth_channels.full$a.590[!(LED_synth_channels.full$wl>=545 & LED_synth_channels.full$wl<=620)]<-NA
LED_synth_channels.full$a.625[!(LED_synth_channels.full$wl>=580 & LED_synth_channels.full$wl<=650)]<-NA
LED_synth_channels.full$a.645[!(LED_synth_channels.full$wl>=590 & LED_synth_channels.full$wl<=675)]<-NA
LED_synth_channels.full$a.660[!(LED_synth_channels.full$wl>=620 & LED_synth_channels.full$wl<=690)]<-NA
LED_synth_channels.full$a.680[!(LED_synth_channels.full$wl>=635 & LED_synth_channels.full$wl<=720)]<-NA
LED_synth_channels.full$a.700[!(LED_synth_channels.full$wl>=655 & LED_synth_channels.full$wl<=755)]<-NA
LED_synth_channels.full$a.740[!(LED_synth_channels.full$wl>=700 & LED_synth_channels.full$wl<=780)]<-NA
#side b
LED_synth_channels.full$b.385[!(LED_synth_channels.full$wl>=370 & LED_synth_channels.full$wl<=435)]<-NA
LED_synth_channels.full$b.395[!(LED_synth_channels.full$wl>=380 & LED_synth_channels.full$wl<=440)]<-NA
LED_synth_channels.full$b.415[!(LED_synth_channels.full$wl>=390 & LED_synth_channels.full$wl<=460)]<-NA
LED_synth_channels.full$b.430[!(LED_synth_channels.full$wl>=395 & LED_synth_channels.full$wl<=480)]<-NA
LED_synth_channels.full$b.450[!(LED_synth_channels.full$wl>=415 & LED_synth_channels.full$wl<=500)]<-NA
LED_synth_channels.full$b.470[!(LED_synth_channels.full$wl>=440 & LED_synth_channels.full$wl<=535)]<-NA
LED_synth_channels.full$b.505[!(LED_synth_channels.full$wl>=466 & LED_synth_channels.full$wl<=560)]<-NA
LED_synth_channels.full$b.525[!(LED_synth_channels.full$wl>=485 & LED_synth_channels.full$wl<=585)]<-NA
LED_synth_channels.full$b.545[!(LED_synth_channels.full$wl>=505 & LED_synth_channels.full$wl<=615)]<-NA
LED_synth_channels.full$b.570[!(LED_synth_channels.full$wl>=540 & LED_synth_channels.full$wl<=600)]<-NA
LED_synth_channels.full$b.590[!(LED_synth_channels.full$wl>=545 & LED_synth_channels.full$wl<=620)]<-NA
LED_synth_channels.full$b.625[!(LED_synth_channels.full$wl>=580 & LED_synth_channels.full$wl<=650)]<-NA
LED_synth_channels.full$b.645[!(LED_synth_channels.full$wl>=590 & LED_synth_channels.full$wl<=675)]<-NA
LED_synth_channels.full$b.660[!(LED_synth_channels.full$wl>=620 & LED_synth_channels.full$wl<=690)]<-NA
LED_synth_channels.full$b.680[!(LED_synth_channels.full$wl>=635 & LED_synth_channels.full$wl<=720)]<-NA
LED_synth_channels.full$b.700[!(LED_synth_channels.full$wl>=655 & LED_synth_channels.full$wl<=755)]<-NA
LED_synth_channels.full$b.740[!(LED_synth_channels.full$wl>=700 & LED_synth_channels.full$wl<=780)]<-NA


#---- Create composite gray spectrum from individual LED channels

#function to scale column using a pwm value looked up
#based on the column name provided
pwm.scale<-function(column, col_name){
  #splits the column name by "." into side (a|b) and a wavelength
  col_name <- strsplit(col_name, "\\.")
  #pulls out side value
  side <- col_name[[1]][1]
  #pulls out wavelength value and makes it numeric
  wl <- as.numeric(col_name[[1]][2])
  #pulls out the correct pwm value from "PWM"
  #depending on "side" and "wl"
  if(side=="a") pwm.value <- PWM$a_pwm[PWM$wl==wl]
  if(side=="b") pwm.value <- PWM$b_pwm[PWM$wl==wl]
  #multiply the column by the pwm value divided by the max pwm value of 4096
  column * (pwm.value/4096)
}

composite.spectra<-
  LED_synth_channels.full %>%
  #narrow down to only the LED channels on side a used in the composite spectra
  select(wl, a.450, a.505, a.525, a.545, a.570, a.590, a.625) %>%
  #scale the spectra from full PWM intensity (4096) down to the level in the composite spectra
  mutate( across( !matches("wl"), ~ pwm.scale(.x, cur_column() ) ) ) %>%
  #take the sum of all the LED channels to create the composite spectra
  mutate(comp_spectra = rowSums(
    data.frame(a.450, a.505, a.525, a.545, a.570, a.590, a.625), na.rm=TRUE)) %>%
  #narrows the composite spectra range to 420-650 nm
  mutate(comp_spectra = case_when(
    wl < 400 ~ NA,
    wl > 660 ~ NA,
    TRUE ~ comp_spectra )) %>%
  #normalize each variable by the mean of the composite spectrum
  mutate(across(!matches("wl"), ~ norm_by_var(.x, comp_spectra)))

#---- outputs the composite spectra to a csv file
write.csv(composite.spectra, file="irradiance/composite gray spectra test.csv", row.names=FALSE)




