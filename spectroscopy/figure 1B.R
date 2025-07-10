#Loads required libraries
library(sfsmisc)
library(dplyr)
library(tidyr)
library(purrr)

#loads custom color functions for looking up color values from 
#a custom color pallette using wavelength and intensity
source("../trajectory analysis/R code/color functions.R")

#Sets the output directory
output_dir<-"../figures/"


#----------------------------------------------------------------------------
#Load final pwm values determined from photographic measurements
#calibrating led channels both between side a and b and among each other
PWM<-readxl::read_excel("../photography/color stimuli/2023-05-10/spectral sweep side comparison v5.xlsx", 
                range = "J94:L111") %>%
  #rename columns
  rename(wl=1, a_pwm=2, b_pwm=3)




#---------------  Imports full intensity spectra (pwm=4096) -----------------

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
    #converts watts to photos
    mutate(spec.id=spec.id*wl*(5050000000000000/1000000)) %>%
    select(last_col())
}

#Generates a list of spectra files
led.spectra.list.full<-
  data.frame(file.name=list.files("irradiance/LED synth channels/full int/", pattern="*.txt")) %>%
  mutate(spec.id=gsub(".txt", "", file.name)) %>%
  mutate(spec.id=gsub(" ", ".", spec.id))

#Creates empty data frame with wl column, used in for loop below
led.spectra.full<-
  data.frame(wl= read.table( 
    paste0("irradiance/LED synth channels/full int/", 
           led.spectra.list.full$file.name[1]) )[,1]
  )

#For loop to extract and concatenate all the spectra
for(i in seq(1, length(led.spectra.list.full[,1]), 1) ) {
  led.spectra.full[,i+1] <- led.spectra.extract.full(led.spectra.list.full$spec.id[i], 
                                                     led.spectra.list.full$file.name[i])
}

#renames the columns to match the list
names(led.spectra.full)<-c("wl",led.spectra.list.full$spec.id)

#narrows the spectra to the extent of the led emmission to avoid incorporating noise
#side a
led.spectra.full$a.385[!(led.spectra.full$wl>=370 & led.spectra.full$wl<=435)]<-NA
led.spectra.full$a.395[!(led.spectra.full$wl>=375 & led.spectra.full$wl<=435)]<-NA
led.spectra.full$a.415[!(led.spectra.full$wl>=385 & led.spectra.full$wl<=465)]<-NA
led.spectra.full$a.430[!(led.spectra.full$wl>=395 & led.spectra.full$wl<=475)]<-NA
led.spectra.full$a.450[!(led.spectra.full$wl>=415 & led.spectra.full$wl<=500)]<-NA
led.spectra.full$a.470[!(led.spectra.full$wl>=440 & led.spectra.full$wl<=535)]<-NA
led.spectra.full$a.505[!(led.spectra.full$wl>=466 & led.spectra.full$wl<=560)]<-NA
led.spectra.full$a.525[!(led.spectra.full$wl>=485 & led.spectra.full$wl<=585)]<-NA
led.spectra.full$a.545[!(led.spectra.full$wl>=505 & led.spectra.full$wl<=615)]<-NA
led.spectra.full$a.570[!(led.spectra.full$wl>=540 & led.spectra.full$wl<=600)]<-NA
led.spectra.full$a.590[!(led.spectra.full$wl>=545 & led.spectra.full$wl<=620)]<-NA
led.spectra.full$a.625[!(led.spectra.full$wl>=580 & led.spectra.full$wl<=650)]<-NA
led.spectra.full$a.645[!(led.spectra.full$wl>=590 & led.spectra.full$wl<=675)]<-NA
led.spectra.full$a.660[!(led.spectra.full$wl>=620 & led.spectra.full$wl<=690)]<-NA
led.spectra.full$a.680[!(led.spectra.full$wl>=635 & led.spectra.full$wl<=720)]<-NA
led.spectra.full$a.700[!(led.spectra.full$wl>=655 & led.spectra.full$wl<=755)]<-NA
led.spectra.full$a.740[!(led.spectra.full$wl>=700 & led.spectra.full$wl<=780)]<-NA
#side b
led.spectra.full$b.385[!(led.spectra.full$wl>=370 & led.spectra.full$wl<=435)]<-NA
led.spectra.full$b.395[!(led.spectra.full$wl>=380 & led.spectra.full$wl<=440)]<-NA
led.spectra.full$b.415[!(led.spectra.full$wl>=390 & led.spectra.full$wl<=460)]<-NA
led.spectra.full$b.430[!(led.spectra.full$wl>=395 & led.spectra.full$wl<=480)]<-NA
led.spectra.full$b.450[!(led.spectra.full$wl>=415 & led.spectra.full$wl<=500)]<-NA
led.spectra.full$b.470[!(led.spectra.full$wl>=440 & led.spectra.full$wl<=535)]<-NA
led.spectra.full$b.505[!(led.spectra.full$wl>=466 & led.spectra.full$wl<=560)]<-NA
led.spectra.full$b.525[!(led.spectra.full$wl>=485 & led.spectra.full$wl<=585)]<-NA
led.spectra.full$b.545[!(led.spectra.full$wl>=505 & led.spectra.full$wl<=615)]<-NA
led.spectra.full$b.570[!(led.spectra.full$wl>=540 & led.spectra.full$wl<=600)]<-NA
led.spectra.full$b.590[!(led.spectra.full$wl>=545 & led.spectra.full$wl<=620)]<-NA
led.spectra.full$b.625[!(led.spectra.full$wl>=580 & led.spectra.full$wl<=650)]<-NA
led.spectra.full$b.645[!(led.spectra.full$wl>=590 & led.spectra.full$wl<=675)]<-NA
led.spectra.full$b.660[!(led.spectra.full$wl>=620 & led.spectra.full$wl<=690)]<-NA
led.spectra.full$b.680[!(led.spectra.full$wl>=635 & led.spectra.full$wl<=720)]<-NA
led.spectra.full$b.700[!(led.spectra.full$wl>=655 & led.spectra.full$wl<=755)]<-NA
led.spectra.full$b.740[!(led.spectra.full$wl>=700 & led.spectra.full$wl<=780)]<-NA




#-----------------  Load transmission data of black tulle  ------------------
black_tulle<-
  read.csv("reflectance and transmission/black tulle tansmission.csv") %>%
  #converts from 0-100 to 0-1
  mutate(across(!matches("wl"),  ~./100))

black_tulle<-
  #multiplies the transmission measurement of 4 layers of black tulle
  #by 3 layers to get the transmission of 7 layers of tulle.
  #This transmission is then interpolated to the wavelength breaks
  #of the LED irradiance measurments
  approx(x=black_tulle$wl, y=black_tulle$tulle_4 * black_tulle$tulle_3, 
         xout=led.spectra.full$wl) %>%
  #convert output to dataframe
  as.data.frame() %>%
  #rename output columns
  rename(wl = 1, tulle_7 = 2)




#---------------------  Calculates predicted spectra -----------------------

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

#create the predicted spectra
led.spectra.pred <-
  #use the full intensity spectra as input
  led.spectra.full %>%
  #scale by the pwm values determined during calibration
  mutate( across( !matches("wl"), ~ pwm.scale(.x, cur_column() ) ) ) %>%
  #intensity of the LEDs is then lowered to account for the transmission
  #of the 7 layers of tulle making up the stimulus bullseye
  mutate( across( !matches("wl"), ~ .x * black_tulle$tulle_7 ) )



#------  output a pdf graph of the predicted a side LED spectra Fig. 1f  ------

#x and y axis limits
xlim<-c(350,800)
ylim<-c(0,2.2e10)

#These cex value set the relative font size relative to pointsize below
cex.axis<-0.75
cex.lab<-1 #mtext is not affected by the cex reduction
#Position of the mtext axis labels in lines, one line 0.2"
mtext.line.x<-1.2
mtext.line.y<-1.6
#tick length
tcl<-(-0.3)
#tick label position, x and y labels need different spacing
mgp.x<-c(3, 0.2, 0)
mgp.y<-c(3, 0.4, 0)
#sets line width for spectra
lwd<-2

#Set figure width and calculate plot sizes and figure height
#All inputs in inches
f_width<-4
f_height<-1.8
m_bottom<-0.3
m_left<-0.4
m_top<-0.1
m_right<-0.1

pdf(file=paste0(output_dir,"figure 1B.pdf"), 
    width=f_width, height=f_height, pointsize = 8)

#Margins
par(mai=c(m_bottom,m_left,m_top,m_right)) #inner margins in inches

plot.new()
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#Axes
axis(1, at=axisTicks(xlim, log=F), mgp=mgp.x, tcl=tcl, cex.axis=cex.axis)
#axis(1, at=seq(xlim[1],xlim[2],5), mgp=mgp.x, tcl=tcl, cex.axis=cex.axis)
axis(2, at=axisTicks(ylim, log=F), mgp=mgp.y, tcl=tcl, cex.axis=cex.axis)
mtext("wavelength (nm)", side=1, line=mtext.line.x)
mtext("photon flux (photons/cm2/s)", side=2, line=mtext.line.y, cex.lab=cex.lab)
box(bty="l")

#draws lines for all a channel led spectra (columns 2 to 18)
#color for each line is determined by the wavelength from a custom color palette
for(i in seq(2, 18, 1) ){
  wl <- names(led.spectra.pred[i]) %>% gsub("[a|b]\\.", "", .)
  lines(led.spectra.pred$wl, led.spectra.pred[,i], lwd=lwd, col=color_lookup(wl,"1.00"))
}

#closes graph pdf
dev.off()
