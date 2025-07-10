#Loads required libraries
library(sfsmisc)
library(dplyr)
library(tidyr)
library(purrr)



#---------------  Calculates baseline spectra -----------------

#Reads in the spectra used a baseline for the intensity of the other spectra
#the brightness of the 570 nm is limits the maximum acheivable isoquantal intensity
baseline<-
  read.table("irradiance/LED synth channels/a 570 baseline.txt") %>%
  #discard extra wavelength columns
  select(V1,V2,V4,V6,V10) %>%
  #renames first column
  rename(wl=V1) %>%
  #takes the mean of all five spectra
  mutate(a.570=rowMeans(.[,-1])) %>%
  select(wl,a.570) %>%
  #converts watts to photos
  mutate(a.570=a.570*wl*(5050000000000000/1000000))

#narrows the spectra to the extent of the led emmission
baseline$a.570[!(baseline$wl>=540 & baseline$wl<=600)]<-NA

#functions to integrate the spectra 
spec.int<-function(x) {
  integrate.xy(subset(baseline$wl, !is.na(x)), subset(x, !is.na(x)))
}
#calculates a base intensity to use 
base.int<-spec.int(baseline$a.570)



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
for(i in seq(1,length(led.spectra.list.full[,1]),1)) {
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


#------ Determine peak wavelength for each LED ----------

peak.wl<-
  led.spectra.full %>% 
  summarise_at(vars(!contains("wl")), which.max) %>%
  #This command pivots the data frame from wide to long taking the names
  #as the wavelengths after removing the leading "x"
  pivot_longer(cols=everything(), names_to = "wl", values_to = "peak_wl") %>%
  #splits "wl" column into two
  separate(wl, into=c("side", "wl")) %>%
  #creates a separate column for sides "a" and "b"
  pivot_wider(names_from=side, values_from=peak_wl) %>%
  #turns wavelength to numeric
  mutate(wl=as.numeric(wl)) %>%
  #looks up peak wavelength from "led.spectra.full$wl"
  mutate(a=led.spectra.full$wl[a],
         b=led.spectra.full$wl[b]) %>%
  #averages the peak wavelengths from "a" and "b"
  mutate(peak_wl=rowMeans(cbind(a,b))) %>%
  #Rounds the peak wavelengths
  mutate(peak_wl=round(peak_wl,0)) %>%
  #drops "a" and "b" columns
  select(!c(a,b))

#------ Integrate spectra to get photon flux ----------

#functions to integrate the spectra 
spec.int<-function(x) {
  integrate.xy(subset(led.spectra.full$wl, !is.na(x)), subset(x, !is.na(x)))
}

led.int.full<-led.spectra.full %>%
  #Summarizes each LED using the integration function from above
  #generates a single row data frame
  summarise_at(vars(!contains("wl")), spec.int) %>%
  #This command pivots the data frame from wide to long taking the names
  #as the wavelengths after removing the leading "x"
  pivot_longer(cols=everything(), names_to = "wl", values_to = "photon_flux") %>%
  #splits "wl" column into two
  separate(wl, into=c("side", "wl")) %>%
  #creates a separate column for sides "a" and "b"
  pivot_wider(names_from=side, values_from=photon_flux) %>%
  #turns wavelength to numeric
  mutate(wl=as.numeric(wl)) %>%
  #adds peak wavelength from above
  left_join(peak.wl, by="wl") %>%
  relocate(peak_wl, .after=wl)


#Writes the output to file
write.csv(led.int.full, file="irradiance/LED synth channels/led integrals full.csv", row.names=F)




#---- Imports intensity spectra measured at the testing intensity (1.0) ----

#Function for extracting, averaging and coverting spectra
led.spectra.extract<-function(spec.id,file.name) {
  #Loads spectra
  read.table(paste0("irradiance/LED synth channels/test int (1.0)/", file.name)) %>%
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
led.spectra.list<-
  data.frame(file.name=list.files("irradiance/LED synth channels/test int (1.0)/", pattern="*.txt")) %>%
  mutate(spec.id=gsub(".txt", "", file.name)) %>%
  mutate(spec.id=gsub(" ", ".", spec.id))

#Creates empty data frame with wl column, used in for loop below
led.spectra<-
  data.frame(wl= read.table( 
    paste0("irradiance/LED synth channels/test int (1.0)/", 
           led.spectra.list$file.name[1]) )[,1]
  )

#For loop to extract and concatenate all the spectra
for(i in seq(1,length(led.spectra.list[,1]),1)) {
  led.spectra[,i+1] <- led.spectra.extract(led.spectra.list$spec.id[i], 
                                              led.spectra.list$file.name[i])
}
#renames the columns to match the list
names(led.spectra)<-c("wl",led.spectra.list$spec.id)

#narrows the spectra to the extent of the led emmission to avoid incorporating noise
#side a
led.spectra$a.385[!(led.spectra$wl>=370 & led.spectra$wl<=435)]<-NA
led.spectra$a.395[!(led.spectra$wl>=375 & led.spectra$wl<=435)]<-NA
led.spectra$a.415[!(led.spectra$wl>=385 & led.spectra$wl<=465)]<-NA
led.spectra$a.430[!(led.spectra$wl>=395 & led.spectra$wl<=475)]<-NA
led.spectra$a.450[!(led.spectra$wl>=415 & led.spectra$wl<=500)]<-NA
led.spectra$a.470[!(led.spectra$wl>=440 & led.spectra$wl<=535)]<-NA
led.spectra$a.505[!(led.spectra$wl>=466 & led.spectra$wl<=560)]<-NA
led.spectra$a.525[!(led.spectra$wl>=485 & led.spectra$wl<=585)]<-NA
led.spectra$a.545[!(led.spectra$wl>=505 & led.spectra$wl<=615)]<-NA
led.spectra$a.570[!(led.spectra$wl>=540 & led.spectra$wl<=600)]<-NA
led.spectra$a.590[!(led.spectra$wl>=545 & led.spectra$wl<=620)]<-NA
led.spectra$a.625[!(led.spectra$wl>=580 & led.spectra$wl<=650)]<-NA
led.spectra$a.645[!(led.spectra$wl>=590 & led.spectra$wl<=675)]<-NA
led.spectra$a.660[!(led.spectra$wl>=620 & led.spectra$wl<=690)]<-NA
led.spectra$a.680[!(led.spectra$wl>=635 & led.spectra$wl<=720)]<-NA
led.spectra$a.700[!(led.spectra$wl>=655 & led.spectra$wl<=755)]<-NA
led.spectra$a.740[!(led.spectra$wl>=700 & led.spectra$wl<=780)]<-NA
#side b
led.spectra$b.385[!(led.spectra$wl>=370 & led.spectra$wl<=435)]<-NA
led.spectra$b.395[!(led.spectra$wl>=380 & led.spectra$wl<=440)]<-NA
led.spectra$b.415[!(led.spectra$wl>=390 & led.spectra$wl<=460)]<-NA
led.spectra$b.430[!(led.spectra$wl>=395 & led.spectra$wl<=480)]<-NA
led.spectra$b.450[!(led.spectra$wl>=415 & led.spectra$wl<=500)]<-NA
led.spectra$b.470[!(led.spectra$wl>=440 & led.spectra$wl<=535)]<-NA
led.spectra$b.505[!(led.spectra$wl>=466 & led.spectra$wl<=560)]<-NA
led.spectra$b.525[!(led.spectra$wl>=485 & led.spectra$wl<=585)]<-NA
led.spectra$b.545[!(led.spectra$wl>=505 & led.spectra$wl<=615)]<-NA
led.spectra$b.570[!(led.spectra$wl>=540 & led.spectra$wl<=600)]<-NA
led.spectra$b.590[!(led.spectra$wl>=545 & led.spectra$wl<=620)]<-NA
led.spectra$b.625[!(led.spectra$wl>=580 & led.spectra$wl<=650)]<-NA
led.spectra$b.645[!(led.spectra$wl>=590 & led.spectra$wl<=675)]<-NA
led.spectra$b.660[!(led.spectra$wl>=620 & led.spectra$wl<=690)]<-NA
led.spectra$b.680[!(led.spectra$wl>=635 & led.spectra$wl<=720)]<-NA
led.spectra$b.700[!(led.spectra$wl>=655 & led.spectra$wl<=755)]<-NA
led.spectra$b.740[!(led.spectra$wl>=700 & led.spectra$wl<=780)]<-NA




#---------------  Calculates fitted spectra  -----------------
#Fits the high intensity data to the intensity of 1.0 spectra
#the testing intensity measurements are very noisy so scaling the full
#intensity down to the level of the 1.0 spectra gives a more conistent estimate
#of the photon flux

#scales each spectra by the ratio of 95% quatiles between the full and 1.0 spectra
led.spectra.fit<-
  led.spectra.full %>%
  mutate(a.385=a.385*(quantile(led.spectra$a.385, p=0.95, na.rm=TRUE)/quantile(a.385, p=0.95, na.rm=TRUE))) %>%
  mutate(a.395=a.395*(quantile(led.spectra$a.395, p=0.95, na.rm=TRUE)/quantile(a.395, p=0.95, na.rm=TRUE))) %>%
  mutate(a.415=a.415*(quantile(led.spectra$a.415, p=0.95, na.rm=TRUE)/quantile(a.415, p=0.95, na.rm=TRUE))) %>%
  mutate(a.430=a.430*(quantile(led.spectra$a.430, p=0.95, na.rm=TRUE)/quantile(a.430, p=0.95, na.rm=TRUE))) %>%
  mutate(a.450=a.450*(quantile(led.spectra$a.450, p=0.95, na.rm=TRUE)/quantile(a.450, p=0.95, na.rm=TRUE))) %>%
  mutate(a.470=a.470*(quantile(led.spectra$a.470, p=0.95, na.rm=TRUE)/quantile(a.470, p=0.95, na.rm=TRUE))) %>%
  mutate(a.505=a.505*(quantile(led.spectra$a.505, p=0.95, na.rm=TRUE)/quantile(a.505, p=0.95, na.rm=TRUE))) %>%
  mutate(a.525=a.525*(quantile(led.spectra$a.525, p=0.95, na.rm=TRUE)/quantile(a.525, p=0.95, na.rm=TRUE))) %>%
  mutate(a.545=a.545*(quantile(led.spectra$a.545, p=0.95, na.rm=TRUE)/quantile(a.545, p=0.95, na.rm=TRUE))) %>%
  mutate(a.570=a.570*(quantile(led.spectra$a.570, p=0.95, na.rm=TRUE)/quantile(a.570, p=0.95, na.rm=TRUE))) %>%
  mutate(a.590=a.590*(quantile(led.spectra$a.590, p=0.95, na.rm=TRUE)/quantile(a.590, p=0.95, na.rm=TRUE))) %>%
  mutate(a.625=a.625*(quantile(led.spectra$a.625, p=0.95, na.rm=TRUE)/quantile(a.625, p=0.95, na.rm=TRUE))) %>%
  mutate(a.645=a.645*(quantile(led.spectra$a.645, p=0.95, na.rm=TRUE)/quantile(a.645, p=0.95, na.rm=TRUE))) %>%
  mutate(a.660=a.660*(quantile(led.spectra$a.660, p=0.95, na.rm=TRUE)/quantile(a.660, p=0.95, na.rm=TRUE))) %>%
  mutate(a.680=a.680*(quantile(led.spectra$a.680, p=0.95, na.rm=TRUE)/quantile(a.680, p=0.95, na.rm=TRUE))) %>%
  mutate(a.700=a.700*(quantile(led.spectra$a.700, p=0.95, na.rm=TRUE)/quantile(a.700, p=0.95, na.rm=TRUE))) %>%
  mutate(a.740=a.740*(quantile(led.spectra$a.740, p=0.95, na.rm=TRUE)/quantile(a.740, p=0.95, na.rm=TRUE))) %>%
  mutate(b.385=b.385*(quantile(led.spectra$b.385, p=0.95, na.rm=TRUE)/quantile(b.385, p=0.95, na.rm=TRUE))) %>%
  mutate(b.395=b.395*(quantile(led.spectra$b.395, p=0.95, na.rm=TRUE)/quantile(b.395, p=0.95, na.rm=TRUE))) %>%
  mutate(b.415=b.415*(quantile(led.spectra$b.415, p=0.95, na.rm=TRUE)/quantile(b.415, p=0.95, na.rm=TRUE))) %>%
  mutate(b.430=b.430*(quantile(led.spectra$b.430, p=0.95, na.rm=TRUE)/quantile(b.430, p=0.95, na.rm=TRUE))) %>%
  mutate(b.450=b.450*(quantile(led.spectra$b.450, p=0.95, na.rm=TRUE)/quantile(b.450, p=0.95, na.rm=TRUE))) %>%
  mutate(b.470=b.470*(quantile(led.spectra$b.470, p=0.95, na.rm=TRUE)/quantile(b.470, p=0.95, na.rm=TRUE))) %>%
  mutate(b.505=b.505*(quantile(led.spectra$b.505, p=0.95, na.rm=TRUE)/quantile(b.505, p=0.95, na.rm=TRUE))) %>%
  mutate(b.525=b.525*(quantile(led.spectra$b.525, p=0.95, na.rm=TRUE)/quantile(b.525, p=0.95, na.rm=TRUE))) %>%
  mutate(b.545=b.545*(quantile(led.spectra$b.545, p=0.95, na.rm=TRUE)/quantile(b.545, p=0.95, na.rm=TRUE))) %>%
  mutate(b.570=b.570*(quantile(led.spectra$b.570, p=0.95, na.rm=TRUE)/quantile(b.570, p=0.95, na.rm=TRUE))) %>%
  mutate(b.590=b.590*(quantile(led.spectra$b.590, p=0.95, na.rm=TRUE)/quantile(b.590, p=0.95, na.rm=TRUE))) %>%
  mutate(b.625=b.625*(quantile(led.spectra$b.625, p=0.95, na.rm=TRUE)/quantile(b.625, p=0.95, na.rm=TRUE))) %>%
  mutate(b.645=b.645*(quantile(led.spectra$b.645, p=0.95, na.rm=TRUE)/quantile(b.645, p=0.95, na.rm=TRUE))) %>%
  mutate(b.660=b.660*(quantile(led.spectra$b.660, p=0.95, na.rm=TRUE)/quantile(b.660, p=0.95, na.rm=TRUE))) %>%
  mutate(b.680=b.680*(quantile(led.spectra$b.680, p=0.95, na.rm=TRUE)/quantile(b.680, p=0.95, na.rm=TRUE))) %>%
  mutate(b.700=b.700*(quantile(led.spectra$b.700, p=0.95, na.rm=TRUE)/quantile(b.700, p=0.95, na.rm=TRUE))) %>%
  mutate(b.740=b.740*(quantile(led.spectra$b.740, p=0.95, na.rm=TRUE)/quantile(b.740, p=0.95, na.rm=TRUE)))




#------ Integrate spectra to get photon flux ----------

#functions to integrate the spectra 
spec.int<-function(x) {
  integrate.xy(subset(led.spectra$wl, !is.na(x)), subset(x, !is.na(x)))
}

led.int.fit<-led.spectra.fit %>%
  #Summarizes each LED using the integration function from above
  #generates a single row data frame
  summarise_at(vars(!contains("wl")), spec.int) %>%
  #This command pivots the data frame from wide to long taking the names
  #as the wavelengths after removing the leading "x"
  pivot_longer(cols=everything(), names_to = "wl", values_to = "photon_flux") %>%
  #splits "wl" column into two
  separate(wl, into=c("side", "wl")) %>%
  #creates a separate column for sides "a" and "b"
  pivot_wider(names_from=side, values_from=photon_flux) %>%
  #Calculates the intensity relative to "base.int"
  mutate(a_rel_to_base=a/base.int) %>%
  mutate(b_rel_to_base=b/base.int)

#Writes the output to file
write.csv(led.int.fit, file="irradiance/LED synth channels/led integrals fit.csv", row.names=F)




#---- Plot comparing the fitted photon flux of side a vs b ------
xlim<-c(350,800)
ylim<-c(0,6.5e10)

#These cex value set the relative font size relative to pointsize below
cex.axis<-1
cex.lab<-1 #mtext is not affected by the cex reduction
cex.leg<-1
#Position of the mtext axis labels in lines, one line 0.2"
mtext.line.x<-2.3
mtext.line.y<-2.2
#tick length
tcl<-(-0.4)
#tick label position, x and y labels need different spacing
mgp.x<-c(3, 0.75, 0)
mgp.y<-c(3, 0.5, 0)
#sets line width for spectra
lwd<-2.5

#Set figure width and calculate plot sizes and figure height
#All inputs in inches
f_width<-8
f_height<-0.4*f_width
m_bottom<-0.7
m_left<-0.7
m_top<-0.2
m_right<-0.2

pdf(file="irradiance/LED synth channels/fitted stimulus intensity.pdf", 
    width=f_width, height=f_height, pointsize = 10)

#sets margins
par(mai=c(0,0,0,0), omi=c(m_bottom,m_left,m_top,m_right))

plot.new()
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#Axes
axis(1, at=axisTicks(xlim, log=F), mgp=mgp.x, tcl=tcl, cex.axis=cex.axis)
#axis(1, at=seq(xlim[1],xlim[2],5), mgp=mgp.x, tcl=tcl, cex.axis=cex.axis)
axis(2, at=axisTicks(ylim, log=F), mgp=mgp.y, tcl=tcl, cex.axis=cex.axis)
mtext("wavelength (nm)", side=1, line=mtext.line.x)
mtext("photon flux (photons/cm2/s)", side=2, line=mtext.line.y, cex.lab=cex.lab)
box(bty="l")

#a spectra
lines(led.spectra.fit$wl, led.spectra.fit$a.385, lwd=lwd, col=hsv(198/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.395, lwd=lwd, col=hsv(190/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.415, lwd=lwd, col=hsv(183/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.430, lwd=lwd, col=hsv(176/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.450, lwd=lwd, col=hsv(169/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.470, lwd=lwd, col=hsv(152/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.505, lwd=lwd, col=hsv(119/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.525, lwd=lwd, col=hsv(102/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.545, lwd=lwd, col=hsv(85/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.570, lwd=lwd, col=hsv(50/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.590, lwd=lwd, col=hsv(32/255,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.625, lwd=lwd, col=hsv(0,1,1))
lines(led.spectra.fit$wl, led.spectra.fit$a.645, lwd=lwd, col=hsv(0,1,0.9))
lines(led.spectra.fit$wl, led.spectra.fit$a.660, lwd=lwd, col=hsv(0,1,0.8))
lines(led.spectra.fit$wl, led.spectra.fit$a.680, lwd=lwd, col=hsv(0,1,0.7))
lines(led.spectra.fit$wl, led.spectra.fit$a.700, lwd=lwd, col=hsv(0,1,0.6))
lines(led.spectra.fit$wl, led.spectra.fit$a.740, lwd=lwd, col=hsv(0,1,0.4))

#b spectra
lines(led.spectra.fit$wl, led.spectra.fit$b.385, lwd=lwd, col=hsv(198/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.395, lwd=lwd, col=hsv(190/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.415, lwd=lwd, col=hsv(183/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.430, lwd=lwd, col=hsv(176/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.450, lwd=lwd, col=hsv(169/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.470, lwd=lwd, col=hsv(152/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.505, lwd=lwd, col=hsv(119/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.525, lwd=lwd, col=hsv(102/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.545, lwd=lwd, col=hsv(85/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.570, lwd=lwd, col=hsv(50/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.590, lwd=lwd, col=hsv(32/255,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.625, lwd=lwd, col=hsv(0,0.5,1))
lines(led.spectra.fit$wl, led.spectra.fit$b.645, lwd=lwd, col=hsv(0,0.5,0.9))
lines(led.spectra.fit$wl, led.spectra.fit$b.660, lwd=lwd, col=hsv(0,0.5,0.8))
lines(led.spectra.fit$wl, led.spectra.fit$b.680, lwd=lwd, col=hsv(0,0.5,0.7))
lines(led.spectra.fit$wl, led.spectra.fit$b.700, lwd=lwd, col=hsv(0,0.5,0.6))
lines(led.spectra.fit$wl, led.spectra.fit$b.740, lwd=lwd, col=hsv(0,0.5,0.4))

#closes graph pdf
dev.off()


