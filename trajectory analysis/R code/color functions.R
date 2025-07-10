#loads dplyr if not already in use
library(dplyr)

#Creates a table of "hex.color"s corresponding to stimulus 
#"color" and "intensity"
color.table<-
  #creates data.frame from all combinations of the two factors
  expand.grid(
    color=factor(c(385,395,415,430,450,470,505,525,545,
                   570,590,625,645,660,680,700,740,"gray")),
    #the format produces the character vector with the correct number of 0s
    intensity=factor(format(seq(0,3,0.1), nsmall=2))
  ) %>%
  #turns the factors back into characters / numbers
  mutate(
    color=as.character(color),
    intensity=as.character(intensity),
  ) %>%
  #Creates hue values for the stimuli colors
  mutate(hue=case_when(
    color==385 ~ 198/255,
    color==395 ~ 191/255,
    color==415 ~ 184/255,
    color==430 ~ 177/255,
    color==450 ~ 170/255,
    color==470 ~ 150/255,
    color==505 ~ 113/255,
    color==525 ~ 85/255,
    color==545 ~ 60/255,
    color==570 ~ 50/255,
    color==590 ~ 22/255,
    color==625 ~ 0/255,
    color==645 ~ 0/255,
    color==660 ~ 0/255,
    color==680 ~ 0/255,
    color==700 ~ 0/255,
    color==740 ~ 0/255,
    TRUE ~ 0
  )) %>%
  #Saturation set to 1, except for gray stimuli
  mutate(saturation=case_when(
    color=="gray" ~ 0,
    #Scales saturation down to 0.5 when intensities are above 1.50
    as.numeric(intensity) > 1.5 ~ 1 - 0.5*(as.numeric(intensity)-1.5)/1.5,
    TRUE ~ 1)) %>%
  #Set a max for the value variable, the max is lower for the far red colors
  mutate(value=case_when(
    color==645 ~ 0.9,
    color==660 ~ 0.8,
    color==680 ~ 0.7,
    color==700 ~ 0.6,
    color==740 ~ 0.4,
    TRUE ~ 1
  )) %>%
  #multiplies the max value by the intensity normalized by max intensity
  #and gamma (using a typical encoding value of 0.45) corrected to account
  #for RGB color intensity increasing non-linearly and the intensity values
  #of the stimuli increasing linearly. When intensities are above 1.5
  #value remains at 1 and saturation decreases.
  mutate(value= case_when(
    as.numeric(intensity) <= 1.5 & color != "gray" ~ 
      value*(as.numeric(intensity)/1.5)^0.45,
    as.numeric(intensity) <= 1.5 & color == "gray" ~ 
      0.95*value*(as.numeric(intensity)/1.5)^0.45,
    as.numeric(intensity) > 1.5 & color != "gray"  ~ 1,
    as.numeric(intensity) > 1.5  & color == "gray" ~ 0.95
  )) %>%
  #Sort by "color" and then "intensity"
  arrange(color, intensity) %>%
  #Add row for the "all", intensity "0' stimuli
  add_row(color="all", intensity="0.00", hue=0, saturation=0, value=0) %>%
  #Changes any possible pure white values to a 95% gray
  mutate(value= case_when(
    value==1 & saturation==0 ~ 0.95,
    TRUE ~ value
  )) %>%
  #Creates hex color values from the HSV variables
  mutate(hex.color=hsv(hue, saturation, value))


#function to lookup up a color bases on the stimulus the color and intensity
#must have a matching row in the "color.table
color_lookup <- function(input_color, #variable name of color info
                        input_intensity) #variable name of intensity info
{
  #ensures "input_color" is a character for use below
  input_color <- as.character(input_color)
  #ensures "input_intensity" is a character number with 2 decimal places
  if( is.numeric(input_intensity) ){
    input_intensity <- round(input_intensity, 2) %>% format(nsmall = 2)
  }
  
  #creates dataframe for join
  color_intensity<-data.frame(color=input_color, 
                              intensity=input_intensity)
  #starts with color table
  color.table %>%
    #drops extra columns
    select(!c(hue, saturation, value)) %>%
    #joins colors to data set
    left_join(color_intensity, .,
              by = c("color" = "color", 
                     "intensity" = "intensity") ) %>%
    #isolates to the output column
    pull(hex.color)
}


#Function to add an alpha value between 0-100%, with 100% being fully opaque
color.add.alpha <- function(hex_color, alpha) {
  #converts percentage to value between 0-255
  alpha<-round(255 * alpha/100, 0)
  #converts hex color to rgb, and transposes it for input to function below
  hex_color %>% col2rgb() %>% t() %>%
    #creates new hex color with the specified alpha
    rgb(alpha=alpha, max=255)
}