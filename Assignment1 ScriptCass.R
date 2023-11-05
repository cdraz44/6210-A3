
##You have to load packages that are to be used onto R studio each time you begin working on a project.  The library() function loads the packages into your session so you can utilize their specific features and functions. 

library(tidyverse)
library(vegan)
library(ggplot2)

##The read_tsv file imports BOLD data
dfBOLD_Oct <- read_tsv(file = "http://www.boldsystems.org/index.php/API_Public/combined?taxon=Octopodidae&format=tsv")


##I want to create a data frame organized by country and BIN uri, as well as count the BINs from each country. I am taking the data from my raw BOLD data and creating a new dataframe. The group_by() function will group the data by both variables which will pipe into the count() function which will give us the number of each BIN uri found in each country.

dfOctBINs.country <- dfBOLD_Oct %>%
  group_by(country, bin_uri) %>%
  count(bin_uri)

##I want to filter this data by removing any non-applicable(NA) values within the 'country' variable by using the !is.na function and saving it to a new dataframe.  This code gave me an error message for a very long time until I found that you had to add a comma after the !is.na function.

 dfOctBINs.country.na.rm <- dfOctBINs.country[!is.na(dfOctBINs.country$country), ]
 
##The pivot_wider function reshapes our data to a format more condusive with the plot that I need to convey my species completeness curve. The names_from function changes the dataframe so that each BIN that was counted in the previous code will be used as a new column with counts of each bin for each country. The values_from function being (n) will specify that each data box will contain a numerical value.
 dfOctBINs.spread.by.country <- pivot_wider(data = dfOctBINs.country.na.rm, names_from = bin_uri, values_from = n)

##This code substitutes all instances of NAs within the dataframe to zero in order to be able to plot the data. This will only affect bin_uri values since we previously removed the NA values from our country variable.
dfOctBINs.spread.by.country[is.na(dfOctBINs.spread.by.country)] <- 0

##This code removes row names and setting the country variable names as row attributes instead since the country variable is a character class.
dfOctBINs.spread.by.country <- dfOctBINs.spread.by.country %>%
  remove_rownames %>%
  column_to_rownames(var = "country")

##Creating the species accumulation curve before plotting.  
OctopodidaeCurve <- specaccum(dfOctBINs.spread.by.country) 

##This will plot my species accumulation curve.  I have included x and y axis labels, a main title, and colouring.  This graph will outline the sampling completeness of the family 'Octopodidae' overall.  If the graph maintains an upward slope, there is the potential for further sampling to yield new species and therefore more BINs in the BOLD data.  If it reaches an asymptote, it means that the family is well samples and further sampling may not yield any further species discovery.
plot(OctopodidaeCurve, xlab = "Countries Sampled", ylab = "BIN Richness", title(main = "Sampling Richness for Family Octopodidae"), ci.col = "Light blue", col="Orange")

##This code gives me raw counts of BINs found in each country so that I can isolate the countries with the highest amount of BINs (species) sampled. 
dfCOUNT <- count(dfBOLD_Oct, country) 

##This code removes NAs from my dfCOUNT countries list.  This allows me to further isolate the countries with the highest amount of BINs sampled.
dfCOUNTNA <- na.omit(dfCOUNT)

##This code sorts the table of counted BINs per country in descending order and isolates the highest 5 countries for further analysis.
 dfSORT <- dfCOUNTNA %>%
   arrange(desc(n)) %>%
   slice(1:5)
 
 ##Here I am indexing by character the values of the five countries with the highest amount of BINs sampled.
 dfSHORT <- dfOctBINs.spread.by.country[c("Mexico", "Japan", "Madagascar", "Australia", "Indonesia"), ]
 
 ##I am loading a colour scheme package that is condusive with those who are colour blind to make my plot more accessible.
 library(viridisLite)
 
 ##Here I am creating a geometric bar plot to visualize the data obtained from isolating and sorting the top BIN values.  
 CountriesAccessible<- ggplot(data = dfSORT, mapping = aes(x = country, y = n, fill = country))  + geom_bar(stat = 'identity') + labs(title = "Species Frequency of Countries With Highest Number of BINs", x = "Country", y = "Number of Species") + scale_fill_viridis_d()+
    vjust = 1.5, colour = "white")+
   theme_bw()
 
 ##
 plot(CountriesAccessible) 


##
library(iNEXT)
 ##
 Matrixswap <- t(dfSHORT)
 ##
 iNEXTOct3 <- iNEXT(Matrixswap, datatype = "abundance")  

 ggiNEXT(x= iNEXTOct3, type = 1, color.var ="Both")
 
 





