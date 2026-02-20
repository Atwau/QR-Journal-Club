
# gapminder package

# loading necessary packages

library(gapminder)
library(tidyverse)
library(ggstats)
library(dplyr)
library(dlookr)
library(ggplot2)
library(ggstatsplot)

options(scipen=999) # for decimal formating of outputs, I just don't like
# scientific notation of output display i.e. 2e-16, 2e+16

#-------------------------------------------------------------------------------
# loading the gapminder dataset

data() # checking for datasets in gapminder package
gm <- gapminder # retrieving gapminder dataset and assigning it to new name gm
View(gm) # visual inspection of entire gapminder dataset
#-------------------------------------------------------------------------------

# exploring the data

?gapminder # overview of gapminder data from documentation
head(gm)
str(gm)
plot(gm)

# summary from the exploration
# discrete variables - year, pop(population)
# continuous variables - lifeexp(life expectancy), gdpPercap(GDP per capita)
# nominal variables- country, continent

#-------------------------------------------------------------------------------

# RESEARCH OBJECTIVE

#---primary objective

# To evaluate if there is a significant difference in life expectancy across
# continents

#---secondary objectives

# To evaluate if there is a correlation between population size and 
# GDP per capita
# To evaluate if there is a correlation between GDP per capita and 
# life expectancy


#-------------------------------------------------------------------------------

# DESCRIPTIVE STATISTICS

# summary per continent

# Population size summary
gm_summary_POP <- gm %>% 
  dplyr::select(country, continent, pop) %>% 
  group_by(continent) %>% 
  summarise(mean_POP=mean(pop),
            meadian_POP=median(pop),
            sd_POP=sd(pop),
            se_POP=sd(pop)/sqrt(n()),
            min_POP=min(pop),
            max_POP=max(pop)) %>% 
  print()
#---error bar graph plot for population usind sd_POP
ggplot(gm_summary_POP, aes(x=continent, y=mean_POP))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin = mean_POP-sd_POP, ymax = mean_POP+sd_POP))+
  theme_bw()+
  labs(title='Error Bar Graph for Population Across Continent',
       x='Continent',
       y='Mean Population Size')

#---error bar graph plot for population using se_POP (standard error of mean)
ggplot(gm_summary_POP, aes(x=continent, y=mean_POP))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin = mean_POP-se_POP, ymax = mean_POP+se_POP))+
  theme_bw()+
  labs(title='Error Bar Graph for Population Across Continent',
       x='Continent',
       y='Mean Population Size')

# Life expectancy summaries
gm_summary_LE <- gm %>% 
  dplyr::select(country, continent, lifeExp) %>% 
  group_by(continent) %>% 
  summarise(mean_LE=mean(lifeExp), 
            median_LE=median(lifeExp), 
            sd_LE=sd(lifeExp), 
            max_LE=max(lifeExp), 
            min_LE=min(lifeExp)) %>% 
  print()

#---error bar graph plot for Life Expectancy
ggplot(gm_summary_LE, aes(x=continent, y=mean_LE))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin = mean_LE-sd_LE, ymax = mean_LE+sd_LE))+
  theme_bw()+
  labs(title='Bar Graph with Error Bars for Life Expectancy Across Continents',
       x='Continent',
       y='Mean Life Expectancy')

#---box plot for Life expectancy across continent
ggplot(gm, aes(x=continent, y=lifeExp))+
  geom_boxplot()

#GDP Per capita summary
gm_summary_GDP <- gm %>% 
  dplyr::select(country, continent, gdpPercap) %>% 
  group_by(continent) %>% 
  summarise(mean_GDP=mean(gdpPercap), 
            median_GDP=median(gdpPercap), 
            sd_GDP=sd(gdpPercap), 
            max_GDP=max(gdpPercap), 
            min_GDP=min(gdpPercap)) %>% 
  print()

#---error bar graph plot for GDP per capita
ggplot(gm_summary_GDP, aes(x=continent, y=mean_GDP))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin = mean_GDP-sd_GDP, ymax = mean_GDP+sd_GDP))+
  theme_bw()+
  labs(title='Error Bar Graph for GDP Per Capita Across Continents',
       x='Continent',
       y='Mean GDP Per Capita')

#---box plot for GDP per capita
ggplot(gm, aes(x=continent, y=gdpPercap))+
  geom_boxplot()+
  labs(title = 'GDP Across Continents',
       x='Continent',
       y='GDP Per Capita')

?geom_boxplot
#-------------------------------------------------------------------------------

# TESTING FOR NORMALITY

#---population size---
#------using normality() from dlookr
pop_norm <- gm %>% dplyr::select(country, continent, pop) %>% 
  group_by(continent) %>% 
  normality() %>% 
  print()



#---GDP per capita---
#------using normality() from dlookr
GDP_norm <- gm %>% dplyr::select(country, continent, gdpPercap) %>% 
  group_by(continent) %>% 
  normality() %>% 
  print()

#---Life expectancy---
#------using normality() from dlookr
LE_norm <- gm %>% dplyr::select(country, continent, lifeExp) %>% 
  group_by(continent) %>% 
  normality() %>% 
  print()


# Summary
# All are not normality distributed except for life expectancy in Oceania


#-------------------------------------------------------------------------------
# INFERENTIAL STATISTICS

LE_KW <- kruskal.test(lifeExp~continent, data = gm) %>% 
  print()


# Conclusion
# There is significant difference in Life Expectancy across continents



