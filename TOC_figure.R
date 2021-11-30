######### TOC-figure ####
#
# Graphical Abstract
# In addition to the Abstract, which will appear at the beginning of the text, the author
# should provide a graphical abstract (TOC). The TOC graphic should capture the reader's
# attention and, in conjunction with the manuscript title, should give the reader a quick visual
# impression of the essence of the paper. It should be in the form of a structure, graph,
# drawing, TEM/SEM/AFM micrograph, or reaction scheme without any added text. Tables or
# spectra are not acceptable. Color is encouraged. Type size of labels, formulas, or numbers
# within the graphic must be legible. The graphic should be submitted as a separate file in
# either TIFF, JPG, Word, or Powerpoint format.

require(tidyverse)
require(ggthemes)

#Ladataan ennustedata
predictions.df <- read.csv2(file = "data/Prediction_vs_rheometer_experiment.csv", 
                            sep = ",",
                            dec = ",")

#Muokataan hieman materiaalien jne. nimia ja valitaan tarvittavat kolumnit

predictions.df <- predictions.df %>% 
  transmute(
    time_min = as.character(time_min),
    material = case_when(
      material == "MED-4735" ~ "Platinum cure",
      material == "MED4-4516" ~ "Peroxide cure"
    ),
    sample = sample,
    temp_C = temp_C,
    conversion = as.character(conversion),
    method = case_when(
      type == "exp" ~ "MDR experiment",
      type == "pred" ~"Kinetic prediction"
    )
    ) %>% 
  mutate(
    time_min = as.numeric(time_min),
    material = as.factor(material),
    temp_C = as.factor(temp_C),
    sample = as.factor(sample),
    method = as.factor(method),
    conversion = as.numeric(conversion)
  )

## VAihdetaan temp_C muuttujaan astemerkki
predictions.df <- predictions.df %>% 
  mutate(
    temp_C = paste(temp_C, " °C", sep = "")
  ) %>% 
  mutate(
    temp_C = as.factor(temp_C)
  )

### Luodaan TOC-kuva

theme_set(theme_clean())

TOC_plot <- predictions.df %>% ggplot(aes(x = time_min, y = conversion, color = sample, lty = method)) +
  geom_line(size = 1.5) +
  facet_grid(material~temp_C) +
  labs(x = "Time (min)", 
       y = "Curing reaction progress (%)", 
       caption = "Cross-linking of silicone rubber") +
  theme(legend.position = "top")


ggsave(filename = "FinalFigures_files/TOC_plot.tiff", 
       device = tiff(), units = "cm",
       height = 5, width = 13, 
       dpi = 300,
       scale = 1.7)
#####################
## TOC-kuva V2
#####################

#Loading the dataset
pred_exp.df <- read.csv(file = "data/Prediction_vs_rheometer_experiment.csv")

#Add celsius mark and make factors factors
pred_exp.df <- pred_exp.df %>%
  mutate(temp2 = if_else(grepl(x = temp_C, pattern = "100"), "100 °C", "120 °C")) %>% 
  mutate(across(c("material", "sample", "temp2", "type"), as.factor))

## first and second derivative for conversion
pred_exp.df <- pred_exp.df %>% 
  group_by(temp2, material, type, sample) %>% 
  mutate(dif_alpha = (lead(conversion) - conversion)/(lead(time_min) - time_min)) %>% # first derivative
  mutate(dif2_alpha = (lead(dif_alpha) - dif_alpha)/(lead(time_min) - time_min)) %>% # second derivative
  ungroup()

#Finding the maxima for the second derivative
#Let's take the mean values as estimates for observed lag, tau_obs
tau_obs <- pred_exp.df %>% 
  filter(type == "exp") %>% 
  group_by(type, temp2, material, sample) %>% 
  summarise(max_2nd = max(dif2_alpha, na.rm=TRUE), 
            time_min = time_min[which.max(dif2_alpha)]) %>% 
  ungroup() %>% 
  group_by(type, temp2, material) %>% 
  summarise(tau_obs = mean(time_min), sd = sd(time_min), n = n()) %>% 
  ungroup() %>% 
  select(-type)

#Calculating values for intrinsic lag, tau_i
tau_i <- pred_exp.df %>% 
  filter(type == "pred") %>% 
  group_by(type, temp2, material, sample) %>% 
  summarise(max_2nd = max(dif2_alpha, na.rm=TRUE), 
            time_min = time_min[which.max(dif2_alpha)]) %>% 
  ungroup() %>% 
  group_by(type, temp2, material) %>% 
  summarise(tau_i = mean(time_min), sd = sd(time_min), n = n()) %>% 
  mutate(tau_i = if_else(material == "MED4-4516", 0, tau_i,),
         sd = if_else(material == "MED4-4516", 0, sd)) %>% 
  ungroup() %>% 
  select(-type)

#Calculating experimental lag, tau_e
lags <- full_join(select(tau_obs, -sd, -n), select(tau_i, -sd, -n)) %>% 
  mutate(tau_e = tau_obs - tau_i)
#Finding the maximum value for G' at 120°C
max_G <- pred_exp.df %>% filter(material == "MED4-4516") %>%
  group_by(material) %>%
  summarise(max_G = max(G., na.rm = TRUE)) %>%
  pull(max_G)

# Correcting for experimental lag
pred_exp.df <- pred_exp.df %>% 
  mutate(tau_e = case_when(
    material == "MED-4735" & temp2 == "100 °C" ~ (lags %>% 
                                                    filter(grepl("4735", material), 
                                                           grepl(100, temp2)) %>% 
                                                    pull(tau_e)),
    material == "MED-4735" & temp2 == "120 °C" ~ (lags %>% 
                                                    filter(grepl("4735", material), 
                                                           grepl(120, temp2)) %>% 
                                                    pull(tau_e)),
    material == "MED4-4516" & temp2 == "100 °C" ~ (lags %>% 
                                                     filter(grepl("4516", material), 
                                                            grepl(100, temp2)) %>% 
                                                     pull(tau_e)),
    material == "MED4-4516" & temp2 == "120 °C" ~ (lags %>% 
                                                     filter(grepl("4516", material), 
                                                            grepl(120, temp2)) %>% 
                                                     pull(tau_e))
  ))

#Finding G_tot
pred_exp.df <- pred_exp.df %>% group_by(type, material, temp2, sample) %>% 
  mutate(G_tot = max(G.)) %>% 
  ungroup()

## Renaming material and method

pred_exp.df <- pred_exp.df %>% 
  mutate(
    material = case_when(
      material == "MED-4735" ~ "Platinum cure",
      material == "MED4-4516" ~ "Peroxide cure"
  ), method = case_when(
      type == "exp" ~ "MDR experiment",
      type == "pred" ~"Kinetic prediction"
  )) %>% 
  mutate(across(.cols = c("material", "method"), as.factor))


#Shifted time and corrected G' for MED4-4516 @ 100°C
#Plotting out the corrected conversion plot

TOC_plot2 <- pred_exp.df %>% mutate(time_min = if_else(type == "exp", time_min - tau_e, time_min)) %>% 
  mutate(conversion = if_else(
    (type == "exp" & temp2 == "100 °C" & material == "Peroxide cure"), 
    true = (G_tot/max_G)*conversion, 
    false = conversion)) %>%
  ggplot(mapping = aes(x = time_min, y = conversion)) +
  geom_line(aes(linetype = method, color = sample), size = 1.5) +
  facet_grid(material ~ temp2) +
  coord_cartesian(xlim = c(0,5)) +
  labs(x = "Time (min)", 
       y = "Curing reaction progress (%)", 
       color = "Sample #", 
       linetype = "Method:",
       caption = "Cross-linking of silicone rubber") +
  theme(legend.position = "top")

ggsave(filename = "FinalFigures_files/TOC_plot.tiff", 
       device = tiff(), units = "cm",
       height = 5, width = 13, 
       dpi = 300,
       scale = 1.7)
