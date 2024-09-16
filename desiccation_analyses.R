#No evidence for the melanin desiccation hypothesis in a larval Lepidopteran
#Sarah Britton and Goggy Davidowitz

# Libraries
library(ggplot2)
library(rlang)
library(ggpubr)
library(car)
library(piecewiseSEM)
library(multcompView)
library(DiagrammeR)
library(dplyr)

#read in data
desiccation_data <- read.csv(file="Desiccation_Data_F22.csv")
image_data <- read.csv(file="Desiccation_Photo_Data_F23.csv")
vapo_data <- read.csv(file="Vapometer_Data_S24.csv")

#combine data sets/ clean data
combo_data <-left_join(desiccation_data, image_data, by = "ID")

desiccation_clean <- combo_data %>% 
  mutate(percent_mass_change = (100*percent_mass_change)) %>% 
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(log_SA = log(SA))%>%
  mutate(log_percent = log(percent_G))

#Settings for plots
treatment_labels = c("Crowded", "Solitary")  

#### Experiment 1 ####

##Figures

#Figure 2A
proportion_plot <- ggplot(desiccation_clean, aes(x=Treatment, y=percent_G)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, alpha=0.6) +
  stat_summary(aes(group=Treatment), fun=mean, shape="diamond", size=0.8) +
  geom_signif(comparison=list(c("C","S")), textsize=7, vjust=0.5, map_signif_level = TRUE)+
  scale_x_discrete(labels=treatment_labels) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Percentage melanic area (%)")
proportion_plot

#Figure 2B
darkness_plot <- ggplot(desiccation_clean, aes(Treatment, darkness_G)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, alpha=0.6) +
  stat_summary(aes(group=Treatment), fun=mean, shape="diamond", size=0.8) +
  geom_signif(comparison=list(c("C","S")), textsize=7, vjust=0.5, map_signif_level = TRUE)+
  scale_x_discrete(labels=treatment_labels) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Darkness")
darkness_plot

#Figure 2 combined 
ggarrange(proportion_plot, darkness_plot, 
          font.label = list(size=12, family="Times New Roman"),labels=c("A", "B"),
          ncol = 2, hjust=-7, align="hv", widths=c(2,2))

#Figure 3A
mass_change_1<-ggplot(desiccation_clean, aes(Treatment, percent_mass_change)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, alpha=0.6) +
  stat_summary(aes(group=Treatment), fun=mean, shape="diamond", size=0.8) +
  geom_signif(comparison=list(c("C","S")), textsize=7, vjust=0.5, map_signif_level = TRUE)+
  scale_x_discrete(labels=treatment_labels) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Percentage water loss (%)")+
  ylim(0,37)
mass_change_1

#Figure 4A
osmo_1 <- ggplot(desiccation_clean, aes(Treatment, delta_osmo)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, alpha=0.6) +
  stat_summary(aes(group=Treatment), fun=mean, shape="diamond", size=0.8) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.7, color="black") +
  geom_signif(comparison=list(c("C","S")), textsize=3, map_signif_level = TRUE)+
  scale_x_discrete(labels=treatment_labels) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Osmolality change")
osmo_1

#Other Osmo
osmo_2 <- ggplot(desiccation_clean, aes(x=before_mass, y=delta_osmo)) + 
  geom_point(aes(color=Treatment), show.legend = FALSE) +
  geom_smooth(method=lm, alpha=0.2) +
  scale_color_manual(values=colors, labels=treatment_labels) +
  geom_hline(yintercept=0, linetype="dashed", size=0.7, color="black") +
  xlim(1,7) +
  xlab("Body mass (g)") + ylab("Osmolality change")
osmo_2

##Statistics

#Summary stats
desiccation_clean %>%
  group_by(Treatment) %>%
  summarize(Percent_mean = mean(percent_G, na.rm = TRUE),
            Percent_sd = sd(percent_G, na.rm = TRUE),
            Darkness_mean = mean(darkness_G, na.rm = TRUE),
            Darkness_sd = sd(darkness_G, na.rm = TRUE),
            Change_mean = mean(percent_mass_change, na.rm = TRUE),
            Change_sd = sd(percent_mass_change, na.rm = TRUE),
            SA_mean = mean(SA, na.rm = TRUE),
            SA_sd = sd(SA, na.rm = TRUE),
            Osmo_mean = mean(delta_osmo, na.rm = TRUE),
            Osmo_sd = sd(delta_osmo, na.rm = TRUE),
            Mass_mean = mean(before_mass, na.rm = TRUE),
            Mass_sd = sd(before_mass, na.rm = TRUE))
as.data.frame

#Treatment and melanin
melanin_mod <- lm(percent_G~ Treatment, data=desiccation_clean)
summary(melanin_mod)
plot(melanin_mod)
ncvTest(melanin_mod)

darkness_mod <- lm(darkness_G ~ Treatment, data=desiccation_clean)
summary(darkness_mod)
plot(darkness_mod)
ncvTest(darkness_mod)

#Body size and melanin
#no size or SA differences between treatments
size_mod <- lm(before_mass ~ Treatment, data=desiccation_clean)
summary(size_mod)
plot(size_mod)
ncvTest(size_mod)

SA_mod <- lm(SA ~ Treatment, data=desiccation_clean)
summary(SA_mod)
plot(SA_mod)
ncvTest(SA_mod)

melanin_mass <- lm(percent_G ~ SA, data=desiccation_clean)
summary(melanin_mass)
plot(melanin_mass)
ncvTest(melanin_mass)

darkness_mass <- lm(darkness_G ~ SA, data=desiccation_clean)
summary(darkness_mass)
ncvTest(darkness_mass)

#Frass
frass_mod <- lm(frass_mass ~ Treatment, data=desiccation_clean)
summary(frass_mod)
qqnorm(residuals(frass_mod)) #normality of residuals
plot(frass_mod)
ncvTest(frass_mod)

#LM- Percent mass change
mass_loss_mod<- lm(percent_mass_change ~ Treatment + SA, data=desiccation_clean)
summary(mass_loss_mod)
plot(mass_loss_mod)
ncvTest(mass_loss_mod)

#SEM- Percent mass change 
sem_model <- psem(
  lm(percent_G ~ Treatment_Binary + SA, data=desiccation_clean, na.action=na.omit),
  lm(percent_mass_change ~ Treatment_Binary + percent_G + SA, desiccation_clean, na.action=na.omit),
  lm(SA~Treatment_Binary, data=desiccation_clean, na.action=na.omit)
)
summary(sem_model)

plot(sem_model, node_attrs = list(
  shape = "rectangle", color = "black",
  fillcolor = "gray"))

#LM- Osmolality change
osmo_change<- lm(log(delta_osmo) ~ Treatment + SA, data=desiccation_clean, na.action=na.omit)
summary(osmo_change)
plot(osmo_change)
ncvTest(osmo_change)

mass_osmo <- lm(delta_osmo ~ before_mass, data=desiccation_clean)
summary(mass_osmo)

osmo <- lm(percent_mass_change ~ delta_osmo, data=desiccation_clean)
summary(osmo)

leveneTest(delta_osmo ~ Treatment, data = desiccation_clean)
summary(variance)

sem_model_2 <- psem(
  lm(percent_G ~ Treatment_Binary + SA, data=desiccation_clean, na.action=na.omit),
  lm(delta_osmo ~ Treatment_Binary + percent_G + SA, desiccation_clean, na.action=na.omit),
  lm(SA~Treatment_Binary, data=desiccation_clean, na.action=na.omit)
)
summary(sem_model_2)

#split up data set for osmolality t-tests
osmo_crowded <- desiccation_clean %>% 
  filter(Treatment == "C" & delta_osmo != 0)

osmo_solo <- desiccation_clean %>% 
  filter(Treatment == "S" & delta_osmo != 0)

#T-tests: is change different from 0?
t.test(osmo_crowded$delta_osmo, mu = 0, alternative = "two.sided")
#p=0.063

t.test(osmo_solo$delta_osmo, mu = 0, alternative = "two.sided")
#p=0.907

#### Experiment 2 ####

#Figure 5A
ggplot(vapo_data, aes(treatment, avg_reading)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, alpha=0.6) +
  stat_summary(aes(group=treatment), fun=mean, shape="diamond", size=0.8) +
  geom_signif(comparison=list(c("crowded","solitary")), textsize=3, map_signif_level = TRUE) +
  scale_x_discrete(labels=treatment_labels) +
  theme_classic(base_size = 20)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + labs(y=expression("Evaporation rate (g/" *m^2*"h)"))

#Summary stats
vapo_data %>%
  group_by(treatment) %>%
  summarize(Percent_mean = mean(area_G, na.rm = TRUE),
            Percent_sd = sd(area_G, na.rm = TRUE),
            Darkness_mean = mean(darkness_G, na.rm = TRUE),
            Darkness_sd = sd(darkness_G, na.rm = TRUE),
            Evap_mean = mean(avg_reading, na.rm = TRUE),
            Evap_sd = sd(avg_reading, na.rm = TRUE),
            Mass_mean = mean(mass, na.rm = TRUE),
            Mass_sd = sd(mass, na.rm = TRUE)) %>%
  as.data.frame

#LM- evaporation                       
vapo_mod<- lm(avg_reading ~ treatment + ambient_humidity + body_temp, data=vapo_data)
summary(vapo_mod)

#Melanin models
vapo_melanin <- lm(area_G~ treatment, data=vapo_data)
summary(vapo_melanin)

vapo_darkness <- lm(darkness_G~ treatment, data=vapo_data)
summary(vapo_darkness)

#SEM- evaporation
sem_model_3 <- psem(
  lm(area_G ~ treatment_binary, data=vapo_data, na.action=na.omit),
  lm(avg_reading ~ treatment_binary + area_G, vapo_data, na.action=na.omit)
)
summary(sem_model_3)



