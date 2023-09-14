### Preamble
calc$ann_final_cdr <- 1 - ((1 - internals$hcdr)^4)
calc$ann_final_cdr_agewise <- as.data.frame(cbind(1900:2099, rowMeans(calc$ann_final_cdr[,1:55]), rowMeans(calc$ann_final_cdr[,56:65]), rowMeans(calc$ann_final_cdr[,66:100]))[51:200,])
colnames(calc$ann_final_cdr_agewise) <- c("Year", "1_C+A", "2_OA", "3_E")

gather(calc$ann_final_cdr_agewise, key=AgeGrp, value=CDR, 2:4) %>%
  ggplot() +
  geom_line(aes(x=Year, y=CDR, color=AgeGrp)) +
  labs(title='Case Detection Rate', subtitle='Scaling factors applied; model internal representation') +
  scale_y_continuous(limits = c(0,NA))

calc$ann_final_kappa <- 1 - ((1 - internals$hkappa)^4)
calc$ann_final_kappa_agewise <- as.data.frame(cbind(1900:2099, rowMeans(calc$ann_final_kappa[,1:55]), rowMeans(calc$ann_final_kappa[,56:65]), rowMeans(calc$ann_final_kappa[,66:100]))[51:200,])
colnames(calc$ann_final_kappa_agewise) <- c("Year", "1_C+A", "2_OA", "3_E")


gather(calc$ann_final_kappa_agewise, key=AgeGrp, value=TIR, 2:4) %>%
  ggplot() +
  geom_line(aes(x=Year, y=TIR, color=AgeGrp)) +
  labs(title='Treatment Initiation Rate', subtitle='Scaling factors applied; model internal representation') +
  scale_y_continuous(limits = c(0,NA))


#
# 1 - ((1 - internals$NTDS_II_n)^4)
# 1 - ((1 - internals$NTDS_II_u)^4)

# internals$NTDS_II_kappa
# internals$NTDS.CDR.II

# ((internals$NTDS_II_u + internals$NTDS_II_n) * internals$NTDS.CDR.II)/(1- internals$NTDS.CDR.II)
# ((internals$NTDS_IN_u + internals$NTDS_IN_n) * internals$NTDS.CDR.IN)/(1- internals$NTDS.CDR.IN)

# 1 - ((1 - internals$NTDS_IN_u)^4)

# calc$ann_final_kappa

# calc$ann_final_kappa_agewise
###

library(cowplot)
library(ggpubr)
# Prev_facet
prev <- calc$ann_ds_prev_norm_ii[,-c(2,3)] %>% gather(key = "AgeGrp", value = "Prevalence", -1) %>%
  ggplot() +
  geom_line(aes(x = Year, y = Prevalence)) +
  facet_wrap(~AgeGrp, nrow = 1) +
  geom_point(data = filter(gather(extdata$lancet_prevxage[,-7], key = "AgeGrp", value = "Prevalence", -c(1, 2)), Type == "M"), aes(x = Year, y = Prevalence)) +
  geom_errorbar(data = spread(gather(extdata$lancet_prevxage[,-7], key = "AgeGrp", value = "Prevalence", -c(1, 2)), key = Type, value = Prevalence), aes(x = Year, ymin = L, ymax = H)) +
  labs(title="AllTB Prevalence") +
  scale_y_continuous(name='Prevalence Rate', limits = c(0,200))


# Inc_facet
inc <- calc$ann_dstb_inc_norm[,-c(2,3)] %>%
  gather(key='AgeGrp', value='Incidence', -1) %>%
  ggplot() + geom_line(aes(x=Year, y=Incidence)) + facet_wrap(~AgeGrp, nrow=1) + labs(title = "AllTB Incidence") +
  scale_y_continuous(name='Incidence Rate')


# Notif_fact
notifR <- filter(extdata$alltb_notifR, AgeGrp!="All")
notif <- calc$ann_alltb_notif_norm[,-c(2,3)] %>%
  gather(key = "AgeGrp", value = "Notifications", 2:4) %>%
  ggplot() +
  geom_line(aes(x=Year, y=Notifications)) +
  facet_wrap(~AgeGrp, nrow = 1) +
  geom_point(data=notifR, aes(x=Year, y=Notifications)) +
  geom_errorbar(data=notifR, aes(x=Year, ymin=Notifications_lo, ymax=Notifications_hi)) +
  labs(title="AllTB Notifications") +
  scale_y_continuous(name='Notification Rate', limits = c(0,NA))

# CDR_facet
cdrf <- calc$ann_final_cdr_agewise %>%
  gather(key = "AgeGrp", value = "CDR", 2:4) %>%
  ggplot() +
  geom_line(aes(x=Year, y=CDR)) +
  facet_wrap(~AgeGrp, nrow=1) +
  labs(title="Case Detection Rate, annual") +
  scale_y_continuous(limits = c(0,NA))

# TIR facet
tirf <- calc$ann_final_kappa_agewise %>%
  gather(key = "AgeGrp", value = "TIRNR", 2:4) %>%
  ggplot() +
  geom_line(aes(x=Year, y=TIRNR)) +
  facet_wrap(~AgeGrp, nrow=1) +
  labs(title="Treatment Initiation Rate (~Notification Proportion), annual") +
  scale_y_continuous(limits = c(0,NA))

# Mortality
mort <- calc$ann_dstb_mort_norm[,-c(2,3)] %>%
  gather(key='AgeGrp', value='Mortality', -1) %>%
  ggplot() + geom_line(aes(x=Year, y=Mortality)) + facet_wrap(~AgeGrp, nrow = 1) + labs(title = "DSTB Mortality", subtitle = "Normalised to Age Groups") +
  scale_y_continuous(name = "Mortality (cases per 100,000)") +
  geom_point(data=filter(extdata$alltb_mortR, AgeGrp!="All"), aes(x=Year, y=Mortality)) +
  geom_errorbar(data=filter(extdata$alltb_mortR, AgeGrp!="All"), aes(x=Year, ymin=Mortality_lo, ymax=Mortality_hi))

# Demography
a <- qplot(calc$wth, rowSums(calc$ann_agewise[, 2:16]), geom = 'line', ylim = c(0,1E06), ylab = "Population (1000s)", xlab='Year')
b <- qplot(calc$wth, rowSums(calc$ann_agewise[, 17:65]), geom = 'line', ylim = c(0,1E06), ylab = "Population (1000s)", xlab='Year')
c <- qplot(calc$wth, rowSums(calc$ann_agewise[, 66:101]), geom = 'line', ylim = c(0,1E06), ylab = "Population (1000s)", xlab='Year')
ggarrange(a,b,c, nrow=1, ncol=3)

# Save plot
master <- ggarrange(prev, inc, mort, cdrf, tirf, notif, ggarrange(a,b,c, nrow=1, ncol=3), nrow=7, ncol=1)
ggsave(master, filename = "master.pdf", width = 14, height = 30, units = "in")

# New Demography

demog <- cbind(calc$wth, rowSums(calc$ann_agewise[, 2:16]), rowSums(calc$ann_agewise[, 17:65]), rowSums(calc$ann_agewise[, 66:101]))
colnames(demog) <- calc$incnames[-c(2,3)]
demog %>%
  gather(key=AgeGrp, value=Population) %>%
  ggplot() +
  geom_line(aes(x=Year, y=Population))


plot(calc$wth,rowSums(calc$ann_alltb_prev_ii[,c(66:101)])/rowSums(calc$ann_agewise[,66:101]))

