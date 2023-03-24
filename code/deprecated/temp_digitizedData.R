


library(tidyverse); library(brms); theme_set(theme_bw())
source("code/002_fnDemography.R")

digitizedMetadata <- read_csv("..\\data\\digitized\\metadata.csv", show_col_types=F)

fig.f <- unique(str_remove(dir("..\\data\\digitized", "_fig.*csv"), "[a-z]?.csv"))

digitized <- fig.f %>% 
  map(~dir("..\\data\\digitized", paste0(.x, "[a-z]?.csv")) %>%
        map_dfr(~read_csv(paste0("..\\data\\digitized\\", .x), show_col_types=F) %>%
                  mutate(reference=str_remove(.x, ".csv"))) %>%
        left_join(., digitizedMetadata, by="reference")) %>%
  setNames(fig.f)

fig_var <- imap_dfr(digitized, 
                    ~tibble(V1=names(.x)[1], V2=names(.x)[2], reference=.y)) %>%
  arrange(V1, V2) %>%
  mutate(vars=paste0(V1, "_", V2))

unique_comparisons <- unique(fig_var$vars)
varCombined.ls <- vector("list", length(unique_comparisons)) %>% 
  setNames(unique_comparisons)
for(i in seq_along(varCombined.ls)) {
  refs <- fig_var$reference[fig_var$vars==unique_comparisons[i]]
  varCombined.ls[[i]] <- do.call('rbind', digitized[refs])
}


varCombined.ls$ageMin_weightStipe %>%
  ggplot(aes(factor(ageMin), weightStipe, fill=factor(depth))) + geom_boxplot()

# Within forest stands, size distribution is more bimodal (F2) vs. at patch limit (L5)
varCombined.ls$lengthStipe_proportion %>%
  ggplot(aes(lengthStipe, ymax=proportion, fill=location, group=reference)) + 
  geom_ribbon(ymin=0, alpha=0.5) + facet_grid(depth~.)


varCombined.ls$day_growthStipe %>%
  mutate(dateMeasured=lubridate::as_date(round(day))) %>%
  ggplot(aes(dateMeasured, growthStipe, colour=ageGroup, group=ageGroup)) + 
  geom_line() + scale_colour_viridis_c()
varCombined.ls$day_growthFrond %>%
  mutate(dateMeasured=lubridate::as_date(round(day))) %>%
  ggplot(aes(dateMeasured, growthFrond, colour=ageGroup, group=ageGroup)) + 
  geom_line() + scale_colour_viridis_c()
varCombined.ls$day_weightPlant %>%
  ggplot(aes(log(day), log(weightPlant), colour=depth, linetype=habitat)) + 
  geom_point() + stat_smooth(method="lm", se=F)
varCombined.ls$day_weightPlant %>%
  ggplot(aes(day, weightPlant, colour=depth, linetype=habitat)) + 
  geom_point()
varCombined.ls$day_areaFrond %>%
  mutate(dateMeasured=lubridate::as_date(round(day))) %>%
  ggplot(aes(dateMeasured, areaFrond, colour=age, group=reference)) + geom_line()

varCombined.ls$day_areaFrond %>%
  group_by(age, location) %>%
  mutate(deltaArea=areaFrond-lag(areaFrond),
         deltaDay=day-lag(day),
         rate=deltaArea/deltaDay,
         midptDay=day+deltaDay/2) %>%
  mutate(dateMeasured=lubridate::as_date(round(midptDay))) %>%
  ggplot(aes(dateMeasured, rate, colour=factor(age), linetype=location)) + 
  geom_line(size=1) + scale_colour_viridis_d("Age") +
  labs(x="Date", y=expression('Growth rate in frond area (cm'^2~'/day)'))
  



varCombined.ls$lengthStipe_weightStipe %>%
  ggplot(aes(log(lengthStipe), log(weightStipe), fill=reference)) + 
  geom_area()



varCombined.ls$weightStipe_weightFrond %>%
  ggplot(aes(weightStipe, weightFrond, colour=habitat)) + geom_point()
varCombined.ls$weightStipe_weightFrond %>%
  ggplot(aes(log(weightStipe), log(weightFrond), colour=habitat)) +
  geom_point() + facet_wrap(~round(age))
varCombined.ls$weightStipe_weightFrond %>%
  ggplot(aes(log(weightStipe), log(weightFrond), colour=habitat=="clearing")) +
  geom_point() + stat_smooth(method="lm")
varCombined.ls$weightStipe_weightFrond %>% filter(habitat=="limit") %>%
  ggplot(aes(log(weightStipe), log(weightFrond), colour=depth, group=depth)) +
  geom_point() + stat_smooth(method="lm")
varCombined.ls$weightStipe_weightFrond %>% filter(habitat!="tbd") %>%
  ggplot(aes(log(weightStipe), log(weightFrond), colour=depth, group=depth)) +
  geom_point() + stat_smooth(method="lm")




lenStipe_wtStipe <- brm(log(weightStipe) ~ log(lengthStipe),
                        data=varCombined.ls$lengthStipe_weightStipe)
wtStipe_wtFrond <- brm(log(weightFrond) ~ log(weightStipe), 
                       data=varCombined.ls$weightStipe_weightFrond)

lenStipe_wtStipe.mx <- as_draws_matrix(lenStipe_wtStipe)
wtStipe_wtFrond.mx <- as_draws_matrix(wtStipe_wtFrond)




plot(NA, NA, xlim=c(4, 7.2), ylim=c(-2,8), xlab="Stipe (ln mm)", ylab="Stipe (ln g)") 
for(i in 1:nrow(lenStipe_wtStipe.mx)) {
  simStipe <- exp(runif(10, 
                        log(min(varCombined.ls$lengthStipe_weightStipe$lengthStipe)),
                        log(max(varCombined.ls$lengthStipe_weightStipe$lengthStipe))))
  simStipeWt <- convertAllometry1D(log(simStipe), lenStipe_wtStipe.mx[i,1:2],
                                   shape="linear", 
                                   resid_err=lenStipe_wtStipe.mx[i,3])
  # points(log(simStipe), simStipeWt,
  #        col=rgb(.109375,0.5625,0.5976562,0.1), pch=1)
  abline(lm(simStipeWt ~ log(simStipe)), col=rgb(.109375,0.5625,0.5976562,0.1))
}
points(log(varCombined.ls$lengthStipe_weightStipe$lengthStipe), 
       log(varCombined.ls$lengthStipe_weightStipe$weightStipe))
abline(a=summary(lenStipe_wtStipe)$fixed[1,1], b=summary(lenStipe_wtStipe)$fixed[2,1])



plot(NA, NA, xlim=c(-3, 6), ylim=c(-4, 8), xlab="Stipe (ln g)", ylab="Frond (ln g)") 
for(i in 1:nrow(wtStipe_wtFrond.mx)) {
  simStipe <- exp(runif(10, 
                        log(min(varCombined.ls$weightStipe_weightFrond$weightStipe)),
                        log(max(varCombined.ls$weightStipe_weightFrond$weightStipe))))
  simFrondWt <- convertAllometry1D(log(simStipe), wtStipe_wtFrond.mx[i,1:2],
                                    shape="linear", 
                                    resid_err=wtStipe_wtFrond.mx[i,3])
    # points(log(simStipe), simFrondLen,
    #      col=rgb(.109375,0.5625,0.5976562,0.3), pch=1)
  abline(lm(simFrondWt ~ log(simStipe)), col=rgb(.109375,0.5625,0.5976562,0.1))
}
points(log(varCombined.ls$weightStipe_weightFrond$weightStipe), 
       log(varCombined.ls$weightStipe_weightFrond$weightFrond))
abline(a=summary(wtStipe_wtFrond)$fixed[1,1], b=summary(wtStipe_wtFrond)$fixed[2,1])



plot(NA, NA, xlim=c(4, 7.2), ylim=c(-2,8), xlab="Stipe (ln mm)", ylab="Frond (ln g)") 
for(i in 1:nrow(wtStipe_wtFrond.mx)) {
  simStipeLen <- exp(runif(100, 
                           log(min(varCombined.ls$lengthStipe_weightStipe$lengthStipe)),
                           log(max(varCombined.ls$lengthStipe_weightStipe$lengthStipe))))
  simStipeWt <- convertAllometry1D(log(simStipeLen), lenStipe_wtStipe.mx[i,1:2],
                                   shape="linear", 
                                   resid_err=lenStipe_wtStipe.mx[i,3])
  simFrondWt <- convertAllometry1D(simStipeWt, wtStipe_wtFrond.mx[i,1:2],
                                   shape="linear", 
                                   resid_err=wtStipe_wtFrond.mx[i,3])
  # points(log(simStipeLen), 
  #        simFrondWt,
  #        col=rgb(.109375,0.5625,0.5976562,0.3), pch=1)
  abline(lm(simFrondWt ~ log(simStipeLen)), col=rgb(.109375,0.5625,0.5976562,0.1))
}


