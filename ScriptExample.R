#script below was run for each Concentration-discharge station pair
library(tidyverse)
library(tidyhydat)
library(EGRET)
library(EGRETci)
library(dplyr)
library(fields)
library(survival)

#Import Q data from HYDAT
Qin<-hy_daily_flows(station_number = "",start_date="#########",end_date="######")
Qout<-hy_daily_flows(station_number = "",start_date="#####",end_date="########")

Qin<-data.frame(date=Qin$Date,Q=Qin$Value)
full_dates <- data.frame(date = seq(min(Qin$date), max(Qin$date), by = "day"))
Qin<-full_dates%>%
  left_join(Qin,by="date")
Qout<-data.frame(date=Qout$Date,Q=Qout$Value)
full_dates <- data.frame(date = seq(min(Qout$date), max(Qout$date), by = "day"))
Qout<-full_dates%>%
  left_join(Qout,by="date")
#gap filled or sum inflows from multiple sites if needed see Supplemental Material in Starling et al. (2026) 

#format Q for WRTDS
CSVFilenameQin <- paste0("./Q/","Qin.csv")
CSVFilenameQout <- paste0("./Q/","Qout.csv")

write.csv(Qin, file=paste0(CSVFilenameQin), row.names= FALSE)
write.csv(Qout, file=paste0(CSVFilenameQout), row.names= FALSE)
filepathQ<-""
Qin_WRTDS<-readUserDaily(filepathQ,filenameQin,qUnit=2)
Qout_WRTDS<-readUserDaily(filepathQ,filenameQout,qUnit=2)

#import concentration files (3 columns: Date (YYYY-MM-DD),comment (< if below detection limit), and value (concentration in mg/L))
filepathC<-""
TPin<-"TPin.csv" #or TN
TPinwrtds<-readUserSample(filepathC,TPin)
TPout<-"TPout.csv" #or TN
TPoutwrtds<-readUserSample(filepathC,TPout)

#write metadata for WRTDS/EGRET - must be run line by line
INFO_TPin<-readNWISInfo("","")
SITE NAME
NA
NA
SITE NAME
Total Phosphorus
TP
TP
mg/L

INFO_TPout<-readNWISInfo("","")
SITE NAME
NA
NA
SITE NAME
Total Phosphorus
TP
TP
mg/L

#calculate base WRTDS
ELIST_Pin<-mergeReport(INFO_TPin,Qin_WRTDS,TPinwrtds)
WRTDS_Pin<-modelEstimation(ELIST_Pin, windowY = 100, windowQ = 2, windowS = 0.5, minNumObs
                               = 30, minNumUncen =20, edgeAdjust = TRUE) 
ELIST_Pout<-mergeReport(INFO_TPout,Qout_WRTDS,TPoutwrtds)
WRTDS_Pout<-modelEstimation(ELIST_Pout, windowY = 100, windowQ = 2, windowS = 0.5, minNumObs
                           = 30, minNumUncen =20, edgeAdjust = TRUE)
#NSE
Sample<-getSample(ELIST_Pin) ##replace sample with each ELIST and run or loop
CrossV<-estCrossVal(2008,2011,Sample,windowY = 100, windowQ = 2, windowS = 0.5, minNumObs
                        = 30, minNumUncen =20, edgeAdjust = TRUE)
observed <- CrossV$ConcAve
modeled <- CrossV$ConcHat
# Calculate the mean of the observed values
mean_observed <- mean(observed)
# Calculate the numerator (sum of squared differences)
numerator <- sum((observed - modeled)^2)
# Calculate the denominator (sum of squared differences from the mean of observed values)
denominator <- sum((observed - mean_observed)^2)
# Calculate Nash-Sutcliffe Efficiency
NSE <- 1 - (numerator / denominator)
# Print the result
print(paste("Nash-Sutcliffe Efficiency (NSE):", round(NSE, 4)))

#RMSE and FBS
flux_error<-errorStats(WRTDS_Pout) 
flux_error<-errorStats(WRTDS_Nout) 
flux_error<-errorStats(WRTDS_Pin) 
flux_error<-errorStats(WRTDS_Nin) 

#Y-window =7 for Long term sites, set minobs proportonal to length of dataset and # of samples, values here used for short term (3 years) sites
#kalman filter and bootstrap
WRTDS_K_Pin <- WRTDSKalman(WRTDS_Pin)
dailyBoot_Pin <- genDailyBoot(WRTDS_K_inP,       # Generate a matrix of daily flux estimates. Rows = days, columns = individual estimates from bootstrapped WRTDS_K (e.g. nboot x nKalman)
                                nBoot   = 20, # The # of bootstrap replicates are performed, the WRTDS model is estimated for each one.
                                nKalman = 5,  # Number of times the Kalman process is performed on each WRTDS model
                                rho     = 0.9) 

annual_boot_Pin <- makeAnnualPI(dailyBoot_Pin,   # Using daily predictions of flux, calculate annual predictions of conc and flux at different percentiles. 
                                   WRTDS_K_Pin,
                                   paStart  = 10,  # Starting month e.g. 10 is water year, 1 is calendar year
                                   paLong   = 12,  # Number of months to consider e.g. 12 for a full year, could be less for seasonal estimates
                                   fluxUnit = 13)   # 3 = kg/day, 13 = kg/year, 8 = MTA

fluxCIinP <- annual_boot_Pin$flux 
concCIinP <- annual_boot_Pin$conc

WRTDS_K_Pout <- WRTDSKalman(WRTDS_Pout)
dailyBoot_Pout<- genDailyBoot(WRTDS_K_outP,       # Generate a matrix of daily flux estimates. Rows = days, columns = individual estimates from bootstrapped WRTDS_K (e.g. nboot x nKalman)
                                 nBoot   = 20, # The # of bootstrap replicates are performed, the WRTDS model is estimated for each one.
                                 nKalman = 5,  # Number of times the Kalman process is performed on each WRTDS model
                                 rho     = 0.9) 

annual_boot_Pout <- makeAnnualPI(dailyBoot_Pout,   # Using daily predictions of flux, calculate annual predictions of conc and flux at different percentiles. 
                                WRTDS_K_Pout,
                                paStart  = 10,  # Starting month e.g. 10 is water year, 1 is calendar year
                                paLong   = 12,  # Number of months to consider e.g. 12 for a full year, could be less for seasonal estimates
                                fluxUnit = 13)   # 3 = kg/day, 13 = kg/year, 8 = MTA

fluxCIoutP <- annual_boot_Pout$flux 
concCIoutP <- annual_boot_Pout$conc
#save results at 25th, 75th and 50th percentiles (alter to 5 and 95 if desired)
p50SITEPL<-data.frame(year=fluxCIinP$DecYear,
                          inflow=fluxCIinP$p50,
                          outflow=fluxCIoutP$p50)
p95SITEPL<-data.frame(year=fluxCIinP$DecYear,
                      inflow=fluxCIinP$p95,
                      outflow=fluxCIoutP$p95)
p75SITEPL<-data.frame(year=fluxCIinP$DecYear,
                      inflow=fluxCIinP$p75,
                      outflow=fluxCIoutP$p75)
p25SITEPL<-data.frame(year=fluxCIinP$DecYear,
                      inflow=fluxCIinP$p25,
                      outflow=fluxCIoutP$p25)
p5SITEPL<-data.frame(year=fluxCIinP$DecYear,
                      inflow=fluxCIinP$p5,
                      outflow=fluxCIoutP$p5)
p50SITEPC<-data.frame(year=concCIinP$DecYear,
                      inflow=concCIinP$p50,
                      outflow=concCIoutP$p50)
p95SITEPC<-data.frame(year=concCIinP$DecYear,
                      inflow=concCIinP$p95,
                      outflow=concCIoutP$p95)
p75SITEPC<-data.frame(year=concCIinP$DecYear,
                      inflow=concCIinP$p75,
                      outflow=concCIoutP$p75)
p25SITEPC<-data.frame(year=concCIinP$DecYear,
                      inflow=concCIinP$p25,
                      outflow=concCIoutP$p25)
p5SITEPC<-data.frame(year=concCIinP$DecYear,
                     inflow=concCIinP$p5,
                     outflow=concCIoutP$p5)
write.csv(p50SITEPL,file = "SITEPL.csv", row.names = FALSE)
write.csv(p50SITEPC,file = "SITEPC.csv", row.names = FALSE)

write.csv(p5SITEPL,file = "SITEPLp5.csv", row.names = FALSE)
write.csv(p5SITEPC,file = "SITEPCp5.csv", row.names = FALSE)

write.csv(p25SITEPL,file = "SITEPLp25.csv", row.names = FALSE)
write.csv(p25SITEPC,file = "SITEPCp25.csv", row.names = FALSE)

write.csv(p75SITEPL,file = "SITEPLp75.csv", row.names = FALSE)
write.csv(p75SITEPC,file = "SITEPCp75.csv", row.names = FALSE)

write.csv(p95SITEPL,file = "SITEPLp95.csv", row.names = FALSE)
write.csv(p95SITEPC,file = "SITEPCp95.csv", row.names = FALSE)

write.csv(dailyBoot_Pout,file="SITEdailyoutP.csv", row.names = FALSE)
write.csv(dailyBoot_Pin,file="SITEdailyinP.csv", row.names = FALSE)
#see README or paper supplement for sites with multiple inflow/outflow stations and weighting
#see corrections file for post-load calculation corrections and resulting Retention/Enrichent factor calculations

#B2 monthly
#must run script from Zhang et al. (2016) first
##Zhang, Q., Harman, C. J., & Ball, W. P. (2016).
##An improved method for interpretation of riverine concentration-discharge
##relationships indicates long-term shifts in reservoir sediment trapping. 
##Geophysical Research Letters, 43(19), 10,215-10,224. https://doi.org/10.1002/2016GL069945

B2_in_P<-modelEstimation1(ELIST_Pin, windowY = 100, windowQ = 2, windowS = 0.5, minNumObs
                             = 30, minNumUncen =20, edgeAdjust = TRUE)
B2_out_P<-modelEstimation1(ELIST_Pout, windowY = 100, windowQ = 2, windowS = 0.5, minNumObs
                          = 30, minNumUncen =20, edgeAdjust = TRUE)
#note set y-window to 7 for long term term sites and adjust minimum obs according to dataset
#pull daily results
inflow<-data.frame(date=B2_in_P[["Daily"]]$Date)
outflow<-data.frame(date=B2_out_P[["Daily"]]$Date)
inflow$B2_P<-B2_in_P[["Daily"]]$ConcDay
outflow$B2_P<-B2_out_P[["Daily"]]$ConcDay
#add water year
water_year <- function(dates) {
  dates <- as.Date(dates)
  yr <- as.numeric(format(dates, "%Y"))
  mo <- as.numeric(format(dates, "%m"))
  wy <- ifelse(mo >= 10, yr + 1, yr)
  return(wy)
}

inflow$wyear<-water_year(inflow$date)
outflow$wyear<-water_year(outflow$date)

#calculate average B2 by water year
B2_annual_in<-inflow%>%
  group_by(wyear)%>%
  summarize(
    avgB2_P=mean(B2_P),
    sdB2_P=sd(B2_P)
  )
B2_annual_out<-outflow%>%
  group_by(wyear)%>%
  summarize(
    avgB2_P=mean(B2_P),
    sdB2_P=sd(B2_P)
  )
#group by month to understand seasonal trends in B2 
B2_month_in_P<-inflow[["Daily"]]%>%
  group_by(Month)%>%
  summarize(
    B2avgP=mean(ConcDay),
    B2sdP=sd(ConcDay)
  )
B2_month_out_P<-outflow[["Daily"]]%>%
  group_by(Month)%>%
  summarize(
    B2avgP=mean(ConcDay),
    B2sdP=sd(ConcDay)
  )

#note - for sites with multiple inflows B2in is the weighted average of monthly avg. inflowing B2 by average annual fraction of inflow from each stream
#note - for Darling B2 of the outflow was corrected for tributary influence as B2_out<- B2_outA$avgB2_P*0.92-B2_tributary$avgB2_P*0.08

monthlyB2<-data.frame(month=B2_month_sher_P,
                    B2in_P=B2_month_in_P$B2avgP,
                    B2out_P=B2_month_out_P$B2avgP
)
annualB2<-data.frame(year=B2_annual_in$,
                     B2in_P=B2_annual_in$avgB2_P,
                     B2out_P=B2_annual_out$avgB2_P
) 
#save results to csv
write.csv(monthlyB2,"B2monthly.csv",row.names = FALSE)
write.csv(annualB2,"B2annual.csv",row.names = FALSE)

############NP molar ratio monthly########################
monthlyPin <- makeMonthPI(dailyBoot_Pin,WRTDS_Pin)
monthlyPout<- makeMonthPI(dailyBoot_Pout,WRTDS_Pout)
#have run both N and P
monthlyNin <- makeMonthPI(dailyBoot_Nin,WRTDS_Nin)
monthlyNout<- makeMonthPI(dailyBoot_Nout,WRTDS_Nout)

decimal_to_date <- function(decimal_year) {
  year <- floor(decimal_year)
  remainder <- decimal_year - year
  start <- as.Date(paste0(year, "-01-01"))
  end <- as.Date(paste0(year + 1, "-01-01"))
  days_in_year <- as.numeric(end - start)
  start + round(remainder * days_in_year)
}
#note for some sites N and P in will be the sum of multiple inflowing loads
monthlyinP <- monthlyPin$flux
monthlyoutP<- monthlyPout$flux
monthlyinN<- monthlyNin$flux
monthlyoutN<- monthlyNout$flux

monthlyinP$date<- decimal_to_date(as.numeric(monthlyinP$DecYear))
monthlyoutP$date<- decimal_to_date(as.numeric(monthlyoutP$DecYear))
monthlyinN$date<- decimal_to_date(as.numeric(monthlyinN$DecYear))
monthlyoutN$date<- decimal_to_date(as.numeric(monthlyoutN$DecYear))

inflowNP<-data.frame(
  date=monthlyinP$date,
  Pin=((monthlyinP$p50*1000)/30.974),
  Nin=((monthlyinN$p50*1000)/14.01)
)
outflowNP<-data.frame(
  date=monthlyMinotP$date,
  Pout=((monthlyoutP$p50*1000)/30.974),
  Nout=((monthlyoutN$p50*1000)/14.01)
)

merge<-inner_join(inflowNP,outflowNP,by="date")

merge$NPin<-merge$Nin/merge$Pin
merge$NPout<-merge$Nout/merge$Pout

write.csv(merge,"NPmontly.csv",row.names = FALSE)

#CVc/CVq Ratio of coefficient of variation of concentration vs coefficient of variation for discharge
#note for sites with multiple inflows Q is the daily sum of inflows and C is the weighted average based on each sites relative contribution to inflow over the study period
Q<-data.frame(Date=Qin$Date,Qin=Qin$Q,Qout=Qout$Q) #note that data should be cut to the overlap of the inflow(s) and outflow concentration record

concCIoutP<-dailyBoot_Pout$conc
concCIoutN<-dailyBoot_Nout$conc
concCIinP<-dailyBoot_Pin$conc
concCIinN<-dailyBoot_Nin$conc

#cut to same temporal range if needed (based on max intersection of all stations)
conc<-data.frame(Date=concCIinN$Date,
                 inP=concCIinP$p50,
                 inN=concCIinN$p50,
                 outP=concCIoutP$p50,
                 outN=concCIoutN$p50)
conc$year<-water_year(conc$Date)
CVC<-conc%>%
  group_by(year)%>%
  summarise(
    CVPin=sd(inP)/mean(inP),
    CVNin=sd(inN)/mean(inN),
    CVPout=sd(outP)/mean(outP),
    CVNout=sd(outN)/mean(outN)
  )

Q$year<-water_year(Q$Date)
CVQ<-Q%>%
  group_by(year)%>%
  summarize(CVqin=sd(inQ)/mean(inQ),
            CVqout=sd(outQ)/mean(inQ))

CV<-merge(CVQ,CVC,by="year")
CV$CVCVPin<-CV$CVPin/CV_GR$CVqin
CV$CVCVNin<-CV$CVNin/CV_GR$CVqin
CV$CVCVPout<-CV$CVPout/CV_GR$CVqout
CV$CVCVNout<-CV$CVNout/CV_GR$CVqout

write.csv(CV,file="CV.csv", row.names = FALSE)

#hydrological metrics of freshet timing and intensity
Qin<-data.frame(year=year(Q$Date),date=Q$Date,Q=Q$WFB) #dataframe of inflow discharge
#add column that is day into the water year
Qin <- Qin %>%
  mutate(
    water_year = if_else(month(date) >= 10, year(date) + 1, year(date)),
    water_day  = yday(date) - yday(as.Date(paste0(year(date), "-10-01"))) + 1
  )
Qin$water_day[Qin$water_day <= 0] <- 
  Qin$water_day[Qin$water_day <= 0] + 
  (365 + leap_year(Qin$date[Qin$water_day <= 0]))


#cummulative volume within each water year
Qin<- Qin %>%
  group_by(water_year) %>%
  arrange(date) %>%
  mutate(cum_volume = cumsum(Q* 86400))
#calculate half of total annual cummulative volume 
sumQ<-Qin%>%
  group_by(water_year)%>%
  filter(water_day == max(water_day)) %>%
  summarize(
    maxV=cum_volume,
    V50=(cum_volume/2)
  )

Q50 <- Qin %>%
  group_by(water_year) %>%
  # join the V50 value for that year
  left_join(sumQ %>% select(water_year, V50), by = "water_year") %>%
  # keep only rows where cum_volume >= V50
  filter(cum_volume >= V50) %>%
  # pick the first occurrence per year
  slice_min(order_by = cum_volume, n = 1) %>%
  ungroup() #this is the timing metric!! (# of days until 50% of water year flow)

#calculate change in cumulative volume each day
Qin<-Qin %>%
  mutate(
    v_after   = lead(cum_volume),  # value from next row (shift up)
    v_before = lag(cum_volume)    # value from previous row (shift down)
  )
Qin$m<-((((Qin$cum_volume-Qin$v_before))+
                ((Qin$v_after-Qin$cum_volume)))/2)
Q_list<-split(Qin,Qin$water_year)

#sort by m (most intensive flow)
sorted_m<-lapply(Q_list, function(df) {
    df %>% arrange(desc(m))
})
sorted_m <-lapply(sorted_m, function(df) {
    df %>% mutate(Volume = Q * 86400)
})

sorted_m <- lapply(sorted_m, function(df) {
    df %>% mutate(v50 = (sum(Volume)/2))
})
sorted_m <- lapply(sorted_m, function(df) {
    df %>% mutate(cVolume_m = cumsum(Volume))
})
##return the row number
sorted_m <- lapply(sorted_m, function(df) {
    df %>% mutate(
      index = which(cVolume_m>= v50)[1]
    )
})#this is the intensity (number of peak flow days to reach 50% of annual flow volume)

#merge and print results
result_df <- do.call(rbind, lapply(sorted_m, function(df) {
  data.frame(
    year = unique(df$water_year),
    m_days = df$index[1],
    V50_value = unique(df$v50)
  )
}))
result_df <- do.call(rbind, lapply(sorted_m, function(df) {
  data.frame(
    year = unique(df$water_year),
    m_days = df$index[1],
    V50 = unique(df$v50)
  )
}))
result_df <- merge(result_df, Q50[, c("year", "water_day")], by = "year", all.x = TRUE)

result_df <- result_df %>%
  rename(timing=water_day,
         intensity=m_days)

getwd()
setwd("D:/")
write.csv(result_df,"CumulativeFlowStats.csv")

#once all sites are run via Cummulative Flow Stats are in one excel file
library(readr)
library(tidyverse)
library(patchwork)

FlowStats <- read_csv("path to CumulativeFlowStats.csv")
FlowStats<- FlowStats %>%
  rename(
    intensity = m_days
  )

P<-FlowStats%>%
  filter(!is.na(Pret))

N<-FlowStats%>%
  filter(!is.na(Nret))

P <- P %>%
  group_by(site) %>%
  mutate(
    volume = case_when(
      V50 >= quantile(V50, 0.75, na.rm = TRUE) ~ "P75",
      V50 <= quantile(V50, 0.25, na.rm = TRUE) ~ "P25",
      TRUE                                     ~ "P5"
    )
  ) %>%
  ungroup()

N <- N %>%
  group_by(site) %>%
  mutate(
    volume = case_when(
      V50 >= quantile(V50, 0.75, na.rm = TRUE) ~ "P75",
      V50 <= quantile(V50, 0.25, na.rm = TRUE) ~ "P25",
      TRUE                                     ~ "P5"
    )
  ) %>%
  ungroup()


Pintensity<-P%>%
  ggplot(aes(intensity,Pret,colour=site))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE)+
  theme_bw()+
  ggtitle("hydrological intensity")+
  ylim(-100, 100)+
  theme(legend.position = "none")+
  ylab("REF TP (%)")

Nintensity<-N%>%
  ggplot(aes(intensity,Nret,colour=site))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE)+
  theme_bw()+
  ggtitle("hydrological intensity")+
  ylim(-100, 100)+
  theme(legend.position = "none")+
  ylab("REF TN (%)")

Ptiming<-P%>%
  ggplot(aes(timing,Pret,colour=site))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE)+
  theme_bw()+
  ggtitle("hydrological timing")+
  ylim(-100, 100)+
  theme(legend.position = "none")+
  ylab("REF TP (%)")

Ntiming<-N%>%
  ggplot(aes(timing,Nret,colour=site))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE)+
  ggtitle("hydrological timing")+
  theme_bw()+
  ylim(-100, 100)+
  theme(legend.position = "none")+
  ylab("REF TN (%)")

Nvolume<-N%>%
  ggplot(aes(volume,Nret,fill=site))+
  geom_boxplot()+
  ggtitle("relative volume of annual flows")+
  theme_bw()+
  ylim(-100,100)+
  ylab("REF TN (%)")

Nvolume
Pvolume<-P%>%
  ggplot(aes(volume,Pret,fill=site))+
  geom_boxplot()+
  ggtitle("relative volume of annual flows")+
  theme_bw()+
  ylim(-100,100)+
  ylab("REF TP (%)")

Pintensity+Ptiming+Nintensity+Ntiming

Pintensity+Ptiming+Pvolume+Nintensity+Ntiming+Nvolume

resultsP<-P%>%
  group_by(site)%>%
  nest()%>%
  mutate(
    model=map(data,~lm(Pret~timing+intensity+volume,data=.))
  )

coef_dfP <- resultsP %>%
  mutate(tidy_model = map(model, broom::tidy)) %>%  # extract coefficients
  unnest(tidy_model) %>%
  select(site, term, estimate) %>%
  mutate(term = recode(term,
                       "(Intercept)" = "intercept",
                       "timing" = "timing",
                       "intensity" = "intensity",
                       "volume" = "volume")) %>%
  pivot_wider(names_from = term, values_from = estimate)

coef_dfP
#ANOVA
resultsP <- P %>%
  group_by(site) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Pret ~ timing + intensity+volume, data = .)),
    anova_tbl = map(model, ~ anova(.x) %>% broom::tidy()),
  ) %>%
  unnest(anova_tbl) %>%
  group_by(site) %>%
  mutate(
    prop_var = sumsq / sum(sumsq)   # proportion variance explained
  ) %>%
  ungroup()

# 2. Site-level summary across all sites
summary_across_sitesP <- resultsP %>%
  #  filter(term != "Residuals") %>%
  group_by(term) %>%
  summarise(
    mean_prop = mean(prop_var, na.rm = TRUE),
    sd_prop   = sd(prop_var, na.rm = TRUE),
    mean_F    = mean(statistic, na.rm = TRUE),
    mean_p    = mean(p.value, na.rm = TRUE),
    n_sites   = n_distinct(site),
    .groups = "drop"
  )

Pres_exp<-data.frame(site=resultsP$site,
                     factor=resultsP$term,
                     p_value=resultsP$p.value,
                     prop_var=resultsP$prop_var)
view(Pres_exp)
getwd()
setwd("D:/MS1scripts/results")
write.csv(Pres_exp,"anova+bydro_wvol_P.csv")

#ANOVA N

resultsN<-N%>%
  group_by(site)%>%
  nest()%>%
  mutate(
    model=map(data,~lm(Nret~timing+intensity+volume,data=.))
  )

coef_dfN <- resultsN %>%
  mutate(tidy_model = map(model, broom::tidy)) %>%  # extract coefficients
  unnest(tidy_model) %>%
  select(site, term, estimate) %>%
  mutate(term = recode(term,
                       "(Intercept)" = "intercept",
                       "timing" = "timing",
                       "intensity" = "intensity",
                       "volume" = "volume")) %>%
  pivot_wider(names_from = term, values_from = estimate)

coef_dfN
#ANOVA
resultsN <- N %>%
  group_by(site) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Nret~ timing + intensity+volume, data = .)),
    anova_tbl = map(model, ~ anova(.x) %>% broom::tidy()),
  ) %>%
  unnest(anova_tbl) %>%
  group_by(site) %>%
  mutate(
    prop_var = sumsq / sum(sumsq)   # proportion variance explained
  ) %>%
  ungroup()

Nres_exp<-data.frame(site=resultsN$site,
                     factor=resultsN$term,
                     p_value=resultsN$p.value,
                     prop_var=resultsN$prop_var)
view(Nres_exp)
getwd()
setwd("D:/MS1scripts/results")
write.csv(Nres_exp,"anova+bydro_wvol_N.csv")

# 2. Site-level summary across all sites
summary_across_sitesN <- resultsN %>%
  #  filter(term != "Residuals") %>%
  group_by(term) %>%
  summarise(
    mean_prop = mean(prop_var, na.rm = TRUE),
    sd_prop   = sd(prop_var, na.rm = TRUE),
    mean_F    = mean(statistic, na.rm = TRUE),
    mean_p    = mean(p.value, na.rm = TRUE),
    n_sites   = n_distinct(site),
    .groups = "drop"
  )
