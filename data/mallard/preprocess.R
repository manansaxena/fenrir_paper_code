# set the path where all the R libraries are stored.
.libPaths(c("random_path", .libPaths()))
library(tidyverse)
library(lubridate)
setwd("random_path_where_code_is_stored/fenrir_paper_code/data/mallard/")

load("mallard_family.RData")

mallard_family$sample_data <- mallard_family$sample_data %>%
  select(X.SampleID,time,Vessel,SampleType,batch)%>%
  group_by(Vessel)%>%
  arrange(time)%>%
  mutate(time_num = as.numeric(time),
         Hour = (time_num - min(time_num))/3600,
         SampleType = ifelse(time_num == "1448031600","Hourly",SampleType),
         SampleType = ifelse(time_num %in% c("1448550000","1448636400","1448722800","1448809200"),"Daily",SampleType),
         Hour_diff = (Hour-mean(Hour))/sd(Hour),
         Hour_sample = ifelse(SampleType == "Hourly",1,0),
         Daily_sample = ifelse(SampleType == "Daily",1,0),
  )%>%
  filter(!duplicated(Hour_diff))%>%
  arrange(Vessel)

order_vector <- mallard_family[["sample_data"]][["X.SampleID"]]

Y_obs <- mallard_family[["otu_table"]][order_vector, ,drop = FALSE]

time_data <- as.POSIXct(mallard_family[["sample_data"]][["time"]], tz = "UTC")

start_time_v1 <- floor_date(min(time_data[1:135]), "hour")
end_time_v1 <- ceiling_date(max(time_data[1:135]), "hour")
all_hours_v1 <- seq(from = start_time_v1, to = end_time_v1, by = "hour")
observed_TT_v1 <- all_hours_v1 %in% time_data[1:135]

start_time_v2 <- floor_date(min(time_data[136:268]), "hour")
end_time_v2 <- ceiling_date(max(time_data[136:368]), "hour")
all_hours_v2 <- seq(from = start_time_v2, to = end_time_v2, by = "hour")
observed_TT_v2 <- all_hours_v2 %in% time_data[136:268]

start_time_v3 <- floor_date(min(time_data[269:401]), "hour")
end_time_v3 <- ceiling_date(max(time_data[269:401]), "hour")
all_hours_v3 <- seq(from = start_time_v3, to = end_time_v3, by = "hour")
observed_TT_v3 <- all_hours_v3 %in% time_data[269:401]

start_time_v4 <- floor_date(min(time_data[402:537]), "hour")
end_time_v4 <- ceiling_date(max(time_data[402:537]), "hour")
all_hours_v4 <- seq(from = start_time_v4, to = end_time_v4, by = "hour")
observed_TT_v4 <- all_hours_v4 %in% time_data[402:537]

observed_TT <- as.vector(c(observed_TT_v1,observed_TT_v2,observed_TT_v3,observed_TT_v4))

data <- list(Y_obs = Y_obs, observed_TT = observed_TT, N_obs_list = c(135,133,133,136))

saveRDS(data, file="data.rds")
write.table(observed_TT, "observed_indices.csv", row.names = FALSE, col.names = FALSE, sep = ",")
write.table(Y_obs, "Y_obs.csv", row.names = FALSE, col.names = FALSE, sep = ",")