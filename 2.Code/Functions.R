#########################################
####              ####
#########################################

# Description
# Insert Description Here

# Author:
# Taylor Spaulding
# tspaulding@esassoc.com

# Date: 

# Read the original file

paths <- list.files(here::here("1.Data","raw"),full.names = TRUE)

read_logger_file <- function(file_path){
  raw <- data.frame(read.delim(file = file_path,
                               skip = 0,
                               header = FALSE))
  colnames(raw) <- "lines"
  raw$commas = lengths(regmatches(raw$lines, gregexpr(",",raw$lines)))
  raw <- raw[raw$commas > 4,]
  headers <- raw$lines[1]
  raw <- raw[-1,]

  ifelse(raw$commas[1] == 5, {
    raw <- tidyr::separate(data = raw,
                           col = lines,
                           sep = ",",
                           into = c(NA,'DateTime','uS','Temp', NA,
                                    NA, NA, NA))
    raw$DateTime <- lubridate::mdy_hms(raw$DateTime, tz = "GMT")
    raw$uS <- as.numeric(raw$uS)
    raw$Temp <- as.numeric(raw$Temp)
    raw <- dplyr::select(raw, DateTime, uS, Temp)
  },
  {
    raw <- tidyr::separate(data = raw,
                           col = lines,
                           sep = ",",
                           into = c('Date','Time','ms','level','Temp', 'uS',
                                    NA, NA, NA))
    raw$Date <- lubridate::as_date(raw$Date, tz = "GMT")
    raw$Time <- lubridate::hms(raw$Time)
    raw$DateTime = raw$Date + raw$DateTime
    raw$uS <- as.numeric(raw$uS)
    raw$Temp <- as.numeric(raw$Temp)
    raw <- dplyr::select(raw, DateTime, uS, Temp)
    })
raw
}

file_list <- list()
for(i in 1:length(paths)){
  file_list[[i]] <- read_logger_file(paths[i])
}

# Remove the start and end data -------------------------------------------
# Data contains a lot of extraneous data at the start and the end. 
# Using "landmarks" begin to hone in on the area of interest

# Calculate Instantaneous slope to get in-water and out-of-water points
inst_slope <- function(data,x_var,y_var){
  data <- data |>
    dplyr::mutate(slope = ({{x_var}}-dplyr::lag({{x_var}}))/
             as.numeric(({{y_var}}-dplyr::lag({{y_var}}))))
  data
}

for(i in 1:length(file_list)){
  file_list[[i]] <- inst_slope(file_list[[i]],uS,DateTime)
}

# Filter out points before and after
filter_extra <- function(data, slope_index){
  data[-c(1:which(data[,slope_index] > 100),
          which(data[,slope_index]<(-100)):length(data[,slope_index])),]
}

filtered_list <- list()
for(i in 1:length(file_list)){
  filtered_list[[i]] <- filter_extra(file_list[[i]],4)
}


# Model the data and clip to curve ---------------------------------------

model_clip <- function(data){
  data$num_DT = as.numeric(data$DateTime)
  model <- mgcv::gam(uS~s(num_DT), data = data)
  fit <- mgcv::predict.gam(model, data)
  data$fit <- fit
  
  data <- inst_slope(data, fit, num_DT)
  
  clip_DT <- data$DateTime[which(data$slope > -0.00002 & 
                                     data$slope < 0.00002)][c(1,3)]
  
  clipped <- data[data$DateTime > clip_DT[1] & 
                      data$DateTime < clip_DT[2],]
  clipped
}

clipped_list <- list()
for(i in 1:length(filtered_list)){
  clipped_list[[i]] <- model_clip(filtered_list[[i]])
}


# Isolate the background curve --------------------------------------------

isolate_slug <- function(data){
  clip_model <- loess(uS ~ num_DT, span = 0.01, data = data)
  clip_fit <- predict(clip_model, newdata = data)
  
  data$fit = clip_fit
  
  data <- inst_slope(data, fit, num_DT)
  
  index <- which(data$slope > 1.8|data$slope < -0.5)
  first_last <- c((min(index)-20):(max(index)+500))
  slug <- data[first_last,]
  background <- data[-first_last,]
  
  print(
    ggplot2::ggplot()+
      ggplot2::geom_point(data = slug, 
                          ggplot2::aes(x = DateTime, y = uS), 
                          color = "blue")+
      ggplot2::geom_point(data = background, 
                          ggplot2::aes(x = DateTime, y = uS), 
                          color = "green")+
      ggplot2::geom_smooth(data = background, 
                           ggplot2::aes(x = DateTime, y = uS),
                  method = "loess",
                  span = 1,
                  size = 1.25,
                  color = "black")+
      ggplot2::geom_smooth(data = slug, 
                  ggplot2::aes(x = DateTime, y = uS),
                  method = "loess",
                  span = 0.01,
                  size = 1.25,
                  color = "red"))
  ggplot2::ggsave(here::here("1.Data",
                             "Outputs",
                             "Figures",
                             paste0("File_",i,".pdf")), 
                  device = "pdf", width = 8, height = 6)
  
  slug_fit <- predict(loess(uS ~ num_DT, 
                            data = slug, 
                            span = 0.01), 
                      newdata = slug)
  baseline_fit <- predict(loess(uS ~ num_DT, 
                                data = background, 
                                span = 0.01), 
                          newdata = slug)
  
  slug$raw_uS = as.numeric(slug_fit)
  slug$background_uS = as.numeric(baseline_fit)
  slug$corrected_fit = slug$raw_uS - slug$background_uS
  slug <- slug[slug$corrected_fit > 0,]
  slug <- dplyr::select(slug,DateTime,corrected_fit)
  slug
}

slug_list <- list()
for(i in 1:length(clipped_list)){
  slug_list[[i]] <- isolate_slug(clipped_list[[i]])
}

# BTC Analysis ----------------------------------------------------
injection_time <- lubridate::mdy_hms("05/10/22 19:05:00", tz = "GMT")
slug_mass <- 63.07


analyze_btc <- function(data, injection_time, slug_mass){
  btc <- data
  btc$time_from_inj = as.numeric(difftime(data$DateTime, 
                                          injection_time, 
                                          units = "secs"))
  btc$conc = btc$corrected_fit/2
  btc$dt = as.numeric(btc$time_from_inj-dplyr::lag(btc$time_from_inj))
  btc$conc_t = ((btc$corrected_fit + lag(btc$corrected_fit))/2)*btc$dt
  
  #set the first value to 0
  btc$conc_t[1] <- 0
  
  sum <- sum(btc$conc_t)
  
  Q_lps <- (slug_mass/(sum/1000))
  Q_cfs <- (Q_lps * 0.035314)
  
  data.frame("Q_lps" = Q_lps, "Q_cfs" = Q_cfs)
}


for(i in length(slug_list)){
  discharge <- analyze_btc(slug_list[[i]], injection_time[i], slug_mass[i])
  data.frame("file" = paths[i])
  output <- dplyr::bind_cols(data.frame("file" = paths[i]), discharge)
  output
}




    
    