# Installing required packages
install.packages(c("ncdf4", "lubridate","tidyverse","readxl"))
# Extracting data from norkyst in order to model

library(ncdf4)
library(lubridate)
library(tidyverse)
library(readxl)

# Extract time series of temperature, salinity and currents from oper. Norkyst from single locations

# Path where main member of operational Norkyst file is located
URL <- "https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be"

# User definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define chosen locations where Norkyst data is extracted, list in vectors lonpos and latpos
lonpos <- numeric(3)
latpos <- numeric(3)
# mfs
lonpos[1] <- 6.83070
latpos[1] <- 58.08520          
#ystestein
lonpos[2] <- 6.87469
latpos[2] <- 58.01121 
# svartskjær roholmen 
lonpos[3] <- 6.85873
latpos[3] <- 58.06199         
# mfs 0 5 50m 17.04 - 29.04
# svartskjær 
# ystestein

# Requested time period (if you need more than ~15 days, do this in chunks to prevent a substantial delay, remember its hourly values)
dnum0 <- as.POSIXct("2026-04-28 08:00:00", tz = "UTC")
dnum1 <- as.POSIXct("2026-04-30 16:00:00", tz = "UTC")

# Choose what depth you want data from (available p.t. is 0, 1, 2, 3, 5, 7, 10, 15, 25, 50, 65, 75, 100, 200 and 300m)
zout <- c(0, 5, 10, 25, 50, 75, 100)

# End user definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Open the NetCDF file
nc <- nc_open(URL)

# Need some static parameters
LON <- ncvar_get(nc, "lon")
LAT <- ncvar_get(nc, "lat")
depth <- ncvar_get(nc, "depth")
TIDER <- ncvar_get(nc, "time") / 86400  # convert seconds to days

# Convert TIDER to POSIXct date-time (origin 1970-01-01)
dnum <- as.POSIXct(TIDER * 86400, origin = "1970-01-01", tz = "UTC")

# Find indices of requested time stamps
T <- which(dnum >= dnum0 & dnum <= dnum1)

# Find nearest grid points to chosen location(s)
xpos <- integer(length(lonpos))
ypos <- integer(length(lonpos))
for (p in seq_along(lonpos)) {
  mindist <- 1000
  for (i in seq_len(dim(LON)[1])) {
    for (j in seq_len(dim(LON)[2])) {
      dist <- sqrt((lonpos[p] - LON[i, j])^2 + (latpos[p] - LAT[i, j])^2)
      if (dist < mindist) {
        xpos[p] <- i
        ypos[p] <- j
        mindist <- dist
      }
    }
  }
}

# Find depth indices for requested depths
K <- integer(0)
for (k in seq_along(zout)) {
  n <- which(depth == zout[k])
  if (length(n) > 0) {
    K <- c(K, n)
  } else {
    message(sprintf("Chosen output depth %d is not available, script continues without that depth level", zout[k]))
  }
}

# Initialize arrays to store data
lenT <- length(T)
lenP <- length(xpos)
lenK <- length(K)

U <- array(NA_real_, dim = c(lenP, lenK, lenT))
V <- array(NA_real_, dim = c(lenP, lenK, lenT))
SALT <- array(NA_real_, dim = c(lenP, lenK, lenT))
TEMP <- array(NA_real_, dim = c(lenP, lenK, lenT))

# Read Norkyst data
for (p in seq_len(lenP)) {
  for (k in seq_len(lenK)) {
    start <- c(xpos[p], ypos[p], K[k], T[1])
    count <- c(1, 1, 1, lenT)
    U[p, k, ] <- ncvar_get(nc, "u_eastward", start = start, count = count)
    V[p, k, ] <- ncvar_get(nc, "v_northward", start = start, count = count)
    SALT[p, k, ] <- ncvar_get(nc, "salinity", start = start, count = count)
    TEMP[p, k, ] <- ncvar_get(nc, "temperature", start = start, count = count)
    message(sprintf("Variables from position no. %d and depth %dm is read from operational Norkyst file", p, zout[k]))
  }
}

nc_close(nc)

# Calculate speed from current vectors
SPD <- array(NA_real_, dim = c(lenP, lenK, lenT))
for (p in seq_len(lenP)) {
  for (k in seq_len(lenK)) {
    SPD[p, k, ] <- sqrt(U[p, k, ]^2 + V[p, k, ]^2)
  }
}

# Vectoring the variables and putting them into corresponding depth data frames

df_0m <- data.frame( 
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 1, ]), as.vector(TEMP[2, 1, ]),  as.vector(TEMP[3, 1, ]) ),
  current = c(as.vector(SPD[1, 1, ]), as.vector(SPD[2, 1, ]), as.vector(SPD[3, 1, ]) ),
  salinity = c(as.vector(SALT[1, 1, ]), as.vector(SALT[2, 1, ]), as.vector(SALT[3, 1, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT))

df_5m <- data.frame(
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 2, ]), as.vector(TEMP[2, 2, ]),  as.vector(TEMP[3, 2, ]) ), 
  current = c(as.vector(SPD[1, 2, ]), as.vector(SPD[2, 2, ]), as.vector(SPD[3, 2, ]) ),
  salinity = c(as.vector(SALT[1, 2, ]), as.vector(SALT[2, 2, ]), as.vector(SALT[3, 2, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT)
)

df_10m <- data.frame(
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 3, ]), as.vector(TEMP[2, 3, ]),  as.vector(TEMP[3, 3, ]) ), 
  current = c(as.vector(SPD[1, 3, ]), as.vector(SPD[2, 3, ]), as.vector(SPD[3, 3, ]) ),
  salinity = c(as.vector(SALT[1, 3, ]), as.vector(SALT[2, 3, ]), as.vector(SALT[3, 3, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT)
)

df_25m <- data.frame(
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 4, ]), as.vector(TEMP[2, 4, ]),  as.vector(TEMP[3, 4, ]) ), 
  current = c(as.vector(SPD[1, 4, ]), as.vector(SPD[2, 4, ]), as.vector(SPD[3, 4, ]) ),
  salinity = c(as.vector(SALT[1, 4, ]), as.vector(SALT[2, 4, ]), as.vector(SALT[3, 3, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT)
)

df_50m <- data.frame(
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 5, ]), as.vector(TEMP[2, 5, ]),  as.vector(TEMP[3, 5, ]) ), 
  current = c(as.vector(SPD[1, 5, ]), as.vector(SPD[2, 5, ]), as.vector(SPD[3, 5, ]) ),
  salinity = c(as.vector(SALT[1, 5, ]), as.vector(SALT[2, 5, ]), as.vector(SALT[3, 5, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT)
)

df_75m <- data.frame(
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 6, ]), as.vector(TEMP[2, 6, ]),  as.vector(TEMP[3, 6, ]) ), 
  current = c(as.vector(SPD[1, 6, ]), as.vector(SPD[2, 6, ]), as.vector(SPD[3, 6, ]) ),
  salinity = c(as.vector(SALT[1, 6, ]), as.vector(SALT[2, 6, ]), as.vector(SALT[3, 6, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT)
)

df_100m <- data.frame(
  datetime = rep(dnum[T], times = 3),
  temperature = c(as.vector(TEMP[1, 7, ]), as.vector(TEMP[2, 7, ]),  as.vector(TEMP[3, 7, ]) ), 
  current = c(as.vector(SPD[1, 7, ]), as.vector(SPD[2, 7, ]), as.vector(SPD[3, 7, ]) ),
  salinity = c(as.vector(SALT[1, 7, ]), as.vector(SALT[2, 7, ]), as.vector(SALT[3, 7, ]) ),
  location = rep(c("mfs", "yst", "sva"), each = lenT)
)

#averaging every vector for each station for the three days

df_0m_mean <- df_0m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)

df_5m_mean <- df_5m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)

df_10m_mean <- df_10m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)

df_25m_mean <- df_25m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)


df_50m_mean <- df_50m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)


df_75m_mean <- df_75m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)

df_100m_mean <- df_100m %>% 
  group_by(location) %>% 
  summarise( mean_current = mean(current, na.rm = TRUE),
             mean_temperature = mean(temperature, na.rm = TRUE),
             mean_salinity = mean(salinity, na.rm = TRUE),)

#combining them into a single data frame

df_list <-  list(
  "0" = df_0m_mean,
  "5" = df_5m_mean,
  "10" = df_10m_mean,
  "25" = df_25m_mean,
  "50" = df_50m_mean,
  "75" = df_75m_mean,
  "100" = df_100m_mean
)

df_modeldata  <- bind_rows(df_list, .id = "depth")


#-----------------------------------------------------

###Analyzing and wrangling protist data

library(tidyverse)

phyto_df <- read_excel("Phyto.xlsx") #or please insert the name of the plankton sample file here if "Phyto.xlsx" is not the name of the given file


#condensing the data set into usable variables

phytowr_df <- phyto_df %>%
  group_by(Station) %>%
  summarise(
    abundance = sum(Chain),     # total individuals each sample. variable went unused due to protist data being purely qualitative,
    richness = n_distinct(Family),           # number of families observed
  )

#creating a column of dates corresponding to the different sampling days and stations
dates <- c(
  MFS   = "2026-04-29",
  SVA     = "2026-04-28",
  YST = "2026-04-30"
)

#combining the wrangled phyto data frame and the dates

plankton_df <- phytowr_df %>% 
  mutate(date = case_when(
    Station == "MFS"   ~ as.Date("2026-04-29"),
    Station == "SVA"     ~ as.Date("2026-04-28"),
    Station == "YST" ~ as.Date("2026-04-30")
  )) %>% 
  select(date, everything()) 



#plotting richness

ggplot(plankton_df, aes(x = Station, y = richness)) +
  geom_col()+
  labs(title = "Protist taxon richness per station",
       x = "Station",
       y = "Richness")+
  scale_color_discrete(labels = c("MFS" = "Midtfjorskj\u00E6r",
                                  "SVA" = "Svartskj\u00E6r",
                                  "YST" = "Ystesteinen"))+
  theme_minimal()


#saving plot 
ggsave("planktonrichness.jpeg", plot = get_last_plot())





#end of script




