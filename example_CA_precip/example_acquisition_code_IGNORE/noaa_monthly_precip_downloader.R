library(rvest)
library(dplyr)
library(purrr)

###Scrape station locations----------------------------------------------------------------------
locations_url  <- "https://www.cnrfc.noaa.gov/monthly_precip_locations.php"
locations_page <- read_html(locations_url)
tables_nodes   <- html_elements(locations_page, "table")

locations_df <- map_dfr(tables_nodes, function(tbl) {
  df <- html_table(tbl, header = FALSE)
  if (nrow(df) < 3) return(NULL)
  colnames(df) <- unlist(df[2, ])
  df <- df[-(1:2), ]
  df
})

#simplify column names
colnames(locations_df)[colnames(locations_df) == "ID"]               <- "Station_ID"
colnames(locations_df)[colnames(locations_df) == "Latitude (°N)"]    <- "Latitude"
colnames(locations_df)[colnames(locations_df) == "Longitude (°W)"]   <- "Longitude"
colnames(locations_df)[colnames(locations_df) == "Elevation (ft)"]   <- "Elevation"

#make sure values are saved as numericals
locations_df$Latitude  <- as.numeric(locations_df$Latitude)
locations_df$Longitude <- -abs(as.numeric(locations_df$Longitude)) #must be converted to negative
locations_df$Elevation <- as.numeric(locations_df$Elevation)

###Scrape monthly precipitation table--------------------------------------------------------------
precip_url  <- "https://www.cnrfc.noaa.gov/monthly_precip_2024.php"
precip_page <- read_html(precip_url)
rows        <- html_elements(precip_page, "table tr")
#extract all rows as vector values
parsed <- lapply(rows, function(row) {
  html_text(html_elements(row, "th, td"), trim = TRUE)
})
#code to keep only rows with data; not headings or subheadings
row_lengths <- sapply(parsed, length)
keep        <- parsed[row_lengths %in% c(15, 17)]
header      <- keep[[1]]
n_cols      <- length(header)
#lose extra columns for those that have them
data_list <- lapply(keep[-1], function(r) {
  length(r) <- n_cols
  r
})

precip_df <- as.data.frame(do.call(rbind, data_list), stringsAsFactors = FALSE)
colnames(precip_df) <- header
colnames(precip_df)[1] <- "Station_ID"

#Convert all columns except Station_ID and Location to numeric
non_numeric_cols <- c("Station_ID", "Location")
for (col in colnames(precip_df)) {
  if (!col %in% non_numeric_cols) {
    precip_df[[col]] <- suppressWarnings(as.numeric(precip_df[[col]]))
  }
}

#Delete rows with NA
numeric_cols   <- setdiff(colnames(precip_df), non_numeric_cols)
has_na         <- apply(precip_df[, numeric_cols], 1, anyNA)
precip_df      <- precip_df[!has_na, ]
precip_df      <- distinct(precip_df)

####Join-----------------------------------------------------------------------------------------------
combined_df <- inner_join(locations_df, precip_df, by = "Station_ID")
combined_df <- combined_df[ , colnames(combined_df) != "Location.y"]
colnames(combined_df)[colnames(combined_df) == "Location.x"] <- "Location"

###Save (optionally)---------------------------------------------------------------------------------
write.csv(combined_df,  "combined_precip_2024.csv",        row.names = FALSE)
write.csv(locations_df, "monthly_precip_stations.csv",     row.names = FALSE)
write.csv(precip_df,    "monthly_precip_2024.csv",         row.names = FALSE)