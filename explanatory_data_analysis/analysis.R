# Load necessary libraries
library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(Matrix)
library(purrr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(igraph)
library(mondate)
library(zoo)
library(ggplot2)
library(dplyr)
library(gridExtra)


#### DATA ####



# Read data from csv and excel file
results <- read.csv("mens_results.csv")
nominated <- read_excel("Qualified_athletes.xlsx",sheet = "Men")


# Add 'sex' column to each dataset
results$Sex <- "M"


head(results)
head(nominated)

# Clean nominated files

nominated<- nominated[!is.na(nominated$QP), ]




# Convert Event.Date to Date format explicitly
results <- results %>%
  mutate(
    Event.Date = dmy(Event.Date),  # Use dmy() to ensure the date format is day/month/year
    Year = year(Event.Date),       # Extract the year from the date
    Quarter = as.yearqtr(Event.Date),  # Extract quarter with year
    Competition_ID = paste(Event.Name, Year, sep = "_")
  )


# Calculate the quarter index starting from Q1 2018
start_quarter <- as.yearqtr("2018 Q1")
print(start_quarter)
results <- results %>%
  mutate(
    time_index = as.numeric((Quarter - start_quarter) * 4) + 1
  )


# Generate a numeric ID for each unique Competition_ID
results <- results %>%
  arrange(Competition_ID) %>%
  mutate(
    Event_ID = as.numeric(factor(Competition_ID, levels = unique(Competition_ID)))
  )



# Filter results based on nominated athletes
data <- results %>%
  filter(NAME %in% nominated$Athlete)

# Filter out rows where the Numeric_ID appears only once
data <- data %>%
  group_by(Sex, Event.Type, Event_ID) %>%
  filter(n() > 1) %>%
  ungroup()


# Sort Place within each Numeric_ID and reassign places starting from 1
# Handle NA values appropriately
data <- data %>%
  group_by(Event.Type, Event_ID) %>%
  arrange(Event_ID, desc(PLACE), na.last = TRUE) %>%
  mutate(
    New_Place = if_else(is.na(PLACE), max(row_number()), row_number())
  ) %>%
  ungroup()


# Count number of races each athlete (Name) participated in
athlete_race_counts <- data %>%
  group_by(NAME) %>%
  summarize(Count_of_races = n_distinct(Event_ID)) %>%
  ungroup()

print(mean(athlete_race_counts$Count_of_races))

number_of_nominated <- n_distinct(nominated$Athlete)

print(number_of_nominated)
uniq<-unique(data$NAME)

# Find the missing athletes
missing_athletes <- setdiff(nominated$Athlete, uniq)
print(missing_athletes)
# 2 athletes nominated in the Olympics but not in the dataset
# Yaseen ABDALLA never ran marathon
# Matthias KYBURZ




#### GRAPHS ####

library(dplyr)
library(igraph)

# Assuming 'data' contains the full dataset with the 'Event.Type' column

# Split the data into four separate datasets

men_half_marathon <- data %>%
  filter(Sex == "M" & grepl("Half Marathon", Event.Type))

men_marathon <- data %>%
  filter(Sex == "M" & grepl("Marathon", Event.Type) & !grepl("Half Marathon", Event.Type))

# Function to create edges from the dataset
create_edges <- function(df) {
  df %>%
    group_by(Event_ID) %>%
    do({
      pairs <- t(combn(unique(.$NAME), 2, simplify = TRUE))
      data.frame(from = pairs[,1], to = pairs[,2])
    }) %>%
    ungroup() %>%
    distinct(from, to)
}

# Create edge lists for each event type
edges_men_half_marathon <- create_edges(men_half_marathon)
edges_men_marathon <- create_edges(men_marathon)

# Store edge lists in a named list
edge_lists <- list(
  "Men Half-Marathon" = edges_men_half_marathon,
  "Men Marathon" = edges_men_marathon
)

# Function to create a graph, plot it, and return the igraph object
plot_graph <- function(edge_list, event_type) {
  g <- graph_from_data_frame(edge_list, directed = FALSE)
  
  # Optionally color nodes by gender
  gender_colors <- ifelse(grepl("Men", event_type), "orange", "pink")
  
  plot(g, vertex.label = V(g)$name, vertex.size = 5, vertex.label.cex = 0.5,
       vertex.color = gender_colors,
       edge.arrow.size = 0.5, edge.curved = TRUE, 
       main = paste("Event type:", event_type))
  
  return(g)  # Return the graph object
}

# Set up a 1x2 plotting layout
par(mfrow = c(1, 2))  # 1x2 layout


# 1. Men's Half-Marathon
plot_graph(edge_lists[["Men Half-Marathon"]], "Men Half-Marathon")

# 2. Men's Marathon
plot_graph(edge_lists[["Men Marathon"]], "Men Marathon")

# Reset layout to default
par(mfrow = c(1, 1))




centrality_summary <- function(graph) {
  # Calculate centrality measures
  degree_centrality <- degree(graph)
  betweenness_centrality <- betweenness(graph, normalized = TRUE)
  closeness_centrality <- closeness(graph, normalized = TRUE)
  
  # Remove NA values and ensure all measures are numeric
  degree_centrality <- degree_centrality[!is.na(degree_centrality)]
  betweenness_centrality <- betweenness_centrality[!is.na(betweenness_centrality)]
  closeness_centrality <- closeness_centrality[!is.na(closeness_centrality)]
  
  if (length(degree_centrality) > 0 && length(betweenness_centrality) > 0 && length(closeness_centrality) > 0) {
    
    # Calculate summary statistics
    summary_stats <- data.frame(
      Metric = c("Degree", "Betweenness", "Closeness"),
      Mean = c(mean(degree_centrality), 
               mean(betweenness_centrality), 
               mean(closeness_centrality)),
      Median = c(median(degree_centrality), 
                 median(betweenness_centrality), 
                 median(closeness_centrality)),
      SD = c(sd(degree_centrality), 
             sd(betweenness_centrality), 
             sd(closeness_centrality))
    )
    
    return(summary_stats)
  } else {
    warning("Centrality measures are not numeric or contain NA values.")
    return(NULL)
  }
}

# Calculate summary statistics for each event type


g_men_half<-plot_graph(edge_lists[["Men Half-Marathon"]], "Men Half-Marathon")

# 2. Men's Marathon
g_men_marathon<-plot_graph(edge_lists[["Men Marathon"]], "Men Marathon")


summary_men_half <- centrality_summary(g_men_half)
summary_men_marathon <- centrality_summary(g_men_marathon)

# Combine the summaries into a single table
combined_summary <- do.call(rbind, list(
  cbind(Event = "Men Half-Marathon", summary_men_half),
  cbind(Event = "Men Marathon", summary_men_marathon)
))

# Display the summary table
print(combined_summary)

# Convert to LaTeX table format
library(xtable)
latex_table <- xtable(combined_summary, caption = "Centrality Measure Summary for Each Event", label = "tab:centrality-summary")
print(latex_table, include.rownames = FALSE)



#### PLOTS ####

# Count number of unique participants (Name) in each race (Numeric_ID)
race_participants <- data %>%
  group_by(Event_ID, Sex) %>%
  summarize(Number_of_people = n_distinct(NAME)) %>%
  ungroup()

# Split the data by gender
race_participants_men <- race_participants %>% filter(Sex == "M")


# Plot histogram of number of people in each race - Men
plot_men_races <- ggplot(race_participants_men, aes(x = Number_of_people)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(x = "Number of People in Race", y = "Count") +
  ggtitle("Distribution of Men's Races by Number of Participants") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



# Count number of races each athlete (Name) participated in
athlete_race_counts <- data %>%
  group_by(NAME, Sex) %>%
  summarize(Count_of_races = n_distinct(Event_ID)) %>%
  ungroup()

# Split the data by gender
athlete_race_counts_men <- athlete_race_counts %>% filter(Sex == "M")


# Plot histogram of count of races each athlete participated in - Men
plot_men_athletes <- ggplot(athlete_race_counts_men, aes(x = reorder(NAME, -Count_of_races), y = Count_of_races)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(x = "Athlete", y = "Count of Races") +
  ggtitle("Number of Races Each Male Athlete Participated In") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5),  # Adjust text size for better readability
        plot.title = element_text(hjust = 0.5)) +
  coord_flip()



# Arrange the four plots into a 1x2 grid
grid.arrange(plot_men_races,  plot_men_athletes, ncol = 2)

#### SAVING RESULTS ####

data_men <- subset(data, Sex == "M")

data_men$competitor_id <- as.numeric(factor(data_men$NAME))

# Create a unique identifier for each athlete
sorted_ids <- data_men %>% distinct(NAME) %>% pull(NAME) %>% sort()
n <- length(sorted_ids)



extract_and_rename <- function(data) {
  # Select the desired columns
  selected_data <- data[c("time_index","NAME","NAT.","BIRTH.DATE" ,"Event.Type", "Event.Date","Event.Location", "Event_ID", "competitor_id", "New_Place")]
  
  # Rename the columns
  colnames(selected_data) <- c("time_index","name", "nationality" , "birth_date","event_type", "event_date","event_location", "event_id", "competitor_id", "place")
  
  return(selected_data)
}

# Assuming your dataset is loaded into a variable named 'mydata'
# Apply the function to the dataset
new_data_men <- extract_and_rename(data_men)


# Saving the dataset as a CSV file
write.csv(new_data_men, "data_men.csv", row.names = FALSE)



