

Calculate percentage:
```
  # Calculate various statistics grouped by zipcode
  statistics <- data %>%
    group_by(zip) %>%
    summarise(
      total_count = n(),
      female_percentage = mean(sex == 2) * 100,
      dual_percentage = mean(dual == 1) * 100,
      mean_age = mean(age),
      death_rate = sum(death) / total_count * 100,
      percentage_race_labelAsian = mean(race_label == "Asian") * 100,
      percentage_race_labelBlack = mean(race_label == "Black") * 100,
      percentage_race_labelOther = mean(race_label == "Other") * 100,
      percentage_race_labelWhite = mean(race_label == "White") * 100,
      percentage_race_labelHispanic = mean(race_label == "Hispanic") * 100,
      percentage_race_labelNorth_American_Native = mean(race_label == "North American Native") * 100
    )
    
```