## 1. MultiQC reports of raw reads

<img width="689" alt="Screen Shot 2023-08-14 at 11 41 43 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/0c98dffa-db4d-4650-ad4c-530d1440f70a">

<img width="700" alt="Screen Shot 2023-08-07 at 11 53 26 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/279d8a32-5159-46f6-bd19-788a8bb407f8">

<img width="682" alt="Screen Shot 2023-08-07 at 11 54 41 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/266abbab-9a57-4b78-903d-60564aaff7d4">

<img width="680" alt="Screen Shot 2023-08-14 at 11 42 15 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7a13d8ad-4ff5-4370-b919-62afa4489372">

<img width="693" alt="Screen Shot 2023-08-14 at 11 43 07 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/4c6544d9-6da9-4a0d-bd48-0e55a1f2b7af">

<img width="695" alt="Screen Shot 2023-08-14 at 11 43 19 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/62cd9b7e-86ee-41aa-b41a-a68ca6a7f75c">

<img width="688" alt="Screen Shot 2023-08-14 at 11 43 35 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/7a139834-e9cd-4902-bd2b-bccfd87480cb">

<img width="692" alt="Screen Shot 2023-08-14 at 11 43 54 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/ca3bc6e4-7311-4e60-b792-a6bad74cc540">


## 2. MultiQC reports of trimmed reads

<img width="699" alt="Screen Shot 2023-08-14 at 11 44 42 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/26cb0769-c9a6-4b7e-b74a-4ba7b3caff5d">

<img width="696" alt="Screen Shot 2023-08-14 at 11 44 58 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/5739b726-2d0b-4eaa-aa1e-040ade9ddb9f">

<img width="690" alt="Screen Shot 2023-08-14 at 11 45 14 AM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/8f2658e5-2e2d-4066-996e-1aa1bcbd3d34">

## 3. Alignment rates

<img width="625" alt="Screen Shot 2023-08-14 at 1 32 07 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/95943c60-6a32-47c2-8ee9-3ed0f464629f">

Code made in R:
```{r}
sequencing_data <- read.csv("../../RNA_extraction_sequencing_data.csv")

sequencing_data %>% select(Species:Treatment,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  ggplot(., aes(x=Sample.ID, y=percent_alignment, fill=Treatment)) + geom_col() + facet_wrap(~Species) + theme_classic()
```

## 4. Number of trimmed reads aligned 

<img width="624" alt="Screen Shot 2023-08-14 at 1 30 10 PM" src="https://github.com/ademerlis/temperaturevariability2023/assets/56000927/9cc119b1-5472-48ab-bb64-98b62b742a1d">

Code made in R:
```{r}
sequencing_data <- read.csv("../../RNA_extraction_sequencing_data.csv")

sequencing_data %>% 
select(Species:Treatment,Raw.Reads..Million.,M.Reads.After.Filtering..Trimmed.reads.,Mapping.Rate..Post.trimming.alignment.) %>% 
  mutate(percent_alignment = gsub("%","", Mapping.Rate..Post.trimming.alignment.)) %>% 
  mutate(percent_alignment = as.numeric(percent_alignment)) %>% 
  mutate(percent_alignment = percent_alignment / 100) %>% 
  mutate(million_reads_aligned = (M.Reads.After.Filtering..Trimmed.reads.)*percent_alignment) %>% 
  ggplot(., aes(x=Sample.ID, y=million_reads_aligned, fill=Treatment)) + geom_col() + facet_wrap(~Species) + theme_classic()
```
