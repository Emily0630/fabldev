library(RecordLinkage)
library(glue)
library(tictoc)
library(fabldev)

data <- RecordLinkage::RLdata10000 %>%
  mutate(unique_id = RecordLinkage::identity.RLdata10000)

duplicates <- nrow(data) * .1

duplicated_ids <- data %>%
  filter(duplicated(unique_id)) %>%
  select(unique_id) %>%
  pull()

duplicated_records <- data %>%
  filter(unique_id %in% duplicated_ids) %>%
  arrange(unique_id) %>%
  mutate(rn = row_number())

duplicated_1 <- duplicated_records %>%
  filter(rn %% 2 == 0)

duplicated_2 <- duplicated_records %>%
  filter(rn %% 2 == 1)

non_duplicated_records <- data %>%
  filter(!(unique_id %in% duplicated_ids)) %>%
  mutate(rn = row_number())

non_duplicated_1 <- non_duplicated_records %>%
  filter(rn %% 2 == 0)

non_duplicated_2 <- non_duplicated_records %>%
  filter(rn %% 2 == 1)

df1 <- rbind(duplicated_1, non_duplicated_1)
df2 <- rbind(duplicated_2, non_duplicated_2)
n1 <- nrow(df1)
n2 <- nrow(df2)
Z_true <- rep(0, nrow(df1))
Z_true[1:duplicates] <- 1:duplicates

fields <- c(1, 3, 5, 6, 7)
types <- c("lv", "lv", "bi", "bi", "bi")

cd <- compare_records(df1, df2, flds = fields, types = types,
                      breaks = c(0, .15))

hash <- hash_comparisons(cd, all_patterns = F)

kappa_vec <- c(.5, .6, .7, .8, .9, 1)
B <- 1000
holdout_size <- 1000
tau <-  1

svi_efficient(hash, k = 1, B = 1000)

chain <- brl_efficient_serge(hash, S = 50, mode = "rejection")
