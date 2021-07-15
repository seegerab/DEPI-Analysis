### Create a mock data frame that models Ally's problem
genotype_col <- c("A", "B", "C", "D", "E", "F")
terminal_go_col <- c(NA, NA, NA, "TERMINAL", "TERMINAL", NA)

data_frame <- as.data.frame(cbind(genotype_col, terminal_go_col))

terminal_go <- list()
j <- 1
for (i in 1:nrow(data_frame)){
  if (!is.na(data_frame[i,2])){
    terminal_go[[j]] <- data_frame[i,1]
    j <- j + 1
  }
}

unlist(terminal_go)

### Alternately:
library(dplyr)
(data_frame%>%
  filter(!is.na(terminal_go_col)))$genotype_col
