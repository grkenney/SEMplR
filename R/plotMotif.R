sem_file_list <- list.files("/users/g/k/gkenney/semplr/testing_sempl/SEMpl/SEMs",
                            pattern = ".sem", full.names = T)
dt <- semMotifBinding(chr = "chr12",
                      pos = 94136009,
                      ref = "G",
                      alt = "C",
                      semFiles = sem_file_list,
                      semBaselines = "/users/g/k/gkenney/semplr/testing_sempl/SEMpl/SEMs/BASELINE/SEMs_baseline_norm.txt")

psm <- plotSemMotifs(dt)


library(tidyr)

sems <- lapply(sem_file_list, fread, header = TRUE)

names(sems) <-
  lapply(sems, \(x) colnames(x)[[1]]) |>
  unlist() |>
  paste0("_", 1:length(sems))

long_sems <- sems$JUN_162 %>%
  pivot_longer(cols=c(A, C, G, T), names_to="AA", values_to="score")

ggplot(long_sems, aes(x=JUN, y=score, color=AA)) +
  geom_hline(yintercept=0, color = "grey", size = 2) +
  geom_hline(yintercept=-0.5, linetype="dashed", color = "grey", size = 2) + # this isn't exactly the right spot
  geom_text(aes(label=AA), size=8) +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Position in JUN Motif") +
  ylab("SNP Effect Score") +
  theme(text = element_text(size = 20))
