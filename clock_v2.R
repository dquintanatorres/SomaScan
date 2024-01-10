library(dplyr)
library(glmnet)
library(ggplot2)
# library(glmmLasso)
library(ggpubr)
library(coefplot)

set.seed(1000)

#### Set colors ----
cols <- list("WT" = "#6DA9B0", "KO" = "#39646A")
 
## Load data ----

leh <- readRDS("mice_lehallier.rds")
uni <- readRDS("mice_uniovi.rds")


# ## Only common genes
common <- intersect(names(leh), names(uni))
leh <- dplyr::select(leh, all_of(common)) %>% 
  mutate(cohort = "Lehallier")
uni <- dplyr::select(uni, all_of(common))


# ## Merge
uni$cohort <- ifelse(uni$cohort %in% c("A", "B"), "Laki", "Zmpste")
merged <- rbind(uni, leh) %>%
  mutate(
    id = NULL,
    dataset = NULL,
    model = NULL
  ) %>%
  group_by(cohort) %>%
  mutate(across(-c(!where(is.numeric), "age"), scale)) %>%
  ungroup()



## Split 
WTpos <- c(1:nrow(merged))[merged$genotype == "WT"]
training_set <- sample(WTpos, size = 0.7 * length(WTpos))
training <- merged[training_set,] %>% dplyr::select(-c(genotype))
testing <- merged[-training_set,] %>% dplyr::select(-c(genotype))

## LASSO
x <- model.matrix( ~ ., data = training[,-1])[, -1]
y <- training$age

cv.out <- cv.glmnet(x, y, alpha = 1, folds = 30, nlambda = 1000)
plot(cv.out)

# best_lambda <- cv.out$lambda.min
best_lambda <- cv.out$lambda.1se

lasso.mod2 <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coeff <- coefficients(lasso.mod2)
clock_genes <- coeff %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  dplyr::filter(s0 != 0,
                rownames(.) != "(Intercept)") %>% 
  rownames() %>% 
  append("age")

newx <- model.matrix( ~ ., data = testing[,-1])[, -1]

clock.pred <- predict(lasso.mod2, newx = newx , type = "response") %>% .[ ,1]

testing2 <- merged[-training_set,] %>% 
  mutate(age_pred = clock.pred,
         cohort = factor(cohort, levels = c("Laki", "Zmpste", "Lehallier")))
  

## f5a ----

f5a <- testing2 %>% 
  dplyr::rename(Group = "genotype") %>% 
  mutate(Group = ifelse(Group == "KO", "Progeroid", "Control"),
         Group = factor(Group, levels = c("Control", "Progeroid")),
         cohort = case_when(cohort == "Laki" ~ "Lmna",
                            cohort == "Zmpste" ~ "Zmpste24", 
                            TRUE ~ "Lehallier"),
         cohort = factor(cohort, levels = c("Lmna", "Zmpste24", "Lehallier"))) %>% 
  ggplot(aes(x = Group, y = age_pred-age, color = Group)) +
  scale_color_manual(values = c(cols$WT, cols$KO)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(0.8), alpha = .4) + 
  ylab("Actual - predicted age") +
  xlab(NULL) +
  stat_compare_means(size = 3, label = "p.signif", label.x.npc = 0.4, show.legend = F) +
  theme_classic() +
  theme(text = element_text(size = 12.5),
        legend.position = "bottom",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) + 
  facet_wrap(~cohort, scales = "free_x") 


  
## f5b ----

f5b <- testing2 %>% 
  dplyr::rename(Group = "genotype") %>% 
  mutate(Group = ifelse(Group == "KO", "Progeroid", "Control"),
         Group = factor(Group, levels = c("Control", "Progeroid")),
         cohort = case_when(cohort == "Laki" ~ "Lmna",
                            cohort == "Zmpste" ~ "Zmpste24", 
                            TRUE ~ "Lehallier"),
         cohort = factor(cohort, levels = c("Lmna", "Zmpste24", "Lehallier"))) %>% 
  ggplot(aes(x = Group, y = age_pred, color = Group)) +
  scale_color_manual(values = c(cols$WT, cols$KO)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(0.8), alpha = .4) + 
  stat_compare_means(size = 3, label = "p.signif", label.x.npc = 0.4, show.legend = F) +
  ylab("Predicted age") +
  xlab(NULL) +
  theme_classic() +
  theme(text = element_text(size = 12.5),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  facet_wrap(~cohort, scales = "free_x") 
  

##f5e ----
temp <- as.matrix(lasso.mod2$beta) 
markers <- rownames(temp)[temp != 0 & rownames(temp) != "cohortLehallier" ]
f5e <- coefplot(lasso.mod2, intercept = F, pointSize = 4, color = "#ECC8AE", coefficients = markers) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(-2, 3), labels = c(-2, "", 0, "", 2, "")) +
  theme(text = element_text(size = 12.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),) +
  ggtitle(NULL)


corrplot::corrplot(cor(merged[clock_genes[-c(1, length(clock_genes))]]), type = "lower", order = "hclust", method = "pie")



## Bootstrap lasso ----

bootstrap_glmnet <- function(df, pos, seed, exclude, split_min, alpha = 1, min=3, predict_type = "response"){
  set.seed(seed)
  
  training_set <- sample(pos, size = 0.7 * length(WTpos))
  
  training_fails <- min(df[training_set,] %>% dplyr::select(any_of(split_min)) %>% group_by_at(split_min) %>% summarise(n = n()) %>% pull(n))  
  if (training_fails < min ) return(NA) # at least min samples in every group
  
  testing_fails <- min(df[-training_set,] %>% dplyr::select(any_of(split_min)) %>% group_by_at(split_min) %>% summarise(n = n()) %>% pull(n))  
  
  if (testing_fails < min) return(NA) # at least min samples in every group
  
  training <- df[training_set,] %>% dplyr::select(-any_of(exclude))
  testing <- df[-training_set,] %>% dplyr::select(-any_of(exclude))
  
  ## LASSO
  x <- model.matrix( ~ ., data = training[,-1])[, -1]
  y <- training$age
  
  cv.out <- cv.glmnet(x, y, alpha = alpha, folds = 20, nlambda = 500)
  
  best_lambda <- cv.out$lambda.min
  
  lasso.mod2 <- glmnet(x, y, alpha = alpha, lambda = best_lambda)
  
  newx <- model.matrix( ~ ., data = testing[,-1])[, -1]
  
  preds <- predict(lasso.mod2, newx = newx , type = predict_type) %>% .[ ,1]
  
  return(list("training" = training_set,
              "model" = lasso.mod2,
              "lambda" = best_lambda,
              "preds" = preds))
}


exclude = c("genotype")
split_min = c("cohort", "genotype")
WTpos <- c(1:nrow(merged))[merged$genotype == "WT"]

l <- list()
for (i in 1:1500) {
  print(paste0("Iteration: ", i, "/1500"))
  l[[i]] <- bootstrap_glmnet(merged, WTpos, i, exclude = exclude, split_min = split_min , alpha = 1)
}

l2 <- l[!is.na(l)]
l3 <- l2[1:200] # We don't need that many points. It was a bit too much. Let's keep "only" 200

summ_l3 <- lapply(l3, function(x){
  print(paste0("Iteration: ", i, "/200"))
  tmp <- merged[-x$training,] %>% 
    mutate(age_pred = x$preds,
           cohort = factor(cohort, levels = c("Laki", "Zmpste", "Lehallier")),
           residual = age_pred - age
           ) %>% 
    group_by(cohort, genotype) %>% 
    summarise(value = mean(residual),
              value2 = mean(age_pred))
  }) %>% data.table::rbindlist()
  


#### Get coeffs ----

slected <- lapply(l3, function(x){
coeff <- coefficients(x$model)
coeff %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  dplyr::filter(s0 != 0,
                rownames(.) != "(Intercept)")
}) 

all_coeff <- lapply(slected, rownames)
coeffs <- names(table(unlist(all_coeff)))[table(unlist(all_coeff))>60]

slected_filtered <- lapply(slected, function(x) rownames(x)[rownames(x)%in%coeffs])




##f5c ----

f5c <- summ_l3 %>%
  dplyr::rename(Group = "genotype") %>% 
  mutate(Group = ifelse(Group == "KO", "Progeroid", "Control"),
         Group = factor(Group, levels = c("Control", "Progeroid")),
         cohort = case_when(cohort == "Laki" ~ "Lmna",
                            cohort == "Zmpste" ~ "Zmpste24", 
                            TRUE ~ "Lehallier"),
         cohort = factor(cohort, levels = c("Lmna", "Zmpste24", "Lehallier"))) %>% 
  ggplot(., aes(Group, y = value, color = Group)) +
  scale_color_manual(values = c(cols$WT, cols$KO)) +
  geom_violin(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(0.8), alpha = .4) + 
  stat_compare_means(size = 3, label = "p.signif", label.x.npc = 0.3, show.legend = F) +
  ylab("Mean actual - predicted age") +
  facet_wrap(~cohort, scales = "free_x") + 
  theme_classic() + 
  xlab(NULL) +
  theme(text = element_text(size = 12.5),
        legend.position = "bottom",
        strip.text = element_text(face = "italic"),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())

##f5d ----

f5d <- summ_l3 %>%
  dplyr::rename(Group = "genotype") %>% 
  filter(cohort != "Lehallier") %>% 
  mutate(Group = ifelse(Group == "KO", "Progeroid", "Control"),
         Group = factor(Group, levels = c("Control", "Progeroid")),
         cohort = case_when(cohort == "Laki" ~ "Lmna",
                            cohort == "Zmpste" ~ "Zmpste24", 
                            TRUE ~ "Lehallier"),
         cohort = factor(cohort, levels = c("Lmna", "Zmpste24", "Lehallier"))) %>% 
  ggplot(., aes(Group, y = value2, color = Group)) +
  scale_color_manual(values = c(cols$WT, cols$KO)) +
  geom_violin(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(0.8), alpha = .4) + 
  stat_compare_means(size = 3, label = "p.signif", label.x.npc = 0.3, show.legend = F) +
  ylab("Mean predicted age") +
  xlab(NULL) +
  facet_wrap(~cohort, scales = "free_x") + 
  theme_classic() + 
  theme(text = element_text(size = 12.5),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.ticks.x = element_blank())


##f5f ----

f5f <- table(unlist(slected_filtered)) %>% 
  as.data.frame() %>% as_tibble() %>% 
  filter(Var1 != "cohortLehallier") %>% 
  arrange(Freq, Var1) %>% 
  mutate(Freq = Freq*100/200,
         Var1 = factor(Var1, levels = Var1)) %>% 
  ggplot(., aes(y = Var1, x = Freq)) +
  geom_segment(aes(x = 0, xend = Freq, yend = Var1), color = "#6A8D73", linetype = 2) +
  geom_point(color = "#ECC8AE", size = 3) + 
  scale_x_continuous(expand = c(0, 4), limits = c(0,100)) +
  theme_bw() + theme(text = element_text(size = 12.5),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor = element_blank(),) +
  ylab("Selected feature") +
  xlab("Frequency (%)")




## Final fig ----
f5boxes <- ggarrange(f5a, f5b, f5c, f5d, labels = c("(a)", "(b)", "(c)", "(d)"), common.legend = T, legend = "bottom")


f5 <- ggarrange(f5boxes, f5e, f5f, nrow = 1, labels = c("", "(e)", "(f)"), widths = c(0.55, 0.2, 0.25))

f5

ggsave("f5.pdf", plot = f5, path = "figures/", device = "pdf", height = 200, width = 180, units = "mm")


pdf("figures/supp.corplot.pdf")
corrplot::corrplot(cor(merged[clock_genes[-c(1, length(clock_genes))]]), type = "lower", order = "hclust", method = "pie")
dev.off()
