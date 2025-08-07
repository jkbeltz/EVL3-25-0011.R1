
plot_pca_tpt = function(df, samps){
  ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint)) +
    geom_point(size = 5) +
    xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 2), " % expl. var.)")) +
    ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 2), " % expl. var.)")) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw(base_size = 15) 
}
load('./orch2020.RData')

#PCA on all samples 
load('./orch2020.RData')
df = afmat
samps = samps %>% rename(Timepoint = tpt)

pca_res <- prcomp(na.omit(t(df)), scale = TRUE, center = TRUE)
pca_data <- as.data.frame(pca_res$x)
pca_data <- cbind(samps, pca_data)
pca_data$Timepoint = as.numeric(pca_data$Timepoint)


ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment)) +
  geom_point(size = 5) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 2), " % expl. var.)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 2), " % expl. var.)")) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 15) 

ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment, label = full.sample.name)) +
  geom_point(size = 5) +
  geom_text(hjust=0, vjust=0) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 2), " % expl. var.)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 2), " % expl. var.)")) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 15)


##removing bloomington - just apple, AT, and LB
load('./orch2020.RData')
df = cbind(samps, t(afmat))
dim(df)
df = df %>% filter(treatment != 'B')
dim(df)
df = df[,-c(1:ncol(samps))]
df = t(df)
samps = samps %>% rename(Timepoint = tpt)
931884042
631884042
samps = samps %>% filter(treatment != 'B')
pca_res <- prcomp(na.omit(t(df)), scale = TRUE, center = TRUE)
pca_data <- as.data.frame(pca_res$x)
pca_data <- cbind(samps, pca_data)
pca_data$Timepoint = as.numeric(pca_data$Timepoint)


ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint, shape = treatment)) +
  geom_point(size = 5) +
  xlab(paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 2), " % expl. var.)")) +
  ylab(paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 2), " % expl. var.)")) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw(base_size = 15)

#Stats
treats = c('A', 'AT', 'LB')
stats = data.frame()
for (treat in treats){
  d = pca_data %>% filter(treatment == treat)
  pval = summary(aov(PC1 ~ Timepoint, data = d))[[1]][1, 5]
  d.s = cbind(treat, pval)
  stats = rbind(stats, d.s)
  
}
stats
A data.frame: 3 × 2
treat	pval
<fct>	<fct>
  A	9.77085020137192e-05
AT	0.00121526821315374
LB	0.000426077951594504