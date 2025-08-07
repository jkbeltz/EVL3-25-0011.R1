
calc_Neff = function (rd, poolCt) 
{
  ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}
##GLM of apple and bloom t1->5 (run on cluster): 


library('dplyr')
library('tidyr')
library('data.table')
library('doMC')



###FUNCTIONS
fit_GLM_ContinuousTime = function (afMatrix, rdMatrix, sampleData, vec, model.vars, poolCt = 100, 
                                   ncores) 
{
  registerDoMC(ncores)
  df = as.data.frame(sampleData[, colnames(sampleData) %in% 
                                  model.vars])
  colnames(df) <- model.vars
  formulaString = paste0(colnames(df), collapse = " + ")
  formulaString = paste0("cts ~ ", formulaString)
  cat("Model Formula is: \n", formulaString, "\n")
  Neff = calc_Neff(rdMatrix, poolCt)
  do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
    if (ix%%10000 == 0) {
      cat("working on site ", ix, "\n")
    }
    cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
    ] * (1 - afMatrix[ix, ]))))
    df$cts = cts
    df = as.data.frame(df)
    model = glm(formulaString, family = "quasibinomial", 
                data = df)
    cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
    results = c(cp[, 1], cp[, 2])
    names(results) = c(paste0("coef.", vec[1], "_", vec[length(vec)]), 
                       paste0("p.", vec[1], "_", vec[length(vec)]))
    return(results)
  }, mc.cores = ncores))
}

#calucluate N effective
calc_Neff = function (rd, poolCt) 
{
  ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}


###Loop through treatments and run glm - save results as RData file
treats = c('A', 'B')
load('./orch2020_filtered.RData')
df.glm.meta = sites
vec = c(1:5)
for (treat in treats){
  load('./orch2020_filtered.RData')
  samps = samps %>% mutate(tpt = as.numeric(tpt))
  df.glm = sites
  afmat = t(afmat)
  eec = t(eec)
  afmat = cbind(samps, afmat)
  eec = cbind(samps, eec)
  afmat = afmat %>%
    filter(tpt %in% vec & treatment == treat)
  eec = eec %>%
    filter(tpt %in% vec & treatment == treat)
  samps = afmat[,1:ncol(samps)]
  afmat = afmat[,-c(1:ncol(samps))]
  afmat = t(afmat)
  afmat = as.data.frame(afmat)
  eec = eec[,-c(1:ncol(samps))]
  eec = t(eec)
  eec = as.data.frame(eec)
  res = fit_GLM_ContinuousTime(afMatrix = afmat ,rdMatrix = eec, vec = vec, sampleData = samps, model.vars = 'tpt', poolCt=100, ncores = 16)
  df.glm = cbind(df.glm, res) %>% dplyr::select(-chrom, -pos)
  names(df.glm) = c(paste0('coef.1_5'),paste0('p.1_5'))
  save(df.glm, file = paste0("../05_GLM/glm.", treat, ".TP1_5.RData"))
}
##Getting significant sites of each treatment through time:
get.sig.sites = function(glm.file, rdata, comps, fdrThreshs, esThreshs){
  chroms = c('2L', '2R', '3L', '3R', 'X')
  load(rdata)
  af.shifts = get_af_shifts(afmat, samps, cage_set = NULL, comparisons = comps)
  #Get fdr 
  load(glm.file)
  FDR = get_glm_FDR.V2(df.glm = df.glm)
  load(glm.file)
  df.glm = cbind(sites, df.glm)
  df.sig = get_sig_sites(df.glm, comparisons = comps, FDR  , af.shifts, fdrThreshs,
                         esThreshs)
  return(df.sig)
}


df.sig.A = get.sig.sites(glm.file = './glm.A.TP1_5.RData', rdata = '../RData/orch2020_filtered.ACages.RData',
                         comps = c('1_5'), fdrThreshs= c(.2, .05, 0.01), esThreshs= c(0.005, 0.02, 0.02))
dim(df.sig.A)
568479
write.csv(df.sig.A, './df.sig.A.1_5.csv', row.names = FALSE)
dim(read.csv('./df.sig.A.1_5.csv'))
dim(read.csv('./df.sig.B.1_5.csv'))
148479
813269
df.sig.A %>% group_by(chrom) %>% count()
A grouped_df: 5 × 2
chrom	n
<chr>	<int>
  2L	5033
2R	2497
3L	2438
3R	3994
X	885
df.sig.B = get.sig.sites(glm.file = './glm.B.TP1_5.RData', rdata = '../RData/orch2020_filtered.BCages.RData',
                         comps = c('1_5'), fdrThreshs= c(.2, .05, 0.01), esThreshs= c(0.005, 0.02, 0.02))
write.csv(df.sig.B, './df.sig.B.1_5.csv', row.names = FALSE)
df.sig.B %>% group_by(chrom) %>% count()
A grouped_df: 5 × 2
chrom	n
<chr>	<int>
  2L	33737
2R	13048
3L	14884
3R	13034
X	6623
write.csv(df.sig.B, './df.sig.B.1_5.csv', row.names = FALSE)
###Plotting significant sites:
df.sig = read.csv('./df.sig.A.1_5.csv')
df.sig = df.sig %>% filter(FDR < 0.1)
df.sig <- df.sig %>% mutate(sigLabel=as.factor(sigLevel),compLabel = as.factor(comparison))

df.sig = df.sig %>% filter(abs(afShift) > 0.005)
p.A = df.sig %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,3.5))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=20),strip.text.y = element_text(size=15),
        axis.text.y=element_text(size=15), axis.text.x=element_text(size=15),
        axis.title = element_text(size = 15)) +
  labs(x="Position (Mb)",y="-log10(FDR-corrected p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=1)) +
  scale_colour_manual(values = c('grey', 'orange', 'red')) +
  theme(strip.text.y = element_blank())
p.A 
ggsave('../Figures/manhattan.A.1_5.png',p.A, height = 10, width = 10)

df.sig = read.csv('./df.sig.B.1_5.csv')
df.sig = df.sig %>% filter(FDR < 0.1)

df.sig <- df.sig %>% mutate(sigLabel=as.factor(sigLevel),compLabel = as.factor(comparison))

df.sig = df.sig %>% filter(abs(afShift) > 0.005)
p.A = df.sig %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,3.5))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=20),strip.text.y = element_text(size=15),
        axis.text.y=element_text(size=15), axis.text.x=element_text(size=15),
        axis.title = element_text(size = 15)) +
  labs(x="Position (Mb)",y="-log10(FDR-corrected p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=1)) +
  scale_colour_manual(values = c('grey', 'orange', 'red')) +
  theme(strip.text.y = element_blank())
p.A 
ggsave('../Figures/manhattan.B.1_5.png',p.A, height = 10, width = 10)

#GLM with all samples through time and treatment*time interaction
setwd('/scratch/groups/dpetrov/MarkB/Orchard2020Data/05_GLM/')

chroms=c("2L","2R","3L","3R","X")
comparisons = c("tpt", "treat", "int")
load('../RData/orch2020_filtered.RData')


library('dplyr')
library('tidyr')
library('data.table')
library('doMC')
setwd('~/dpetrov/MarkB/Orchard2020Data/')

##Functions

fit_GLM_Treatment = function (afMatrix, rdMatrix, sampleData, model.vars = "treatment", 
                              poolCt = 100, ncores) 
{
  registerDoMC(ncores)
  df = as.data.frame(sampleData[, colnames(sampleData) %in% 
                                  model.vars])
  colnames(df) <- model.vars
  formulaString = paste0(colnames(df), collapse = " + ")
  formulaString = paste0("cts ~ ", formulaString)
  cat("Model Formula is: \n", formulaString, "\n")
  Neff = calc_Neff(rdMatrix, poolCt)
  do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
    if (ix%%10000 == 0) {
      cat("working on site ", ix, "\n")
    }
    cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
    ] * (1 - afMatrix[ix, ]))))
    df$cts = cts
    df = as.data.frame(df)
    model = glm(formulaString, family = "quasibinomial", 
                data = df)
    cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
    results = c(cp[, 1], cp[, 2])
    return(results)
  }, mc.cores = ncores))
}


calc_Neff = function (rd, poolCt) 
{
  ((poolCt * 2 * rd) - 1)/(poolCt * 2 + rd)
}



treats = c('A', 'B')
load('./RData/orch2020.RData')

df.glm = sites

load('./RData/orch2020.RData')
afmat = cbind(samps, t(afmat))
samps = samps %>% mutate(treatment = as.factor(treatment))
eec = cbind(samps, t(eec))
afmat = afmat %>% filter(treatment %in% treats)
eec = eec %>% filter(treatment %in% treats)
samps = afmat[,1:ncol(samps)]
afmat = afmat[,-c(1:ncol(samps))]
afmat = as.data.frame(t(afmat))
eec = eec[,-c(1:ncol(samps))]
eec = as.data.frame(t(eec))
res = fit_GLM_Treatment(afMatrix = afmat, rdMatrix = eec, sampleData = samps, model.vars = 'treatment',
                        poolCt = 100, ncores = 16)
res = as.data.frame(res)
names(res) = c(paste0("coef.", treats[1], "_", treats[2]), paste0("pval.", treats[1], "_", treats[2]))
df.glm = cbind(df.glm, res)



save(df.glm, file = "./05_GLM/glm.AppleVBloom.RData")
#t1->5 
load('./glm.AppleVBloom.T1_5.RData')
head(df.glm)
A data.frame: 6 × 8
chrom	pos	coef.tpt	p.tpt	coef.treat	p.treat	coef.int	p.int
<chr>	<int>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>
  1	2L	5390	0.05468918	0.02319510	0.05428193	0.6156774	-0.008347045	0.8004225
2	2L	5399	-0.04029232	0.56180627	0.03181246	0.9183557	-0.029034305	0.7692138
3	2L	5465	0.08627595	0.05629994	0.06883283	0.7411048	-0.002863101	0.9630850
4	2L	5495	-0.04166446	0.30888387	0.16116501	0.3599644	-0.009430046	0.8653824
5	2L	5598	0.10688694	0.01652663	0.11726918	0.5676004	-0.019007255	0.7530812
6	2L	5762	0.04155147	0.08738923	0.03418010	0.7556382	0.005090452	0.8796194
fdr.1_5 = as.data.frame(get_glm_FDR.V2(df.glm))
write.csv(fdr.1_5, './fdr.AppleVBloom.T1_5.csv', row.names = FALSE)
fdr.1_5 = read.csv('./fdr.AppleVBloom.T1_5.csv')
fdr.1_5 = read.csv('./fdr.AppleVBloom.T1_5.csv')
load('../RData/orch2020_filtered.RData')

fdr.1_5 = cbind(sites, fdr.1_5)
head(fdr.1_5)
A data.frame: 6 × 5
chrom	pos	p.tpt	p.treat	p.int
<chr>	<int>	<dbl>	<dbl>	<dbl>
  1	2L	5390	0.3159935	0.9681455	0.9871637
2	2L	5399	0.8674358	0.9953415	0.9844266
3	2L	5465	0.4266383	0.9812863	0.9979160
4	2L	5495	0.7290291	0.9275304	0.9918014
5	2L	5598	0.2819813	0.9617595	0.9832273
6	2L	5762	0.4920755	0.9824860	0.9926615
fdr.1_5 = fdr.1_5 %>% dplyr::select(chrom, pos, p.int) %>% rename(FDR = p.int)
fdr.1_5$comparison = 'Treatment_Interaction'
fdr.1_5 <- fdr.1_5 %>% mutate(compLabel = as.factor(comparison))
fdr.1_5 = fdr.1_5 %>% mutate(sigLabel = if_else(FDR < 0.1, 'sig', 'non-sig')) %>% 
  mutate(sigLabel=as.factor(sigLabel))
fdr.1_5 %>% filter(FDR < 0.1)
A data.frame: 14 × 6
chrom	pos	FDR	comparison	compLabel	sigLabel
<chr>	<int>	<dbl>	<chr>	<fct>	<fct>
  3L	14360958	0.05246537	Treatment_Interaction	Treatment_Interaction	sig
3L	16911479	0.05246537	Treatment_Interaction	Treatment_Interaction	sig
3L	16911491	0.05246537	Treatment_Interaction	Treatment_Interaction	sig
3L	16911529	0.02518424	Treatment_Interaction	Treatment_Interaction	sig
3R	12959498	0.02759900	Treatment_Interaction	Treatment_Interaction	sig
3R	12960332	0.05246537	Treatment_Interaction	Treatment_Interaction	sig
3R	12960369	0.05246537	Treatment_Interaction	Treatment_Interaction	sig
3R	12961464	0.02518424	Treatment_Interaction	Treatment_Interaction	sig
3R	12961492	0.02518424	Treatment_Interaction	Treatment_Interaction	sig
3R	12961501	0.02518424	Treatment_Interaction	Treatment_Interaction	sig
3R	12962391	0.04858069	Treatment_Interaction	Treatment_Interaction	sig
3R	12963490	0.03376865	Treatment_Interaction	Treatment_Interaction	sig
3R	12964265	0.07288027	Treatment_Interaction	Treatment_Interaction	sig
3R	12998860	0.05246537	Treatment_Interaction	Treatment_Interaction	sig
p.int = fdr.1_5 %>% ggplot(aes(x=pos/1000000))+ theme_minimal() + lims(y=c(0,2))  +
  facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
  theme(legend.position = 'none',strip.text = element_text(size=20),strip.text.y = element_text(size=15),
        axis.text.y=element_text(size=15), axis.text.x=element_text(size=15),
        axis.title = element_text(size = 15)) +
  labs(x="Position (Mb)",y="-log10(FDR-corrected p-value)") +
  suppressWarnings(geom_point(aes(y=-log10(FDR),color=sigLabel),alpha=.8,size=1)) +
  scale_colour_manual(values = c('dark grey', 'red')) +
  theme(strip.text.y = element_blank())
p.int
ggsave('../Figures/manhattan.Interaction.1_5.png',p.int, height = 10, width = 10)
source("/home/users/mcbitter/OrchardProject/Code/Orchard2020Scripts/orch2020_functions.R")

fit_GLM_ContinuousTime.Treatment
function (afMatrix, rdMatrix, sampleData, poolCt = 100, ncores) 
{
  registerDoMC(ncores)
  df = as.data.frame(sampleData[, colnames(sampleData) %in% 
                                  c("tpt", "treatment")])
  Neff = calc_Neff(rdMatrix, poolCt)
  do.call(rbind, mclapply(1:nrow(afMatrix), function(ix) {
    if (ix%%10000 == 0) {
      cat("working on site ", ix, "\n")
    }
    cts = cbind(t(round(Neff[ix, ] * afMatrix[ix, ])), t(round(Neff[ix, 
    ] * (1 - afMatrix[ix, ]))))
    df$cts = cts
    df = as.data.frame(df)
    model = glm(cts ~ tpt + treatment + treatment * tpt, 
                family = "quasibinomial", data = df)
    cp = summary(model)$coefficients[-1, c(1, 4), drop = FALSE]
    results = c(cp[1, 1], cp[1, 2], cp[2, 1], cp[2, 2], cp[3, 
                                                           1], cp[3, 2])
    names(results) = c("coef.tpt", "p.tpt", "coef.treat", 
                       "p.treat", "coef.int", "p.int")
    return(results)
  }, mc.cores = ncores))
}
###Get trajectories of interaction-significant alleles:
int.alleles = fdr.1_5 %>% filter(FDR < 0.1)
int.alleles = int.alleles %>% mutate(snp = paste0(chrom, pos))
##apple trjaectory generation
load('../RData/orch2020_filtered.ACages.RData', verbose = TRUE)
df = cbind(sites, afmat)
df = df %>% mutate(snp = paste0(chrom, pos))
df = df %>% filter(snp %in% int.alleles$snp) #filter frequencies to sites of interest
df = df %>% dplyr::select(-snp)

shifts.int.apple = get_af_shifts(afmat = df[,-c(1:2)], samps = samps, comparison = c('1_5'))

shifts.int.apple = shifts.int.apple %>% mutate(sign.shift = sign(dAF.1_5))
Loading objects:
  sites
samps
afmat
eec
load('../RData/orch2020_filtered.ACages.RData', verbose = TRUE)
df = cbind(sites, afmat)
df = df %>% mutate(snp = paste0(chrom, pos))
df = df %>% filter(snp %in% int.alleles$snp) #filter frequencies to sites of interest
df = df %>% dplyr::select(-snp)

shifts.int.apple = get_af_shifts(afmat = df[,-c(1:2)], samps = samps, comparison = c('1_5'))

shifts.int.apple = shifts.int.apple %>% mutate(sign.shift = sign(dAF.1_5))


traj.apple = get_af_shifts(df %>% dplyr::select(-chrom, -pos), samps,cage_set = NULL, comparisons = c('1_2', '2_3','3_4','4_5'))

#Generate trajectories:
shifts.apple = matrix(nrow = nrow(df), ncol = 5)
for (row in 1:nrow(traj.apple)){
  shifts = as.vector(as.numeric(traj.apple[row, ]))
  shifts = cumsum(shifts)
  shifts = c(0, shifts)
  shifts.apple = rbind(shifts.apple, shifts)
  
}

shifts.apple = as.data.frame(na.omit(shifts.apple))
names(shifts.apple) = c('t1', 't2', 't3','t4','t5')
shifts.apple = cbind((df %>% dplyr::select(chrom, pos)), shifts.apple)
head(shifts.apple)
##alternative trjaectory generation
load('../RData/orch2020_filtered.BCages.RData', verbose = TRUE)
df = cbind(sites, afmat)
df = df %>% mutate(snp = paste0(chrom, pos))
df = df %>% filter(snp %in% int.alleles$snp) #filter frequencies to sites of interest
df = df %>% dplyr::select(-snp)

shifts.int.bloom = get_af_shifts(afmat = df[,-c(1:2)], samps = samps, comparison = c('1_5'))

shifts.int.bloom = shifts.int.bloom %>% mutate(sign.shift = sign(dAF.1_5))


traj.bloom = get_af_shifts(df %>% dplyr::select(-chrom, -pos), samps,cage_set = NULL, comparisons = c('1_2', '2_3','3_4','4_5'))

#Generate trajectories:
shifts.bloom = matrix(nrow = nrow(df), ncol = 5)
for (row in 1:nrow(traj.bloom)){
  shifts = as.vector(as.numeric(traj.bloom[row, ]))
  shifts = cumsum(shifts)
  shifts = c(0, shifts)
  shifts.bloom = rbind(shifts.bloom, shifts)
  
}

shifts.bloom = as.data.frame(na.omit(shifts.bloom))
names(shifts.bloom) = c('t1', 't2', 't3','t4','t5')
shifts.bloom = cbind((df %>% dplyr::select(chrom, pos)), shifts.bloom)
head(shifts.bloom)
Loading objects:
  sites
samps
afmat
eec
A data.frame: 6 × 7
chrom	pos	t1	t2	t3	t4	t5
<chr>	<int>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>
  1152903	3L	14360958	0	-0.013690333	0.005148333	0.02389650	0.03435983
1199457	3L	16911479	0	0.008593833	0.023486417	0.02539950	0.05255395
1199459	3L	16911491	0	0.008593833	0.023486417	0.02539950	0.05255395
1199462	3L	16911529	0	0.012551767	0.030256467	0.03073792	0.04169615
1431787	3R	12959498	0	0.012651900	0.032971767	0.03653415	0.03740252
1431816	3R	12960332	0	0.003453817	0.029801383	0.03765867	0.03251978
shifts.apple = shifts.apple %>% gather(3:7, key = tpt, value = frequency)
shifts.bloom = shifts.bloom %>% gather(3:7, key = tpt, value = frequency)
shifts.apple$treatment = 'Apple'
shifts.bloom$treatment = 'Bloomington'
shifts.apple = shifts.apple %>% mutate(snp = paste0(chrom, pos))
shifts.bloom = shifts.bloom %>% mutate(snp = paste0(chrom, pos))
shifts.meta = rbind(shifts.apple, shifts.bloom)
write.csv(shifts.meta, './shifts.InteractionSites.csv', row.names = FALSE)
ggplot(shifts.apple, aes(x = tpt, y = frequency, group = snp, colour = treatment)) +
  geom_line(size = 0.8, alpha = 0.8) +
  theme_bw(base_size = 35) +
  scale_colour_manual(values = c('grey')) +
  theme(legend.position = 'none')

ggplot(shifts.bloom, aes(x = tpt, y = frequency, group = snp, colour = treatment)) +
  geom_line(size = 0.8, alpha = 0.8) +
  theme_bw(base_size = 35) +
  scale_colour_manual(values = c('grey')) +
  theme(legend.position = 'none')

shifts.meta = shifts.meta %>% mutate(tpt = case_when(
  tpt == 't1' ~ 1,
  tpt == 't2' ~ 2,
  tpt == 't3' ~ 3,
  tpt == 't4' ~ 4,
  tpt == 't5' ~ 5,))
shifts.meta = read.csv('./shifts.InteractionSites.csv')
ggplot(shifts.meta, aes(x = tpt, y = frequency, color = treatment, group = interaction(snp, treatment))) +
  geom_line() +
  labs(
    x = "Time Point",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c('green', 'orange'))

##Testing apple sites in bloom and vice versa:
setwd('~/dpetrov/MarkB/Orchard2020Data/AppleBloomMSRevisions/')
df.a = read.csv('./df.sig.A.1_5.csv')
df.b = read.csv('./df.sig.B.1_5.csv')
#apple sites in bloom
df.a = read.csv('./df.sig.Matched.A.1_5.csv')
df.a = df.a %>% mutate(sign.shift = sign(afShift))
df.a = df.a %>% filter(sigLevel > 1)
target.a.sites = df.a %>% dplyr::select(chrom, pos)
matched.a.sites = df.a %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)
cages = c(1:6)
for (c in cages){
  load(paste0('../RData/RData_ByCage/orch2020.B',c, '.RData'))
  shifts = get_af_shifts(afmat = afmat.lo, samps = samps.lo, cage_set = NULL, comparison = c('1_5'))
  shifts = cbind(sites, shifts)
  shifts.target = left_join(target.a.sites, shifts)
  shifts.matched = left_join(matched.a.sites, shifts)
  shifts.target = shifts.target %>% mutate(dAF.1_5.target = dAF.1_5 * df.a$sign.shift )
  shifts.matched = shifts.matched %>% rename(dAF.1_5.matched = dAF.1_5) %>% rename(pos.matched = pos)
  shifts = cbind(shifts.target %>% dplyr::select(chrom, pos,dAF.1_5.target ), 
                 shifts.matched %>% dplyr::select(pos.matched, dAF.1_5.matched))
  shifts$cage = c
  shifts$id_treatment = 'Apple'
  shifts$test_treatment = 'Bloom'
  write.csv(shifts, paste0('./Apple.Sites.BloomCage',c,'Shifts.csv'), row.names = FALSE)
}


#bloom sites in apple
df.b = read.csv('./df.sig.Matched.B.1_5.csv')
df.b = df.b %>% mutate(sign.shift = sign(afShift))
df.b = df.b %>% filter(sigLevel > 1)
target.b.sites = df.b %>% dplyr::select(chrom, pos)
matched.b.sites = df.b %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)
cages = c(1:6)
for (c in cages){
  load(paste0('../RData/RData_ByCage/orch2020.A',c, '.RData'))
  shifts = get_af_shifts(afmat = afmat.lo, samps = samps.lo, cage_set = NULL, comparison = c('1_5'))
  shifts = cbind(sites, shifts)
  shifts.target = left_join(target.b.sites, shifts)
  shifts.matched = left_join(matched.b.sites, shifts)
  shifts.target = shifts.target %>% mutate(dAF.1_5.target = dAF.1_5 * df.b$sign.shift )
  shifts.matched = shifts.matched %>% rename(dAF.1_5.matched = dAF.1_5) %>% rename(pos.matched = pos)
  shifts = cbind(shifts.target %>% dplyr::select(chrom, pos,dAF.1_5.target ), 
                 shifts.matched %>% dplyr::select(pos.matched, dAF.1_5.matched))
  shifts$cage = c
  shifts$id_treatment = 'Bloom'
  shifts$test_treatment = 'Apple'
  write.csv(shifts, paste0('./Bloom.Sites.AppleCage',c,'Shifts.csv'), row.names = FALSE)
}
files = list.files(pattern = 'Shifts')
head(read.csv('Bloom.Sites.AppleCage5Shifts.csv'))
head(read.csv('Bloom.Sites.AppleCage4Shifts.csv'))
A data.frame: 6 × 8
chrom	pos	dAF.1_5.target	pos.matched	dAF.1_5.matched	cage	id_treatment	test_treatment
<fct>	<int>	<dbl>	<int>	<dbl>	<int>	<fct>	<fct>
  1	2L	292130	-0.0547930	10028226	0.0361430	5	Bloom	Apple
2	2L	377371	0.0017838	3930889	0.0331138	5	Bloom	Apple
3	2L	402435	0.0319020	13430904	0.0167100	5	Bloom	Apple
4	2L	402773	0.0654060	1863964	0.0297180	5	Bloom	Apple
5	2L	404524	0.0704920	9265285	0.0174520	5	Bloom	Apple
6	2L	495458	0.0356240	3274387	0.0192070	5	Bloom	Apple
A data.frame: 6 × 8
chrom	pos	dAF.1_5.target	pos.matched	dAF.1_5.matched	cage	id_treatment	test_treatment
<fct>	<int>	<dbl>	<int>	<dbl>	<int>	<fct>	<fct>
  1	2L	292130	0.0310090	10028226	0.0754150	4	Bloom	Apple
2	2L	377371	-0.0006615	3930889	-0.0091686	4	Bloom	Apple
3	2L	402435	0.0478490	13430904	0.0616750	4	Bloom	Apple
4	2L	402773	0.0605320	1863964	0.0180540	4	Bloom	Apple
5	2L	404524	0.0561930	9265285	0.0225170	4	Bloom	Apple
6	2L	495458	0.1429870	3274387	0.0947270	4	Bloom	Apple
#genome wide stats:
data.stats.genome = data.frame()
files = list.files(pattern = 'Shifts')
for (f in files){
  d = read.csv(f)
  id.treat = strsplit(f, '[.]')[[1]][1]
  test.treat = strsplit(strsplit(f, '[.]')[[1]][3], 'Cage' )[[1]][1]
  test.cage = as.character(d %>% distinct(cage))
  chrom = 'Genome'
  median.target = median(d$dAF.1_5.target)
  median.matched = median(d$dAF.1_5.matched)
  pval = t.test(d$dAF.1_5.target, d$dAF.1_5.matched)$p.val
  stats = as.data.frame(cbind(id.treat, test.treat, test.cage,chrom, median.target, median.matched, pval))
  data.stats.genome = rbind(data.stats.genome, stats)
}
#chromsome-level stats:
data.stats.chr = data.frame()
files = list.files(pattern = 'Shifts')
for (f in files){
  d = read.csv(f)
  id.treat = strsplit(f, '[.]')[[1]][1]
  test.treat = strsplit(strsplit(f, '[.]')[[1]][3], 'Cage' )[[1]][1]
  test.cage = as.character(d %>% distinct(cage))
  stats.cage = data.frame()
  for (chr in as.character(unique(d$chrom))){    
    d.c = d %>% filter(chrom == chr)
    c = chr
    median.target = median(d.c$dAF.1_5.target)
    median.matched = median(d.c$dAF.1_5.matched)
    pval = t.test(d.c$dAF.1_5.target, d.c$dAF.1_5.matched)$p.val
    stats = as.data.frame(cbind(id.treat, test.treat, test.cage, c, median.target, median.matched, pval))
    stats.cage = rbind(stats.cage, stats)    
  }
  data.stats.chr = rbind(data.stats.chr, stats.cage)
}
data.stats.chr = data.stats.chr %>% rename(chrom = c)
stats.meta = rbind(data.stats.genome, data.stats.chr)
stats.meta = unique(stats.meta)
stats.meta = stats.meta %>% mutate(median.target = as.numeric(as.character(median.target)),
                                   median.matched = as.numeric(as.character(median.matched)),
                                   pval = as.numeric(as.character(pval)))

stats.meta = stats.meta %>% rowwise() %>% mutate(FDR = p.adjust(pval))
stats.meta = stats.meta %>% mutate(sig = if_else(FDR < 0.05 & abs(median.target) > abs(median.matched), 'sig', 'non-sig'))
stats.meta = stats.meta %>% rowwise() %>%
  mutate(dynamic = case_when(sig == 'sig' & abs(median.target > 0) ~ 'parallel',
                             sig == 'sig' & abs(median.target < 0) ~ 'fluctuating',
                             sig == 'non-sig' & abs(median.target > 0) ~ 'non-sig',))
write.csv(stats.meta, './GLM.TreatmentSiteTest.Stats.csv', row.names = FALSE)
stats.meta = read.csv('./GLM.TreatmentSiteTest.Stats.csv')
df.stats = read.csv('./Code_AppleVBloom/GLM.TreatmentSiteTest.Stats.csv')
##add N sites to stats data:
df.a = read.csv('./df.sig.A.1_5.csv')
df.a = df.a %>% filter(sigLevel > 1)
nrow(df.a)
df.a %>% group_by(chrom) %>% count
df.b = read.csv('./df.sig.B.1_5.csv')
df.b = df.b %>% filter(sigLevel > 1)
nrow(df.b)
df.b %>% group_by(chrom) %>% count
2203
A grouped_df: 5 × 2
chrom	n
<fct>	<int>
  2L	938
2R	278
3L	372
3R	504
X	111
29303
A grouped_df: 5 × 2
chrom	n
<fct>	<int>
  2L	14925
2R	4374
3L	4871
3R	3660
X	1473
df.stats = df.stats %>% mutate(N.SNPs = case_when(id.treat == 'Apple' & chrom == 'Genome' ~ 2203,
                                                  id.treat == 'Apple' & chrom == '2L' ~ 938,
                                                  id.treat == 'Apple' & chrom == '2R' ~ 278,
                                                  id.treat == 'Apple' & chrom == '3L' ~ 372,
                                                  id.treat == 'Apple' & chrom == '3R' ~ 504,
                                                  id.treat == 'Apple' & chrom == 'X' ~ 111,
                                                  id.treat == 'Bloom' & chrom == 'Genome' ~ 29303,
                                                  id.treat == 'Bloom' & chrom == '2L' ~ 14925,
                                                  id.treat == 'Bloom' & chrom == '2R' ~ 4374,
                                                  id.treat == 'Bloom' & chrom == '3L' ~ 4871,
                                                  id.treat == 'Bloom' & chrom == '3R' ~ 3660,
                                                  id.treat == 'Bloom' & chrom == 'X' ~ 1473))
df.stats = df.stats %>% dplyr::select(id.treat, test.treat, N.SNPs, test.cage, chrom, median.target, median.matched, FDR)
tail(df.stats)
A data.frame: 6 × 8
id.treat	test.treat	N.SNPs	test.cage	chrom	median.target	median.matched	FDR
<fct>	<fct>	<dbl>	<int>	<fct>	<dbl>	<dbl>	<dbl>
  67	Bloom	Apple	1473	5	X	0.01957270	0.00257900	1.305942e-59
68	Bloom	Apple	14925	6	2L	0.02090100	0.00401750	0.000000e+00
69	Bloom	Apple	4374	6	2R	0.00858000	0.00112075	5.947830e-31
70	Bloom	Apple	4871	6	3L	0.01552017	0.00003200	3.918344e-127
71	Bloom	Apple	3660	6	3R	0.01861150	0.00331650	3.965674e-101
72	Bloom	Apple	1473	6	X	0.01494900	0.00052080	2.907537e-54
write.csv(df.stats, './T.Test.Stats.csv', row.names = FALSE)
stats.meta = stats.meta %>% gather(median.target, median.matched, key = type, value = median)
head(stats.meta)
A data.frame: 6 × 10
id.treat	test.treat	test.cage	chrom	pval	FDR	sig	dynamic	type	median
<fct>	<fct>	<int>	<fct>	<dbl>	<dbl>	<fct>	<fct>	<chr>	<dbl>
  1	Apple	Bloom	1	Genome	1.445478e-88	1.445478e-88	sig	parallel	median.target	0.0248795
2	Apple	Bloom	2	Genome	1.946528e-152	1.946528e-152	sig	parallel	median.target	0.0285690
3	Apple	Bloom	3	Genome	7.323402e-180	7.323402e-180	sig	parallel	median.target	0.0350248
4	Apple	Bloom	4	Genome	4.498558e-71	4.498558e-71	sig	parallel	median.target	0.0300626
5	Apple	Bloom	5	Genome	4.104603e-122	4.104603e-122	sig	parallel	median.target	0.0237203
6	Apple	Bloom	6	Genome	1.187323e-111	1.187323e-111	sig	parallel	median.target	0.0288159
stats.meta = stats.meta %>% rowwise() %>%
  mutate(SNPs = strsplit(as.character(type), '[.]')[[1]][2]) %>%
  mutate(Color = case_when(SNPs == 'matched' ~ 'black',
                           SNPs == 'target' & dynamic == 'parallel' ~ 'red',
                           SNPs == 'target'& dynamic == 'fluctuating' ~ 'blue',
                           SNPs == 'target'& dynamic == 'non-sig' ~ 'dark grey'))
stats.meta$chrom = factor(stats.meta$chrom, levels = c('Genome', '2L', '2R', '3L', '3R', 'X'))
stats.meta = stats.meta %>% mutate(SNPs = if_else(SNPs == 'matched', 'Matched Control', 'Target'))
p = ggplot(stats.meta %>% filter(id.treat == 'Apple'), aes(x = SNPs, y = median, shape = SNPs, colour = Color))+
  geom_jitter(width = 0.2, size = 4) +
  scale_colour_identity() + 
  facet_wrap(~chrom) +
  scale_shape_manual(values = c(4, 19)) +
  ylim(c(-0.02, 0.07)) +
  scale_x_discrete(limits = c('Target', 'Matched Control')) +
  theme_bw(base_size = 20) +
  xlab('')
p
ggsave('../Figures/AppleSites.BloomCageTest.pdf', p, height = 10, width = 15)
Warning message:
  “Removed 2 rows containing missing values (`geom_point()`).”
Warning message:
  “Removed 2 rows containing missing values (`geom_point()`).”

p = ggplot(stats.meta %>% filter(id.treat == 'Bloom'), aes(x = SNPs, y = median, shape = SNPs, colour = Color))+
  geom_jitter(width = 0.2, size = 4) +
  scale_colour_identity() + 
  facet_wrap(~chrom) +
  ylim(c(-0.02, 0.07)) +
  scale_shape_manual(values = c(4, 19)) +
  scale_x_discrete(limits = c('Target', 'Matched Control')) +
  theme_bw(base_size = 20) +
  xlab('')
p
ggsave('../Figures/BloomSites.AppleCageTest.pdf', p, height = 10, width = 15)
Warning message:
  “Removed 3 rows containing missing values (`geom_point()`).”
Warning message:
  “Removed 3 rows containing missing values (`geom_point()`).”

Allele frequency correlations
setwd('../AppleBloomMSRevisions/')
load('./orch2020_filtered.RData')
##segregate data into Apple and Bloomington to average frequencies per treatment
df = cbind(samps, t(afmat))
df.a = df %>% filter(treatment == "A")
df.b = df %>% filter(treatment == "B")

samps.a = df.a[,1:ncol(samps)]
afmat.a = df.a[,-c(1:ncol(samps))]
afmat.a = as.data.frame(t(afmat.a ))


samps.b = df.b[,1:ncol(samps)]
afmat.b = df.b[,-c(1:ncol(samps))]
afmat.b = as.data.frame(t(afmat.b ))
#Also get sites significant in each treatment for assessement
df.sig.a = read.csv('./df.sig.Matched.A.1_5.csv')
df.sig.a = df.sig.a %>% filter(sigLevel > 1)

df.sig.b = read.csv('./df.sig.Matched.B.1_5.csv')
df.sig.b = df.sig.b %>% filter(sigLevel > 1)

sites.A = df.sig.a %>% dplyr::select(chrom, pos)
sites.B = df.sig.b %>% dplyr::select(chrom, pos)

sites.A.matched = df.sig.a %>% dplyr::select(chrom, pos.matched) %>% rename(pos = pos.matched)
sites.B.matched = df.sig.b %>% dplyr::select(chrom, pos.matched)  %>% rename(pos = pos.matched)
#Get allele frequency shifts genome wide and for each set of set of SNPs, by treatment
shifts.genome.a = get_af_shifts(afmat = afmat.a, samps = samps.a, comparisons = c('1_5'))
shifts.genome.a = cbind(sites, shifts.genome.a)
shifts.genome.b = get_af_shifts(afmat = afmat.b, samps = samps.b, comparisons = c('1_5'))
shifts.genome.b = cbind(sites, shifts.genome.b)

shifts.ASites.a = left_join(sites.A, shifts.genome.a)
shifts.ASites.a.matched = left_join(sites.A.matched, shifts.genome.a)

shifts.ASites.b = left_join(sites.A, shifts.genome.b)
shifts.ASites.b.matched = left_join(sites.A.matched, shifts.genome.b)

shifts.BSites.a = left_join(sites.B, shifts.genome.a)
shifts.BSites.a.matched = left_join(sites.B.matched, shifts.genome.a)

shifts.BSites.b = left_join(sites.B, shifts.genome.b)
shifts.BSites.b.matched = left_join(sites.B.matched, shifts.genome.b)
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
Joining with `by = join_by(chrom, pos)`
#quick look at point estimates of correlations between treatments
cor.test(shifts.genome.a$dAF.1_5, shifts.genome.b$dAF.1_5 )$estimate[[1]][1]

cor.test(shifts.ASites.a$dAF.1_5, shifts.ASites.b$dAF.1_5 )$estimate[[1]][1]
cor.test(shifts.ASites.a.matched$dAF.1_5, shifts.ASites.b.matched$dAF.1_5 )$estimate[[1]][1]

cor.test(shifts.BSites.a$dAF.1_5, shifts.BSites.b$dAF.1_5 )$estimate[[1]][1]
cor.test(shifts.BSites.a.matched$dAF.1_5, shifts.BSites.b.matched$dAF.1_5 )$estimate[[1]][1]
0.418777805356415
0.863919610744922
0.430926912004359
0.823382116270171
0.419951363917199
##Compute genome-wide correlations and correlations by chromosome between treatments:
chroms = c( '2L', '2R', '3L', '3R', 'X')

all.sites = data.frame()
for (chr in chroms){
  a.shifts = shifts.genome.a %>% filter(chrom == chr)
  b.shifts = shifts.genome.b %>% filter(chrom == chr)
  cor = cor.test(a.shifts$dAF.1_5, b.shifts$dAF.1_5)$estimate
  d = cbind(chr, cor)
  all.sites = rbind(all.sites, d)
  
  
}
all.gen = as.data.frame(cbind('Genome', cor.test(shifts.genome.a$dAF.1_5, shifts.genome.b$dAF.1_5 )$estimate[[1]][1]))
names(all.gen) = c('chr', 'cor')
all.sites = rbind(all.sites, all.gen)
head(all.sites)
A data.frame: 6 × 2
chr	cor
<fct>	<fct>
  cor	2L	0.586217398438326
cor1	2R	0.447392295901473
cor2	3L	0.281969757852913
cor3	3R	0.38686307202988
cor4	X	0.198579750786953
1	Genome	0.418777805356415
#permutations of genome-wide and by chrom
all.sites.perm = data.frame()
for (i in 1:100){
  for (chr in chroms){
    a.shifts = shifts.genome.a %>% filter(chrom == chr)
    a.shifts = a.shifts[sample(nrow(a.shifts)),]    #permute a shifts
    b.shifts = shifts.genome.b %>% filter(chrom == chr)
    cor = cor.test(a.shifts$dAF.1_5, b.shifts$dAF.1_5)$estimate
    d = cbind(i, chr, cor)
    all.sites.perm = rbind(all.sites.perm, d)
  }
  
}
all.sites.perm = all.sites.perm %>% mutate(cor = as.numeric(as.character(cor)))
all.sites.perm.genome = data.frame()
for (i in 1:100){
  chr = 'Genome'
  a.shifts = shifts.genome.a[sample(nrow(shifts.genome.a)),]    #permute a shifts
  b.shifts = shifts.genome.b 
  cor = cor.test(a.shifts$dAF.1_5, b.shifts$dAF.1_5)$estimate
  d = cbind(i, chr, cor)
  all.sites.perm.genome = rbind(all.sites.perm.genome, d)
}

all.sites.perm.genome = all.sites.perm.genome %>% mutate(cor = as.numeric(as.character(cor)))
all.sites.perm = rbind(all.sites.perm, all.sites.perm.genome)
##Plot chromosome and genome observations relative to null distribution
ggplot(all.sites.perm %>% filter(chr == '2L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.02, 0.99)

ggplot(all.sites.perm %>% filter(chr == '2L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == '2L'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('2L') +
  xlab('Pearson Correlation Coefficient' )

#plots:
p.2l = ggplot(all.sites.perm %>% filter(chr == '2L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == '2L'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('2L') +
  xlab('Pearson Correlation Coefficient' )

p.2r = ggplot(all.sites.perm %>% filter(chr == '2R'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == '2R'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('2R') +
  xlab('Pearson Correlation Coefficient' )

p.3l = ggplot(all.sites.perm %>% filter(chr == '3L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == '3L'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('3L') +
  xlab('Pearson Correlation Coefficient' )

p.3r = ggplot(all.sites.perm %>% filter(chr == '3R'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == '3R'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('3R') +
  xlab('Pearson Correlation Coefficient' )

p.x = ggplot(all.sites.perm %>% filter(chr == 'X'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == 'X'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('X') +
  xlab('Pearson Correlation Coefficient' )

p.gen = ggplot(all.sites.perm %>% filter(chr == 'Genome'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(-0.04, 0.99) +
  geom_vline(xintercept = as.numeric(as.character((all.sites %>% filter(chr == 'Genome'))$cor)), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('Genome') +
  xlab('Pearson Correlation Coefficient' )
library(cowplot)
p= plot_grid(p.gen, p.2l, p.2r, p.3l, p.3r, p.x,  nrow = 1) 
p
ggsave('./AllSites.Correlations.pdf', height = 14, width = 32)

#are observed correlations significant relative to these permutations?
##get 99.9% quantile 

perm.results = all.sites.perm %>% 
  group_by(chr) %>%  
  summarise(quantile = scales::percent(c(0.999)),
            cor = quantile(cor, c(0.999))) %>% rename(null.cor = cor) %>% mutate(null.cor = as.numeric(as.character(null.cor)))
all.sites.res = left_join(perm.results, all.sites)
all.sites.res %>% mutate(significant = if_else(as.numeric(cor) > null.cor, 'yes', 'no')) %>% dplyr::select(-quantile)
Joining with `by = join_by(chr)`
A tibble: 6 × 4
chr	null.cor	cor	significant
<fct>	<dbl>	<fct>	<chr>
  2L	0.003272145	0.586217398438326	yes
2R	0.003609942	0.447392295901473	yes
3L	0.002746358	0.281969757852913	yes
3R	0.002991775	0.38686307202988	yes
X	0.005500572	0.198579750786953	yes
Genome	0.001397263	0.418777805356415	yes
head(all.sites.res)
A tibble: 6 × 4
chr	quantile	null.cor	cor
<fct>	<chr>	<dbl>	<fct>
  2L	100%	0.003272145	0.586217398438326
2R	100%	0.003609942	0.447392295901473
3L	100%	0.002746358	0.281969757852913
3R	100%	0.002991775	0.38686307202988
X	100%	0.005500572	0.198579750786953
Genome	100%	0.001397263	0.418777805356415
##Plot chromosome and genome observations relative to null distribution
##Now assess correlations of Apple and Bloomington sites, relative to random samples from the genome
#get number of
n.sites.A = sites.A %>% group_by(chrom) %>% count
n.sites.B = sites.B %>% group_by(chrom) %>% count
#Apple sites by chrom:
chroms = c( '2L', '2R', '3L', '3R', 'X')

apple.sites = data.frame()
for (chr in chroms){
  a.shifts = shifts.ASites.a %>% filter(chrom == chr)
  b.shifts = shifts.ASites.b %>% filter(chrom == chr)
  cor = cor.test(a.shifts$dAF.1_5, b.shifts$dAF.1_5)$estimate
  d = cbind(chr, cor)
  apple.sites = rbind(apple.sites, d)
  
  
}
apple.gen = as.data.frame(cbind('Genome', cor.test(shifts.ASites.a$dAF.1_5, shifts.ASites.b$dAF.1_5 )$estimate[[1]][1]))
names(apple.gen) = c('chr', 'cor')
apple.sites = rbind(apple.sites, apple.gen)
apple.sites
A data.frame: 6 × 2
chr	cor
<fct>	<fct>
  cor	2L	0.929629733374436
cor1	2R	0.899586220176079
cor2	3L	0.686901015524171
cor3	3R	0.819187946435207
cor4	X	0.639913831304598
1	Genome	0.863919610744922
##random sample of sites from genome of same length as apple sites

apple.sites.null = data.frame()
for (i in 1:100){
  for (chr in chroms){
    n.sites = n.sites.A %>% filter(chrom == chr) %>% pull(n)
    a.shifts = shifts.genome.a %>% filter(chrom == chr) %>% sample_n(n.sites)
    sites.null = a.shifts %>% dplyr::select(chrom, pos)
    b.shifts = left_join(sites.null, shifts.genome.b)
    b.shifts = na.omit(b.shifts)
    cor = cor.test(a.shifts$dAF.1_5, b.shifts$dAF.1_5)$estimate
    d = cbind(i, chr, cor)
    apple.sites.null = rbind(apple.sites.null, d)
  }
}

#generate statistics - what quantile do the observed values fall in relative to the random sample distributions?
apple.sites.null = apple.sites.null %>% mutate(cor = as.numeric(as.character(cor)))
apple.sites  = apple.sites %>% mutate(cor = as.numeric(as.character(cor)))
apple.sites$type = 'Observed'
apple.sites.null$type = 'Random'
#emprical pval by chrom
sig.df = data.frame()
for (chrom in chroms){
  corr = apple.sites %>% filter(chr == chrom) %>% pull(cor)
  null.dist = (apple.sites.null %>% filter(chr == chrom))$cor
  p.val = 1 - ecdf(null.dist)(corr)
  sig.df = rbind(sig.df, p.val)
  
}
names(sig.df) = c('p.val')
#genome-wide empirical pval
corr = apple.sites %>% filter(chr == 'Genome') %>% pull(cor)
null.dist = apple.sites.null$cor
p.val = 1 - ecdf(null.dist)(corr)
genome.p = 1 - ecdf(null.dist)(corr)
sig.df = rbind(sig.df, genome.p)

apple.sites = cbind(apple.sites, sig.df)

head(apple.sites.null)
A data.frame: 6 × 4
i	chr	cor	type
<fct>	<fct>	<dbl>	<chr>
  cor	1	2L	0.54860675	Random
cor1	1	2R	0.48287352	Random
cor2	1	3L	0.34467865	Random
cor3	1	3R	0.31926601	Random
cor4	1	X	0.05784908	Random
cor5	2	2L	0.54914491	Random
#plots:
a.2l = ggplot(apple.sites.null %>% filter(chr == '2L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = apple.sites %>% filter(chr == '2L') %>% pull(cor), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('2L') +
  xlab('Pearson Correlation Coefficient' )

a.2r = ggplot(apple.sites.null %>% filter(chr == '2R'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = apple.sites %>% filter(chr == '2R') %>% pull(cor), color = 'red' ) + 
  guides(fill="none")  +
  ggtitle('2R') +
  xlab('Pearson Correlation Coefficient' )

a.3l = ggplot(apple.sites.null %>% filter(chr == '3L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = apple.sites %>% filter(chr == '3L') %>% pull(cor), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('3L') +
  xlab('Pearson Correlation Coefficient' )

a.3r = ggplot(apple.sites.null %>% filter(chr == '3R'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 23) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = apple.sites %>% filter(chr == '3R') %>% pull(cor), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('3R') +
  xlab('Pearson Correlation Coefficient' )

a.x = ggplot(apple.sites.null %>% filter(chr == 'X'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = apple.sites %>% filter(chr == 'X') %>% pull(cor), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('X') +
  xlab('Pearson Correlation Coefficient' )

a.gen = ggplot(apple.sites.null, aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = apple.sites %>% filter(chr == 'Genome') %>% pull(cor), color = 'red' ) +
  guides(fill="none")  +
  ggtitle('Genome') +
  xlab('Pearson Correlation Coefficient' )
a.gen
Warning message:
  “Removed 6 rows containing non-finite values (`stat_density()`).”

library(cowplot)
p= plot_grid(a.gen, a.2l, a.2r, a.3l, a.3r, a.x,  nrow = 1) 
p
ggsave('./ASites.Correlations.pdf', height = 14, width = 32)
Warning message:
  “Removed 5 rows containing non-finite values (`stat_density()`).”
Warning message:
  “Removed 5 rows containing non-finite values (`stat_density()`).”

# now look at Bloomington sites by chrom:
chroms = c( '2L', '2R', '3L', '3R', 'X')

bloom.sites = data.frame()
for (chr in chroms){
  a.shifts = shifts.BSites.a %>% filter(chrom == chr)
  b.shifts = shifts.BSites.b %>% filter(chrom == chr)
  cor = cor.test(a.shifts$dAF.1_5, b.shifts$dAF.1_5)$estimate
  d = cbind(chr, cor)
  bloom.sites = rbind(bloom.sites, d)
  
  
}
bloom.gen = as.data.frame(cbind('Genome', cor.test(shifts.BSites.a$dAF.1_5, shifts.BSites.b$dAF.1_5 )$estimate[[1]][1]))
names(bloom.gen) = c('chr', 'cor')
bloom.sites = rbind(bloom.sites, bloom.gen)
bloom.sites
A data.frame: 6 × 2
chr	cor
<fct>	<fct>
  cor	2L	0.878072765007211
cor1	2R	0.843933007636718
cor2	3L	0.651665639073897
cor3	3R	0.737474213459068
cor4	X	0.708586720338408
1	Genome	0.823382116270171
##random sample of sites from genome of same length as apple sites

bloom.sites.null = data.frame()
for (i in 1:100){
  for (chr in chroms){
    n.sites = n.sites.B %>% filter(chrom == chr) %>% pull(n)
    b.shifts = shifts.genome.b %>% filter(chrom == chr) %>% sample_n(n.sites)
    sites.null = b.shifts %>% dplyr::select(chrom, pos)
    a.shifts = left_join(sites.null, shifts.genome.a)
    a.shifts = na.omit(a.shifts)
    cor = cor.test(b.shifts$dAF.1_5, a.shifts$dAF.1_5)$estimate
    d = cbind(i, chr, cor)
    bloom.sites.null = rbind(bloom.sites.null, d)
  }
}

head(bloom.sites.null)
A data.frame: 6 × 3
i	chr	cor
<fct>	<fct>	<fct>
  cor	1	2L	0.585806152219517
cor1	1	2R	0.439225341892664
cor2	1	3L	0.278572242304196
cor3	1	3R	0.381340788247217
cor4	1	X	0.225503264650847
cor5	2	2L	0.582339931617618
#generate statistics - what quantile do the observed values fall in relative to the random sample distributions?
bloom.sites.null = bloom.sites.null %>% mutate(cor = as.numeric(as.character(cor)))
bloom.sites  = bloom.sites %>% mutate(cor = as.numeric(as.character(cor)))
bloom.sites$type = 'Observed'
bloom.sites.null$type = 'Random'

#emprical pval by chrom
sig.df = data.frame()
for (chrom in chroms){
  corr = apple.sites %>% filter(chr == chrom) %>% pull(cor)
  null.dist = (bloom.sites.null %>% filter(chr == chrom))$cor
  p.val = 1 - ecdf(null.dist)(corr)
  sig.df = rbind(sig.df, p.val)
  
}
names(sig.df) = c('p.val')
#genome-wide empirical pval
corr = bloom.sites %>% filter(chr == 'Genome') %>% pull(cor)
null.dist = bloom.sites.null$cor
p.val = 1 - ecdf(null.dist)(corr)
genome.p = 1 - ecdf(null.dist)(corr)
sig.df = rbind(sig.df, genome.p)

bloom.sites = cbind(bloom.sites, sig.df)
#plots:
#plots:
b.2l = ggplot(bloom.sites.null %>% filter(chr == '2L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 15) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == '2L') %>% pull(cor), color = 'red' ) +
  guides(fill="none") 
b.2l

b.2l = ggplot(bloom.sites.null %>% filter(chr == '2L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == '2L') %>% pull(cor), color = 'red' ) +
  ggtitle('2L') +
  guides(fill="none")+
  xlab('Pearson Correlation Coefficient' )

b.2r = ggplot(bloom.sites.null %>% filter(chr == '2R'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == '2R') %>% pull(cor), color = 'red' ) +
  ggtitle('2R') +
  guides(fill="none") +
  xlab('Pearson Correlation Coefficient' )

b.3l = ggplot(bloom.sites.null %>% filter(chr == '3L'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == '3L') %>% pull(cor), color = 'red' ) +
  ggtitle('3L') +
  guides(fill="none") +
  xlab('Pearson Correlation Coefficient' )


b.3r = ggplot(bloom.sites.null %>% filter(chr == '3R'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == '3R') %>% pull(cor), color = 'red' ) +
  ggtitle('3R') +
  guides(fill="none") +
  xlab('Pearson Correlation Coefficient' )


b.X = ggplot(bloom.sites.null %>% filter(chr == 'X'), aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == 'X') %>% pull(cor), color = 'red' ) +
  ggtitle('X') +
  guides(fill="none") +
  xlab('Pearson Correlation Coefficient' )


b.gen = ggplot(bloom.sites.null, aes(x = cor, fill = 'grey')) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c('grey')) +
  theme_bw(base_size = 23) +
  xlim(0, 0.99) +
  geom_vline(xintercept = bloom.sites %>% filter(chr == 'Genome') %>% pull(cor), color = 'red' ) +
  ggtitle('Genome') +
  guides(fill="none") +
  xlab('Pearson Correlation Coefficient' )
p = plot_grid(b.gen, b.2l, b.2r, b.3l, b.3r, b.X, nrow = 1)
ggsave('./BSites.Correlations.pdf', p, height = 14, width = 32)
