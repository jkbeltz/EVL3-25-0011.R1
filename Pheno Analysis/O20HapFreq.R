hapfreq=read_csv("HapFreqs.csv")
View(hapfreq)

hapfreq.sum <- hapfreq %>% 
  group_by(treatment, tpt, replicate) %>% 
  dplyr::summarize(tpt = tpt,
          treatment = treatment,
          count1 = length(Haplotype[which(Hap.Freq > .01)]),
          count05 = length(Haplotype[which(Hap.Freq > .015)]))
hapfreq.sum=unique(hapfreq.sum)
View(hapfreq.sum)
