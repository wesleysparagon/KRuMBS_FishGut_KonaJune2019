#Analyze the distribution of abundances of ASVs across the various subsections. Use the results of this to pick out ASVs that peak in abundance in the various sections.
mean.gut_subsection=function(x,y) { #function of inputs x where x=row in relabund.t.fish.F4.F8.cull that corresponds to an ASV column and y=character string corresponding to Gut_SubSection
  means=mean(x[relabund.t.fish.F4.F8.cull$Gut_SubSection==y]) 
  return(means)}
relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection=cbind(relabund.t.fish.F4.F8.cull[,colnames(relabund.t.fish.F4.F8.cull) %in% rownames(mixed.mod.sig.otu.Gut_SubSection)],relabund.t.fish.F4.F8.cull[,790:799]) #make a new relabund df that only has the significant gut subsection otus.
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$ST.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"ST")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$PC.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"PC")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$GI_1.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"GI_1")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$GI_2.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"GI_2")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$GI_3.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"GI_3")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$GI_4.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"GI_4")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$HG_1.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"HG_1")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$HG_2.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"HG_2")
}
for (i in 1:length(rownames(mixed.mod.sig.otu.Gut_SubSection))) {
  mixed.mod.sig.otu.Gut_SubSection$HG_3.mean[i]=mean.gut_subsection(relabund.t.fish.F4.F8.cull.sig.otu.gut_subsection[,i],"HG_3")
}
write.csv(mixed.mod.sig.otu.Gut_SubSection,file="mixed.mod.sig.otu.Gut_SubSection.csv") #export as CSV







#run post hoc testing to identify which Gut_SubSections are significantly different from each other
summary(glht(mixed.mod.relabund.t.fish.F4.F8.cull.trans$Otu00004, linfct = mcp(Gut_SubSection = "Tukey")), test = adjusted("holm")) 






#visualize differentially abundant ASVs
otu.4.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00004,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+
  ylim(0,.85)+
  theme_minimal()+
  theme(legend.position="none")+ #remove legend, crop a legend out and combine it with multiple plots
  #remove legend, add 1 in later in power point for multiple plots.
  ylab("Relative Abundance")+
  ggtitle("Pasteurellaceae")

#significance label code if necessary: geom_text(x=1,y=.8,label="A")

#Highly abundant in the stomach while not found in any other gut sections at high abundance. Facultative anaerobe suggests potentially microaerophillic environment.

summary(glht(mixed.mod.relabund.t.fish.F4.F8.cull.trans$Otu00309, linfct = mcp(Gut_SubSection = "Tukey")), test = adjusted("holm")) #run post hoc testing to identify which Gut_SubSections are significantly different from each other

otu.309.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00309,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+
  ylim(0,.015)+
  theme_minimal()+
  theme(legend.position="none")+ #remove legend, crop a legend out and combine it with multiple plots in power point. This can be done more eloquently by merging otu relabund column data that are to be visualized (and ONLY that data) into a df and plotting with facet().
  ylab("Relative Abundance")+
  ggtitle("Rivularia_PCC-7116 (Cyanobacteria)")

#significance label code if necessary: ?

plot_grid(otu.4.relabund.plot,otu.309.relabund.plot) #plot these 2 graphs together

#Cyanobacteria found growing colonially on submerged stones (among other locations). Potentially a food source/incidental food source of grazing Nenue? Reduction in relative abundance through the gut suggests it is being digested.

ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00634,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Hyphomonas")

#Trend potentially driven by outlier. However, Hyphomonas is still important since it is a genus of Rhodobacterales that is a obligate aerobe, lending creedence to the hypothesis that the Nenue stomach is aerobic.

#Now let's visualize ASVs that are abundant in the early GI (PC-GI2)

summary(glht(mixed.mod.relabund.t.fish.F4.F8.cull.trans$Otu00015, linfct = mcp(Gut_SubSection = "Tukey")), test = adjusted("holm")) #run post hoc testing to identify which Gut_SubSections are significantly different from each other

otu.15.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00015,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="none")+
  ylim(0,.4)+
  ylab("Relative Abundance")+
  ggtitle("Peptostreptococcaceae")
  
#significance code if necessary:annotate(geom="text",x=2:9,y=c(.05,.375,.275,.075,.075,.05,.05,.05),label=c("B/","A","A/B","B/","B/","B/","B/","B/"))

#Anaerobic Clostridia fermenters. Suggests potentially anaerobic early GI with fermentation as a dominant process.

summary(glht(mixed.mod.relabund.t.fish.F4.F8.cull.trans$Otu00003, linfct = mcp(Gut_SubSection = "Tukey")), test = adjusted("holm")) #run post hoc testing to identify which Gut_SubSections are significantly different from each other

otu.3.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00003,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Vibrionaceae")

#significance code if necessary:annotate(geom="text",x=1:7,y=c(.175,.175,.225,.35,.125,.05,.05),label=c("A/","A/","A/","A","A/","A/","A/"))

plot_grid(otu.15.relabund.plot,otu.3.relabund.plot) #plot these 2 graphs together

#Vibrionaceae peaking in the early GI (GI_1-2). Facultative anaerobes, fermenters.

#Now let's visualize ASVs that are abundant in the early GI (GI3-4).

#No longer going to run post hocs until I have a better idea how to display them.

otu.5.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00005,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Turicibacter (Erysipelotrichaceae)")

otu.33.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00033,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Brevinema (Brevinemataceae)")

otu.81.relabund.plot=
  ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00081,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Clostridiaceae") #Generally anaerobic

plot_grid(otu.5.relabund.plot,otu.33.relabund.plot,otu.81.relabund.plot)

#Now let's visualize ASVs that are abundant in the hindgut (HG1-3)

otu.19.relabund.plot=
ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00019,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Desulfovibrio")

#Sulfate reducers higher in HG. Obligate anaerobes but can actually be aerobically tolerant.

otu.25.relabund.plot=
ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00025,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Ruminococcaceae")

#Obligate anaerobes.

otu.40.relabund.plot=
ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00040,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Akkermansia")

#chemoorganotroph obligate anaerobes.

otu.60.relabund.plot=
ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00062,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Clostridiales vadinBB60 group")

#don't include this one for now
ggplot(relabund.t.fish.F4.F8.cull,aes(y=Otu00047,x=Gut_SubSection,fill=Gut_SubSection))+
  scale_fill_manual(values=viridis(n=9,option="D"))+
  geom_boxplot()+ #Ask Craig if I should remove outliers using geom_boxplot(outlier.shape=NA)
  theme_minimal()+
  theme(legend.position="none")+
  ylab("Relative Abundance")+
  ggtitle("Alistipes (Rikenellaceae)")

plot_grid(otu.19.relabund.plot,otu.25.relabund.plot,otu.40.relabund.plot,otu.60.relabund.plot)




#compare number of significant ASVs in full model vs. models without species
#now run a mixed model on the transformed ASV data. Have relative abundance be a function of Species, Gut_SubSection, and Fish_Number be a random effect. Interaction term has been removed since it was not significant in the multivariate analysis. Could re-run model with Species also as a random effect.
mixed.mod=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Species+Gut_SubSection+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod) #apply mixed.mod only to columns containing ASV data. Save as new object.

#Now let's inspect the results of the mixed models.
mixed.mod.singularity=lapply(mixed.mod.relabund.t.fish.F4.F8.cull.trans,FUN=isSingular) #test each mixed model for singularity, save in list. After inspecting the list it looks like most models are not singular.
anova.kw=function(x) { #make a function that performs anova (Kenward-Roger)
  results=anova(x,ddf="Kenward-Roger")
  return(results)
}
mixed.mod.results=lapply(mixed.mod.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod.pvals=as.data.frame(t(as.data.frame(mixed.mod.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).
#Work up p values
mixed.mod.pvals$Gut_SubSection.padjust=p.adjust(mixed.mod.pvals$Gut_SubSection,method="BH") #Adjust p values for multiple comparisons
mixed.mod.pvals$Species.padjust=p.adjust(mixed.mod.pvals$Species,method="BH") #Adjust p values for multiple comparisons
rownames(mixed.mod.pvals)=colnames(relabund.t.fish.F4.F8.cull.trans)[1:692] #Change rownames to correspond with tested OTU names.
mixed.mod.sig.otu.Gut_SubSection=mixed.mod.pvals[mixed.mod.pvals$Gut_SubSection.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_SubSection pvals. Save as new df.

mixed.mod1=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Gut_SubSection+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod1.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod1) #apply mixed.mod only to columns containing ASV data. Save as new object.
mixed.mod1.results=lapply(mixed.mod1.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod1.pvals=as.data.frame(t(as.data.frame(mixed.mod1.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).
#Work up p values
colnames(mixed.mod1.pvals)="Gut_SubSection"
mixed.mod1.pvals$Gut_SubSection.padjust=p.adjust(mixed.mod1.pvals$Gut_SubSection,method="BH") #Adjust p values for multiple comparisons
mixed.mod1.sig.otu.Gut_SubSection=mixed.mod1.pvals[mixed.mod1.pvals$Gut_SubSection.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_SubSection pvals. Save as new df.

mixed.mod2=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Gut_SubSection+(1|Species/Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod2.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod2) #apply mixed.mod only to columns containing ASV data. Save as new object.
mixed.mod2.results=lapply(mixed.mod2.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod2.pvals=as.data.frame(t(as.data.frame(mixed.mod2.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).
#Work up p values
colnames(mixed.mod2.pvals)="Gut_SubSection"
mixed.mod2.pvals$Gut_SubSection.padjust=p.adjust(mixed.mod2.pvals$Gut_SubSection,method="BH") #Adjust p values for multiple comparisons
mixed.mod2.sig.otu.Gut_SubSection=mixed.mod2.pvals[mixed.mod2.pvals$Gut_SubSection.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_SubSection pvals. Save as new df.

#now check to see if the padjusted pvals are equal in the 3 models
mixed.mod.pvals$Gut_SubSection.padjust==mixed.mod1.pvals$Gut_SubSection.padjust #no values are the same
mixed.mod.pvals$Gut_SubSection.padjust==mixed.mod2.pvals$Gut_SubSection.padjust #no values are the same
mixed.mod1.pvals$Gut_SubSection.padjust==mixed.mod2.pvals$Gut_SubSection.padjust #no values are the same

#it looks like the full model has 285 significant ASVs, whereas the reduced models (both of them) have 283 significant ASVs.

