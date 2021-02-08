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




#Here is all the code to analyze Gut_Zone, which was originally called Gut_Zone but then had it's name changed

permanova.gut.section=adonis.II(weighted.unifrac.dist~Species*Gut_Zone,by="margin",strata=metadata.fish.1$Fish_Number,data=metadata.fish.1,permutations=999)
permanova.gut.section #Gut_Subsection and Species are both significant, no significant interaction.

###ERRORS IN THIS PLOT
plot(nmds.weighted.unifrac.2,type="n") #Set up the NMDS plot environment
points(nmds.weighted.unifrac.2, #plot each sample as a point
       cex=2,
       pch=pointvec[metadata.fish.2.unifrac$Species], #point shape corresponds to fish Species
       col=colvec[metadata.fish.2.unifrac$Gut_Zone]) #point color corresponds to Gut_Zone
ordiellipse(nmds.weighted.unifrac.2,metadata.fish.2.unifrac$Gut_,kind="sd",draw="polygon",alpha=c(0,0),lwd=2.5,border=colvec) #add SD ellipses, colored by Gut_Zone
legend(.2,-.22,legend = c("Stomach","Fore/midgut","Hindgut","K. cinerascens","K. hawaiiensis","K. vaigiensis"),cex=1,y.intersp=.25,bg="transparent",bty="o",col=c(colvec,"black","black","black"),pch=c(16,16,16,15,16,17),pt.cex=2) #Add legend. Save image as 8.5x11 PDF (landscape)

#Calculate means and SE of the sobs and shannon evenness alpha diversity metrics for each Gut_Zone and save as a df
alpha.diversity.2.mean.Gut_Zone=aggregate(alpha.diversity.fish.2$sobs,by=list(alpha.diversity.fish.2$Gut_Zone),FUN=mean)$x #aggregate (mean) sobs by Gut_Zone, save in new df.
alpha.diversity.2.mean.Gut_Zone=as.data.frame(alpha.diversity.2.mean.Gut_Zone)
colnames(alpha.diversity.2.mean.Gut_Zone)[1]="sobs" #rename column 1 to sobs.
alpha.diversity.2.mean.Gut_Zone$SE_sobs=aggregate(alpha.diversity.fish.2$sobs,by=list(alpha.diversity.fish.2$Gut_Zone),FUN=stderror)$x #aggregate (stderror) sobs by Gut_Zone, save as a new column in df.
#repeat this for shannoneven
alpha.diversity.2.mean.Gut_Zone$shannoneven=aggregate(alpha.diversity.fish.2$shannoneven,by=list(alpha.diversity.fish.2$Gut_Zone),FUN=mean)$x
alpha.diversity.2.mean.Gut_Zone$SE_shannoneven=aggregate(alpha.diversity.fish.2$shannoneven,by=list(alpha.diversity.fish.2$Gut_Zone),FUN=stderror)$x
#now add in a Gut_Zone column
alpha.diversity.2.mean.Gut_Zone$Gut_Zone=aggregate(alpha.diversity.fish.2$shannoneven,by=list(alpha.diversity.fish.2$Gut_Zone),FUN=stderror)$Group.1

#Model the sobs data in response to Gut_seciton, fish species, and fish ID
adiv.mod.gut.section.sobs=lmer(sobs~Gut_Zone+Species+(1|Fish_Number),data=alpha.diversity.fish.2)
summary(adiv.mod.gut.section.sobs)
anova(adiv.mod.gut.section.sobs,ddf="Kenward-Roger") #Gut_Zone is significant
lsmeans(adiv.mod.gut.section.sobs, pairwise ~ Gut_Zone, adjust="tukey") #HG and ST are not significantly different from each other, both are significantly different from GI.

#Model the shannoneven data in response to Gut_seciton, fish species, and fish ID
adiv.mod.gut.section.shannoneven=lmer(shannoneven~Gut_Zone+Species+(1|Fish_Number),data=alpha.diversity.fish.2)
summary(adiv.mod.gut.section.shannoneven)
anova(adiv.mod.gut.section.shannoneven,ddf="Kenward-Roger") #Gut_Zone is significant
lsmeans(adiv.mod.gut.section.shannoneven, pairwise ~ Gut_Zone, adjust="tukey")
#HG is significantly different from ST and GI. ST and GI are not significantly different from each other.

#Graph alpha diversity indices as a function of Gut_Zone
sobs.gut.section.plot=ggplot(alpha.diversity.2.mean.Gut_Zone,aes(x=Gut_Zone,y=sobs,fill=Gut_Zone))+
  ylim(0,675)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=viridis(n=5))+
  geom_errorbar(aes(ymin=sobs-SE_sobs,ymax=sobs+SE_sobs),width=.5)+
  ggtitle("Observed Sequences")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_text(x=1:3,y=c(650,275,550),label=c("A","B","A"))

#now repeat for shannoneven
shannoneven.gut.section.plot=ggplot(alpha.diversity.2.mean.Gut_Zone,aes(x=Gut_Zone,y=shannoneven,fill=Gut_Zone))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=viridis(n=3))+
  geom_errorbar(aes(ymin=shannoneven-SE_shannoneven,ymax=shannoneven+SE_shannoneven),width=.5)+
  ggtitle("Shannon Eveness")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_text(x=1:3,y=c(.625,.625,.84),label=c("A","A","B"))+
  ylim(0,.85)

plot_grid(sobs.gut.section.plot,shannoneven.gut.section.plot) #use the plot_grid function to plot the above 4 plots (par(mfrow=())) doesnt work with ggplot2. Export as PNG (800 x 455).

dispersion.gut.section=betadisper(weighted.unifrac.2.dist,type="centroid",group=metadata.fish.2.unifrac$Gut_Zone)
dispersion.gut.section
hist(dispersion.gut.section$distances) #asses the distribution of the dispersion data
shapiro.test(dispersion.gut.section$distances) #data is not normal and needs to be transformed.
dispersion.gut.section$sqrt.distances=sqrt(dispersion.gut.section$distances) #square root transform the distance data
hist(dispersion.gut.section$sqrt.distances) #asses distribution of sqrt transformed dispersion data
shapiro.test(dispersion.gut.section$sqrt.distances) #the data appears (barely) normal
dispersion.gut.section1=as.data.frame(dispersion.gut.section[[7]]) #make a new df with just the sqrt distances.
dispersion.gut.section1$distances=dispersion.gut.section[[3]] #add the original distance values as well.
colnames(dispersion.gut.section1)[1]="sqrt.distances"
dispersion.gut.section1$samples=rownames(dispersion.gut.section1) #add in metadata
dispersion.gut.section1$Gut_Zone=metadata.fish.2.unifrac$Gut_Zone
dispersion.gut.section1$Gut_SubSection=metadata.fish.2.unifrac$Gut_SubSection
dispersion.gut.section1$Species=metadata.fish.2.unifrac$Species
dispersion.gut.section1$Fish_Number=metadata.fish.2.unifrac$Fish_Number
dispersion.gut.section2=dispersion.gut.section1[c(-1:-2,-4:-7),] #remove F1-F3 GI and HG samples and save as new df.
#convert dispersion.gut.section2 into df of means and SE so it can be graphed as a barplot.
dispersion.gut.section3=as.data.frame(aggregate(dispersion.gut.section2$distances,by=list(dispersion.gut.section2$Gut_Zone),FUN=mean))
dispersion.gut.section3$se=aggregate(dispersion.gut.section2$distances,by=list(dispersion.gut.section2$Gut_Zone),FUN=stderror)[[2]]
colnames(dispersion.gut.section3)=c("Gut_Zone","distance_to_centroid","se")
#do the same for gut subsection
dispersion.gut.subsection4=as.data.frame(aggregate(dispersion.gut.section2$distances,by=list(dispersion.gut.section2$Gut_SubSection),FUN=mean))
dispersion.gut.subsection4$se=aggregate(dispersion.gut.section2$distances,by=list(dispersion.gut.section2$Gut_SubSection),FUN=stderror)[[2]]
colnames(dispersion.gut.subsection4)=c("Gut_SubSection","distance_to_centroid","se")

#extra NMDS script
#repeat for the PC culled nmds
nmds.scores.2=as.data.frame(scores(nmds.weighted.unifrac.F4.F8)) #extract NMDS scores and save in new df
nmds.scores.2$Gut_SubSection=metadata.fish.F4.F8.unifrac$Gut_SubSection #add the Fish_GutSubsection to the new df
nmds.scores.2$Ranked_Distance_from_Stomach=metadata.fish.F4.F8.unifrac$Ranked_Gut_SubSection_Distance_from_Stomach #Add Ranked_Gut_SubSection_Distance_from_Stomach to df, give it a shorter name
nmds.scores.2$Species=metadata.fish.F4.F8.unifrac$Species #Add Species to the new df.
levels(nmds.scores.2$Species)=c("Kyphosus hawaiiensis","Kyphosus vaigiensis","Kyphosus cinerascens","")
nmds.scores.2$Fish_Number=metadata.fish.F4.F8.unifrac$Fish_Number
nmds.scores.2$Gut_Section=metadata.fish.F4.F8.unifrac$Gut_Section
nmds.scores.2$Sample_Type=c(rep("16S",times=14),rep("Metabolomics",times=23))

#test different color gradients
colvectest=c("#7b5804",viridis(n=7,option="D",begin=.15))
colvectest2=c("mediumpurple1",viridis(n=7,option="D",begin=.15))
pointvectest=c(22,21,24)

ggplot(nmds.scores.2,aes(x=NMDS1,y=NMDS2,shape=Species,color=Gut_SubSection,size=2))+ #plot nmds.scores, specify x and y axis, shape corresponding to species, color corresponding to Ranked_Distance_from_Stomach
  scale_color_manual(values=colvec2)+ #set gradient scale to desired viridis scale
  geom_point(data=nmds.scores.2,aes(x=NMDS1,y=NMDS2))+ #add points
  theme_classic()+
  ggtitle("new gradient1")#export graphs as 8.5x11 PDF (landscape)

plot(nmds.weighted.unifrac.F4.F8,type="n") #Set up the NMDS plot environment
points(nmds.weighted.unifrac.F4.F8, #plot each sample as a point
       cex=2,
       pch=pointvec[metadata.fish.F4.F8.unifrac$Species], #point shape corresponds to fish Species
       col=colvectest2[as.factor(metadata.fish.F4.F8.unifrac$Gut_SubSection)]) 
text(-.5,.275,label="Stress=.14")
ordiellipse(nmds.weighted.unifrac.F4.F8,metadata.fish.F4.F8.unifrac$Gut_Zone,kind="sd",draw="polygon",alpha=c(0,0),lwd=2.5,border=c("mediumpurple1",viridis(n=3,option="D")[2:3])) #add SD ellipses, colored by Gut_Section. For now I will not use this.

plot(nmds.weighted.unifrac.F4.F8,type="n") #Set up the NMDS plot environment
points(nmds.weighted.unifrac.F4.F8, #plot each sample as a point
       cex=2,
       pch=pointvectest[metadata.fish.F4.F8.unifrac$Species], #point shape corresponds to fish Species
       bg=colvec2[as.factor(metadata.fish.F4.F8.unifrac$Gut_SubSection)],
       col=colvec[as.factor(metadata.fish.F4.F8.unifrac$Gut_Zone)],
       lwd=2.5)

ggplot(nmds.scores.2,aes(x=NMDS1,y=NMDS2,shape=Species,color=Gut_SubSection,size=2))+ #plot nmds.scores, specify x and y axis, shape corresponding to species, color corresponding to Ranked_Distance_from_Stomach
  scale_color_manual(values=colvec2)+ #set gradient scale to desired viridis scale
  geom_point(data=nmds.scores.2,aes(x=NMDS1,y=NMDS2))+ #add points
  theme_classic()+
  ggtitle("old gradient")

ggplot(alpha.diversity.mean.gut_subsection,aes(x=Gut_SubSection,y=sobs,fill=Gut_SubSection))+
  ylim(0,670)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=viridis(n=8))+
  geom_errorbar(aes(ymin=sobs-SE_sobs,ymax=sobs+SE_sobs),width=.5)+
  ggtitle("A")+
  theme_classic()+
  theme(legend.position = "none")+
  ggtitle("old gradient")+
  geom_text(x=1:8,y=c(655,275,240,275,440,625,595,570),aes(label=c("A","B","B","B","B/C","A/C","A/C","A/C")))

ggplot(alpha.diversity.mean.gut_subsection,aes(x=Gut_SubSection,y=sobs,fill=Gut_SubSection))+
  ylim(0,670)+
  geom_bar(stat="identity")+
  scale_fill_manual(values=colvectest)+
  geom_errorbar(aes(ymin=sobs-SE_sobs,ymax=sobs+SE_sobs),width=.5)+
  ggtitle("A")+
  theme_classic()+
  theme(legend.position = "none")+
  ggtitle("new gradient1")+
  geom_text(x=1:8,y=c(655,275,240,275,440,625,595,570),aes(label=c("A","B","B","B","B/C","A/C","A/C","A/C")))

plot(nmds.weighted.unifrac.F4.F8,type="n") #Set up the NMDS plot environment
points(nmds.weighted.unifrac.F4.F8, #plot each sample as a point
       cex=2,
       pch=pointvec[metadata.fish.F4.F8.unifrac$Species], #point shape corresponds to fish Species
       col=colvec2[as.factor(metadata.fish.F4.F8.unifrac$Gut_SubSection)]) 
text(-.5,.275,label="Stress=.14")
ordiellipse(nmds.weighted.unifrac.F4.F8,metadata.fish.F4.F8.unifrac$Gut_Zone,kind="sd",draw="polygon",alpha=c(0,0),lwd=2.5,border=c("mediumpurple1",viridis(n=3,option="D")[2:3])) #add SD ellipses, colored by Gut_Section.
ordiarrows(nmds.weighted.unifrac.F4.F8,group=metadata.fish.F4.F8.unifrac$Fish_Number,order.by=metadata.fish.F4.F8.unifrac$Gut_SubSection,lwd=2,col="black")



#this was testing color gradient, but now it is for generating the fill scale for the legen for the procrustes figure
ggplot(nmds.scores.2,aes(x=NMDS1,y=NMDS2))+ #plot nmds.scores, specify x and y axis, shape corresponding to species, color corresponding to Ranked_Distance_from_Stomach
  geom_point(data=nmds.scores.2,aes(x=NMDS1,y=NMDS2,shape=Sample_Type,size=2))+ #add points
  scale_shape_manual(values = c(21, 16))+
  theme_classic()+
  ggtitle("new gradient1")#export graphs as 8.5x11 PDF (landscape)

#graph the relative abundance of significant ASVs across the gut subsections as stacked barcharts. Graph within family,fill by family ID but keep individual ASVs seperated by black lines. Lable ASV# on each portion of the barchart, with the size of the label proportional to the relabund of that ASV in that sample.
#proteo.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Proteobacteria",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Proteobacteria",]$mean_relabund))*15,position = position_stack(vjust = 0.5))+
#ggtitle(label="Proteobacteria") 

#unclassified.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Bacteria_unclassified",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Bacteria_unclassified",]$mean_relabund))*30,position = position_stack(vjust = 0.5))+
#ggtitle(label="Unclassified_Phylum") #save as landscape PDF

#bacteroidetes.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Bacteroidetes",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
# geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Bacteroidetes",]$mean_relabund))*25,position = position_stack(vjust = 0.5))+
#ggtitle(label="Bacteroidetes") #save as landscape PDF

#tenericutes.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Tenericutes",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Tenericutes",]$mean_relabund))*25,position = position_stack(vjust = 0.5))+
# ggtitle(label="Tenericutes") #save as landscape PDF

#actino.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Actinobacteria",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Actinobacteria",]$mean_relabund))*80,position = position_stack(vjust = 0.5))+
#ggtitle(label="Actinobacteria") #save as landscape PDF

#firmicutes.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Firmicutes",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
# geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Firmicutes",]$mean_relabund))*15,position = position_stack(vjust = 0.5))+
# ggtitle(label="Firmicutes")

#cyano.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Cyanobacteria",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
# geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Cyanobacteria",]$mean_relabund))*75,position = position_stack(vjust = 0.5))+
# ggtitle(label="Cyanobacteria")

#lentis.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Lentisphaerae",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
# geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Lentisphaerae",]$mean_relabund))*50,position = position_stack(vjust = 0.5))+
#ggtitle(label="Lentisphaerae")

#parabasalia.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Parabasalia",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Parabasalia",]$mean_relabund))*40,position = position_stack(vjust = 0.5))+
#ggtitle(label="Parabasalia")

#spirochaetes.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Spirochaetes",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Spirochaetes",]$mean_relabund))*40,position = position_stack(vjust = 0.5))+
#ggtitle(label="Spirochaetes")

#synergist.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Synergistetes",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
# geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Synergistetes",]$mean_relabund))*200,position = position_stack(vjust = 0.5))+
#ggtitle(label="Synergistetes")

#verruco.sig.otus=ggplot(mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Verrucomicrobia",],aes(x=Gut_SubSection,y=mean_relabund,fill=Family))+
#geom_bar(stat="identity",col="black")+
#geom_text(aes(label=OTU),size=sqrt((mixed.mod.sig.otu.Gut_SubSection.longformat[mixed.mod.sig.otu.Gut_SubSection.longformat$Phylum=="Verrucomicrobia",]$mean_relabund))*30,position = position_stack(vjust = 0.5))+
# ggtitle(label="Verrucomicrobia")
