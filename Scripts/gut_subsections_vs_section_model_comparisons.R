#now run a mixed model on the transformed ASV data. Have relative abundance be a function of Species, Gut_SubSection, and Fish_Number be a random effect. Interaction term has been removed since it was not significant in the multivariate analysis. Could re-run model with Species also as a random effect.
mixed.mod=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Species+Gut_Section+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod1.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod) #apply mixed.mod1 only to columns containing ASV data. Save as new object.

#Now let's inspect the results of the mixed models.
mixed.mod1.singularity=lapply(mixed.mod1.relabund.t.fish.F4.F8.cull.trans,FUN=isSingular) #test each mixed model for singularity, save in list. After inspecting the list it looks like most models are not singular.
anova.kw=function(x) { #make a function that performs anova (Kenward-Roger)
  results=anova(x,ddf="Kenward-Roger")
  return(results)
}
mixed.mod1.results=lapply(mixed.mod1.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod1.pvals=as.data.frame(t(as.data.frame(mixed.mod1.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod1.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).

#Work up p values. #NEED TO REDO PVALUE ADJUSTMENT, COMBINING PVALS FOR GUT_SUBSECTION AND SPCIES INTO 1 VECTOR PRIOR TO ADJUSTMENT?
mixed.mod1.pvals$Gut_Section.padjust=p.adjust(mixed.mod1.pvals$Gut_Section,method="BH") #Adjust p values for multiple comparisons
mixed.mod1.pvals$Species.padjust=p.adjust(mixed.mod1.pvals$Species,method="BH") #Adjust p values for multiple comparisons
mixed.mod1.pvals$Gut_Section.padjust1=p.adjust(c(mixed.mod1.pvals$Gut_Section,mixed.mod1.pvals$Species),method="BH")[1:692] #Adjust p values for multiple comparisons
mixed.mod1.pvals$Species.padjust1=p.adjust(c(mixed.mod1.pvals$Gut_Section,mixed.mod1.pvals$Species),method="BH")[693:1384] #Adjust p values for multiple comparisons

rownames(mixed.mod1.pvals)=colnames(relabund.t.fish.F4.F8.cull.trans)[1:692] #Change rownames to correspond with tested OTU names.
mixed.mod1.pvals=cbind(mixed.mod1.pvals,taxonomy.fish.F4.F8.cull,relabund.fish.F4.F8.cull) #merge mixed.mod.sig.otu.Gut_Section with the associated OTU taxonomies.
mixed.mod.sig.otu.Gut_Section=cbind(mixed.mod.sig.otu.Gut_Section,taxonomy.fish.F4.F8.cull[rownames(taxonomy.fish.F4.F8.cull) %in% rownames(mixed.mod.sig.otu.Gut_Section),],relabund.fish.F4.F8.cull[rownames(relabund.fish.F4.F8.cull) %in% rownames(mixed.mod.sig.otu.Gut_Section),]) #merge mixed.mod.sig.otu.Gut_Section with the associated OTU taxonomies.
mixed.mod.sig.otu.Species=mixed.mod1.pvals[mixed.mod1.pvals$Species.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_Section pvals. Save as new df. It appears there are no significant Species pvals after padjustment. Will need to look into why this is. Perhaps if using the full sample set that includes F1-F3 GI/HG samples there would be enough replication to find significant species pvals.
#create an output zscore df and calculate z score for each otu using the scale function.
mixed.mod.sig.otu.Gut_Section.trans=cbind(mixed.mod.sig.otu.Gut_Section[,1:10],asin(sqrt(mixed.mod.sig.otu.Gut_Section[,-1:-10]))) #make a asinsqrt transfomred mixed mod sig otu gut subsection df. Normalizing ASVs is necessary for z scoring.
mixed.mod.sig.otu.Gut_Section_zscore=as.data.frame(t(mixed.mod.sig.otu.Gut_Section.trans[,-1:-10]))

for (i in 1:ncol(mixed.mod.sig.otu.Gut_Section_zscore)) { #calculate the zscore
  mixed.mod.sig.otu.Gut_Section_zscore[,i]=scale(mixed.mod.sig.otu.Gut_Section_zscore[,i],center=TRUE,scale=TRUE)
}

#now compare the OTUs that made it into the Section model but not the Subsection model
mixed.mod.sig.otu.Gut_Section_new_ASVs_zscore=mixed.mod.sig.otu.Gut_Section_zscore[,colnames(mixed.mod.sig.otu.Gut_Section_zscore) %in% colnames(mixed.mod.sig.otu.Gut_SubSection_zscore)=="FALSE"]
mixed.mod.sig.otu.Gut_Section_new_ASVS=mixed.mod.sig.otu.Gut_Section[rownames(mixed.mod.sig.otu.Gut_Section) %in% rownames(mixed.mod.sig.otu.Gut_SubSection)=="FALSE",]

write.csv(mixed.mod.sig.otu.Gut_Section_new_ASVs_zscore,file="mixed.mod.sig.otu.Gut_Section_new_ASVs_zscore.csv")
write.csv(mixed.mod.sig.otu.Gut_Section_new_ASVS,file="mixed.mod.sig.otu.Gut_Section_new_ASVS.csv")

#subset the relabund.t.fish.F4.F8.cull df to select only the significant ASVs (by gut subsection).
relabund.t.fish.F4.F8.cull.sig.section=relabund.t.fish.F4.F8.cull[,colnames(relabund.t.fish.F4.F8.cull) %in% rownames(mixed.mod.sig.otu.Gut_Section)]

cull.otu(relabund.t.fish.F4.F8.cull.sig.section,1,.01,.01)
relabund.t.fish.F4.F8.cull.sig.section.cull=relabund.df.cull
ncol(relabund.t.fish.F4.F8.cull.sig.section.cull) #59 significant ASVs remain after culling at the threshold: .01 relabund (1%) in at least 1 sample.
relabund.t.fish.F4.F8.cull.sig.section.cull=cbind(relabund.t.fish.F4.F8.cull.sig.section.cull,metadata.fish.F4.F8.noPC.noST.unifrac) #add in metadata

#subset the mixed.mod.sig.otu.Gut_SubSection df for only these 59 ASvs.
mixed.mod.sig.otu.Gut_Section.cull=mixed.mod.sig.otu.Gut_Section[rownames(mixed.mod.sig.otu.Gut_Section) %in% colnames(relabund.t.fish.F4.F8.cull.sig.section.cull),]

#Updated colnames
colnames(relabund.t.fish.F4.F8.cull.sig.section.cull)[1:90]=paste(mixed.mod.sig.otu.Gut_Section.cull$Genus,rownames(mixed.mod.sig.otu.Gut_Section.cull),sep="_") #change colnames to correspond to ASV Family_ASVnumber

#subset the mixed.mod.sig.otu.Gut_SubSection_zscore df to select only the culled significant ASVs.
mixed.mod.sig.otu.Gut_Section_zscore.cull=mixed.mod.sig.otu.Gut_Section_zscore[,colnames(mixed.mod.sig.otu.Gut_Section_zscore) %in% rownames(mixed.mod.sig.otu.Gut_Section.cull)]




#NOW RUN THE NESTED MODEL
#now run a mixed model on the transformed ASV data. Have relative abundance be a function of Species, Gut_SubSection, and Fish_Number be a random effect. Interaction term has been removed since it was not significant in the multivariate analysis. Could re-run model with Species also as a random effect.
mixed.mod.nest=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Gut_Section/Gut_SubSection+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod.nest.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod.nest) #apply mixed.mod1 only to columns containing ASV data. Save as new object.

#Now let's inspect the results of the mixed models.
mixed.mod.nest.singularity=lapply(mixed.mod.nest.relabund.t.fish.F4.F8.cull.trans,FUN=isSingular) #test each mixed model for singularity, save in list. After inspecting the list it looks like most models are not singular.
anova.kw=function(x) { #make a function that performs anova (Kenward-Roger)
  results=anova(x,ddf="Kenward-Roger")
  return(results)
}
mixed.mod.nest.results=lapply(mixed.mod.nest.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod.nest.pvals=as.data.frame(t(as.data.frame(mixed.mod.nest.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod.nest.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).

#Work up p values. #NEED TO REDO PVALUE ADJUSTMENT, COMBINING PVALS FOR GUT_SUBSECTION AND SPCIES INTO 1 VECTOR PRIOR TO ADJUSTMENT?
mixed.mod.nest.pvals$Gut_Section.padjust=p.adjust(mixed.mod.nest.pvals$Gut_Section,method="BH") #Adjust p values for multiple comparisons
mixed.mod.nest.pvals$Gut_Section.padjust1=p.adjust(c(mixed.mod.nest.pvals$Gut_Section,mixed.mod.nest.pvals$`Gut_Section:Gut_SubSection`),method="BH")[1:692] #Adjust p values for multiple comparisons
mixed.mod.nest.pvals$Gut_Section.Subsection.padjust=p.adjust(mixed.mod.nest.pvals$`Gut_Section:Gut_SubSection`,method="BH") #Adjust p values for multiple comparisons
mixed.mod.nest.pvals$Gut_Section.Subsection.padjust1=p.adjust(c(mixed.mod.nest.pvals$Gut_Section,mixed.mod.nest.pvals$`Gut_Section:Gut_SubSection`),method="BH")[693:1384]
rownames(mixed.mod.nest.pvals)=colnames(relabund.t.fish.F4.F8.cull.trans)[1:692] #Change rownames to correspond with tested OTU names.
mixed.mod.nest.sig.otu.Gut_Section=mixed.mod.nest.pvals[mixed.mod.nest.pvals$Gut_Section.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_Section pvals. Save as new df.
mixed.mod.nest.sig.otu.Gut_Section=cbind(mixed.mod.nest.sig.otu.Gut_Section,taxonomy.fish.F4.F8.cull[rownames(taxonomy.fish.F4.F8.cull) %in% rownames(mixed.mod.nest.sig.otu.Gut_Section),],relabund.fish.F4.F8.cull[rownames(relabund.fish.F4.F8.cull) %in% rownames(mixed.mod.nest.sig.otu.Gut_Section),]) #merge mixed.mod.nest.sig.otu.Gut_Section with the associated OTU taxonomies.
mixed.mod.nest.sig.otu.Species=mixed.mod.nest.pvals[mixed.mod.nest.pvals$Species.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_Section pvals. Save as new df. It appears there are no significant Species pvals after padjustment. Will need to look into why this is. Perhaps if using the full sample set that includes F1-F3 GI/HG samples there would be enough replication to find significant species pvals.
#create an output zscore df and calculate z score for each otu using the scale function.
mixed.mod.nest.sig.otu.Gut_Section.trans=cbind(mixed.mod.nest.sig.otu.Gut_Section[,1:12],asin(sqrt(mixed.mod.nest.sig.otu.Gut_Section[,-1:-12]))) #make a asinsqrt transfomred mixed mod sig otu gut subsection df. Normalizing ASVs is necessary for z scoring.
mixed.mod.nest.sig.otu.Gut_Section_zscore=as.data.frame(t(mixed.mod.nest.sig.otu.Gut_Section.trans[,-1:-12]))

for (i in 1:ncol(mixed.mod.nest.sig.otu.Gut_Section_zscore)) { #calculate the zscore
  mixed.mod.nest.sig.otu.Gut_Section_zscore[,i]=scale(mixed.mod.nest.sig.otu.Gut_Section_zscore[,i],center=TRUE,scale=TRUE)
}







#NOW RUN AN UPDATED SUBSECTION MODEL WITHOUT SPECIES
#now run a mixed model on the transformed ASV data. Have relative abundance be a function of Species, Gut_SubSection, and Fish_Number be a random effect. Interaction term has been removed since it was not significant in the multivariate analysis. Could re-run model with Species also as a random effect.
mixed.mod2=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Gut_SubSection+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod2.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod2) #apply mixed.mod2 only to columns containing ASV data. Save as new object.

#Now let's inspect the results of the mixed models.
mixed.mod2.singularity=lapply(mixed.mod2.relabund.t.fish.F4.F8.cull.trans,FUN=isSingular) #test each mixed model for singularity, save in list. After inspecting the list it looks like most models are not singular.
anova.kw=function(x) { #make a function that performs anova (Kenward-Roger)
  results=anova(x,ddf="Kenward-Roger")
  return(results)
}
mixed.mod2.results=lapply(mixed.mod2.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod2.pvals=as.data.frame(t(as.data.frame(mixed.mod2.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod2.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).

#Work up p values. #NEED TO REDO PVALUE ADJUSTMENT, COMBINING PVALS FOR GUT_SUBSECTION AND SPCIES INTO 1 VECTOR PRIOR TO ADJUSTMENT?
colnames(mixed.mod2.pvals)="Gut_SubSection"
mixed.mod2.pvals$Gut_SubSection.padjust=p.adjust(mixed.mod2.pvals$Gut_SubSection,method="BH") #Adjust p values for multiple comparisons
rownames(mixed.mod2.pvals)=colnames(relabund.t.fish.F4.F8.cull.trans)[1:692] #Change rownames to correspond with tested OTU names.
mixed.mod2.sig.otu.Gut_SubSection=mixed.mod2.pvals[mixed.mod2.pvals$Gut_SubSection.padjust<=.05,] #subset the mixed.mod2.pvals df for significant (p<=.05) Gut_SubSection pvals. Save as new df.
mixed.mod2.sig.otu.Gut_SubSection=cbind(mixed.mod2.sig.otu.Gut_SubSection,taxonomy.fish.F4.F8.cull[rownames(taxonomy.fish.F4.F8.cull) %in% rownames(mixed.mod2.sig.otu.Gut_SubSection),],relabund.fish.F4.F8.cull[rownames(relabund.fish.F4.F8.cull) %in% rownames(mixed.mod2.sig.otu.Gut_SubSection),]) #merge mixed.mod2.sig.otu.Gut_SubSection with the associated OTU taxonomies.
#create an output zscore df and calculate z score for each otu using the scale function.
mixed.mod2.sig.otu.Gut_SubSection.trans=cbind(mixed.mod2.sig.otu.Gut_SubSection[,1:8],asin(sqrt(mixed.mod2.sig.otu.Gut_SubSection[,-1:-8]))) #make a asinsqrt transfomred mixed mod sig otu gut subsection df. Normalizing ASVs is necessary for z scoring.
mixed.mod2.sig.otu.Gut_SubSection_zscore=as.data.frame(t(mixed.mod2.sig.otu.Gut_SubSection.trans[,-1:-8]))

for (i in 1:ncol(mixed.mod2.sig.otu.Gut_SubSection_zscore)) { #calculate the zscore
  mixed.mod2.sig.otu.Gut_SubSection_zscore[,i]=scale(mixed.mod2.sig.otu.Gut_SubSection_zscore[,i],center=TRUE,scale=TRUE)
}

#now compare the OTUs that made it into the Section model but not the Subsection model
mixed.mod2.sig.otu.Gut_SubSection_new_ASVs_zscore=mixed.mod2.sig.otu.Gut_SubSection_zscore[,colnames(mixed.mod2.sig.otu.Gut_SubSection_zscore) %in% colnames(mixed.mod.nest.sig.otu.Gut_Section_zscore)=="FALSE"]
mixed.mod2.sig.otu.Gut_SubSection_new_ASVS=mixed.mod2.sig.otu.Gut_SubSection[rownames(mixed.mod2.sig.otu.Gut_SubSection) %in% rownames(mixed.mod.nest.sig.otu.Gut_Section_zscore)=="FALSE",]

write.csv(mixed.mod2.sig.otu.Gut_SubSection_new_ASVs_zscore,file="mixed.mod2.sig.otu.Gut_SubSection_new_ASVs_zscore.csv")
write.csv(mixed.mod2.sig.otu.Gut_SubSection_new_ASVS,file="mixed.mod2.sig.otu.Gut_SubSection_new_ASVS.csv")







#NOW RUN a mode with just gut section
#now run a mixed model on the transformed ASV data. Have relative abundance be a function of Species, Gut_SubSection, and Fish_Number be a random effect. Interaction term has been removed since it was not significant in the multivariate analysis. Could re-run model with Species also as a random effect.
mixed.mod.section2=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Gut_Section+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod.section2.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod.section2) #apply mixed.mod1 only to columns containing ASV data. Save as new object.

#Now let's inspect the results of the mixed models.
mixed.mod.section2.singularity=lapply(mixed.mod.section2.relabund.t.fish.F4.F8.cull.trans,FUN=isSingular) #test each mixed model for singularity, save in list. After inspecting the list it looks like most models are not singular.
anova.kw=function(x) { #make a function that performs anova (Kenward-Roger)
  results=anova(x,ddf="Kenward-Roger")
  return(results)
}
mixed.mod.section2.results=lapply(mixed.mod.section2.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod.section2.pvals=as.data.frame(t(as.data.frame(mixed.mod.section2.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod.section2.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).
colnames(mixed.mod.section2.pvals)="Gut_Section"

#Work up p values. #NEED TO REDO PVALUE ADJUSTMENT, COMBINING PVALS FOR GUT_SUBSECTION AND SPCIES INTO 1 VECTOR PRIOR TO ADJUSTMENT?
mixed.mod.section2.pvals$Gut_Section.padjust1=p.adjust(mixed.mod.section2.pvals$Gut_Section,method="BH") #Adjust p values for multiple comparisons
rownames(mixed.mod.section2.pvals)=colnames(relabund.t.fish.F4.F8.cull.trans)[1:692] #Change rownames to correspond with tested OTU names.
mixed.mod.section2.sig.otu.Gut_Section=mixed.mod.section2.pvals[mixed.mod.section2.pvals$Gut_Section.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_Section pvals. Save as new df.
mixed.mod.section2.sig.otu.Gut_Section=cbind(mixed.mod.section2.sig.otu.Gut_Section,taxonomy.fish.F4.F8.cull[rownames(taxonomy.fish.F4.F8.cull) %in% rownames(mixed.mod.section2.sig.otu.Gut_Section),],relabund.fish.F4.F8.cull[rownames(relabund.fish.F4.F8.cull) %in% rownames(mixed.mod.section2.sig.otu.Gut_Section),]) #merge mixed.mod.section2.sig.otu.Gut_Section with the associated OTU taxonomies.
#create an output zscore df and calculate z score for each otu using the scale function.
mixed.mod.section2.sig.otu.Gut_Section.trans=cbind(mixed.mod.section2.sig.otu.Gut_Section[,1:12],asin(sqrt(mixed.mod.section2.sig.otu.Gut_Section[,-1:-12]))) #make a asinsqrt transfomred mixed mod sig otu gut subsection df. Normalizing ASVs is necessary for z scoring.
mixed.mod.section2.sig.otu.Gut_Section_zscore=as.data.frame(t(mixed.mod.section2.sig.otu.Gut_Section.trans[,-1:-12]))

for (i in 1:ncol(mixed.mod.section2.sig.otu.Gut_Section_zscore)) { #calculate the zscore
  mixed.mod.section2.sig.otu.Gut_Section_zscore[,i]=scale(mixed.mod.section2.sig.otu.Gut_Section_zscore[,i],center=TRUE,scale=TRUE)
}






#NOW RUN THE updated NESTED MODEL with species
#now run a mixed model on the transformed ASV data. Have relative abundance be a function of Species, Gut_SubSection, and Fish_Number be a random effect. Interaction term has been removed since it was not significant in the multivariate analysis. Could re-run model with Species also as a random effect.
mixed.mod.nest1=function(x) { #create a function that performs the desired mixed model
  mod=lmer(x~Gut_Section/Gut_SubSection+Species+(1|Fish_Number),data=relabund.t.fish.F4.F8.cull.trans)
  return(mod)
}
mixed.mod.nest1.relabund.t.fish.F4.F8.cull.trans=apply(relabund.t.fish.F4.F8.cull.trans[,-693:-702],2,FUN=mixed.mod.nest1) #apply mixed.mod1 only to columns containing ASV data. Save as new object.

#Now let's inspect the results of the mixed models.
mixed.mod.nest1.singularity=lapply(mixed.mod.nest1.relabund.t.fish.F4.F8.cull.trans,FUN=isSingular) #test each mixed model for singularity, save in list. After inspecting the list it looks like most models are not singular.
anova.kw=function(x) { #make a function that performs anova (Kenward-Roger)
  results=anova(x,ddf="Kenward-Roger")
  return(results)
}
mixed.mod.nest1.results=lapply(mixed.mod.nest1.relabund.t.fish.F4.F8.cull.trans,FUN=anova.kw) #perform an anova(ddf="Kenward-Roger") on all mixed models, save output.
mixed.mod.nest1.pvals=as.data.frame(t(as.data.frame(mixed.mod.nest1.results))[seq(from=6,to=4152,by=6),]) #extract the pvas from the mixed.mod.nest1.results list. Do this by converting list to df, transposing it, and then extracting every 6th column (which correspond to the pvals from each model).

#Work up p values. #NEED TO REDO PVALUE ADJUSTMENT, COMBINING PVALS FOR GUT_SUBSECTION AND SPCIES INTO 1 VECTOR PRIOR TO ADJUSTMENT?
mixed.mod.nest1.pvals$Gut_Section.padjust=p.adjust(mixed.mod.nest1.pvals$Gut_Section,method="BH") #Adjust p values for multiple comparisons
mixed.mod.nest1.pvals$Gut_Section.padjust1=p.adjust(c(mixed.mod.nest1.pvals$Gut_Section,mixed.mod.nest1.pvals$`Gut_Section:Gut_SubSection`),method="BH")[1:692] #Adjust p values for multiple comparisons
mixed.mod.nest1.pvals$Gut_Section.Subsection.padjust=p.adjust(mixed.mod.nest1.pvals$`Gut_Section:Gut_SubSection`,method="BH") #Adjust p values for multiple comparisons
mixed.mod.nest1.pvals$Gut_Section.Subsection.padjust1=p.adjust(c(mixed.mod.nest1.pvals$Gut_Section,mixed.mod.nest1.pvals$`Gut_Section:Gut_SubSection`),method="BH")[693:1384]
mixed.mod.nest1.pvals$Gut_Section.Subsection.padjust=p.adjust(mixed.mod.nest1.pvals$`Gut_Section:Gut_SubSection`,method="BH") #Adjust p values for multiple comparisons
mixed.mod.nest1.pvals$Gut_Section.Subsection.padjust1=p.adjust(c(mixed.mod.nest1.pvals$Gut_Section,mixed.mod.nest1.pvals$`Gut_Section:Gut_SubSection`),method="BH")[693:1384]
rownames(mixed.mod.nest1.pvals)=colnames(relabund.t.fish.F4.F8.cull.trans)[1:692] #Change rownames to correspond with tested OTU names.
mixed.mod.nest1.pvals=cbind(mixed.mod.nest1.pvals,taxonomy.fish.F4.F8.cull)
mixed.mod.nest1.sig.otu.Gut_Section=mixed.mod.nest1.pvals[mixed.mod.nest1.pvals$Gut_Section.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_Section pvals. Save as new df.
mixed.mod.nest1.sig.otu.Gut_Section=cbind(mixed.mod.nest1.sig.otu.Gut_Section,taxonomy.fish.F4.F8.cull[rownames(taxonomy.fish.F4.F8.cull) %in% rownames(mixed.mod.nest1.sig.otu.Gut_Section),],relabund.fish.F4.F8.cull[rownames(relabund.fish.F4.F8.cull) %in% rownames(mixed.mod.nest1.sig.otu.Gut_Section),]) #merge mixed.mod.nest1.sig.otu.Gut_Section with the associated OTU taxonomies.
mixed.mod.nest1.sig.otu.Species=mixed.mod.nest1.pvals[mixed.mod.nest1.pvals$Species.padjust<=.05,] #subset the mixed.mod.pvals df for significant (p<=.05) Gut_Section pvals. Save as new df. It appears there are no significant Species pvals after padjustment. Will need to look into why this is. Perhaps if using the full sample set that includes F1-F3 GI/HG samples there would be enough replication to find significant species pvals.
#create an output zscore df and calculate z score for each otu using the scale function.
mixed.mod.nest1.sig.otu.Gut_Section.trans=cbind(mixed.mod.nest1.sig.otu.Gut_Section[,1:12],asin(sqrt(mixed.mod.nest1.sig.otu.Gut_Section[,-1:-12]))) #make a asinsqrt transfomred mixed mod sig otu gut subsection df. Normalizing ASVs is necessary for z scoring.
mixed.mod.nest1.sig.otu.Gut_Section_zscore=as.data.frame(t(mixed.mod.nest1.sig.otu.Gut_Section.trans[,-1:-12]))

for (i in 1:ncol(mixed.mod.nest1.sig.otu.Gut_Section_zscore)) { #calculate the zscore
  mixed.mod.nest1.sig.otu.Gut_Section_zscore[,i]=scale(mixed.mod.nest1.sig.otu.Gut_Section_zscore[,i],center=TRUE,scale=TRUE)
}







