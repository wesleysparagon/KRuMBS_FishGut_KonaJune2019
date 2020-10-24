#read in dfs
cluster.dat=read.csv("mixed.mod.sig.otu.Gut_SubSection_zscore.cull.v1.csv")
relabund=read.csv(file="relabund.t.fish.F4.F8.cull.sig.subsection.cull.csv")
pvals.dat=read.csv(file="mixed.mod.sig.otu.Gut_SubSection.cull.csv")

#merge cluster.dat and pvals.dat
merged=merge(cluster.dat,pvals.dat,by.x="Column.1",by.y="X")

#remove unecessaey columns
merged=merged[,c(1,35:46)]
merged$Genus_ASVNumber=paste(merged$Genus.x,merged$Column.1,sep="_")

#avergae relabund data by gut subsection
relabund2=aggregate(relabund,by=list(relabund$Gut_SubSection),FUN=mean)
rownames(relabund2)=relabund2$Group.1
relabund2=as.data.frame(t(relabund2))
relabund2$ASV=rownames(relabund2)

merged2=merge(relabund2,merged,by.x="ASV",by.y="Genus_ASVNumber")

#write.csv(merged2,file="Table X4.csv")
