#### set-up#####

library(pavo)
library(onewaytests)
library(mclust)
library(RColorBrewer)
library(ComplexHeatmap)
library(diptest)

options(scipen = 999)
rawdat <- readRDS("/Users/rs155/Dropbox/Attila/rawdat.rds")

# rename specimens that now has lsumz number
colnames(rawdat) <- gsub(x = colnames(rawdat), pattern = "000002", replacement = "191084")
colnames(rawdat) <- gsub(x = colnames(rawdat), pattern = "000001", replacement = "191085")

# filter to cis-andean ssp
rawdat <- subset(x = rawdat, subset = c("_uro_", "_spa_"))
rawdat <- procspec(rspecdata = rawdat, opt = "smooth", fixneg = "zero", span = 0.05)

# exclude specimens w unknown sex
rawdat <- rawdat[, grep(x = colnames(rawdat), pattern = "_u_", invert = T)]

# exclude juves
rawdat <- subset(x = rawdat, subset = c("LSU175476", "LSU178366", "LSU178367", "LSU052174"), invert = T)

# rename to plumage patches
names(rawdat) <- gsub(x = names(rawdat), pattern = "_01$", replacement = "_crow")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_02$", replacement = "_crow")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_03$", replacement = "_crow")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_04$", replacement = "_back")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_05$", replacement = "_back")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_06$", replacement = "_back")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_07$", replacement = "_rump")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_08$", replacement = "_rump")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_09$", replacement = "_rump")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_10$", replacement = "_tail")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_11$", replacement = "_tail")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_12$", replacement = "_tail")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_13$", replacement = "_bell")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_14$", replacement = "_bell")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_15$", replacement = "_bell")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_16$", replacement = "_thro")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_17$", replacement = "_thro")
names(rawdat) <- gsub(x = names(rawdat), pattern = "_18$", replacement = "_thro")

indspecs <- aggspec(rawdat, by = names(rawdat), trim = FALSE)

# separate by patch
crow <- subset(x = indspecs, subset = "crow")
back <- subset(x = indspecs, subset = "back")
rump <- subset(x = indspecs, subset = "rump")
tail <- subset(x = indspecs, subset = "tail")
bell <- subset(x = indspecs, subset = "bell")
thro <- subset(x = indspecs, subset = "thro")

######## PCAs###########

# center and PCA each patch
pcawithsex <- T

if (pcawithsex == T) {
  crow <- procspec(crow, opt = "center")
  crow$wl <- NULL
  crow[nrow(crow) + 1, ] <- substr(x = colnames(crow), start = 17, stop = 17) # add a row for sex
  crow[nrow(crow), ] <- gsub(x = crow[nrow(crow), ], pattern = "f", replacement = 0) # convert sex to numeric
  crow[nrow(crow), ] <- gsub(x = crow[nrow(crow), ], pattern = "m", replacement = 1)
  rownames(crow)[nrow(crow)] <- "sex"
  crow <- as.data.frame(sapply(crow, as.numeric))
  crow <- crow[, colSums(is.na(crow)) == 0] # remove nas
  crow.pca <- prcomp(t(crow), scale = TRUE)

  back <- procspec(back, opt = "center")
  back$wl <- NULL
  back[nrow(back) + 1, ] <- substr(x = colnames(back), start = 17, stop = 17) # add a row for sex
  back[nrow(back), ] <- gsub(x = back[nrow(back), ], pattern = "f", replacement = 0) # convert sex to numeric
  back[nrow(back), ] <- gsub(x = back[nrow(back), ], pattern = "m", replacement = 1)
  rownames(back)[nrow(back)] <- "sex"
  back <- as.data.frame(sapply(back, as.numeric))
  back <- back[, colSums(is.na(back)) == 0] # remove nas
  back.pca <- prcomp(t(back), scale = TRUE)


  rump <- procspec(rump, opt = "center")
  rump$wl <- NULL
  rump[nrow(rump) + 1, ] <- substr(x = colnames(rump), start = 17, stop = 17) # add a row for sex
  rump[nrow(rump), ] <- gsub(x = rump[nrow(rump), ], pattern = "f", replacement = 0) # convert sex to numeric
  rump[nrow(rump), ] <- gsub(x = rump[nrow(rump), ], pattern = "m", replacement = 1)
  rownames(rump)[nrow(rump)] <- "sex"
  rump <- as.data.frame(sapply(rump, as.numeric))
  rump <- rump[, colSums(is.na(rump)) == 0] # remove nas
  rump.pca <- prcomp(t(rump), scale = TRUE)


  tail <- procspec(tail, opt = "center")
  tail$wl <- NULL
  tail[nrow(tail) + 1, ] <- substr(x = colnames(tail), start = 17, stop = 17) # add a row for sex
  tail[nrow(tail), ] <- gsub(x = tail[nrow(tail), ], pattern = "f", replacement = 0) # convert sex to numeric
  tail[nrow(tail), ] <- gsub(x = tail[nrow(tail), ], pattern = "m", replacement = 1)
  rownames(tail)[nrow(tail)] <- "sex"
  tail <- as.data.frame(sapply(tail, as.numeric))
  tail <- tail[, colSums(is.na(tail)) == 0] # remove nas
  tail.pca <- prcomp(t(tail), scale = TRUE)


  bell <- procspec(bell, opt = "center")
  bell$wl <- NULL
  bell[nrow(bell) + 1, ] <- substr(x = colnames(bell), start = 17, stop = 17) # add a row for sex
  bell[nrow(bell), ] <- gsub(x = bell[nrow(bell), ], pattern = "f", replacement = 0) # convert sex to numeric
  bell[nrow(bell), ] <- gsub(x = bell[nrow(bell), ], pattern = "m", replacement = 1)
  rownames(bell)[nrow(bell)] <- "sex"
  bell <- as.data.frame(sapply(bell, as.numeric))
  bell <- bell[, colSums(is.na(bell)) == 0] # remove nas
  bell.pca <- prcomp(t(bell), scale = TRUE)

  thro <- procspec(thro, opt = "center")
  thro$wl <- NULL
  thro[nrow(thro) + 1, ] <- substr(x = colnames(thro), start = 17, stop = 17) # add a row for sex
  thro[nrow(thro), ] <- gsub(x = thro[nrow(thro), ], pattern = "f", replacement = 0) # convert sex to numeric
  thro[nrow(thro), ] <- gsub(x = thro[nrow(thro), ], pattern = "m", replacement = 1)
  rownames(thro)[nrow(thro)] <- "sex"
  thro <- as.data.frame(sapply(thro, as.numeric))
  thro <- thro[, colSums(is.na(thro)) == 0] # remove nas
  thro.pca <- prcomp(t(thro), scale = TRUE)
}

######### FIGURE 2###########
pdf(file = "/Users/rs155/Dropbox/Attila/2022/FIG2.pdf", height = 6, width = 10)

layout(mat = matrix(data = 1:6, nrow = 2, ncol = 3))

rownames(crow.pca$rotation) <- c(rawdat$wl, "sex")
plot(y = crow.pca$rotation[, 1], x = rownames(crow.pca$rotation), ylab = "PC1 loading", xlab = "Wavelength (nm)", main = "Crown")
points(x = 702, y = crow.pca$rotation[402, 1], pch = 8)
pvar <- round(summary(crow.pca)$importance[2, 1], 2)
mtext(text = "Sex", side = 1, line = 3, at = 702, cex = 0.65)
text(labels = paste0("Proportion of variance explained=\n", pvar), x = 450, y = 0.05, cex = 0.85)

rownames(back.pca$rotation) <- c(rawdat$wl, "sex")
plot(y = back.pca$rotation[, 1], x = rownames(back.pca$rotation), ylab = "PC1 loading", xlab = "Wavelength (nm)", main = "Back")
points(x = 702, y = back.pca$rotation[402, 1], pch = 8)
pvar <- round(summary(back.pca)$importance[2, 1], 2)
mtext(text = "Sex", side = 1, line = 3, at = 702, cex = 0.65)
text(labels = paste0("Proportion of variance explained=\n", pvar), x = 450, y = 0.05, cex = 0.85)

rownames(rump.pca$rotation) <- c(rawdat$wl, "sex")
plot(y = rump.pca$rotation[, 1], x = rownames(rump.pca$rotation), ylab = "PC1 loading", xlab = "Wavelength (nm)", main = "Rump")
points(x = 702, y = rump.pca$rotation[402, 1], pch = 8)
pvar <- round(summary(rump.pca)$importance[2, 1], 2)
mtext(text = "Sex", side = 1, line = 3, at = 702, cex = 0.65)
text(labels = paste0("Proportion of \nvariance explained=\n", pvar), x = 600, y = 0.04, cex = 0.85)

rownames(tail.pca$rotation) <- c(rawdat$wl, "sex")
plot(y = tail.pca$rotation[, 1], x = rownames(tail.pca$rotation), ylab = "PC1 loading", xlab = "Wavelength (nm)", main = "Tail")
points(x = 702, y = tail.pca$rotation[402, 1], pch = 8)
pvar <- round(summary(tail.pca)$importance[2, 1], 2)
mtext(text = "Sex", side = 1, line = 3, at = 702, cex = 0.65)
text(labels = paste0("Proportion of variance explained=\n", pvar), x = 450, y = 0.05, cex = 0.85)

rownames(bell.pca$rotation) <- c(rawdat$wl, "sex")
plot(y = bell.pca$rotation[, 1], x = rownames(bell.pca$rotation), ylab = "PC1 loading", xlab = "Wavelength (nm)", main = "Belly")
points(x = 702, y = bell.pca$rotation[402, 1], pch = 8)
pvar <- round(summary(bell.pca)$importance[2, 1], 2)
mtext(text = "Sex", side = 1, line = 3, at = 702, cex = 0.65)
text(labels = paste0("Proportion of variance explained=\n", pvar), x = 420, y = 0.05, cex = 0.85)

rownames(thro.pca$rotation) <- c(rawdat$wl, "sex")
plot(y = thro.pca$rotation[, 1], x = rownames(thro.pca$rotation), ylab = "PC1 loading", xlab = "Wavelength (nm)", main = "Throat")
points(x = 702, y = thro.pca$rotation[402, 1], pch = 8)
pvar <- round(summary(thro.pca)$importance[2, 1], 2)
mtext(text = "Sex", side = 1, line = 3, at = 702, cex = 0.65)
text(labels = paste0("Proportion of variance explained=\n", pvar), x = 420, y = 0.05, cex = 0.85)
dev.off()

####################### PCAs cont'd########

if (pcawithsex == F) {
  crow <- procspec(crow, opt = "center")
  crow$wl <- NULL
  crow.pca <- prcomp(t(crow), scale = TRUE)

  back <- procspec(back, opt = "center")
  back$wl <- NULL
  back.pca <- prcomp(t(back), scale = TRUE)


  rump <- procspec(rump, opt = "center")
  rump$wl <- NULL
  rump.pca <- prcomp(t(rump), scale = TRUE)


  tail <- procspec(tail, opt = "center")
  tail$wl <- NULL
  tail.pca <- prcomp(t(tail), scale = TRUE)


  bell <- procspec(bell, opt = "center")
  bell$wl <- NULL
  bell.pca <- prcomp(t(bell), scale = TRUE)

  thro <- procspec(thro, opt = "center")
  thro$wl <- NULL
  thro.pca <- prcomp(t(thro), scale = TRUE)
}

####### McLust##########

# put first 5 PCs of each patch in a dataframe
crow5 <- as.data.frame(crow.pca$x[, 1:5])
row.names(crow5) <- substr(x = row.names(crow5), start = 1, stop = 27)
colnames(crow5) <- paste0(colnames(crow5), ".crow")

back5 <- as.data.frame(back.pca$x[, 1:5])
row.names(back5) <- substr(x = row.names(back5), start = 1, stop = 27)
colnames(back5) <- paste0(colnames(back5), ".back")

rump5 <- as.data.frame(rump.pca$x[, 1:5])
row.names(rump5) <- substr(x = row.names(rump5), start = 1, stop = 27)
colnames(rump5) <- paste0(colnames(rump5), ".rump")

tail5 <- as.data.frame(tail.pca$x[, 1:5])
row.names(tail5) <- substr(x = row.names(tail5), start = 1, stop = 27)
colnames(tail5) <- paste0(colnames(tail5), ".tail")

bell5 <- as.data.frame(bell.pca$x[, 1:5])
row.names(bell5) <- substr(x = row.names(bell5), start = 1, stop = 27)
colnames(bell5) <- paste0(colnames(bell5), ".bell")

thro5 <- as.data.frame(thro.pca$x[, 1:5])
row.names(thro5) <- substr(x = row.names(thro5), start = 1, stop = 27)
colnames(thro5) <- paste0(colnames(thro5), ".thro")

mergedpcs <- merge(x = crow5, y = back5, by = 0)
mergedpcs <- merge(x = mergedpcs, y = rump5, by.x = "Row.names", by.y = 0)
mergedpcs <- merge(x = mergedpcs, y = tail5, by.x = "Row.names", by.y = 0)
mergedpcs <- merge(x = mergedpcs, y = bell5, by.x = "Row.names", by.y = 0)
mergedpcs <- merge(x = mergedpcs, y = thro5, by.x = "Row.names", by.y = 0)

row.names(mergedpcs) <- mergedpcs$Row.names
mergedpcs$Row.names <- NULL

fitted_bic <- mclustBIC(mergedpcs)

######### FIGURE 3###########
pdf(file = "/Users/rs155/Dropbox/Attila/2022/FIG3.pdf", height = 6, width = 10)
plot(x = fitted_bic, legendArgs = list(bg = "white", bty = "o"))
dev.off()

clusters <- Mclust(data = mergedpcs, modelNames = "EVI") # fit best model

# histogram of assignment probs
# hist(clusters$z[,1],breaks=20,main='p(green/grey cluster)')

# import specimen metadata, including visual classification of morphs
metadat <- read.csv("/Users/rs155/Dropbox/Attila/2022/ATTILA MASTER SHEET _ECC.csv")
metadat$eamon.visual <- gsub(x = metadat$eamon.visual, pattern = "yellow-brown", replacement = "yellow brown")

# classif=as.data.frame(clusters$classification) #classification from mclust
# colnames(classif)='cluster'

# metadat=merge(x=metadat,y=classif,by.x='NAME',by.y=0)

# contingency table of visual morphs by clusters
# conti.table=table(metadat[,c(4,30)])


### chroma metrics#######

summ.crow <- summary(object = subset(x = indspecs, subset = "crow"), subset = c("S1U", "S1B", "S1G", "S1R"))
row.names(summ.crow) <- gsub(x = row.names(summ.crow), pattern = "_Reflection_crow", replacement = "")
colnames(summ.crow) <- paste0(colnames(summ.crow), ".crow")

summ.back <- summary(object = subset(x = indspecs, subset = "back"), subset = c("S1U", "S1B", "S1G", "S1R"))
row.names(summ.back) <- gsub(x = row.names(summ.back), pattern = "_Reflection_back", replacement = "")
colnames(summ.back) <- paste0(colnames(summ.back), ".back")

summ.rump <- summary(object = subset(x = indspecs, subset = "rump"), subset = c("S1U", "S1B", "S1G", "S1R"))
row.names(summ.rump) <- gsub(x = row.names(summ.rump), pattern = "_Reflection_rump", replacement = "")
colnames(summ.rump) <- paste0(colnames(summ.rump), ".rump")

summ.tail <- summary(object = subset(x = indspecs, subset = "tail"), subset = c("S1U", "S1B", "S1G", "S1R"))
row.names(summ.tail) <- gsub(x = row.names(summ.tail), pattern = "_Reflection_tail", replacement = "")
colnames(summ.tail) <- paste0(colnames(summ.tail), ".tail")

summ.bell <- summary(object = subset(x = indspecs, subset = "bell"), subset = c("S1U", "S1B", "S1G", "S1R"))
row.names(summ.bell) <- gsub(x = row.names(summ.bell), pattern = "_Reflection_bell", replacement = "")
colnames(summ.bell) <- paste0(colnames(summ.bell), ".bell")

summ.thro <- summary(object = subset(x = indspecs, subset = "thro"), subset = c("S1U", "S1B", "S1G", "S1R"))
row.names(summ.thro) <- gsub(x = row.names(summ.thro), pattern = "_Reflection_thro", replacement = "")
colnames(summ.thro) <- paste0(colnames(summ.thro), ".thro")

# merge all the above to metadat

metadat <- merge(x = metadat, y = summ.crow, by.x = "NAME", by.y = 0)
metadat <- merge(x = metadat, y = summ.back, by.x = "NAME", by.y = 0)
metadat <- merge(x = metadat, y = summ.rump, by.x = "NAME", by.y = 0)
metadat <- merge(x = metadat, y = summ.tail, by.x = "NAME", by.y = 0)
metadat <- merge(x = metadat, y = summ.bell, by.x = "NAME", by.y = 0)
metadat <- merge(x = metadat, y = summ.thro, by.x = "NAME", by.y = 0)

#### dip test#####
dip <- vector()
for (i in colnames(metadat)[31:54]) {
  p <- dip.test(metadat[, i])$p.value
  dip[i] <- p
}

#### visual models and heatmaps####
crow <- subset(x = indspecs, subset = "thro")

# add cluster to colnames
for (n in 2:length(colnames(crow))) {
  name <- colnames(crow)[n]
  name <- substr(name, 1, 27)
  cluster <- metadat[metadat$NAME == name, "eamon.visual"]
  colnames(crow)[n] <- substr(name, 1, 27)
  colnames(crow)[n] <- paste0(colnames(crow)[n], "_", cluster)
}

crow.model <- vismodel(rspecdata = crow, visual = "avg.v")
crow.dist <- coldist(modeldata = crow.model)

crow.dist$patch1 <- substr(crow.dist$patch1, 1, 27)
crow.dist$patch2 <- substr(crow.dist$patch2, 1, 27)


# create crowmatrix and populate it
crow.model$shortnames <- rownames(crow.model)
rownames(crow.model) <- substr(rownames(crow.model), 1, 27)
crow.model <- merge(crow.model, metadat, by.x = 0, by.y = 1)
crow.model$sex <- substr(rownames(crow.model), 17, 17)
crow.model$Row.names <- paste0(substr(crow.model$Row.names, 1, 15), substr(crow.model$Row.names, 18, 27))
crow.model$Row.names <- paste(crow.model$Row.names, crow.model$Sex, crow.model$eamon.visual, sep = "_")
crow.model <- crow.model[order(crow.model$Sex), ]
crow.model <- crow.model[order(crow.model$eamon.visual), ]

# adjust order of morphs
order <- c("", "gray", "green", "yellow brown", "brown", "rufous")

crow.model <- crow.model[order(match(table = order, x = crow.model$eamon.visual)), ]
crow.matrix <- as.data.frame(matrix(ncol = nrow(crow.model), nrow = nrow(crow.model)))

colnames(crow.matrix) <- crow.model$shortnames
rownames(crow.matrix) <- crow.model$shortnames

for (c in 1:nrow(crow.dist)) {
  sp1 <- crow.dist[c, 1]
  sp2 <- crow.dist[c, 2]
  sp1 <- grep(pattern = sp1, x = colnames(crow.matrix))
  sp2 <- grep(pattern = sp2, x = colnames(crow.matrix))
  jnd <- crow.dist[c, 3]
  crow.matrix[sp1, sp2] <- jnd
  crow.matrix[sp2, sp1] <- jnd
}


crow.matrix <- crow.matrix[nchar(rownames(crow.matrix)) > 29, nchar(colnames(crow.matrix)) > 29] # exclude ones with ? morph

colnames(crow.matrix) <- paste0(substr(x = colnames(crow.matrix), start = 17, stop = 18), substr(x = colnames(crow.matrix), start = 29, stop = 50))
colnames(crow.matrix) <- gsub(x = colnames(crow.matrix), pattern = "f_", replacement = "F_")
colnames(crow.matrix) <- gsub(x = colnames(crow.matrix), pattern = "m_", replacement = "M_")
colnames(crow.matrix) <- paste(substr(x = rownames(crow.matrix), start = 22, stop = 27), colnames(crow.matrix), sep = "_")
rownames(crow.matrix) <- colnames(crow.matrix)

pdf(file = "/Users/rs155/Dropbox/Attila/2022/HM_thro_w_labels.pdf", width = 10, family = "sans")


hm <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(crow.matrix),
  col = brewer.pal(7, "YlOrRd"),
  cluster_columns = F,
  cluster_rows = F,
  name = "Color distance (JNDs)",
  show_row_names = T,
  show_column_names = T
)

draw(hm)

dev.off()


##### JOINT PCA#######

# process patch dataframe
for (patch in c("crow", "back", "rump", "tail", "bell", "thro")) {
  p <- get(patch)
  p <- subset(x = indspecs, subset = patch)
  rownames(p) <- paste0(rawdat$wl, patch)
  p$wl <- NULL
  p <- as.data.frame(t(p))
  rownames(p) <- substr(x = rownames(p), start = 1, stop = 27)
  assign(patch, p)
}

# merge patches
uber.pca.dat <- merge(crow, back, by = 0)
uber.pca.dat <- merge(uber.pca.dat, rump, by.x = 1, by.y = 0)
uber.pca.dat <- merge(uber.pca.dat, tail, by.x = 1, by.y = 0)
uber.pca.dat <- merge(uber.pca.dat, bell, by.x = 1, by.y = 0)
uber.pca.dat <- merge(uber.pca.dat, thro, by.x = 1, by.y = 0)
rownames(uber.pca.dat) <- uber.pca.dat$Row.names
uber.pca.dat$Row.names <- NULL

# run PCA
uber.pca <- prcomp(uber.pca.dat, scale = TRUE)
summ.uber.pca <- summary(uber.pca)

# plot
dat <- merge(uber.pca$x, metadat, by.x = 0, by.y = 1)
dat$color <- dat$eamon.visual
dat$color <- gsub(x = dat$color, pattern = "green", replacement = "chartreuse4")
dat$color <- gsub(x = dat$color, pattern = "rufous", replacement = "chocolate3")
dat$color <- gsub(x = dat$color, pattern = "yellow brown", replacement = "darkolivegreen2")
dat$color <- gsub(x = dat$color, pattern = "brown", replacement = "burlywood4")
dat <- dat[dat$color != "", ]

### FIG 4 #####
pdf("/Users/rs155/Dropbox/Attila/2022/uber_pca_fig.pdf")
xlab <- "PC1 (33.0% of variance)"
ylab <- "PC2 (19.2% of variance)"
plot(dat$PC1, dat$PC2, col = dat$color, pch = 16, xlab = xlab, ylab = ylab, xlim = c(-100, 100), cex = 2.5)
legend(pch = 16, x = "topright", legend = c("green", "rufous", "yellow brown", "brown", "grey"), col = c("chartreuse4", "chocolate3", "darkolivegreen2", "burlywood4", "grey"), title = "Morphs", cex = 0.7)
dev.off()
