library(bestNormalize)
library(BoutrosLab.plotting.general)
library(tidyverse)

COMBINED_DATA <- read.csv('/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/CLINGEN_HG38_2021/CLINGEN_LOSS_DATA/COMBINED_DATA_LOSS.txt',sep="\t",header=FALSE)
colnames(COMBINED_DATA) <- c("PathogenicList", "PARAMETERS", "X")
ordernorm_obj <- orderNorm(COMBINED_DATA$X)
COMBINED_DATA$Normalized <- ordernorm_obj$x.t
COMBINED_DATA$EVENTS <- vapply(strsplit(as.character(COMBINED_DATA$PARAMETERS),"_"), `[`, 1, FUN.VALUE=character(1))
COMBINED_DATA$EVENTS <- as.numeric(as.character(COMBINED_DATA$EVENTS))
COMBINED_DATA$CLASS <- gsub('[[:digit:]]+_', '', COMBINED_DATA$PARAMETERS)
COMBINED_DATA$CLASS <- factor(COMBINED_DATA$CLASS, levels = c("B","USLB","VUS","B_USLB","USLB_VUS","B_USLB_VUS"))
COMBINED_DATA <- with(COMBINED_DATA, COMBINED_DATA[order(EVENTS, CLASS),])
supp.labs <- character()
for( i in unique(COMBINED_DATA$EVENTS)) {
supp.labs <- c(supp.labs, paste0("EVENTS=",i))}
names(supp.labs) <- unique(COMBINED_DATA$EVENTS)
cnts <- as.data.frame(tail(count(COMBINED_DATA, PathogenicList, EVENTS)[,3],length(unique(COMBINED_DATA$EVENTS))))[,1]

COMBINED_DATA_AVG <- aggregate(Normalized~PARAMETERS, COMBINED_DATA, mean)
COMBINED_DATA_AVG$EVENTS <- vapply(strsplit(as.character(COMBINED_DATA_AVG$PARAMETERS),"_"), `[`, 1, FUN.VALUE=character(1))
COMBINED_DATA_AVG$EVENTS <- as.numeric(as.character(COMBINED_DATA_AVG$EVENTS))
COMBINED_DATA_AVG$CLASS <- gsub('[[:digit:]]+_', '', COMBINED_DATA_AVG$PARAMETERS)
COMBINED_DATA_AVG$CLASS <- factor(COMBINED_DATA_AVG$CLASS, levels = c("B","USLB","VUS","B_USLB","USLB_VUS","B_USLB_VUS"))
COMBINED_DATA_AVG <- with(COMBINED_DATA_AVG, COMBINED_DATA_AVG[order(EVENTS, CLASS),])
COMBINED_DATA_AVG$COLOUR <- "NA"
COMBINED_DATA_AVG$COLOUR[COMBINED_DATA_AVG$Normalized < 0] <- "firebrick1"
COMBINED_DATA_AVG$COLOUR[COMBINED_DATA_AVG$Normalized >= 0] <- "dodgerblue1"

COMBINED_DATA_PathList <- aggregate(Normalized~PathogenicList, COMBINED_DATA, mean)
COMBINED_DATA_PathList$COLOUR <- "NA"
COMBINED_DATA_PathList$COLOUR[COMBINED_DATA_PathList$Normalized < 0] <- "firebrick1"
COMBINED_DATA_PathList$COLOUR[COMBINED_DATA_PathList$Normalized >= 0] <- "dodgerblue1"

MATRIX <- as.data.frame(pivot_wider(COMBINED_DATA[,c(1,2,4)], names_from=PathogenicList, values_from=Normalized))
rownames(MATRIX) <- MATRIX[,1]
MATRIX <- MATRIX[,-1]
write.csv(MATRIX,'/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/CLINGEN_HG38_2021/CLINGEN_LOSS_DATA/MATRIX_LOSS.csv')
MATRIX2 <- MATRIX
MATRIX2$AVERAGE <- rowMeans(MATRIX2)
BestParam <- cbind(as.data.frame(apply(MATRIX2, 2, which.max)), as.data.frame(apply(MATRIX2,2,max)))
colnames(BestParam) <- c("ROW","SCORE")
BestParam$PARAMETER <- rownames(MATRIX2)[BestParam$ROW]
BestParam$ROW <- NULL
write.csv(BestParam,'/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/CLINGEN_HG38_2021/CLINGEN_LOSS_DATA/LIST_SCORES_LOSS.csv')

key_data <- data.frame(x <- seq(-50,50,1))
key <- create.heatmap(
x = key_data,
clustering.method = 'none',
scale.data = FALSE,
colour.scheme = c("firebrick1", "white", "dodgerblue1"),
print.colour.key = FALSE,
yaxis.tck = 0,
xat = c(1,101), xaxis.lab = c('Negative', 'Positive'), xaxis.rot = 0, xaxis.cex = 1,
same.as.matrix = TRUE
)
HeatMap <- create.heatmap(
x = MATRIX, clustering.method = 'none',
colour.scheme = c('firebrick1','white','dodgerblue1'), print.colour.key = FALSE,
xaxis.lab = rownames(MATRIX), xaxis.cex = 0.75, xaxis.rot = 90, axis.xlab.padding = 1, xaxis.tck = 0.5,
yaxis.lab = colnames(MATRIX), yaxis.cex = 0.75, yaxis.tck = 0.5,
grid.row = TRUE, row.colour = 'black', grid.col=TRUE, force.grid.col = TRUE,
col.lines = c(cumsum(cnts)+0.5),
top.padding = 2, bottom.padding = 2, left.padding = 2, right.padding = 2
)
PathList <- create.barplot(
data = COMBINED_DATA_PathList, formula = PathogenicList~Normalized,
plot.horizontal = TRUE,
col = COMBINED_DATA_PathList$COLOUR,
yaxis.tck = 0, ylab.axis.padding = 0, ylab.label = '', yaxis.lab = rep('',14),
xlab.label = "Average", xlab.cex = 1, xaxis.cex = 0.75,
add.grid=TRUE, grid.lwd = 1, grid.col = "grey",
xlimits = c(round(min(COMBINED_DATA_PathList$Normalized)/0.25)*0.25-0.2,round(max(COMBINED_DATA_PathList$Normalized)/0.25)*0.25+0.2),
xgrid.at = seq(round(min(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,0.25),
xat = seq(round(min(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,0.25),
ygrid.at = c(seq(1,ncol(MATRIX))+0.5)
)
ParamSim <- create.barplot(
data = COMBINED_DATA_AVG, formula = Normalized~PARAMETERS,
col = COMBINED_DATA_AVG$COLOUR,
sample.order = COMBINED_DATA_AVG$PARAMETERS,
xaxis.tck = 0, xlab.axis.padding = 0, xlab.label = '', xaxis.rot=90, xaxis.cex=0.5,
ylab.label = "Average", ylab.cex = 1, yaxis.cex = 0.75,
add.grid=TRUE, grid.lwd = 1, grid.col = "grey",
ylimits = c(round(min(COMBINED_DATA_AVG$Normalized)/0.25)*0.25-0.2,round(max(COMBINED_DATA_AVG$Normalized)/0.25)*0.25+0.2),
ygrid.at = seq(round(min(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,0.25),
yat = seq(round(min(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,0.25),
xgrid.at = c(cumsum(cnts)+0.5)
)
pdf(file='/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/CLINGEN_HG38_2021/CLINGEN_LOSS_DATA/CLINGEN_LOSS.pdf', width=25, height=10)
create.multiplot(
plot.objects = list(key, HeatMap, PathList, ParamSim),
plot.layout = c(2, 3), layout.skip = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
panel.heights = c(0.25, 1, 0.05), panel.widths = c(1, 0.095),
ylab.padding = 15, x.spacing = 5, y.spacing=0.25,
xaxis.lab = list(
c(-3,3),
rownames(MATRIX),
seq(round(min(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,0.25),
NULL),
xaxis.cex = 0.75, xaxis.fontface = 1, xaxis.rot = 90,
yaxis.lab = list(
NULL,
sub('[_][^_]+$', '', colnames(MATRIX)),
BestParam$PARAMETER,
seq(round(min(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,0.25)),
yaxis.cex = 0.75, yaxis.fontface = 1,
xaxis.top.tck.lab = supp.labs, xat.top = c(cumsum(cnts)-cnts/2+0.5), xlab.top.label = " ", xlab.top.cex = 5.5, main.key.padding = 2,
main=paste0('CLINGEN LOSS (HG38): ', COMBINED_DATA_AVG$PARAMETERS[which(COMBINED_DATA_AVG$Normalized == max(COMBINED_DATA_AVG$Normalized))][1], " (",round(max(COMBINED_DATA_AVG$Normalized),3),")"),
main.cex = 1.25,
xlab.label = c("","","","Benign-Ex Score","","","",""), xlab.cex = 1, xlab.to.xaxis.padding = -5, bottom.padding = 1
)
dev.off()
#sessionInfo()
