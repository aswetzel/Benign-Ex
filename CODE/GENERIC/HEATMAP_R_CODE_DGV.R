library(bestNormalize)
library(BoutrosLab.plotting.general)
library(tidyverse)

COMBINED_DATA <- read.csv('/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/DGV_HG38_2020/DGV_LOSS_DATA/COMBINED_DATA_LOSS.txt',sep="\t",header=FALSE)
colnames(COMBINED_DATA) <- c("PathogenicList", "PARAMETERS", "X")
ordernorm_obj <- orderNorm(COMBINED_DATA$X)
COMBINED_DATA$Normalized <- ordernorm_obj$x.t
COMBINED_DATA$EVENTS <- vapply(strsplit(as.character(COMBINED_DATA$PARAMETERS),"_"), `[`, 1, FUN.VALUE=character(1))
COMBINED_DATA$EVENTS <- as.numeric(as.character(COMBINED_DATA$EVENTS))
COMBINED_DATA$N <- vapply(strsplit(as.character(COMBINED_DATA$PARAMETERS),"_"), `[`, 2, FUN.VALUE=character(1))
COMBINED_DATA$AF <- vapply(strsplit(as.character(COMBINED_DATA$PARAMETERS),"_"), `[`, 3, FUN.VALUE=character(1))
COMBINED_DATA$MTHD <- vapply(strsplit(as.character(COMBINED_DATA$PARAMETERS),"_"), `[`, 4, FUN.VALUE=character(1))
COMBINED_DATA$YR <- vapply(strsplit(as.character(COMBINED_DATA$PARAMETERS),"_"), `[`, 5, FUN.VALUE=character(1))
COMBINED_DATA$MORE <- paste(COMBINED_DATA$MTHD,COMBINED_DATA$YR,sep = "_")
COMBINED_DATA <- with(COMBINED_DATA, COMBINED_DATA[order(EVENTS, N, MTHD, YR, N, AF),])
cnts <- as.data.frame(tail(count(COMBINED_DATA, PathogenicList, EVENTS)[,3],length(unique(COMBINED_DATA$EVENTS))))[,1]

COMBINED_DATA_AVG <- aggregate(Normalized~PARAMETERS, COMBINED_DATA, mean)
COMBINED_DATA_AVG$EVENTS <- vapply(strsplit(as.character(COMBINED_DATA_AVG$PARAMETERS),"_"), `[`, 1, FUN.VALUE=character(1))
COMBINED_DATA_AVG$EVENTS <- as.numeric(as.character(COMBINED_DATA_AVG$EVENTS))
COMBINED_DATA_AVG <- arrange(COMBINED_DATA_AVG, EVENTS)
COMBINED_DATA_AVG$COLOUR <- "NA"
COMBINED_DATA_AVG$COLOUR[COMBINED_DATA_AVG$Normalized < 0] <- "firebrick1"
COMBINED_DATA_AVG$COLOUR[COMBINED_DATA_AVG$Normalized >= 0] <- "dodgerblue1"
COMBINED_DATA_AVG$N <- vapply(strsplit(as.character(COMBINED_DATA_AVG$PARAMETERS),"_"), `[`, 2, FUN.VALUE=character(1))
COMBINED_DATA_AVG$AF <- vapply(strsplit(as.character(COMBINED_DATA_AVG$PARAMETERS),"_"), `[`, 3, FUN.VALUE=character(1))
COMBINED_DATA_AVG$MTHD <- vapply(strsplit(as.character(COMBINED_DATA_AVG$PARAMETERS),"_"), `[`, 4, FUN.VALUE=character(1))
COMBINED_DATA_AVG$YR <- vapply(strsplit(as.character(COMBINED_DATA_AVG$PARAMETERS),"_"), `[`, 5, FUN.VALUE=character(1))
COMBINED_DATA_AVG <- with(COMBINED_DATA_AVG, COMBINED_DATA_AVG[order(EVENTS, N, MTHD, YR, AF),])
COMBINED_DATA_PathList <- aggregate(Normalized~PathogenicList, COMBINED_DATA, mean)
COMBINED_DATA_PathList$COLOUR <- "NA"
COMBINED_DATA_PathList$COLOUR[COMBINED_DATA_PathList$Normalized < 0] <- "firebrick1"
COMBINED_DATA_PathList$COLOUR[COMBINED_DATA_PathList$Normalized >= 0] <- "dodgerblue1"

MATRIX <- as.data.frame(pivot_wider(COMBINED_DATA[,c(1,2,4)], names_from=PathogenicList, values_from=Normalized))
rownames(MATRIX) <- MATRIX[,1]
MATRIX <- MATRIX[,-1]
MATRIX2 <- MATRIX
MATRIX2$AVERAGE <- rowMeans(MATRIX2)
write.csv(MATRIX2,'/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/DGV_HG38_2020/DGV_LOSS_DATA/MATRIX_LOSS.csv')
BestParam <- cbind(as.data.frame(apply(MATRIX, 2, which.max)), as.data.frame(apply(MATRIX,2,max)))
colnames(BestParam) <- c("ROW","SCORE")
BestParam$PARAMETER <- rownames(MATRIX)[BestParam$ROW]
BestParam$ROW <- NULL
write.csv(BestParam,'/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/DGV_HG38_2020/DGV_LOSS_DATA/LIST_SCORES_LOSS.csv')

AF <- vapply(strsplit(as.character(rownames(MATRIX)),"_"), `[`, 3, FUN.VALUE=character(1))
if( length(AF)<500 ) {
AF2 <- c(rep("",25), AF)[match(AF, c(0.002, 0.003, 0.004, 0.006, 0.007, 0.008, 0.009, 0.02, 0.03, 0.04, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 2,3,4,5, AF))]
} else {
AF2 <- c(rep("",28), AF)[match(AF, c(0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2,3,4,5, AF))]}
a <- strsplit(rownames(MATRIX),split="_")
b <- c("EVENTS=","N=","METHOD=","YEAR=")
rownme <- NULL
N <- NULL
YRS <- NULL
MTHDS <- NULL
EVNTS <- NULL
for(value in a) {
if (length(value)==4) { value <- c(value,"NA") }
EVNTS <- c(EVNTS, value[1])
N <- c(N, value[2])
MTHDS <- c(MTHDS, value[4])
YRS <- c(YRS, value[5]) }
EVNTS <- unique(EVNTS)
N <- unique(N)
YRS <- unique(YRS)
YRS <- YRS[c(length(YRS),1:length(YRS)-1)]
MTHDS <- unique(MTHDS)
if (length(N)==1 & length(MTHDS)==1 & length(YRS)==1) {
for(value in a) {
addme <- paste(paste(b,value[-c(3)],sep=""),collapse="\n")
rownme <- c(rownme,addme)
}
rownme <- unique(rownme) } else {
for(value in a) {
c <- c(value[1],paste(N,collapse="/"),paste(MTHDS,collapse="/"),paste(YRS,collapse="/"))
addme <- paste(paste0(b,c,collapse="\n"),collapse="")
rownme <- c(rownme,addme)
}}
rownme <- unique(rownme)

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
pdf(file='/home/darbro/Desktop/BENIGNEX/OUTPUT/SMALLESTRUN/DGV_HG38_2020/DGV_LOSS_DATA/DGV_LOSS.pdf', width=25, height=10)
create.multiplot(
plot.objects = list(key, HeatMap, PathList, ParamSim),
plot.layout = c(2, 3), layout.skip = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
panel.heights = c(0.25, 1, 0.05), panel.widths = c(1, 0.095),
ylab.padding = 15, x.spacing = 5, y.spacing=0.25,
xaxis.lab = list(
c(-3,3),
AF2,
seq(round(min(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_PathList$Normalized)/0.25)*0.25,0.25),
NULL),
xaxis.cex = 0.75, xaxis.fontface = 1, xaxis.rot = 90,
yaxis.lab = list(
NULL,
sub('[_][^_]+$', '', colnames(MATRIX)),
BestParam$PARAMETER,
seq(round(min(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,round(max(COMBINED_DATA_AVG$Normalized)/0.25)*0.25,0.25)),
yaxis.cex = 0.75, yaxis.fontface = 1,
xaxis.top.tck.lab = rownme, xat.top = c(cumsum(cnts)-cnts/2+0.5), xlab.top.label = " ", xlab.top.cex = 5.5, main.key.padding = 2,
main=paste0('DGV LOSS (HG38): ', COMBINED_DATA_AVG$PARAMETERS[which(COMBINED_DATA_AVG$Normalized == max(COMBINED_DATA_AVG$Normalized))][1], " (",round(max(COMBINED_DATA_AVG$Normalized),3),")"),
main.cex = 1.25,
xlab.label = c("","","","Benign-Ex Score","","","",""), xlab.cex = 1, xlab.to.xaxis.padding = -5
)
dev.off()
#sessionInfo()

