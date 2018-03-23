######################################################
#                   DESCRIPTION

# The goal of this script it to create density plot from the average
# frequency of generated txt files by the program
# random_fasta_comparsion_high_N_low.py
setwd("/home/nicolas/PycharmProjects/2018_03_02_Difference_between_Up_N_Down_random_sequence")

list_file = dir("result/")


pdf("result/density.pdf", onefile=TRUE, paper="a4")
for(file in list_file){
  print(file)
  file = paste("result/",file,sep="")
  file_list = strsplit(file, "_")[[1]]
  patron = file_list[7]
  unit = file_list[2]
  high = as.numeric(strsplit(file_list[9], ":")[[1]][2])
  low = as.numeric(gsub(".tsv", "", strsplit(file_list[10], ":")[[1]][2]))
  dt <- read.csv(file, sep="\t", h=T)
  dt<- dt[0:100, ]
  names(dt) <- c("V1", "V2", "V3")
  dt$V1 <- as.numeric(as.vector(dt$V1))
  dt$V2 <- as.numeric(as.vector(dt$V2))
  xrange = c(min(c(dt[,1], dt[,2])),  max(c(dt[,1], dt[,2])))
  a = density(dt[,1])
  b = density(dt[,2])
  yrange=c(0, max(c(a$y,b$y)))
  plot(density(dt[,1]), xlim=xrange, ylim=yrange, col="red", xlab=paste("frequency of ", unit), main=paste("distribution of mean frequency of ",unit,"\n in the enriched and impovrished fasta files\nGenesis : ", patron, " high freq : ", high, " low-freq : ", low))
  lines(density(dt[,2]), col="blue")
  abline(v=(high * 100), lwd=2, col="red", lty=3)
  abline(v=(low * 100), lwd=2, col="blue", lty=3)

  lg1 = paste("density of mean frequency of ", unit, " in low fasta")
  lg2 = paste("fixed frequency of ", unit, " for low fasta")
  lg3 = paste("density of mean frequency of ", unit, " in high fasta")
  lg4 = paste("fixed frequency of ", unit, " for high fasta")
  legend("topright", col=c("blue", "blue", "red", "red"), lty=c(1,3,1,3), lwd=c(1,2,1,2), legend=c(lg1, lg2, lg3, lg4), cex=0.5)

}
dev.off()
