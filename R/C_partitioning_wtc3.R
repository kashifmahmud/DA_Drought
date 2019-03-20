################ Figure 7 #####################
# Calculate total C partitioning for individual treatments
#-------------------------------------------------------------------------------------
# cbPalette = c("gray", "orange", "skyblue", "green3", "yellow3", "#0072B2", "#D55E00")
cbPalette = c("gray", "orange", "skyblue", "green3", "#009E73", "yellow3", "#0072B2", "#D55E00")
# vol_group <- list(c(5,10,15),c(20,25,35),1000)
# treat.group = as.factor(c("ambient drought","ambient watered","elevated drought","elevated watered")) # Assign all treatments
Ct.group = data.frame(matrix(ncol = 9, nrow = length(treat.group)))
names(Ct.group) = c("GPP","Rg","Rm","Cs","Cr","Cw","Cf","Cflit","Crlit")

param.summary = result[[2]]
data.summary = result[[3]]
output.summary = result[[4]]
Cstorage.data = result[[7]]
Cstorage.data = Cstorage.data[with(Cstorage.data, order(Date,treatment)), 1:3]
names(Cstorage.data)[1] = "Cs"

data.part = data.all[c("Date","Treatment","GPP")]
names(data.part)[2] = "treatment"
data.part$Date = as.Date(data.part$Date)

Rm.data = subset(output.summary, variable %in% "Rm")
names(Rm.data)[3] = "Rm"

data.part = merge(data.part, Rm.data[,c(1,3,4)], by=c("Date","treatment"))
data.part = merge(data.part, Cstorage.data, by=c("Date","treatment"))

cpool = as.factor(c("Mleaf.modelled","Mwood.modelled","Mroot.modelled","Mlit.modelled","Sleaf.modelled"))
for (i in 1:length(cpool)) {
  cpool.data = subset(output.summary, variable %in% cpool[i])
  cpool.data = cpool.data[, c("Date","value","treatment")]
  cpool.data = cpool.data[with(cpool.data, order(Date,treatment)), ]
  names(cpool.data)[2] = as.character(cpool[i])
  data.part = merge(data.part,cpool.data, all = TRUE)
}

# combine data and parameters together
var = as.factor(c("k","Y","af","as","ar","sf","sr"))
if (no.param.par.var < 5) {
  for (i in 1:length(var)) {
    param = subset(param.summary, variable %in% var[i])
    param = param[, c("Date","Parameter","treatment")]
    param = param[with(param, order(Date,treatment)), ]
    names(param)[2] = as.character(var[i])
    data.part = merge(data.part,param, all = TRUE)
  }
} else {
  for (i in 1:length(var)) {
    for (j in 1:length(unique(param.summary$Date))) {
      if (j < length(unique(param.summary$Date))) {
        param = subset(param.summary, variable %in% var[i] & Date >= unique(param.summary$Date)[j] & Date < unique(param.summary$Date)[j+1])
      } else {
        param = subset(param.summary, variable %in% var[i] & Date >= unique(param.summary$Date)[j])
      }
      param = param[, c("Date","Parameter","treatment")]
      param = param[with(param, order(treatment)), ]
      names(param)[2] = as.character(var[i])
      if (j < length(unique(param.summary$Date))) {
        param.part.set = merge (subset(data.part, Date >= unique(param.summary$Date)[j] & Date < unique(param.summary$Date)[j+1]), param[,c(2:3)], by="treatment")
      } else {
        param.part.set = merge (subset(data.part, Date >= unique(param.summary$Date)[j]), param[,c(2:3)], by="treatment")
      }
      
      if (j == 1) {
        param.part = param.part.set
      } else {
        param.part = rbind(param.part,param.part.set)
      }
    }
    data.part = merge(data.part,param.part, all = TRUE)
  }
}

data.part = data.part[with(data.part, order(treatment,Date)), ]
# data.set = data.part[data.part$Date <= as.Date("2013-05-21"), ]
for (v in 1:length(treat.group)) {
  Ct.group$GPP[v] = sum ( data.part$GPP[which(data.part$treatment %in% treat.group[v])] )
  Ct.group$Rm[v] = sum ( data.part$Rm[which(data.part$treatment %in% treat.group[v])] )
  Ct.group$Cs[v] = data.part$Cs[which(data.part$treatment %in% treat.group[v] & data.part$Date %in% as.Date("2014-02-12"))]
  Ct.group[v, c(5:7)] = data.part[which(data.part$treatment %in% treat.group[v] & data.part$Date %in% as.Date("2014-02-12")), 8:6] - data.part[which(data.part$treatment %in% treat.group[v] & data.part$Date %in% as.Date("2013-09-17")), 8:6]
  Ct.group$Cflit[v] = data.part[which(data.part$treatment %in% treat.group[v] & data.part$Date %in% as.Date("2014-02-12")), 9]
  Ct.group$Crlit[v] = sum ( data.part$sr [which(data.part$treatment %in% treat.group[v])] * data.part$Mroot.modelled [which(data.part$treatment %in% treat.group[v])])
  Ct.group$Rg[v] = Ct.group$GPP[v] - sum(Ct.group[v,c(3:9)])
}

Ct.fraction.group = Ct.group[, c(2:9)]
Ct.fraction.group[,] = Ct.fraction.group[,] / Ct.group[, 1] * 100
# row.names(Ct.fraction.group) <- c("amb-dry","amb-wet","warm-dry","warm-wet")
row.names(Ct.fraction.group) <- c("ambient","warmed")

Ct.group$treatment = treat.group
colnames(Ct.group) <- c("GPP (g C)", "Rg (g C)", "Rm (g C)", "Cs (g C)", "Cr (g C)", "Cw (g C)", "Cf (g C)", "Cflit (g C)", "Crlit (g C)", "Treatment")
Ct.group = Ct.group[,c(10,1,2,3,4,7,5,6,8,9)]
write.csv(Ct.group, file = "output/C_partitioning_wtc3.csv", row.names = FALSE)

png("output/Figure_7_C_partitioning_wtc3.png", units="px", width=1200, height=1000, res=200)
par(mfrow = c(1, 1), mar=c(5, 4, 2, 6))
# bb = barplot(as.matrix(t(Ct.fraction.group)), ylim=c(0, 107), ylab = "C Partitioning (%)", xlab = "Treatments (Container size)",  
#         col = rainbow(20),legend = colnames(Ct.fraction.group), 
#         args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))
bb = barplot(as.matrix(t(Ct.fraction.group)), ylim=c(0, 107), ylab = "C Partitioning (%)", xlab = "Container size (L))",  
             col = cbPalette[1:8],legend = c(expression(R[g]),expression(R["m,tot"]),expression(C[s]),expression(C["s,r"]),expression(C["s,w"]),expression(C["s,f"]),expression(C["f,lit"]),expression(C["r,lit"])), 
             args.legend = list(x = "topright", bty = "n", inset=c(-0.18, 0)))
# text( bb, Ct.fraction.group[,1]-3, labels = round(Ct.group[,3],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]-4, labels = round(Ct.group[,4],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]-1, labels = round(Ct.group[,5],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]-3, labels = round(Ct.group[,6],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]-3, labels = round(Ct.group[,7],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]+Ct.fraction.group[,6]-2, labels = round(Ct.group[,8],1), cex=.9)
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]+Ct.fraction.group[,6]+Ct.fraction.group[,7]-1, labels = round(Ct.group[,9],1), cex=.9)
text( bb, rowSums(Ct.fraction.group)+0.5, labels = round(Ct.group[,2],1), pos = 3, cex=1, col="red")

dev.off()

