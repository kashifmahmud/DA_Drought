# Matching C balance of the entire experiment considering C inputs and outputs
C.balance = data.frame(matrix(ncol = 11, nrow = length(treat.group)))
names(C.balance) = c("Treatment","GPP","Ra","Rm_root","Rg_root","Cs","Clit_foliage","Cn","Clit_root","C.output","C.imbalance")
C.balance$Treatment = treat.group
  
for (v in 1:length(treat.group)) {
  data.set = subset(data.all,(Treatment %in% treat.group[v]))
  # data.set[nrow(data.set),c(10:17)] = data.set[nrow(data.set)-1,c(10:17)]
  data.set[,c("LM","WM","RM","litter")] = na.spline(data.set[,c("LM","WM","RM","litter")])
  # plot(data.set$Date, data.set$LM)
  
  C.balance$GPP[v] = sum(data.set$GPP) 
  C.balance$Ra[v] = sum(data.set$Ra) 
  C.balance$Cs[v] = (data.set$LM[nrow(data.set)] - data.set$LM[1]) + (data.set$WM[nrow(data.set)] - data.set$WM[1]) + (data.set$RM[nrow(data.set)] - data.set$RM[1])
  C.balance$Rm_root[v] = sum(data.set$Rd.fineroot.mean*data.set$RM*data.set$FRratio + data.set$Rd.intermediateroot.mean*data.set$RM*data.set$IRratio + 
    data.set$Rd.coarseroot.mean*data.set$RM*data.set$CRratio + data.set$Rd.boleroot.mean*data.set$RM*data.set$BRratio)
  C.balance$Rg_root[v] = 0.3 * (data.set$RM[nrow(data.set)] - data.set$RM[1])
  C.balance$Clit_foliage[v] = data.set$litter[nrow(data.set)] - data.set$litter[1]
  C.balance$Cn[v] = data.set$TNC_tot[max(which(complete.cases(data.set$TNC_tot)))] - data.set$TNC_tot[min(which(complete.cases(data.set$TNC_tot)))]
  C.balance$Clit_root[v] = 0.1 * C.balance$Clit_foliage[v]
    
  C.balance$C.output[v] = C.balance$Ra[v] + C.balance$Cs[v] + C.balance$Rm_root[v] + C.balance$Rg_root[v] + C.balance$Clit_foliage[v] + C.balance$Cn[v] + C.balance$Clit_root[v]
  C.balance$C.imbalance[v] = C.balance$GPP[v] - C.balance$C.output[v]
}

C.balance.fraction = C.balance[, c(3:9)]
C.balance.fraction[,] = C.balance.fraction[,] / C.balance[,2] * 100
row.names(C.balance.fraction) <- treat.group
# C.balance.fraction = abs(C.balance.fraction)
  
C.balance = C.balance[,-c(10,11)]
colnames(C.balance) <- c("Treatment", "GPP (g C)", "Ra (g C)", "Rm_root (g C)", "Rg_root (g C)", "Cs (g C)", "Clit_foliage (g C)", "Cn (g C)", "Clit_root (g C)")
# C.balance = C.balance[,c(10,1,2,3,4,7,5,6,8,9)]
write.csv(C.balance, file = "output/C_partitioning_wtc3.csv", row.names = FALSE)


cbPalette = c("gray", "orange", "skyblue", "green3", "#009E73", "yellow3", "#0072B2", "#D55E00")
png("output/Figure_1a_C_balance_wtc3.png", units="px", width=1200, height=1000, res=200)
par(mfrow = c(1, 1), mar=c(5, 4, 2, 6))
# bb = barplot(as.matrix(t(Ct.fraction.group)), ylim=c(0, 107), ylab = "C Partitioning (%)", xlab = "Treatments (Container size)",  
#         col = rainbow(20),legend = colnames(Ct.fraction.group), 
#         args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))
bb = barplot(as.matrix(t(C.balance.fraction)), ylim=c(0, 110), ylab = "C Partitioning (%)", xlab = "Container size (L))",  
             col = cbPalette[1:7],legend = c(expression(R[a]),expression(R["m,root"]),expression(R["g,root"]),expression(C[s]),
                                             expression(C["lit,foliage"]),expression(C[n]),expression(C["lit,root"])), 
             args.legend = list(x = "topright", bty = "n", inset=c(-0.22, 0)))
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]+Ct.fraction.group[,6]+Ct.fraction.group[,7]-1, labels = round(Ct.group[,9],1), cex=.9)
text( bb, rowSums(C.balance.fraction)+0.5, labels = round(C.balance[,2],1), pos = 3, cex=1, col="red")

dev.off()
