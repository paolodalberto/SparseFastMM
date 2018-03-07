R=read.csv("results.txt", sep=",")

R$Complexity = R$N*R$D
R$NPerf      = R$Complexity/R$Time


Ps = unique(R$P)

M = max(R$NPerf)
m = 0
s = log2(unique(R$Complexity))
w = s*0
w[1]=M
w = s*0
w[1]=M


plot(s,w, type='n', main="BLAH",xlab="size", ylab = "GFLOPS")
i = 1
col   = c("red", "blue", "black","darkred", "darkblue","gray")
for (p in Ps) {


    T = subset(R, P==p)	   
    M = max(log2(T$NPerf))
    m = 0
    s = unique(T$N)
    w = s*0
    w[1]=M
    w = s*0
    w[1]=M

    if (FILE) png(paste(p,"processors",".png",sep=""))
   
    plot(s,w, type='n', main=paste(p,"processors" ),xlab="N", ylab = "log2(N*D/time)")
    i = 1
    print(p)
    for (d in unique(T$D)) { 
        TD = subset(T, D==d)
    	str(TD)
        lines(TD$N,
	     log2(TD$NPerf), col=col[i], lty=2, lwd=2)    
        i = i +1
    
    }
    legend("bottomright", legend=paste("degree",unique(T$D)),col=col,lty=2,lwd=2)
    if (FILE) dev.off()
    else X11()
}


#types = unique(R$TYPE)
#algs  = unique(R$NAME)
#col   = c("red", "blue", "black")
#if (FALSE) {
#for (t in types[1]) {
#    print(t);
#    
#    data = subset(R, TYPE==t)
#    M = max(data$HOT)
#    m = 0
#    s = unique(data$SIZE)
#    w = s*0
#    w[1]=M
#    plot(s,w, type='n', main=t,xlab="size", ylab = "GFLOPS")
#    i = 1
#    for (a in algs) {
#       d = subset(data, NAME==a)
#       lines(d$SIZE,d$HOT,type='l',col=col[i],lty=2, lwd=2)
#       i = i + 1
#
#   }
#   legend("bottomleft", legend=algs,col=col,lty=2,lwd=2)
#   X11()
#}
#}
#
#	    
#    M = max(R$HOT)
#    m = 0
#    s = unique(data$SIZE)
#    w = s*0
#    w[1]=M
#    png("Epyc7601.png")
#    plot(s,w, type='n', main="EPYC 7601 32-Core Processor 2.2GHz OpenBLAS GEMM",xlab="size", ylab = "GFLOPS")
#    i = 1
#    label =list()
#    
#    for (a in algs) {
#    	for (t in types) {
#	      
# 	      d = subset(R, NAME==a & TYPE==t)
#	      print(paste(a,t))
#	      str(d)
#       	      lines(d$SIZE,d$HOT,type='l',col=col[i],lty=2, lwd=2)
#	      label[[length(label)+1]] = paste(a,t)
#	      			
#      	      }						    
#	i= i+1
#   }
#   legend("bottomright", legend=algs,col=col,lty=2,lwd=2)
#   dev.off()



