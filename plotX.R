plotX=function(X,tcm,title="",xl="",yl="",sp=1){
I=dim(X)[1]  
color=rainbow(length(table(tcm)))
plot(X[1,], type="l", col=color[tcm[1]],ylim=c(min(X),max(X)),
     ylab=yl,xlab=xl,main=title,lwd=3)
for (i in 1:I){
  lines(X[i,], col=color[tcm[i]],lwd=sp)
}
}