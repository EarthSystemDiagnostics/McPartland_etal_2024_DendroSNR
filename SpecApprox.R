SpecApprox<-function(spec, xout=NULL,...){
  if(is.null(xout)){
    #frq.bnds<-log10(range(spec$freq))
    frq.bnds<-c(1/1000,1/2) #bin output by regular frequency bands
    xout<-seq(frq.bnds[1],frq.bnds[2], 0.001)
  }
  int.spec<-approx(x = spec$freq,y = spec$spec,xout = xout,...)
  xnt.spec<-approx(x=spec$freq, y=spec$dof, xout = xout,...)
  #znt.spec<-approx(x=spec$freq, y=spec$Scales, xout = xout,...)
  rtrn.spec<-list(freq=int.spec$x,
                  spec=int.spec$y,
                  dofs=xnt.spec$y)
  #scales=znt.spec$y)
  class(rtrn.spec)<-"spec"
  return(rtrn.spec)
}