######################################################
spd=function(v1,v2,b)
{
  if (v1<0 | v2<0 | b<0) stop("v1, v2, b should be positive integers.")
  main_plot = matrix(0,b,v1)
  for (i in 1:b) main_plot[i,] = sample(1:v1,v1)
  sub_plot = matrix(0,b,v2)
  for (i in 1:b) sub_plot[i,] = sample(1:v2,v2)
  design = matrix("",b,v1)
  alt.layout = matrix(0,b*v1*v2,3)
  colnames(alt.layout) = c("block","mp trt","sp trt")
  rownames(alt.layout)=NULL
  cnt=0
  for(i in 1:v1)
  {
    random_order=sample(1:b,b)
    for(j in 1:b)
    {
      position=which(main_plot[j,]==i)
      temp=paste('(')
      for(u in 1:v2) 
      {
        sub.trt=sub_plot[random_order[j],u]  
        temp=paste(temp,as.character(sub.trt))
        cnt=cnt+1
        alt.layout[cnt,1]=j
        alt.layout[cnt,2]=i
        alt.layout[cnt,3]=sub.trt
      }    
      temp=paste(temp,')')
      design[j,position]=paste(i,temp)
    }
  }
  alt.layout=alt.layout[order(alt.layout[,1],alt.layout[,2],alt.layout[,3]),]
  parameters=c(v1,b,v2)
  names(parameters)=c("v1","b","v2")
  design = as.data.frame(design)
  colnames(design) <- paste0("Main Plot-", 1:v1)
  rownames(design)=paste0('Block-',1:b)
  output = list(parameters = parameters, design = design, column.layout = alt.layout)
  return(output)
}
ispd.cmis=function(v1,v2,b,k2)
{
  if (v1<0 | v2<0 | b<0 | k2<0) stop("v1,v2,b,k2 should be positive integers.")
  if(!is.wholenumber(b*k2/v2)) stop("b*k2 should be exactly divisible by v2.")
  if(k2 >= v2) stop("k2 should be less than v2")
  main_plot = matrix(0,b,v1)
  for (i in 1:b) main_plot[i,] = sample(1:v1,v1)
  x = ibd(v2,b,k2)
  sub_plot = x$design
  design = matrix("",b,v1)
  alt.layout = matrix(0,b*v1*k2,3)
  colnames(alt.layout) = c("block","mp trt","sp trt")
  rownames(alt.layout)=NULL
  cnt=0
  for(i in 1:v1)
  {
    random_order=sample(1:b,b)
    for(j in 1:b)
    {
      position=which(main_plot[j,]==i)
      temp=paste('(')
      for(u in 1:k2) 
      {
        sub.trt=sub_plot[random_order[j],u]  
        temp=paste(temp,as.character(sub.trt))
        cnt=cnt+1
        alt.layout[cnt,1]=j
        alt.layout[cnt,2]=i
        alt.layout[cnt,3]=sub.trt
      }    
      temp=paste(temp,')')
      design[j,position]=paste(i,temp)
    }
  }
  alt.layout=alt.layout[order(alt.layout[,1],alt.layout[,2],alt.layout[,3]),]
  parameters=c(v1,b,v2,k2)
  names(parameters)=c("v1","b","v2","k2")
  design=as.data.frame(design)
  colnames(design) <- paste0("Main Plot-", 1:v1)
  rownames(design)=paste0('Block-',1:b)
  output = list(parameters = parameters, design = design, column.layout = alt.layout)
  return(output)
}
######################################################
ispd.imcs=function(v1,b,k1,v2)
{
  if (v1<0 | v2<0 | b<0 | k1<0) stop("v1,v2,b,k1 should be positive integers.")
  if(!is.wholenumber(b*k1/v1)) stop("b*k1 should be exactly divisible by v1.")
  if(k1 >= v1) stop("k1 should be less than v1")
  layout_design=NULL
  x=ibd(v1,b,k1)
  design=x$design
  sub_plot =matrix(0,b,v2)
  for (i in 1:b) sub_plot[i,]=sample(1:v2,v2)
  for(i1 in 1:v1)
  {
    a=0 
    random_order=sample(1:b,b)
    for(j in 1:b)
    {
      position=which(design[j,]==i1)
      if (length(position)>0) 
      {
        a=a+1
        temp='('
        for(u in 1:v2) 
        {
            temp=paste(temp,(as.character(sub_plot[random_order[a],u])))
            layout_design=rbind(layout_design,c(j,i1,u)) 
        }  
        temp=paste(temp,')')
        design[j,position]=paste(i1,temp)
      }
    }
  }
  layout_design=layout_design[order(layout_design[,1],layout_design[,2],layout_design[,3]),]
  colnames(layout_design)=c("block","mp","sp")
  col1=v2*b*k1
  rownames(layout_design)=c(1:col1)
  design = as.data.frame(design)
  colnames(design) <- paste0("Main Plot-", 1:k1)
  rownames(design)=paste0('Block-',1:b)
  parameters=paste("v1=",v1,",","b=",b,",","v2=",v2,",","k1=",k1)
  parameters=c(v1,b,k1,v2)
  names(parameters)=c("v1","b","k1","v2")
  output = list(parameters = parameters, design = design, column.layout = layout_design)
  return(output)
}
######################################################
ispd.imis=function(v1,b,k1,v2,k2)
{
  if (v1<0 | v2<0 | b<0 | k1<0 | k2<0) stop("v1,v2,b,k1,k2 should be positive integers.")
  r = b*k1/v1
  if(!is.wholenumber(r)) stop("b*k1 should be exactly divisible by v1.")
  if(k1 >= v1) stop("k1 should be less than v1")
  t = r*k2/v2
  if(!is.wholenumber(t)) stop("b*k1*k2 should be exactly divisible by v1*v2.")
  layout_design=NULL
  x=ibd(v1,b,k1)
  y=ibd(v2,r,k2)
  design=x$design
  sub_plot=y$design
  if(all(diag(y$conc.mat)==t))
  {
    for(i in 1:v1)
    {
      a=0 
      random_order=sample(1:r,r)
      for(j in 1:b)
      {
        position=which(design[j,]==i)
        if (length(position)>0) 
        {
          a=a+1
          temp='('
          for(u in 1:k2) 
          {
            sptrt = sub_plot[random_order[a],u]
            temp=paste(temp,as.character(sptrt))
            layout_design=rbind(layout_design,c(j,i,sptrt)) 
          }
          temp=paste(temp,')')
          design[j,position]=paste(i,temp)
        }
      }
    }
    layout_design=layout_design[order(layout_design[,1],layout_design[,2],layout_design[,3]),]
    colnames(layout_design)=c("block","mp","sp")
    col1=v1*v2*t
    design = as.data.frame(design)
    colnames(design) = paste0("Main Plot-", 1:k1)
    rownames(design)=paste0("Block-",1:b)
    parameters=c(v1,b,k1,v2,k2)
    names(parameters)=c("v1","b","k1","v2","k2")
	  output = list(parameters = parameters, design = design, column.layout = layout_design)
	  return(output)
  } else return("Try again")
}
###################################################################################
ispd = function(v1,v2,b,k1 = NULL,k2 = NULL)
{
  if(missing(k1) & missing(k2)) output = spd(v1,v2,b)
  if(missing(k1) & !missing(k2)) output = ispd.cmis(v1,v2,b,k2)
  if(!missing(k1) & missing(k2)) output = ispd.imcs(v1,b,k1,v2)
  if(!missing(k1) & !missing(k2)) output = ispd.imis(v1,b,k1,v2,k2)
  return(output) 
}
######################################################
aov.ispd.cmis = function(obs, block, mp, sp, y) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,mp))
  m=ncol(x2)
  blockmp=noquote(paste0(block,mp))
  x3=as.matrix(table(obs,blockmp))
  k=colSums(x3)[1] 
  x4=as.matrix(table(obs,sp))
  s=ncol(x4)
  mpsp=noquote(paste0(mp,sp))
  x5=as.matrix(table(obs,mpsp))
  r1=b*m*k/s 
  r=b*k/s 
  x=cbind(one,x1,x2,x3,x4,x5)
  xpx=t(x)%*%x
  y...=t(one)%*%y
  yB..=t(x1)%*%y
  y.M.=t(x2)%*%y
  yBM.=t(x3)%*%y
  blockSS=t(yB..)%*%yB../(m*k)-y...^2/(b*m*k)
  mpSS=t(y.M.)%*%y.M./(b*k)-y...^2/(b*m*k)
  blockmpSS=(t(yBM.)%*%yBM./k-y...^2/(b*m*k))-blockSS-mpSS
  y..S=t(x4)%*%y
  y..S_curl=y..S[1:(s-1)]
  N3=t(x4)%*%x3
  N3_curl=N3[1:(s-1),]
  spSS=(t(y..S_curl-N3_curl%*%yBM./k)%*%solve(diag(s-1)-N3_curl%*%t(N3_curl)/(r1*k))%*%(y..S_curl-(N3_curl%*%yBM.)/k))/r1
  y.MS=t(x5)%*%y
  y.MS_curl=y.MS[-(seq(s,m*s,s))]
  N6=t(x5)%*%x3
  N6_curl=N6[-(seq(s,m*s,s)),]
  mpspSS=t(y.MS_curl-N6_curl%*%yBM./k)%*%solve(r*diag(m*(s-1))-N6_curl%*%t(N6_curl)/(r1*k))%*%(y.MS_curl-N6_curl%*%yBM./k)-spSS
  resSS=t(y)%*%y-t(yBM.)%*%yBM./k-t(y.MS_curl-N6_curl%*%yBM./k)%*%solve(r*diag(m*(s-1))-N6_curl%*%t(N6_curl)/(r1*k))%*%(y.MS_curl-N6_curl%*%yBM./k)                                
  mpspF=(mpspSS/((m-1)*(s-1)))/(resSS/(b*m*k-b*m-m*s+m))
  p.value.mpsp=pf(mpspF, df1=(m-1)*(s-1), df2=b*m*k-b*m-m*s+m, ncp=0, lower.tail = FALSE)
  spF=(spSS/(s-1))/(resSS/(b*m*k-b*m-m*s+m))
  p.value.sp=pf(spF, df1=(s-1), df2=b*m*k-b*m-m*s+m, ncp=0, lower.tail = FALSE)
  mpF=(mpSS/(m-1))/(blockmpSS/((b-1)*(m-1)))
  p.value.mp=pf(mpF, df1=(m-1), df2=(b-1)*(m-1), ncp=0, lower.tail = FALSE)
  blockF=(blockSS/(b-1))/(blockmpSS/((b-1)*(m-1)))
  p.value.block=pf(blockF, df1=(b-1), df2=(b-1)*(m-1), ncp=0, lower.tail = FALSE)
  df=c(b-1,m-1,(b-1)*(m-1),s-1,(m-1)*(s-1),b*m*k-b*m-m*s+m,b*m*k-1)
  totalSS=t(y)%*%y-y...^2/(b*m*k)
  ss=c(blockSS,mpSS,blockmpSS,spSS,mpspSS,resSS,totalSS)
  ms=ss/df
  F.value=c(blockF,mpF,NA,spF,mpspF,NA,NA)
  p.value=c(p.value.block,p.value.mp,NA,p.value.sp,p.value.mpsp,NA,NA)
  anova.table=cbind(df,ss,ms,F.value,p.value)
  rownames(anova.table)=c("Block","A","Error(A)","B","A x B","Error(B)","Total")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
}
#######################################################
aov.ispd.imcs = function(obs, block, mp, sp, y) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,mp))
  m=ncol(x2)
  blockmp=noquote(paste0(block,mp))
  x3=as.matrix(table(obs,blockmp))
  mptbyblockN=as.matrix(table(mp,block))
  k1=length(which(mptbyblockN[,1]>0)) 
  x4=as.matrix(table(obs,sp))
  s=ncol(x4)
  mpsp=noquote(paste0(mp,sp))
  x5=as.matrix(table(obs,mpsp))
  x=cbind(one,x1,x2,x3,x4,x5)
  xpx=t(x)%*%x
  y...=t(one)%*%y
  yB..=t(x1)%*%y
  y.M.=t(x2)%*%y
  yBM.=t(x3)%*%y
  class(yBM.)
  blockSS=t(yB..)%*%yB../(s*k1)-y...^2/(b*k1*s)
  N1=t(x2)%*%x1/s
  N1_curl=N1[1:(m-1),]
  r=rowSums(N1)
  R=diag(r)
  y.M._curl=y.M.[1:(m-1)]
  mpSS=t(y.M._curl-N1_curl%*%yB../k1)%*%solve(R[1:(m-1),1:(m-1)]-N1_curl%*%t(N1_curl)/k1)%*%(y.M._curl-N1_curl%*%yB../k1)/s
  blockmpSS=t(yBM.)%*%yBM./s-y...^2/(b*k1*s)-blockSS-mpSS
  y..S=t(x4)%*%y
  y..S_curl=y..S[1:(s-1)]
  Q_s=(y..S_curl - (y.../s)*matrix(1,(s-1),1))
  C_s=b*k1*(diag(s-1)-matrix(1,s-1,s-1)/s)
  spSS=t(Q_s)%*%solve(C_s)%*%Q_s
  y.MS=t(x5)%*%y
  y.MS_curl=y.MS[-seq(s,m*s,s)]
  N4=t(x5)%*%x3
  N4_curl=N4[-seq(s,m*s,s),]
  Q_MS=(y.MS_curl - N4_curl%*%yBM./s)
  C_MS=R %x% diag(s-1) - N4_curl%*%t(N4_curl)/s
  mpspSS=t(Q_MS)%*%solve(C_MS)%*%Q_MS - spSS
  resSS=t(y)%*%y-t(yBM.)%*%yBM./s-mpspSS-spSS
  mpspF=(mpspSS/((m-1)*(s-1)))/(resSS/(b*k1*s-b*k1-m*s+m))
  p.value.mpsp=pf(mpspF, df1=(m-1)*(s-1), df2=b*k1*s-b*k1-m*s+m, ncp=0, lower.tail = FALSE)
  spF=(spSS/(s-1))/(resSS/(b*k1*s-b*k1-m*s+m))
  p.value.sp=pf(spF, df1=(s-1), df2=b*k1*s-b*k1-m*s+m, ncp=0, lower.tail = FALSE)
  mpF=(mpSS/(m-1))/(blockmpSS/(b*k1-b-m+1))
  p.value.mp=pf(mpF, df1=(m-1), df2=b*k1-b-m+1, ncp=0, lower.tail = FALSE)
  blockF=(blockSS/(b-1))/(blockmpSS/(b*k1-b-m+1))
  p.value.block=pf(blockF, df1=(b-1), df2=b*k1-b-m+1, ncp=0, lower.tail = FALSE)
  df=c(b-1,m-1,b*k1-b-m+1,s-1,(m-1)*(s-1),b*k1*s-b*k1-m*s+m,b*k1*s-1)
  totalSS=t(y)%*%y-y...^2/(b*k1*s)
  ss=c(blockSS,mpSS,blockmpSS,spSS,mpspSS,resSS,totalSS)
  ms=ss/df
  F=c(blockF,mpF,NA,spF,mpspF,NA,NA)
  p.value=c(p.value.block,p.value.mp,NA,p.value.sp,p.value.mpsp,NA,NA)
  anova.table=cbind(df,ss,ms,F, p.value)
  rownames(anova.table)=c("Block","A","Error(A)","B","A x B","Error(B)","Total")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
}  
###################################################################
aov.ispd.imis = function(obs, block, mp, sp, y) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,mp))
  m=ncol(x2)
  blockmp=noquote(paste0(block,mp))
  x3=as.matrix(table(obs,blockmp))
  mptbyblockN=as.matrix(table(mp,block))
  k1=length(which(mptbyblockN[,1]>0)) 
  x4=as.matrix(table(obs,sp))
  s=ncol(x4)
  k2=colSums(x3)[1] 
  mpsp=noquote(paste0(mp,sp))
  x5=as.matrix(table(obs,mpsp))
  x=cbind(one,x1,x2,x3,x4,x5)
  xpx=t(x)%*%x
  y...=t(one)%*%y
  yB..=t(x1)%*%y
  y.M.=t(x2)%*%y
  yBM.=t(x3)%*%y
  blockSS=t(yB..)%*%yB../(k2*k1)-y...^2/(b*k1*k2)
  N1=t(x2)%*%x1/k2
  N1_curl=N1[1:(m-1),]
  r=rowSums(N1)
  R=diag(r)
  y.M._curl=y.M.[1:(m-1)]
  mpSS=t(y.M._curl-N1_curl%*%yB../k1)%*%solve(R[1:(m-1),1:(m-1)]-N1_curl%*%t(N1_curl)/k1)%*%(y.M._curl-N1_curl%*%yB../k1)/k2
  blockmpSS=t(yBM.)%*%yBM./k2-y...^2/(b*k1*k2)-blockSS-mpSS
  y..S=t(x4)%*%y
  y..S_curl=y..S[1:(s-1)]
  N4=t(x4)%*%x3
  N4_curl=N4[1:(s-1),]
  Q_s=(y..S_curl - N4_curl%*%yBM./k2)
  t=r[1]*k2/s
  C_s=m*t*(diag(s-1))-N4_curl%*%t(N4_curl)/k2
  spSS=t(Q_s)%*%solve(C_s)%*%Q_s
  y.MS=t(x5)%*%y
  y.MS_curl=y.MS[-seq(s,m*s,s)]
  N6=t(x5)%*%x3
  N6_curl=N6[-seq(s,m*s,s),]
  Q_MS=(y.MS_curl - N6_curl%*%yBM./k2)
  C_MS=t*diag(m*(s-1)) - N6_curl%*%t(N6_curl)/k2
  mpspSS=t(Q_MS)%*%solve(C_MS)%*%Q_MS - spSS
  resSS=t(y)%*%y-t(yBM.)%*%yBM./k2-mpspSS-spSS
  mpspF=(mpspSS/((m-1)*(s-1)))/(resSS/(b*k1*k2-b*k1-m*s+m))
  p.value.mpsp=pf(mpspF, df1=(m-1)*(s-1), df2=b*k1*k2-b*k1-m*s+m, ncp=0, lower.tail = FALSE)
  spF=(spSS/(s-1))/(resSS/(b*k1*k2-b*k1-m*s+m))
  p.value.sp=pf(spF, df1=(s-1), df2=b*k1*k2-b*k1-m*s+m, ncp=0, lower.tail = FALSE)
  mpF=(mpSS/(m-1))/(blockmpSS/(b*k1-b-m+1))
  p.value.mp=pf(mpF, df1=(m-1), df2=(b*k1-b-m+1), ncp=0, lower.tail = FALSE)
  blockF=(blockSS/(b-1))/(blockmpSS/(b*k1-b-m+1))
  p.value.block=pf(blockF, df1=(b-1), df2=(b*k1-b-m+1), ncp=0, lower.tail = FALSE)
  df=c(b-1,m-1,(b*k1-b-m+1),s-1,(m-1)*(s-1),b*k1*k2-b*k1-m*s+m,b*k1*k2-1)
  totalSS=t(y)%*%y-y...^2/(b*k1*k2)
  ss=c(blockSS,mpSS,blockmpSS,spSS,mpspSS,resSS,totalSS)
  ms=ss/df
  F=c(blockF,mpF,NA,spF,mpspF,NA,NA)
  p.value=c(p.value.block,p.value.mp,NA,p.value.sp,p.value.mpsp,NA,NA)
  anova.table=cbind(df,ss,ms,F, p.value)
  rownames(anova.table)=c("Block","A","Error(A)","B","A x B","Error(B)","Total")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
}  

aov.ispd = function(obs, block, mp, sp, y, incomplete.block = FALSE, incomplete.mp = TRUE)
{
  if(incomplete.block == FALSE & incomplete.mp == TRUE) output = aov.ispd.cmis(obs, block, mp, sp, y)
  if(incomplete.block == TRUE & incomplete.mp == FALSE) output = aov.ispd.imcs(obs, block, mp, sp, y)
  if(incomplete.block == TRUE & incomplete.mp == TRUE) output = aov.ispd.imis(obs, block, mp, sp, y)
  if(incomplete.block == FALSE & incomplete.mp == FALSE) output = "Please use of aov() function or sp.plot() function of agricolae package"
  return(output) 
}