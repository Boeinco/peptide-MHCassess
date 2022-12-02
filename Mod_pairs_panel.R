merge_hlathena$HLAthena <- rowSums(merge_hlathena[,2:9], na.rm=TRUE)
merge_netmhcpan$netMHCpan <- rowSums(merge_netmhcpan[,2:9], na.rm=TRUE)
merge_mhcnuggets$MHCnuggets <- rowSums(merge_mhcnuggets[,2:9], na.rm=TRUE)
merge_mhcflurry$MHCflurry <- rowSums(merge_mhcflurry[,2:9], na.rm=TRUE)

merge_mhcnuggets$Allele <- gsub(":", "", merge_mhcnuggets$Allele)
merge_human_viral <- merge(merge_hlathena[,c(1,11)], merge_mhcflurry[,c(1,12)], by="Allele", all=TRUE)
merge_human_viral <- merge(merge_human_viral, merge_mhcnuggets[,c(1,12)], by="Allele", all=TRUE)
merge_human_viral <- merge(merge_human_viral, merge_netmhcpan[,c(1,11)], by="Allele", all=TRUE)
merge_human_viral <- merge_human_viral[merge_human_viral$Allele %in% list_alleles$Allele,]

max_binders <- merge(hlathena_group_table[,c(1,14)], merge_mhcflurry[,c(1,14)], by="Allele", all=TRUE)
max_binders <- merge(max_binders, merge_mhcnuggets[,c(1,14)], by="Allele", all=TRUE)
max_binders <- merge(max_binders, netmhcpan_group_table[,c(1,14)], by="Allele", all=TRUE)
max_binders[is.na(max_binders)] <- 0
max_binders$tool_avg <- rowMeans(max_binders[,2:5])
options(scipen=999)
write.csv(max_binders, "max_binders_human_viral.csv", row.names = FALSE, quote=FALSE)
max_binders <- read.csv("max_binders_human_viral.csv")
max_binders <- max_binders[max_binders$Allele %in% list_alleles$Allele,]
temp_max <- max_binders[,c(1,6)]
temp_max$tool_avg <- temp_max$tool_avg*100
write.csv(temp_max, "avg_percent_binders_human_viral.csv", row.names = FALSE, quote=FALSE)


colnames(merge_human_viral)[2] <- "HLAthena"
merge_human_viral[is.na(merge_human_viral)] <- 0
merge_human_viral <- merge_human_viral[,c(1,2,5,3,4)]
#HLAcol <- ifelse(grepl("A", merge_alleles$Allele), '#1B9E77', ifelse(grepl("B", merge_alleles$Allele), '#D95F02', '#7570B3'))
HLAcol <- ifelse(grepl("A", merge_human_viral$Allele), '#1B9E77', ifelse(grepl("B", merge_human_viral$Allele), '#D95F02', '#7570B3'))
library(ggplot2)
library(psych)

setwd("/Users/nguyenau/Documents/Coding/PeptideAlleleBinding/Results/netmhcpan_full")
merge_alleles <- read.csv("merge_alleles_random_0.5.csv")
merge_alleles <- read.csv("merge_alleles_random_0.6.csv")
merge_alleles <- read.csv("merge_alleles_random_0.7.csv")

HLAcol <- ifelse(grepl("A", merge_alleles$Allele), '#1B9E77', ifelse(grepl("B", merge_alleles$Allele), '#D95F02', '#7570B3'))
merge_alleles[is.na(merge_alleles)] <- 0

pairs.panels(merge_alleles[,-1], 
             method = "spearman", # correlation method
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

pairs.panels(merge_human_viral[,-1], 
             method = "spearman", # correlation method
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


"pairs.panels2" <- function (x, smooth = TRUE, scale = FALSE, density=TRUE,ellipses=TRUE,digits = 2, method="pearson",pch = 20,lm=FALSE,cor=TRUE,jiggle=FALSE,factor=2,hist.col="cyan",show.points=TRUE,rug=TRUE, breaks="Sturges", cex.cor = 1 ,wt=NULL,smoother=FALSE,stars=FALSE,ci=FALSE,alpha=.05,...)   #combines a splom, histograms, and correlations
  {
   
    "panel.hist.density" <- function(x, ...)
      {
        usr <- par("usr")
        par(usr = c(0, max(x), 0, 1.5) ) #consider max(x) instead of global max (~38000) 
        nx <- length(x)
        x.sort <- order(x)
        rect(0, 1:nx/nx, x[x.sort], (1:nx-1)/nx, border=NA, col = HLAcol[x.sort], ...)
      }
    
    
    
    "panel.cor" <-
      function(x, y, prefix="",...)  {
        
        usr <- par("usr"); on.exit(par("usr"))
        par(usr = c(0, 1, 0, 1))
        if(is.null(wt)) { r  <- cor(x, y,use="pairwise",method=method)} else {
          r <- cor.wt(data.frame(x,y),w=wt[,c(1:2)])$r[1,2]}
        txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(stars) {pval <- r.test(sum(!is.na(x*y)),r)$p
        symp <- symnum(pval, corr = FALSE,cutpoints = c(0,  .001,.01,.05, 1),
                       symbols = c("***","**","*"," "),legend=FALSE)
        txt <- paste0(txt,symp)}
        cex <- cex.cor*0.8/(max(strwidth("0.12***"),strwidth(txt)))
        if(scale)  {cex1 <- cex  * abs(r)
        if(cex1 < .25) cex1 <- .25 #otherwise they just vanish
        text(0.5, 0.5, txt, cex = cex1) } else {
          text(0.5, 0.5, txt,cex=cex)}
      }
    
    "panel.smoother" <- 
      function (x, y,pch = par("pch"), 
                col.smooth = "red", span = 2/3, iter = 3, ...) 
      {
        # usr <- par("usr"); on.exit(par(usr))
        #  par(usr = c(usr[1]-abs(.05*usr[1]) ,usr[2]+ abs(.05*usr[2])  , usr[3],usr[4]) )     #doensn't affect the axis correctly
        xm <- mean(x,na.rm=TRUE)
        ym <- mean(y,na.rm=TRUE)
        xs <- sd(x,na.rm=TRUE)
        ys <- sd(y,na.rm=TRUE)
        r = cor(x, y,use="pairwise",method=method)
        if(jiggle) { x <- jitter(x,factor=factor)
        y <- jitter(y,factor=factor)}
        if(smoother) {smoothScatter(x,y,add=TRUE, nrpoints=0)} else {if(show.points)  points(x, y, pch = pch, ...)}
        
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) {   
          if(smooth & ci) {   lml <- loess(y~x ,degree=1,family="symmetric") 
          tempx <- data.frame(x = seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length.out=47))
          pred <-  predict(lml,newdata=tempx,se=TRUE ) 
          
          if(ci) {  upperci <- pred$fit + confid*pred$se.fit
          lowerci <- pred$fit - confid*pred$se.fit 
          polygon(c(tempx$x,rev(tempx$x)),c(lowerci,rev(upperci)),col=adjustcolor("light grey", alpha.f=0.8), border=NA)
          }
          lines(tempx$x,pred$fit,  col = col.smooth, ...)   #this is the loess fit
          }  else {if(smooth)  lines(stats::lowess(x[ok],y[ok],f=span,iter=iter),col=col.smooth) }}
        if(ellipses)  draw.ellipse(xm,ym,xs,ys,r,col.smooth=col.smooth,...)  #this just draws the ellipse 
      }
    
    "panel.lm" <- 
      function (x, y,  pch = par("pch"), 
                col.lm = "red",  ...) 
      {   ymin <- min(y)
      ymax <- max(y)
      xmin <- min(x)
      xmax <- max(x)
      ylim <- c(min(ymin,xmin),max(ymax,xmax))
      xlim <- ylim
      if(jiggle) { x <- jitter(x,factor=factor)
      y <- jitter(y,factor=factor)}
      if(smoother) {smoothScatter(x,y,add=TRUE, nrpoints=0)} else {if(show.points) {points(x, y, pch = pch,ylim = ylim, xlim= xlim, ...)}}# if(show.points) points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)
      ok <- is.finite(x) & is.finite(y)
      if (any(ok)) {
        lml <- lm(y ~ x)  
        
        
        if(ci) {
          tempx <- data.frame(x = seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length.out=47))
          pred <-  predict.lm(lml,newdata=tempx,se.fit=TRUE)  #from Julian Martins 
          upperci <- pred$fit + confid*pred$se.fit
          lowerci <- pred$fit - confid*pred$se.fit
          polygon(c(tempx$x,rev(tempx$x)),c(lowerci,rev(upperci)),col=adjustcolor("light grey", alpha.f=0.8), border=NA)
        }
        if(ellipses) {
          xm <- mean(x,na.rm=TRUE)
          ym <- mean(y,na.rm=TRUE)
          xs <- sd(x,na.rm=TRUE)
          ys <- sd(y,na.rm=TRUE)
          r = cor(x, y,use="pairwise",method=method)
          draw.ellipse(xm,ym,xs,ys,r,col.smooth=col.lm,...)   #just draw the ellipse
        }
        abline(lml, col = col.lm, ...)
      }
      }
    
    
    "draw.ellipse" <-  function(x=0,y=0,xs=1,ys=1,r=0,col.smooth,add=TRUE,segments=51,...) {
      #based upon John Fox's ellipse functions
      angles <- (0:segments) * 2 * pi/segments
      unit.circle <- cbind(cos(angles), sin(angles))
      if(!is.na(r)) {
        if (abs(r)>0 )theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
        
        shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
        ellipse <- unit.circle %*% shape 
        ellipse[,1] <- ellipse[,1]*xs + x
        ellipse[,2] <- ellipse[,2]*ys + y
        if(show.points) points(x,y,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
        lines(ellipse, ...)   }    
    }
    
    "panel.ellipse" <-
      function (x, y,   pch = par("pch"), 
                col.smooth = "red", ...) 
      { segments=51
      usr <- par("usr"); on.exit(par("usr"))
      par(usr = c(usr[1]-abs(.05*usr[1]) ,usr[2]+ abs(.05*usr[2])  , 0, 1.5) ) 
      xm <- mean(x,na.rm=TRUE)
      ym <- mean(y,na.rm=TRUE)
      xs <- sd(x,na.rm=TRUE)
      ys <- sd(y,na.rm=TRUE)
      r = cor(x, y,use="pairwise",method=method)
      if(jiggle) { x <- jitter(x,factor=factor)
      y <- jitter(y,factor=factor)}
      if(smoother) {smoothScatter(x,y,add=TRUE, nrpoints=0)} else {if(show.points) {points(x, y, pch = pch, ...)}}
      
      angles <- (0:segments) * 2 * pi/segments
      unit.circle <- cbind(cos(angles), sin(angles))
      if(!is.na(r)) {
        if (abs(r)>0 ) theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
        
        shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
        ellipse <- unit.circle %*% shape 
        ellipse[,1] <- ellipse[,1]*xs + xm
        ellipse[,2] <- ellipse[,2]*ys + ym
        points(xm,ym,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
        if(ellipses) lines(ellipse, ...) 
      }    
      }

    old.par <- par(no.readonly = TRUE) # save default, for resetting... 
    on.exit(par(old.par))     #and when we quit the function, restore to original values
    
    
    if(missing(cex.cor)) cex.cor <- 1   #this allows us to scale the points separately from the correlations 
    
    for(i in 1:ncol(x)) {  #treat character data as numeric
      if(is.character(x[[i]] ))  { x[[i]] <- as.numeric(as.factor(x[[i]]) )
      colnames(x)[i] <- paste(colnames(x)[i],"*",sep="")}
    }
    n.obs <- nrow(x)     
    confid <- qt(1-alpha/2,n.obs-2)   #used in finding confidence intervals for regressions and loess
    
    if(!lm) { #the basic default is here
      if(cor) {
        pairs2(x, diag.panel = panel.hist.density, upper.panel = panel.cor
              , lower.panel = panel.smoother, pch=pch, ...)} else {
                pairs2(x, diag.panel = panel.hist.density, upper.panel = panel.smoother, lower.panel = panel.smoother, pch=pch, ...)} 
      
    } else { #lm is TRUE
      if(!cor)  { #this case does not show the correlations, but rather shows the regression lines above and below the diagonal
        pairs2(x, diag.panel = panel.hist.density, upper.panel = panel.lm, lower.panel = panel.lm, pch=pch, ...)   
      } else {  #the normal case is to show the regressions below and the rs above
        pairs2(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.lm,pch=pch,  ...)   
        
      }
    }
    
  }   #end of pairs.panels 
###


"histo" <- function(x,breaks="Sturges", ...) {  
  tax <- table(x)
  if(length(tax) < 11) {breaks <- as.numeric(names(tax))
  y <- tax/max(tax)
  interbreak <- min(diff(breaks))*(length(tax)-1)/21
  rect(breaks-interbreak,0,breaks + interbreak,y)
  } else {
    
    h <- hist(x,breaks=breaks)
  }}

