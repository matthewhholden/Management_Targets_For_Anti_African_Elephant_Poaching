##### Plots in PNAS Manuscript
##############################



# Fig 1A - Baseline bar plot ####
    #plot Elephant dynamics under 3 management options [baseline parameterization: legal trade vs. no trade]
      col.width = 3.42
    
      pdf('ElephantEquilibria.pdf', width = col.width, height = .75*col.width)
        lw = 1.3; txts = .8; ticks = 1; xMax =20; yMax=k; ca=.6;
        labs1 = c('status-quo', 'legal trade', ' ')
        labs2 = c(' ', ' ', 'legal trade funds\nextra enforcement')
        par( mfrow = c(1,1), oma = c(0,0,0,0), mar=c(1.65,2,.75,.35) )
        bp = barplot(c(out, outNFE, outFE), yaxt = 'n', ylim = c(0,.22) )
        mtext(side=2, 
              expression(paste("Equilibrium elephant denisty [ind / km"^"2","]")), 
              cex=txts, padj=-1.1)
        axis(side=2, seq(0,.2,by=.05), cex.axis=ca, padj =2.6, tck=-0.015)
        abline(h=N0, lty = 1, lwd = lw)
        abline(h=k, lty = 2, lwd = lw)
        text(x=bp, y= 0, labels=labs1, xpd = TRUE,  adj = c(0.5,2.2), cex=txts)
        text(x=bp, y= 0, labels=labs2, xpd = TRUE,  adj = c(0.5,1.2), cex=txts)
        text(3.1,N0+.009,  'current density', cex = txts)
        text(3.05,k+.01, 'carrying capacity', cex = txts)
        box()
     dev.off()
    
    
# Fig 1B - Equilibrium density vs. shift in demand curve ####
  
    pdf( file = 'demandShift_baseline.pdf', width = 8, height = 5)
      par( mfrow=c(1,1), oma = c(10,4,2.3,1), mar=c(1,1,0,0) )
      lw = 3; sizeS=2; ca=1.3; ca2=1.1; 
      plot(shift, N.eq.demand, 
           xlim = c(0, max(shift)), ylim = c(0,.2),
           xlab='', ylab='', 
           xaxt='n', yaxt='n',
           type='l', pch =1, lty=1,  col=1, lwd=lw, cex= sizeS)
      axis(side = 1, at = seq(0,p0/b,300), cex.axis=ca)
      axis(side = 2, at = seq(0,k,.05), cex.axis=ca)
      mtext(side=1, 'Demand reduction ( elephants / day )', cex=ca, padj = 2.5)
      mtext(side=2, expression(paste("Equilibrium density (elephants/km"^"2",")") ),
            cex=ca, padj = -1.6, outer=TRUE)
      mtext(side=1, 'Demand [elephants/day]', outer=TRUE, cex = ca2, padj=14)
      mtext(side=2, 'Price [USD]', outer=TRUE, cex = ca2, adj=-1.08, padj=-3.5)
    dev.off()
    
    #plot demand curves to go below the x-axis
      for(i in 1:5){
        PlotDemandCurve(p0, p, QD, 300*(i-1), b, pa=i)
      }
      
      PlotDemandCurve = function(p0, pBaseline, QDBaseline, shift, b, pa){
        lw=15; cs = 2.5; ca = 3
        p1 = p0 - shift*b
        xlim = c(0,p0/b)
        ylim = c(0,p0)
        x1 = c(0, p0/b)
        x2 = c(0, p1/b)
        filename = paste('SupplyDemandShift', ceiling(shift),'.png')
        pdf(filename)
        par( mfrow = c(1,1), oma = c(2,2,2,1), mar=c(1,1,1,1), xaxs="i" )
        plot(x1, p0-b*x1, type = 'l', col = 'red', lty = 2, lwd=lw, axes = FALSE)
        points(QDBaseline, pBaseline, col='black',lwd=3, cex=6, pch=16)
        lines(x2 , p1 - b*x2, col = 'blue', lty = 1, lwd=lw)
        y = 1000
        arrows( (p0-y)/b - 70, y,
                (p1-y)/b + 70, y, lwd = lw)
        atx = seq(0,1200, by = 400)
        aty = seq(0,24000, by = 8000)
        axis(side =1, at = atx, labels = prettyNum(atx), cex.axis =1.5)
        if(pa==1){
          axis(side =2, at = aty, labels = prettyNum(aty, big.mark=','), cex.axis =ca)
        } 
        box()
        dev.off()
      }
    
# Fig 1C - Baseline Both shift and punish ####
  
  
  library(fields)
  col.width = 3.42
  
  pdf('Baseline_Punish_shift2.pdf', height = 3.0, width = col.width)
    par(mfrow=c(1,1), oma = c(1.5,1.6,0,1), mar=c(.5,.5,.5, 2), las=0 )
    lw = 2; sizeS=2; ca=.6; cxt = .8; ca2=.5;
    image(N.eq,            
          xlab='',  ylab='', 
          xaxt='n', yaxt='n',
          col = gray(seq(0,1, by =1/np)) )
    box()
    #plot the line where Equilibrium density = N0
    point1 = which(  (N.eq[1,] - N0)^2 == min( (N.eq[1,] - N0)^2 ) )/np
    point2 = which(  (N.eq[,1] - N0)^2 == min( (N.eq[,1] - N0)^2 ) )/np
    
    point1k = which(  (N.eq[1,] - k)^2 == min( (N.eq[1,] - k)^2 ) )[1]/np
    point2k = which(  (N.eq[,1] - k)^2 == min( (N.eq[,1] - k)^2 ) )[1]/np
    
    
    lines( c(0, point2), c(point1,0), lwd = lw )
    lines( c(0, point2k), c(point1k,0), lwd = lw, lty=2 )
    
    slope = point1/(-point2)
    
    demand.shift.arrow = 200
    arrows(0, (s0*pc)/pmax,  demand.shift.arrow / (p0/b) , (s0*pc)/pmax, col='yellow',lwd=lw )
    arrows(demand.shift.arrow / (p0/b),     (s0*pc)/pmax,  
           demand.shift.arrow / (p0/b) ,    point1 + demand.shift.arrow/(p0/b)*slope, 
           col='yellow',lwd=lw )
    
    #add baseline point
    points(0, (s0*pc)/pmax , col = 'black', cex=2, pch =19, xpd=TRUE)
    points(0, (s0*pc)/pmax , col = 'yellow', cex=1.5, pch =19, xpd=TRUE)
    
    #add descriptions
    text(0.3,0.03, 'extinction', col = 'white', cex = ca)
    text(.7,.72, 'carrying capacity', col='black', cex = ca)
    
    #axes
    yMarks = seq(.005, .02, by =.005)    
    xMarks = seq(0, 1200, by=300)
    axis(side = 1, at = xMarks/shift[np], labels = round(xMarks, digits = 0),
         cex.axis=ca, tck=-0.015, padj = -3)
    axis(side = 2, at = yMarks/punish.prob[np], labels = round(yMarks, digits = 4), 
         cex.axis=ca, tck=-0.015, padj = 2.5)
    mtext(side=2, expression('Punishment probability'), 
          cex=cxt, padj = -1.4)
    mtext(side=1, expression(paste("Demand reduction ( elephants / day )") ), 
          cex=cxt, padj=1.2)
    
    #legend  
    image.plot(N.eq, legend.only=TRUE, zlim = c(0,.2),
               col = gray(seq(0,1, by =1/np)), 
               #smallplot= c(.90,.93,0,1), 
               add=TRUE, 
               legend.mar=1.5, axes=F,
               axis.args=list(cex.axis=ca2, tck = -0.15 , line=-.75),
               legend.args=list(text=expression(paste("Equilibrium elephant density,   ind / km"^"2") ),
                                side=4, line=1.15, cex=.7, srt=-90)
    )
    pos.nums = 1.04
  dev.off()

# Fig 2A Effect of slope of demand Curve on legal trade #### 
  # Run 'Elasticity_EquilibriumAnalysis.R' to generate the data for this figure
    pdf( file = 'bSensitivityEFF.pdf', width = 8, height = 5)
      par( mfrow=c(1,1), oma = c(10,4,1,1), mar=c(1,1,0,0) )
      lw = 3; sizeS=2; ca=1.15; ca2=1.1;
      plot(log10(bSeq), prop.SQ, 
           xlim = c(1,5), ylim=c(0,.46),
           xlab='', ylab='', xaxt='n', yaxt='n', yaxs='i',
           type='b', pch =1, lty=1, col=1, lwd=lw, cex= sizeS)
      lines(log10(bSeq), prop.M, type='b', pch=2, lty=2, col =2, lwd=lw, cex= sizeS)
      lines(log10(bSeq), prop.FE, type='b', pch=0, lty=3, col=3, lwd=lw, cex= sizeS)
      legend('topleft', c('status-quo enforcement', 'status-quo enforcement + legal trade', 
                          'legal trade funds extra enforcement'),
             lty = c(1,2,3), col = c( 1,2,3), pch=c(1,2,0), bty = 'n', lwd=lw, cex=ca )
      options(scipen=99) 
      xMarks = 10^(0:5)
      axis(side = 1, at = log10(xMarks), labels = prettyNum(xMarks, big.mark = ","), cex.axis=ca)
      axis(side = 2, at = seq(0,1,.2), cex.axis=ca)
      mtext(side=1, expression('Slope of demand curve'), cex=ca, padj = 2.4)
      mtext(side=2, 'Fraction of scenarios with equilibrium\ndensity at or above current levels', 
            cex=ca, padj = -0.8, outer=TRUE)
      mtext(side=1, 'Demand [elephants/day]', outer=TRUE, cex = ca2, padj=14)
      mtext(side=2, 'Price [USD]', outer=TRUE, cex = ca2, adj=-.8, padj=-3.5)
    dev.off()
    
  
# Fig 2B - Sensitivity both shift and punish ####
    
    
    #load data from 
    
    library(fields)
    col.width = 3.42
    pdf('Sensitivity_Punish_shift.pdf', height = 3, width = col.width)
    
      par(mfrow=c(1,1), oma = c(1.3,1.4,0,2.0), mar=c(.5,.5,.5, 2), las=0 )
      lw = 2; sizeS=2; ca=.6; cxt = .8; ca2=.5;
      #heatmap
        colfunc <- colorRampPalette(c("black", "#FFFDFE"))#lavenderblush"))
        
        image(shift, punish.prob, prop,          
              xlab='', ylab='', 
              xaxt='n', yaxt='n', 
              xaxs="i", yaxs="i",
              col = colfunc(np))#gray(seq(0,1, by =1/np)) )
      
      #axes 
        yMarks = seq(.02, .08, by =.02)    
        xMarks = seq(0, 1, by=.25)
        
        axis(side = 1, at = xMarks, labels = round(xMarks, digits = 2), 
             cex.axis=ca, tck=-0.015, padj = -3)
        axis(side = 2, at = yMarks, labels = round(yMarks, digits = 4), 
             cex.axis=ca, tck=-0.015, padj = 2.5)
        mtext(side=2, expression('Punishment probability'), 
              cex=cxt, padj = -1.4)
        mtext(side=1, expression(paste("Demand reduction ( % max demand )") ), 
              cex=cxt, padj=1.2)
        box()
      
      #plot the line where Equilibrium density = N0
        eps1 = .95; # tolerance for number of scenarios
        eps2 = .8;
        eps3 = .5;
        
        crit.prop.eps1 = apply(prop, 1, 
                               function(x) which( (x- eps1)^2 == min( (x - eps1)^2 ) )[1] ) 
        crit.prop.eps2 = apply(prop, 1, 
                               function(x) which( (x- eps2)^2 == min( (x - eps2)^2 ) )[1] ) 
        crit.prop.eps3 = apply(prop, 1, 
                               function(x) which( (x- eps3)^2 == min( (x - eps3)^2 ) )[1] ) 
        
        lines(shift, punish.prob[crit.prop.eps1], lwd = lw , lty=1)
        lines(shift, punish.prob[crit.prop.eps2], lwd = lw , lty=2)
        lines(shift, punish.prob[crit.prop.eps3], lwd = lw , lty=3)
        
        points(0, s0*pc, col = 'black', cex=1.15, pch =19, xpd = TRUE)
        points(0, s0*pc, col = 'yellow', cex=1.0, pch =19, xpd = TRUE)
        
        text(.45, .052, "95%", cex=ca)
        text(.45, .0215, "80%", cex=ca)
        text(.45, .0096, "50%", cex=ca)

      #legend
        text(1.30, .045, 
             labels = paste("Fraction of scenarios with\nelephant density > current estimate"), 
             xpd = NA, srt = -90, cex = cxt)

        image.plot(prop, legend.only=TRUE, zlim = c(0,1),
                   col = colfunc(np), 
                   add=TRUE, legend.mar=1.5, axes=F,
                   axis.args=list(cex.axis=ca2, tck = -0.1 , line=-.7)
        )
        
    dev.off()
# Fig S2 - Equilibrium elephant density vs punishment probability ####

  pdf( file = 'PunishProb.pdf', width = 8, height = 5)
    par(mfrow=c(1,1), oma = c(5,7,1.1,1.1), mar=c(1,1,1,1) )
    lw = 3; sizeS=2; ca=1.5; cxt = 2.5
    plot(punish.prob, N.eq.punish, 
         xlim = c(0, max(punish.prob)), ylim = c(0,.2),
         xlab='', ylab='', 
         xaxt='n', yaxt='n',
         type='l', pch =1, lty=1,  col=1, lwd=lw, cex= sizeS)
    axis(side = 1, at = seq(0,.02,.005), cex.axis=ca)
    axis(side = 2, at = seq(.05,.2,.05), cex.axis=ca)
    mtext(side=1, expression('Punishment probability'), cex=cxt, padj = 2)
    mtext(side=2, expression(paste("Equilibrium density ( elephants / km"^" 2"," )") ), 
          cex=cxt, padj=-1.5)
    ind.baseline = which( (pc*s0-punish.prob)^2 == min( (pc*s0-punish.prob)^2 ) )
    points(punish.prob[ind.baseline], N.eq.punish[ind.baseline], col = 'red', pch=19, cex=3 )
  dev.off()
  
    