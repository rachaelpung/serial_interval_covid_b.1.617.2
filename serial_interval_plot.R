library(grid)

# colour
CONFIG = list() 
# CONFIG$cols = c('#00468B','#ED0000','#42B540','#0099B4','#925E9F') # BMJ
# CONFIG$cols = c('#00468BFF','#FDAF91FF','#42B540FF','#0099B4FF','#925E9FFF') # Lancet Dark
# CONFIG$cols = c('#9CD6E9', '#E0A6A1', '#C9E0B4', '#039FC6', '#9FADD4') # Lancet Pastel - blue, pink, green, turquiose, purple
CONFIG$cols = c('#E64B35FF', '#4DBBD5FF', '#00A087FF', '#3C5488FF', '#F39B7FFF', 
                '#8491B4FF', '#91D1C2FF', '#DC0000FF', '#7E6148FF', '#B09C85FF') # Nature


lightup = function(c, alpha)
{
  z=col2rgb(c)/255
  return(rgb(z[1],z[2],z[3],alpha))  # 0.125 for col3
}
CONFIG$colsLight1 = c();for(i in seq(CONFIG$cols))CONFIG$colsLight1[i] = lightup(CONFIG$cols[i], alpha = 0.6)
CONFIG$colsLight2 = c();for(i in seq(CONFIG$cols))CONFIG$colsLight2[i] = lightup(CONFIG$cols[i], alpha = 0.4)
CONFIG$colsLight3 = c();for(i in seq(CONFIG$cols))CONFIG$colsLight3[i] = lightup(CONFIG$cols[i], alpha = 0.25)

rm(lightup,i)

# fig 1 main text
paneller=function(row=1,column=1)
{

  if(column==1){xlm=c(-3.5,15.5); ylm=c(0,.35)} 
  if(column==2){xlm=c(-3,15); ylm=c(0,.35)}
  if(column==3){xlm=c(-3,15); ylm=c(0,1)}
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=1))
  pushViewport(plotViewport(c(3,4,1,1),xscale=xlm,yscale=ylm))
  grid.rect()
  grid.text(c('A','B','C')[column],unit(-3,'lines'),unit(1,'npc')+unit(0,'lines'))
  
  
  if(column ==1){
    
    data=casePair.VOC
    data=data[,.N, by= .(VARIANT_WGS,SERIAL_INT)]
    
    data_617=data[VARIANT_WGS == 'B.1.617.2']
    xv_1 = data_617[,SERIAL_INT]
    yv_1 = data_617[,N]/sum(data_617[,N])
    
    for(i in 1:nrow(data_617)){
      grid.polygon(xv_1[i]+c(-0.5,-0.5,0.5,0.5),
                   c(0,  yv_1[i],  yv_1[i], 0), default.units = 'native', gp=gpar(col='white',fill=CONFIG$colsLight2[3]))
    }
    

  }
  
  if(column == 2){
    
    data=distSerialInt[, .(SERIAL_INTERVAL, PMF_PRE_LD_MEAN, PMF_PRE_LD_LOWER_CI, PMF_PRE_LD_UPPER_CI)]
    
    rm = which(data$SERIAL_INTERVAL %in% c(-3,10:15))
    xv=data[[1]]
    yv_mean=data[[2]]
    yv_lower=data[[3]]
    yv_upper=data[[4]]
    
    grid.polygon(c(xv,rev(xv)),c(yv_lower,rev(yv_upper)),default.units = 'native',gp=gpar(col=NA,fill=CONFIG$colsLight2[4]))
    grid.lines(xv,yv_mean,default.units = 'native',gp=gpar(col=CONFIG$cols[4]))
    grid.points(xv[-rm],yv_mean[-rm],default.units = 'native',gp=gpar(col=CONFIG$cols[4],cex=0.5),pch=16)
    
  }
  
  if(column == 3){
    
    data=distSerialInt[, .(SERIAL_INTERVAL, CDF_PRE_LD_MEAN, CDF_PRE_LD_LOWER_CI, CDF_PRE_LD_UPPER_CI)]
    
    rm = which(data$SERIAL_INTERVAL %in% c(-3,9:15))
    xv=data[[1]]
    yv_mean=data[[2]]
    yv_lower=data[[3]]
    yv_upper=data[[4]]
    
    grid.polygon(c(xv,rev(xv)),c(yv_lower,rev(yv_upper)),default.units = 'native',gp=gpar(col=NA,fill=CONFIG$colsLight2[4]))
    grid.lines(xv,yv_mean,default.units = 'native',gp=gpar(col=CONFIG$cols[4]))
    grid.points(xv[-rm],yv_mean[-rm],default.units = 'native',gp=gpar(col=CONFIG$cols[4],cex=0.5),pch=16)
    
    
    data = casePair.VOC[,.N, by=.(SERIAL_INT)][order(SERIAL_INT)]

    data = rbindlist(list(data, data.table(SERIAL_INT=-3:15, N=NA)), use.names = T, fill= T)
    data = data[order(SERIAL_INT,N,na.last=F)]
    data[is.na(N),N:=0]
    data[,CUM_N:=cumsum(N)]
    data[,CUM_PROP:=CUM_N/sum(N)]
    
    xv = data[[1]]
    yv = data[[4]]
    rm = which(data$SERIAL_INT %in% c(-3,11:15))
  
    
    grid.lines(xv,yv,default.units = 'native',gp=gpar(col=CONFIG$cols[3]))
    grid.points(xv[-rm],yv[-rm],default.units = 'native',gp=gpar(col=CONFIG$cols[3],cex=0.3),pch=4)
    
  }
  
  
  
  
  grid.xaxis(at=seq(0,15,5),label=seq(0,15,5))
  if(column==3)grid.yaxis(at=seq(0,1,0.2),label=seq(0,100,20))
  if(column%in%c(1,2))grid.yaxis(at=seq(0,0.35,0.05),label=seq(0,35,5))

  grid.text('Time (d)',y=unit(-3,'lines'))
  if(column%in%c(1,2))grid.text('Probability (%)',x=unit(-3,'lines'),rot=90)
  if(column==3)grid.text('Cumulative probability (%)',x=unit(-3,'lines'),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('serial interval.png',height=7.5,width=22.5,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(1,0,0,0)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=3)))
paneller(1,1)
paneller(1,2)
paneller(1,3)
popViewport()
popViewport()
dev.off()

rm(paneller)


# supp fig 1
paneller=function(row=1,column=1)
{
  xlm=c(-0.5,20)
  ylm=c(-5,20)
  
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=1))
  pushViewport(plotViewport(c(3,4,1,1),xscale=xlm,yscale=ylm))
  grid.rect()
  grid.text(c('A','B')[column],unit(-3,'lines'),unit(1,'npc')+unit(0,'lines'))
  
  
  if(row==1 & column==1) data=casePair.VOC
  if(row==1 & column==2) data=casePair.preLD
  
  
  if(column==1){
  
    n = data[VARIANT_WGS == 'B.1.617.2',.N]
    xv_1 = data[VARIANT_WGS == 'B.1.617.2', ONSET_TO_ISOLATION_INFECTOR] + runif(n, min=-0.5, max=0.5)
    yv_1 = data[VARIANT_WGS == 'B.1.617.2', SERIAL_INT] + runif(n, min=-0.5, max=0.5)
    grid.points(xv_1,yv_1,default.units = 'native',gp=gpar(col=CONFIG$cols[3],cex=0.3),pch=4)
    
    # n = data[VARIANT_WGS == 'Unknown',.N]
    # xv_2 = data[VARIANT_WGS == 'Unknown', ONSET_TO_ISOLATION_INFECTOR] + runif(n, min=-0.5, max=0.5)
    # yv_2 = data[VARIANT_WGS == 'Unknown', SERIAL_INT] + runif(n, min=-0.5, max=0.5)
    # grid.points(xv_2,yv_2,default.units = 'native',gp=gpar(col=CONFIG$colsLight2[3],cex=0.3),pch=4)
    
  }
  
  
  
  if(column==2){
    
    n = data[,.N]
    xv = data[, ONSET_TO_ISOLATION_INFECTOR] + runif(n, min=-0.5, max=0.5)
    yv = data[, SERIAL_INT] + runif(n, min=-0.5, max=0.5)
    grid.points(xv,yv,default.units = 'native',gp=gpar(col=CONFIG$cols[4],cex=0.3),pch=16)
    
  }
  
  
  
  
  grid.xaxis(at=seq(0,20,5),label=seq(0,20,5))
  grid.yaxis(at=seq(-5,20,5),label=seq(-5,20,5))
  grid.text('Time onset to isolation in index case (d)',y=unit(-3,'lines'))
  if(column==1)grid.text('Serial interval (d)',x=unit(-3,'lines'),rot=90)
  
  popViewport()
  popViewport()
  
}

png('onset_isolation.png',height=7.5,width=15,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(1,0,0,0)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))
paneller(1,1)
paneller(1,2)
popViewport()
popViewport()
dev.off()

rm(paneller)


# supp fig 2 barplot
paneller=function(row=1,column=1)
{
  xlm=c(-3.5,17.5)
  ylm=c(0,12) 
  
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=1))
  pushViewport(plotViewport(c(3,4,1,1),xscale=xlm,yscale=ylm))
  grid.rect()
  grid.text(c('A','B')[column],unit(-3,'lines'),unit(1,'npc')+unit(0,'lines'))
  
  
  if(column ==1){
    
    data=casePair.VOC
    data=data[,.N, by= .(VARIANT_WGS,SERIAL_INT)]
    
    data_617=data[VARIANT_WGS == 'B.1.617.2']
    xv_1 = data_617[,SERIAL_INT]
    yv_1 = data_617[,N]
    
    for(i in 1:nrow(data_617)){
      grid.polygon(xv_1[i]+c(-0.5,-0.5,0.5,0.5),
                   c(0,  yv_1[i],  yv_1[i], 0), default.units = 'native', gp=gpar(col='white',fill=CONFIG$colsLight2[3]))
      
    }
   
    # data_unknown = data[VARIANT_WGS == 'Unknown']
    # data_unknown[data_617, LWR := i.N, on=c(SERIAL_INT='SERIAL_INT')]
    # data_unknown[is.na(LWR), LWR:=0]
    # data_unknown[,UPP:=LWR+N]
    # xv_2 = data_unknown[,SERIAL_INT]
    # yv_2lwr = data_unknown[,LWR]
    # yv_2upp = data_unknown[,UPP]
    # 
    # 
    # for(i in 1:nrow(data_unknown)){
    #   grid.polygon(xv_2[i]+c(-0.5,-0.5,0.5,0.5),
    #                c(yv_2lwr[i], yv_2upp[i],  yv_2upp[i], yv_2lwr[i]), default.units = 'native', gp=gpar(col='white',fill=CONFIG$colsLight2[3]))
    #   
    # }
    
  }
  
  
  if(column != 1){
    
    data = casePair.preLD
    data=data[,.N, by= .(SERIAL_INT)]
    
    xv_1 = data[,SERIAL_INT]
    yv_1 = data[,N]
    
    for(i in 1:nrow(data)){
      grid.polygon(xv_1[i]+c(-0.5,-0.5,0.5,0.5),
                   c(0,  yv_1[i],  yv_1[i], 0), default.units = 'native', gp=gpar(col='white',fill=CONFIG$colsLight2[4]))
      
    }
    
  }
  
  
  
  grid.xaxis(at=seq(0,15,5),label=seq(0,15,5))
  grid.yaxis(at=seq(0,12,2),label=seq(0,12,2))
  grid.text('Time (d)',y=unit(-3,'lines'))
  if(column==1)grid.text('Cases (n)',x=unit(-3,'lines'),rot=90)
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  popViewport()
  popViewport()
  
}


png('serial_interval_non_adjust.png',height=7.5,width=15,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(1,0,0,0)))
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2)))
paneller(1,1)
paneller(1,2)
popViewport()
popViewport()
dev.off()

rm(paneller)


