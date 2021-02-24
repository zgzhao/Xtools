##** Limma结果下载UI
output$saveLimma <- downloadHandler(
            filename = "results.limma.csv",
            content = function(file){
                write.csv(calResult()$limma, file, row.names=FALSE, sep='\t')
            })
        
        ##** 归一化结果图
        output$FigNorm <- renderPlot({
            if(is.null(setData())) return(NULL)
            datax <- calResult()
            xraw <- exprs(datax$d.raw)
            xnorm <- exprs(datax$d.norm)
            par(mar=c(4, 6, 6, 1))
            cols <- brewer.pal(nrow(xraw), "Spectral")
            plot(xraw, xnorm, pch = 20, cex=2, main = "DeltaCt Normalization", cex.main=2, 
                 col = cols, xlab='Raw data', ylab='Normalized', cex.lab=2, cex.axis=1.5)
        })
        
        ##** 聚类分析图1
        output$FigClust1 <- renderPlot({
            if(is.null(setData())) return(NULL)
            datax <- calResult()
            plotCtCor(datax$d.raw, main = "Raw Ct")
        })
        
        ##** 聚类分析图2
        output$FigClust2 <- renderPlot({
            if(is.null(setData())) return(NULL)
            datax <- calResult()
            plotCtCor(datax$d.norm, main = "ddCt-normed Ct")
        })
        
        ##** 设置相对表达量绘图数据与图形参数
        getRelPars <- reactive({
            if(is.null(setData())) return(NULL)
            
            df <- calResult()$ddct
            glabs <- input$geneNameRel
            df <- df[df$Gene==glabs, c("Sample", "FC", "FC.sd")]
            xlabs <- df$Sample
            datax <- cbind(t(df$FC), t(df$FC.sd))
            
            colnames(datax) <- NULL
            rownames(datax) <- NULL
            n <- ncol(datax)/2
            y.cols <- 1:n
            sdu.cols <- (n + 1):(2 * n)
            sdd.cols <- sdu.cols
            
            w <- 600
            h <- 400 * input$imageh
            
            
            if (input$graycol) cols <- 'gray80'
            else {
                if (input$cols=='heat.colors') cols <- heat.colors(n)
                else if (input$cols=='cm.colors') cols <- cm.colors(n)
                else if (input$cols=='terrain.colors') cols <- terrain.colors(n)
                else if (input$cols=='RainBow') cols <- rainbow(n)
                else if (input$cols=='黑白') cols <- rev(gray(0:n/n))
                else cols <- brewer.pal(n, input$cols)
                cols <- cols[1:n]
            }
            
            xsrt <- input$xsrt
            xadj <- c(input$xmgp1, input$xmgp2)
            
            ## 是否用断点
            min.range <- 10
            max.fold <- 5
            if(input$brkset=='常规') {
                min.range <- 10000000
                max.fold <- 10000000
            }
            
            if(input$boxlwd==0) brklwd <- 1
            else brklwd <- input$boxlwd
            
            list(datax=datax, y.cols=y.cols, sdu.cols=sdu.cols, sdd.cols=sdd.cols, glabs=glabs, xlabs=xlabs, 
                 w=w, h=h, cols=cols, xsrt=xsrt, xadj=xadj, brklwd=brklwd, 
                 min.range=min.range, max.fold=max.fold)
        })

        ##** cal results
        calResult <- reactive({
            if(is.null(setData())) return(NULL)
            info <- setData()
            df <- read.table(info$file, header=TRUE, stringsAsFactor=FALSE, sep=info$sep)
            genes <- unique(df$gene)
            refs <- input$Refs
            targets <- genes[!genes %in% refs]
            updateRadioButtons(session, 'figName', choices=targets, selected=targets[1])
            
            data.raw   <- readDataHTqPCR(files=info$file, header=TRUE, sep=info$sep)
            result.ht <- calHTqPCR(dt, refs, norm.method="deltaCt", verbose=FALSE)
            
            ##*** ddCt 计算表达量和误差
            coreData <- data.frame(matrix(nrow = nrow(df), ncol = 4))
            colnames(coreData) <- c("Sample", "Detector", "Ct", "Platename")
            
            ddct.raw <- new("InputFrame")
            ddct.raw@coreData <- coreData
            ddct.raw@coreData$Sample <- df$treatment
            ddct.raw@coreData$Ct <- df$Ct
            ddct.raw@coreData$Detector <- df$gene
            ddct.raw@coreData$Platename <- rep(1, nrow(df))
            ddct.raw@files <- info$file
            mock <- input$Mock
            
            ##*** 相对表达量
            result.dd <- ddCtExpression(ddct.raw, calibrationSample = mock, housekeepingGene = refs)
            result.ddct <- elist(result.dd)[, c(2, 1, 3, 4)]
            colnames(result.ddct) <- c('Sample', 'Gene', 'FC', 'FC.sd')
            result.ddct <- result.ddct[result.ddct$Gene %in% targets, ]
            
            ##*** 绝对表达量
            result.abs <- elist(result.dd)[, c(2, 1, 7, 8)]
            colnames(result.abs) <- c('Sample', 'Gene', 'expr', 'expr.sd')
            x1 <- result.abs$expr
            x2 <- result.abs$expr.sd
            result.abs$expr <- 2^(-x1) * 1000
            result.abs$expr.sd <- abs(result.abs$expr - 2^(-(x1 + x2)) * 1000)
            result.abs <- result.abs[result.abs$Gene %in% targets, ]
            
            ##*** 更新UI
            updateRadioButtons(session, 'geneNameRel', choices=targets, selected=targets[1])            
            updateCheckboxGroupInput(session, 'geneNameAbs', choices=targets, selected=targets)
            list(d.raw=data.raw, d.norm=d.norm, limma=result.ht, ddct=result.ddct, absExpr=result.abs)
        })
        
        ##** 输出相对表达量
        output$resultDdct <- renderTable({
            if(is.null(setData())) return(NULL)
            calResult()$ddct
        })
        
        ##** 相对表达量下载UI
        output$saveDdct <- downloadHandler(
            filename = "results.ddct.csv",
            content = function(file){
                write.csv(calResult()$ddct, file, row.names=FALSE, sep='\t')
            }
            )
        
        ##** 输出绝对表达量
        output$resultAbs <- renderTable({
            if(is.null(setData())) return(NULL)
            calResult()$absExpr
        })
        ##** 绝对表达量下载UI
        output$saveAbs <- downloadHandler(
            filename = "results.exprs.csv",
            content = function(file){
                write.csv(calResult()$absExpr, file, row.names=FALSE, sep='\t')
            })
        
        ##** 绘制相对表达量图片
        output$FigRel <- renderImage(
        {
            xx <- getRelPars()
            attach(xx)
            
            outfile <- tempfile(fileext='.png')
            png(outfile, width=w, height=h)
            if(is.null(xx)){
                dev.off()
                return(NULL)
            }
            
            par(mar=c(input$side1,input$side2,input$side3,input$side4),
                mgp=c(input$mgp2, input$mgp1, 0),
                cex.axis=input$cexlab * 0.8)
            
            ## 绘图
            xx <- gap.barplot(datax, y.cols=y.cols, sdu.cols=sdu.cols, sdd.cols=sdd.cols, yext=input$yext, 
                              col=cols, min.range=min.range, max.fold=max.fold,
                              box.lwd=input$boxlwd, box.col=input$boxcol,
                              error.cex=input$errcex, error.lwd=input$errlwd, error.col=input$errcol, 
                              brk.type=input$brktype, ratio=input$brkratio, gap.width=input$brkwidth,
                              brk.srt = input$brksrt, brk.size = input$brksize, brk.lwd = brklwd,
                              brk.col=input$boxcol)
            
            ## 标题
            title(ylab="Relative Expression", cex.lab=input$cexlab)
            
            axis(1, at=xx, labels=NA)
            text(xx, 0, labels=xlabs, srt=xsrt, xpd=TRUE, adj=xadj, cex=input$cexlabx)
            legpos <- switch(input$legpos, 'topleft', 'top', 'topright', 'left', 'right', 'bottomleft', 'bottom', 'bottomright')
            legend(legpos, legend=glabs, box.col=NA, bg="transparent", cex=input$legcex)
            
            dev.off()
            
            list(src = outfile,
                 contentType = 'image/png',
                 width = w,
                 height = h,
                 alt = "No figures or plotting error!")
        }, deleteFile = TRUE)
        
        ##** 相对表达量图片下载UI
        output$saveFigRel <- downloadHandler(
            filename = function(){
                paste("figure", input$geneNameRel, "tiff", sep=".")
            }, 
            content = function(file) {
                xx <- getRelPars()
                attach(xx)
                tiff(file, width=w, height=h)
                
                if(!is.null(xx)) {
                    par(mar=c(input$side1,input$side2,input$side3,input$side4),
                        mgp=c(input$mgp2, input$mgp1, 0),
                        cex.axis=input$cexlab * 0.8)
                    ## 绘图
                    xx <- gap.barplot(datax, y.cols=y.cols, sdu.cols=sdu.cols, sdd.cols=sdd.cols, yext=input$yext, 
                                      col=cols, min.range=min.range, max.fold=max.fold,
                                      box.lwd=input$boxlwd, box.col=input$boxcol,
                                      error.cex=input$errcex, error.lwd=input$errlwd, error.col=input$errcol, 
                                      brk.type=input$brktype, ratio=input$brkratio, gap.width=input$brkwidth,
                                      brk.srt = input$brksrt, brk.size = input$brksize, brk.lwd = brklwd,
                                      brk.col=input$boxcol)
                    
                    ## 标题
                    title(ylab="Relative Expression", cex.lab=input$cexlab)
                    
                    axis(1, at=xx, labels=NA)
                    text(xx, 0, labels=xlabs, srt=xsrt, xpd=TRUE, adj=xadj, cex=input$cexlabx)
                    
                    legpos <- switch(input$legpos, 'topleft', 'top', 'topright', 'left', 'right', 'bottomleft', 'bottom', 'bottomright')
                    legend(legpos, legend=glabs, box.col=NA, bg="transparent", cex=input$legcex)
                    
                }            
                dev.off()
            }
            )
        
        ##** 绝对表达量绘图数据与图形参数
        getAbsPars <- reactive({
            if(is.null(setData())) return(NULL)
            df <- calResult()$absExpr
            df <- df[df$Gene %in% input$geneNameAbs, ]
            meanx <- df[, c("Gene", "Sample", "expr")]
            sdx <- df[, c("Gene", "Sample", "expr.sd")]
            
            meanx <- dcast(meanx, Gene~Sample)
            sdx <- dcast(sdx, Gene~Sample)
            ## meanx <- dcast(meanx, Sample~Gene)
            ## sdx <- dcast(sdx, Sample~Gene)
            
            xlabs <- meanx[, 1]
            meanx <- meanx[, -1]
            sdx <- sdx[, -1]
            glabs <- colnames(meanx)
            
            df <- cbind(meanx, sdx)
            m <- ncol(df)/2
            y.cols <- 1:m
            sdu.cols <- (m + 1):(2 * m)
            sdd.cols <- sdu.cols
            n <- nrow(df)
            
            w <- 600
            h <- 400 * input$imageh
            
            nn <- m
            if(nn < 4) nn <- 4
            if (input$graycol) cols <- 'gray80'
            else {
                if (input$cols=='heat.colors') cols <- heat.colors(nn)
                else if (input$cols=='cm.colors') cols <- cm.colors(nn)
                else if (input$cols=='terrain.colors') cols <- terrain.colors(nn)
                else if (input$cols=='RainBow') cols <- rainbow(nn)
                else if (input$cols=='黑白') cols <- rev(gray(0:nn/nn))
                else cols <- brewer.pal(nn, input$cols)
                cols <- cols[1:m]
            }
            
            xsrt <- input$xsrt
            xadj <- c(input$xmgp1, input$xmgp2)
            
            ## 是否用断点
            min.range <- 10
            max.fold <- 5
            if(input$brkset=='常规') {
                min.range <- 10000000
                max.fold <- 10000000
            }
            
            if(input$boxlwd==0) brklwd <- 1
            else brklwd <- input$boxlwd
            
            list(datax=df, y.cols=y.cols, sdu.cols=sdu.cols, sdd.cols=sdd.cols, glabs=glabs, xlabs=xlabs, 
                 w=w, h=h, cols=cols, xsrt=xsrt, xadj=xadj, brklwd=brklwd, 
                 min.range=min.range, max.fold=max.fold)
        })
        
        ##** 绘制绝对表达量图片
        output$FigAbs <- renderImage(
        {
            ##*** 设置参数
            xx <- getAbsPars()
            attach(xx)
            
            outfile <- tempfile(fileext='.png')
            png(outfile, width=w, height=h)
            if(is.null(xx)){
                dev.off()
                return(NULL)
            }
            
            par(mar=c(input$side1,input$side2,input$side3,input$side4),
                mgp=c(input$mgp2, input$mgp1, 0),
                cex.axis=input$cexlab * 0.8)
            
            ##*** 绘图
            xx <- gap.barplot(datax, y.cols=y.cols, sdu.cols=sdu.cols, sdd.cols=sdd.cols, yext=input$yext, 
                              col=cols, min.range=min.range, max.fold=max.fold,
                              box.lwd=input$boxlwd, box.col=input$boxcol,
                              error.cex=input$errcex, error.lwd=input$errlwd, error.col=input$errcol, 
                              brk.type=input$brktype, ratio=input$brkratio, gap.width=input$brkwidth,
                              brk.srt = input$brksrt, brk.size = input$brksize, brk.lwd = brklwd,
                              brk.col=input$boxcol)
            
            ##*** 标题
            title(ylab="Relative Expression", cex.lab=input$cexlab)
            
            if(nrow(xx) > 1) xx <- colMeans(xx)
            axis(1, at=xx, labels=NA)
            text(xx, 0, labels=xlabs, srt=xsrt, xpd=TRUE, adj=xadj, cex=input$cexlabx)
            legpos <- switch(input$legpos, 'topleft', 'top', 'topright', 'left', 'right', 'bottomleft', 'bottom', 'bottomright')
            legend(legpos, legend=glabs, fill=cols, box.col=NA, bg="transparent", cex=input$legcex)
            dev.off()
            
            list(src = outfile,
                 contentType = 'image/png',
                 width = w,
                 height = h,
                 alt = "No figures or plotting error!")
        }, deleteFile = TRUE)


library(HTqPCR)
library(ddCt)
library(zQPCR)
ff <- "/home/zhao/A.Labwork/0.BR.n.Cold.Acclimation/qPCR/mxp.files/BZR1-targets.2014.7-8/Ct.results20140820.csv"
dt <- read.qPCRtable(ff, sep=",")
ddct.raw <- dt$data.ddct
mock <- "Col-0@NA@a"
refs <- "ACT2"
result.dd <- ddCtExpression(ddct.raw, calibrationSample = mock, housekeepingGene = refs)
result.ddct <-
str(elist(result.dd))
[, c(2, 1, 3, 4)]
colnames(result.ddct) <- c('Sample', 'Gene', 'FC', 'FC.sd')
            result.ddct <- result.ddct[result.ddct$Gene %in% targets, ]
