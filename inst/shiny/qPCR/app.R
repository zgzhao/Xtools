library('Xtools')
library('RColorBrewer')
library('HTqPCR')
library('ddCt')
library('reshape2')

shinyApp(
    ##* UI section
    ui = fluidPage(
    ##** header
    responsive=FALSE,
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "http://ppmb.lzu.edu.cn/css/shiny.css")
        ),
    includeHTML("http://ppmb.lzu.edu.cn/includes/shiny.header.html"),
    
    HTML("<div class='content'>"),
    titlePanel("qPCR数据处理工具"),
    ##** Layout
    sidebarLayout(
        ##*** sidebar panel
        sidebarPanel(
            tabsetPanel(
                tabPanel("数据设置", 
                         fileInput("dtfile", h4("上传数据"), 
                                   accept = c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                                   'text/plain', '.csv', '.txt')),
                         radioButtons("Sep", h4("数据分隔符"), choices=c('制表符', '逗号', '空格')), 
                         checkboxGroupInput("Refs", h4("内参基因"), choices=''), 
                         radioButtons("Mock", h4("对照样品（仅相对表达量使用）"), choices='')
                         ),
                tabPanel("图形参数",
                         tabsetPanel(type="pills", 
                                     tabPanel("位置",
                                              sliderInput("legpos", h5("图例位置"), min=1, max = 8, value=1, step = 1),
                                              sliderInput("imageh", h5("图像高度"), min=0.1, max = 3, value=1, step = 0.1),
                                              sliderInput("side1", h5("下边距"), min=0.5, max = 10, value=4, step = 0.5), 
                                              sliderInput("side2", h5("左边距"), min=0.5, max = 10, value=6, step = 0.5), 
                                              sliderInput("side3", h5("上边距"), min=0.5, max = 10, value=2, step = 0.5), 
                                              sliderInput("side4", h5("右边距"), min=0.5, max = 10, value=1, step = 0.5),
                                              sliderInput("mgp1", h5("Y轴刻度标签距离"), min=0, max = 5, value=0.5, step = 0.1),
                                              sliderInput("mgp2", h5("Y轴标题距离"), min=1, max = 10, value=3, step = 0.2),
                                              sliderInput("yext", h5("Y轴扩展"), min=1, max = 3, value=1.05, step = 0.05),
                                              br()
                                              ),
                                     tabPanel("文字", 
                                              sliderInput("cexlabx", h5("X轴标签字体大小"), min=0, max = 5, value=1.5, step = 0.1),
                                              sliderInput("xsrt", h5("X标签文字旋转"), min=0, max = 180, value=0, step = 15),
                                              sliderInput("xmgp1", h5("X标签水平位置"), min=-2, max = 2, value=1, step = 0.01),
                                              sliderInput("xmgp2", h5("X标签垂直位置"), min=-5, max = 5, value=2, step = 0.1),
                                              br(), 
                                              sliderInput("cexlab", h5("Y轴字体大小"), min=0, max = 5, value=1.8, step = 0.1),
                                              br(), 
                                              sliderInput("legcex", h5("图例大小"), min=0.8, max = 3, value=1.5, step = 0.1),
                                              br()
                                              ), 
                                     tabPanel("线条与颜色",
                                              h5('图形颜色'), 
                                              checkboxInput('graycol', '单一灰度颜色', value=FALSE), 
                                              selectInput('cols', '选择色彩系列',
                                                          choices=c('Set1', 'Set2', 'Set3', 'Paired', 'Pastel1',
                                                          'Pastel2', 'Dark2', 'Accent', 'Spectral', 'RdYlGn',
                                                          'RdYlBu', 'RdGy', 'RdBu', 'PuOr', 'PRGn', 'PiYG', 'BrBG',
                                                          'RainBow', 'heat.colors', 'cm.colors', 'terrain.colors', '黑白'),
                                                          selected='Paired', multiple=FALSE, width='120px'),
                                              sliderInput("boxlwd", h5("边框线宽度"), min=0, max = 6, value=1, step = 1),
                                              sliderInput("boxcol", h5("边框线颜色"), min=1, max = 657, value=1, step = 1),
                                              sliderInput("errcex", h5("误差线长度"), min=0, max = 6, value=1, step = 0.1),
                                              sliderInput("errlwd", h5("误差线粗细"), min=0, max = 6, value=1, step = 0.1),
                                              sliderInput("errcol", h5("误差线颜色"), min=1, max = 657, value=1, step = 1),
                                              br()
                                              )
                                     )
                         ),
                tabPanel("柱形图断点", 
                         radioButtons("brkset", h5("柱形图类型"), choices=c('常规', '断裂'), selected='断裂'),
                         conditionalPanel(
                             condition="input.brkset=='断裂'", 
                             radioButtons("brktype", h5("断点类型"), choices=c('normal', 'zigzag')), 
                             sliderInput("brkratio", h5("上下高度比"), min=0.3, max = 3, value=1.4, step = 0.1),
                             sliderInput("brkwidth", h5("断裂宽度"), min=0.1, max = 2, value=1, step = 0.1),
                             sliderInput("brksrt", h5("断裂线旋转"), min=0, max = 180, value=135, step = 5),
                             sliderInput("brksize", h5("断裂线长度"), min=0, max = 5, value=0.6, step = 0.1)
                             ),
                         br()
                         )
                )
            ),
        ##*** main panel
        mainPanel(
            tabsetPanel(
                tabPanel("数据查看",
                         tabsetPanel(type='pills', 
                                     tabPanel("原始数据", tableOutput("datax")),
                                     tabPanel("Ct归一化效果", plotOutput("FigNorm", width="600px", height="600px")), 
                                     tabPanel("样品聚类图",
                                              plotOutput("FigClust1", width="400px", height="400px"), 
                                              plotOutput("FigClust2", width="400px", height="400px"))
                                     )
                         ), 
                tabPanel("相对表达量",
                         tabsetPanel(type='pills', 
                                     tabPanel("表格",
                                              conditionalPanel(
                                                  condition="output.resultDdct", 
                                                  downloadButton('saveDdct', '保存为 CSV 文件'), br(), br()
                                                  ), 
                                              tableOutput("resultDdct")
                                              ),
                                     tabPanel("柱形图",
                                              radioButtons("geneNameRel", "选择基因", choices="", inline=TRUE),
                                              conditionalPanel(
                                                  condition="output.resultDdct",
                                                  imageOutput("FigRel", height='100%')
                                                  ), 
                                              conditionalPanel(
                                                  condition="output.FigRel",
                                                  br(),
                                                  downloadButton('saveFigRel', '下载和保存本图片')
                                                  )
                                              ),
                                     tabPanel("样本间比较",
                                              radioButtons("compType", h4("双样本比较类型"), choices=c("单对照", "自定义"), inline=TRUE),
                                              HTML(
                                                  '<label for="Comps">“单对照”只需输入一个样品名称，“自定义”的输入格式为：T1:CK, T2:CK</label><input id="Comps" type="text" value="" style="width: 95%" />'
                                                  ),
                                              br(), 
                                              conditionalPanel(
                                                  condition="output.resultLimma",
                                                  downloadButton('saveLimma', '保存分析结果')
                                                  ),
                                              br(), 
                                              tableOutput("resultLimma")
                                              )                                     
                                     )
                         ), ## tabPanel 相对表达量
                tabPanel("绝对表达量",
                         tabsetPanel(type='pills', 
                                     tabPanel("表格",
                                              conditionalPanel(
                                                  condition="output.resultAbs", 
                                                  downloadButton('saveAbs', '保存为 CSV 文件'), br(), br()
                                                  ), 
                                              tableOutput("resultAbs")
                                              ),
                                     tabPanel("图形",
                                              checkboxGroupInput("geneNameAbs", "选择图形中包含的基因", choices="", inline=TRUE),
                                              conditionalPanel(
                                                  condition="output.resultAbs",
                                                  imageOutput("FigAbs", height='100%')
                                                  )
                                              )
                                     )
                         ), ## tabPanel 绝对表达量
                tabPanel("方法说明",
                         HTML('
<h4><li>本工具使用<font color="red">Bioconductor</font>高通量qPCR数据分析软件包<font color="red">HTqPCR</font>和<font color="red">ddCt</font>工作</h4>
<h4><li><font color="red">相对表达量</font>使用标准的ddCt方法计算，即内参基因和对照样品都要使用。内参基因在所有样品中的表达量都归一化，且进一步将对照样品的所有基因的表达量都归一化</h4>
<h4><li><font color="red">绝对表达量</font>仅使用内参基因，内参基因在所有样品中的表达量都归一化，对其他基因不做针对某样品的归一化处理</h4>'),
                         hr(),
                         hr(), 
                         HTML("<h4 style='color: #1E90FF'>请从实验室FTP下载 /R语言/data/qPCR.demo.csv 文件查看格式和练习操作</h4>")
                         )
                )
            )
        ), 
    ##* end UI
    HTML("</div>"), 
    includeHTML("http://ppmb.lzu.edu.cn/includes/shiny.footer.html")
    ),
    ##* server section
    server = function(input, output, session) {
        ##** set files
        setData <- reactive({
            f <- input$dtfile
            if(is.null(f)) return(NULL)
            else filex <- f$datapath
            sep <- '\t'
            if (input$Sep=='逗号') sep <- ','
            if (input$Sep=='空格') sep <- ' '
            
            xx <- read.table(filex, header=TRUE, stringsAsFactor=FALSE, sep=sep)
            
            ## 内参基因设置选项
            genes <- unique(xx$gene)
            refs <- grep('^(ACT|UBQ|REF|RF).*', genes, ignore.case=TRUE)
            if (length(refs)==0) refs <- genes[length(genes)]
            else refs <- genes[refs]
            updateCheckboxGroupInput(session, 'Refs', choices=genes, selected=refs)
            
            ## 对照样品设置选项
            treatments <- unique(xx$treatment)
            ck <- grep('^(CK|Control|Mock).*', treatments, ignore.case=TRUE)
            if (length(ck)==0) ck <- treatments[1]
            else ck <- treatments[ck[1]]
            updateRadioButtons(session, 'Mock', choices=treatments, selected=ck)
            if(length(treatments) > 4 | max(nchar(treatments)) > 6) {
                updateSliderInput(session, 'side1', value=8)
                updateSliderInput(session, 'xsrt', value=30)
            }
            return(list(file=filex, sep=sep, samples=treatments))
        })
        ##** output data
        output$datax <- renderTable({
            if(is.null(setData())) return(NULL)
            info <- setData()
            read.table(info$file, header=TRUE, stringsAsFactor=FALSE, sep=info$sep)
        })
        ##** cal results
        calResult <- reactive({
            if(is.null(setData())) return(NULL)
            info <- setData()
            df <- read.table(info$file, header=TRUE, stringsAsFactor=FALSE, sep=info$sep)
            ##*** 参数设置
            infox <- unique(paste(df$treatment, df$sample))
            ## 样品处理类型
            treatments <- gsub("(.*) .*", "\\1", infox)
            ## 样品名称
            samples <- gsub(".* (.*)", "\\1", infox)
            ## 样品数量
            nsample <- length(unique(samples))
            ## 基因（features）名称
            genes <- unique(df$gene)
            ## 基因（features）数量
            ngene <- length(genes)
            ## 内参基因
            refs <- input$Refs
            targets <- genes[!genes %in% refs]
            updateRadioButtons(session, 'figName', choices=targets, selected=targets[1])
            
            data.raw <- readCtData(files=info$file, header=TRUE, sep=info$sep,
                                   n.features=ngene, 
                                   column.info = list(position=1, Ct=2, feature=5, type=6),
                                   n.data=length(samples), samples=samples)
            
            sel <- ! duplicated(paste(treatments, samples))
            pData(data.raw) <- data.frame(sample=samples[sel], treatment=treatments[sel])
            
            d.norm <- normalizeCtData(data.raw, norm = "deltaCt", deltaCt.genes = refs, verbose=FALSE)
            design <- model.matrix(~0 + treatments)
            colnames(design) <- gsub('treatments', '', colnames(design))
            
            ##*** 样本间比较
            vss <- gsub('^\\s*(.*)\\s*$', '\\1', input$Comps)
            if(grepl("^$", vss)) {
                vss <- NULL
            } else {
                if(input$compType=='自定义')
                {
                    if (!grepl('\\w+:\\w+', vss)) vss <- NULL
                    else {
                        xx <- strsplit(vss, '[:, ]+')[[1]]
                        if(length(xx)%%2 > 0) vss <- NULL
                        if (any(! xx %in% treatments)) vss <- NULL
                    }
                    if(!is.null(vss)) {
                        vss <- strsplit(vss, '[, ]+')[[1]]
                        vss <- gsub(':', '-', vss)
                    }
                } else {
                    trs <- unique(df$treatment)
                    if(vss %in% trs){
                        trs <- trs[trs != vss]
                        vss <- paste(trs, vss, sep='-')
                    }
                }
            }
            
            
            result.ht <- NULL
            if(!is.null(vss)) {
                contrasts <- makeContrasts(contrasts=vss, levels=design)
                qDE.limma <- limmaCtData(d.norm, design = design, contrasts = contrasts, ndups = 1,
                                         spacing = 1)
                for(i in 1:length(vss)) 
                    result.ht <- rbind(result.ht, qDE.limma[[i]][targets, c(7, 6, 3, 4, 5)])
                
                GENE <- rep(targets, times=length(vss))
                COMPARE <- rep(gsub('-', '/', vss), each=length(targets))
                result.ht <- cbind(COMPARE, GENE, result.ht)
                rownames(result.ht) <- NULL
            }
            
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
            }
            )
        ##** 输出Limma分析结果
        output$resultLimma <- renderTable({
            if(is.null(setData())) return(NULL)
            calResult()$limma
        })
        
        ##** Limma结果下载UI
        output$saveLimma <- downloadHandler(
            filename = "results.limma.csv",
            content = function(file){
                write.csv(calResult()$limma, file, row.names=FALSE, sep='\t')
            }
            )
        
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
    ##* end server
    }
    )

