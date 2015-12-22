require('Xtools')
shinyApp(
    ##* UI section ==========================================
    ui=fluidPage(
        ##** header
        responsive=FALSE,
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "http://ppmb.lzu.edu.cn/css/shiny.css")
        ),
        includeHTML("http://ppmb.lzu.edu.cn/includes/shiny.header.html"),
        
        HTML("<div class='content'>"),
        h1("低温半致死温度计算"),
        ##** Layout
        sidebarLayout(
            ##*** sidebarPanel
            sidebarPanel(
                tabsetPanel(type='pills', 
                            tabPanel("数据选项", 
                                     checkboxInput('demo', HTML('<span style="color:#FF0000">使用演示数据</span>')),
                                     conditionalPanel(
                                         condition = '!input.demo', 
                                         fileInput("dtfile", "上传电导率结果文件：",
                                                   accept = c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                                                       'text/plain', '.csv', '.txt')),
                                         radioButtons("sep", h5("数据分隔符"), choices=c('制表符', '逗号', '空格'), inline=TRUE)
                                     ),
                                     sliderInput("Temps", h5("温度范围"), min=-25, max =10, value=c(-18, 2), step = 1),
                                     checkboxGroupInput("hideCols", h5("选择图中不显示的样品"), choices='', inline=TRUE)
                            ),
                            tabPanel("图形设置", 
                                     selectInput('cols', h5('图形颜色'),
                                                 choices=c('Set1', 'Set2', 'Set3', 'Paired', 'Pastel1',
                                                     'Pastel2', 'Dark2', 'Accent', 'Spectral', 'RdYlGn',
                                                     'RdYlBu', 'RdGy', 'RdBu', 'PuOr', 'PRGn', 'PiYG', 'BrBG',
                                                     'RainBow', 'heat.colors', 'cm.colors', 'terrain.colors', 'Gray'),
                                                 selected='Dark2', multiple=FALSE, width='120px'),
                                     radioButtons("seplot", h5("置信区间"), choices=c('显示', '隐藏'), inline=TRUE),
                                     radioButtons("pointplot", h5("数据点"), choices=c('显示', '隐藏'), inline=TRUE),
                                     br(), 
                                     fluidRow(
                                         div(numericInput("cexlab", "坐标轴字体大小", value=1, min=0.2, max = 5, step = 0.1),
                                             class="col-md-6"), 
                                         div(numericInput("cexpoint", "数据点大小", value=1, min=0.2, max = 5, step = 0.1),
                                             class="col-md-6"), 
                                         div(numericInput("cexleg", "图例大小", min=0.8, max = 3, value=1.2, step = 0.1),
                                             class="col-md-6"),                                          
                                         div(numericInput("lwd", "线宽", min=0.1, max = 10, value=1, step = 0.1), 
                                             class="col-md-6")
                                     ), 
                                     sliderInput('alpha', '阴影透明度', min=0.1, max=0.5, value=0.2, step=0.05)
                            )
                )
            ),
            ##*** mainPanel
            mainPanel(
                ##**** section1
                HTML('<input id="status" type="number" value="1" hidden/>'),
                conditionalPanel(
                    condition = 'input.status==0',
                    div(style="font-size:16px; margin-top:40px;", 
                        HTML('<h4 style="color:#FF4500">数据读取错误，或没有上传文件。<br/>',
                             '如果你已上传数据，请选择合适的数据分隔符。</h4>', 
                             '<ul style="color:#006400">',
                             '<li>数据的第一列必需为温度且列名称为“TMP”</li>',
                             '<li>缺失、无效值请用“NA”表示</li>', 
                             '<li>文件需保存为UTF-8编码，<br/>如果你不知道如何设置，',
                             '<span style="color:#FF0000">请使用演示数据保存为录入模板。</span></li>',
                             '</ul>'
                        )
                    )
                ),
                ##**** section2
                conditionalPanel(
                    condition='input.status==1', 
                    tabsetPanel(
                        ##***** tabpan1
                        tabPanel("原始数据", 
                                 conditionalPanel(
                                     condition = 'input.demo',
                                     HTML('<h4>以下为演示数据，你可以'),
                                     downloadLink('template', '保存为数据录入模板'),
                                     HTML('</h4>')
                                 ),
                                 tableOutput("datax"),
                                 style="width:100%; max-height:600px;"
                        ), 
                        ##***** tabpan2
                        tabPanel("图",
                                 imageOutput("lt50Curve", width="600px", height="400px"),
                                 br(),
                                 imageOutput("bars",  width="600px", height="400px")
                        ), 
                        ##***** tabpan3
                        tabPanel("表",
                                 h4("CI95和CI99分别表示95%和99%置信区间，L为区间下限，U为上限。"),
                                 p("95%置信区间不重叠可认为有显著差异，即：p<0.05"),
                                 p("99%置信区间不重叠可认为有极显著差异，即：p<0.01"),
                                 downloadButton('saveTable', '保存为 CSV 文件'), br(), br(),
                                 tableOutput("lt50s")
                        ),
                        ##***** tabpan4
                        tabPanel("显著性检验",
                                 br(), 
                                 checkboxInput('uplt50', '上传LT50数据'), 
                                 tableOutput('sigresult'), 
                                 conditionalPanel(
                                     condition = 'input.uplt50', 
                                     fileInput("ltfile", "",
                                               accept = c('text/csv', 'text/comma-separated-values',
                                                   'text/tab-separated-values',
                                                   'text/plain', '.csv', '.txt')
                                     )
                                 )
                        )
                    )
                )
            )
        ), 
        HTML("</div>"), 
        includeHTML("http://ppmb.lzu.edu.cn/includes/shiny.footer.html")
    ), 
    ##* server section
    server = function(input, output, session) {
        ##** getData
        getData <- reactive({
            if(input$demo) {
                data(lt50.demo)
                df <- lt50.demo
            } else {
                filex <- input$dtfile
                if(is.null(filex)) df <- NULL
                else {
                    filex <- filex$datapath
                    sep <- '\t'
                    if(input$sep=='空格') sep <- ' '
                    if(input$sep=='逗号') sep <- ','
                    df <- read.table(filex, header=TRUE, sep=sep)
                    if(ncol(df) < 2 | is.null(df$TMP)) df <- NULL
                }
            }
            df
        })

        observe({
            df <- getData()
            if(is.null(df)) {
                namex <- ''
                status <- 0
            } else {
                status <- 1
                namex <- colnames(df)[-1]
            }
            updateCheckboxGroupInput(session, "hideCols", choices = namex)
            updateNumericInput(session, "status", value = status)
        })
        
        ##** getSel
        getSel <- reactive({
            df <- getData()
            namex <- colnames(df)
            sel <- ! namex %in% input$hideCols
            x <- 1:length(namex)
            x <- x[sel]
            x[ - 1] - 1
        })
        
        ##** View original data
        output$datax <- renderTable({
            getData()
        })        
        ##** cal LT50s
        calLT50s <- reactive({
            df <- getData()
            if(is.null(df)) return(NULL)
            rt <- c(input$Temps[1], input$Temps[2])            
            getLT50s(df, rt, all.info=TRUE)
        })

        ##** 显著性分析
        data4sig <- reactive({
            if(!input$uplt50) {
                if(is.null(getData)) dt <- NULL
                else dt <- t(calLT50s()$vals)
            } else {
                filex <- input$ltfile
                if(is.null(filex)) dt <- NULL
                else {
                    filex <- filex$datapath
                    seps <- c(',', '\t', ' ')
                    dt <- NULL
                    for(sep in seps) {
                        try(df <- read.table(filex, header=TRUE, sep=sep))
                        if(ncol(df) > 5) dt <- as.matrix(df[, 2:6])
                    }
                }
            }
            dt
        })

        calsig <- reactive({
            res <- data4sig()
            if(is.null(res)) return(NULL)
            m <- nrow(res)
            sigtable <- NULL
            for(i in 1:m) {
                ir1 <- res[i, 2:3]
                ir2 <- res[i, 4:5]
                aaa <- NULL
                for(j in 1:m) {
                    ir3 <- res[j, 2:3]
                    ir4 <- res[j, 4:5]
                    xx <- ''
                    if(!isOverlap(ir2, ir4)) xx <- '**'
                    else if(!isOverlap(ir1, ir3)) xx <- '*'
                    aaa <- c(aaa, xx)
                }
                sigtable <- rbind(sigtable, aaa)
            }
            rownames(sigtable) <- rownames(res)
            colnames(sigtable) <- rownames(res)
            sigtable
        })
        output$sigresult <- renderTable({
            calsig()
        })
        
        ##** output LT50s table
        output$lt50s <- renderTable({
            if(is.null(getData)) return(NULL)
            xx <- calLT50s()$vals
            t(xx)
        })
        
        ##** 分析结果表格保存
        output$saveTable <- downloadHandler(
            filename = "result.LT50s.csv",
            content = function(file){
                if(is.null(getData)) xx <- NULL
                else xx <- calLT50s()$vals
                write.csv(t(xx), file, sep='\t')
            }
        )
        
        ##** 绘图与显示结果
        observe({
            if(is.null(getData)) return(NULL)
            ##*** 数据处理
            info <- calLT50s()
            m <- length(info$preds)
            if(is.null(m)) m <- 1
            if(input$cols=="Gray") cols <- rep('gray40', m)
            else if (input$cols=='heat.colors') cols <- heat.colors(m)
            else if (input$cols=='cm.colors') cols <- cm.colors(m)
            else if (input$cols=='terrain.colors') cols <- terrain.colors(m)
            else if (input$cols=='RainBow') cols <- rainbow(m)
            else cols <- brewer.pal(m, input$cols)
            
            sels <- getSel()
            output$lt50Curve <- renderImage({
                image1 <- tempfile(fileext='.jpg')
                jpeg(filename=image1, width=600, height=400, quality=100)
                se.plot <- input$seplot=='显示'
                point.plot <- input$pointplot=='显示'
                xlim <- c(input$Temps[1], input$Temps[2])
                par(mar=c(4, 8, 1, 0.5))                
                plot(info, show.cols=sels,
                     point.plot=point.plot, se.plot=se.plot, cex.point=input$cexpoint,
                     col=cols, xlim=xlim, lwd=input$lwd, 
                     cex.axis=input$cexlab * 0.8, se.alpha=input$alpha, 
                     cex.lab=input$cexlab, cex.legend =input$cexleg)
                dev.off()
                
                list(src=image1)
            })
            
            ##*** 柱形图
            output$bars <- renderImage({
                image2 <- tempfile(fileext='.jpg')
                jpeg(image2, width=600, height=400, quality=100)
                lt50s <- round(info$vals, 2)
                lt50s <- lt50s
                meanx <- unlist(lt50s[1, sels])
                sd1 <- unlist(lt50s[2, sels])
                sd2 <- unlist(lt50s[3, sels])
                
                par(mar=c(1, 8, 6, 0.5))
                xx <- barplot(meanx, ylim=c(min(lt50s[1:3, ])*1.1, 0), col=cols[sels], 
                              axes=F, ylab=expression(paste('LT50 (', degree, 'C)')), cex.lab=input$cexlab)
                axis(2, cex.axis=input$cexlab)
                axis(3, at=xx, labels=rep("", length(xx)))
                text(xx, 0, labels=colnames(lt50s)[sels], cex=input$cexlab * 1.1, adj=c(0.8, -1), srt=-15, xpd=TRUE)
                text(xx, lt50s[1, sels], labels=lt50s[1, sels], adj=c(1, 1.2))
                errorbar(x=xx, y=meanx, sd.lwr=meanx-sd1, sd.upr=sd2-meanx)
                box()
                dev.off()
                list(src=image2)
            })
        })

        output$template <- downloadHandler(
            filename = 'template.EL.4.LT50.csv',
            content = function(file) {
                data(lt50.demo)
                write.csv(lt50.demo, file, row.names=FALSE)
            }
        )
        ##*** server end
    })
