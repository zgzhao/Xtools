library('HTqPCR')
library('ddCt')
library('reshape2')
library('RColorBrewer')
library('Xtools')

shinyApp(
    ## UI section
    ui = fluidPage(
        ## header
        includeCSS("www/style.css"),
        includeHTML("www/header.html"),

        HTML("<div class='content'>"),
        titlePanel("RT-qPCR data analysis"),
        ## Layout
        sidebarLayout(
            ## sidebar panel
            sidebarPanel(
                downloadLink('template', 'Download Ct template file here.'),
                fileInput("dtfile", h4("Upload file"), accept = c('text/csv', 'text/comma-separated-values', 'text/tab-separated-values',
                                                                  'text/plain', '.csv', '.txt')),
                radioButtons("Sep", h4("Data seperator (in file)"), choices=c('comma', 'tab', 'space'), inline=TRUE),
                checkboxGroupInput("Refs", h4("Reference gene"), choices=''),
                radioButtons("Mock", h4("Control sample（Relative expression only）"), choices='')
            ),
            ## main panel
            mainPanel(
                tabsetPanel(
                    tabPanel("DataView",
                             tableOutput("datax")
                             ),
                    tabPanel("ddCt",
                             tabsetPanel(type='pills',
                                         tabPanel("Expression (relative)",
                                                  conditionalPanel(
                                                      condition="output.ddctRelative",
                                                      downloadButton('saveRelative', 'download result'), br(), br()
                                                  ), 
                                                  tableOutput("ddctRelative")
                                                  ),
                                         tabPanel("Expression (obsolute)",
                                                  conditionalPanel(
                                                      condition="output.ddctAbsolute",
                                                      downloadButton('saveAbsolute', 'download result'), br(), br()
                                                  ), 
                                                  tableOutput("ddctAbsolute")
                                                  )
                                         )
                             ),
                    tabPanel("HTqPCR",
                             tabsetPanel(type='pills',
                                         tabPanel("Ct normalization", plotOutput("FigNorm", width="600px", height="600px")),
                                         tabPanel("Sample comparison",
                                                  radioButtons("compType", h4("two sample comparison"), choices=c("Single control", "Free input"), inline=TRUE),
                                                  HTML('<label for="Comps">Single control: input a sample name (copy from left panel).</br>',
                                                       'Free input: input comparisons like this：T1:CK, T2:CK</label><input id="Comps" type="text" value="" style="width: 95%" />'
                                                       ),
                                                  br(),
                                                  conditionalPanel(
                                                      condition="output.resultLimma",
                                                      downloadButton('saveLimma', 'download result')
                                                  ),
                                                  br(),
                                                  tableOutput("resultLimma")
                                                  )
                                         )
                             ),
                    tabPanel("Data format", includeHTML("www/data.format.html"))
                )
            )
        ), 
        ##* end UI
        HTML("</div>"),
        includeHTML("www/footer.html")
    ),
    ##* server section
    server = function(input, output, session) {
    ##** set files
    setData <- reactive({
        f <- input$dtfile
        if(is.null(f)) return(NULL)
        else filex <- f$datapath
        sep <- '\t'
        if (input$Sep=='comma') sep <- ','
        if (input$Sep=='space') sep <- ' '

        xx <- read.table(filex, header=TRUE, stringsAsFactor=FALSE, sep=sep)
        ## 内参基因设置选项
        genes <- unique(xx$gene)
        refs <- grep('^(ACT|UBQ|REF|RF).*', genes, ignore.case=TRUE)
        if (length(refs)==0) refs <- genes[length(genes)]
        else refs <- genes[refs]
        updateCheckboxGroupInput(session, 'Refs', choices=genes, selected=refs)

        ## 设置qPCR数据
        htset <- NULL
        ddset <- NULL
        xdata <- read.qPCRtable(filex, sep=sep)
        if(!is.null(xdata)) {
            htset <- xdata$data.ht
            ddset <- xdata$data.ddct
            samples <- unique(as.character(pData(htset)$Sample))
            updateRadioButtons(session, 'Mock', choices=samples, selected=samples[1])
            if(length(samples) > 4 | max(nchar(samples)) > 6) {
                updateSliderInput(session, 'side1', value=8)
                updateSliderInput(session, 'xsrt', value=30)
            }
        }

        return(list(data=xx, data.ht=htset, data.ddct=ddset))
    })

    ##*** 原始数据输出
    output$datax <- renderTable({
        info <- setData()
        if(is.null(info)) NULL else setData()$data
    })

    ##*** HTqPCR：数据分析
    resultHTqPCR <- reactive({
        info <- setData()
        if(is.null(info)) return(NULL)
        refs <- input$Refs
        ctype <- if( input$compType == "Single control") 1 else 2
        comps <- input$Comps
        res <- calHTqPCR(info$data.ht, refs, comps, ctype)
        return(res)
    })

    ##*** Limma (HTqPCR)
    calLimma <- reactive({
        info <- resultHTqPCR()
        if(is.null(info)) return(NULL)
        if(is.null(info$result)) return(NULL)
        cts <- info$contrast
        cns <- names(cts)
        ##data.frame(cns)
        res <- NULL
        for (i in 1:length(cts)) {
            aa <- info$result[[i]]
            aa$Label <- rep(cns[i], nrow(aa))
            aa <- aa[, c("Label", "genes", "FC", "t.test", "p.value", "adj.p.value")]
            res <- rbind(res, aa)
        }
        res <- res[order(res$Label, res$genes), ]
        res[! res$genes %in% input$Refs, ]
    })
    ##** HTqPCR：输出Limma结果
    output$resultLimma <- renderTable({
        dt <- calLimma()
        if(is.null(dt)) return(NULL)
        dt
    })

    ##** Limma结果下载UI
    output$saveLimma <- downloadHandler(
        filename = "results.expr.limma.csv",
        content = function(file){
        write.csv(calLimma(), file, row.names=FALSE, sep='\t')
    })

    ##** 归一化结果图
    output$FigNorm <- renderPlot({
        info1 <- setData()
        info2 <- resultHTqPCR()
        if(is.null(info1) || is.null(info2)) return(NULL)
        plot(exprs(info1$data.ht), exprs(info2$d.norm), pch = 20, main = "DeltaCt normalisation",
             col = rep(brewer.pal(6, "Spectral")))
    })

    ##** ddCt计算
    resultDdct <- reactive({
        info <- setData()
        if(is.null(info)) return(NULL)
        ddct.raw <- info$data.ddct
        if(is.null(ddct.raw)) return(NULL)
        mock <- input$Mock
        refs <- input$Refs

        result.dd <- ddCtExpression(ddct.raw, calibrationSample = mock, housekeepingGene = refs)
        result.dd <- elist(result.dd)
        expr.rel <- result.dd[, c("Sample", "Detector", "exprs", "level.err")]
        colnames(expr.rel) <- c('Sample', 'Gene', 'FC', 'FC.sd')
        expr.rel <- expr.rel[! expr.rel$Gene %in% refs, ]

        expr.abs <- result.dd[, c("Sample", "Detector", "dCt", "dCt.error")]
        colnames(expr.abs) <- c('Sample', 'Gene', 'expr', 'expr.sd')
        aa <- expr.abs$expr
        bb <- expr.abs$expr.sd
        expr.abs$expr <- 2^(- aa) * 1000
        expr.abs$expr.sd <- abs(expr.abs$expr - 2^(-(aa + bb)) * 1000)
        expr.abs <- expr.abs[! expr.abs$Gene %in% refs, ]

        ##updateRadioButtons(session, 'geneNameRel', choices=targets, selected=targets[1])
        ##updateCheckboxGroupInput(session, 'geneNameAbs', choices=targets, selected=targets)
        list(relative=expr.rel, absolute=expr.abs)
    })

    ##** 输出相对表达量
    output$ddctRelative <- renderTable({
        info <- resultDdct()
        if(is.null(info)) return(NULL)
        info$relative
    })

    ##** 输出绝对表达量
    output$ddctAbsolute <- renderTable({
        info <- resultDdct()
        if(is.null(info)) return(NULL)
        info$absolute
    })

    ##** 相对表达量下载UI
    output$saveRelative <- downloadHandler(
        filename = "results.expr.rel.csv",
        content = function(file){
        write.csv(resultDdct()$relative, file, row.names=FALSE, sep='\t')
    })
    ##** 绝对表达量下载UI
    output$saveAbsolute <- downloadHandler(
        filename = "results.expr.abs.csv",
        content = function(file){
        write.csv(resultDdct()$absolute, file, row.names=FALSE, sep='\t')
    })
    output$template <- downloadHandler(
        filename = 'template.Ct.data.csv',
        content = function(file) {
        data(Ct.demo)
        write.csv(Ct.demo, file, row.names=FALSE)
    }
    )
})
