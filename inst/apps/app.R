library(shiny)
library(shinythemes)
library(plotly)
library(shinyalert)
library(dplyr)
library(FBN)
library(purrr)
library(edgeR)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(vctrs)

# ui ----
# Where all the user interface decisions are
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  navbarPage("Normalisation and visualisation of TraDIS data", id = "tab",
             #tags$script(HTML("var header = $('.navbar> .container-fluid');
            #           header.append('<div style=\"float:right;color:white\"><h5>Last updated: 8 Nov 2022</h4></div>');
            #           console.log(header)")),
             # tab 1 -------
             ## normalisation ----
             tabPanel("Normalisation", value = "1",
                      sidebarPanel(width = 3,
                                   conditionalPanel(condition="input.tab==1 & input.normtab==2", fileInput("uploadrc", "Upload your traDIS read files here", buttonLabel = "Upload...", multiple = TRUE)),
                                   conditionalPanel(condition="input.tab==1 & input.normtab==1", fileInput("uploadfc", "Upload your traDIS output file(s) here", buttonLabel = "Upload...", multiple = TRUE)),
                                   conditionalPanel(condition="input.tab==1 & input.normtab==1", uiOutput("datasetsUIfc")),
                                   conditionalPanel(condition="input.tab==1 && input.normtab==2", uiOutput("controlUIrc")),
                                   conditionalPanel(condition="input.tab==1 && input.normtab==2", actionButton("controlselect", label = "Select")),
                                   br(),
                                   conditionalPanel(condition="input.tab==1 && input.normtab==2", numericInput("window", "Sliding window size", value = 500, step = 50, min = 100, max = 1000)),
             ),
                      mainPanel(
                        tabsetPanel(id = "normtab",
                                    tabPanel("Diagnostics", value = 1,
                                             br(),
                                             h4("Upload your Bio::TraDIS output files to determine whether chromosomal bias is affecting your data"),
                                             p("If the overall trend of your fold changes does not match the red line, your data needs normalising."),
                                             plotOutput("diag_fc")),
                                    tabPanel("Correction", value = 2,
                                             h4("Upload your Bio::TraDIS read files to correct the chromosomal bias affecting your data"),
                                             p("This requires two control files and two condition files."),
                                             plotOutput("corrected_plot"),
                                             actionButton("download_attempt", label = "Write to file"))))),
             
             # overarching tab 2 -----
             ## analysis ----
             tabPanel("Visualisation", value = 2,
                      sidebarPanel(width = 3,
                                   conditionalPanel(condition="input.tab==2", fileInput("upload", "Upload your traDIS output file(s) here", buttonLabel = "Upload...", multiple = TRUE)),
                                   conditionalPanel(condition="input.tab==2 && input.subtab==3",fileInput("eggnog_upload", "Upload your eggnog_mapper file here", buttonLabel = "Upload...", multiple = TRUE)),
                                   conditionalPanel(condition="input.tab==2 & input.subtab==1 | input.subtab==2", uiOutput("datasetsUI")),
                                   numericInput("fc", "Absolute value fold change cut-off", value = 1, step = 0.5),
                                   numericInput("sig", "Significance cut-off", value = 0.05, step = 0.01, min = 0, max = 1),
                                   ),
                      mainPanel(
                        tabsetPanel(id = "subtab",
                                    tabPanel("Volcano plot", value = 1,
                                             plotlyOutput("plotlyVolcano"),
                                             DT::dataTableOutput("sigdf")),
                                    tabPanel("MA plot", value = 2,
                                             br(),
                                             p("The MA plot shows average log counts versus fold change. The lower the counts, the higher the variability in fold change.
                                                The points should tighten as the average counts increase."),
                                             plotOutput("maplot")),
                                    tabPanel("Pathways", value = 3,
                                             h4(tagList("This tab requires an eggnog mapper file. Please go to the ", 
                                                        a("eggnog mapper site", href = "http://eggnog-mapper.embl.de/"),
                                                        " to perform this analysis.")),
                                             downloadButton("download_plot", "Download plot"),
                                             plotOutput("cog"),
                                             br(),
                                             downloadButton("download_cogtable", "Download table"),
                                             DT::dataTableOutput("cogtable")),
                                    tabPanel("Filtering", value = 14,
                                             h4("Use this tab to filter and/or download your data"),
                                             downloadButton("filter_download"),
                                             DT::dataTableOutput("filter"))))),
             
  # overarching tab 3 ----
  # information
  tabPanel("App feedback",
           h4("Reporting bugs"),
           p("Email Geri at geraldine.sullivan@hdr.mq.edu.au if there are 
                errors in the app."),
           h4("General feedback"),
           p("Feedback for functionality or additional capabilities
                is also welcome."))))

server <- function(input, output, session) {
  
  output$datasetsUIfc <- renderUI({
    selectizeInput("datasetsnorm", "Select dataset for visualising:",
                   choices = gsub(".csv", "", input$uploadfc$name))
  })
  
  output$datasetsUI <- renderUI({
    selectizeInput("datasets", "Select dataset for visualising:",
                   choices = gsub(".csv", "", input$upload$name))
  })
  
  output$controlUIrc <- renderUI({
    selectizeInput("controlrc", "Select which condition is your control:",
                   choices = unique(gsub("_[0-9].tradis.gene.insert.sites.csv", "", input$uploadrc$name)))
  })
  
  output$diag_fc <- renderPlot({
      req(input$uploadfc)
      num <- grep(value = FALSE, pattern = input$datasetsnorm, x = input$uploadfc$name)
      data <- read.csv(input$uploadfc[[num, "datapath"]])
      data$obs <- 1:nrow(data)
      ggplot(data, aes(x = obs, y = logFC)) +
        geom_point(cex = 0.5) +
        geom_hline(yintercept = 0, col = "red") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size = 16)) +
        labs(x = "Locus", y = "Log2 Fold Change", title = paste("Fold change by locus scatterplot - ", 
                                                                gsub(pattern = ".csv", "", input$uploadfc$name[[num]])))
  })
  
  observeEvent(input$subtab, {
    if (input$subtab == "3") {
      updateNumericInput(inputId = "fc", value = 0, session = session)
      updateNumericInput(inputId = "sig", value = 1, session = session)
    } else if (input$subtab == "2") {
      updateNumericInput(inputId = "fc", value = 0, session = session)
    } else {
      updateNumericInput(inputId = "fc", value = 1, session = session)
      updateNumericInput(inputId = "sig", value = 0.05, session = session)
    }
  })
  
  correction <- reactive({
    req(input$uploadrc)
    req(input$controlselect)
    conds <- unique(gsub("_[0-9].tradis.gene.insert.sites.csv", "", input$uploadrc$name))
    ctrl <- as.character(input$controlrc)
    condition <<- subset(conds, conds != ctrl)
    myfiles <- purrr::map(input$uploadrc$datapath, read.delim) %>%
      purrr::set_names(input$uploadrc$name)
    showNotification(paste("Joining read count files together"), duration = 8, type = "message")
    joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag") #join together by locus tag
    filenames <- input$uploadrc$name %>% 
      gsub(pattern = ".tradis_gene_insert_sites.csv", replacement = "") #make file names to rename columns later
    rc <- joined %>% select(contains(c("locus_tag", "read_count"))) #extract only read counts
    locus_tags <- rc$locus_tag # handy for later
    names <- character(0) # rename columns
    for (i in 1:length(filenames)){
      readcount <- paste0(filenames[i])
      names <- append(names, readcount)
      rm(readcount)
    }
    colnames(rc)[2:ncol(rc)] <- names #rename columns
    rc$observation <- 1:nrow(rc) #make observation column, 1-4488
    rc <- rc[c(1,ncol(rc),2:(ncol(rc)-1))] #move observation column to front
    norm_counts <- as.data.frame(1:nrow(rc)) #make data frame for normalised counts and offset
    rownames(norm_counts) <- rc$locus_tag #naming rows - locus tags
    showNotification(paste("Normalising read counts"), duration = 8, type = "message")
    for (i in 3:ncol(rc)){ #normalise read counts using ratios of averages, window size of 20
      calc <- rc[,c(1, 2, i)]
      calc[(nrow(calc)+1):(nrow(calc)+1000),] <- calc[1:1000,]
      calc$keep <- 1:nrow(calc)
      #calc$pred <- medianFilter(inputData = calc[,3], windowSize = 500)
      calc$pred <- medianFilter(inputData = calc[,3], windowSize = input$window)
      calc <- calc[!calc$keep > nrow(rc),]
      calc$ratio <- calc$pred/mean(calc$pred)
      calc$norm <- as.integer(round(calc[,3]/calc$ratio))
      norm_counts[,i-1] <- calc$norm #each column is one condition replicate
    }
    offset <- (log(norm_counts[,2:5] + 0.01) - log(rc[,3:6] + 0.01))
    eff.lib <- calcNormFactors(norm_counts[,2:5]) * colSums(norm_counts[,2:5])
    offset <- sweep(offset, 2, log(eff.lib), "-")
    colnames(norm_counts)[2:5] <- colnames(offset) <- colnames(rc)[3:ncol(rc)] # rename columns
    norm_counts <- norm_counts[,-1] #remove observation column
    rownames(rc) <- rc$locus_tag #name readcount columns
    rc <- rc[,-c(1:2)] #remove descriptive columns, just keeping row names
    #data.table(offset)
    showNotification(paste("Removing low-count genes"), duration = 8, type = "message")
    norm_counts <- norm_counts[apply(apply(rc, 1, ">", 10), 2, any),]
    offset <- offset[apply(apply(rc, 1, ">", 10), 2, any),]
    locus_info <- joined[apply(apply(rc, 1, ">", 10), 2, any),c(1:2,11)]
    locus_info <<- locus_info
    rc <- rc[apply(apply(rc, 1, ">", 10), 2, any),]
    conds_edgeR <- unique(gsub("_[0-9]$", replacement = "", x = filenames)) 
    group <- as.factor(rep(c(conds_edgeR), each = 2))
    group <- relevel(group, ref=ctrl)
    #group <- relevel(group, ref="MH")
    design <- model.matrix(~0+group) # create design
    d <- DGEList(counts = rc, group=group)
    plotMDS.DGEList(d, labels=group)
    d <- calcNormFactors(d)
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d)
    de.tgw <- exactTest(d,pair=c(ctrl, condition))
    #de.tgw <- exactTest(d,pair=c("MH", "Cip"))
    tags_before <- data.frame(de.tgw$table)
    tags_before$col <- ifelse(tags_before$PValue < 0.05, "yes", "no")
    tags_before <<- tags_before
    contrast <- makeContrasts(contrasts = paste0("group", condition, " - group", ctrl), levels = design)
    y <- DGEList(counts=rc, group=group)
    y <- scaleOffset(y, -as.matrix(offset))
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design, robust=TRUE)
    lrt <- glmLRT(fit, contrast=contrast)
    tags_after <- lrt$table
    tags_after$col <- ifelse(tags_after$PValue < 0.05, "yes", "no")
    tags_after <<- tags_after
    rm(myfiles)
    cond <<- condition
  })
  
  output$corrected_plot <- renderPlot({
    req(input$uploadrc)
    req(input$controlselect)
    correction()
    before <- ggplot(tags_before, aes(x = 1:nrow(tags_before), y = logFC, col = col)) +
      geom_point(size = 0.5) +
      theme_classic() +
      theme(text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, size = 14)) +
      labs(title = paste0("Uncorrected fold change plot by locus - ", cond),
           x = "Locus", y = "Log2 fold change") +
      scale_color_manual(values = c("no" = "black", "yes" = "red"), guide = "none")
    after <- ggplot(tags_after, aes(x = 1:nrow(tags_after), y = logFC, col = col)) +
      geom_point(size = 0.5) +
      theme_classic() +
      theme(text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, size = 14)) +
      labs(title = paste0("Fold change plot by locus - ", cond),
           x = "Locus", y = NULL, col = "Significant?") +
      scale_color_manual(values = c("no" = "black", "yes" = "red"))
    before+after
  })
  
  observeEvent(input$download_attempt, {
    correction()
    tags <- tags_after[,c(1:4)]
    diff <- cbind(locus_info, tags)
    write.table(diff,file=paste0(cond, "_corrected.csv"), 
                append=FALSE, quote=TRUE, sep=",", row.names=FALSE, 
                col.names=c("locus_tag","gene_name","function","logFC",
                            "logCPM","PValue","q.value"))
    showNotification(paste0(condition, "_corrected.csv has been saved to ", getwd()))
    shinyalert(title = "Success",
               text = paste0(condition, "_corrected.csv has been saved to ", getwd()))
  })
  
  dataproc <- reactive({
    req(input$upload)
    num <- grep(value = FALSE, pattern = input$datasets, x = input$upload$name)
    read.csv(input$upload[[num, "datapath"]])
  })
  
  output$plotlyVolcano <- renderPlotly({
    #data <- read.csv(file = "files/logfcs/Cip_uncorrected.csv")
    data <- dataproc()
    data$group <- "NotSignificant"
    data[which(data$q.value < input$sig & abs(data$logFC) > input$fc ),"group"] <- "Significant&FoldChange"
    #data[which(data$q.value < 0.05 & abs(data$logFC) > 1 ),"group"] <- "Significant&FoldChange"
    
    p <- ggplot(data, aes_string(x = "logFC", y = "-log10(q.value)", text = "gene_name")) +
      geom_point(aes(color = group)) +
      theme_classic() +
      #geom_vline(xintercept = c(-1, 1), linetype = "dotted", size = 0.3) +
      geom_vline(xintercept = c(input$fc, -input$fc), linetype = "dotted", size = 0.3) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", size = 0.3) +
      scale_color_manual(values = c("NotSignificant"="gray70", "Significant&FoldChange"="red")) +
      labs(x = "Log2 Fold Change", y = "-Log2 Q Value") +
      theme(legend.position = "none",
            text = element_text(size = 16))
    ggplotly(p, tooltip = c("gene_name"))
  })
  
  output$sigdf <- DT::renderDataTable({
    data <- dataproc()
    filt <- subset(data, abs(data$logFC) > input$fc & data$q.value < input$sig)
    filt <- filt %>% dplyr::select("locus_tag", "gene_name", "function.", "logFC", "q.value")
    filt <- filt[order(filt$logFC, decreasing = FALSE),]
    filt$logFC <- round(filt$logFC, digits = 2)
    filt$q.value <- round(as.numeric(filt$q.value), digits = 3)
    DT::datatable(filt, rownames = FALSE, options = list(pageLength = 25)) %>%
      DT::formatStyle('logFC', backgroundColor = DT::styleInterval(c(0), c('#7ec384', '#6aa5a9')))
  })
  
  output$maplot <- renderPlot({
    #data <- read.csv(file = "files/logfcs/TpSu.csv")
    data <- dataproc()
    ggplot(data, aes(x = logCPM, y = logFC)) +
      geom_point(size = 0.5) +
      theme_classic() +
      labs(title = "MA plot") +
      theme(text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5))
  })
  
  
}

shinyApp(ui = ui, server = server)


