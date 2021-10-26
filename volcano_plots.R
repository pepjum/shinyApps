library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)
library(ggplot2)
library(ggrepel)
#library(readr)
library(dplyr)
options(ggrepel.max.overlaps = Inf)
ui<-dashboardPage(
     skin="yellow",
     dashboardHeader(title = "Programa de Inmunología e Inmunoterapia ",
     titleWidth = 450),
     dashboardSidebar(

        tags$head(
            tags$style(HTML("
                        @import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
                        body {
                              background-color: white;
                              color: black;
                            }
                            .shiny-input-container {
                     color: #474747;
                                }"))

            ),
        nav_links <- tags$ul(
        tags$li(
            tags$a(href = "/welcome/", "Welcome"),
        ),
        tags$li(
            tags$a(href = "/preprocess/", "Pre-process RNA-Seq PE data"),
        ),
        tags$li(
            tags$a(href = "/ideal/", "Perform RNA-Seq DEg analysis"),
        ),
        tags$li(
            tags$a(href = "/gating/", "Perform gating FC "),
        ),
        tags$li(
            tags$a(href = "/prognosis/", "Prognosis Marker Analysis"),
        ),
        tags$li(
            tags$a(href = "/tsne/", "tSNE analysis"),
        ),
        tags$li(
            tags$a(href = "/volcano/", "Volcano Plot"),
        ),



        )  #links navegadores
     ),
     dashboardBody(
            tags$head(
            tags$style(HTML("
                         @import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
                      body {
                              background-color: white;
                              color: black;
                            }
                            .shiny-input-container {
                     color: #474747;
                                }"))

            ),

              box(
                width = 12,
                title = "Volcano Plot", status = "danger", solidHeader = TRUE,
                h5("Here you can achive a Volcano Plot of your data. Make sure that you select the correct organism before running."),
                h5(strong("WARNING: ")," This tool needs as input the output matrix of the pipeline for study DEg with edgeR package. If it is not possible, enter a matrix with at least the ENSG column, the logFC and the P.Value named on this way."),
                h5("You can download an example of the format by clicking in the following button:"),
                fluidRow(
                    column(1, offset = 7,
                    downloadButton("Example_matrixVolcano", "Example edgeR matrix",style = "background-color:#FFFFFF; color:#009846; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"))),
                    tags$hr(),
                    column(1, offset = 7,
                    downloadButton("Example_minimal", "Minimal example matrix",style = "background-color:#FFFFFF; color:#009846; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"))
                    )
                ),


                
                fileInput("upload2", NULL, width = 500, buttonLabel = "Upload your matrix data...", multiple = FALSE, accept=".txt"),
#                fileInput("upload2", NULL, width = 500, buttonLabel = "Upload your matrix tumor data...", multiple = FALSE, accept=".txt"),

                #sliderInput("perplexity", "Perplexity:", 30, 50, 40, width = 120),
                fluidRow(column(2, tags$br(  
                   selectInput("ColorPalette_Volcano", "Color Groups:",
                   choices = c("Blue-Red","Green-Red","Black-Red"), width="500px"))),
                   column(4, tags$br(
                   selectInput("Species_Volcano", "Human or Mouse data?:",
                   choices = c("Human","Mouse"), width="200px"))),
                   column(6, tags$br( 
                   sliderInput("Nogenes_Volcano", "Plot the name of top N most differentially expressed genes:", 0,50, 10, width="500px"),
                   )),
                ),
                textInput("Title_Volcano","Enter the main title of the plot"),  
                 

  # sidebarPanel
                actionButton("goButton_Volcano", "Draw VolcanoPlot!", style = "background-color:#FA1414; color:#FFFFFF", icon= icon("fas fa-pencil") , lib = "font-awesome" ),

                downloadButton("download_Volcano", "Download VolcanoPlot", class="btn-success"),

            #    actionButton("plot_surv", "Run Analysis!", style = "background-color:#FA1414; color:#FFFFFF"),

            ),
            verbatimTextOutput("path"),
            verbatimTextOutput("files"),
            verbatimTextOutput("geneExpression"),

            # br(),
            # br(),
            # fluidRow(
            #       splitLayout(style = "border: 1px solid silver:", cellWidths = c(200,200,200), 
            #           plotOutput(outputId = "distPlot",width = "30%"),
            #           plotOutput(outputId = "distPlot2",width = "30%"),
            #           plotOutput(outputId = "distPlot3",width = "30%")
            #                           )
            #               )
            plotOutput(outputId = "VolcanoPlot",width ="500px", height = "500px"),




     ) # dashboardbody
)


server = function(input, output, session){
    options(shiny.maxRequestSize=90000000*1024^2)
    session$onSessionEnded(stopApp)


     Matrix_dataVolcano <- reactive({
         req(input$upload2)
     })   # DATA SELECTED VolcanoPlot

     ColV <- reactive({
         req(input$ColorPalette_Volcano)
     })   # ColorVolcano

     SpecV <- reactive({
         req(input$Species_Volcano)
     })   # SpeciesVolcano

     NgenesV <- reactive({
         req(input$Nogenes_Volcano)
     })   # NgenesVolcano

     MainTitleVolcano <- reactive({
         req(input$Title_Volcano)
     })   # TitleVolcano

    

    matrix_example<-read.table("/media/inmuno/data_projects/example_Volcano.txt", header=T)
    minimal_example<-read.table("/media/inmuno/data_projects/minimal_example.txt", header=T)
   # matrix_example$ENSEMBL<-row.names(matrix_example)

     

    output$Example_matrixVolcano <- downloadHandler(
    filename = function() {
      "Example_matrix.txt"
    },
    content = function(file) {
     	    
      write.table(matrix_example, file, col.names=T, row.names=F, quote=F, sep="\t")
    }) # page 1

    output$Example_minimal <- downloadHandler(
    filename = function() {
      "minimal_example.txt"
    },
    content = function(file) {
     	    
      write.table(minimal_example, file, col.names=T, row.names=F, quote=F, sep="\t")
    }) # page 1


     observeEvent(input$goButton_Volcano, {

       validate(
           need(Matrix_dataVolcano() !="", "Please upload a file!")
       )

       withProgress(
           message = "Creating Volcano Plot...",
           detail = "Created by JGonzalezGom",
           value = 0,
           {

              incProgress(0.1)

              dataVolcano<-read.table(Matrix_dataVolcano()$datapath, header=T, sep="\t") 
             if(ncol(dataVolcano)> 5){

               
              #datamerged<-rbind(dataset,dataset2)
                    FCth<-1
                    Bth<-0
            

                    if(SpecV()=="Mouse"){
                    
                    library(org.Mm.eg.db)

                        ensids<-rownames(dataVolcano)
                        ensids<-paste(lapply(strsplit(paste(ensids),"\\."),"[",1))
                        ensids<-unique(ensids)

                        symbols <- AnnotationDbi::select(
                        org.Mm.eg.db,
                        columns = c("SYMBOL"),
                        keys = ensids,
                        keytype = "ENSEMBL"
                        )


                        dataVolcano$ENSG<-paste(lapply(strsplit(rownames(dataVolcano),"\\."),"[",1))
                        DataVolcano_annot<-merge(dataVolcano,symbols, by.x="ENSG", by.y="ENSEMBL")
                    
                    }else if(SpecV()=="Human"){

                        library(org.Hs.eg.db)
                        
                        ensids<-rownames(dataVolcano)
                        ensids<-paste(lapply(strsplit(paste(ensids),"\\."),"[",1))
                        ensids<-unique(ensids)

                        symbols <- AnnotationDbi::select(
                        org.Hs.eg.db,
                        columns = c("SYMBOL"),
                        keys = ensids,
                        keytype = "ENSEMBL"
                        )


                        dataVolcano$ENSG<-paste(lapply(strsplit(rownames(dataVolcano),"\\."),"[",1))
                        DataVolcano_annot<-merge(dataVolcano,symbols, by.x="ENSG", by.y="ENSEMBL")
                    
                    }

                        dataset<-DataVolcano_annot
                        dataset<-unique(dataset)

                        dataset <- dataset %>% 
                        mutate(
                            Expression = case_when(logFC >= FCth & B > Bth ~ "Up-regulated",
                                                logFC <= -FCth & B > Bth ~ "Down-regulated",
                                                TRUE ~ "Unchanged")
                            )



                        y_int<-(-1)*log10(max(dataset$P.Value[dataset$B> Bth]))

                        if(ColV()=="Blue-Red"){
                            colvalors<-c("#00008b", "gray50", "firebrick3")
                        }else if(ColV()=="Green-Red"){
                            colvalors<-c("#7CAE00", "gray50", "firebrick3")
                        }else if(ColV()=="Black-Red"){
                            colvalors<-c("#000000", "gray50", "firebrick3")
                        }    

                        p2 <- ggplot(dataset, aes(logFC, -log(P.Value,10))) +
                        geom_point(aes(color = Expression), size = 2/5) +
                        xlab(expression("log"[2]*"FC")) + 
                        ylab(expression("-log"[10]*"P.Value")) +
                        scale_color_manual(values = colvalors) +
                        guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
                        geom_vline(xintercept=FCth) + geom_vline(xintercept=-(FCth))+
                        geom_hline(yintercept=y_int)+
                        ggtitle(MainTitleVolcano())+ theme(plot.title=element_text( size=8, face='bold'))


                            if(NgenesV()>0){
                            top <- NgenesV()

                            top_genes <- bind_rows(
                            dataset %>% 
                                filter(Expression == 'Up-regulated') %>% 
                                arrange(P.Value, desc(abs(logFC))) %>% 
                                head(top),
                            dataset %>% 
                                filter(Expression == 'Down-regulated') %>% 
                                arrange(P.Value, desc(abs(logFC))) %>% 
                                head(top)
                            )

                            p3<-p2+ geom_label_repel(data = top_genes,
                                            mapping = aes(logFC, -log(P.Value,10), label = SYMBOL),
                                            size = 2)
                            }
            }else if(ncol(dataVolcano)==3){   ### objects not provided by edgeR

                 #dataVolcano$ENSG<-paste(lapply(strsplit(paste(dataVolcano$ENSG),"\\."),"[",1))
                 rownames(dataVolcano)<-paste(dataVolcano$ENSG)
                 FCth<-1

                    if(SpecV()=="Mouse"){
                    
                    library(org.Mm.eg.db)

                        ensids<-rownames(dataVolcano)
                        ensids<-paste(lapply(strsplit(paste(ensids),"\\."),"[",1))
                        ensids<-unique(ensids)

                        symbols <- AnnotationDbi::select(
                        org.Mm.eg.db,
                        columns = c("SYMBOL"),
                        keys = ensids,
                        keytype = "ENSEMBL"
                        )


                        dataVolcano$ENSG<-paste(ensids)
                        DataVolcano_annot<-merge(dataVolcano,symbols, by.x="ENSG", by.y="ENSEMBL")
                    
                    }else if(SpecV()=="Human"){

                        library(org.Hs.eg.db)
                        
                        ensids<-rownames(dataVolcano)
                        ensids<-paste(lapply(strsplit(paste(ensids),"\\."),"[",1))
                        ensids<-unique(ensids)

                        symbols <- AnnotationDbi::select(
                        org.Hs.eg.db,
                        columns = c("SYMBOL"),
                        keys = ensids,
                        keytype = "ENSEMBL"
                        )


                        dataVolcano$ENSG<-paste(ensids)
                        DataVolcano_annot<-merge(dataVolcano,symbols, by.x="ENSG", by.y="ENSEMBL")
                    
                    }

                dataset<-DataVolcano_annot
                dataset<-unique(dataset)

                dataset <- dataset %>% 
                mutate(
                    Expression = case_when(logFC >= FCth  ~ "Up-regulated",
                                        logFC <= -FCth   ~ "Down-regulated",
                                        TRUE ~ "Unchanged")
                    )



                #y_int<-(-1)*log10(max(dataset$adj.P.Value[dataset$B> Bth]))

                if(ColV()=="Blue-Red"){
                    colvalors<-c("#00008b", "gray50", "firebrick3")
                }else if(ColV()=="Green-Red"){
                    colvalors<-c("#7CAE00", "gray50", "firebrick3")
                }else if(ColV()=="Black-Red"){
                    colvalors<-c("#000000", "gray50", "firebrick3")
                }    

                p2 <- ggplot(dataset, aes(logFC, -log(P.Value,10))) +
                geom_point(aes(color = Expression), size = 2/5) +
                xlab(expression("log"[2]*"FC")) + 
                ylab(expression("-log"[10]*"P.Value")) +
                scale_color_manual(values = colvalors) +
                guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
                geom_vline(xintercept=FCth) + geom_vline(xintercept=-(FCth))+
                #geom_hline(yintercept=y_int)+
                ggtitle(MainTitleVolcano())+ theme(plot.title=element_text( size=8, face='bold'))

             
                if(NgenesV()>0){
                top <- NgenesV()

                top_genes <- bind_rows(
                dataset %>% 
                    filter(Expression == 'Up-regulated') %>% 
                    arrange(P.Value, desc(abs(logFC))) %>% 
                    head(top),
                dataset %>% 
                    filter(Expression == 'Down-regulated') %>% 
                    arrange(P.Value, desc(abs(logFC))) %>% 
                    head(top)
                )

                p3<-p2+ geom_label_repel(data = top_genes,
                                mapping = aes(logFC, -log(P.Value,10), label = SYMBOL),
                                size = 2)
                }

            } # objects not provided by edgeR END 

             
              incProgress(0.89)

           }
          ) #end withProgress
            
            if(NgenesV()>0){    
            output$VolcanoPlot <- renderPlot({p3})

            
            pdf("test2.pdf")
            print(p3)            
            dev.off()
            }else{

            output$VolcanoPlot <- renderPlot({p2})

            pdf("test2.pdf")
            print(p2)            
            dev.off()


            }

      #  output$path<-renderText({wd})

        output$download_Volcano<-downloadHandler(
            filename <- function(){'VolcanoPlot.pdf'},

            content = function(file){
              file.copy("test2.pdf", file)

            }
        )




       })  # collect data Volcano END


       

       session$onSessionEnded(function() {
         system(paste("rm -f ", "test.pdf"))
       })

}

shinyApp(ui = ui, server = server)
