library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)
library(Rtsne)
library(ggplot2)

#library(readr)


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
                title = "t-SNE analysis", status = "danger", solidHeader = TRUE,
                h5("Here you can view a t-SNE analysis of your data. Penultimate column must be the type (i.e: NK cells or CD8). Last column of the matrix must be the condition (i.e: treated or untreated) for each row. "),


                
                fileInput("upload1", NULL, width = 500, buttonLabel = "Upload your merged data...", multiple = FALSE, accept=".txt"),
#                fileInput("upload2", NULL, width = 500, buttonLabel = "Upload your matrix tumor data...", multiple = FALSE, accept=".txt"),

                sliderInput("perplexity", "Perplexity:", 30, 50, 40, width = 120),
                 # sidebarPanel(width=4,
                 #  selectInput("YEARS", "Limit study (in years):",
                 #  choices = c("5","6","7","8","9","10","11","12","13","14","15")),
                 #
                 #  ),
                 

  # sidebarPanel
                actionButton("goButton_Tsne", "Run tSNE!", style = "background-color:#FA1414; color:#FFFFFF"),

                downloadButton("download_tsne", "Download tSNE", class="btn-success"),

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
            h4("here you can view only the main t-SNE plot. However, in the file that you can download you will see 1 plot for each condition and the main plot"),
            plotOutput(outputId = "distPlot",width = "50%"),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
            br(),
      
                    fluidRow(
             column(6,plotOutput(outputId="distPlot2", width="400px",height="400px")),  
             column(6,plotOutput(outputId="distPlot3", width="400px",height="400px"))
                )    
        #    plotOutput(outputId = "distPlot",width = "50%"),   
        #   plotOutput(outputId = "distPlot2",width = "50%"),
        #   plotOutput(outputId = "distPlot3",width = "50%"),
     ) # dashboardbody
)


server = function(input, output, session){
    options(shiny.maxRequestSize=90000000*1024^2)
    session$onSessionEnded(stopApp)


     Matrix_data <- reactive({
         req(input$upload1)
     })   # DATA SELECTED

     #output$files<-renderTable(data())
     perp <- reactive({
         req(input$perplexity)
     })   # GENE SELECTED

     observeEvent(input$goButton_Tsne, {

       validate(
           need(Matrix_data() !="", "Please upload a file!")
       )

       withProgress(
           message = "Running tSNE analyisis...",
           detail = "This step might take a while",
           value = 0,
           {

              incProgress(0.1)

              datamerged<-read.table(Matrix_data()$datapath, header=T, sep="\t") 
               
              #datamerged<-rbind(dataset,dataset2)
              
              columns<-ncol(datamerged)
              tsne <- Rtsne(as.matrix(datamerged[,1:(columns-2)]), check_duplicates = FALSE, pca = FALSE, perplexity=perp(), theta=0.5, dims=2)

             # tsne <- Rtsne(as.matrix(normalize_input(t(dataset[,1:(columns-1)]))), check_duplicates = FALSE, pca = FALSE, perplexity=perp, theta=0.5, dims=2)

             # tsne2 <- Rtsne(as.matrix(dataset[,4:(columns-1)]), check_duplicates = FALSE, pca = FALSE, perplexity=perp, theta=0.5, dims=2)
              # each row is an observation, each column is a variable
              tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = as.factor(datamerged[,columns-1]))
              incProgress(0.89)

           }
          )
            plot_list<-list()

            p<-ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col)) + ggtitle("Global")
            
            output$distPlot <- renderPlot({p})

            k<-length(table(datamerged[,ncol(datamerged)]))
            
            counter<-0
            for (i in 1:k) {
            if(i==1){   
            coord_end<-as.numeric(paste(table(datamerged[,ncol(datamerged)])[i]))
            pl = ggplot(tsne_plot[1:coord_end,]) + geom_point(aes(x=x, y=y, color=col)) + ggtitle(names(table(datamerged[,ncol(datamerged)])[i]))
            counter<-coord_end
            }else{
            dim_end<-as.numeric(paste(table(datamerged[,ncol(datamerged)])[i]))
            coord_end<-counter+dim_end 
            pl = ggplot(tsne_plot[counter+1:coord_end,]) + geom_point(aes(x=x, y=y, color=col)) + ggtitle(names(table(datamerged[,ncol(datamerged)])[i]))
            counter<-counter+dim_end    
            } 
            plot_list[[i]] = pl
            }

            k<-length(plot_list)
            
            plot_list[[k+1]]<-p
              # end of withProgress run analysis
            #aqui
            wd<-getwd()
  #      tmpdir<-tempdir()

            pdf("test.pdf")
            for (i in 1:length(plot_list)) {
            print(plot_list[[i]])
            }            
            dev.off()


      #  output$path<-renderText({wd})

        output$download_tsne<-downloadHandler(
            filename <- function(){'tSNE.pdf'},

            content = function(file){
              file.copy("test.pdf", file)

            }
        )




       })  # collect data

       session$onSessionEnded(function() {
         system(paste("rm -f ", "test.pdf"))
       })

}

shinyApp(ui = ui, server = server)
