
options(shiny.maxRequestSize=90000000000000000000000000*1024^2)

library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)

ui<-dashboardPage(

     skin="yellow",
     dashboardHeader(title = "Programa de Inmunología e Inmunoterapia",
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
            tags$a(href = "/gating/", "Automated gating FC "),
        ),
        tags$li(
            tags$a(href = "/prognosis/", "Prognosis Marker analysis"),
        ),
        tags$li(
            tags$a(href = "/tsne/", "tSNE analysis"),
        ),

        )
      ),  #links navegadores
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

        h5("Here you can carry out the processing of fastq files from an RNA-seq experiment until the generation of the count matrix. This page is optimized for the alignment of PE libraries in the absence of new updates."),
         box(  # FASTQC REPORT FRONTEND
                width = 12,
                title = "Step 1. Quality Control with FastQC tool", status = "danger", solidHeader = TRUE,
                h6("Select your FASTQ files. You must have a forward file '_1' and a reverse file '_2' for each sample. Don´t click in run button until files be fully loaded. Temporary paths to files will appear when they are fully loaded."),
                fluidRow(
                    column(1, offset = 7,
                    actionButton("FQ_help", "FastQC Help",style = "background-color:#FFFFFF; color:#6A93DE; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"), onclick ="window.open('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/', '_blank')")),
                    tags$br()
                ),
                fileInput("upload", NULL, width = 500, buttonLabel = "Select your forward sequences file...", placeholder="Only fastq.gz files", multiple = FALSE, accept=".fastq.gz"),
                fileInput("upload2", NULL, width = 500, buttonLabel = "Select your reverse sequences file...", placeholder="Only fastq.gz files", multiple = FALSE, accept=".fastq.gz"),

                # sidebarPanel(width=5,
                
                #  selectInput("dataset_FQ", "Select the project where you want to store the data:",
                #  choices = c("Intrust", "Linterna", "Independent")),

                 sidebarPanel(width=5,
                 uiOutput('dataset_FQ')), 


              


                fluidRow(
                    column(1, offset = 7,
                    actionButton("goButtonFQ", "Run FastQC!", style = "background-color:#FA1414; color:#FFFFFF", icon= icon("fas fa-running") , lib = "font-awesome"  ),
                    tags$br(),
                    downloadButton("demo", "Download report", class = "btn btn-success"),
                    tags$br()
                    )

                ),
                tableOutput("files"),
                tableOutput("files2"),
		             h4("Cli command example"),
            	  verbatimTextOutput("path2"),
                verbatimTextOutput("datafile")


                #verbatimTextOutput("upload", placeholder = TRUE)
            ),   #### FASTQC REPORT BOX FRONTEND

         box(    #### TRIMMOMATIC BOX FRONTEND
                width = 12,
                title = "Step 2. Trimming low-quality reads and adapters with Trimmomatic ", status = "info", solidHeader = TRUE,
                h6("Select the FASTQ files you want to trim"),
                fluidRow(
                    column(1, offset = 7,
                    actionButton("TR_help", "Trimmomatic Help",style = "background-color:#FFFFFF; color:#6A93DE; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"), onclick ="window.open('http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf', '_blank')")),
                    tags$br()
                ),

                #fileInput("upload3", NULL, width = 500, buttonLabel = "Select your forward sequences file...", placeholder="Only fastq.gz files", multiple = FALSE, accept=".fastq.gz"),
                #fileInput("upload4", NULL, width = 500, buttonLabel = "Select your reverse sequences file...", placeholder="Only fastq.gz files", multiple = FALSE, accept=".fastq.gz"),

                sidebarPanel(width=6,
                 selectInput("upload3", "Select your forward file:", choices= ""),
                ),
                sidebarPanel(width=6,
                 selectInput("upload4", "Select your reverse file:", choices= ""),

                ),

                fluidRow(column(2, tags$br(sliderInput("slider_sliding_min", "Scan the read with a X-base wide sliding window:", 10, 50, 15, width = 120)),
                        #tags$br()
                        ),
                         column(4, tags$br(sliderInput("slider_sliding_max", "Cutting when the average quality per
                          base drops below:", 10, 50, 25, width = 120)),
                         #tags$br()
                        ),
                        column(2, tags$br(sliderInput("slider_MINLEN", "Drop reads which are less than X bases long:", 10, 50, 30, width = 120)),
                        #        tags$br()
                        ),
                ),

                fluidRow(column(4, tags$br(sidebarPanel(width=20,
                        # selectInput("dataset_TRIM", "Select the project where you want to store the data:",
                        # choices = c("Intrust", "Linterna","Independent")))),
                        uiOutput("dataset_TRIM")))
                        ),
                        column(1, offset = 3,

                        actionButton("goButtonTrim", "Run Trimmomatic!", style = "background-color:#FA1414; color:#FFFFFF", icon= icon("fas fa-running") , lib = "font-awesome"   ),
                        tags$br(),
                        #downloadButton("download_TRIM", "Download output files", class = "btn btn-success"),
                    )
                ),

                tableOutput("files3"),
                tableOutput("files4"),

                h4("Command launched"),
                #verbatimTextOutput("filetrim"),
                verbatimTextOutput("pathtrim"),
                #verbatimTextOutput("name_trim_exit"),
                #tableOutput("outVar_P1_TRIM")


                #verbatimTextOutput("upload", placeholder = TRUE)
            ),  ### TRIMMOMATIC BOX FRONTEND

            box(
                width = 12,
                title = "Step 3. Align your data with STAR aligner ", status = "warning", solidHeader = TRUE,
                h6("Select the FASTQ files you want to align. This process may take a long time..."),
                fluidRow(
                    column(1, offset = 7,
                    actionButton("STAR_help", "STAR Help",style = "background-color:#FFFFFF; color:#6A93DE; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"), onclick ="window.open('https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf', '_blank')")),
                    tags$br()
                ),

                #fileInput("upload5", NULL, width = 500, buttonLabel = "Select your forward sequences file...", placeholder="Only fastq.gz files", multiple = FALSE, accept=".fastq.gz"),
                #fileInput("upload6", NULL, width = 500, buttonLabel = "Select your reverse sequences file...", placeholder="Only fastq.gz files", multiple = FALSE, accept=".fastq.gz"),

                sidebarPanel(width=6,
                 selectInput("upload5", "Select your forward file:", choices= ""),
                ),
                sidebarPanel(width=6,
                 selectInput("upload6", "Select your reverse file:", choices= ""),

                ),
                sidebarPanel(width=6,
                 selectInput("trimmed", "Are you trimmed your sample before?", choices=c("YES","NO")),

                ),




                fluidRow(column(2, tags$br(selectInput("dataset_species", "Select the Genome Reference:",
                        choices = c("Human hg38", "Mouse mm39")))),
                        #tags$br()

                        column(4, tags$br(sidebarPanel(width=20,
                        # selectInput("dataset_ALIGN", "Select the project where you want to store the data:",
                        # choices = c("Intrust", "Linterna","None")))),
                            #tags$br()
                            uiOutput("dataset_ALIGN")))
                        )
                ),

                fluidRow(
                    column(1, offset = 7,
                    actionButton("goButtonALIGN", "Run STAR!", style = "background-color:#FA1414; color:#FFFFFF", icon= icon("fas fa-running") , lib = "font-awesome"   ),
                    tags$br(),
                #    downloadButton("download_ALIGN", "Download output alignment files", class = "btn btn-success"),

                )
                ),

                tableOutput("files5"),
                tableOutput("files6"),
                verbatimTextOutput("pathalign"),

                #verbatimTextOutput("upload", placeholder = TRUE)
            ),       #### TRIMMOMATIC BOX FRONTEND
            box(
                width = 12,
                title = "Step 4. Gene expression quantification ", status = "primary", solidHeader = TRUE,
                h6("Select the BAM file created by STAR"),
                fluidRow(
                    column(1, offset = 7,
                    actionButton("FC_help", "featureCounts Help",style = "background-color:#FFFFFF; color:#6A93DE; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"), onclick ="window.open('http://bioinf.wehi.edu.au/featureCounts/', '_blank')")),
                    tags$br()
                ),

                #fileInput("upload7", NULL, width = 500, buttonLabel = "Select your BAM file...", placeholder="Only .BAM files", multiple = FALSE, accept=".BAM"),
                sidebarPanel(width=6,
                 selectInput("upload7", "Select your aligned BAM file:", choices= ""),

                ),



                fluidRow(column(2, tags$br(selectInput("dataset_species_FC", "Select the Genome Reference:",
                        choices = c("Human hg38", "Mouse mm39")))),
                        #tags$br()

                        column(4, tags$br(sidebarPanel(width=20,
                        # selectInput("dataset_FC", "Select the project where you want to store the data:",
                        # choices = c("Intrust", "Linterna","Independent")))),
                        uiOutput("dataset_FC")))
                            #tags$br()
                        ),
                        # column(6, tags$br(selectInput("dataset_species_FC", "Select the Genome Reference:",
                        #         choices = c("Human hg38", "Mouse mm39")))),

                ),

                fluidRow(
                    column(1, offset = 7,
                    actionButton("goButtonFC", "Run featureCounts!", style = "background-color:#FA1414; color:#FFFFFF", icon= icon("fas fa-running") , lib = "font-awesome"  ),
                    tags$br(),
                    downloadButton("download_Features", "Download fC output files", class = "btn btn-success"),

                )
                ),

                tableOutput("files7"),

                verbatimTextOutput("featurecounts"),
                verbatimTextOutput("NAMEFILE_FC"),
            )       #### featureCounts


     ) # dashboardbody

)


server = function(input, output, session){
      options(shiny.maxRequestSize=9000000000000000000000000000*1024^2)
      session$onSessionEnded(stopApp)


        data <- reactive({
        req(input$upload)
    })   # Fichero PE forward fastQC

        datab <- reactive({
        req(input$upload2)
    })  #fichero PE reverse FastQC

     myprojects = list.dirs("/media/inmuno/data_projects", recursive=FALSE)
     myprojects<-basename(myprojects)

     output$dataset_FQ <- renderUI({
       selectInput(inputId = 'dataset_FQ', label = 'Choose project to store the data',
                   choices = myprojects)
     })

     output$dataset_ALIGN <- renderUI({
       selectInput(inputId = 'dataset_ALIGN', label = 'Choose project to store the data',
                   choices = myprojects)
     })

     output$dataset_TRIM <- renderUI({
       selectInput(inputId = 'dataset_TRIM', label = 'Choose project to store the data',
                   choices = myprojects)
     })

     output$dataset_FC <- renderUI({
       selectInput(inputId = 'dataset_FC', label = 'Choose project to store the data',
                   choices = myprojects)
     })


        outVar_P1_TRIM = reactive({
          mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_TRIM, "/RNA_seq_data/FASTQ_Files/"), pattern=".gz")
          mydata<-as.list(mydata)
          #return(mydata)
          #names(mydata)
    })  # fichero PE trimmomatic

        observe({
            updateSelectInput(session, "upload3",
            choices = outVar_P1_TRIM()
    )})

        outVar_P2_TRIM = reactive({
          mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_TRIM, "/RNA_seq_data/FASTQ_Files/"), pattern=".gz")
          mydata<-as.list(mydata)
          #return(mydata)
          #names(mydata)
    })  # fichero PE trimmomatic

        observe({
            updateSelectInput(session, "upload4",
            choices = outVar_P2_TRIM()
    )})


        datac <- reactive({
        req(input$upload3)
        datac<-paste0("/media/inmuno/data_projects/", input$dataset_TRIM, "/RNA_seq_data/FASTQ_Files/", input$upload3)
    }) #Fichero PE forward Trimmomatic

        datad <- reactive({
        req(input$upload4)
        datad<-paste0("/media/inmuno/data_projects/", input$dataset_TRIM, "/RNA_seq_data/FASTQ_Files/", input$upload4)
    }) #Fichero PE reverse Trimmomatic


    outVar_P1_ALIGN = reactive({
      if(input$trimmed=="YES"){
      mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_ALIGN, "/RNA_seq_data/FASTQC_Trimmed/"), pattern="T.fastq.gz")
      }else{
        mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_ALIGN, "/RNA_seq_data/FASTQ_Files/"), pattern=".gz")
      }
      })  # fichero PE1  align
    observe({
        updateSelectInput(session, "upload5",
        choices = outVar_P1_ALIGN()
        )})

        datae <- reactive({
        req(input$upload5)
        #tmp<-paste(lapply(strsplit(paste(input$upload5),"_"),"[",1))
        datae<-paste0("/media/inmuno/data_projects/", input$dataset_ALIGN, "/RNA_seq_data/FASTQC_Trimmed/", input$upload5)

    }) #Fichero PE forward STAR

    outVar_P2_ALIGN = reactive({
      if(input$trimmed=="YES"){
      mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_ALIGN, "/RNA_seq_data/FASTQC_Trimmed/"), pattern="T.fastq.gz")
      }else{
        mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_ALIGN, "/RNA_seq_data/FASTQ_Files/"), pattern=".gz")
      }
      })  # fichero PE1  align
    observe({
        updateSelectInput(session, "upload6",
        choices = outVar_P2_ALIGN()
        )})

        dataf <- reactive({
        req(input$upload6)
        #tmp<-paste(lapply(strsplit(paste(input$upload5),"_"),"[",1))
        dataf<-paste0("/media/inmuno/data_projects/", input$dataset_ALIGN, "/RNA_seq_data/FASTQC_Trimmed/", input$upload6)

    }) #Fichero PE reverse STAR


    outVar_FC = reactive({

      mydata = list.files(paste0("/media/inmuno/data_projects/",input$dataset_FC, "/RNA_seq_data/STAR_results/"), pattern=".bam")
      mydata<-as.list(mydata)
      #return(mydata)
      #names(mydata)
      })  # fichero PE1  align
    observe({
        updateSelectInput(session, "upload7",
        choices = outVar_FC()
        )})

        datag <- reactive({
        req(input$upload7)
        #tmp<-paste(lapply(strsplit(paste(input$upload5),"_"),"[",1))
        datag<-paste0("/media/inmuno/data_projects/", input$dataset_FC, "/RNA_seq_data/STAR_results/", input$upload7)

    }) #Fichero BAM FC


      output$files<-renderTable(data())  # show routes to temporary folders
      output$files2<-renderTable(datab())
      output$files3<-renderTable(datac())
      output$files4<-renderTable(datad())
      output$files5<-renderTable(datae())
      output$files6<-renderTable(dataf())
      output$files7<-renderTable(datag())
      #output$outVar_P1_TRIM<-renderTable(outVar_P1_TRIM())
#     browser()



        #FASTQC BOX
      observeEvent(input$goButtonFQ, {
          validate(
              need(data()$datapath !="", "Please upload your files before!"),
              need(datab()$datapath !="", "Please upload your files before!")
          )
        withProgress(
            message = "Running FastQC on your data...",
            detail = "This step might take a while",
            value = 0,
            {
            a<-reactive({input$dataset_FQ})
            b<-reactive({data()$name})
      	    c<-reactive({datab()$name})

            ProjectOutput_FQ<-paste0("/media/inmuno/data_projects/",a(),"/RNA_seq_data/FASTQC_Report/",b(),"/")
            ProjectOutput_FQ2<-paste0("/media/inmuno/data_projects/",a(),"/RNA_seq_data/FASTQC_Report/",c(),"/")

            dir.create(ProjectOutput_FQ)
	          dir.create(ProjectOutput_FQ2)

            #cat("FASTQC ", data()$datapath, "\n")
            a<-reactive({input$dataset_FQ})
            b<-reactive({data()$name})
      	    c<-reactive({datab()$name})

            saveFile_PATH_FQ<-paste0("/media/inmuno/data_projects/",a(),"/RNA_seq_data/FASTQ_Files/")
            dir.create(saveFile_PATH_FQ)
            saveFile_Name<-paste0(saveFile_PATH_FQ, b())
            saveFile_Name2<-paste0(saveFile_PATH_FQ, c())

            incProgress(0.1)
            y<-paste0("fastqc -o ", ProjectOutput_FQ ," ", data()$datapath)
            output$path2<-renderText({y})
            system(y)

            if(a()!="None"){
              yy<-paste0("mv ", data()$datapath," ", saveFile_Name)
              system(yy)
            }

            incProgress(0.5)
            #cat("FASTQC ", datab()$datapath, "\n")
            y2<-paste0("fastqc -o ", ProjectOutput_FQ2 ," ", datab()$datapath)
            output$path3<-renderText({y2})
            system(y2)

            if(a()!="None"){
              yy2<-paste0("mv ", datab()$datapath," ", saveFile_Name2)
              system(yy2)
            }


            incProgress(0.89)
            }


        )
        # session$onSessionEnded(function() {
        #   if(a()=="None"){
        #   system(paste("rm -rf ", ProjectOutput_FQ))
        #   system(paste("rm -rf ", ProjectOutput_FQ2))
        #
        #   }
        # })

        wd<-getwd()

        output$demo<-downloadHandler(

                filename <- function(){'data.zip'},
                content = function(file){


                  tmpdir<-tempdir()
                  # rename file 1
                  filesN<-list.files(ProjectOutput_FQ, pattern="html")
                  filesN<-paste0(ProjectOutput_FQ, filesN)
                  z<-paste0("mv ",filesN," ",paste0(ProjectOutput_FQ, data()$name,".html"))
                  system(z)
                  filesZ<-list.files(ProjectOutput_FQ, pattern=".html")
                  filesZ<-paste0(ProjectOutput_FQ, filesZ)
                  file.copy(filesZ, tmpdir)
                  #rename file 2

                  filesQ<-list.files(ProjectOutput_FQ2, pattern="html")
                  filesQ<-paste0(ProjectOutput_FQ2, filesQ)
                  z<-paste0("mv ",filesQ," ", paste0(ProjectOutput_FQ2, datab()$name,".html"))
                  system(z)
                  filesM<-list.files(ProjectOutput_FQ2, pattern=".html")
                  filesM<-paste0(ProjectOutput_FQ2, filesM)
                  file.copy(filesM, tmpdir)
                  # temporary path
                  setwd(tmpdir)
                  testdata<-as.character(list.files(getwd(), pattern=".html"))
                  print(testdata)
                  zip(file, files=testdata)
                  setwd(wd)
                  unlink(tmpdir)
                  },
                # content <- function(file) {
                #   file.copy(page, file)
                #   },
              contentType="application/zip"
            )

      }) #LAUNCH FASTQC




      ### TRIMMOMATIC BOX
      observeEvent(input$goButtonTrim, {


        withProgress(
            message = "Running Trimmomatic on your data...",
            detail = "This step might take a while",
            value = 0,
            {


            d<-reactive({input$dataset_TRIM})

            g<-reactive({input$slider_sliding_min})
            h<-reactive({input$slider_sliding_max})
            i<-reactive({input$slider_MINLEN})
            #output$filetrim<-renderText({e()})

            NAMEFILE_tmp<-paste(lapply(strsplit(paste(datac()),"\\."),"[",1))
            NAMEFILE_tmp2<-basename(NAMEFILE_tmp)
            NAMEFILE<-paste(lapply(strsplit(paste(NAMEFILE_tmp2),"_"),"[",1))

            ProjectOutput_TRIM<-paste0("/media/inmuno/data_projects/",d(),"/RNA_seq_data/FASTQC_Trimmed/")
            dir.create(ProjectOutput_TRIM)
            #output$filetrim<-renderText({datac()})



            #output$namefile_trim<-renderText({NAMEFILE})
            incProgress(0.1)
            y<-paste0("java -jar /usr/share/java/trimmomatic-0.39.jar PE -phred33 ", datac()," ", datad()," ", paste0(ProjectOutput_TRIM,NAMEFILE,"_1_T.fastq.gz "), paste0(ProjectOutput_TRIM,NAMEFILE,"_1_T_Unpaired.fastq.gz "), paste0(ProjectOutput_TRIM,NAMEFILE,"_2_T.fastq.gz "), paste0(ProjectOutput_TRIM,NAMEFILE,"_2_T_Unpaired.fastq.gz "),"ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE-2.fa/:2:30:10 " ,"SLIDINGWINDOW:",g(),":",h()," MINLEN:",i())

            output$pathtrim<-renderText({y})
            system(y)
            incProgress(0.89)

            }) #LAUNCH TRIMMOMATIC


        })   # TRIMMOMATIC BOX




    #
    #
    # #   ### ALIGNEMENT BOX
    #
      observeEvent(input$goButtonALIGN, {

       withProgress(
            message = "Aligning with STAR your data...",
            detail = "This step might take a long time",
            value = 0,
        {
          d<-reactive({input$dataset_ALIGN})

          NAMEFILE_tmp3<-paste(lapply(strsplit(paste(datae()),"\\."),"[",1))
          NAMEFILE_tmp4<-basename(NAMEFILE_tmp3)
          NAMEFILE_ALIGN<-paste(lapply(strsplit(paste(NAMEFILE_tmp4),"_"),"[",1))


         ProjectOutput_ALIGN<-paste0("/media/inmuno/data_projects/",d(),"/RNA_seq_data/STAR_results/")
         dir.create(ProjectOutput_ALIGN)


         incProgress(0.1)
         if(input$dataset_species =="Human hg38"){
         j<-paste0("STAR --genomeDir /media/inmuno/references/STAR_human_hg38/ --readFilesCommand zcat --runThreadN 6 --readFilesIn ", datae() ," ", dataf()," ","--outFileNamePrefix ", paste0(ProjectOutput_ALIGN,NAMEFILE_ALIGN)," ","--outSAMtype BAM ", "SortedByCoordinate ", "--outSAMunmapped Within --outSAMattributes Standard")
       }else{
         j<-paste0("STAR --genomeDir /media/inmuno/references/STAR_mouse_mm39/ --readFilesCommand zcat --runThreadN 6 --readFilesIn ", datae() ," ", dataf()," ","--outFileNamePrefix ", paste0(ProjectOutput_ALIGN,NAMEFILE_ALIGN)," ","--outSAMtype BAM ", "SortedByCoordinate ", "--outSAMunmapped Within --outSAMattributes Standard")
       }
         #system(j)
         output$pathalign<-renderText({j})
         incProgress(0.89)
        }
        )
      })   #LAUNCH STAR 
    #

    #
    #   ### featureCounts BOX
            # Launch featureCounts
      observeEvent(input$goButtonFC, {

       withProgress(
            message = "Extracting counts from your data...",
            detail = "This step might take a while",
            value = 0,
        {
         d<-reactive({input$dataset_FC})

         NAMEFILE_FC_tmp<-paste(lapply(strsplit(paste(datag()),"\\."),"[", 1))
         NAMEFILE_FC<-basename(NAMEFILE_FC_tmp)
         output$NAMEFILE_FC<-renderText({NAMEFILE_FC})


         ProjectOutput_FC<-paste0("/media/inmuno/data_projects/",d(),"/RNA_seq_data/featureCounts_results/")
         dir.create(ProjectOutput_FC)

         incProgress(0.1)
         if(input$dataset_species_FC=="Human hg38"){
         y<-paste0("featureCounts -p -t exon -g gene_id -a $PATH_TO_GTF_HUMAN -o ", paste0(ProjectOutput_FC, NAMEFILE_FC,"_counts.txt " ), datag())
       }else{
         y<-paste0("featureCounts -p -t exon -g gene_id -a $PATH_TO_GTF_MOUSE -o ", paste0(ProjectOutput_FC, NAMEFILE_FC,"_counts.txt " ), datag())
       }
         #system(y)
         output$featurecounts<-renderText({y})
         incProgress(0.89)
        }
        )

        wd<-getwd()

        output$download_Features<-downloadHandler(
            filename <- function(){'Counts.zip'},

            content = function(file){

                tmpdir<-tempdir()

                all_files<-list.files(ProjectOutput_FC, pattern=NAMEFILE_FC)
                all_files<-paste0(ProjectOutput_FC, all_files)

                file.copy(all_files, tmpdir)
                # temporary path
                setwd(tmpdir)
                testdata<-as.character(list.files(getwd(), pattern=".txt"))
                print(testdata)
                zip(file, files=testdata)
                setwd(wd)
                unlink(tmpdir)
            },
        )
        session$onSessionEnded(function() {
          if(d()=="Independent"){
          system("find /media/inmuno/data_projects/Independent/ ! -type d -delete")

          }
        })   # delete files if query is of an independent project once files are downloaded

      })   #LAUNCH featureCounts


}


shinyApp(ui = ui, server = server)
