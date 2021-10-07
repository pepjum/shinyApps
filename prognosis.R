library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)
library(TCGAbiolinks)
library("SummarizedExperiment")
library(dplyr)
library(survival)
library(ggplot2)
library(survminer)
library(edgeR)
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
            tags$a(href = "/ideal/", "Perform RNA-Seq DE analysis"),
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
                title = "Prognosis marker analysis", status = "danger", solidHeader = TRUE,
                h5("Check the TCGA database to see if the expression of your gene has a prognosis correlate. Press the GDC button if you need some explanation about TCGA Study Abbreviations."),
                fluidRow(
                    column(1, offset = 7,
                    actionButton("GDC_help", "GDC abbreviations",style = "background-color:#FFFFFF; color:#6A93DE; border-color:#BEBEBE; border-style:none; border-width:1px; border-radius:81%; font-size:18px;", icon = icon("th"), onclick ="window.open('https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations/', '_blank')")),
                    tags$br()
                ),


                sidebarPanel(width=3,
                 selectInput("TISSUE", "Choose a Study:",
                 choices = c("LAML", "ACC","BLCA","LGG","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THYM","THCA","UCS","UCEC","UVM")),

                 ),
                 sidebarPanel(width=4,
                  textInput("GENE", "Enter ENSG code of your Gene:", ""),
                 ),
                 # sidebarPanel(width=4,
                 #  selectInput("YEARS", "Limit study (in years):",
                 #  choices = c("5","6","7","8","9","10","11","12","13","14","15")),
                 #
                 #  ),
                 sidebarPanel(width=4,
                  selectInput("metrics", "Expression cut-off for group delimitation (LOW-HIGH):",
                  choices=c("optimal","1st quartile and 3rd quartile", "1st quartile and mean", "1st quartile and median", "mean", "median")),
                 ),

                 sliderInput("YEARS", "Limit study to (in years):", 5, 50, 25, width = 250),

  # sidebarPanel
                actionButton("goButton_flow", "Run analysis!", style = "background-color:#FA1414; color:#FFFFFF"),

                downloadButton("download_kaplan", "Download Kaplan Meier", class="btn-success"),

            #    actionButton("plot_surv", "Run Analysis!", style = "background-color:#FA1414; color:#FFFFFF"),

            ),
            verbatimTextOutput("path"),

            verbatimTextOutput("geneExpression"),

            plotOutput(outputId = "distPlot",width = "90%"),

     ) # dashboardbody
)


server = function(input, output, session){
    options(shiny.maxRequestSize=90000000*1024^2)
    session$onSessionEnded(stopApp)


    TCGA_data <- reactive({
        req(input$TISSUE)
    })   # DATA SELECTED

     #output$files<-renderTable(data())
     ENSG <- reactive({
         req(input$GENE)
     })   # GENE SELECTED

     metrics <- reactive({
         req(input$metrics)
     })   # metrics SELECTED

     years <- reactive({
         req(input$YEARS)
     })   # metrics SELECTED


     cancer_types<-reactive({paste0("TCGA-",TCGA_data())})

     output$query<-renderText({cancer_types()})
     output$metrics<-renderText({metrics()})
     output$months<-renderText({years()*12})
     observeEvent(input$goButton_flow, {

       validate(
           need(ENSG() !="", "Please INTRODUCE A VALID ENSG!")
       )

       withProgress(
           message = "Fetching expression matrix...",
           detail = "This step might take a while",
           value = 0,
           {
            cancer_types<-reactive({paste0("TCGA-",TCGA_data())})

            kk<-cancer_types()

            if(file.exists(paste0("/media/inmuno/TCGA_database/",kk,".Rdata"))){
              #setwd("/media/inmuno/data/TCGA_database/")

              incProgress(0.1)

              cancer_df<-get(load(paste0("/media/inmuno/TCGA_database/",cancer_types(),".Rdata")))

              incProgress(0.89)

            }else{

           incProgress(0.1)

           #setwd("/media/inmuno/data/TCGA_database/")

           #output$exists<-renderText({"No entra"})

           query_target <- GDCquery(project = cancer_types() ,data.category="Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",workflow.type = "HTSeq - Counts",
                          experimental.strategy = "RNA-Seq")
           #incProgress(0.3)
            GDCdownload(query_target)
           #incProgress(0.5)
            datas<-GDCprepare(query_target)
           #incProgress(0.7)
            cancer_df<-assay(datas)
            save(cancer_df, file=paste0("/media/inmuno/TCGA_database/",cancer_types(),".Rdata"))
           incProgress(0.89)
          }
          # system("rm -r /media/inmuno/data/TCGA_database/GDCData/*")
           })


       withProgress(
           message = "Fetching clinical data...",
           detail = "This step might take a while",
           value = 0,
           {

            incProgress(0.1)
           #
            clinical <- GDCquery_clinic(project = cancer_types(), type = "clinical")
            #save(clinical, file="/media/inmuno/data/test_clinical_ago21.R")

           incProgress(0.89)

           })


           withProgress(

               message = "Normalizing matrix...",
               detail = "This step might take a while",
               value = 0,
               {

               incProgress(0.1)

               d0<-DGEList(cancer_df)

               d0<-calcNormFactors(d0)

               cutoff<-5
               drop<-which(apply(cpm(d0),1,max) < cutoff)
               d<-d0[-drop,]
               if(ENSG() %in% names(drop)){

                 output$geneExpression<-renderText({"This gene has no expression here"})
                 break
               }

               incProgress(0.5)
               y<-voom(d, design=NULL, plot=F)

               cancer_df<-y$E
               incProgress(0.89)
               })



           withProgress(
               message = "Running analysis...",
               detail = "This step might take a while",
               value = 0,
               {

               incProgress(0.1)
               headers<-colnames(cancer_df)
               headers_samples<-paste(lapply(strsplit(paste(headers),"-"),"[", 4))
               tumor_samples<-which(as.numeric(substr(headers_samples,1,2))<10)
               Tissue_tumor<-cancer_df[,tumor_samples]

               Gene_data<-Tissue_tumor[c(ENSG()),]

               #save(Gene_data, file="/media/inmuno/data/test_ago17.R")

               samples<-names(Gene_data)
               samples_cleaned<-paste(lapply(strsplit(paste(samples),"-0"),"[",1))

               names(Gene_data)<-samples_cleaned

               #median(Gene_data)
               #[1] 56.79155


               if(metrics()=="optimal"){

                 Gene_data_df<-as.data.frame(Gene_data)
                 Gene_data_df$sample<-names(Gene_data)
                 names(Gene_data_df)[1]<-"Expression"
                 Gene_data_df$sample<-paste(lapply(strsplit(paste(Gene_data_df$sample),"-0"),"[",1))

                 clinical<-as.data.frame(clinical)
                 clinical<-clinical[,c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up")]
                 clinical$time<-ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
                 clinical<-clinical[,c("submitter_id","vital_status","time")]
                 clinical<-unique(clinical)

                 Gene_data_annot<-merge(Gene_data_df, clinical, by.x="sample", by.y="submitter_id")
                 Gene_data_annot$Status<-ifelse(Gene_data_annot$vital_status=="Dead", 1, 0)
                 Gene_data_annot$time<-as.numeric(Gene_data_annot$time)
                 Gene_data_annot<-na.omit(Gene_data_annot)
                 Gene_data_annot$Status<-as.numeric(Gene_data_annot$Status)
                 Gene_data_annot$time<-(Gene_data_annot$time)/30
                 months<-as.numeric(years())*12
                 Gene_data_annot<-Gene_data_annot[which(Gene_data_annot$time <= months),]

                 res.cut<-surv_cutpoint(Gene_data_annot, time="time", event="Status", variables=c("Expression"))
                 res.cat<-surv_categorize(res.cut)
                 data_prepared<-res.cat
                 model_fit<-survfit(Surv(time, Status)~Expression, data=res.cat)

                 GENE<-ENSG()
                 CANCER<-TCGA_data()
                 p<-ggsurvplot(model_fit, res.cat,
                 pval = TRUE, conf.int = TRUE,
                 title=paste0(GENE," expression in ",CANCER),
                 legend.title="Expression",
                 xlab="Months",
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                 legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"))

                 output$distPlot <- renderPlot({p})


                 incProgress(0.89)


               }else if(metrics()=="1st quartile and 3rd quartile"){

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[2])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[5])

                 GROUP_LOW<-Gene_data[GROUP_LOW_index ]
                 GROUP_LOW_df<-as.data.frame(GROUP_LOW)
                 GROUP_LOW_df$sample<-names(GROUP_LOW)
                 names(GROUP_LOW_df)[1]<-"Expression"

                 output_low_df<-data.frame()
                 for(sample in unique(GROUP_LOW_df$sample)){
                   subset_df<-GROUP_LOW_df[which(GROUP_LOW_df$sample ==sample),]
                   kk<-subset_df[1,]
                   output_low_df<-rbind(output_low_df,kk)
                 }
                 rownames(output_low_df)<-paste(output_low_df$sample)

                 GROUP_HIGH<-Gene_data[GROUP_HIGH_index ]
                 GROUP_HIGH_df<-as.data.frame(GROUP_HIGH)
                 GROUP_HIGH_df$sample<-names(GROUP_HIGH)
                 names(GROUP_HIGH_df)[1]<-"Expression"

                 output_high_df<-data.frame()
                 for(sample in unique(GROUP_HIGH_df$sample)){
                   subset_df<-GROUP_HIGH_df[which(GROUP_HIGH_df$sample ==sample),]
                   subset_df<-subset_df[1,]
                   output_high_df<-rbind(output_high_df,subset_df)
                 }
                 rownames(output_high_df)<-paste(output_high_df$sample)

                 clinical<-as.data.frame(clinical)
                 clinical<-clinical[,c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up")]
                 clinical$time<-ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
                 clinical<-clinical[,c("submitter_id","vital_status","time")]
                 clinical<-unique(clinical)

                 output_low_df$sample<-paste(lapply(strsplit(paste(output_low_df$sample),"-01"),"[",1))
                 output_high_df$sample<-paste(lapply(strsplit(paste(output_high_df$sample),"-01"),"[",1))

                 GL_annot<-merge(output_low_df, clinical, by.x="sample", by.y="submitter_id")
                 GL_annot$group<-"Low"
                 GL_annot<-unique(GL_annot)

                 GH_annot<-merge(output_high_df, clinical, by.x="sample", by.y="submitter_id")
                 GH_annot$group<-"High"
                 GH_annot<-unique(GH_annot)

                 GH_annot$Status<-ifelse(GH_annot$vital_status=="Dead", 1, 0)
                 GH_annot$group<-as.factor(GH_annot$group)
                 GH_annot$time<-as.numeric(GH_annot$time)
                 GH_annot<-na.omit(GH_annot)
                 GH_annot$Status<-as.numeric(GH_annot$Status)

                 GL_annot$Status<-ifelse(GL_annot$vital_status=="Dead", 1, 0)
                 GL_annot$group<-as.factor(GL_annot$group)
                 GL_annot$time<-as.numeric(GL_annot$time)
                 GL_annot<-na.omit(GL_annot)
                 GL_annot$Status<-as.numeric(GL_annot$Status)

                 data_prepared<-rbind(GH_annot,GL_annot)
                 data_prepared$time<-(data_prepared$time)/30
                 months<-as.numeric(years())*12

                 data_prepared<-data_prepared[which(data_prepared$time <= months),]
                 model_fit<-survfit(Surv(time, Status) ~ group, data=data_prepared)

                 GENE<-ENSG()
                 CANCER<-TCGA_data()
                 p<-ggsurvplot(model_fit, data_prepared,
                 pval = TRUE, conf.int = TRUE,
                 title=paste0(GENE," expression in ",CANCER),
                 legend.title="Expression",
                 xlab="Months",
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                 legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"))

                 output$distPlot <- renderPlot({p})

                 incProgress(0.89)


               }else if(metrics()=="median"){

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[3])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[3])

                 GROUP_LOW<-Gene_data[GROUP_LOW_index ]
                 GROUP_LOW_df<-as.data.frame(GROUP_LOW)
                 GROUP_LOW_df$sample<-names(GROUP_LOW)
                 names(GROUP_LOW_df)[1]<-"Expression"

                 output_low_df<-data.frame()
                 for(sample in unique(GROUP_LOW_df$sample)){
                   subset_df<-GROUP_LOW_df[which(GROUP_LOW_df$sample ==sample),]
                   kk<-subset_df[1,]
                   output_low_df<-rbind(output_low_df,kk)
                 }
                 rownames(output_low_df)<-paste(output_low_df$sample)

                 GROUP_HIGH<-Gene_data[GROUP_HIGH_index ]
                 GROUP_HIGH_df<-as.data.frame(GROUP_HIGH)
                 GROUP_HIGH_df$sample<-names(GROUP_HIGH)
                 names(GROUP_HIGH_df)[1]<-"Expression"

                 output_high_df<-data.frame()
                 for(sample in unique(GROUP_HIGH_df$sample)){
                   subset_df<-GROUP_HIGH_df[which(GROUP_HIGH_df$sample ==sample),]
                   subset_df<-subset_df[1,]
                   output_high_df<-rbind(output_high_df,subset_df)
                 }
                 rownames(output_high_df)<-paste(output_high_df$sample)

                 clinical<-as.data.frame(clinical)
                 clinical<-clinical[,c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up")]
                 clinical$time<-ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
                 clinical<-clinical[,c("submitter_id","vital_status","time")]
                 clinical<-unique(clinical)

                 output_low_df$sample<-paste(lapply(strsplit(paste(output_low_df$sample),"-01"),"[",1))
                 output_high_df$sample<-paste(lapply(strsplit(paste(output_high_df$sample),"-01"),"[",1))

                 GL_annot<-merge(output_low_df, clinical, by.x="sample", by.y="submitter_id")
                 GL_annot$group<-"Low"
                 GL_annot<-unique(GL_annot)
                 #save(GL_annot, file="/media/inmuno/data/GL_annot.R")

                 GH_annot<-merge(output_high_df, clinical, by.x="sample", by.y="submitter_id")
                 GH_annot$group<-"High"
                 GH_annot<-unique(GH_annot)

                 GH_annot$Status<-ifelse(GH_annot$vital_status=="Dead", 1, 0)
                 GH_annot$group<-as.factor(GH_annot$group)
                 GH_annot$time<-as.numeric(GH_annot$time)
                 GH_annot<-na.omit(GH_annot)
                 GH_annot$Status<-as.numeric(GH_annot$Status)

                 GL_annot$Status<-ifelse(GL_annot$vital_status=="Dead", 1, 0)
                 GL_annot$group<-as.factor(GL_annot$group)
                 GL_annot$time<-as.numeric(GL_annot$time)
                 GL_annot<-na.omit(GL_annot)
                 GL_annot$Status<-as.numeric(GL_annot$Status)

                 data_prepared<-rbind(GH_annot,GL_annot)
                 data_prepared$time<-(data_prepared$time)/30
                 months<-as.numeric(years())*12

                 data_prepared<-data_prepared[which(data_prepared$time <= months),]
                 model_fit<-survfit(Surv(time, Status) ~ group, data=data_prepared)

                 GENE<-ENSG()
                 CANCER<-TCGA_data()
                 p<-ggsurvplot(model_fit, data_prepared,
                 pval = TRUE, conf.int = TRUE,
                 title=paste0(GENE," expression in ",CANCER),
                 legend.title="Expression",
                 xlab="Months",
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                 legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"))

                 output$distPlot <- renderPlot({p})

                 incProgress(0.89)


               }else if(metrics()=="mean"){

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[4])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[4])

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[3])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[3])

                 GROUP_LOW<-Gene_data[GROUP_LOW_index ]
                 GROUP_LOW_df<-as.data.frame(GROUP_LOW)
                 GROUP_LOW_df$sample<-names(GROUP_LOW)
                 names(GROUP_LOW_df)[1]<-"Expression"

                 output_low_df<-data.frame()
                 for(sample in unique(GROUP_LOW_df$sample)){
                   subset_df<-GROUP_LOW_df[which(GROUP_LOW_df$sample ==sample),]
                   kk<-subset_df[1,]
                   output_low_df<-rbind(output_low_df,kk)
                 }
                 rownames(output_low_df)<-paste(output_low_df$sample)

                 GROUP_HIGH<-Gene_data[GROUP_HIGH_index ]
                 GROUP_HIGH_df<-as.data.frame(GROUP_HIGH)
                 GROUP_HIGH_df$sample<-names(GROUP_HIGH)
                 names(GROUP_HIGH_df)[1]<-"Expression"

                 output_high_df<-data.frame()
                 for(sample in unique(GROUP_HIGH_df$sample)){
                   subset_df<-GROUP_HIGH_df[which(GROUP_HIGH_df$sample ==sample),]
                   subset_df<-subset_df[1,]
                   output_high_df<-rbind(output_high_df,subset_df)
                 }
                 rownames(output_high_df)<-paste(output_high_df$sample)

                 clinical<-as.data.frame(clinical)
                 clinical<-clinical[,c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up")]
                 clinical$time<-ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
                 clinical<-clinical[,c("submitter_id","vital_status","time")]
                 clinical<-unique(clinical)

                 output_low_df$sample<-paste(lapply(strsplit(paste(output_low_df$sample),"-01"),"[",1))
                 output_high_df$sample<-paste(lapply(strsplit(paste(output_high_df$sample),"-01"),"[",1))

                 GL_annot<-merge(output_low_df, clinical, by.x="sample", by.y="submitter_id")
                 GL_annot$group<-"Low"
                 GL_annot<-unique(GL_annot)
                 #save(GL_annot, file="/media/inmuno/data/GL_annot.R")

                 GH_annot<-merge(output_high_df, clinical, by.x="sample", by.y="submitter_id")
                 GH_annot$group<-"High"
                 GH_annot<-unique(GH_annot)

                 GH_annot$Status<-ifelse(GH_annot$vital_status=="Dead", 1, 0)
                 GH_annot$group<-as.factor(GH_annot$group)
                 GH_annot$time<-as.numeric(GH_annot$time)
                 GH_annot<-na.omit(GH_annot)
                 GH_annot$Status<-as.numeric(GH_annot$Status)

                 GL_annot$Status<-ifelse(GL_annot$vital_status=="Dead", 1, 0)
                 GL_annot$group<-as.factor(GL_annot$group)
                 GL_annot$time<-as.numeric(GL_annot$time)
                 GL_annot<-na.omit(GL_annot)
                 GL_annot$Status<-as.numeric(GL_annot$Status)

                 data_prepared<-rbind(GH_annot,GL_annot)
                 data_prepared$time<-(data_prepared$time)/30
                 months<-as.numeric(years())*12

                 data_prepared<-data_prepared[which(data_prepared$time <= months),]
                 model_fit<-survfit(Surv(time, Status) ~ group, data=data_prepared)

                 GENE<-ENSG()
                 CANCER<-TCGA_data()
                 p<-ggsurvplot(model_fit, data_prepared,
                 pval = TRUE, conf.int = TRUE,
                 title=paste0(GENE," expression in ",CANCER),
                 legend.title="Expression",
                 xlab="Months",
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                 legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"))

                 output$distPlot <- renderPlot({p})

                 incProgress(0.89)


               }else if(metrics()=="1st quartile and median"){

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[2])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[3])

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[3])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[3])

                 GROUP_LOW<-Gene_data[GROUP_LOW_index ]
                 GROUP_LOW_df<-as.data.frame(GROUP_LOW)
                 GROUP_LOW_df$sample<-names(GROUP_LOW)
                 names(GROUP_LOW_df)[1]<-"Expression"

                 output_low_df<-data.frame()
                 for(sample in unique(GROUP_LOW_df$sample)){
                   subset_df<-GROUP_LOW_df[which(GROUP_LOW_df$sample ==sample),]
                   kk<-subset_df[1,]
                   output_low_df<-rbind(output_low_df,kk)
                 }
                 rownames(output_low_df)<-paste(output_low_df$sample)

                 GROUP_HIGH<-Gene_data[GROUP_HIGH_index ]
                 GROUP_HIGH_df<-as.data.frame(GROUP_HIGH)
                 GROUP_HIGH_df$sample<-names(GROUP_HIGH)
                 names(GROUP_HIGH_df)[1]<-"Expression"

                 output_high_df<-data.frame()
                 for(sample in unique(GROUP_HIGH_df$sample)){
                   subset_df<-GROUP_HIGH_df[which(GROUP_HIGH_df$sample ==sample),]
                   subset_df<-subset_df[1,]
                   output_high_df<-rbind(output_high_df,subset_df)
                 }
                 rownames(output_high_df)<-paste(output_high_df$sample)

                 clinical<-as.data.frame(clinical)
                 clinical<-clinical[,c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up")]
                 clinical$time<-ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
                 clinical<-clinical[,c("submitter_id","vital_status","time")]
                 clinical<-unique(clinical)

                 output_low_df$sample<-paste(lapply(strsplit(paste(output_low_df$sample),"-01"),"[",1))
                 output_high_df$sample<-paste(lapply(strsplit(paste(output_high_df$sample),"-01"),"[",1))

                 GL_annot<-merge(output_low_df, clinical, by.x="sample", by.y="submitter_id")
                 GL_annot$group<-"Low"
                 GL_annot<-unique(GL_annot)
                 #save(GL_annot, file="/media/inmuno/data/GL_annot.R")

                 GH_annot<-merge(output_high_df, clinical, by.x="sample", by.y="submitter_id")
                 GH_annot$group<-"High"
                 GH_annot<-unique(GH_annot)

                 GH_annot$Status<-ifelse(GH_annot$vital_status=="Dead", 1, 0)
                 GH_annot$group<-as.factor(GH_annot$group)
                 GH_annot$time<-as.numeric(GH_annot$time)
                 GH_annot<-na.omit(GH_annot)
                 GH_annot$Status<-as.numeric(GH_annot$Status)

                 GL_annot$Status<-ifelse(GL_annot$vital_status=="Dead", 1, 0)
                 GL_annot$group<-as.factor(GL_annot$group)
                 GL_annot$time<-as.numeric(GL_annot$time)
                 GL_annot<-na.omit(GL_annot)
                 GL_annot$Status<-as.numeric(GL_annot$Status)

                 data_prepared<-rbind(GH_annot,GL_annot)
                 data_prepared$time<-(data_prepared$time)/30
                 months<-as.numeric(years())*12

                 data_prepared<-data_prepared[which(data_prepared$time <= months),]
                 model_fit<-survfit(Surv(time, Status) ~ group, data=data_prepared)

                 GENE<-ENSG()
                 CANCER<-TCGA_data()
                 p<-ggsurvplot(model_fit, data_prepared,
                 pval = TRUE, conf.int = TRUE,
                 title=paste0(GENE," expression in ",CANCER),
                 legend.title="Expression",
                 xlab="Months",
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                 legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"))

                 output$distPlot <- renderPlot({p})

                 incProgress(0.89)


               }else if(metrics()=="1st quartile and mean"){

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[2])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[4])

                 GROUP_LOW_index<-which(Gene_data < summary(Gene_data)[3])
                 GROUP_HIGH_index<-which(Gene_data >= summary(Gene_data)[3])

                 GROUP_LOW<-Gene_data[GROUP_LOW_index ]
                 GROUP_LOW_df<-as.data.frame(GROUP_LOW)
                 GROUP_LOW_df$sample<-names(GROUP_LOW)
                 names(GROUP_LOW_df)[1]<-"Expression"

                 output_low_df<-data.frame()
                 for(sample in unique(GROUP_LOW_df$sample)){
                   subset_df<-GROUP_LOW_df[which(GROUP_LOW_df$sample ==sample),]
                   kk<-subset_df[1,]
                   output_low_df<-rbind(output_low_df,kk)
                 }
                 rownames(output_low_df)<-paste(output_low_df$sample)

                 GROUP_HIGH<-Gene_data[GROUP_HIGH_index ]
                 GROUP_HIGH_df<-as.data.frame(GROUP_HIGH)
                 GROUP_HIGH_df$sample<-names(GROUP_HIGH)
                 names(GROUP_HIGH_df)[1]<-"Expression"

                 output_high_df<-data.frame()
                 for(sample in unique(GROUP_HIGH_df$sample)){
                   subset_df<-GROUP_HIGH_df[which(GROUP_HIGH_df$sample ==sample),]
                   subset_df<-subset_df[1,]
                   output_high_df<-rbind(output_high_df,subset_df)
                 }
                 rownames(output_high_df)<-paste(output_high_df$sample)

                 clinical<-as.data.frame(clinical)
                 clinical<-clinical[,c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up")]
                 clinical$time<-ifelse(is.na(clinical$days_to_death), clinical$days_to_last_follow_up, clinical$days_to_death)
                 clinical<-clinical[,c("submitter_id","vital_status","time")]
                 clinical<-unique(clinical)

                 output_low_df$sample<-paste(lapply(strsplit(paste(output_low_df$sample),"-01"),"[",1))
                 output_high_df$sample<-paste(lapply(strsplit(paste(output_high_df$sample),"-01"),"[",1))

                 GL_annot<-merge(output_low_df, clinical, by.x="sample", by.y="submitter_id")
                 GL_annot$group<-"Low"
                 GL_annot<-unique(GL_annot)
                 #save(GL_annot, file="/media/inmuno/data/GL_annot.R")

                 GH_annot<-merge(output_high_df, clinical, by.x="sample", by.y="submitter_id")
                 GH_annot$group<-"High"
                 GH_annot<-unique(GH_annot)

                 GH_annot$Status<-ifelse(GH_annot$vital_status=="Dead", 1, 0)
                 GH_annot$group<-as.factor(GH_annot$group)
                 GH_annot$time<-as.numeric(GH_annot$time)
                 GH_annot<-na.omit(GH_annot)
                 GH_annot$Status<-as.numeric(GH_annot$Status)

                 GL_annot$Status<-ifelse(GL_annot$vital_status=="Dead", 1, 0)
                 GL_annot$group<-as.factor(GL_annot$group)/
                 GL_annot$time<-as.numeric(GL_annot$time)
                 GL_annot<-na.omit(GL_annot)
                 GL_annot$Status<-as.numeric(GL_annot$Status)

                 data_prepared<-rbind(GH_annot,GL_annot)
                 data_prepared$time<-(data_prepared$time)/30
                 months<-as.numeric(years())*12

                 data_prepared<-data_prepared[which(data_prepared$time <= months),]
                 model_fit<-survfit(Surv(time, Status) ~ group, data=data_prepared)

                 GENE<-ENSG()
                 CANCER<-TCGA_data()
                 p<-ggsurvplot(model_fit, data_prepared,
                 pval = TRUE, conf.int = TRUE,
                 title=paste0(GENE," expression in ",CANCER),
                 legend.title="Expression",
                 xlab="Months",
                 risk.table = TRUE, # Add risk table
                 risk.table.col = "strata", # Change risk table color by groups
                 linetype = "strata", # Change line type by groups
                 ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                 legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
                 palette = c("#E7B800", "#2E9FDF"))

                 output$distPlot <- renderPlot({p})

                 incProgress(0.89)


               }   # end of if condition


               })   # end of withProgress run analysis
        #aqui
        pdf("test.pdf", onefile=FALSE)
        print(ggsurvplot(model_fit, data_prepared,
        pval = TRUE, conf.int = TRUE,
        title=paste0(GENE," expression in ",CANCER),
        legend.title="Expression",
        xlab="Months",
        risk.table = TRUE, # Add risk table
        risk.table.col = "strata", # Change risk table color by groups
        linetype = "strata", # Change line type by groups
        ggtheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10)), # Change ggplot2 theme
        palette = c("#E7B800", "#2E9FDF"))
        )
        dev.off()


        

        output$download_kaplan<-downloadHandler(
            filename <- function(){'Kaplan-meier.pdf'},

            content = function(file){
              file.copy("test.pdf", file)

            }
        )




       })  # collect data

       session$onSessionEnded(function() {
         system(paste("rm -f ", "test.pdf"))
       })
        wd<-getwd()
  #      tmpdir<-tempdir()
        output$path<-renderText({wd})

}

shinyApp(ui = ui, server = server)
