library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(sodium)
#library(readxl)
library(shinythemes)


# Main login screen
loginpage <- div(id = "loginpage", style = "width: 500px; max-width: 100%; margin: 0 auto; padding: 20px;",
                 wellPanel(
                   tags$h2("LOGIN", class = "text-center", style = "padding-top: 0;color:#333; font-weight:600;"),
                   textInput("userName", placeholder="Username", label = tagList(icon("user"), "Username")),
                   passwordInput("passwd", placeholder="Password", label = tagList(icon("unlock-alt"), "Password")),
                   br(),
                   div(
                     style = "text-align: center;",
                     actionButton("login", "SIGN IN", style = "color: white; background-color:#3c8dbc;
                                 padding: 10px 15px; width: 150px; cursor: pointer;
                                 font-size: 18px; font-weight: 600;"),
                     shinyjs::hidden(
                       div(id = "nomatch",
                           tags$p("Oops! Incorrect username or password!",
                                  style = "color: red; font-weight: 600; 
                                            padding-top: 5px;font-size:16px;", 
                                  class = "text-center"))),
                     br(),
                     br(),
                     tags$code("Username: myuser  Password: mypass"),
                     br(),
                     tags$code("Username: myuser1  Password: mypass1")
                   ))
)

credentials = data.frame(
  username_id = c("myuser", "myuser1"),
  passod   = sapply(c("mypass", "mypass1"),password_store),
  permission  = c("basic", "advanced"), 
  stringsAsFactors = F
)

header <- dashboardHeader( title = "Programa de Inmunología e Inmunoterapia", uiOutput("logoutbtn"),titleWidth=450 )

sidebar <- dashboardSidebar(uiOutput("sidebarpanel"),             tags$head(
            tags$style(HTML("
                        @import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
                        body {
                              background-color: white;
                              color: black;
                            }
                            .shiny-input-container {
                     color: #474747;
                                }"))

            )
) 
body <- dashboardBody(shinyjs::useShinyjs(), uiOutput("body"))
ui<-dashboardPage(header, sidebar, body, skin = "yellow")

server <- function(input, output, session) {
  
  login = FALSE
  USER <- reactiveValues(login = login)
  options(shiny.maxRequestSize=90000000*1024^2)
    #session$onSessionEnded(stopApp)

    dataInput <-reactive({get(load(paste0("/media/inmuno/data_projects/",input$dataset,".Rda")))})
    #dataInput <-reactive(get(load(paste0("/Users/jgonzalez.69/Documents/WEB_APP/",input$dataset,".Rda"))))


    #session$onSessionEnded(stopApp)


    observe({ 
    if (USER$login == FALSE) {
      if (!is.null(input$login)) {
        if (input$login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          if(length(which(credentials$username_id==Username))==1) { 
            pasmatch  <- credentials["passod"][which(credentials$username_id==Username),]
            pasverify <- password_verify(pasmatch, Password)
            if(pasverify) {
              USER$login <- TRUE
            } else {
              shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
            }
          } else {
            shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
            shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
          }
        } 
      }
    }    
  })
  
   output$logoutbtn <- renderUI({
    req(USER$login)
    actionButton("logout", "Logout",
    style = "background-color:#0000ff;
                      color:#FFFFFF;
                      border-color:#BEBEBE;
                      border-style:solid;
                      border-width:3px;
                      border-radius:14%;
                      font-size: 70%;" 
   , icon= icon("fas fa-sign-out-alt") , lib = "font-awesome" )    
  })

  
  
  output$sidebarpanel <- renderUI({
    if (USER$login == TRUE ){ 
      sidebarMenu(  # aqui meter el menu de pasar de páginas 
               nav_links <- tags$ul(
            tags$li(
                tags$a(href = "/", "Welcome"),
            ),
            tags$li(
                tags$a(href = "/preprocess/", "Pre-process RNA-Seq PE data"),
            ),
            tags$li(
                tags$a(href = "/ideal/", "Perform RNA-Seq DEg analysis"),
            ),
            tags$li(
                tags$a(href = "/gating", "Automated gating FC "),
            ),
            tags$li(
                tags$a(href = "/prognosis/", "Prognosis Marker analysis"),
            ),
            tags$li(
                tags$a(href = "/tsne/", "tSNE analysis"),
            ),


            ),
            sidebarPanel(width=10,
                 uiOutput('dataset')), 
              #    selectInput("dataset", "Choose a Project:",
              #    #choices = c("Intrust", "Linterna")),
              #     choices= "", selected="Intrust"),  
              # style="color:#253529"), # sidebarPanel  #links navegadores
              # br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              br(),
              box(
                width=10,
                title="Download Project", solidHeader=TRUE, style="font-size: 50%",
                downloadButton("download", "Download .tsv",
                style = "background-color:#13963F; color:#FFFFFF" )
                ),
               br(),
               br(),
               br(),
               br()
              #  box(  # INFO ACTUALIZATION
              #   width = 12,
              #   title = "Modify Project data", status = "danger", solidHeader = TRUE,                   
              #   fileInput("upload", "Select xls file", width = 350, buttonLabel = "Upload...", placeholder="Only .xlsx files", multiple = FALSE, accept=".xlsx"),
              #   actionButton("goButtonActu", "Actualize data!", style = "background-color:#13963F; color:#FFFFFF" ), 
              #   tableOutput("files"),
              #   verbatimTextOutput("name_actualization")
   
              #   )

              #  box(  # INFO ACTUALIZATION
              #   width = 12,
              #   title = "Create New Project", status = "warning", solidHeader = TRUE,                   
              #   fileInput("upload_new", "Project xls file", width = 350, buttonLabel = "Upload...", placeholder="Only .xlsx files", multiple = FALSE, accept=".xlsx"),
              #   textInput("nameproject","Enter project name"),
              #   actionButton("goButtonCreate", "Create Project!", style = "background-color:#13963F; color:#FFFFFF" ), 
              #   #tableOutput("files")   
              #   )

                

                

 
 
      )
    }
  })
  
 output$body <- renderUI({
    if (USER$login == TRUE ) {

    br()
    br()
    br()
    br()
    br()
    br()
    br()
    br()
    br()
    br()
    br()
    br()
    br()


      tabsetPanel(
      tabPanel("Main Information",
      
      # First tab
            tabItem(tabName ="dashboard", class = "active",
                    fluidRow(
                      box(width = 12, dataTableOutput('mytable')),
                      # box(  # INFO ACTUALIZATION
                      # width = 12,
                      # title = "Modify Project data", status = "danger", solidHeader = TRUE,                   
                      # fileInput("upload", "Select xls file", width = 350, buttonLabel = "Upload...", placeholder="Only .xlsx files", multiple = FALSE, accept=".xlsx"),
                      # actionButton("goButtonActu", "Actualize data!", style = "background-color:#13963F; color:#FFFFFF" ), 
                      # tableOutput("files")   
                      # ),

                      
                      box(  # INFO ACTUALIZATION
                      width = 12,
                      title = "Create New Project", status = "warning", solidHeader = TRUE,
                      h5("Please, make sure that the Project name and the file have the same name. Example: 'Foo.xlsx' for the name file and 'Foo' for the `project you want to create."),                   
                      fileInput("upload_new", "Project xlsx file", width = 350, buttonLabel = "Upload...", placeholder="Only .xlsx files", multiple = FALSE, accept=".xlsx"),
                      textInput("nameproject","Enter project name"),
                      actionButton("goButtonCreate", "Create Project!", style = "background-color:#13963F; color:#FFFFFF", icon= icon("fas fa-running") , lib = "font-awesome"  ),
                      tableOutput("files2"),
                      verbatimTextOutput("a"),
                      verbatimTextOutput("b"),
                      #verbatimTextOutput("chorizo")
                      

                      ),
                      box(  # INFO MODIFICATION
                      width = 12,
                      title = "Modify Project data", status = "danger", solidHeader = TRUE,                   
                      fileInput("upload", "Select xlsx file", width = 350, buttonLabel = "Upload...", placeholder="Only .xlsx files", multiple = FALSE, accept=".xlsx"),
                      actionButton("goButtonActu", "Actualize data!", style = "background-color:#13963F; color:#FFFFFF", icon= icon("fas fa-running") , lib = "font-awesome" ), 
                      h5("You are going to modify the following project:"),
                      verbatimTextOutput("name_actualization"),
                      tableOutput("files"),
                      )


                    ))
         
      ),# fin del primer tab
      tabPanel("Visualization", h2("under construction")), 
    ) 
    }
    else {
      loginpage
    }
  })
  
  output$mytable <-  DT::renderDataTable({
    dataInput()
  })
  
    #   output$mytable = DT::renderDataTable({
    	    
    # dataInput()
    # })
    output$download <- downloadHandler(
    filename = function() {
      paste0(input$dataset, ".tsv")
    },
    content = function(file) {
     	    
      vroom::vroom_write(dataInput(), file)
    }) # page 1

    output$download <- downloadHandler(
    filename = function() {
      paste0(input$dataset, ".tsv")
    },
    content = function(file) {
     	    
      vroom::vroom_write(dataInput(), file)
    }) # page 1

    data <- reactive({
        req(input$upload)
    })   # Modification file

    datab <- reactive({
        req(input$upload_new)
    })   # New project file

    nameproject <- reactive({
        req(input$nameproject)
    })   # Actualization file

    name <- reactive({
        req(input$dataset)
    })   # Actualization file


    output$files<-renderTable(data())

    output$files2<-renderTable(datab())

    output$name_actualization<-renderText({name()})

    myprojects = list.dirs("/media/inmuno/data_projects", recursive=FALSE)
    myprojects<-basename(myprojects)
    myprojects<-myprojects[! myprojects %in% c("Independent")]
      #myprojects<-as.list(myprojects)
     
    
     output$dataset <- renderUI({
      selectInput(inputId = 'dataset', label = 'choose project',
                  choices = myprojects)
    })


    name <- reactive({
        req(input$dataset)
    })   # Actualization file


    observeEvent(input$goButtonActu, {
    validate(
        need(data()$datapath !="", "Please upload your files before!")
    )
    
    #browser()
    fileIN<-read_xlsx(data()$datapath, sheet=1, col_names=T, na="",skip=0)
    #browser()
    save(fileIN, file=paste0("/media/inmuno/data_projects/",name(),".Rda"))
    }
    )
    a<-reactive({paste0("mkdir ",paste0("/media/inmuno/data_projects/",nameproject(),"/"))})
    b<-reactive({paste0('rsync -a -f"+ */" -f"- *" /media/inmuno/data_projects/Intrust/ ', paste0("/media/inmuno/data_projects/",nameproject(),"/"))})
    
    output$a<-renderText({b()})
    output$b<-renderText({a()})
    observeEvent(input$goButtonCreate, {
    validate(
        need(datab()$datapath !="", "Please upload your files before!")
    )
    
    
    #browser()
    system(a())
    
    system(b())
    fileINP<-read_xlsx(datab()$datapath, sheet=1, col_names=T, na="",skip=0)
    #browser()
    save(fileINP, file=paste0("/media/inmuno/data_projects/",nameproject(),".Rda"))
    }
    )

  observeEvent(input$logout, {
    
      session$reload()
})

}

shinyApp(ui = ui, server = server)
