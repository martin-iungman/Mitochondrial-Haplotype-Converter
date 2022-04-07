#cargo paquetes (o instalo si no estan)
packages<-c("shiny", "tidyverse","readxl","writexl","XML", "shinyjs")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

bold = function(x) strong(x, .noWS = "outside")
ital = function(x) em(x, .noWS = "outside")
link = function(s, href = s) a(s, href = href, .noWS = "outside")

registered_mutations<-read_delim("Mutaciones_mito.csv", delim=";")%>%names()

#Funcion para transformar una mutacion
transf_indiv=function(ht, type=type){
  #crea vector vacio
  vector<-c()
  #como el seqscape incluye a veces muchas mutaciones dentro de un unico elemento, hay que identificar cada uno y analizarlo separado:
  #si hay una sola (por ej 73A>G), i solo va a tomar el valor 1. 
  for (i in seq(1,str_length(str_extract(ht, "(?<=(ins|>|del))[a-zA-Z]+")))) {
    #me quedo con la posicion
    pos1<-str_extract(ht, "[:digit:]+")
    #la transformo en numerico. si es unica, 290+1-1=290. si es doble, 290+2-1=291 (para la segunda posicion)
    pos<-as.numeric(pos1)+i-1
    #cual es la base original?? 
    orig<-str_extract(ht, "(?<=([:digit:]|del))[:upper:]+")%>%str_sub(start=i, end=i)
    #cual es la mutada??
    mut<-str_extract(ht, "(?<=(>|ins))[a-zA-Z]+")%>%str_sub(start=i, end=i)%>%toupper()
    #si se selecciona opcion "Nueva"
    if (type=="Nueva") {
      #Si hay una insersion:
      if (str_count(ht,"ins")==1) {
        mutacion<-str_glue("-",pos1,".",i,mut)
      }
      #Si hay delecion:
      if (str_count(ht,"del")==1) {
        mutacion<-str_glue(orig, pos, "DEL")
      }
      #Si hay sustitucion (>)
      if(str_count(ht,">")==1) {
        mutacion<-str_glue(orig,pos,mut)
      }
      
    } else { #para la opcion "Old"
      if (length(grep("ins",ht))==1) {
        mutacion<-str_glue(pos1,".",i,mut)
      }
      if (length(grep("del",ht))==1) {
        mutacion<-str_glue("DEL", pos, orig)
      }
      if(str_count(ht,">")==1) {
        mutacion<-str_glue(pos,mut)
      }
    }
    #agregar al vector la mutacion
    vector[length(vector)+1]<-mutacion
  }
  return(vector)
  
}


#Funcion para transformar todo el haplotipo (predeterminado: nomenclatura nueva)
transf_ht<-function(datapath, type){
  # si el archivo que se sube es un .xlsx o si es un .xlm
  if(grepl("xlsx$",datapath)){
    df<-read_xlsx(datapath)
  } else if (grepl("xml$",datapath)){
    df<-xmlParse(datapath)%>%xmlToDataFrame()
  }
  # separar el id en numero de muestra (sampleid) y rango 
  id<-df$Specimen[!is.na(df$Specimen)]%>%unique()%>%strsplit(.,split="_")
  sampleid<-id[[1]][1]
  rango<-id[[1]][2]%>%gsub("\\(|\\)","",.)
  # Selecciono las mutaciones
  ht<-df$Base_Change[!is.na(df$Base_Change)]
  if(!is_empty(grep("\\[", ht))){ht=ht[-grep("\\[", ht)]}
  #aplico la transformacion a todas las mutaciones y genero vector
  vector=unlist(map(ht, ~transf_indiv(.x,type)))
  #armo un dataframe de salida
  output<-as_tibble_row(c(Muestra=sampleid, Rango=rango, mut=vector))
  return(output)
}

check_mutations<-function(vector){
  unidentified<-vector[!vector%in%registered_mutations]
  unidentified<-unidentified[!grepl("-(309|573).",unidentified)]
  if(!is_empty(unidentified)){
    print(paste0("Atencion! Mutaciones no registradas en la base de datos: ",paste0(unidentified, collapse=", ")))}else{
      print("OK!")}
}

ui <- fluidPage(
  useShinyjs(),  # Set up shinyjs
  
  tags$head(
    tags$style(type = "text/css", "
      .body {font-size: small}
      .well {padding-top: 10px;}
      .selectize-dropdown {width: 250px !important;}
      .fa-check { font-size:xx-large; color:Lime}
      
  ")),
  
  # Application title
  h2(id = "title-h2", bold("SeqConverter:  "), "Conversion de haplotipos desde reporte de SeqScape"),
  tags$style(HTML("#title-h2 {background-color: gray; color: white; padding: 15px}")),
  
  p(bold("Objetivo: "),
    "Extraer los haplotipos de ADN mitocondrial del reporte de SeqScape, transformando a las nomenclaturas aceptadas internacionalmente de ADNmt."),
  
  p("En primer lugar exportar el/los reporte/s generado en SeqScape como .XML. En caso de tratarse de la region control completa, y por lo tanto, haber sido generada en dos proyectos distintos, subir ambos archivos a la aplicacion"),
  p("Automaticamente se generara una tabla de Excel descargable que se mostrara en pantalla, con el haplotipo de la muestra en la nomenclatura que se le indique"),
  p("La nomenclatura referida como 'Nueva' refiere a las recomendaciones de ", link("Parson et al. (2014)", "https://pubmed.ncbi.nlm.nih.gov/25117402/")),
  p("La referida como vieja es la que requiere el sistema de comparacion actual del BNDG"),
  
  p(bold("Chequear mutaciones:  "), "Esta funcion notificara si alguna de las mutaciones no estuviera registrada previamente en el archivo adjunto a la app, realizado a partir de 2700 muestras del BNDG. Para su uso correcto, la nomenclatura debe figurar como Nueva"),
  
  sidebarPanel(
  #Inputs"
  fileInput("file1", "Subí un archivo", accept = c(".xlsx","xls",".xml")),
  fileInput("file2", "Subí otro archivo", accept = c(".xlsx","xls",".xml")),
  selectInput("type", "Convertir a nomenclatura:", choices = c("Nueva", "Vieja")),
  actionButton("check", "Chequear Mutaciones"),
  downloadButton("download")),
  mainPanel(
  tableOutput("tabla")
))

server <- function(input, output, session) {
  #Transformar haplotipos de las tablas (o si hay una sola, solo esa), y pegarlas en dos filas de una tabla  
  data<-reactive({
    req(input$file1)
    row1<-transf_ht(input$file1$datapath, input$type)
    if(!is_empty(input$file2)){
      row2<-transf_ht(input$file2$datapath,input$type)
      bind_rows(row1,row2)
    }else{
      row1
    }
  })
  

  
  observeEvent(input$check, 
              #showNotification("HOLA")
              
              #showNotification( check_mutations(unlist(data())), type="warning")
              #showNotification(paste(registered_mutations, collapse=","))
              #check_mutations(as.character(unlist(data()))[-c(1,2)])
              showModal(modalDialog(check_mutations(as.character(unlist(data()[,-c(1,2)])[!is.na(unlist(data()[,-c(1,2)]))]))))
  )
  
  #mostrar tabla
  output$tabla<-renderTable({data()})
  #descargar en xlsx. El nombre predeterminado es el nombre del file1 junto al type
  output$download <- downloadHandler(
    filename = function() {
      paste0(str_replace(input$file1, "\\..*",""),"_", input$type,".xlsx")
    },
    content = function(file) {
      write_xlsx(data(), file)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
