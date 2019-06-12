#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/
#

#library(DT)
source("helpers.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(HTML('
  <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=UA-18253988-3"></script>
    <script>
    window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag("js", new Date());
  
  gtag("config", "UA-18253988-3");
  </script>')),
  
  shinyjs::useShinyjs(),
  #img(src='images/SingleR_logo.jpg', align = "left"),
  titlePanel(windowTitle='SingleRbrowseR',title=div(img(src="SingleR_logo.jpg"), HTML('<span style="text-decoration: underline;">Single</span>','-cell <span style="text-decoration: underline;">R</span>ecognition of cell types'))),
  tabsetPanel(
    tabPanel(title='Data sets',
             fluidRow(
               column(12,
                      h5('Important note: the SingeR browseR has been upgraded. It now uses a new SingleR object. Please refer to https://github.com/dviraran/SingleR for more information.'),
                column(12,
                      fileInput(inputId="mydata", label = "Import Data (SingleR object)",multiple = TRUE,accept=c('.rds')))
               ),width="100%")
    ),
    tabPanel(title='Analysis',
             fluidRow(
               sidebarLayout(
                 sidebarPanel(
                   selectInput('data1','Data set:',NULL),
                   selectInput('ref_data','Refererence set:',NULL),
                   sliderInput("dot.size", "Point size:", min=1, max=25, value=3),
                   radioButtons('data2','Color by:',analysis_sets, selected = 'SingleR annotations'),
                   checkboxInput('by_clusters','Aggregate by clusters'),
                   
                   conditionalPanel(condition = "input.by_clusters == true",
                                    selectizeInput('data_clusters',"Choose clustering:",
                                                   choices=NULL,
                                                   selected=NULL,
                                                   multiple=F)
                   ),
                   
                   conditionalPanel(condition = "input.data2 == 'Gene'",
                                    selectizeInput('data6',"Choose gene:",
                                                   choices=NULL,
                                                   selected=NULL,
                                                   multiple=T),
                                    helpText('Multiple genes will plot counts of non-zero expression.')
                                    
                   ),
                   
                   conditionalPanel(condition = "input.data2 == 'Meta.data'",
                                    selectizeInput('data7',"Choose meta data:",
                                                   choices=NULL,
                                                   selected=NULL,
                                                   multiple=F)
                   ),
                   
                   conditionalPanel(condition = "input.data2 == 'Other'",
                                    selectizeInput('data8',"Choose other data:",
                                                   choices=NULL,
                                                   selected=NULL,
                                                   multiple=F)
                   ),
                   
                   verbatimTextOutput('about_text'),
                   checkboxInput('show_heatmap','Show scores heatmap?'),
                   helpText('Use the toolbar on the top of the plot to save. Click ESC to reset plot.')
                   
                   
                 ,width=3),
                 mainPanel(
                   canvasXpressOutput("scatter_plot",height=610),
                   shinysky::busyIndicator(wait = 1500)
                 
                  # uiOutput('main_panel') 
                   #tabsetPanel(id = "inTabset",
                   #            tabPanel("Plot",                    
                   #                      canvasXpressOutput("scatter_plot",height=610),
                   #                     shinysky::busyIndicator(wait = 1000)),
                   #            tabPanel("panel2", h2("This is the second panel."))
                   #)
                  ,width=9)),
               fluidRow(
                 conditionalPanel(condition = "input.show_heatmap == true",
                                  column(10,
                                         helpText('Use shift to select multiple cells in the scatter plot to show in the heatmap.'),
                                         plotOutput("response_plot",height=500)
                                         
                                  ), column(2,
                                            #   conditionalPanel('is_clicked',
                                            br(),br(),br(),
                                            sliderInput('levelsShow','Levels show:', min=2, max=200, value=40),
                                            checkboxInput('order_by_cluster','Order by cluster?'),
                                            conditionalPanel(condition = "input.order_by_cluster == true",
                                                             selectizeInput('heatmap_cluster',"Annotate by:",
                                                                            choices=NULL,
                                                                            selected=NULL,
                                                                            multiple=F)
                                            ),
                                            
                                            #actionButton('prevBtn','Previous'),br(),br(),
                                            #actionButton('nextBtn','Next'),
                                            downloadButton('save_heatmap', 'Save heatmap'),
                                            helpText('Use the lasso in the plot tool bar to choose single-cells for the heatmap.'),
                                            
                                            uiOutput("selectedLabels")
                                  )
                 ),
                 
                 conditionalPanel(condition = "input.data2 == 'Gene'",
                                  plotOutput("gene_box_plot2",width='100%',height=400),
                                  
                                  plotOutput("gene_box_plot",width='100%',height=400)
                 )
                 
               )
               
             )
    ),
    
    tabPanel(title='Cluster',
             sidebarLayout(
               sidebarPanel(
                 uiOutput("header_panel5"),
                 selectInput('ref_data_cluster','Refererence set:',NULL),
                 sliderInput("dot.size.cluster", "Point size:", min=1, max=25, value=3),
                 sliderInput("cluster_num", "Number of clusters:", min=2, max=30, value=10),
                 actionButton('recluster_btn',"Cluster", icon("rocket")),
                 downloadButton('downloadClusterPDF', 'Download CSV'),
                 br(),br(),
                 plotOutput("cluster_dend",height=250)
               ),
               mainPanel(
                 canvasXpressOutput("cluster_plot",height=600),
                 shinysky::busyIndicator(wait = 1500))
             )
    ),
    tabPanel(title='Proportions',
             fluidRow(
               sidebarLayout(
                 sidebarPanel(
                   selectInput('ref_data_prop','Refererence set:',NULL),
                   selectInput('data_prop','Group by:',c('Batches','Clusters')),
                   sliderInput("prop.thres", "Thresold (%):", min=0, max=50, value=3)
                   ,width="100%"),
                 
                 mainPanel(column(1),column(10,
                                            plotOutput("prop_plot",height=600),
                                            dataTableOutput('prop_table'),
                                            downloadButton('downloadTable', 'Download table',width='100%'))
                           ,width="100%"))
               
             )),
    
    tabPanel(title='Differential analysis',
             
             fluidRow(
               uiOutput("header_panel4"),
               column(9,offset=1,h4("Differential expression (DE) analysis between two sets.",align='center'),
                      h5("In each set choose cells by using the three types of filters: by SingleR annotations,
                    by identities and by clusters. The bottom box shows the number of cells in the set
                    after accounting for the filters. Set 2 may be left empty, than all cells not in Set 1
                   are considered for Set 2. Once the Sets are chosen, click Run DE analysis. May take a
                   while to run.", align = "center")),
               column(4,wellPanel(
                 h4('Set 1',align = "center"),
                 column(8,selectInput('de_set1_ref_data1','SingleR reference:',NULL)),
                 selectizeInput('de_set1_labels1','','',multiple = TRUE),
                 column(8, selectInput('de_set1_group1',' & Group 1:',c())),
                 selectizeInput('de_set1_labels2','','',multiple = TRUE),
                 column(8, selectInput('de_set1_group2','& Group 2:',c())),
                 selectizeInput('de_set1_labels3','','',multiple = TRUE),
                 verbatimTextOutput('de_set1_num_cells')
               )),
               column(4,wellPanel(
                 h4('Set 2',align = "center"),
                 column(8,selectInput('de_set2_ref_data1','SingleR reference:',NULL)),
                 selectizeInput('de_set2_labels1','','',multiple = TRUE),
                 column(8, selectInput('de_set2_group1',' & Group 1:',c())),
                 selectizeInput('de_set2_labels2','','',multiple = TRUE),
                 column(8, selectInput('de_set2_group2','& Group 2:',c())),
                 selectizeInput('de_set2_labels3','','',multiple = TRUE),
                 verbatimTextOutput('de_set2_num_cells')
               )),
               column(4,  wellPanel(   
                 h4('Parameters',align = "center"),
                 sliderInput("min_diff_pct", "Minimum difference percentage:", min=0, max=100, value=5),
                 sliderInput("min_pct", "Minimum percentage:", min=0, max=100, value=25),
                 sliderInput("de.num.genes", "Num. of genes:", min=1, max=50, value=20),
                 selectInput('data_method','Method:',c('Wilcoxon test','Bimod test','ROC test','t-Test','Tobit-censoring',"Poisson test","Negative-binomial test","MAST"),
                             selected='Wilcoxon test'),
                 downloadButton('downloadData2', 'Download'),
                 actionButton('run2','Run DE analysis')
                 
               )),
               plotOutput("diff_plot",height=800,width='100%'),
               shinysky::busyIndicator(wait = 1500)
             ) 
             
             
             
    )
  )
)


server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=1000*1024^2)
  
  #outputOptions(output, "header_panel2", suspendWhenHidden = FALSE)
  
  obj.subset <- reactiveValues()
  obj <- reactiveValues()
  singler_diff <- reactiveValues()

  clicked_cells <- reactiveValues(count=1)
  
  observeEvent(singler.upload(),{
    
  })
  
  singler.upload = eventReactive(input$mydata,{
    print('mydata')
    data = list()
    for(i in 1:length(input$mydata[,1])){
      data[[i]] = readRDS(input$mydata[[i, 'datapath']])
    }
    
    titles = unlist(lapply(data,function(x) x@project.name))
    updateSelectInput(session, "data1",
                      label = 'Data set:',
                      choices = c(titles),
                      selected = NULL)
    
    return(data)
  })
  
  observeEvent(input$data1, {
    print('data1')
    temp = singler.upload()[[1]]
    obj$singler = temp
    N = nrow(obj$singler@xy)
    updateSliderInput(session,"dot.size", "Point size:", min=1, max=25, value=min(max((18.5-log2(N))*2,1),25))
    updateSliderInput(session,"dot.size.cluster", "Point size:", min=1, max=25, value=min(max((18.5-log2(N))*2,1),25))
    
    isolate({

      updateSelectInput(session, "ref_data",
                        label = 'Reference set:',
                        choices =  colnames(obj$singler@labels),
                        selected =  colnames(obj$singler@labels)[1])
      
      updateSelectInput(session, "ref_data_cluster",
                        label = 'Reference set:',
                        choices =  colnames(obj$singler@labels),
                        selected =  colnames(obj$singler@labels)[1])
      
      updateSelectInput(session, "ref_data_prop",
                        label = 'Reference set:',
                        choices =  colnames(obj$singler@labels),
                        selected =  colnames(obj$singler@labels)[1])
      
      if (length(obj$singler@other)==0) {
        updateRadioButtons(session,'data2','Color by:',analysis_sets[-4], selected = 'SingleR annotations')
      }
      g = sort(rownames(obj$singler@expr))
      obj$gene = g[1]
      updateSelectizeInput(session,inputId='data6','Choose gene:',choices=g,selected=obj$gene,server=T)
      updateSelectizeInput(session,inputId='data7','Choose meta data:',choices=c(colnames(obj$singler@ident),colnames(obj$singler@clusters)),selected=colnames(obj$singler@ident)[1], server = TRUE)
      updateSelectizeInput(session,inputId='data8','Choose other data:',choices=colnames(obj$singler@other),selected=colnames(obj$singler@other)[1], server = TRUE)
      updateSelectizeInput(session,inputId='data_prop','Group by:',choices=c(colnames(obj$singler@ident),colnames(obj$singler@clusters)),selected=colnames(obj$singler@ident)[1], server = TRUE)
      updateSelectizeInput(session,inputId='data_clusters','Choose clustering:',choices=c(colnames(obj$singler@clusters)),selected=colnames(obj$singler@clusters)[1], server = TRUE)
      updateSelectizeInput(session,inputId='heatmap_cluster','Annotate by:',choices=c(colnames(obj$singler@clusters),colnames(obj$singler@ident)),selected=colnames(obj$singler@clusters)[1], server = TRUE)
      
      updateSelectInput(session, "de_set1_ref_data1",label = 'SingleR reference:',choices =  colnames(obj$singler@labels),selected =  colnames(obj$singler@labels)[1])
      updateSelectInput(session, "de_set1_group1",label = '& Group 1:',choices = colnames(obj$singler@ident),selected=colnames(obj$singler@ident)[1])
      updateSelectInput(session, "de_set1_group2",label = '& Group 2:',choices = colnames(obj$singler@clusters),selected=colnames(obj$singler@clusters)[1])
      
      updateSelectInput(session, "de_set2_ref_data1",label = 'SingleR reference:',choices =  colnames(obj$singler@labels),selected =  colnames(obj$singler@labels)[1])
      updateSelectInput(session, "de_set2_group1",label = '& Group 1:',choices = colnames(obj$singler@ident),selected=colnames(obj$singler@ident)[1])
      updateSelectInput(session, "de_set2_group2",label = '& Group 2:',choices = colnames(obj$singler@clusters),selected=colnames(obj$singler@clusters)[1])
      
      
      updateCheckboxInput(session,inputId='show_heatmap','Show scores heatmap?',value=F)
      
      if (length(levels(obj$singler@clusters))==1) {
        shinyjs::disable("by_clusters")
      } else {
        shinyjs::enable("by_clusters")
      }
      shinyjs::disable("prevBtn")
      shinyjs::disable("nextBtn")
      
    })
    
    if(input$de_set1_ref_data1 %in% colnames(obj$singler@labels)) {
      if(input$de_set1_ref_data1!="")
        updateSelectizeInput(session, "de_set1_labels1",label = '',choices =  sort(unique(obj$singler@labels[,input$de_set1_ref_data1])),selected = NULL)
    }
    print(input$de_set1_group1)
    
    if(input$de_set1_group1 %in% colnames(obj$singler@ident))
      updateSelectizeInput(session, "de_set1_labels2",label = '',choices =  sort(unique(obj$singler@ident[,input$de_set1_group1])),selected = NULL)
    if(input$de_set1_group2 %in% colnames(obj$singler@clusters))
      updateSelectizeInput(session, "de_set1_labels3",label = '',choices =  sort(unique(obj$singler@clusters[,input$de_set1_group2])),selected = NULL)
    if(input$de_set2_ref_data1 %in% colnames(obj$singler@labels)) {
      if(input$de_set2_ref_data1!="")
        updateSelectizeInput(session, "de_set2_labels1",label = '',choices =  sort(unique(obj$singler@labels[,input$de_set2_ref_data1])),selected = NULL)
    }
    if(input$de_set2_group1 %in% colnames(obj$singler@ident))
      updateSelectizeInput(session, "de_set2_labels2",label = '',choices =  sort(unique(obj$singler@ident[,input$de_set2_group1])),selected = NULL)
    if(input$de_set2_group2 %in% colnames(obj$singler@clusters))
      updateSelectizeInput(session, "de_set2_labels3",label = '',choices =  sort(unique(obj$singler@clusters[,input$de_set2_group2])),selected = NULL)
    
    
  }) 
  
  output$main_panel = renderUI({
    if (FALSE) {
      if (input$data2 == 'Gene' && length(input$data6)>1) {
        output$mytabs = renderUI({
          nTabs = length(input$data6)
          myTabs = lapply(seq_len(nTabs), function(i) {
            tabPanel(input$data6[i],
                     canvasXpressOutput(paste0('scatter_plot',i), width = "100%", height = "610px")
            )
          })
          do.call(tabsetPanel, myTabs)
        })
        
      } else {
        canvasXpressOutput('scatter_plot', width = "100%", height = "610px")
      }
    }
    
    canvasXpressOutput('scatter_plot', width = "100%", height = "610px")
    
  })
  
  output$scatter_plot1 <- renderCanvasXpress({
    print('scatter_plot1')
    createScatterPlot(obj,input,input$data6[1])
  })
  
  output$scatter_plot2 <- renderCanvasXpress({
    print('scatter_plot2')
    createScatterPlot(obj,input,input$data6[2])
  })
  
  output$scatter_plot <- renderCanvasXpress({
    print('scatter_plot')
    createScatterPlot(obj,input,NULL)
  })
  
  output$selectedLabels <- renderUI({
    selectInput("selected", "", rownames(obj$singler@labels), selectize = FALSE, multiple = TRUE, selected = rownames(obj$singler@labels))
  })
  
  output$response_plot <- renderPlot({
    print('response_plot')
    sel = input$selected
    order_by_cluster = input$order_by_cluster
    heatmap_cluster = input$heatmap_cluster
    n.show = input$levelsShow
    isolate({
      if (length(sel)>1) {
        cluster.by = NULL
        if (order_by_cluster == T) {
          print(heatmap_cluster)
          if (heatmap_cluster %in% colnames(obj$singler@clusters)) {
            cluster.by = obj$singler@clusters[,heatmap_cluster,drop=F]
          } else {
            cluster.by = obj$singler@ident[,heatmap_cluster,drop=F]
          }
        }
        annot = cbind(obj$singler@ident,obj$singler@clusters)
        colnames(annot) = c(colnames(obj$singler@ident),colnames(obj$singler@clusters))
        rownames(annot) = rownames(obj$singler@ident)
        obj$filename = paste0(obj$singler@project.name,'_heatmap.pdf')
        obj$plot = singlerDrawHeatmap(obj$singler@scores[[input$ref_data]],sel,annot=annot,cluster.by=cluster.by,n.show = n.show)
      } else {
        
        #singlerDrawBoxPlot2(scores,obj$reference,clicked_cells$cell_id,n.show = input$levelsShow,is.main.types=input$ref_data_main)
        
        if (FALSE) {
          if (vNext$doPlot == TRUE) {
            clicked_cells$top_labels[[clicked_cells$count+1]] = FineTuningStep(obj$singler.use,obj$expr,obj$reference,types,clicked_cells$top_labels[clicked_cells$count][[1]]$labels,clicked_cells$cell_id)
            if (length(clicked_cells$top_labels[clicked_cells$count][[1]]$labels)>1) {
              singlerDrawBoxPlot(obj$singler.use, obj$reference,types, clicked_cells$top_labels[clicked_cells$count][[1]]$r,clicked_cells$cell_id,labels.use=clicked_cells$top_labels[clicked_cells$count][[1]]$labels,n.show = input$levelsShow)
            }
            vNext$doPlot = FALSE
            
          }
          if (vPrev$doPlot == TRUE) {
            if (clicked_cells$count>1) {
              singlerDrawBoxPlot(obj$singler.use, obj$reference,types, clicked_cells$top_labels[clicked_cells$count][[1]]$r,clicked_cells$cell_id, labels.use=clicked_cells$top_labels[clicked_cells$count][[1]]$labels,n.show = input$levelsShow)
            }
            vPrev$doPlot = FALSE
          }
        }
      }
    })
  })
  
  output$save_heatmap <- downloadHandler(
    filename = function() {
      obj$filename
    },
    content = function(file) {
      save_pheatmap_pdf(obj$plot,file ,width=10,height=8)
    }
  )
  
  output$about_text <- renderText({
    paste(apply(cbind(names(obj$singler@meta.data),obj$singler@meta.data),1,
                FUN=function(x) paste0(x[1],": ",x[2])),collapse = '\n')
  })
  
  output$prop_plot <- renderPlot({
    print('prop_plot')
    if (input$data_prop %in% colnames(obj$singler@clusters)) {
      tbl = table(obj$singler@clusters[,input$data_prop],obj$singler@labels[,input$ref_data_prop])
    } else {
      tbl = table(obj$singler@ident[,input$data_prop],obj$singler@labels[,input$ref_data_prop])
      
    }
    
    obj$table = tbl
    tbl = tbl/rowSums(tbl)
    obj$filename = paste0('SingleR_prop_',obj$singler@project.name)
    
    A = colMaxs(tbl)>input$prop.thres/100
    df = melt(tbl[,A])
    colnames(df)= c('Group','Cell_types','Percentage')
    corrplot(method='color',tbl[,A],is.corr=F,tl.srt=45,number.cex=1,addCoef.col = "gray")
    #ggplot(data=df, aes(x=Cell_types, y=Percentage, fill=Group))+geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  })
  
  output$prop_table <- renderDataTable({
    print('prop_table')
    #isolate({
    tbl = obj$table
    tbl = cbind(rownames(tbl),tbl,apply(tbl,1,sum))
    colnames(tbl)[dim(tbl)[2]]='Sum'
    tbl
    #})
  })
  
  
  observeEvent(input$recluster_btn, {
    print('recluster_btn')
    isolate({
      scores = t(scale(t(obj$singler@scores[[input$ref_data]]^3)))
      hc = hclust(dist(scores,method='euclidean'),method='ward.D2')
      cl = cutree(hc,k=input$cluster_num)
      obj$singler.clusters  = list(hc=hc,cl=factor(cl))
      #obj$singler@clusters$`SingleR cluster` = obj$singler.clusters$cl 
      #print(obj$singler@clusters)
      #updateSelectizeInput(session,inputId='data_clusters','Choose clustering:',choices=c(colnames(obj$singler@clusters)),selected=colnames(obj$singler@clusters)[1], server = TRUE)
      
    })
  })
  
  output$downloadClusterPDF <- downloadHandler(
    filename = function() {
      paste0(obj$filename,'_clustering.csv')
    },
    content = function(file) {
      write.csv(obj$df_cluster, file)
    }
  )
  
  output$cluster_dend <- renderPlot({
    print('cluster_dend')
    cluster_num = input$cluster_num
    if(!is.null(obj$singler.clusters)) {
      isolate({
        dend = as.dendrogram(obj$singler.clusters$hc)
        dend =  color_branches(dend,input$cluster_num,col=colors[1:input$cluster_num])
        #branches_color(dend) = colors[obj$singler[[objAll$id]]$SingleR.single$clusters$cl][order.dendrogram(dend)]
        dend = set(dend,'branches_lwd',4)
        d = dend
        obj$singler.clusters$cl = cutree(d,cluster_num)
        n = length(labels(d))
        labels(d) = rep('',n)
        par(mar=c(1,2,1,0.5))
        plot(d)
      })
    }
    
  })
  
  output$cluster_plot <- renderCanvasXpress({
    print('cluster_plot')
    shinyjs::enable("recluster_btn")
    shinyjs::enable("cluster_num")
    
    if(!is.null(obj$singler.clusters)) {
      clusters = data.frame(Clusters=obj$singler.clusters$cl)
      dot.size = input$dot.size.cluster
      canvasXpress(obj$singler@xy,varAnnot=clusters,graphType='Scatter2D',
                   colorBy='Clusters',dataPointSize=dot.size,outlineWidth=dot.size/40,
                   stringVariableFactors=list('Clusters'),sizeByShowLegend=F,xAxisMinorTicks=F,yAxisMinorTicks=F,
                   xAxisMajorTicks=F,yAxisMajorTicks=F,xAxisTitle='tSNE1',yAxisTitle='tSNE2',
                   colors=colors[1:input$cluster_num])
      
    }
    #p
  })
  
  output$gene_box_plot <- renderPlot({
    print('gene_box_plot')
    if (input$data2!='Gene') return()
    ref.name = gsub('\\.main','',tolower(input$ref_data))
    if (ref.name=='mouse-rnaseq' || ref.name=='mousernaseq') {
      ref.name = 'mouse.rnaseq'
    }
    ref = eval(parse(text = ref.name))
    if (!is.null(ref)) {
      if (is.null(input$data6[1])) {
        gene = obj$gene
      } else {
        gene = input$data6
      }
      if (sum(gene %in% rownames(ref$data))==length(gene)) {
        g = ref$data[gene,]
        
        if (length(gene)>1) {
          df = data.frame(t(g),main_types=(ref$main_types))
          df = melt(df,'main_types')
          gene = paste0(gene,collapse='+')
          
          ggplot(df,aes(x=main_types, y=value)) + geom_violin(aes(color=main_types,fill=main_types)) + 
            geom_boxplot(width=0.1, outlier.size = 0) + geom_jitter(shape=16, size=1,color='grey') + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(paste(gene,'in reference data set'))+
            ylab('Expression')+theme(legend.position="none")+facet_wrap(~variable, scale="free")
        }  else {
          df = data.frame(x=g,types=(ref$main_types),main_types=(ref$main_types))
          colnames(df) = c('x','types','main_types')
          ggplot(df,aes(x=main_types, y=x)) + geom_violin(aes(color=main_types,fill=main_types)) + 
            geom_boxplot(width=0.1, outlier.size = 0) + geom_jitter(shape=16, size=1,color='grey') + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(paste(gene,'in reference data set'))+
            ylab('Expression')+theme(legend.position="none")
        }
      }
    }
  })
  
  
  output$gene_box_plot2 <- renderPlot({
    print('gene_box_plot2')
    if (input$data2!='Gene') return()
    
    if (is.null(input$data6[1])) {
      gene = obj$gene
    } else {
      gene = as.character(input$data6)
    }
    g = obj$singler@expr[gene,]
    
    if (!empty(obj$singler@labels.clusters)) {
      clusters = clusters.map.values(obj$singler@clusters[,input$data_clusters],obj$singler@labels.clusters[,input$ref_data,drop=F])
    } else {
      clusters = obj$singler@clusters[,input$data_clusters]
    }
    
    if (length(gene)>1) {
      df = data.frame(x=colSums(g>0),Clusters=clusters)
      colnames(df) = c('x','Clusters')
      df$x = factor(df$x)
      gene = paste0(gene,collapse='+')
      
      ggplot(df,aes(x=Clusters, y=x)) + geom_col(aes(fill=x)) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(paste(gene,'in clusters'))+
        ylab('Expression')
    }  else {
      
      df = data.frame(x=g,Clusters=clusters)
      colnames(df) = c('x','Clusters')
      ggplot(df,aes(x=Clusters, y=x)) + geom_violin(aes(color=Clusters,fill=Clusters)) + 
        geom_boxplot(width=0.1, outlier.size = 0) + geom_jitter(shape=16, size=0.5,color='grey') + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle(paste(gene,'in clusters'))+
        ylab('Expression')+theme(legend.position="none")
    }
  })
  
  
  
  #------------- diff tab -----------#
  
  observeEvent(input$de_set1_ref_data1, {
    print('de_set1_ref_data1s')
    if (input$de_set1_ref_data1=="") return()
    updateSelectizeInput(session, "de_set1_labels1",label = '',choices =  sort(unique(obj$singler@labels[,input$de_set1_ref_data1])),selected = NULL)
  })
  
  observeEvent(input$de_set1_group1, {
    if (input$de_set1_group1=="") return()
    updateSelectizeInput(session, "de_set1_labels2",label = '',choices =  sort(unique(obj$singler@ident[,input$de_set1_group1])),selected = NULL)
  })
  
  observeEvent(input$de_set1_group2, {
    if (input$de_set1_group2=="") return()
    updateSelectizeInput(session, "de_set1_labels3",label = '',choices =  sort(unique(obj$singler@clusters[,input$de_set1_group2])),selected = NULL)
  })
  
  observeEvent(input$de_set2_ref_data1, {
    if (input$de_set2_ref_data1=="") return()
    updateSelectizeInput(session, "de_set2_labels1",label = '',choices =  sort(unique(obj$singler@labels[,input$de_set2_ref_data1])),selected = NULL)
  })
  
  observeEvent(input$de_set2_group1, {
    if (input$de_set2_group1=="") return()
    print(input$de_set2_group1)
    updateSelectizeInput(session, "de_set2_labels2",label = '',choices =  sort(unique(obj$singler@ident[,input$de_set2_group1])),selected = NULL)
  })
  
  observeEvent(input$de_set2_group2, {
    if (input$de_set2_group2=="") return()
    updateSelectizeInput(session, "de_set2_labels3",label = '',choices =  sort(unique(obj$singler@clusters[,input$de_set2_group2])),selected = NULL)
  })
  
  output$de_set1_num_cells <- renderText({
    A = rep(TRUE,nrow(obj$singler@labels)); B = A; C = A
    
    if (length(input$de_set1_labels1)>0)
      A = obj$singler@labels[,input$de_set1_ref_data1] %in% input$de_set1_labels1
    if (length(input$de_set1_labels2)>0)
      B =  obj$singler@ident[,input$de_set1_group1] %in% input$de_set1_labels2
    if (length(input$de_set1_labels3)>0)
      C =  obj$singler@clusters[,input$de_set1_group2] %in% input$de_set1_labels3
    
    obj$de.set1 = A&B&C
    sprintf('Set 1 number of cells: %d', sum(A&B&C))
  })
  
  output$de_set2_num_cells <- renderText({
    A = rep(TRUE,nrow(obj$singler@labels)); B = A; C = A
    
    if (length(input$de_set2_labels1)>0)
      A = obj$singler@labels[,input$de_set2_ref_data1] %in% input$de_set2_labels1
    if (length(input$de_set2_labels2)>0)
      B =  obj$singler@ident[,input$de_set2_group1] %in% input$de_set2_labels2
    if (length(input$de_set2_labels3)>0)
      C =  obj$singler@clusters[,input$de_set2_group2] %in% input$de_set2_labels3
    
    obj$de.set2 = A&B&C
    sprintf('Set 1 number of cells: %d', sum(A&B&C))
  })
  
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$run2, {
    v$doPlot <- input$run2
  })
  
  output$diff_plot <- renderPlot({
    print('diff_plot')
    if (v$doPlot == FALSE) return()
    isolate({
      obj$dm = Identify.Diff.Markers(obj$singler@expr,obj$de.set1,obj$de.set2,input$data_method,input$min_pct,input$min_diff_pct)
      shinyjs::enable("downloadData2")
      top.genes = NULL
      #  if(is.null(input$obj$de.set1)) {
      #    obj$dm$data %>% group_by(cluster) %>% top_n(ceiling(input$de.num.genes/length(unique(dm$data$cluster))), avg_diff) -> top.genes
      #  }
      Make.HeatMap(obj$singler@expr,obj$dm,input$de.num.genes,top.genes$gene)
      
    })
  })
  
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0('SingleR_Diff_',obj$singler@project.name,'.csv')
    },
    content = function(file) {
      write.csv(obj$dm$data, file)
    }
  )
  
}




shinyApp(ui = ui, server = server)
