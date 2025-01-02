#increasing max file size upload
options(shiny.maxRequestSize = 30 * 1024**2)


# Importing required libraries
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this
library(tidyverse)
library(patchwork)
library(pheatmap)
library(DT)
library(DESeq2)
library(biomaRt)
library(fgsea)


# Defining the UI page 
ui <- fluidPage(
  titlePanel("Analyze multiple processes in one all-encompassing app!"),
  
  # Main panel with analysis tabs and their respective sub-tabs
  mainPanel(
    tabsetPanel(
      id = "analyses",
      type = "tabs", 
      
      # Samples tab with its own sidebar
      tabPanel(
        "Samples", 
        sidebarLayout(
          sidebarPanel(
            h4("Sidebar for Samples"),
            #inputs Samples sidebar
            fileInput("meta_file", "upload metadata csv", accept = c(".csv",".tsv")), 
            actionButton("sum_exp","Summarize"))
          , 
          mainPanel(
            tabsetPanel(
              tabPanel("Summary", tableOutput("meta_sum")), 
              tabPanel("Table", dataTableOutput("meta_dt")), 
              tabPanel("Plots", plotOutput("dt_plts"))
            )
          )
        )
      ),
      
      # Counts tab
      tabPanel(
        "Counts", 
        sidebarLayout(
          sidebarPanel(
            h4("Sidebar for Counts"),
            #input file 
            fileInput("norm_counts", "upload normalized counts csv", accept = c(".csv",".tsv",".txt")), 
            #slider to include genes with atl x% variance 
            sliderInput("Variance",
                        "Include genes with at least this % of variance",
                        min = 0,
                        max = 100,
                        value = 50),
            
            #slider to inclide genes w at least x% non-zero
            sliderInput("non_zero",
                        "Include genes with more than number of non-zero samples",
                        min = 0,
                        max = 36, 
                        value = 18),
            
            #drop down for choosing which PCA principle components to plot
            radioButtons("PC_1_choice",
                         "Pick 1st Principle component to plot",
                         choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",  "PC9", "PC10"),
                         selected = "PC1"),
            radioButtons("PC_2_choice",
                         "Pick 2nd Principle component to plot",
                         choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",  "PC9", "PC10"),
                         selected = "PC2"),
            
            
            #action button for counts
            actionButton("do_counts", "run filtering")
          ),
          mainPanel(
            #sub taba for counts tab outputs
            tabsetPanel(
              tabPanel("Summary", tableOutput("counts_summary_table")), 
              tabPanel("Scatter Plot", plotOutput("counts_scatter_plot")), 
              tabPanel("Heatmap", plotOutput("counts_heatmap")), 
              tabPanel("PCA", plotOutput("counts_pca"))
            )
          )
        )
      ),
      
      # DE tab
      tabPanel(
        "DE", 
        sidebarLayout(
          sidebarPanel(
            h4("Sidebar for DE"),
            # inputs for the DE sidebar below
            #file upload
            fileInput("DE_counts", "upload unnormalized counts csv", accept = c(".csv",".tsv",".txt")),
            #radio buttons for users to select which two groups they want to analyse DE for 
            radioButtons('DE_g1',
                         "Select 1st experimental group to compare",
                         c("0 hr explanted dissociated adult cardiomyocytes",
                           "24 hr explanted dissociated adult cardiomyocytes",
                           "48 hr explanted dissociated adult cardiomyocytes",
                           "72 hr explanted dissociated adult cardiomyocytes",
                           "0 day postnatal ventricular myocardium",
                           "4 day post-partum ventricular myocardium",
                           "7 day post-partum ventricular myocardium",
                           "Adult ventricular myocardium",
                           "1 day post-sham surgery ventricular myocardium",
                           "1 day post-resection surgery ventricular myocardium",
                           "7 days post-resection surgery ventricular myocardium",
                           "7 days post-sham surgery ventricular myocardium",
                           "0 day postnatal myocardium (isolated CMs)",
                           "4 day postnatal myocardium (isolated CMs)",
                           "7 days post-sham surgery myocardium (isolated CMs)",
                           "7 days post-resection myocardium (isolated CMs)"
                         )),
            radioButtons('DE_g2',
                         "Select 2nd experimental group to compare",
                         c("0 hr explanted dissociated adult cardiomyocytes",
                           "24 hr explanted dissociated adult cardiomyocytes",
                           "48 hr explanted dissociated adult cardiomyocytes",
                           "72 hr explanted dissociated adult cardiomyocytes",
                           "0 day postnatal ventricular myocardium",
                           "4 day post-partum ventricular myocardium",
                           "7 day post-partum ventricular myocardium",
                           "Adult ventricular myocardium",
                           "1 day post-sham surgery ventricular myocardium",
                           "1 day post-resection surgery ventricular myocardium",
                           "7 days post-resection surgery ventricular myocardium",
                           "7 days post-sham surgery ventricular myocardium",
                           "0 day postnatal myocardium (isolated CMs)",
                           "4 day postnatal myocardium (isolated CMs)",
                           "7 days post-sham surgery myocardium (isolated CMs)",
                           "7 days post-resection myocardium (isolated CMs)"
                         )),    
            "A volcano plot can be generated with ",
            strong("log2 fold-change"),
            " on the x-axis and ",
            strong("p-adjusted"),
            " on the y-axis.",
            #radio buttons to choose volcano plot x axis
            radioButtons("x_axis",
                         "Choose the column for the x-axis: ",
                         c("baseMean" = "baseMean",
                           "log2FoldChange"= "log2FoldChange",
                           "lfcSC" = "lfcSC",
                           "stat" = "stat",
                           "pvalue" = "pvalue",
                           "padj" = "padj")
                         
            ),
            # radio buttons to choose volcano plot y axis
            radioButtons("y_axis",
                         "Choose the column for the y-axis",
                         c("baseMean" = "baseMean",
                           "log2FoldChange"= "log2FoldChange",
                           "lfcSC" = 'lfcSC',
                           "stat" = "stat",
                           "pvalue" = "pvalue",
                           "padj" = "padj")
            )
            ,
           
            #Determining colors
            colourInput("base", "Base point color: ", value="skyblue"),
            colourInput("highlight", "Highlight point color: ", value="orange"),
            
            #slider input to determinep -adjusted threshold 
            sliderInput("mag",
                        "Select the magnitude of the p adjusted threshold",
                        min = -30,
                        max = 0,
                        value = -15),
            #action button to run analysis
            actionButton("do_DE", "run differential expression")
          ),
          mainPanel(
            #sub DE tabs for output
            tabsetPanel(
              tabPanel("DE results", dataTableOutput("Deseq_res")), 
              tabPanel("Volcano Plot",plotOutput("Volcano_plt"))
            )
          )
        )
      ),
      
      # Gene enrichment analysis tab
      tabPanel(
        "Gene Enrichment Analysis", 
        sidebarLayout(
          sidebarPanel(
            h4("Please run DE analysis prior to runninig Gene Enrichment Analysis"),
            #input file 
            fileInput("mouse_gmt", "upload mice geneset file", accept = c(".gmt")), 
            #slider for padj val cutoff
            sliderInput("padj_lim_ge",
                        "Adjusted P-value threshold",
                        min = -30,
                        max = 0,
                        value = -5),
            
            #drop down for NES pathway selection 
            radioButtons("nes_type",
                         "Pick types of NES pathway to download",
                         choices = c("Positive", "Negative", "All"),
                         selected = "All"),
            
            
            #action
            actionButton("do_fgsea", "run fgsea results"),
            downloadButton("download_data", "Download results")
          ),
          mainPanel(
            tabsetPanel(
              #output for FGSEA tab analysis
              tabPanel("Barplot", plotOutput("barplt_top_pathways")),
              tabPanel("Table", dataTableOutput("tbl_nes_pathways"), DTOutput("fgsea_results_file")),
              tabPanel("scatter plot", plotOutput("scatter_plt_nes"))
              
              
            )
          )
        )
      )
    )
    
  )
)


server <- function(input, output, session){


  #Samples tab logic

  #load in metadata file
  load_metadata <- reactive({
    req(input$sum_exp) #req action button
    req(input$meta_file) # req file upload
    meta_df <- read_delim(input$meta_file$datapath)
    print("Metadata File:")
    print(head(meta_df))
    return(meta_df)
  })
  
  #return meta_df datatable 
  
  output$meta_dt <- renderDataTable({
    return(load_metadata())
  })
  
  meta_sum_table <- function(meta_data_df) {
    # Extract column names
    columns <- colnames(meta_data_df)
    
    # Ensure `types` is a character vector, will run into errors otherwise
    types <- sapply(meta_data_df, function(col) paste(class(col), collapse = ", "))
    
    # Compute summaries
    #for numeric columns compute mean +/- SD, for other complile list of unique values if possible 
    summaries <- sapply(meta_data_df, function(col) {
      if (is.numeric(col)) {
        paste0(round(mean(col, na.rm = TRUE), 1), " +/- ", round(sd(col, na.rm = TRUE), 1))
      } else if (is.factor(col) || is.character(col)) {
        paste(head(unique(col), 5), collapse = ", ") 
      } else {
        "Not summarized"
      }
    })
    
    # Combine summary into a data frame and return
    sum_table_meta <- data.frame(
      columns = columns,
      type = types,
      summary = summaries,
      stringsAsFactors = FALSE
    )
    
    return(sum_table_meta)
  }
  
  
  #call summarize fucntion  
  call_sum_meta <- reactive({
    meta_df <- load_metadata()
    return(meta_sum_table(meta_df))
  })
  #output summarized table
  output$meta_sum <- renderTable({
    return(call_sum_meta())
  })
  
  #plots tab of summarize tab
  
  plot_hist <- function(meta_data_df) {
    # Identify numeric columns
    numeric_df <- meta_data_df %>% dplyr::select(c("AvgSpotLen", "Bases", "Bytes"))
    
    #create plot for each numeric column and return combined plot
    plt1 <- numeric_df %>% ggplot(aes(x=AvgSpotLen)) +
      geom_histogram(color="darkblue", fill = "skyblue", alpha=0.7) +
      labs (title = "Histogram of Avg Spot Length", 
            y = "Frequency",
            x = "Average Spot Length") +
      theme_minimal()
    
    plt2 <- numeric_df %>% ggplot(aes(x=Bases)) +
      geom_histogram(color="darkblue", fill = "skyblue", alpha=0.7) +
      labs (title = "Histogram of Bases", 
            y = "Frequency",
            x = "Number of Bases") +
      theme_minimal()
    
    plt3 <- numeric_df %>% ggplot(aes(x=Bytes)) +
      geom_histogram(color="darkblue", fill = "skyblue", alpha=0.7) +
      labs (title = "Histogram of Bytes", 
            y = "Frequency",
            x = "Bytes") +
      theme_minimal()
    
    return(plt1 + plt2 + plt3)
  }
  
  #call histogram fucntion  
  call_hist <- reactive({
    meta_df <- load_metadata()
    return(plot_hist(meta_df))
  })
  #output summarized table
  output$dt_plts <- renderPlot({
    return(call_hist())
  })
  
  #Counts tab logic
  
  # Reactive expression to load data after action button is clicked
  load_norm_counts <- reactive({
    req(input$do_counts)  # Wait for counts action button to be pressed
    req(input$norm_counts)  # Ensure a file is uploaded
    counts_df <- read_delim(input$norm_counts$datapath)
    #make sure all genes counts are numeric and return dataframe
    counts_df <- counts_df %>%
      mutate(across(-c(Gene, GeneID, Coordinates), as.numeric))
    print("load_norm_counts executed")
    print(head(counts_df))
    return(counts_df)
  })
  
  # Function to filter data based on variance and non-zero sample count
  filter_counts <- function(counts_df, var_lim, nonzero_lim) {
    
    
    # Add variance values
    filtered_df <- counts_df %>%
      rowwise() %>%
      mutate(variance = var(c_across(-c(Gene, GeneID, Coordinates))))
    
    # Calculate variance threshold percentile
    var_threshold <- quantile(filtered_df$variance, var_lim / 100.0)
    
    # Apply filters andn return
    filtered_df <- filtered_df %>%
      filter(variance >= var_threshold) %>%
      rowwise() %>%
      filter(sum(c_across(-c(Gene, GeneID, Coordinates, variance)) != 0) >= nonzero_lim)
    print("executed filtered_counts funtion, returning filtered_df now")
    print(head(filtered_df))
    return(filtered_df)
  }
  
  # Function to create a summary DataFrame
  create_summary <- function(filtered_df, counts_df) {
    sum_df <- data.frame(
      Number_of_samples = ncol(counts_df) - 3,  # -3 for metadata columns
      Total_genes = nrow(counts_df),
      Number_passed = nrow(filtered_df),
      Percent_passed = (nrow(filtered_df) / nrow(counts_df)) * 100,
      Number_failed = nrow(counts_df) - nrow(filtered_df),
      Percent_failed = 100 - (nrow(filtered_df) / nrow(counts_df)) * 100
    )
    print("create_summary() function executed")
    print(head(sum_df))
    return(sum_df)
    
  }
  
  # Reactive expression to filter data
  filtered_data <- reactive({
    req(input$Variance, input$non_zero)  # Ensure inputs are available
    counts_df <- load_norm_counts()
    print("verified all inputs available in filtered_data reactive exp, calling filtered_counts function")
    return(filter_counts(counts_df, input$Variance, input$non_zero))
    
  })
  
  # Reactive expression to create the summary DataFrame
  summary_table <- reactive({
    counts_df <- load_norm_counts()
    filtered_df <- filtered_data()
    print("in reactive summary_table function, calling create_summary now")
    return(create_summary(filtered_df, counts_df))
    
  })
  
  # Output: Render the counts summary table
  output$counts_summary_table <- renderTable({
    summary_table()
  })
  
  # function that makes scatter plot 
  create_filtered_scatter <- function(filtered_df, counts_df){
    #creating 2 plots: subplots
    #genes passing filter: darker color, genes failing: lighter color
    
    df_to_plt <- counts_df %>% dplyr::select(Gene) %>%
      mutate(
        median_counts = apply(counts_df[, !(names(counts_df) %in% c("Gene", "GeneID", "Coordinates"))], 1, median),
        variance = apply(counts_df[, !(names(counts_df) %in% c("Gene", "GeneID", "Coordinates"))], 1, var),
        non_zero_samples = rowSums(counts_df[, !(names(counts_df) %in% c("Gene", "GeneID", "Coordinates", "variance"))] != 0),
        pass_filter = Gene %in% filtered_df$Gene,
        pass_filter = factor(pass_filter, levels = c(TRUE, FALSE))
      )
    
    #1st plot: median count vs variance 
    p1 <- df_to_plt %>% ggplot(aes(x= log(median_counts), y= log(variance), color = pass_filter)) + 
      geom_point() + scale_color_manual(values = c("TRUE" = "blue4" , "FALSE" = "cornflowerblue" ), 
                                        labels = c("Passed filers", "Failed")) +
      labs(x="log of median gene count over samples", y="variance of Gene") +
      theme_minimal()
    
    
    #2nd plot: median count vs # of zeros
    p2 <- df_to_plt %>% ggplot(aes(x=log(median_counts), y= non_zero_samples, color = pass_filter)) + 
      geom_point() + scale_color_manual(values = c("TRUE" = "darkolivegreen" , "FALSE" = "chartreuse3" ), 
                                        labels = c("Passed filers", "Failed")) +
      labs(x="Log of Median gene count over samples", y="number of samples gene appears in") +
      theme_minimal()
    
    
    final_plt <- (p1 + p2)
    print("scatter plots created, returning to be outputted")
    return(final_plt)
  }
  
  # call function to make scatter plot 
  scatter_plt <- reactive({
    counts_df <- load_norm_counts()
    filtered_df <- filtered_data()
    print("in Reactive scatter plot funtion, calling scatter plot creator")
    return(create_filtered_scatter(filtered_df, counts_df))
  })
  
  
  #Output
  output$counts_scatter_plot <- renderPlot({
    scatter_plt()
  })
  
  #Logic for HEATMAP tab 
  
  #function to make heatmap
  create_heatmap <- function(filtered_df){
    data_heatmap <- filtered_df %>% dplyr::select(-Gene, -GeneID, -Coordinates, -variance)
    data_heatmap <- data_heatmap %>% mutate(across(everything(), as.numeric)) %>%
      replace(is.na(.), 0) 
    data_heatmap <- log1p(data_heatmap)
    data_heatmap <- as.matrix(data_heatmap)
    my_colors <- colorRampPalette(c("purple", "green")) 
    plt <- pheatmap(data_heatmap, color = my_colors(10),
                    legend_labels = "Expression Level (log1p transformed)")
    print("in function to create heatmap, returning heatmap")
    return(plt)
  }
  
  #function call to make heatmap 
  heatmap_counts <- reactive({
    filtered_df <- filtered_data()
    print("in reactive heatmap fucntion")
    return(create_heatmap(filtered_df))
  })
  
  #output heatmap
  output$counts_heatmap <- renderPlot({
    heatmap_counts()
  })
  
  
  #logic to make PCA 
                    
  #function to make PCA
  make_pca <- function(filtered_df, pc1_choice, pc2_choice) {
    #running PCA
    filtered_df <- filtered_df %>% dplyr::select(-Gene, -GeneID, -Coordinates, -variance) %>% t()
    pca_obj <- prcomp(filtered_df)
    pca_res <- pca_obj$x %>% as_tibble(rownames = "sample")
    pca_vars <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2))
    pc1_index <- as.numeric(gsub("PC", "", pc1_choice))
    pc2_index <- as.numeric(gsub("PC", "", pc2_choice))
    
    #plotting and returning PCA
    plt <- ggplot(data = pca_res, 
                  aes(x = !!sym(pc1_choice), 
                      y = !!sym(pc2_choice), 
                      color = sample)) +
      geom_point() +
      labs(
        x = paste(pc1_choice, ":", round(pca_vars[pc1_index] * 100), "% variance"),
        y = paste(pc2_choice, ":", round(pca_vars[pc2_index] * 100), "% variance")
      )
    return(plt)
  }
  
  # Reactive function to call make_pca
  pca_counts <- reactive({
    filtered_df <- filtered_data()
    pc1_choice <- input$PC_1_choice
    pc2_choice <- input$PC_2_choice
    print("in reactive pca function")
    return(make_pca(filtered_df, pc1_choice, pc2_choice))
  })
  
  # Render PCA plot
  output$counts_pca <- renderPlot({
    pca_counts()
  })
  
  
  #DE tab
  #load data after action button is pressed
  load_DE_counts <- reactive({
    req(input$do_DE)  # Wait for DE counts action button to be pressed
    req(input$DE_counts)  # Ensure a file is uploaded
    DE_df <- read_delim(input$DE_counts$datapath)
    #make sure all genes counts are numeric 
    DEs_df <- DE_df %>%
      mutate(across(-c(Gene, GeneID, Coordinates), as.numeric))  
    
    print("load_DE_counts executed")
    print(head(DE_df))
    return(DE_df)
  })
  
  #Logic to calculate deseq2 results
  calc_des_seq <- function(count_data, exp1, exp2){
    gene_names <- count_data$Gene
    rownames(count_data) <- gene_names
    count_data <- count_data %>% dplyr::select(-Gene, -Coordinates, -GeneID)
    print("count_data output")
    print(head(count_data))
    
    #creating coldata 
    coldata <- data.frame(sample = c("ex0hr_1",
                                     "ex0hr_2",
                                     "ex24hr_1",
                                     "ex24hr_2",
                                     "ex48hr_1",
                                     "ex48hr_2",
                                     "ex72hr_1",
                                     "ex72hr_2",
                                     "vP0_1",
                                     "vP0_2",
                                     "vP4_1",
                                     "vP4_2",
                                     "vP7_1",
                                     "vP7_2",
                                     "vAd_1",
                                     "vAd_2",
                                     "vD1S_1",
                                     "vD1S_2",
                                     "vD1R_1",
                                     "vD1R_2",
                                     "vD7R_1",
                                     "vD7R_2",
                                     "vD7S_1",
                                     "vD7S_2",
                                     "iP0_1",
                                     "iP0_2",
                                     "iP0_3",
                                     "iP4_1",
                                     "iP4_2",
                                     "iP4_3",
                                     "iD7S_1",
                                     "iD7S_2",
                                     "iD7S_3",
                                     "iD7R_1",
                                     "iD7R_2",
                                     "iD7R_3"
    ),
    condition = c("0 hr explanted dissociated adult cardiomyocytes",
                  "0 hr explanted dissociated adult cardiomyocytes",
                  "24 hr explanted dissociated adult cardiomyocytes",
                  "24 hr explanted dissociated adult cardiomyocytes",
                  "48 hr explanted dissociated adult cardiomyocytes",
                  "48 hr explanted dissociated adult cardiomyocytes",
                  "72 hr explanted dissociated adult cardiomyocytes",
                  "72 hr explanted dissociated adult cardiomyocytes",
                  "0 day postnatal ventricular myocardium",
                  "0 day postnatal ventricular myocardium",
                  "4 day post-partum ventricular myocardium",
                  "4 day post-partum ventricular myocardium",
                  "7 day post-partum ventricular myocardium",
                  "7 day post-partum ventricular myocardium",
                  "Adult ventricular myocardium",
                  "Adult ventricular myocardium",
                  "1 day post-sham surgery ventricular myocardium",
                  "1 day post-sham surgery ventricular myocardium",
                  "1 day post-resection surgery ventricular myocardium",
                  "1 day post-resection surgery ventricular myocardium",
                  "7 days post-resection surgery ventricular myocardium",
                  "7 days post-resection surgery ventricular myocardium",
                  "7 days post-sham surgery ventricular myocardium",
                  "7 days post-sham surgery ventricular myocardium",
                  "0 day postnatal myocardium (isolated CMs)",
                  "0 day postnatal myocardium (isolated CMs)",
                  "0 day postnatal myocardium (isolated CMs)",
                  "4 day postnatal myocardium (isolated CMs)",
                  "4 day postnatal myocardium (isolated CMs)",
                  "4 day postnatal myocardium (isolated CMs)",
                  "7 days post-sham surgery myocardium (isolated CMs)",
                  "7 days post-sham surgery myocardium (isolated CMs)",
                  "7 days post-sham surgery myocardium (isolated CMs)",
                  "7 days post-resection myocardium (isolated CMs)",
                  "7 days post-resection myocardium (isolated CMs)",
                  "7 days post-resection myocardium (isolated CMs)"
    ))
    # filtering col data and counts
    coldata <- coldata[coldata$condition==exp1 | coldata$condition==exp2,]
    coldata$condition <- as.factor(coldata$condition)
    count_data <- count_data %>% dplyr::select(all_of(coldata$sample))
    
    print("filtered count data")
    print(head(count_data))
    print("final colData")
    print(head(coldata))
    counts_mat <- (as.matrix(count_data))
    counts_mat <- round(counts_mat)
    row.names(counts_mat) <- row.names(count_data)
    
    #do deseq analysis
    dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                  colData = coldata,
                                  design = ~ condition
    )
    dds <- DESeq(dds)
    dds_res <- results(dds) %>% as.tibble() %>% mutate(Gene=gene_names)
    print(head(dds_res))
    
    return(dds_res)
  }
  
  #calling calculation function 
  call_deseq <- reactive({
    DE_counts_df <- load_DE_counts()
    print("DE counts file path: ")
    print(input$DE_counts)
    exp1 <- input$DE_g1
    exp2 <- input$DE_g2
    return(calc_des_seq(DE_counts_df, exp1, exp2))
  })
  
  #function to send result to output variable
  output$Deseq_res <- renderDataTable({
    call_deseq()
  })
  
  #volcano plot logic
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    # add proper color encoding to df
    dataf <- dataf %>%
      mutate(Color = ifelse(!is.na(.data[[y_name]]) & padj < 1*10**(slider), color2, color1))
    
    # Create the plot
    plt <- ggplot(dataf, aes(x = .data[[x_name]], y = -log10(.data[[y_name]]), color = Color)) +
      geom_point() +
      theme_minimal() +
      scale_color_manual(
        values = c(color1, color2),
        name = paste0( "padj < ", slider), 
        labels = c("Not Significant", "Significant")
      )
    return(plt)
  }
  
  #call function to make volcano plot and output plot
  output$Volcano_plt <- renderPlot({volcano_plot(call_deseq(),
                                                 input$x_axis,
                                                 input$y_axis,
                                                 input$mag,
                                                 input$base,
                                                 input$highlight)})
  
  #FGSEA logic 
  
  #run fgsea
  
  run_fgsea <- function(gmt_file_path, counts_df, dds_res){
    

    #making ranked list
      
    id2gene <- counts_df %>% dplyr::select(c("Gene", "GeneID"))
    joined_res <- dds_res %>% left_join( id2gene, by='Gene')
    joined_res <- joined_res %>% arrange(desc(log2FoldChange))
    print("joined_res:")
    print(head(joined_res))
    vec <- setNames(joined_res$log2FoldChange, joined_res$GeneID)
    vec <- vec[!is.na(vec)]
    
    #running FGSEA
    pathways <- gmtPathways(gmt_file_path)
    fgseaRes <- fgsea(pathways, vec, 0, 150)
    fgseaRes <- fgseaRes %>% as_tibble() %>% arrange(padj)
    
    print("fgseaRes: ")
    print(head( fgseaRes))
    return(fgseaRes)
  }

  
  #table logic, main call fgsea function 
  
  call_fgsea2 <- reactive({
    counts_df <- load_DE_counts()
    p_adj_lim <- input$padj_lim_ge
    dds_res <-  call_deseq()
    #label res
    dds_res <- dds_res %>% arrange(padj)
    dds_res <- as_tibble(dds_res)
    dds_res <- dds_res %>% mutate(status = case_when( padj < p_adj_lim & log2FoldChange > 0 ~ "Positive",
                                                      padj < p_adj_lim & log2FoldChange <0 ~ "Negative",
                                                      padj >= p_adj_lim ~ "NS"))
    #filter dd_res
    if (input$nes_type != "All"){
      dds_res <- dds_res %>% filter(status == input$nes_type)
    }
    print("DDS_res post filterinf for fgsea")
    
    gmt_file_path <- "/projectnb/bf591-r/students/djoshi89/final-project-s-joshid/Dhruvi_Joshi_final/mh.all.v2024.1.Mm.symbols.gmt"
    #call fgsea function
    fgsea_res <- run_fgsea(gmt_file_path, counts_df, dds_res)
    
    return(fgsea_res)
    
  })

  #output pathways  
  output$tbl_nes_pathways <- renderDataTable({
    fgsea_res  <- call_fgsea2()
    p_adj_lim <- input$padj_lim_ge
    fgsea_res <- fgsea_res %>% filter(padj <= 10**(as.integer(p_adj_lim)))
    return(fgsea_res)
  })
  
  #download logic
  output$download_data <- downloadHandler(
    filename = function() {
      paste("fgsea_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the data to a CSV file
      p_adj_lim <- input$padj_lim_ge
      fgsea_res <- call_fgsea2() %>% filter(padj <= 10**(as.integer(p_adj_lim)))
      fgsea_results_clean <- fgsea_res  %>%
        mutate(across(where(is.list), ~ sapply(., toString)))
      write.csv( fgsea_results_clean , file, row.names = FALSE)
    }
  )
    
  #Barplot logic
  # filter fgsea pathways by indicated p-ajusted value
  top_pathways <- function(fgsea_results, num_paths){
    p_adj_lim <- input$padj_lim_ge
    fgsea_res <- call_fgsea2() %>% filter(padj <= 10**(as.integer(p_adj_lim)))
    top_10 <- fgsea_results %>% arrange(desc(NES)) %>% 
      filter(row_number() <= num_paths/2) %>%
      pull(pathway)
    bottom_10 <- fgsea_results %>% arrange(NES) %>% 
      filter(row_number() <= num_paths/2) %>%
      pull(pathway)
    fgsea_results <- fgsea_results %>% filter (pathway %in% c(top_10, bottom_10)) %>% arrange(desc(as.numeric(NES)))
    plt <- fgsea_results %>% ggplot(aes(x=NES, y=pathway)) +
      geom_bar(stat="identity", aes(fill = NES > 0)) + scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
      theme_minimal() + labs(title = "fgsea results for Hallmark MSigDB genes")
    return(plt)
  }
  

  
  output$barplt_top_pathways <- renderPlot({
  
    return(top_pathways(call_fgsea2(), 10))
  })
  
  #scatter plot logic 
  #runs fgsea without filtering lby padj, joind dseq log fols change values, and makes scatter plot to return
  make_fgsea_plt <- reactive({

    fgsea_res <- call_fgsea2()
    p_adj_lim <- input$padj_lim_ge
    fgsea_res <- fgsea_res %>%
      mutate(threshold = ifelse(padj <= 10^(as.integer(p_adj_lim)), "lightgreen", "darkgreen"))
    print( "FGSEA 1 :")
    print(fgsea_res)
    plt <- fgsea_res %>%
      ggplot(aes(x = NES, y = -log10(pval), color = threshold)) +  
      geom_point() +                                              
      scale_color_identity(name = "Threshold",                    
                           breaks = c("lightgreen", "darkgreen"),
                           labels = c("Significant", "Not Significant"),
                           guide = "legend") +                    
      theme_minimal() +                                        
      labs(x = "Normalized Enrichment Score (NES)",               
           y = "-log10(P-value)", 
           title = "FGSEA Results")
    return(plt)
    
  })
  
  
  output$scatter_plt_nes <- renderPlot({
    return(make_fgsea_plt())
  })
    
}

# Run the application
shinyApp(ui = ui, server = server)
