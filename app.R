# app.R

# Load required packages
suppressPackageStartupMessages({
  library(BayesTraitR)
  library(ape)
  library(ratematrix)
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(DT)
  library(colourpicker)

})

# ---------- Helpers ----------
evo_for_file <- function(evo) if (tolower(evo) == "vr") "VR" else tolower(evo)

# Build log file paths: REGION_d1_REGMODEL-EVOMODEL-REPLICATE.txt.Log.txt
build_log_paths <- function(region, regmodel, evo_file, tags_vec) {
  file.path("data", sprintf("%s_d1_%s-%s-%s.txt.Log.txt", region, regmodel, evo_file, tags_vec))}

# Choose the input data file for predBTlog based on model and x-axis selection
input_data_for_log <- function(log_path, regmodel, x_choice) {
  if (regmodel %in% c("m1","m2","m3","m4","m5")) {
    sub("\\.txt\\.Log\\.txt$", ".txt", log_path)
  } else {
    # m6–m9 depend on x-axis choice:
    if (!is.null(x_choice) && identical(x_choice, "rest of brain")) {
      sub("\\.txt\\.Log\\.txt$", "-BODmean.txt", log_path)
    } else {
      sub("\\.txt\\.Log\\.txt$", "-ROBmean.txt", log_path)
    }
  }
}

# Always read the ORIGINAL data file as ".txt" (for alignment & x variable); strict
read_table_strict <- function(p) {
  if (!file.exists(p)) stop(sprintf("Missing required data file: %s", p))
  utils::read.table(p, header = TRUE, sep = "", quote = "\"'", comment.char = "", check.names = FALSE)
}

# Decide x column name
pick_x_col <- function(df, regmodel, x_choice) {
  if (regmodel %in% c("m6","m7","m8","m9")) {
    if (!is.null(x_choice) && x_choice == "rest of brain") "rest_of_brain" else "BOD"
  } else {
    if ("rest_of_brain" %in% names(df)) "rest_of_brain"
    else if ("BOD" %in% names(df)) "BOD"
    else names(df)[3]
  }
}

# ---------- Constants / mappings ----------
REGION_LEVELS <- c("MOB","SEP","PAL","STR","AMY","HIP","SCH","NEO","MED","MES","CEREB","DIE")


# Regression model choices: names shown to user, values used internally
REG_MODEL_CHOICES <- c(
  "m1: body"                          = "m1",
  "m2: body + group + group:body"                  = "m2",
  "m3: body + body²"                  = "m3",
  "m4: ROB"                           = "m4",
  "m5: ROB + ROB²"                    = "m5",
  "m6: ROB + body"                    = "m6",
  "m7: ROB + body + body²"            = "m7",
  "m8: ROB + ROB² + body"             = "m8",
  "m9: ROB + ROB² + body + body²"     = "m9"
)

EVO_LEVELS <- c("glob","fab","bm","vr")

# Consistent colours per replicate (001..006) across the app
REPLICATE_LEVELS <- sprintf("%03d", 1:6)
REPLICATE_COLORS <- c(
  "001" = "#1b9e77",
  "002" = "#d95f02",
  "003" = "#7570b3",
  "004" = "#e7298a",
  "005" = "#66a61e",
  "006" = "#e6ab02"
)

# Group colours for m2 plots
GROUP_COLORS <- c(
  "Eulipotyphla"     = "#1b9e77",  # green
  "Afrotheria"       = "#d95f02",  # orange
  "Marsupials"       = "#7f7f7f",  # grey
  "Chiroptera_mega"  = "#1f78b4",  # dark blue
  "Chiroptera_micro" = "#a6cee3",  # light blue
  "Euarchontoglires"  = "#e7298a"   # pink
)

# ---------- UI ----------
ui <- navbarPage(
  title = "Brain Region Results",
  id = "main_nav",

  # ---- Tab 1: Marginal Likelihoods ----
  tabPanel(
    title = "Marginal Likelihoods",
    sidebarLayout(
      sidebarPanel(
        
        selectInput(
          "ml_region_filter",
          "Region filter",
          choices = c("all regions", REGION_LEVELS),
          selected = "all regions"
        ),
        
        radioButtons(
          "ml_plot_type",
          label = "Plot Type",
          choices = c("evolutionary model", "regression model"),
          selected = "evolutionary model"
        ),
        helpText("Data source: data/ResultsTable.RDS (if present)")
      ),
      mainPanel(
        plotOutput("ml_plot", height = "500px")

      )
    )
  ),

  # ---- Tab 2: ratematrix ----
  tabPanel(
    title = "ratematrix",
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          "rm_results",
          "Results",
          choices = c("ungrouped", "grouped"),
          selected = "ungrouped"
        ),
        tags$small("Loads: data/ratematrix.RDS or data/ratematrix-grouped.RDS.")
      ),
      mainPanel(
        fluidRow(
          column(
            12,
            div(
              style = "display:flex;align-items:center;gap:12px;margin:6px 0 10px;",
              actionButton("rm3_prev", "◀ Prev"),
              htmlOutput("rm_corr_label"),
              actionButton("rm3_next", "Next ▶")
            )
          )
        ),

        fluidRow(
          column(12, plotOutput("rm_plot3", height = "420px"))
        ),
        fluidRow(
          column(12, style = "height:500px:", plotOutput("rm_plot2", height = "100%"))
        ),


        fluidRow(
          column(
            12,
            div(
              style = "height:1500px;",
              plotOutput("rm_plot1", height = "100%")
            )
          )

        ),

        verbatimTextOutput("rm_status")

      )
    )
  ),

  # ---- Tab 3: regions ----
  tabPanel(
    title = "regions",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "reg_region",
          "Brain region",
          choices = REGION_LEVELS,
          selected = REGION_LEVELS[1]
        ),
        selectInput(
          "reg_evo_model",
          "Evolutionary model",
          choices = c(
            "fabric model"                         = "fab",
            "fabric model with global trend"       = "glob",
            "Brownian motion"                      = "bm",
            "variable rates (VR)"                  = "vr"
          ),
          selected = "glob"   # keep the internal value here
        ),
        selectInput(
          "reg_regr_model",
          "Regression model",
          choices  = REG_MODEL_CHOICES,
          selected = "m1"
        ),
        
        hr(),
        h4("Predicted relationship"),
        selectInput(
          "pred_dataset",
          "Dataset",
          choices = c(paste0("replicate ", 1:6), "all replicates"),
          selected = "replicate 1"
        ),
        # Summary option tied to selected dataset; disabled for "all replicates"
        conditionalPanel(
          "input.pred_dataset != 'all replicates'",
          checkboxInput("show_summary", "Show regression parameter summary table", value = FALSE)
        ),
        conditionalPanel(
          "input.pred_dataset == 'all replicates'",
          tags$em("Parameter summaries are available only when a single replicate is selected.")
        ),
        hr(),
        conditionalPanel(
          "input.reg_evo_model != 'bm'",
          checkboxInput("show_rate_summ", "View/download rate summaries", value = FALSE),
          conditionalPanel(
            "input.show_rate_summ",
            downloadButton("download_rates_file", "Download rate summaries"),
            uiOutput("rate_reminder")  # note under the download link
          ),
          conditionalPanel(
            "input.reg_evo_model == 'vr' && input.pred_dataset != 'all replicates'",
            tags$hr(),
            checkboxInput("vr_color_branches", "Colour tree by VR rates", TRUE),
            fluidRow(
              column(6, colourpicker::colourInput("vr_col_start", "Start colour", value = "#0000ff")),
              column(6, colourpicker::colourInput("vr_col_end",   "End colour",   value = "#FF0000"))
            ),
            checkboxInput("vr_use_mid", "Use intermediate colour", FALSE),
            conditionalPanel(
              "input.vr_use_mid",
              colourpicker::colourInput("vr_col_mid", "Intermediate colour", value = "#FFFFFF")
            )
          )
        )
      ),
      mainPanel(
        # 1) Summary table (optional) — tied to pred_dataset
        conditionalPanel(
          "input.show_summary",
          h4("Regression parameter summary"),
          uiOutput("summary_status"),
          DTOutput("summary_table"),
          checkboxInput("enable_summary_download", "Download summary table", value = FALSE),
          conditionalPanel(
            "input.enable_summary_download",
            downloadButton("download_summary_csv", "Download summary as CSV")
          )
        ),


        # 2) Predicted relationship plot
        h4("Average predicted relationship"),
        uiOutput("reg_formula"),

        # ----- Controls moved here (above the plot) -----
        div(
          style = "display:flex;align-items:center;gap:18px;margin:6px 0 6px;",
          checkboxInput("show_post_lines", "Add posterior sample lines", value = TRUE, width = "auto"),
          conditionalPanel(
            "input.show_post_lines",
            sliderInput("n_post_lines", "Number of posterior lines", min = 10, max = 300, value = 50, step = 10, width = "320px")
          )
        ),

        # The plot (now a bit taller)
        plotOutput("pred_plot", height = "460px"),

        # X-axis choice for models m6–m9 (shown below the plot)
        conditionalPanel(
          "input.reg_regr_model == 'm6' || input.reg_regr_model == 'm7' || input.reg_regr_model == 'm8' || input.reg_regr_model == 'm9'",
          radioButtons(
            "pred_xvar", "X-axis",
            choices = c("body", "rest of brain"),
            selected = "body",
            inline = TRUE
          )
        ),

        # Download controls (below x-axis)
        checkboxInput("enable_pred_download", "Download plotted data", value = FALSE),
        conditionalPanel(
          "input.enable_pred_download",
          downloadButton("download_pred_csv", "Download predictions as CSV")
        ),



        # 3) Tree
        conditionalPanel(
          "input.reg_evo_model != 'bm'",
          h4("Phylogenetic tree"),

          plotOutput("tree_plot", height = "1000px"),
          downloadButton("download_tree_file", "Download tree (.nex)")
        ),

        # Legend for VR colouring (shown only when VR + single replicate + colouring enabled)
        conditionalPanel(
          "input.reg_evo_model == 'vr' && input.pred_dataset != 'all replicates' && input.vr_color_branches",
          plotOutput("tree_legend", height = "160px")
        ),


        # --- Second tree: magnitude/frequency highlighting (VR + single replicate) ---
        conditionalPanel(
          "input.reg_evo_model == 'vr' && input.pred_dataset != 'all replicates'",
          tags$hr(),
          h4("VR tree (magnitude / frequency highlighting)"),
          sliderInput("alt_mag",  "magnitude", min = 1, max = 10, value = 2, step = 1, width = "300px"),
          sliderInput("alt_freq", "frequency", min = 0, max = 100, value = 95, step = 5, width = "300px"),
          plotOutput("tree_plot_alt", height = "1000px")
        )

      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {


  # ---- Tab 1: Marginal Likelihoods ----
  res_df <- reactive({
    p <- "data/ResultsTable.RDS"
    if (file.exists(p)) readRDS(p) else NULL
  })

  output$ml_plot <- renderPlot({
    df <- res_df()
    validate(need(!is.null(df), "No dataset available."))
    
    # --- Region filtering (case-insensitive column name) ---
    reg_col <- names(df)[tolower(names(df)) == "region"][1]
    if (!is.na(reg_col) && !identical(input$ml_region_filter, "all regions")) {
      df <- df[df[[reg_col]] == input$ml_region_filter, , drop = FALSE]
    }
    validate(need(nrow(df) > 0, "No rows for the selected region."))
    
    # --- Existing plot logic (unchanged) ---
    if (input$ml_plot_type == "evolutionary model") {
      x_var <- "form"; group_var <- "form"; fill_var <- "form"; facet_var <- "model"; color_var <- "form"
    } else {
      x_var <- "model"; group_var <- "model"; fill_var <- "form"; facet_var <- "form"; color_var <- "model"
    }
    
    ggplot(
      df,
      aes_string(x = x_var, y = "mL", group = group_var, fill = fill_var, color = color_var)
    ) +
      geom_point(size = 3, position = position_jitter(width = 0.07, height = 0)) +
      facet_wrap(as.formula(paste0("~", facet_var)), scales = "free_y") +
      theme_bw() +
      labs(x = x_var, color = color_var, fill = fill_var, y = "mL")
  }, height = 1000)
  

  # ---- Tab 2: ratematrix ----
  ratematrix_obj <- reactive({
    path <- if (input$rm_results == "grouped") "data/ratematrix-grouped.RDS" else "data/ratematrix.RDS"
    if (file.exists(path)) readRDS(path) else NULL
  })

  # Some reactives required to switch between plots

  # correlation ID value (for ratematrix)
  corr_idx <- reactiveVal(1)  # 1..12

  observeEvent(input$rm3_prev, {
    i <- corr_idx()
    corr_idx(ifelse(i <= 1, 12, i - 1))
  })

  observeEvent(input$rm3_next, {
    i <- corr_idx()
    corr_idx(ifelse(i >= 12, 1, i + 1))
  })

  output$rm_plot3_dbg <- renderText({
    paste("Active region:", REGION_LEVELS[corr_idx()])
  })

  # Add a brain region anchor for the plot switch
  output$rm_corr_label <- renderUI({
    anchor <- REGION_LEVELS[corr_idx()]
    HTML(sprintf("<b>Region group:</b> %s-*", anchor))
  })

  output$rm_plot1 <- renderPlot({
    # Load: grouped/ungrouped RM output
    path <- file.path("data",
                      if (input$rm_results == "grouped") "grouped-RM_output.RDS" else "ungrouped-RM_output.RDS")
    validate(need(file.exists(path), paste0("File not found: ", path)))

    obj <- readRDS(path)
    RM_out_mcc <- if (is.list(obj) && !is.null(obj$RM_out_mcc)) obj$RM_out_mcc else obj
    validate(need(!is.null(RM_out_mcc), "RM_out_mcc object not found in the RDS."))

    par(mar = c(4, 4, 2, 2))
    ratematrix::plotRatematrix(RM_out_mcc, show.zero = TRUE, colors = "pink")
  }, height = 1500, res = 96)


  output$rm_plot2 <- renderPlot({
    # Load data file
    path <- file.path("data",
                      if (input$rm_results == "grouped") "grouped-var_mcc.RDS" else "ungrouped-var_mcc.RDS")
    validate(need(file.exists(path), paste0("File not found: ", path)))


    obj <- readRDS(path)
    # The RDS may contain a list with $sigmas_mcc or be the data.frame itself
    sigmas_mcc <- if (is.list(obj) && !is.null(obj$sigmas_mcc)) obj$sigmas_mcc else obj

    validate(need(is.data.frame(sigmas_mcc), "Could not find a data.frame `sigmas_mcc` in the RDS."))
    validate(need(all(c("Rate", "Region") %in% names(sigmas_mcc)),
                  "The `sigmas_mcc` data must contain columns 'Rate' and 'Region'."))

    sigmas_mcc$Rate = log10(sigmas_mcc$Rate)

    ggplot(data = sigmas_mcc, aes(x = Rate, group = Region, fill = Region)) +
      geom_density(adjust = .5, alpha = 0.5) +
      scale_x_continuous(expand = expansion(mult = 0.02)) +
      scale_y_continuous(expand = expansion(mult = 0.02)) +
      labs(x = "Variance", y = NULL) +
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Remove panel background
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.5)
      )
  }, height = 500, res = 120)

  # Load once and normalize names/casing for corr data (used by Plot 3)
  corr_all <- reactive({
    path <- file.path("data",
                      if (input$rm_results == "grouped") "grouped-corr_mcc.RDS" else "ungrouped-corr_mcc.RDS")

    req(file.exists(path))
    obj <- readRDS(path)
    df <- if (is.list(obj) && !is.null(obj$corr_mcc)) obj$corr_mcc else obj
    stopifnot(is.data.frame(df))
    names(df) <- tolower(names(df))
    stopifnot(all(c("corr","region") %in% names(df)))
    df$region <- toupper(trimws(as.character(df$region)))
    df
  })

  # Global x-range across ALL 12 plots
  corr_xlim <- reactive({
    df <- corr_all()
    range(df$corr[is.finite(df$corr)], na.rm = TRUE)
  })

  output$rm_plot3 <- renderPlot({
    df_all <- corr_all()
    anchor_name <- toupper(REGION_LEVELS[corr_idx()])   # ticker picks 1..12 by name
    prefix <- sub("^\\s*([^\\-]+)\\s*-.*$", "\\1", df_all$region)
    subset_df <- df_all[prefix == anchor_name, , drop = FALSE]
    validate(need(nrow(subset_df) > 0,
                  paste0("No rows beginning with '", anchor_name, "-'. Check 'region' strings.")))

    ggplot(data = subset_df, aes(x = corr, group = region, fill = region)) +
      geom_density(adjust = 1.5, alpha = 0.5) +
      coord_cartesian(xlim = corr_xlim()) +                               # <— fixed global limits
      scale_x_continuous(expand = expansion(mult = 0.02)) +
      scale_y_continuous(expand = expansion(mult = 0.02)) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line = element_line(color = "black", linewidth = 0.5)
      ) +
      labs(
        title = paste0("Covariances for regions starting with ", anchor_name),
        x = "Covariance", y = NULL
      )
  }, res = 120)




  # ---- Tab 3: regions ----
  # Auto-disable summary when "all replicates" is selected
  observe({
    if (identical(input$pred_dataset, "all replicates") && isTRUE(input$show_summary)) {
      updateCheckboxInput(session, "show_summary", value = FALSE)
    }
  })

  # -------- Regression formula (display under the title) --------
  output$reg_formula <- renderUI({
    # Find the label (name) corresponding to the selected value
    lab <- names(REG_MODEL_CHOICES)[match(input$reg_regr_model, REG_MODEL_CHOICES)]
    if (is.na(lab)) lab <- "(unknown model)"
    # Strip the leading "m#: " if you only want the formula part:
    pretty_lab <- sub("^m[1-9]:\\s*", "", lab)
    tags$p(tags$strong("Regression formula: "), pretty_lab)
  })
  
  

  # 2) Predicted relationship plot (strict: sort + cbind medians and posterior)

  # Minimal group map from fullData (Species_TTOL + Group)
  full_groups <- reactive({
    path <- file.path("data", "fullData.RDS")
    if (!file.exists(path)) return(NULL)
    df <- readRDS(path)
    out <- df[, c("Species_TTOL", "Group"), drop = FALSE]
    out$Group <- as.character(out$Group)
    out
  })

  # Prediction data reactive
  pred_data <- reactive({

    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    evo_file <- evo_for_file(input$reg_evo_model)

    # Which replicates?
    target_tags <- if (tolower(input$pred_dataset) == "all replicates") {
      sprintf("%03d", 1:6)
    } else {
      sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
    }

    # Build paths
    log_paths <- build_log_paths(region, regmodel, evo_file, target_tags)
    input_paths <- vapply(log_paths, input_data_for_log, FUN.VALUE = character(1),
                          regmodel = regmodel, x_choice = input$pred_xvar)
    original_path <- sub("\\.txt\\.Log\\.txt$", ".txt", log_paths[1])  # original .txt is the same for all reps

    # --- Read original dataset (STRICT) and sort by Species_TTOL ---
    df_org <- tryCatch(read_table_strict(original_path), error = function(e) NULL)
    if (is.null(df_org)) return(NULL)
    if (!"Species_TTOL" %in% names(df_org)) return(NULL)
    df_org <- df_org[order(df_org$Species_TTOL), , drop = FALSE]

    # Prepare containers for medians and (optional) posterior samples
    med_cols <- list()
    post_cols <- list()

    for (i in seq_along(log_paths)) {
      log_p <- log_paths[i]
      inp_p <- input_paths[i]
      if (!file.exists(log_p) || !file.exists(inp_p)) next

      # Predictions matrix from BayesTraitR
      pred_mat <- tryCatch({
        BayesTraitR::predBTlog(input = inp_p, output = log_p)
      }, error = function(e) NULL)
      if (is.null(pred_mat)) next
      pred_mat <- as.matrix(pred_mat)

      # --- Sort predictions by rownames alphabetically ---
      rn <- rownames(pred_mat)
      if (is.null(rn)) next
      ord_pred <- order(rn)
      pred_mat <- pred_mat[ord_pred, , drop = FALSE]
      rn <- rn[ord_pred]

      # Align by common, alphabetically-sorted species
      common_sp <- intersect(df_org$Species_TTOL, rn)
      if (length(common_sp) == 0) next
      common_sp <- sort(common_sp)

      df_use <- df_org[df_org$Species_TTOL %in% common_sp, , drop = FALSE]
      df_use <- df_use[order(df_use$Species_TTOL), , drop = FALSE]
      pred_mat <- pred_mat[match(common_sp, rn), , drop = FALSE]

      # Row medians
      med <- apply(pred_mat, 1, stats::median, na.rm = TRUE)
      med_cols[[i]] <- med
      names(med_cols)[i] <- paste0("pred_med_", target_tags[i])

      # Posterior columns (random k)
      if (isTRUE(input$show_post_lines) && ncol(pred_mat) > 0) {
        k <- max(1, min(ncol(pred_mat), input$n_post_lines))
        set.seed(20251015 + as.integer(target_tags[i]))
        sel <- sample(seq_len(ncol(pred_mat)), k)
        pc <- as.data.frame(pred_mat[, sel, drop = FALSE], check.names = FALSE)
        colnames(pc) <- paste0("post_", target_tags[i], "_", seq_along(sel))
        post_cols[[i]] <- pc
        names(post_cols)[i] <- target_tags[i]
      } else {
        post_cols[[i]] <- NULL
        names(post_cols)[i] <- target_tags[i]
      }

      # Keep original aligned to the first usable replicate
      df_org <- df_use
    }

    if (length(Filter(Negate(is.null), med_cols)) == 0) return(NULL)

    # cbind medians and posterior columns to the aligned original df
    for (nm in names(med_cols)) df_org[[nm]] <- med_cols[[nm]]
    for (i in seq_along(post_cols)) {
      pc <- post_cols[[i]]
      if (!is.null(pc)) df_org <- cbind(df_org, pc)
    }

    # Decide x column name
    x_col <- pick_x_col(df_org, regmodel, input$pred_xvar)

    # Attach Group only for m2 (body~grp)
    if (identical(input$reg_regr_model, "m2")) {
      grp <- full_groups()
      if (!is.null(grp)) {
        idx <- match(df_org$Species_TTOL, grp$Species_TTOL)
        df_org$Group <- grp$Group[idx]
        df_org$Group <- factor(df_org$Group)
      }
    }

    list(df = df_org,x_col = x_col)
  })

  # --- REPLACE load_summary_for_inset() with this exact version ---
  load_summary_for_inset <- reactive({
    # Only meaningful for a single replicate
    if (identical(input$pred_dataset, "all replicates")) return(NULL)
    
    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    evo_file <- if (tolower(input$reg_evo_model) == "vr") "VR" else tolower(input$reg_evo_model)
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
    path <- file.path("data", sprintf("%s_d1_%s-%s-%s-regSum.RDS", region, regmodel, evo_file, tag))
    if (!file.exists(path)) return(NULL)
    
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    
    # Keep exactly as stored; don't coerce. Rownames must include "Median", "pMCMC".
    if (is.data.frame(obj)) return(obj)
    if (is.matrix(obj)) return(as.data.frame(obj, stringsAsFactors = FALSE))
    as.data.frame(obj, stringsAsFactors = FALSE)
  })
  

  extract_model_inset <- function(df, model_code) {
    if (is.null(df) || !NROW(df)) return(NULL)
    
    # Helper: exact row/column lookup → numeric
    get_num <- function(row_lab, col_lab) {
      if (!(row_lab %in% rownames(df))) return(NA_real_)
      if (!(col_lab %in% colnames(df))) return(NA_real_)
      suppressWarnings(as.numeric(df[row_lab, col_lab]))
    }
    
    # Formatting: 4 dp for coefficients, 3 dp for p-values (as per your example)
    fmt4 <- function(x) if (is.na(x)) "NA" else format(round(x, 4), nsmall = 4, trim = TRUE, scientific = FALSE)
    fmt3 <- function(x) if (is.na(x)) "NA" else format(round(x, 3), nsmall = 3, trim = TRUE, scientific = FALSE)
    
    # Common pulls
    R2 <- get_num("Median", "R^2")
    B1 <- get_num("Median", "Beta 1"); P1 <- get_num("pMCMC", "Beta 1")
    B2 <- get_num("Median", "Beta 2"); P2 <- get_num("pMCMC", "Beta 2")
    B3 <- get_num("Median", "Beta 3"); P3 <- get_num("pMCMC", "Beta 3")
    B4 <- get_num("Median", "Beta 4"); P4 <- get_num("pMCMC", "Beta 4")
    
    txt <- switch(
      model_code,
      "m1" = paste0(
        "median slope = ", fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m2" = paste0(
        "R² = ", fmt4(R2), "\n",
        "(group diffs coming soon)"
      ),
      "m3" = paste0(
        "body slope = ",        fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "body quadratic = ",    fmt4(B2), " [pₓ = ", fmt3(P2), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m4" = paste0(
        "ROB slope = ",         fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m5" = paste0(
        "ROB slope = ",         fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "ROB quadratic = ",     fmt4(B2), " [pₓ = ", fmt3(P2), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m6" = paste0(
        "ROB slope = ",         fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "body slope = ",        fmt4(B2), " [pₓ = ", fmt3(P2), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m7" = paste0(
        "ROB slope = ",         fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "body slope = ",        fmt4(B2), " [pₓ = ", fmt3(P2), "]\n",
        "body quadratic = ",    fmt4(B3), " [pₓ = ", fmt3(P3), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m8" = paste0(
        "ROB slope = ",         fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "ROB quadratic = ",     fmt4(B2), " [pₓ = ", fmt3(P2), "]\n",
        "body slope = ",        fmt4(B3), " [pₓ = ", fmt3(P3), "]\n",
        "R² = ", fmt4(R2)
      ),
      "m9" = paste0(
        "ROB slope = ",         fmt4(B1), " [pₓ = ", fmt3(P1), "]\n",
        "ROB quadratic = ",     fmt4(B2), " [pₓ = ", fmt3(P2), "]\n",
        "body slope = ",        fmt4(B3), " [pₓ = ", fmt3(P3), "]\n",
        "body quadratic = ",    fmt4(B4), " [pₓ = ", fmt3(P4), "]\n",
        "R² = ", fmt4(R2)
      ),
      NULL
    )
    
    list(text = txt)
  }
  
  

  output$pred_plot <- renderPlot({
    pd <- pred_data()
    validate(need(!is.null(pd), "No dataset available."))

    df    <- pd$df
    x_col <- pd$x_col

    prettify_x <- function(xn) if (xn == "rest_of_brain") "Rest of brain" else if (xn == "BOD") "Body" else xn
    x_lab <- prettify_x(x_col)

    # ----- Special plotting for m2: colour by Group (ignore replicate colours) -----
    if (identical(input$reg_regr_model, "m2")) {
      validate(need("Group" %in% names(df), "Group labels not found (check data/fullData.RDS)."))

      # Posterior spaghetti (all replicates, but coloured by Group)
      p <- ggplot() + theme_bw()

      if (isTRUE(input$show_post_lines)) {
        post_cols <- grep("^post_[0-9]{3}_", names(df), value = TRUE)
        if (length(post_cols)) {
          post_long <- tidyr::pivot_longer(
            df[, c(x_col, "Group", post_cols), drop = FALSE],
            cols = tidyselect::all_of(post_cols),
            names_to = "sample", values_to = "y"
          )
          p <- p + geom_line(
            data = post_long,
            mapping = aes(x = .data[[x_col]], y = y,
                          color = Group, group = interaction(Group, sample)),
            alpha = 0.18, linewidth = 0.3, show.legend = TRUE
          )
        }
      }

      # Median lines (one per replicate per group if multiple reps are present)
      med_cols <- grep("^pred_med_[0-9]{3}$", names(df), value = TRUE)
      if (length(med_cols)) {
        med_long <- tidyr::pivot_longer(
          df[, c(x_col, "Group", med_cols), drop = FALSE],
          cols = tidyselect::all_of(med_cols),
          names_to = "rep", values_to = "y"
        )
        p <- p + geom_line(
          data = med_long,
          mapping = aes(x = .data[[x_col]], y = y, color = Group,
                        group = interaction(Group, rep)),
          linewidth = 1.2, show.legend = TRUE
        )
      }

      # Use pretty colours per group (only for groups present)
      groups_present <- intersect(names(GROUP_COLORS), unique(as.character(df$Group)))
      p <- p + scale_color_manual(
        name   = "Group",
        values = GROUP_COLORS[groups_present],
        breaks = groups_present,
        drop   = FALSE
      )

      # ---- Inset parameter text (top-right) for m2 ----
      txt <- NULL
      if (!identical(input$pred_dataset, "all replicates")) {
        inset_df <- load_summary_for_inset()
        ins <- extract_model_inset(inset_df, "m2")
        if (!is.null(ins) && !is.null(ins$text)) {
          # ins$text already contains: "R² = ...\n(group diffs coming soon)"
          txt <- ins$text
        }
      } else {
        # All replicates selected: we can’t pull a single R²
        txt <- "R² (select a single replicate)\n(group diffs coming soon)"
      }

      if (!is.null(txt)) {
        p <- p + annotate(
          "text",
          x = Inf, y = Inf,
          label = txt,
          hjust = 1.02, vjust = 1.2,
          size = 3.6,
          lineheight = 1.05
        )
      }

      return(p + labs(x = x_lab, y = "Predicted (median across posterior)"))

    }

    # ----- Default plotting (non-m2): your original replicate-colour version -----
    p <- ggplot() + theme_bw()

    # Posterior spaghetti (if present)
    if (isTRUE(input$show_post_lines)) {
      post_all <- grep("^post_[0-9]{3}_", names(df), value = TRUE)
      if (length(post_all)) {
        reps_in_df <- sort(unique(sub("^post_([0-9]{3})_.*$", "\\1", post_all)))
        for (rep_id in reps_in_df) {
          rep_cols <- grep(paste0("^post_", rep_id, "_"), names(df), value = TRUE)
          if (!length(rep_cols)) next
          post_long <- tidyr::pivot_longer(
            df[, c(x_col, rep_cols), drop = FALSE],
            cols = tidyselect::all_of(rep_cols),
            names_to = "sample", values_to = "y"
          )
          post_long$rep <- rep_id
          p <- p + geom_line(
            data = post_long,
            mapping = aes(x = .data[[x_col]], y = y,
                          group = interaction(rep, sample), color = rep),
            alpha = 0.20, linewidth = 0.3, show.legend = FALSE
          )
        }
      }
    }

    # Median lines for ALL available replicates (together)
    med_cols <- grep("^pred_med_[0-9]{3}$", names(df), value = TRUE)
    if (length(med_cols)) {
      med_long <- tidyr::pivot_longer(
        df[, c(x_col, med_cols), drop = FALSE],
        cols = tidyselect::all_of(med_cols),
        names_to = "rep", values_to = "y"
      )
      med_long$rep <- sub("^pred_med_", "", med_long$rep)
      p <- p + geom_line(
        data = med_long,
        mapping = aes(x = .data[[x_col]], y = y, color = rep, group = rep),
        linewidth = 1.2
      )
    }

    # ---- Inset parameter text (top-right) for m1/m3–m9 ----
    txt <- NULL
    if (!identical(input$pred_dataset, "all replicates")) {
      inset_df <- load_summary_for_inset()
      ins <- extract_model_inset(inset_df, input$reg_regr_model)
      if (!is.null(ins) && !is.null(ins$text)) txt <- ins$text
    }

    if (!is.null(txt)) {
      p <- p + annotate(
        "text",
        x = Inf, y = Inf,
        label = txt,
        hjust = 1.02, vjust = 1.2,
        size = 3.6,
        lineheight = 1.05
      )
    }

    p + scale_color_manual(
      name   = "Replicate",
      breaks = REPLICATE_LEVELS,
      values = REPLICATE_COLORS
    ) +
      labs(x = x_lab, y = "Predicted (median across posterior)")
  })


  # Download the plotted predictions (original + medians + posterior samples)
  output$download_pred_csv <- downloadHandler(
    filename = function() {
      region   <- input$reg_region
      evo_file <- evo_for_file(input$reg_evo_model)
      regmodel <- input$reg_regr_model
      dataset  <- gsub("\\s+", "", tolower(input$pred_dataset))
      paste0("predictions_", region, "_", regmodel, "_", evo_file, "_", dataset, ".csv")
    },
    content = function(file) {
      pd <- pred_data()
      if (is.null(pd) || is.null(pd$df) || !NROW(pd$df)) {
        stop("No plotted data available to download for the current selection.")
      }
      df_out <- pd$df
      if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_csv(df_out, file)
      } else {
        utils::write.csv(df_out, file, row.names = FALSE)
      }
    }
  )

  # -------- Regression parameter summary (single replicate only) --------
  .evo_to_filename <- function(evo) {
    if (tolower(evo) == "vr") "VR" else tolower(evo)
  }

  summary_paths <- reactive({
    req(input$pred_dataset != "all replicates")
    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    evo_file <- .evo_to_filename(input$reg_evo_model)
    tag <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
    file.path("data", sprintf("%s_d1_%s-%s-%s-regSum.RDS", region, regmodel, evo_file, tag))
  })

  summary_data <- reactive({
    req(isTRUE(input$show_summary))
    path <- summary_paths()
    if (!file.exists(path)) return(NULL)
    
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    
    # Coerce to data.frame while preserving row names if present
    df <- if (is.data.frame(obj)) obj else as.data.frame(obj, stringsAsFactors = FALSE)
    
    if (nrow(df) == 0) return(df)
    
    # --- Keep top row verbatim; process the rest numerically ---
    top_row <- df[1, , drop = FALSE]
    rest    <- if (nrow(df) >= 2) df[-1, , drop = FALSE] else df[0, , drop = FALSE]
    
    if (nrow(rest)) {
      # Best-effort numeric typing on the 'rest' only
      if (requireNamespace("readr", quietly = TRUE)) {
        suppressWarnings({
          rest <- readr::type_convert(rest, guess_integer = TRUE, trim_ws = TRUE,
                                      locale = readr::locale(decimal_mark = "."))
        })
      } else {
        # Fallback: coerce character columns that are mostly numeric
        is_mostly_numeric <- function(x) {
          if (!is.character(x)) return(FALSE)
          sup <- suppressWarnings(as.numeric(x))
          sum(!is.na(sup)) >= (length(x) - 1L)
        }
        char_cols <- vapply(rest, is.character, logical(1))
        mostly_num_cols <- names(rest)[char_cols & vapply(rest, is_mostly_numeric, logical(1))]
        for (nm in mostly_num_cols) {
          rest[[nm]] <- suppressWarnings(as.numeric(rest[[nm]]))
        }
      }
      
      # Round numeric columns (rows 2..n only)
      num_cols <- vapply(rest, is.numeric, logical(1))
      rest[num_cols] <- lapply(rest[num_cols], function(x) round(x, 4))
      
      # For display in DT (and to keep top row intact), convert the 'rest' to strings
      # with consistent formatting for numerics.
      fmt_num <- function(x) ifelse(is.na(x), NA_character_,
                                    formatC(x, digits = 4, format = "fg", drop0trailing = TRUE))
      for (nm in names(rest)) {
        if (is.numeric(rest[[nm]])) {
          rest[[nm]] <- fmt_num(rest[[nm]])
        } else {
          # keep as character as-is
          rest[[nm]] <- as.character(rest[[nm]])
        }
      }
    }
    
    # Bind top row back on (everything is character now, preserving the text in row 1)
    out <- rbind(top_row, rest)
    rownames(out) <- rownames(df)  # preserve any original row names if present
    out
  })
  

  output$summary_status <- renderUI({
    req(isTRUE(input$show_summary))
    dat  <- summary_data()
    path <- try(summary_paths(), silent = TRUE)
    if (!is.null(dat)) return(NULL)
    tags$div(
      class = "alert alert-warning",
      tags$strong("No summary file found for the current selection."),
      tags$p("Tried path:"),
      tags$code(if (inherits(path, "try-error")) "(unknown path)" else path)
    )
  })
  
  output$summary_table <- renderDT({
    req(isTRUE(input$show_summary))
    dat <- summary_data()
    req(!is.null(dat))
    datatable(dat, rownames = TRUE, options = list(pageLength = 15, scrollX = TRUE))
  })
  
  output$download_summary_csv <- downloadHandler(
    filename = function() {
      region   <- input$reg_region
      evo_file <- .evo_to_filename(input$reg_evo_model)
      regmodel <- input$reg_regr_model
      tag <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
      paste0("regSummary_", region, "_", regmodel, "_", evo_file, "_", tag, ".csv")
    },
    content = function(file) {
      req(isTRUE(input$show_summary))
      dat <- summary_data()
      req(!is.null(dat))
      if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_csv(dat, file)
      } else {
        utils::write.csv(dat, file, row.names = TRUE)
      }
    }
  )

  # -------- Rate summaries: simple download (no sliders), evo-specific paths --------
  output$rate_reminder <- renderUI({
    req(input$reg_evo_model != "bm", isTRUE(input$show_rate_summ))
    evo <- tolower(input$reg_evo_model)
    if (evo %in% c("glob","fab")) {
      tags$div(style = "margin-top:6px;",
               tags$em("Merged Fabric results are identical for all six replicates."))
    } else if (evo == "vr" && identical(input$pred_dataset, "all replicates")) {
      tags$div(class = "alert alert-warning",
               "Select a single replicate to download VR rate summary.")
    } else {
      NULL
    }
  })

  output$download_rates_file <- downloadHandler(
    filename = function() {
      evo <- tolower(input$reg_evo_model)
      region   <- input$reg_region
      regmodel <- input$reg_regr_model
      if (evo %in% c("glob","fab")) {
        # Use the same name as the source CSV
        paste0(region, "_d1_", regmodel, "-", evo, "-MergedFabricResults.csv")
      } else if (evo == "vr") {
        validate(need(input$pred_dataset != "all replicates",
                      "Select a single replicate to download VR rate summary."))
        tag <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
        # We output CSV (converted from RDS)
        paste0(region, "_d1_", regmodel, "-VR-", tag, "-rateSummary.csv")
      } else {
        "rate_summary.csv"
      }
    },
    content = function(file) {
      evo <- tolower(input$reg_evo_model)
      region   <- input$reg_region
      regmodel <- input$reg_regr_model

      if (evo %in% c("glob","fab")) {
        src <- file.path("data", sprintf("%s_d1_%s-%s-MergedFabricResults.csv", region, regmodel, evo))
        validate(need(file.exists(src), paste0("Rate summary file not found: ", src)))
        file.copy(src, file, overwrite = TRUE)

      } else if (evo == "vr") {
        validate(need(input$pred_dataset != "all replicates",
                      "Select a single replicate to download VR rate summary."))
        tag <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
        src <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
        validate(need(file.exists(src), paste0("Rate summary file not found: ", src)))
        obj <- readRDS(src)
        df  <- as.data.frame(obj$ratesummary)
        utils::write.csv(df, file, row.names = FALSE)

      } else {
        stop("Unsupported evolutionary model.")
      }
    }
  )

  # -------- Tree plot + download (single replicate for glob/fab/vr) --------

  # Read the VR rateSummary RDS once and prepare:
  # - tree with branch lengths scaled by medianscalar[-1]
  # - rate table without the first row (rt2)
  # - posterior (for frequency thresholding in alt plot)
  vr_data <- reactive({
    req(tolower(input$reg_evo_model) == "vr", input$pred_dataset != "all replicates")
    
    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
    
    rds_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
    validate(need(file.exists(rds_path), paste0("VR rate summary not found: ", rds_path)))
    
    obj <- readRDS(rds_path)
    
    # Extract components
    validate(need(!is.null(obj$tree), "VR RDS does not contain a $tree object."))
    validate(need(!is.null(obj$ratesummary), "VR RDS does not contain a $ratesummary table."))
    
    tr <- obj$tree
    rt <- as.data.frame(obj$ratesummary, stringsAsFactors = FALSE)
    
    validate(need("medianscalar" %in% names(rt), "ratesummary must contain 'medianscalar'."))
    
    # Prepare the -1 version of the rates (rows map to edges)
    rt2 <- rt[-1, , drop = FALSE]
    
    # Scale branch lengths by medianscalar[-1]
    nE <- nrow(tr$edge)
    ms <- rt$medianscalar[-1]
    
    if (length(ms) < nE) {
      stop(sprintf("Not enough medianscalar values (%d) to match edges (%d).", length(ms), nE))
    }
    if (length(ms) > nE) {
      ms <- ms[seq_len(nE)]  # truncate if an extra row sneaks in
    }
    validate(need(all(is.finite(ms) & ms > 0), "All medianscalar[-1] must be positive and finite."))
    
    tr$edge.length <- tr$edge.length * ms
    
    list(
      tree      = tr,
      rates2    = rt2,
      posterior = obj$posterior  # used in the alt (magnitude/frequency) plot
    )
  })
  

  # --- REPLACE the old fab_rates() reactive with this ---
  fab_merged <- reactive({
    req(tolower(input$reg_evo_model) %in% c("fab","glob"))
    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    evo      <- tolower(input$reg_evo_model)
    
    # One merged file across replicates, e.g. MED_d1_m7-glob-MergedFabricResults.csv
    csv_path <- file.path("data", sprintf("%s_d1_%s-%s-MergedFabricResults.csv", region, regmodel, evo))
    validate(need(file.exists(csv_path), paste0("Merged Fabric results not found: ", csv_path)))
    
    # Keep default check.names=TRUE so we get "Mean..Beta...BL..NZ" etc.
    df <- utils::read.csv(csv_path, header = TRUE)
    validate(need(nrow(df) >= 2, "Merged Fabric file has insufficient rows."))
    
    # Drop first row as specified
    df <- df[-1, , drop = FALSE]
    df
  })
  
  vr_edge_colors <- reactive({
    req(isTRUE(input$vr_color_branches))
    rt <- vr_rate_tbl()
    validate(need("medianscalar" %in% names(rt), "rates table must contain 'medianscalar'."))

    # 1) log10 range of medianscalar (guard against non-positive)
    ms_rounded <- round(rt$medianscalar, 2)
    validate(need(all(ms_rounded > 0, na.rm = TRUE),
                  "All 'medianscalar' values must be > 0 for log10 scaling."))
    log_vals <- log10(ms_rounded)
    rng <- range(log_vals, finite = TRUE, na.rm = TRUE)

    # 2) sequence by 0.01
    grid <- seq(rng[1], rng[2], by = 0.01)
    n <- length(grid)
    validate(need(n >= 2, "Insufficient range to build colour grid."))

    # 3) palette from user colours
    cols_in <- if (isTRUE(input$vr_use_mid)) {
      c(input$vr_col_start, input$vr_col_mid, input$vr_col_end)
    } else {
      c(input$vr_col_start, input$vr_col_end)
    }
    pal <- grDevices::colorRampPalette(cols_in)
    lut <- pal(n)  # look-up table for grid points

    # 4) map each (rounded) medianscalar to nearest grid colour (on log10 scale)
    idx <- vapply(log_vals, function(v) {
      which.min(abs(grid - v))
    }, integer(1L))
    edge_cols_full <- lut[pmax(1, pmin(n, idx))]

    # 5/6) drop the first value per spec (align with edge ordering expectation)
    # (We leave length checks to the plotting call; this matches your stated requirement.)
    edge_cols_full[-1]
  })

  vr_palette_params <- reactive({
    req(tolower(input$reg_evo_model) == "vr",
        input$pred_dataset != "all replicates",
        isTRUE(input$vr_color_branches))
    
    vd <- vr_data()
    rt2 <- vd$rates2
    validate(need("medianscalar" %in% names(rt2), "rates table must contain 'medianscalar'."))
    
    ms <- rt2$medianscalar
    validate(need(all(ms > 0, na.rm = TRUE), "All 'medianscalar' must be > 0 for log10 scaling."))
    ms_rounded <- round(ms, 2)
    
    log_vals <- log10(ms_rounded)
    rng  <- range(log_vals, finite = TRUE, na.rm = TRUE)
    grid <- seq(rng[1], rng[2], by = 0.01)
    
    cols_in <- if (isTRUE(input$vr_use_mid)) {
      c(input$vr_col_start, input$vr_col_mid, input$vr_col_end)
    } else {
      c(input$vr_col_start, input$vr_col_end)
    }
    lut <- grDevices::colorRampPalette(cols_in)(length(grid))
    
    list(
      grid   = grid,
      lut    = lut,
      min_s  = min(ms_rounded, na.rm = TRUE),
      med_s  = stats::median(ms_rounded, na.rm = TRUE),
      max_s  = max(ms_rounded, na.rm = TRUE)
    )
  })
  output$tree_plot <- renderPlot({
    tryCatch({
      req(input$reg_evo_model != "bm")
      validate(need(input$pred_dataset != "all replicates",
                    "Tree plotting is only available for a single replicate."))
      
      evo      <- tolower(input$reg_evo_model)
      region   <- input$reg_region
      regmodel <- input$reg_regr_model
      tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
      
      # Tree path (unchanged)
      tree_path <- if (evo %in% c("glob","fab")) {
        file.path("data", "regions_tree.trees")
      } else if (evo == "vr") {
        file.path("data", sprintf("%s_d1_%s-VR-%s-medTree.trees", region, regmodel, tag))
      } else NA_character_
      
      validate(need(!is.na(tree_path) && file.exists(tree_path),
                    paste0("Tree file not found: ", tree_path)))
      tr <- ape::read.nexus(tree_path)
      
      # ---- FAB/GLOB colouring from merged CSV ----
      if (evo %in% c("fab","glob")) {
        fabsum <- fab_merged()  # from reactive above
        nE <- nrow(tr$edge)
        brcols <- rep("lightgrey", nE)
        
        # Safety: clamp any index set to available edges
        clamp <- function(idx) idx[idx >= 1 & idx <= nE]
        
        # 1) purple for rows where No.Sig.Nodes == 6
        if ("No.Sig.Nodes" %in% names(fabsum)) {
          idx_node <- clamp(which(fabsum$No.Sig.Nodes == 6))
          if (length(idx_node)) brcols[idx_node] <- "purple"
        }
        
        # 2) split where No.Sig.Betas == 6 by sign of Mean..Beta...BL..NZ
        if (all(c("No.Sig.Betas","Mean..Beta...BL..NZ") %in% names(fabsum))) {
          idx_beta <- which(fabsum$No.Sig.Betas == 6)
          if (length(idx_beta)) {
            vals <- fabsum$Mean..Beta...BL..NZ[idx_beta]
            idx_neg <- clamp(idx_beta[which(vals < 0)])
            idx_pos <- clamp(idx_beta[which(vals > 0)])
            if (length(idx_neg)) brcols[idx_neg] <- "blue"
            if (length(idx_pos)) brcols[idx_pos] <- "red"
          }
        }
        
        par(mar = c(1,1,1,1))
        ape::plot.phylo(tr, edge.color = brcols, edge.width = 2, show.tip.label = FALSE)
        legend("topleft",
               legend = c("node scalar", "positive beta", "negative beta"),
               col    = c("purple", "red", "blue"),
               lwd    = 3, bty = "n", cex = 0.9)
        
        # Little note as requested 🙂
        mtext("Magnitudes coming soon… but Jo needs to sleep.", side = 1, line = -1.5, adj = 1, cex = 0.8)
        
        return(invisible(NULL))
      }
      
      # ---- existing VR branch (unchanged) ----
      edge_cols <- NULL
      if (evo == "vr" && isTRUE(input$vr_color_branches)) {
        rate_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
        validate(need(file.exists(rate_path), paste0("VR rate summary not found: ", rate_path)))
        rt <- as.data.frame(readRDS(rate_path)$ratesummary)
        validate(need("medianscalar" %in% names(rt), "rates table must contain 'medianscalar'."))
        ms <- rt$medianscalar
        validate(need(all(ms > 0, na.rm = TRUE), "All 'medianscalar' must be > 0 to color by log10."))
        pal <- vr_palette_params()
        log_vals <- log10(round(ms, 2))
        idx <- findInterval(log_vals, pal$grid, all.inside = TRUE)
        cols_full <- pal$lut[idx]
        nE <- nrow(tr$edge)
        edge_cols <- if (length(cols_full) == nE + 1L) cols_full[-1L]
        else if (length(cols_full) == nE) cols_full
        else rep(cols_full, length.out = nE)
      }
      
      par(mar = c(1,1,1,1))
      if (!is.null(edge_cols)) ape::plot.phylo(tr, edge.color = edge_cols, edge.width = 2, show.tip.label = FALSE)
      else                     ape::plot.phylo(tr,                 edge.width = 2, show.tip.label = FALSE)
      
      invisible(NULL)
    }, error = function(e) {
      plot.new(); text(0.5, 0.5, labels = paste("Tree plot error:\n", e$message), cex = 0.9)
    })
  }, res = 96, height = 1000)
  
  
  output$tree_legend <- renderPlot({
    # Only for VR + single replicate
    req(tolower(input$reg_evo_model) == "vr", input$pred_dataset != "all replicates")
    
    # Pull the post-processed VR data (rates table after dropping first row)
    vd  <- vr_data()
    rt2 <- vd$rates2
    validate(need("medianscalar" %in% names(rt2), "Rates table missing 'medianscalar'."))
    
    ms <- rt2$medianscalar
    ms <- ms[is.finite(ms) & ms > 0]
    validate(need(length(ms) > 0, "No positive medianscalar values to show."))
    
    # Build colour grid in log10 space (same palette inputs as your tree)
    log_vals <- log10(round(ms, 2))
    rng    <- range(log_vals, na.rm = TRUE)
    grid_x <- seq(rng[1], rng[2], by = 0.01)
    
    cols_in <- if (isTRUE(input$vr_use_mid)) {
      c(input$vr_col_start, input$vr_col_mid, input$vr_col_end)
    } else {
      c(input$vr_col_start, input$vr_col_end)
    }
    lut <- grDevices::colorRampPalette(cols_in)(length(grid_x))
    
    # Tick labels in UN-logged space: min, median, 1, max (clamped to range)
    ticks_val <- c(min(ms), stats::median(ms), 1, max(ms))
    ticks_lab <- formatC(ticks_val, format = "fg", digits = 3, drop0trailing = TRUE)
    ticks_at  <- log10(ticks_val)
    ticks_at  <- pmax(min(grid_x), pmin(max(grid_x), ticks_at))  # clamp to gradient range
    
    # ---- GRID DRAWING ----
    grid::grid.newpage()
    vp <- grid::viewport(
      x = 0.5, y = 0.5,
      width  = grid::unit(0.92, "npc"),
      height = grid::unit(0.80, "npc"),
      just = c("center", "center"),
      xscale = range(grid_x), yscale = c(0, 1)
    )
    grid::pushViewport(vp)
    
    # Gradient bar
    grad <- matrix(lut, nrow = 1)
    grid::grid.raster(
      grad,
      x = grid::unit(0, "npc"),
      y = grid::unit(0.60, "npc"),
      width  = grid::unit(1, "npc"),
      height = grid::unit(0.38, "npc"),
      just = c("left", "center"),
      interpolate = TRUE
    )
    
    # Ticks + labels
    for (i in seq_along(ticks_at)) {
      grid::grid.lines(
        x = grid::unit(c(ticks_at[i], ticks_at[i]), "native"),
        y = grid::unit(c(0.42, 0.78), "npc"),
        gp = grid::gpar(col = "black", lwd = 1)
      )
      grid::grid.text(
        ticks_lab[i],
        x = grid::unit(ticks_at[i], "native"),
        y = grid::unit(0.22, "npc"),
        just = c("center", "center"),
        gp = grid::gpar(cex = 0.95)
      )
    }
    
    # Title
    grid::grid.text(
      "medianscalar",
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.08, "npc"),
      gp = grid::gpar(cex = 0.95)
    )
    
    grid::popViewport()
  })
  


  # ---- Tree file download (keep only filename + content here) ----
  output$download_tree_file <- downloadHandler(
    filename = function() {
      validate(
        need(input$reg_evo_model != "bm", "Tree download is not available for the 'bm' model."),
        need(input$pred_dataset != "all replicates",
             "Tree download is only available for a single replicate.")
      )

      evo <- tolower(input$reg_evo_model)
      if (evo == "vr") {
        region   <- input$reg_region
        regmodel <- input$reg_regr_model
        tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
        tree_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-medTree.trees", region, regmodel, tag))
      } else if (evo %in% c("glob","fab")) {
        tree_path <- file.path("data", "regions_tree.trees")
      } else {
        validate(need(FALSE, "Unsupported evolutionary model for tree download."))
      }

      validate(need(file.exists(tree_path), paste0("Tree file not found: ", tree_path)))
      sub("^data/", "", tree_path)  # suggested name = actual filename without "data/"
    },
    content = function(file) {
      evo <- tolower(input$reg_evo_model)
      if (evo == "vr") {
        region   <- input$reg_region
        regmodel <- input$reg_regr_model
        tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
        tree_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-medTree.trees", region, regmodel, tag))
      } else {
        tree_path <- file.path("data", "regions_tree.trees")
      }
      tr <- ape::read.nexus(tree_path)
      ape::write.nexus(tr, file = file)
    }
  )

  # ---- Second tree: magnitude / frequency highlighting ----
  output$tree_plot_alt <- renderPlot({
    # Only for VR and a single replicate
    req(tolower(input$reg_evo_model) == "vr", input$pred_dataset != "all replicates")
    
    vd <- vr_data()
    tr  <- vd$tree
    rt2 <- vd$rates2
    
    # Total iterations (for frequency threshold)
    denom_total <- if (!is.null(vd$posterior)) ncol(vd$posterior) - 1 else NA_integer_
    validate(need(is.finite(denom_total) && denom_total > 0, "Could not determine posterior iteration count."))
    
    # Edges & base colours
    nE <- nrow(tr$edge)
    selcols <- rep("lightgrey", nE)
    
    # Thresholds
    mag  <- input$alt_mag
    freq <- input$alt_freq / 100
    
    # Hits using rates2 (already -1)
    validate(need(all(c("medianscalar","n_scaled") %in% names(rt2)), "rates table must contain 'medianscalar' and 'n_scaled'."))
    hi_hits  <- which(rt2$medianscalar >  mag & (rt2$n_scaled / denom_total) > freq)
    low_hits <- which(rt2$medianscalar < (1 / mag) & (rt2$n_scaled / denom_total) > freq)
    
    # Clamp to available edges
    hi_hits  <- hi_hits[hi_hits <= nE]
    low_hits <- low_hits[low_hits <= nE]
    
    # Apply colours
    if (length(hi_hits))  selcols[hi_hits]  <- "red"
    if (length(low_hits)) selcols[low_hits] <- "blue"
    
    # Plot
    par(mar = c(1,1,1,1))
    ape::plot.phylo(tr, edge.color = selcols, show.tip.label = FALSE)
    invisible(NULL)
  }, res = 96, height = 1000)
  

}

shinyApp(ui, server)
