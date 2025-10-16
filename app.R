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

MODEL_LEVELS <- paste0("m", 1:9)
MODEL_TO_FORM <- c(
  m1 = "body",
  m2 = "body~grp",
  m3 = "body-Q",
  m4 = "ROB",
  m5 = "ROB-Q",
  m6 = "ROB+body",
  m7 = "ROB+body-Q",
  m8 = "ROB-Q+body",
  m9 = "ROB-Q+body-Q"
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
          choices = EVO_LEVELS,
          selected = "glob"
        ),
        selectInput(
          "reg_regr_model",
          "Regression model",
          choices = MODEL_LEVELS,
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

          # FAB/GLOB: Beta probability slider (only when single replicate)
          conditionalPanel(
            "(input.reg_evo_model == 'fab' || input.reg_evo_model == 'glob') && input.pred_dataset != 'all replicates'",
            sliderInput("fab_beta", "Beta probability", min = 0, max = 1, value = 0.95, step = 0.01)
          ),

          plotOutput("tree_plot", height = "1000px"),
          downloadButton("download_tree_file", "Download tree (.nex)")
        ),

        # Legend for FAB/GLOB colouring (only when fab/glob + single replicate)
        conditionalPanel(
          "(input.reg_evo_model == 'fab' || input.reg_evo_model == 'glob') && input.pred_dataset != 'all replicates'",
          plotOutput("fab_legend", height = "120px")
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
        title = paste0("Correlations for regions starting with ", anchor_name),
        x = "Correlation", y = NULL
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
    form <- MODEL_TO_FORM[[input$reg_regr_model]]
    if (is.null(form)) form <- "(unknown model)"
    tags$p(tags$strong("Regression formula: "), form)
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

    df <- if (is.data.frame(obj)) obj else as.data.frame(obj, stringsAsFactors = FALSE)

    # Make best effort to numeric-ify
    if (requireNamespace("readr", quietly = TRUE)) {
      suppressWarnings({
        df <- readr::type_convert(df, guess_integer = TRUE, trim_ws = TRUE, locale = readr::locale(decimal_mark = "."))
      })
    }
    df
  })

  extract_model_inset <- function(df, model_code) {
    if (is.null(df)) return(NULL)

    # Normalise names to lower, strip spaces
    cn <- tolower(gsub("\\s+", "", names(df)))

    # Two common shapes:
    #  (i) wide: columns like "beta1median","beta1pmcmc","r2median"
    #  (ii) long: a "parameter" column + statistic columns such as "median","pmcmc"
    get_stat <- function(df, beta_idx = NULL, stat = c("median","pmcmc","r2")) {
      stat <- match.arg(stat)

      # Attempt wide first
      wide_hit <- function(key) {
        j <- match(key, cn)
        if (!is.na(j)) suppressWarnings(as.numeric(df[[j]])) else NULL
      }

      # R2 (wide)
      if (stat == "r2") {
        out <- wide_hit("r2") %||% wide_hit("r2median") %||% wide_hit("medianr2")
        if (!is.null(out)) return(out[1])
      }

      # Beta (wide)
      if (!is.null(beta_idx)) {
        base <- paste0("beta", beta_idx)
        out  <- wide_hit(paste0(base,"median"))
        if (is.null(out)) out <- wide_hit(base)
        if (stat == "pmcmc") {
          out2 <- wide_hit(paste0(base,"pmcmc"))
          if (!is.null(out2)) return(out2[1])
        }
        if (!is.null(out)) return(out[1])
      }

      # Long format: look for a 'parameter' column and typical stat columns
      pcol <- which(tolower(names(df)) %in% c("parameter","param","name"))[1]
      if (!is.na(pcol)) {
        par <- as.character(df[[pcol]])
        nms <- tolower(names(df))
        medcol <- which(nms %in% c("median","med","mean"))[1]
        pmcol  <- which(nms %in% c("pmcmc","p.mcmc","pmc","p"))[1]
        r2row  <- grep("^r2$", tolower(gsub("\\s+","",par)))
        if (stat == "r2" && !is.na(medcol) && length(r2row)) return(suppressWarnings(as.numeric(df[r2row[1], medcol])))
        if (!is.null(beta_idx)) {
          brow <- grep(paste0("^beta\\s*",beta_idx,"$"), tolower(gsub("\\s+","",par)))
          if (length(brow)) {
            if (stat == "median" && !is.na(medcol)) return(suppressWarnings(as.numeric(df[brow[1], medcol])))
            if (stat == "pmcmc"  && !is.na(pmcol))  return(suppressWarnings(as.numeric(df[brow[1], pmcol])))
          }
        }
      }
      NA_real_
    }

    # Build text per model
    fmt <- function(x, digits=3) if (is.na(x)) "NA" else formatC(x, digits=digits, format="fg", drop0trailing=TRUE)

    txt <- switch(
      model_code,
      "m1" = {
        b1  <- get_stat(df, 1, "median"); p1 <- get_stat(df, 1, "pmcmc"); r2 <- get_stat(df, stat="r2")
        paste0("median slope = ", fmt(b1), " [pₓ=", fmt(p1), "]\nR² = ", fmt(r2))
      },
      "m2" = {
        r2 <- get_stat(df, stat="r2")
        paste0("R² = ", fmt(r2), "\n(group diffs coming soon)")
      },
      "m3" = {
        b1 <- get_stat(df,1,"median"); p1 <- get_stat(df,1,"pmcmc")
        b2 <- get_stat(df,2,"median"); p2 <- get_stat(df,2,"pmcmc")
        r2 <- get_stat(df, stat="r2")
        paste0("body slope = ", fmt(b1), " [pₓ=", fmt(p1), "]\n",
               "body quadratic = ", fmt(b2), " [pₓ=", fmt(p2), "]\n",
               "R² = ", fmt(r2))
      },
      "m4" = {
        b1 <- get_stat(df,1,"median"); r2 <- get_stat(df, stat="r2")
        paste0("ROB slope = ", fmt(b1), "\nR² = ", fmt(r2))
      },
      "m5" = {
        b1 <- get_stat(df,1,"median"); b2 <- get_stat(df,2,"median"); r2 <- get_stat(df, stat="r2")
        paste0("ROB slope = ", fmt(b1), "\nROB quadratic = ", fmt(b2), "\nR² = ", fmt(r2))
      },
      "m6" = {
        b1 <- get_stat(df,1,"median"); b2 <- get_stat(df,2,"median"); r2 <- get_stat(df, stat="r2")
        paste0("ROB slope = ", fmt(b1), "\nbody slope = ", fmt(b2), "\nR² = ", fmt(r2))
      },
      "m7" = {
        b1 <- get_stat(df,1,"median"); b2 <- get_stat(df,2,"median"); b3 <- get_stat(df,3,"median"); r2 <- get_stat(df, stat="r2")
        paste0("ROB slope = ", fmt(b1), "\nbody slope = ", fmt(b2), "\nbody quadratic = ", fmt(b3), "\nR² = ", fmt(r2))
      },
      "m8" = {
        b1 <- get_stat(df,1,"median"); b2 <- get_stat(df,2,"median"); b3 <- get_stat(df,3,"median"); r2 <- get_stat(df, stat="r2")
        paste0("ROB slope = ", fmt(b1), "\nROB quadratic = ", fmt(b2), "\nbody slope = ", fmt(b3), "\nR² = ", fmt(r2))
      },
      "m9" = {
        b1 <- get_stat(df,1,"median"); b2 <- get_stat(df,2,"median"); b3 <- get_stat(df,3,"median"); b4 <- get_stat(df,4,"median"); r2 <- get_stat(df, stat="r2")
        paste0("ROB slope = ", fmt(b1), "\nROB quadratic = ", fmt(b2), "\nbody slope = ", fmt(b3), "\nbody quadratic = ", fmt(b4), "\nR² = ", fmt(r2))
      },
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
    if (is.matrix(obj)) {
      df <- as.data.frame(obj, stringsAsFactors = FALSE)
    } else if (is.data.frame(obj)) {
      df <- obj
    } else {
      df <- as.data.frame(obj, stringsAsFactors = FALSE)
    }

    # 1) Use readr::type_convert if available
    if (requireNamespace("readr", quietly = TRUE)) {
      suppressWarnings({
        df <- readr::type_convert(df, guess_integer = TRUE, trim_ws = TRUE, locale = readr::locale(decimal_mark = "."))
      })
    }

    # 2) For character columns that are mostly numeric, coerce to numeric
    is_mostly_numeric <- function(x) {
      if (!is.character(x)) return(FALSE)
      sup <- suppressWarnings(as.numeric(x))
      sum(!is.na(sup)) >= (length(x) - 1L)
    }
    char_cols <- vapply(df, is.character, logical(1))
    mostly_num_cols <- names(df)[char_cols & vapply(df, is_mostly_numeric, logical(1))]
    for (nm in mostly_num_cols) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    }

    # 3) Round all numeric columns to 4 dp
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))

    df
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
    dt <- datatable(dat, rownames = TRUE, options = list(pageLength = 15, scrollX = TRUE))
    # Format all numeric cols to 4 dp for consistent display
    num_cols <- which(vapply(dat, is.numeric, logical(1)))
    if (length(num_cols)) {
      dt <- DT::formatRound(dt, columns = num_cols, digits = 4)
    }
    dt
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

  vr_rate_tbl <- reactive({
    req(tolower(input$reg_evo_model) == "vr", input$pred_dataset != "all replicates")
    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
    src <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
    validate(need(file.exists(src), paste0("VR rate summary not found: ", src)))
    as.data.frame(readRDS(src))
  })

  fab_rates <- reactive({
    req(tolower(input$reg_evo_model) %in% c("fab","glob"), input$pred_dataset != "all replicates")
    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    evo      <- tolower(input$reg_evo_model)
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))

    csv_path <- file.path("data", sprintf("%s_d1_%s-%s-%s-FabPP.csv", region, regmodel, evo, tag))
    validate(need(file.exists(csv_path), paste0("FabPP file not found: ", csv_path)))

    # Preserve original column names (with spaces/symbols)
    df <- utils::read.csv(csv_path, check.names = FALSE)
    validate(need(nrow(df) >= 2, "FabPP file has insufficient rows."))

    # 1) Drop first row per spec
    df <- df[-1, , drop = FALSE]

    # If first column looks like edge indices, use as rownames
    if (!is.null(df[[1]]) && all(grepl("^[0-9]+$", as.character(df[[1]])))) {
      rownames(df) <- as.character(df[[1]])
    }
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

    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))
    rate_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
    validate(need(file.exists(rate_path), paste0("VR rate summary not found: ", rate_path)))

    obj <- readRDS(rate_path)
    rt  <- if (is.data.frame(obj)) obj else if (is.list(obj) && !is.null(obj$ratesummary)) obj$ratesummary else obj
    rt  <- as.data.frame(rt, stringsAsFactors = FALSE)
    validate(need("medianscalar" %in% names(rt), "rates table must contain 'medianscalar'."))

    ms <- rt$medianscalar
    validate(need(all(ms > 0, na.rm = TRUE), "All 'medianscalar' must be > 0 for log10 scaling."))
    ms_rounded <- round(ms, 2)

    log_vals <- log10(ms_rounded)
    rng <- range(log_vals, finite = TRUE, na.rm = TRUE)
    grid <- seq(rng[1], rng[2], by = 0.01)

    cols_in <- if (isTRUE(input$vr_use_mid)) {
      c(input$vr_col_start, input$vr_col_mid, input$vr_col_end)
    } else {
      c(input$vr_col_start, input$vr_col_end)
    }
    lut <- grDevices::colorRampPalette(cols_in)(length(grid))  # same length used by findInterval mapping

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

      evo <- tolower(input$reg_evo_model)
      region   <- input$reg_region
      regmodel <- input$reg_regr_model
      tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))

      tree_path <- if (evo %in% c("glob","fab")) {
        file.path("data", "regions_tree.trees")
      } else if (evo == "vr") {
        file.path("data", sprintf("%s_d1_%s-VR-%s-medTree.trees", region, regmodel, tag))
      } else NA_character_

      validate(need(!is.na(tree_path) && file.exists(tree_path),
                    paste0("Tree file not found: ", tree_path)))

      tr <- ape::read.nexus(tree_path)

      # ... inside tryCatch, after you computed evo/region/regmodel/tag and read tr ...

      edge_cols <- NULL

      if (evo %in% c("glob","fab")) {
        # ---- FAB/GLOB: scale edge lengths and colour by Beta effects ----
        rt <- fab_rates()
        nE <- nrow(tr$edge)

        # Get "Median Non 1" for edge-length scaling
        col_mn1 <- "Median Non 1"
        validate(need(col_mn1 %in% names(rt),
                      paste0("FabPP table missing column '", col_mn1, "'.")))
        scale_vec <- rep(1, nE)

        # Map rows to edges by rownames if possible; otherwise fallback by position
        idx_map <- suppressWarnings(as.integer(rownames(rt)))
        if (all(is.finite(idx_map))) {
          idx_use <- pmin(nE, pmax(1, idx_map))
          scale_vec[idx_use] <- as.numeric(rt[[col_mn1]])
        } else if (nrow(rt) == nE) {
          scale_vec <- as.numeric(rt[[col_mn1]])
        } else {
          k <- min(nE, nrow(rt))
          scale_vec[seq_len(k)] <- as.numeric(rt[[col_mn1]][seq_len(k)])
        }

        # Apply edge-length scaling
        tr$edge.length <- tr$edge.length * scale_vec

        # --- Build colour vector from Mean(Beta*BL) with β probability filter ---
        selcols <- rep("lightgrey", nE)

        # Columns we need (names preserved by check.names=FALSE)
        col_p_lt0 <- if ("P < 0" %in% names(rt)) "P < 0" else "P<0"
        col_p_gt0 <- if ("P > 0" %in% names(rt)) "P > 0" else "P>0"
        col_mean  <- "Mean (Beta * BL) NZ"
        validate(need(col_mean %in% names(rt),
                      paste0("FabPP table missing column '", col_mean, "'.")))
        validate(need(all(c(col_p_lt0, col_p_gt0) %in% names(rt)),
                      "FabPP table missing 'P < 0' and/or 'P > 0' columns."))

        beta_thresh <- input$fab_beta
        prob_ok <- (rt[[col_p_lt0]] >= beta_thresh) | (rt[[col_p_gt0]] >= beta_thresh)

        vals <- as.numeric(rt[[col_mean]])
        sel_idx_rows <- which(prob_ok)
        vals_ok <- vals[sel_idx_rows]

        if (length(vals_ok) && any(vals_ok != 0, na.rm = TRUE)) {
          max_abs <- max(abs(vals_ok), na.rm = TRUE)
          n_grad <- 200
          blue_pal <- grDevices::colorRampPalette(c("lightblue", "#08306B"))(n_grad)  # neg
          red_pal  <- grDevices::colorRampPalette(c("mistyrose", "#7F0000"))(n_grad)  # pos

          # Edge index mapping for selected rows
          if (all(is.finite(idx_map))) {
            edge_idx <- pmin(nE, pmax(1, idx_map[sel_idx_rows]))
          } else {
            # fallback: assume row order aligns with edges
            edge_idx <- pmin(nE, pmax(1, sel_idx_rows))
          }

          for (k in seq_along(sel_idx_rows)) {
            v <- vals_ok[k]
            if (is.na(v) || v == 0) next
            id <- max(1, min(n_grad, ceiling(abs(v) / max_abs * n_grad)))
            selcols[edge_idx[k]] <- if (v < 0) blue_pal[id] else red_pal[id]
          }
        }

        edge_cols <- selcols
      }

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
      if (!is.null(edge_cols)) ape::plot.phylo(tr, edge.color = edge_cols, edge.width = 2, show.tip.label = F)
      else                     ape::plot.phylo(tr,                edge.width = 2, show.tip.label = F)

      invisible(NULL)
    }, error = function(e) {
      # draw the error in the plot area (prevents ggplot from being called)
      plot.new(); text(0.5, 0.5, labels = paste("Tree plot error:\n", e$message), cex = 0.9)
    })
  }, res = 96, height = 1000)

  output$tree_legend <- renderPlot({
    # Only for VR + single replicate
    req(tolower(input$reg_evo_model) == "vr", input$pred_dataset != "all replicates")

    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))

    # Load the rates and extract medianscalar (>0)
    rate_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
    validate(need(file.exists(rate_path), paste0("VR rate summary not found: ", rate_path)))

    obj <- readRDS(rate_path)
    rt  <- if (is.data.frame(obj)) obj else if (is.list(obj) && !is.null(obj$ratesummary)) obj$ratesummary else obj
    rt  <- as.data.frame(rt, stringsAsFactors = FALSE)
    validate(need("medianscalar" %in% names(rt), "Rates table missing 'medianscalar'."))

    ms <- rt$medianscalar
    ms <- ms[is.finite(ms) & ms > 0]
    validate(need(length(ms) > 0, "No positive medianscalar values to show."))

    # Build colour grid in log10 space (same palette inputs as your tree)
    log_vals <- log10(round(ms, 2))
    rng  <- range(log_vals, na.rm = TRUE)
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
    ticks_at  <- pmax(min(grid_x), pmin(max(grid_x), ticks_at))  # clamp

    # ---- GRID DRAWING (no base par/plot.new) ----
    grid::grid.newpage()

    # Inset the drawing area a bit (so labels don’t get clipped)
    vp <- grid::viewport(
      x = 0.5, y = 0.5,
      width  = grid::unit(0.92, "npc"),   # was 1 npc
      height = grid::unit(0.80, "npc"),   # was 1 npc
      just = c("center", "center"),
      xscale = range(grid_x), yscale = c(0, 1)
    )
    grid::pushViewport(vp)

    # Draw full-width gradient (now inside the inset viewport)
    grad <- matrix(lut, nrow = 1)
    grid::grid.raster(
      grad,
      x = grid::unit(0, "npc"),
      y = grid::unit(0.60, "npc"),        # was 0.55
      width  = grid::unit(1, "npc"),
      height = grid::unit(0.38, "npc"),   # was 0.35
      just = c("left", "center"),
      interpolate = TRUE
    )

    # Tick marks & labels (lifted up a bit)
    for (i in seq_along(ticks_at)) {
      grid::grid.lines(
        x = grid::unit(c(ticks_at[i], ticks_at[i]), "native"),
        y = grid::unit(c(0.42, 0.78), "npc"),      # was c(0.35, 0.70)
        gp = grid::gpar(col = "black", lwd = 1)
      )
      grid::grid.text(
        ticks_lab[i],
        x = grid::unit(ticks_at[i], "native"),
        y = grid::unit(0.22, "npc"),               # was 0.15
        just = c("center", "center"),
        gp = grid::gpar(cex = 0.95)
      )
    }

    # Axis title (lifted a touch)
    grid::grid.text(
      "medianscalar",
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.08, "npc"),                 # was 0.02
      gp = grid::gpar(cex = 0.95)
    )

    grid::popViewport()
  })

  output$fab_legend <- renderPlot({
    req(tolower(input$reg_evo_model) %in% c("fab","glob"), input$pred_dataset != "all replicates")
    rt <- fab_rates()

    # Use the same filter as the plot
    col_p_lt0 <- if ("P < 0" %in% names(rt)) "P < 0" else "P<0"
    col_p_gt0 <- if ("P > 0" %in% names(rt)) "P > 0" else "P>0"
    col_mean  <- "Mean (Beta * BL) NZ"
    validate(need(all(c(col_mean, col_p_lt0, col_p_gt0) %in% names(rt)),
                  "FabPP table missing required columns for legend."))

    beta_thresh <- input$fab_beta
    prob_ok <- (rt[[col_p_lt0]] >= beta_thresh) | (rt[[col_p_gt0]] >= beta_thresh)
    vals <- as.numeric(rt[[col_mean]])
    vals_ok <- vals[prob_ok]
    validate(need(length(vals_ok) > 0, "No edges pass the selected β probability."))

    max_abs <- max(abs(vals_ok), na.rm = TRUE)
    n_grad  <- 256
    blue_pal <- grDevices::colorRampPalette(c("lightblue", "#08306B"))(n_grad)  # neg
    red_pal  <- grDevices::colorRampPalette(c("mistyrose", "#7F0000"))(n_grad)  # pos

    # Draw with grid to avoid margin issues
    grid::grid.newpage()
    vp <- grid::viewport(x = 0.5, y = 0.5, width = grid::unit(0.95, "npc"), height = grid::unit(0.8, "npc"))
    grid::pushViewport(vp)

    # Left half (negative), right half (positive)
    grid::grid.raster(matrix(blue_pal, nrow = 1),
                      x = grid::unit(0.0, "npc"), y = grid::unit(0.6, "npc"),
                      width = grid::unit(0.49, "npc"), height = grid::unit(0.35, "npc"),
                      just = c("left","center"), interpolate = TRUE)
    grid::grid.raster(matrix(red_pal, nrow = 1),
                      x = grid::unit(0.51, "npc"), y = grid::unit(0.6, "npc"),
                      width = grid::unit(0.49, "npc"), height = grid::unit(0.35, "npc"),
                      just = c("left","center"), interpolate = TRUE)

    # Labels
    lab_neg <- paste0("−", formatC(max_abs, format = "fg", digits = 3))
    lab_pos <- formatC(max_abs, format = "fg", digits = 3)

    grid::grid.text(lab_neg, x = grid::unit(0.01, "npc"), y = grid::unit(0.20, "npc"), just = c("left","center"))
    grid::grid.text("0",      x = grid::unit(0.50, "npc"), y = grid::unit(0.20, "npc"), just = c("center","center"))
    grid::grid.text(lab_pos,  x = grid::unit(0.99, "npc"), y = grid::unit(0.20, "npc"), just = c("right","center"))
    grid::grid.text("Mean(Beta × BL)  (negative → blue, positive → red)",
                    x = grid::unit(0.5, "npc"), y = grid::unit(0.05, "npc"), just = c("center","center"))

    grid::popViewport()
  }, height = 120)

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

    region   <- input$reg_region
    regmodel <- input$reg_regr_model
    tag      <- sprintf("%03d", as.integer(gsub("[^0-9]", "", input$pred_dataset)))

    # Read tree
    tree_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-medTree.trees", region, regmodel, tag))
    validate(need(file.exists(tree_path), paste0("Tree file not found: ", tree_path)))
    tr <- ape::read.nexus(tree_path)

    # Read VR rate summary
    rate_path <- file.path("data", sprintf("%s_d1_%s-VR-%s-rateSummary.RDS", region, regmodel, tag))
    validate(need(file.exists(rate_path), paste0("VR rate summary not found: ", rate_path)))
    obj <- readRDS(rate_path)
    rt  <- obj$ratesummary

    # Total iterations (for frequency threshold)
    denom_total <- ncol(obj$posterior) - 1

    # Edges & base colours
    nE <- nrow(tr$edge)
    selcols <- rep("lightgrey", nE)

    # Remove first row of rate table
    rt <- rt[-1, , drop = FALSE]

    # Thresholds
    mag  <- input$alt_mag
    freq <- input$alt_freq / 100

    # Hits
    hi_hits  <- which(rt$medianscalar >  mag & rt$n_scaled / denom_total > freq)
    low_hits <- which(rt$medianscalar < (1 / mag) & rt$n_scaled / denom_total > freq)

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
