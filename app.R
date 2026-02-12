library(shiny)
library(tidyverse)
library(readxl)
library(cowplot)
library(ggprism)
library(rstatix)
library(tidyplots)
library(zip)
library(colourpicker)
library(emmeans)

# Anchor defaults (exactly from your FINAL SCRIPT)
MY_DEFAULT_COLORS <- c("#737373", "#00BF7D", "#76008F", "#FF8800")

ui <- fluidPage(
  titlePanel("Plot your flow cytometry data from Excel: FAST!"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "1. Upload Excel File", accept = c(".xls", ".xlsx")),
      uiOutput("sheet_selector"),
      tags$hr(),
      
      # NEW: choose ONE grouping column
      uiOutput("grouping_selector"),
      actionButton("reset_grouping", "Reset grouping selection", class = "btn-default btn-sm"),
      tags$hr(),
      
      tags$hr(),
      tags$b("Statistics"),
      selectInput(
        "ref_group",
        "Compare all groups vs:",
        choices = c("All-by-all (default)" = "__ALL__"),
        selected = "__ALL__"
      ),
      
      uiOutput("var_selector"),
      actionButton("select_all_vars", "Select all variables", class = "btn-default btn-sm"),
      actionButton("deselect_all_vars", "Deselect all variables", class = "btn-default btn-sm"),
      checkboxInput("show_legend", "Display Legend", value = FALSE),
      
      tags$hr(),
      tags$b("Preview"),
      selectInput("preview_var", "Preview variable:", choices = character(0)),
      
      tags$hr(),
      tags$b("Dimensions & Aesthetics"),
      sliderInput("bar_width", "Width per Bar (cm)", min = 0.3, max = 1.5, value = 0.6, step = 0.1),
      sliderInput("plot_height", "Plot Height (cm)", min = 3, max = 15, value = 6),
      
      textInput(
        "colors",
        "Anchor default colors (reference)",
        value = paste(MY_DEFAULT_COLORS, collapse = ", ")
      ),
      
      checkboxInput("use_custom_group_colors", "Custom group colors (override defaults)", value = FALSE),
      uiOutput("group_color_pickers"),
      
      tags$hr(),
      downloadButton("download_combined_pdf", "Download Selected (PDF + Legend ZIP)", class = "btn-success"),
      tags$br(), tags$br(),
      downloadButton("download_all_png", "Download Selected (PNG ZIP + Legend)", class = "btn-primary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Preview", plotOutput("big_plot", height = "1200px", width = '1200px')),
        tabPanel("Data Preview", tableOutput("table_preview")),
        tabPanel("Help & Formatting",
                 fluidRow(
                   column(12,
                          h3("Excel Formatting Guide"),
                          tableOutput("example_table"),
                          p("Select ONE grouping column. Numeric columns are plotted. Non significant comparison are hidden.")
                   )
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Draw a plot/grob inside a smaller canvas to create padding (prevents PDF cropping)
  pad_right <- function(p, right_frac = 0.04, left_frac = 0.00, top_frac = 0.00, bottom_frac = 0.00) {
    # right_frac = 0.04 means keep 4% blank space on the right
    cowplot::ggdraw() +
      cowplot::draw_plot(
        p,
        x = left_frac,
        y = bottom_frac,
        width  = 1 - left_frac - right_frac,
        height = 1 - top_frac - bottom_frac,
        hjust = 0, vjust = 0
      )
  }
  
  is_id_like <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) return(TRUE)
    
    # High uniqueness → likely an ID
    uniq_ratio <- length(unique(x)) / length(x)
    
    # Mostly numeric strings → likely an ID
    numeric_ratio <- mean(!is.na(suppressWarnings(as.numeric(x))))
    
    uniq_ratio > 0.8 || numeric_ratio > 0.8
  }
  
  output$example_table <- renderTable({
    tibble(ID = "S1", Genotype = "WT", Condition = "Ctrl",
           B_Cells = "15.5%", T_Cells = 45.2)
  })
  
  output$sheet_selector <- renderUI({
    req(input$file)
    sheet_names <- excel_sheets(input$file$datapath)
    selectInput("sheet_select", "2. Select Sheet:", choices = sheet_names)
  })
  
  # Read raw sheet once
  raw_df <- reactive({
    req(input$file, input$sheet_select)
    read_excel(input$file$datapath, sheet = input$sheet_select) %>%
      filter(!is.na(.[[1]]))
  })
  
  numeric_cols_from_raw <- reactive({
    req(raw_df())
    df <- raw_df()
    
    # TRUE numeric columns as read from Excel
    raw_num <- names(df)[sapply(df, is.numeric)]
    
    # Also treat columns that are "mostly numeric strings" as numeric (common with % and commas)
    mostly_numeric <- names(df)[sapply(df, function(x) {
      x <- as.character(x)
      x <- stringr::str_replace_all(x, "%", "")
      x <- stringr::str_replace_all(x, ",", ".")
      x <- stringr::str_trim(x)
      x <- x[!is.na(x) & x != ""]
      if (length(x) == 0) return(FALSE)
      mean(!is.na(suppressWarnings(as.numeric(x)))) >= 0.8
    })]
    
    unique(c(raw_num, mostly_numeric))
  })
  
  # Detect candidate grouping columns (non-numeric)
  non_numeric_cols <- reactive({
    req(raw_df())
    df <- raw_df()
    
    is_numeric_like_col <- function(x, threshold = 0.8) {
      x <- as.character(x)
      x <- stringr::str_replace_all(x, "\u00A0", " ")
      x <- stringr::str_squish(x)
      x <- x[!is.na(x) & x != ""]
      if (length(x) == 0) return(FALSE)
      
      # Remove common numeric decorations
      x_clean <- x %>%
        stringr::str_replace_all("%", "") %>%
        stringr::str_replace_all(",", ".") %>%
        stringr::str_replace_all(" ", "")
      
      # If most values become numeric, it's numeric-like (so NOT a grouping column)
      frac_numeric <- mean(!is.na(suppressWarnings(as.numeric(x_clean))))
      frac_numeric >= threshold
    }
    
    # candidates: character/factor columns ONLY
    char_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    
    # keep only those that are NOT numeric-like
    keep <- char_cols[!sapply(df[char_cols], is_numeric_like_col)]
    
    keep
  })
  
  # ============================================================================
  # EDIT 1) Replace your clean_groups() with the "nuclear" cleaner
  # ============================================================================
  clean_groups <- function(x) {
    x <- as.character(x)
    x <- gsub("[\u200b\u200c\u200d\ufeff\u00A0]", "", x, perl = TRUE) # remove invisibles + NBSP
    x <- stringr::str_trim(x)
    
    keep <- !is.na(x) &
      x != "" &
      x != "NA" &
      x != "_" &
      x != "-" &
      grepl("\\S", x, perl = TRUE)
    
    x <- x[keep]
    unique(x)
  }
  
  # UI to choose ONE grouping column (no collapsing/pasting)
  output$grouping_selector <- renderUI({
    req(raw_df())
    df <- raw_df()
    nn <- non_numeric_cols()
    
    if (length(nn) == 0) {
      return(helpText("No character columns found for grouping. Using Group = 'All'."))
    }
    
    good_group_cols <- nn[!sapply(df[nn], is_id_like)]
    default_sel <- if (length(good_group_cols) > 0) good_group_cols[1] else nn[1]
    
    checkboxGroupInput(
      "group_cols",
      "Grouping columns (choose one or more):",
      choices = nn,
      selected = default_sel
    )
  })
  
  observeEvent(input$reset_grouping, {
    req(raw_df())
    df <- raw_df()
    nn <- non_numeric_cols()
    
    good_group_cols <- nn[!sapply(df[nn], is_id_like)]
    default_sel <- if (length(good_group_cols) > 0) good_group_cols[1] else nn[1]
    
    updateCheckboxGroupInput(session, "group_cols", selected = default_sel)
  })
  
  # ============================================================================
  # EDIT 2) In processed_data(): aggressively clean the Group COLUMN ITSELF
  #         and remove "ghost" rows BEFORE factor/droplevels
  # ============================================================================
  processed_data <- reactive({
    req(raw_df())
    df <- raw_df()
    
    group_cols <- input$group_cols
    if (is.null(group_cols)) group_cols <- character(0)
    group_cols <- intersect(group_cols, names(df))
    
    if (length(group_cols) == 0) {
      df <- df %>% mutate(Group = "All")
    } else if (length(group_cols) == 1) {
      df <- df %>% mutate(Group = as.character(.data[[group_cols[1]]]))
    } else {
      df <- df %>% unite("Group", all_of(group_cols), sep = "_", remove = FALSE, na.rm = TRUE)
    }
    
    # --- APPLY NUCLEAR CLEANING TO THE COLUMN + FILTER ROWS ---
    df$Group <- as.character(df$Group)
    df$Group <- gsub("[\u200b\u200c\u200d\ufeff\u00A0]", "", df$Group, perl = TRUE)
    df$Group <- stringr::str_trim(df$Group)
    
    df <- df %>%
      filter(
        !is.na(Group),
        Group != "",
        Group != "NA",
        !Group %in% c("_", "__", "-", "--"),
        grepl("\\S", Group, perl = TRUE)
      )
    
    df$Group <- droplevels(factor(df$Group))
    # ---------------------------------------------------------
    
    protected <- unique(c("Group", group_cols))
    candidate_cols <- setdiff(names(df), protected)
    
    df <- df %>%
      mutate(across(all_of(candidate_cols), ~{
        x <- as.character(.x)
        x <- str_replace_all(x, "%", "")
        x <- str_replace_all(x, ",", ".")
        suppressWarnings(as.numeric(x))
      }))
    
    num_cols <- names(df)[sapply(df, is.numeric)]
    if (length(num_cols) > 0) {
      df <- df %>% filter(if_any(all_of(num_cols), ~ !is.na(.x)))
    }
    
    df
  })
  
  observeEvent(processed_data(), {
    df_proc <- processed_data()
    allowed <- numeric_cols_from_raw()
    new_choices <- intersect(allowed, names(df_proc)[sapply(df_proc, is.numeric)])
    
    old_sel <- isolate(input$selected_vars)
    if (is.null(old_sel)) old_sel <- character(0)
    kept <- intersect(old_sel, new_choices)
    if (length(kept) == 0 && length(new_choices) > 0) kept <- new_choices[1]
    
    updateCheckboxGroupInput(
      session,
      "selected_vars",
      choices = new_choices,
      selected = kept
    )
  }, ignoreInit = TRUE)
  
  observeEvent(groups_in_data(), {
    grps <- as.character(groups_in_data())
    grps <- grps[!is.na(grps) & grps != ""]
    if (length(grps) == 0) return()
    
    cur <- isolate(input$ref_group)
    if (is.null(cur)) cur <- "__ALL__"
    if (cur != "__ALL__" && !(cur %in% grps)) cur <- "__ALL__"
    
    updateSelectInput(
      session,
      "ref_group",
      choices = c("All-by-all (default)" = "__ALL__", setNames(grps, grps)),
      selected = cur
    )
  }, ignoreInit = TRUE)
  
  output$var_selector <- renderUI({
    checkboxGroupInput(
      "selected_vars",
      "3. Select Variables:",
      choices = character(0),
      selected = character(0)
    )
  })
  
  observeEvent(input$select_all_vars, {
    req(processed_data())
    df <- processed_data()
    cols <- names(df)[sapply(df, is.numeric)]
    updateCheckboxGroupInput(session, "selected_vars", selected = cols)
  })
  
  observeEvent(input$deselect_all_vars, {
    updateCheckboxGroupInput(session, "selected_vars", selected = character(0))
  })
  
  observeEvent(list(processed_data(), input$selected_vars), {
    req(processed_data())
    cols <- input$selected_vars
    if (is.null(cols) || length(cols) == 0) return()
    cur <- isolate(input$preview_var)
    if (is.null(cur) || !(cur %in% cols)) cur <- cols[[1]]
    updateSelectInput(session, "preview_var", choices = cols, selected = cur)
  }, ignoreInit = TRUE)
  
  to_gg <- function(tp_obj) {
    if (is.ggplot(tp_obj)) return(tp_obj)
    plot(tp_obj)
  }
  
  groups_in_data <- reactive({
    req(processed_data())
    clean_groups(processed_data()$Group)
  })
  
  default_palette_n <- function(n) {
    anchors <- MY_DEFAULT_COLORS
    if (n <= length(anchors)) return(anchors[seq_len(n)])
    extra_n <- n - length(anchors)
    extra <- grDevices::hcl.colors(extra_n, palette = "Dark 3")
    c(anchors, extra)
  }
  
  default_palette_vec <- reactive({
    req(groups_in_data())
    default_palette_n(length(groups_in_data()))
  })
  
  output$group_color_pickers <- renderUI({
    req(groups_in_data())
    if (!isTRUE(input$use_custom_group_colors)) return(NULL)
    
    grps <- as.character(groups_in_data())
    base_pal <- default_palette_vec()
    
    tagList(
      tags$hr(),
      tags$b("Custom colors per group"),
      lapply(seq_along(grps), function(i) {
        gid <- paste0("col_", make.names(grps[i]))
        colourInput(
          inputId = gid,
          label   = grps[i],
          value   = base_pal[i],
          allowTransparent = FALSE
        )
      })
    )
  })
  
  group_palette <- reactive({
    req(groups_in_data())
    grps <- levels(processed_data()$Group)
    n <- length(grps)
    base_pal <- default_palette_n(n)
    
    if (isTRUE(input$use_custom_group_colors)) {
      cols <- vapply(grps, function(g) {
        id <- paste0("col_", make.names(g))
        val <- input[[id]]
        if (is.null(val) || val == "") NA_character_ else val
      }, character(1))
      cols[is.na(cols)] <- base_pal[is.na(cols)]
      names(cols) <- grps
      cols
    } else {
      names(base_pal) <- grps
      base_pal
    }
  })
  
  # ============================================================================
  # EDIT 3) Replace legend_plot(): build a manual legend (no ggplot legend key)
  #         This removes the "extra empty box" permanently.
  # ============================================================================
  # --- helper: compute a good width/height for a vertical legend (prevents cropping) ---
  legend_dims_vertical <- function(groups, text_pt = 11, key_cm = 0.42) {
    g <- clean_groups(groups)
    if (length(g) == 0) return(list(w_cm = 6, h_cm = 3))
    
    max_label <- g[which.max(nchar(g))]
    
    # open null device for safe measurement (prevents device errors)
    grDevices::pdf(NULL)
    
    # set font size for measurement
    grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = text_pt)))
    
    txt_w_in <- grid::convertWidth(
      grid::stringWidth(max_label),
      unitTo = "in",
      valueOnly = TRUE
    )
    
    grid::popViewport()
    grDevices::dev.off()
    
    txt_w_cm <- txt_w_in * 2.54
    
    # total width = left pad + key + gap + text + right pad
    w_cm <- 0.8 + key_cm + 0.6 + txt_w_cm + 1.2
    
    # height per row + padding
    row_h_cm <- 0.60
    h_cm <- 0.8 + row_h_cm * length(g) + 0.8
    
    list(w_cm = w_cm, h_cm = h_cm)
  }
  
  # --- manual vertical legend: square keys + no ggplot legend (so no “extra box”) ---
  legend_plot_vertical <- function(groups, pal_named,
                                   key_cm = 0.42,
                                   text_pt = 8) {
    groups <- clean_groups(groups)
    if (length(groups) == 0) stop("No valid groups found after cleaning.")
    
    cols <- pal_named[groups]
    if (any(is.na(cols))) {
      fb <- default_palette_n(length(groups))
      cols[is.na(cols)] <- fb[is.na(cols)]
    }
    cols <- unname(cols)
    
    df_leg <- tibble(
      label = as.character(groups),
      fill  = cols,
      y     = rev(seq_along(groups))
    ) %>%
      mutate(
        x_key  = 0,
        x_text = key_cm/2 + 0.30
      )
    
    max_chars <- max(nchar(df_leg$label))
    x_right <- df_leg$x_text[1] + max_chars * 0.18  # a bit more generous than 0.14
    
    ggplot(df_leg, aes(y = y)) +
      geom_tile(aes(x = x_key, fill = fill),
                width = key_cm, height = key_cm, color = "black") +
      geom_text(aes(x = x_text, label = label),
                hjust = 0, vjust = 0.5, color = "black") +
      scale_fill_identity() +
      scale_x_continuous(
        limits = c(-key_cm, x_right),
        expand = expansion(mult = c(0, 0))
      ) +
      coord_cartesian(clip = "off") +        # <-- IMPORTANT
      theme_void() +
      theme(
        text = element_text(size = as.numeric(text_pt)),
        plot.margin = margin(t = 8, r = 120, b = 8, l = 8, unit = "pt")  # <-- more right room
      )
  }
  
  # ============================================================================
  
  # ---- plot_generator unchanged (kept as-is) ----
  plot_generator <- function(var_name, force_no_legend = FALSE) {
    df <- processed_data()
    pal_named <- group_palette()
    n_groups <- length(levels(df$Group))
    fixed_width <- n_groups * input$bar_width
    
    tidyplots_options(width = fixed_width, height = input$plot_height, unit = "cm")
    
    df %>%
      select(Group, valor = all_of(var_name)) %>%
      mutate(Group = droplevels(Group)) %>%
      tidyplot(x = Group, y = valor, color = Group) %>%
      add(stat_summary(fun = mean, geom = "bar", color = "black", width = 0.8, linewidth = 0.3) ) %>%
      add_sem_errorbar(color = "black", linewidth = 0.3) %>%
      add(geom_jitter(position = position_jitter(width = 0.15, height = 0),
                      size = 1.75, color = "black", fill = "white", shape = 21)) %>%
      adjust_colors(unname(pal_named)) %>%
      add(coord_cartesian( clip = "off")) %>%
      adjust_y_axis_title(var_name) %>%
      adjust_x_axis_title("") %>%
      add_title(word(var_name, 1)) %>%
      add(scale_y_continuous(
        guide = "prism_offset_minor", 
        expand = expansion(mult = c(0, 0.075))
      )) %>% 
      add(theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        axis.ticks.y = element_line(linewidth = 0.3),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        
        # IMPORTANT: extra right margin prevents clipping of p-values etc.
        plot.margin = margin(t = 10, r = 60, b = 10, l = 40, unit = "pt"),
        
        # IMPORTANT: never allow per-plot legend in exports
        legend.position = if (force_no_legend) "none" else if (isTRUE(input$show_legend)) "right" else "none"
      )) -> p
    
    # ... keep the stats section exactly as you already have it ...
    
    p <- p %>%
      add_reference_lines(y = 0, linetype = "solid", linewidth = 0.3) %>%
      add(scale_y_continuous(
        limits = c(0, NA),
        guide  = "prism_offset_minor",
        expand = expansion(mult = c(0, 0.075))
      )) %>%
      add(theme = theme(
        axis.line.y = element_line(colour = "black", linewidth = 0.3),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
      ))
    
    return(p)
  }
  
  output$big_plot <- renderPlot({
    req(input$preview_var)
    print(plot_generator(input$preview_var))
  }, res = 300)
  
  output$table_preview <- renderTable({
    req(processed_data())
    head(processed_data() %>% select(any_of(input$group_cols), Group, everything()), 10)
  })
  
  output$download_combined_pdf <- downloadHandler(
    filename = function() { paste0("PDF_and_Legend_", input$sheet_select, ".zip") },
    content = function(file) {
      req(input$selected_vars)
      
      tmpdir <- tempfile("pdfpkg_")
      dir.create(tmpdir, recursive = TRUE)
      
      # ----- A) Build plots with legends FORCED OFF -----
      plot_list <- lapply(input$selected_vars, function(v) to_gg(plot_generator(v, force_no_legend = TRUE)))
      num_plots <- length(plot_list)
      num_cols  <- ceiling(sqrt(num_plots))
      num_rows  <- ceiling(num_plots / num_cols)
      
      grid_plots <- plot_grid(plotlist = plot_list, ncol = num_cols, align = "hv", axis = "tblr")
      
      # ----- B) Build ONE legend (top) -----
      groups    <- unique(processed_data()$Group)
      pal_named <- group_palette()
      
      leg_plot <- legend_plot_vertical(groups, pal_named, key_cm = 0.42, text_pt = 11)
      
      # put legend on top (as one block) then plots below
      final_pdf <- plot_grid(
        leg_plot,
        grid_plots,
        ncol = 1,
        rel_heights = c(0.18, 0.82)
      )
      
      # ----- C) Save combined PDF with extra right padding to avoid page-edge clipping -----
      n_groups <- length(unique(processed_data()$Group))
      single_w <- (n_groups * input$bar_width) + 4
      single_h <- input$plot_height + 3
      
      plots_w_cm <- num_cols * single_w
      plots_h_cm <- num_rows * single_h
      
      combined_path <- file.path(tmpdir, paste0("Combined_", input$sheet_select, ".pdf"))
      
      ggsave(
        combined_path, final_pdf,
        width  = plots_w_cm + 2,   # <-- extra width buffer
        height = plots_h_cm + 4,   # <-- extra height buffer for legend
        units = "cm", limitsize = FALSE, device = cairo_pdf, bg = "transparent"
      )
      
      # ----- D) Export legend ALSO as separate PDF (like PNG workflow) -----
      legend_pdf_path <- file.path(tmpdir, paste0("Legend_", input$sheet_select, ".pdf"))
      dims <- legend_dims_vertical(groups, text_pt = 11, key_cm = 0.42)
      
      ggsave(
        legend_pdf_path, leg_plot,
        width = dims$w_cm, height = dims$h_cm, units = "cm",
        bg = "transparent", device = cairo_pdf, limitsize = FALSE
      )
      
      # zip
      zip::zipr(zipfile = file, files = list.files(tmpdir, full.names = TRUE))
    }
  )
  
  output$download_all_png <- downloadHandler(
    filename = function() { paste0("PNGs_and_Legend_", input$sheet_select, ".zip") },
    content = function(file) {
      req(input$selected_vars)
      
      tmpdir <- tempfile("pngs_")
      dir.create(tmpdir, recursive = TRUE)
      
      n_groups <- length(unique(processed_data()$Group))
      w_cm <- (n_groups * input$bar_width) + (if (input$show_legend) 7 else 4)
      h_cm <- input$plot_height + 3
      
      for (v in input$selected_vars) {
        safe <- gsub("[^A-Za-z0-9_\\-]", "_", v)
        out_png <- file.path(tmpdir, paste0(safe, ".png"))
        p <- to_gg(plot_generator(v))
        ggsave(out_png, plot = p, width = w_cm, height = h_cm, units = "cm", dpi = 600, bg = "transparent")
      }
      
      groups <- unique(processed_data()$Group)
      pal_named <- group_palette()
      
      lp <- legend_plot_vertical(groups, pal_named, key_cm = 0.42, text_pt = 11)
      dims <- legend_dims_vertical(groups)
      
      legend_png_path <- file.path(tmpdir, "Legend.png")
      
      ggsave(
        legend_png_path, lp,
        width = dims$w_cm, height = dims$h_cm, units = "cm",
        dpi = 600, bg = "transparent", limitsize = FALSE
      )
      
      zip::zipr(zipfile = file, files = list.files(tmpdir, full.names = TRUE))
    }
  )
}

shinyApp(ui = ui, server = server)