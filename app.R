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
  
  clean_groups <- function(x) {
    x <- as.character(x)
    x <- str_replace_all(x, "\u00A0", " ")
    x <- str_squish(x)
    x <- x[!is.na(x) & x != "" & x != "NA"]
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
    
    # Identify non-ID-like grouping candidates
    good_group_cols <- nn[!sapply(df[nn], is_id_like)]
    
    # Default selection logic:
    # 1) first non-ID-like column if available
    # 2) otherwise fall back to first character column
    default_sel <- if (length(good_group_cols) > 0) {
      good_group_cols[1]
    } else {
      nn[1]
    }
    
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
    default_sel <- if (length(good_group_cols) > 0) {
      good_group_cols[1]
    } else {
      nn[1]
    }
    
    updateCheckboxGroupInput(session, "group_cols", selected = default_sel)
  })
  
  # Processed data: keep character columns as-is; Group comes ONLY from selected column
  processed_data <- reactive({
    req(raw_df())
    df <- raw_df()
    
    group_cols <- input$group_cols
    if (is.null(group_cols)) group_cols <- character(0)
    group_cols <- intersect(group_cols, names(df))  # safety
    
    # Build Group AFTER user selection:
    if (length(group_cols) == 0) {
      df <- df %>% mutate(Group = "All")
    } else if (length(group_cols) == 1) {
      df <- df %>% mutate(Group = as.character(.data[[group_cols[1]]]))
    } else {
      # paste selected columns with "_"
      df <- df %>%
        unite("Group", all_of(group_cols), sep = "_", remove = FALSE, na.rm = TRUE)
    }
    
    # Clean Group
    df <- df %>%
      mutate(
        Group = as.character(Group),
        Group = str_replace_all(Group, "\u00A0", " "),
        Group = str_squish(Group)
      ) %>%
      filter(!is.na(Group) & Group != "" & Group != "NA")
    
    # Convert all columns NOT involved in grouping (and not Group) to numeric where possible
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
    
    df$Group <- droplevels(as.factor(df$Group))
    df
  })
  
  observeEvent(processed_data(), {
    df_proc <- processed_data()
    
    # numeric columns allowed based on raw excel (prevents char columns creeping in)
    allowed <- numeric_cols_from_raw()
    
    # only keep columns that exist in processed_data and are numeric there
    new_choices <- intersect(allowed, names(df_proc)[sapply(df_proc, is.numeric)])
    
    old_sel <- isolate(input$selected_vars)
    if (is.null(old_sel)) old_sel <- character(0)
    
    kept <- intersect(old_sel, new_choices)
    
    # default = first numeric column only (your preference)
    if (length(kept) == 0 && length(new_choices) > 0) {
      kept <- new_choices[1]
    }
    
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
    
    # preserve current selection if possible
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
    # render the control once; choices/selection will be managed by observers below
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
    updateCheckboxGroupInput(
      session,
      "selected_vars",
      selected = character(0)
    )
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
  
  # Default palette: exactly N colors, never recycles
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
  
  legend_plot <- function(groups, pal_named, horizontal = TRUE) {
    groups <- clean_groups(groups)
    pal <- pal_named[groups]
    if (any(is.na(pal))) pal[is.na(pal)] <- default_palette_n(length(groups))[is.na(pal)]
    pal <- unname(pal); names(pal) <- groups
    
    df_leg <- tibble(Group = factor(groups, levels = groups), x = 1, y = 1)
    
    ggplot(df_leg, aes(x = x, y = y, fill = Group)) +
      geom_point(shape = 22, size = 6, color = "black") +
      scale_fill_manual(values = pal, breaks = groups, limits = groups, drop = TRUE, na.translate = FALSE) +
      theme_void() +
      theme(
        legend.title = element_blank(),
        legend.position = if (horizontal) "top" else "right",
        legend.direction = if (horizontal) "horizontal" else "vertical",
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0, 0, 0, 0),
        legend.justification = "center"
      ) +
      guides(fill = guide_legend(
        nrow = if (horizontal && length(groups) > 3) 2 else 1,
        byrow = TRUE
      ))
  }
  
  plot_generator <- function(var_name) {
    df <- processed_data()
    pal_named <- group_palette()
    n_groups <- length(levels(df$Group))
    fixed_width <- n_groups * input$bar_width
    
    tidyplots_options(width = fixed_width, height = input$plot_height, unit = "cm")
    
    df %>%
      select(Group, valor = all_of(var_name)) %>%
      mutate(Group = droplevels(Group)) %>%
      tidyplot(x = Group, y = valor, color = Group) %>%
      add_reference_lines(y = 0, linetype = "solid", linewidth = 0.3) %>%
      add(stat_summary(fun = mean, geom = "bar", color = "black", width = 0.8, linewidth = 0.3) ) %>%
      add_sem_errorbar(color = "black", linewidth = 0.3) %>%
      add(geom_jitter(position = position_jitter(width = 0.15, height = 0),
                      size = 1.75, color = "black", fill = "white", shape = 21)) %>%
      adjust_colors(unname(pal_named)) %>%
      # add_test_pvalue(
      #   method = "tukey_hsd",
      #   hide_info = TRUE,
      #   hide.ns = TRUE,
      #   step.increase = 0.15,
      #   label = "{p.adj.signif}"
      # ) %>% 
      add(scale_y_continuous(guide = "prism_offset_minor",
                             expand = expansion(mult = c(0, 0.075)))) %>%
      adjust_y_axis_title(var_name) %>%
      adjust_x_axis_title("") %>%
      add_title(word(var_name, 1)) %>%
      add(theme = theme(
        axis.line = element_line(colour = "black", linewidth = 0.3),
        axis.ticks.y = element_line(linewidth = 0.3),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 40, unit = "pt"),
        legend.position = if (isTRUE(input$show_legend)) "right" else "none"
      )) -> p
    
    # --- Stats: all-by-all by default, or vs selected reference group ---
    ref <- input$ref_group
    if (is.null(ref)) ref <- "__ALL__"
    
    n_groups <- nlevels(df$Group)
    
    if (n_groups >= 3) {
      
      if (ref == "__ALL__") {
        
        # ALL-BY-ALL (Tukey)
        p <- p %>% add_test_pvalue(
          method = "tukey_hsd",
          hide_info = TRUE,
          hide.ns = T,
          step.increase = 0.15,
          label = "{p.adj.signif}"
        )
        
      } else {
        
        # MANY-vs-ONE using emmeans contrasts (treatment vs control)
        # IMPORTANT: emmeans_test is supported by geom_pwc
        p <- p %>% add_test_pvalue(
          method = "emmeans_test",
          ref.group = ref,
          hide_info = TRUE,
          hide.ns = T,
          step.increase = 0.15,
          label = "{p.adj.signif}"
        )
        
      }
      
    } else if (n_groups == 2) {
      
      # Two groups: t-test (has p, not p.adj)
      p <- p %>% add_test_pvalue(
        method = "t_test",
        hide_info = TRUE,
        hide.ns = T,
        label = "{p.signif}"
      )
    }
    
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
      
      plot_list <- lapply(input$selected_vars, function(v) to_gg(plot_generator(v)))
      num_plots <- length(plot_list)
      num_cols <- ceiling(sqrt(num_plots))
      num_rows <- ceiling(num_plots / num_cols)
      
      grid_plots <- plot_grid(plotlist = plot_list, ncol = num_cols, align = "hv", axis = "tblr")
      
      final_pdf <- grid_plots
      if (isTRUE(input$show_legend)) {
        groups <- unique(processed_data()$Group)
        pal_named <- group_palette()
        lp <- legend_plot(groups, pal_named, horizontal = TRUE)
        leg_grob <- cowplot::get_legend(lp)
        final_pdf <- plot_grid(NULL, leg_grob, NULL, grid_plots, ncol = 1,
                               rel_heights = c(0.03, 0.12, 0.03, 0.82))
      }
      
      n_groups <- length(unique(processed_data()$Group))
      single_w <- (n_groups * input$bar_width) + 4
      single_h <- input$plot_height + 3
      
      combined_path <- file.path(tmpdir, paste0("Combined_", input$sheet_select, ".pdf"))
      ggsave(combined_path, final_pdf,
             width = num_cols * single_w,
             height = (num_rows * single_h) + (if (input$show_legend) 4 else 0),
             units = "cm", limitsize = FALSE, device = cairo_pdf)
      
      groups <- unique(processed_data()$Group)
      pal_named <- group_palette()
      legend_pdf_path <- file.path(tmpdir, paste0("Legend_", input$sheet_select, ".pdf"))
      lp <- legend_plot(groups, pal_named, horizontal = TRUE)
      
      leg_w <- 22
      leg_h <- if (length(groups) > 3) 2.2 else 1.4
      
      ggsave(legend_pdf_path, lp,
             width = leg_w, height = leg_h, units = "cm",
             bg = "transparent", device = cairo_pdf, limitsize = FALSE)
      
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
        ggsave(out_png, plot = p, width = w_cm, height = h_cm, units = "cm", dpi = 600, bg = "white")
      }
      
      groups <- unique(processed_data()$Group)
      pal_named <- group_palette()
      lp <- legend_plot(groups, pal_named, horizontal = TRUE)
      
      legend_png_path <- file.path(tmpdir, "Legend.png")
      leg_w <- 22
      leg_h <- if (length(groups) > 3) 2.2 else 1.4
      
      ggsave(legend_png_path, lp,
             width = leg_w, height = leg_h, units = "cm",
             dpi = 600, bg = "transparent")
      
      zip::zipr(zipfile = file, files = list.files(tmpdir, full.names = TRUE))
    }
  )
}

shinyApp(ui = ui, server = server)