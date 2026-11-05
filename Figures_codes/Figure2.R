figure2_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(paste("Missing R packages:", paste(missing, collapse = ", ")), call. = FALSE)
  }
}

figure2_check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

figure2_header_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  gsub("[^a-z0-9]", "", tolower(trimws(x)))
}

figure2_read <- function(workbook, sheet = NULL) {
  figure2_require("readxl")
  figure2_check(length(workbook) == 1L && file.exists(workbook),
    "Workbook not found.")

  sheets <- readxl::excel_sheets(workbook)
  if (!is.null(sheet)) {
    if (is.numeric(sheet)) {
      figure2_check(length(sheet) == 1L && sheet %in% seq_along(sheets),
        "Invalid sheet index.")
      sheet <- sheets[[sheet]]
    }
    figure2_check(length(sheet) == 1L && sheet %in% sheets, "Sheet not found.")
    sheets <- sheet
  }

  required <- c("poshg38", "ref", "alt", "variantclass", "sge")
  candidates <- list()
  for (current_sheet in sheets) {
    preview <- readxl::read_excel(workbook, sheet = current_sheet,
      col_names = FALSE, col_types = "text",
      range = readxl::cell_limits(c(1L, 1L), c(20L, NA)),
      .name_repair = "minimal", progress = FALSE)
    if (!nrow(preview)) next
    hits <- which(apply(as.matrix(preview), 1L, function(row) {
      keys <- figure2_header_key(row)
      all(required %in% keys) &&
        any(grepl("^varcallfunctionalscore", keys)) &&
        any(grepl("^varcallfunctionalcategory", keys))
    }))
    for (header_row in hits) {
      candidates[[length(candidates) + 1L]] <- list(
        sheet = current_sheet, row = header_row)
    }
  }
  figure2_check(length(candidates) == 1L,
    "A unique variant-table header is required; specify sheet if necessary.")

  selected <- candidates[[1L]]
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = selected$sheet,
    range = readxl::cell_limits(c(selected$row, 1L), c(NA, NA)),
    col_types = "text", .name_repair = "minimal", progress = FALSE),
    stringsAsFactors = FALSE, check.names = FALSE)
  keys <- figure2_header_key(names(raw))
  figure2_check(!anyDuplicated(keys), "Duplicate normalized column names.")
  column <- function(key) {
    index <- which(keys == key)
    figure2_check(length(index) == 1L, paste("Missing or duplicate column:", key))
    raw[[index]]
  }
  starts_with <- function(prefix) {
    index <- which(startsWith(keys, prefix))
    figure2_check(length(index) == 1L,
      paste("Expected one column matching:", prefix))
    raw[[index]]
  }

  data <- data.frame(
    POS_hg38 = trimws(as.character(column("poshg38"))),
    REF = trimws(as.character(column("ref"))),
    ALT = trimws(as.character(column("alt"))),
    variant_class = trimws(as.character(column("variantclass"))),
    SGE = trimws(as.character(column("sge"))),
    functional_score = trimws(as.character(starts_with("varcallfunctionalscore"))),
    functional_category = trimws(as.character(starts_with("varcallfunctionalcategory"))),
    stringsAsFactors = FALSE
  )
  data[data %in% c("", ".", "NA", "N/A")] <- NA
  data <- data[!is.na(data$POS_hg38) & !is.na(data$REF) & !is.na(data$ALT), ,
    drop = FALSE]
  figure2_check(nrow(data) > 0L, "No variant rows found in the selected sheet.")
  data$POS_hg38 <- suppressWarnings(as.numeric(data$POS_hg38))
  data$functional_score <- suppressWarnings(as.numeric(data$functional_score))
  figure2_check(all(is.finite(data$POS_hg38) & data$POS_hg38 > 0),
    "Invalid genomic positions.")
  figure2_check(!anyNA(data$functional_category),
    "Missing functional category in variant data.")
  valid_categories <- c("P strong", "P moderate", "P supporting", "VUS",
    "B supporting", "B moderate", "B strong")
  figure2_check(all(data$functional_category %in% valid_categories),
    "Missing or unrecognized VarCall functional categories.")
  data
}

figure2_data <- function(workbook, sheet = NULL) {
  raw <- figure2_read(workbook, sheet)
  raw$uPOS <- paste(raw$POS_hg38, raw$REF, raw$ALT, sep = "_")
  raw$GRCh38Location <- raw$POS_hg38
  raw$EventType <- raw$variant_class
  event_key <- tolower(raw$EventType)
  raw$EventType[event_key %in% c("missense", "start codon missense",
    "nonsense_lost", "missense/start-codon-lost",
    "missense/stop-codon-lost")] <- "Missense"
  raw$EventType[event_key %in% c("5' utr", "3' utr", "utr")] <- "UTR"

  event_levels <- c("Canonical splice", "Intronic", "Missense", "Silent",
    "Nonsense", "UTR")
  event_display_order <- c("Silent", "Intronic", "Missense",
    "Canonical splice", "Nonsense", "UTR")
  variants <- raw[raw$EventType %in% event_levels, , drop = FALSE]
  figure2_check(nrow(variants) > 0L, "No supported variant classes found.")
  variants$EventType <- factor(variants$EventType, levels = event_levels)
  variants$Eta <- variants$functional_score
  variants$Exon <- factor(variants$SGE, levels = paste0("E", 1:10))
  category_levels <- c("P strong", "P moderate", "P supporting", "VUS",
    "B supporting", "B moderate", "B strong")
  variants$functional_prediction <- factor(variants$functional_category,
    levels = category_levels)

  figure2_check(all(paste0("E", 1:10) %in% as.character(variants$Exon)),
    "Expected assay exons E1-E10.")
  figure2_check(any(is.finite(variants$Eta)),
    "No numeric functional scores found.")
  event_counts <- as.data.frame(dplyr::count(variants, EventType,
    functional_prediction, name = "n"))
  event_counts <- dplyr::group_by(event_counts, EventType)
  event_counts <- dplyr::mutate(event_counts,
    percent = n / sum(n), total = sum(n),
    category_label = paste0(EventType, ", n=", total))
  event_counts <- dplyr::ungroup(event_counts)
  event_totals <- event_counts[!duplicated(event_counts$EventType),
    c("EventType", "total")]
  present_events <- event_display_order[
    event_display_order %in% as.character(event_totals$EventType)]
  label_levels <- paste0(present_events, ", n=",
    event_totals$total[match(present_events, as.character(event_totals$EventType))])
  event_counts$category_label <- factor(event_counts$category_label,
    levels = label_levels)

  exon_counts <- as.data.frame(dplyr::count(variants, Exon, EventType,
    functional_prediction, name = "n"))
  exon_counts <- dplyr::group_by(exon_counts, Exon, EventType)
  exon_counts <- dplyr::mutate(exon_counts, percent = n / sum(n))
  exon_counts <- dplyr::ungroup(exon_counts)
  exon_counts$EventTypeLabel <- gsub("Canonical splice", "Canonical\nsplice",
    as.character(exon_counts$EventType), fixed = TRUE)
  facet_levels <- rev(gsub("Canonical splice", "Canonical\nsplice",
    event_display_order, fixed = TRUE))
  exon_counts$EventTypeLabel <- factor(exon_counts$EventTypeLabel,
    levels = facet_levels)

  list(
    variants = variants,
    event_counts = event_counts,
    exon_counts = exon_counts
  )
}

figure2_export_plot <- function(plot, output_dir, stem, width, height, dpi) {
  tiff_path <- file.path(output_dir, paste0(stem, "_", dpi, "dpi.tiff"))
  ggplot2::ggsave(tiff_path, plot = plot, width = width, height = height,
    units = "in", dpi = dpi, device = ragg::agg_tiff,
    compression = "lzw", limitsize = FALSE, bg = "white")
  tiff_path
}

make_figure2 <- function(workbook, output_dir = "Figure2_output",
                         sheet = NULL, dpi = 600L, font_family = "Arial") {
  figure2_require(c("readxl", "ggplot2", "dplyr",
    "ggh4x", "ggbreak", "ragg"))
  figure2_check(length(dpi) == 1L && is.finite(dpi) && dpi >= 72 &&
    dpi == as.integer(dpi), "Invalid DPI.")
  figure2_check(length(font_family) == 1L && !is.na(font_family),
    "Invalid font family.")
  result <- figure2_data(workbook, sheet)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figure2_check(dir.exists(output_dir),
    "Output directory could not be created.")

  event_colors <- c("#a3a500", "#0059a3", "#4f9c9c", "#33BBEE",
    "#db6d00", "#800A01")
  names(event_colors) <- levels(result$variants$EventType)
  category_colors <- c("#800000", "#FFA500", "#FFD700", "gray",
    "#00b4c5", "#0073e6", "#0059a3")
  names(category_colors) <- levels(result$variants$functional_prediction)

  plot_2a <- ggplot2::ggplot(result$variants,
    ggplot2::aes(x = Eta, fill = EventType)) +
    ggplot2::geom_histogram(binwidth = 0.025, position = "stack") +
    ggplot2::scale_fill_manual(values = event_colors) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "cm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(.35, .7),
      legend.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 18, face = "bold"),
      axis.text = ggplot2::element_text(size = 16, face = "bold"),
      legend.text = ggplot2::element_text(size = 16, face = "bold")
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Model based functional score", y = "Number of variants")

  exon_positions <- result$variants |>
    dplyr::group_by(Exon) |>
    dplyr::summarise(start = min(GRCh38Location), end = max(GRCh38Location),
      .groups = "drop") |>
    dplyr::arrange(start)
  figure2_check(nrow(exon_positions) == 10L,
    "Expected one coordinate range for each assay exon E1-E10.")
  figure2_check(all(head(exon_positions$end, -1L) <
    tail(exon_positions$start, -1L)),
    "Assay exon genomic ranges overlap or are out of order.")
  gap <- 20
  break_positions <- as.vector(rbind(exon_positions$end[1:9],
    exon_positions$start[2:10]))
  annotation_positions <- c(exon_positions$start[1] + gap,
    as.vector(rbind(exon_positions$end[1:9] - gap,
      exon_positions$start[2:10] + gap)),
    exon_positions$end[10] - gap)
  annotation_labels <- c(exon_positions$start[1],
    break_positions, exon_positions$end[10])
  exon_label_positions <- exon_positions$end
  exon_labels <- as.character(exon_positions$Exon)
  plot_2b <- ggplot2::ggplot(result$variants |>
    dplyr::arrange(GRCh38Location),
    ggplot2::aes(x = GRCh38Location, y = Eta, colour = EventType)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_manual(values = event_colors) +
    ggplot2::scale_x_reverse() +
    ggplot2::annotate("text", x = annotation_positions, y = 0,
      label = annotation_labels, angle = 90, hjust = 0, size = 8) +
    ggplot2::annotate("text", x = exon_label_positions, y = 1.5,
      label = exon_labels, hjust = 0, size = 8) +
    ggplot2::coord_cartesian(ylim = c(0, 1.5), clip = "off") +
    ggplot2::labs(x = "GRCh38 Genomic Coordinate", y = "Functional score") +
    ggplot2::theme_classic(base_family = font_family) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 22.4, face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(size = 24),
      axis.title.x = ggplot2::element_text(size = 32),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 18, face = "bold"),
      legend.margin = ggplot2::margin(0, 0, 0, -11),
      legend.box.margin = ggplot2::margin(0, 0, 0, -10),
      panel.border = ggplot2::element_rect(colour = "black",
        linewidth = 0.5, fill = NA)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 5))) +
    ggbreak::scale_x_break(breaks = break_positions, ticklabels = NULL,
      expand = c(0.01, 0.01))

  plot_2c <- ggplot2::ggplot(result$event_counts,
    ggplot2::aes(x = category_label, y = percent * 100,
      fill = functional_prediction)) +
    ggplot2::geom_bar(stat = "identity", width = 0.8) +
    ggplot2::scale_fill_manual(values = category_colors, drop = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Percentage %") +
    ggplot2::theme_classic(base_family = font_family) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 18, face = "bold"),
      axis.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 18, face = "bold")
    )

  plot_2d <- ggplot2::ggplot(result$exon_counts,
    ggplot2::aes(x = Exon, y = percent * 100,
      fill = functional_prediction)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5) +
    ggplot2::scale_fill_manual(values = category_colors, drop = FALSE) +
    ggplot2::facet_wrap(~EventTypeLabel, ncol = 1, strip.position = "right",
      scales = "free_x") +
    ggplot2::labs(x = "Exon", y = "Percentage %") +
    ggplot2::theme_classic(base_family = font_family) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(size = 15, angle = 0),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 18, face = "bold")
    ) +
    ggh4x::facetted_pos_scales(x = rep(list(ggplot2::scale_x_discrete(
      limits = paste0("E", 1:10))), 6L))

  plots <- list(
    Figure2A_eta_hist = list(plot = plot_2a, width = 5, height = 5),
    Figure2B_eta_byexon_POS = list(plot = plot_2b, width = 20, height = 5),
    Figure2C_barplot_percentage_all = list(plot = plot_2c, width = 9, height = 5),
    Figure2D_barplot_percentage_byexon = list(plot = plot_2d, width = 10, height = 5)
  )
  files <- unlist(lapply(names(plots), function(stem) {
    spec <- plots[[stem]]
    figure2_export_plot(spec$plot, output_dir, stem,
      spec$width, spec$height, dpi)
  }), use.names = TRUE)
  csv_files <- c(
    variants = file.path(output_dir, "Figure2_variants.csv"),
    event_counts = file.path(output_dir, "Figure2C_counts.csv"),
    exon_counts = file.path(output_dir, "Figure2D_counts.csv")
  )
  for (name in names(csv_files)) {
    utils::write.csv(result[[name]], csv_files[[name]],
      row.names = FALSE, na = "")
  }
  invisible(list(files = c(files, csv_files), data = result))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  workbook <- if (length(args)) args[[1L]] else if (interactive()) {
    file.choose()
  } else {
    stop("Workbook path required.", call. = FALSE)
  }
  output_dir <- if (length(args) >= 2L) args[[2L]] else "Figure2_output"
  dpi <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 600L
  sheet <- if (length(args) >= 4L) args[[4L]] else NULL
  make_figure2(workbook, output_dir = output_dir, dpi = dpi, sheet = sheet)
}
