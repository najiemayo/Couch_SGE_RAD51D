figure4a_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(paste("Missing R packages:", paste(missing, collapse = ", ")), call. = FALSE)
  }
}

figure4a_check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

figure4a_header_key <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  gsub("[^a-z0-9]", "", tolower(trimws(x)))
}

figure4a_data <- function(workbook, sheet = NULL) {
  figure4a_require("readxl")
  figure4a_check(length(workbook) == 1L && file.exists(workbook),
    "Workbook not found.")
  sheets <- readxl::excel_sheets(workbook)
  if (!is.null(sheet)) {
    if (is.numeric(sheet)) {
      figure4a_check(length(sheet) == 1L && sheet %in% seq_along(sheets),
        "Invalid sheet index.")
      sheet <- sheets[[sheet]]
    }
    figure4a_check(length(sheet) == 1L && sheet %in% sheets, "Sheet not found.")
    sheets <- sheet
  }

  required <- c("variantclass", "aa")
  candidates <- list()
  for (current_sheet in sheets) {
    preview <- readxl::read_excel(workbook, sheet = current_sheet,
      col_names = FALSE, col_types = "text",
      range = readxl::cell_limits(c(1L, 1L), c(20L, NA)),
      .name_repair = "minimal", progress = FALSE)
    if (!nrow(preview)) next
    hits <- which(apply(as.matrix(preview), 1L, function(row) {
      keys <- figure4a_header_key(row)
      all(required %in% keys) &&
        any(grepl("^varcallfunctionalcategory", keys))
    }))
    for (header_row in hits) {
      candidates[[length(candidates) + 1L]] <- list(
        sheet = current_sheet, row = header_row)
    }
  }
  figure4a_check(length(candidates) == 1L,
    "A unique variant-table header is required; specify sheet if necessary.")

  selected <- candidates[[1L]]
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = selected$sheet,
    range = readxl::cell_limits(c(selected$row, 1L), c(NA, NA)),
    col_types = "text", .name_repair = "minimal", progress = FALSE),
    stringsAsFactors = FALSE, check.names = FALSE)
  keys <- figure4a_header_key(names(raw))
  figure4a_check(!anyDuplicated(keys), "Duplicate normalized column names.")
  column <- function(key) {
    index <- which(keys == key)
    figure4a_check(length(index) == 1L, paste("Missing or duplicate column:", key))
    raw[[index]]
  }
  category_index <- which(grepl("^varcallfunctionalcategory", keys))
  figure4a_check(length(category_index) == 1L,
    "Expected one VarCall functional-category column.")
  data <- data.frame(
    AA = trimws(as.character(column("aa"))),
    variant_class = trimws(as.character(column("variantclass"))),
    functional_category = trimws(as.character(raw[[category_index]])),
    stringsAsFactors = FALSE
  )
  data[data %in% c("", ".", "NA", "N/A")] <- NA
  event_key <- tolower(data$variant_class)
  event_type <- data$variant_class
  event_type[event_key %in% c("missense", "start codon missense",
    "nonsense_lost", "missense/start-codon-lost",
    "missense/stop-codon-lost")] <- "Missense"
  event_type[event_key %in% c("5' utr", "3' utr", "utr")] <- "UTR"
  supported_events <- c("Canonical splice", "Intronic", "Missense", "Silent",
    "Nonsense", "UTR")
  keep <- !is.na(event_type) & event_type %in% supported_events
  data <- data[keep, , drop = FALSE]
  figure4a_check(nrow(data) > 0L, "No supported variant classes found.")

  category_levels <- c("P strong", "P moderate", "P supporting", "VUS",
    "B supporting", "B moderate", "B strong")
  figure4a_check(all(!is.na(data$functional_category) &
    data$functional_category %in% category_levels),
    "Missing or unrecognized VarCall functional categories.")
  data$AApos <- suppressWarnings(as.numeric(
    sub("^[^0-9]*([0-9]+).*", "\\1", data$AA)))
  data$AltAA <- substr(data$AA, nchar(data$AA), nchar(data$AA))
  data$AltAA[data$AltAA == "r"] <- "stop"
  aa_levels <- c("A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N",
    "Q", "C", "G", "P", "R", "H", "K", "D", "E", "stop")
  data <- data[!is.na(data$AA) & data$AA != "." & data$AltAA %in% aa_levels &
    is.finite(data$AApos), , drop = FALSE]
  figure4a_check(nrow(data) > 0L, "No amino-acid substitutions found.")
  data$AltAA <- factor(data$AltAA, levels = aa_levels)
  data$functional_prediction <- factor(data$functional_category,
    levels = category_levels)

  heatmap <- data |>
    dplyr::arrange(functional_prediction) |>
    dplyr::group_by(AApos, AltAA) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(AApos, AltAA, functional_prediction)
  max_position <- max(heatmap$AApos)
  gaps <- setdiff(seq_len(max_position), unique(heatmap$AApos))
  if (length(gaps)) {
    gap_rows <- expand.grid(AApos = gaps, AltAA = aa_levels,
      stringsAsFactors = FALSE)
    gap_rows$AltAA <- factor(gap_rows$AltAA, levels = aa_levels)
    gap_rows$functional_prediction <- factor("gap",
      levels = c(category_levels, "gap"))
    heatmap$functional_prediction <- factor(heatmap$functional_prediction,
      levels = c(category_levels, "gap"))
    heatmap <- rbind(heatmap, gap_rows)
  }
  heatmap$AApos <- factor(heatmap$AApos,
    levels = as.character(seq_len(max_position)))
  heatmap$functional_prediction <- factor(heatmap$functional_prediction,
    levels = c(category_levels, "gap"))
  dplyr::arrange(heatmap, AApos)
}

make_figure4a <- function(workbook, output_dir = "Figure4A_output",
                          sheet = NULL, dpi = 600L,
                          font_family = "Calibri") {
  figure4a_require(c("readxl", "dplyr", "ggplot2", "ragg"))
  figure4a_check(length(dpi) == 1L && is.finite(dpi) && dpi >= 72 &&
    dpi == as.integer(dpi), "Invalid DPI.")
  figure4a_check(length(font_family) == 1L && !is.na(font_family),
    "Invalid font family.")

  heatmap <- figure4a_data(workbook, sheet)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figure4a_check(dir.exists(output_dir), "Output directory could not be created.")

  aa_colors <- c("#B22222", "#FFA500", "#FFD700", "gray", "#00FFFF",
    "#66ddff", "#00b4c5", "white")
  names(aa_colors) <- levels(heatmap$functional_prediction)
  desired_levels <- levels(heatmap$functional_prediction)
  desired_levels <- desired_levels[desired_levels != "gap"]
  observed_positions <- sort(unique(as.integer(as.character(heatmap$AApos))))
  heatmap_breaks <- observed_positions[
    observed_positions == 1L | observed_positions %% 50L == 0L |
      observed_positions == max(observed_positions)]

  plot <- ggplot2::ggplot(heatmap,
    ggplot2::aes(x = AApos, y = AltAA, fill = functional_prediction)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = aa_colors, breaks = desired_levels) +
    ggplot2::scale_x_discrete(position = "bottom",
      name = "Amino Acid Position", expand = c(0, 0), breaks = heatmap_breaks) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(y = "Amino Acid Substitution") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      axis.text = ggplot2::element_text(size = 15, face = "bold"),
      axis.title = ggplot2::element_text(size = 16.5, face = "bold"),
      legend.position = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 18, face = "bold"),
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 0.5,
        margin = ggplot2::margin(r = 4)),
      axis.text.x = ggplot2::element_text(vjust = 1),
      plot.margin = ggplot2::margin(t = 5.5, r = 14, b = 5.5, l = 5.5)
    )

  tiff_path <- file.path(output_dir, paste0("Figure4A_", dpi, "dpi.tiff"))
  csv_path <- file.path(output_dir, "Figure4A_heatmap.csv")
  ggplot2::ggsave(tiff_path, plot = plot, width = 16, height = 5,
    units = "in", dpi = dpi, device = ragg::agg_tiff,
    compression = "lzw", limitsize = FALSE)
  utils::write.csv(heatmap, csv_path, row.names = FALSE, na = "")
  invisible(list(files = c(tiff = tiff_path, heatmap = csv_path),
    data = heatmap))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  workbook <- if (length(args)) args[[1L]] else if (interactive()) {
    file.choose()
  } else {
    stop("Workbook path required.", call. = FALSE)
  }
  output_dir <- if (length(args) >= 2L) args[[2L]] else "Figure4A_output"
  dpi <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 600L
  sheet <- if (length(args) >= 4L) args[[4L]] else NULL
  make_figure4a(workbook, output_dir = output_dir, dpi = dpi, sheet = sheet)
}
