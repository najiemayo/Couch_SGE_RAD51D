#!/usr/bin/env Rscript
# RAD51D Supplementary Figure S3A -- refined, no-leader-line R implementation.
# Version 1.0.0
#
# The supplied refined TIFF was generated in Python. This is an R translation
# of its workbook-processing rules and layout, not that original renderer.
# Fonts and rendering engines can cause small visual differences.
#
# Input: the unmodified Supplementary Data 2 XLSX workbook. Images and archived
# CSVs are NOT plotting inputs. CSV references are used only for validation.
# A red marker means at least one assayed missense SNV at the residue has the
# explicit category "P strong". Black means no such hit, NOT benignity.
# No VarCall refitting, posterior threshold, SpliceAI filter or new structural
# nomination is performed. Annotations are inherited from workbook columns U/V.
#
# R console:
#   source("scripts/panel_a.R")
#   make_panel_a("data/Supplementary_Data_2_updatedVarCall092826.xlsx",
#                reference_dir = "validation/reference")
# Terminal:
#   Rscript scripts/panel_a.R --workbook data/WORKBOOK.xlsx --outdir output
#
# Required: readxl, ragg, systemfonts, digest; svglite for SVG export.
# Installation is explicit; this script does not install packages.
# Documentation: https://readxl.tidyverse.org/reference/read_excel.html
#                https://ragg.r-lib.org/reference/agg_tiff.html
#                https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grobWidth.html

PANEL_A_VERSION <- "1.0.0"
PANEL_A_CATEGORIES <- c("P strong", "P moderate", "P supporting", "VUS",
                         "B supporting", "B moderate", "B strong")
PANEL_A_ROWS <- data.frame(
  row_index = seq_len(6L),
  source = c(rep("Greenhough et al. (2023)", 4L), rep("Rawal et al. (2023)", 2L)),
  structural_context = c("Interaction with ATP", "Interaction with XRCC2",
    "Interaction with RAD51C", "Monomer stability", "ATP binding/hydrolysis", "ssDNA binding"),
  marker_shape = c("v", "^", "s", "o", "D", "v"), stringsAsFactors = FALSE
)

panel_a_assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
panel_a_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages. Run: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}
panel_a_clean <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""; trimws(x)
}
panel_a_missing_annotation <- function(x) x %in% c("", "NA", "N/A", "None", ".")
panel_a_join <- function(x) paste(x, collapse = "; ")
panel_a_write_csv <- function(x, file) {
  old <- options(digits = 17); on.exit(options(old), add = TRUE)
  utils::write.csv(x, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

# Read from A1 to retain actual Excel row numbers, including blank rows.
read_panel_a_workbook <- function(workbook, sheet = NULL) {
  panel_a_require("readxl")
  panel_a_assert(file.exists(workbook), paste("Workbook not found:", workbook))
  sheets <- readxl::excel_sheets(workbook)
  if (is.null(sheet)) {
    if (length(sheets) != 1L) stop("Multiple sheets; supply sheet explicitly.", call. = FALSE)
    sheet <- sheets[[1L]]
  }
  panel_a_assert(is.character(sheet) && length(sheet) == 1L && sheet %in% sheets,
                 "sheet must be an existing worksheet name.")
  raw <- as.data.frame(readxl::read_excel(
    workbook, sheet = sheet, range = readxl::cell_limits(c(1, 1), c(NA, NA)),
    col_types = "text", trim_ws = TRUE, na = "", .name_repair = "minimal"
  ), stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- trimws(names(raw))
  required <- c("Annotation_hg38", "c.", "p.", "AA", "variant class", "SGE",
    "Greenhouph Nature 2023", "Rawal Nature 2023",
    "VarCall Posterior probability of pathogenicity",
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
  for (name in required) panel_a_assert(sum(names(raw) == name) == 1L,
                                        paste("Missing/duplicate workbook column:", name))
  category <- grep("^VarCall functional category", names(raw), value = TRUE)
  panel_a_assert(length(category) == 1L, "Expected one explicit VarCall functional-category column.")
  raw$workbook_row <- seq_len(nrow(raw)) + 1L
  list(data = raw, sheet = sheet, category_column = category)
}

extract_panel_a_snvs <- function(raw, category_column) {
  m <- raw[tolower(panel_a_clean(raw[["variant class"]])) == "missense", , drop = FALSE]
  panel_a_assert(nrow(m) > 0L, "No missense-labeled SNVs found.")
  get <- function(name) panel_a_clean(m[[name]])
  for (key in c("Annotation_hg38", "c.")) {
    z <- get(key)
    panel_a_assert(all(nzchar(z)), paste("Missing SNV identifier in", key))
    panel_a_assert(!anyDuplicated(z), paste("Duplicate SNV identifiers in", key,
                                           "-- resolve rather than count twice."))
  }
  aa <- get("AA"); cdna <- get("c.")
  pattern <- "^([A-Z]|Ter)([0-9]+)([A-Z?]|Ter)$"
  panel_a_assert(all(grepl(pattern, aa)), "Unrecognized protein annotation in AA column.")
  ref <- sub(pattern, "\\1", aa)
  position <- as.integer(sub(pattern, "\\2", aa))
  panel_a_assert(all(position >= 1L & position <= 329L), "Expected RAD51D positions 1-329.")
  cp <- "^c\\.([0-9]+)[ACGT]>[ACGT]$"
  panel_a_assert(all(grepl(cp, cdna)), "Unrecognized coding-SNV notation.")
  cpos <- as.integer(sub(cp, "\\1", cdna))
  panel_a_assert(all((cpos + 2L) %/% 3L == position), "cDNA and protein positions disagree.")
  for (p in unique(position)) panel_a_assert(length(unique(ref[position == p])) == 1L,
                                            paste("Conflicting reference amino acids at", p))
  category <- get(category_column)
  panel_a_assert(all(category %in% PANEL_A_CATEGORIES),
                 "Missing/unrecognized explicit VarCall category; do not infer a replacement.")
  posterior <- suppressWarnings(as.numeric(get("VarCall Posterior probability of pathogenicity")))
  panel_a_assert(all(is.finite(posterior) & posterior >= 0 & posterior <= 1),
                 "Missing/invalid posterior value; inspect the workbook. Used only for audit.")
  splice_names <- c("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
  splice <- vapply(splice_names, function(name) {
    text <- get(name); missing <- panel_a_missing_annotation(text)
    numeric <- suppressWarnings(as.numeric(text))
    panel_a_assert(all(missing | (is.finite(numeric) & numeric >= 0 & numeric <= 1)),
                   paste("Invalid SpliceAI value in", name))
    numeric[missing] <- NA_real_; numeric
  }, numeric(nrow(m)))
  splice <- matrix(splice, nrow = nrow(m), ncol = length(splice_names))
  splice_max <- apply(splice, 1L, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE))
  # Distinct SNVs are retained even when their protein changes coincide.
  data.frame(
    workbook_row = m$workbook_row, variant_hg38 = get("Annotation_hg38"),
    cDNA_change = cdna, protein_HGVS = get("p."), protein_change = aa,
    reference_AA = ref, position = position, variant_class = get("variant class"),
    exon_target = get("SGE"), VarCall_functional_category = category,
    VarCall_posterior = posterior, is_P_strong = as.integer(category == "P strong"),
    maximum_SpliceAI_delta = splice_max,
    Greenhough_workbook_annotation = get("Greenhouph Nature 2023"),
    Rawal_workbook_annotation = get("Rawal Nature 2023"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

# A residue's structural annotations are unioned across its assayed SNVs.
# Its P_strong flag considers ALL assayed missense SNVs at that residue, not
# only the particular SNVs whose annotation cells are populated.
build_panel_a_tables <- function(snvs) {
  plotted <- list(); annotated <- list()
  for (position in sort(unique(snvs$position))) {
    z <- snvs[snvs$position == position, , drop = FALSE]
    g <- sort(unique(z$Greenhough_workbook_annotation[
      !panel_a_missing_annotation(z$Greenhough_workbook_annotation)]), method = "radix")
    r <- sort(unique(z$Rawal_workbook_annotation[
      !panel_a_missing_annotation(z$Rawal_workbook_annotation)]), method = "radix")
    if (!length(g) && !length(r)) next
    hit <- z[z$is_P_strong == 1L, , drop = FALSE]
    annot_rows <- z$workbook_row[!panel_a_missing_annotation(z$Greenhough_workbook_annotation) |
                                 !panel_a_missing_annotation(z$Rawal_workbook_annotation)]
    base <- data.frame(
      residue = paste0(z$reference_AA[[1L]], position), position = position,
      reference_AA = z$reference_AA[[1L]], n_assayed_missense_SNVs = nrow(z),
      n_P_strong_SNVs = nrow(hit), P_strong_hit = as.integer(nrow(hit) > 0L),
      point_color = if (nrow(hit)) "#FF0000" else "#000000",
      P_strong_cDNA_changes = panel_a_join(hit$cDNA_change),
      P_strong_protein_changes = panel_a_join(hit$protein_change),
      maximum_missense_posterior = max(z$VarCall_posterior),
      Greenhough_workbook_annotation = panel_a_join(g),
      Rawal_workbook_annotation = panel_a_join(r),
      annotation_workbook_rows = panel_a_join(annot_rows),
      P_strong_workbook_rows = panel_a_join(hit$workbook_row),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    row_ids <- integer()
    # These rules deliberately match the executed Python source exactly.
    # ATP matches both normal wording and I311's "RAD51D (ATP, Non-bonded)".
    # XRCC2 in an ATP-ligand annotation alone does not imply a protein contact.
    if (any(grepl("ATP", g, fixed = TRUE))) row_ids <- c(row_ids, 1L)
    if (any(startsWith(g, "XRCC2 ("))) row_ids <- c(row_ids, 2L)
    if (any(grepl("RAD51C (", g, fixed = TRUE))) row_ids <- c(row_ids, 3L)
    if (any(grepl("Predicted monomer instability", g, fixed = TRUE))) row_ids <- c(row_ids, 4L)
    if (any(r %in% c("ATP binding", "Mg ion coordination", "ATPase site"))) row_ids <- c(row_ids, 5L)
    if ("ssDNA binding" %in% r) row_ids <- c(row_ids, 6L)
    for (i in row_ids) plotted[[length(plotted) + 1L]] <- cbind(PANEL_A_ROWS[i, , drop = FALSE], base)
    base$in_displayed_six_contexts <- as.integer(length(row_ids) > 0L)
    base$displayed_contexts <- panel_a_join(PANEL_A_ROWS$structural_context[row_ids])
    annotated[[length(annotated) + 1L]] <- base
  }
  panel_a_assert(length(plotted) > 0L, "No residues fall in the six displayed structural contexts.")
  plotted <- do.call(rbind, plotted)
  plotted <- plotted[order(plotted$row_index, plotted$position), , drop = FALSE]
  rownames(plotted) <- NULL
  annotated <- do.call(rbind, annotated); rownames(annotated) <- NULL
  omitted <- annotated[annotated$in_displayed_six_contexts == 0L, , drop = FALSE]
  relevant <- snvs[snvs$position %in% plotted$position, , drop = FALSE]
  counts <- do.call(rbind, lapply(seq_len(6L), function(i) {
    z <- plotted[plotted$row_index == i, , drop = FALSE]
    data.frame(source = PANEL_A_ROWS$source[i], context = PANEL_A_ROWS$structural_context[i],
      n_positions = nrow(z), n_P_strong_hit_positions = sum(z$P_strong_hit),
      P_strong_hit_residues = panel_a_join(z$residue[z$P_strong_hit == 1L]),
      stringsAsFactors = FALSE)
  }))
  panel_a_assert(!anyDuplicated(paste(plotted$row_index, plotted$position)),
                 "Duplicate residue/context markers detected.")
  list(snvs = snvs, plotted = plotted, annotated = annotated, omitted = omitted,
       relevant_snvs = relevant, counts = counts)
}

analyse_panel_a <- function(workbook, sheet = NULL) {
  input <- read_panel_a_workbook(workbook, sheet)
  out <- build_panel_a_tables(extract_panel_a_snvs(input$data, input$category_column))
  out$source_sheet <- input$sheet; out$category_column <- input$category_column
  out
}

panel_a_table_files <- c(
  plotted = "PanelA_plotted_residue_contexts.csv",
  annotated = "PanelA_all_structurally_annotated_positions.csv",
  relevant_snvs = "PanelA_SNV_audit.csv",
  omitted = "PanelA_structural_contexts_not_displayed.csv",
  counts = "PanelA_context_counts.csv"
)

# Frozen CSVs are tests, not an alternative source for plotting.
check_panel_a_reference <- function(result, reference_dir, tolerance = 1e-12) {
  reports <- character()
  for (key in names(panel_a_table_files)) {
    path <- file.path(reference_dir, panel_a_table_files[[key]])
    panel_a_assert(file.exists(path), paste("Missing reference file:", path))
    ref <- utils::read.csv(path, colClasses = "character", check.names = FALSE,
      na.strings = "", fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
    got <- result[[key]]
    panel_a_assert(identical(names(got), names(ref)), paste("Reference column mismatch:", key))
    panel_a_assert(nrow(got) == nrow(ref), paste("Reference row-count mismatch:", key))
    for (name in names(got)) {
      if (is.numeric(got[[name]])) {
        expected <- suppressWarnings(as.numeric(ref[[name]])); actual <- got[[name]]
        parsed <- is.na(ref[[name]]) | !nzchar(ref[[name]]) | !is.na(expected)
        panel_a_assert(all(parsed), paste("Non-numeric reference:", key, name))
        ok <- (is.na(actual) & is.na(expected)) |
          (!is.na(actual) & !is.na(expected) & abs(actual - expected) <= tolerance * pmax(1, abs(expected)))
        ok[is.na(ok)] <- FALSE
      } else {
        ok <- panel_a_clean(got[[name]]) == panel_a_clean(ref[[name]])
      }
      panel_a_assert(all(ok), paste("Reference mismatch in", key, name,
        "at rows", paste(head(which(!ok), 8L), collapse = ", ")))
    }
    reports <- c(reports, sprintf("PASS %s: %d rows, %d columns", basename(path), nrow(got), ncol(got)))
  }
  reports
}

panel_a_style <- function(residue_size_pt = 9.5, row_size_pt = 11.5,
                          tick_size_pt = 11.5, title_size_pt = 13,
                          legend_size_pt = 10, dpi = 600) {
  sizes <- c(residue_size_pt, row_size_pt, tick_size_pt, title_size_pt, legend_size_pt)
  panel_a_assert(all(is.finite(sizes) & sizes > 0), "Font sizes must be positive.")
  panel_a_assert(length(dpi) == 1L && is.finite(dpi) && dpi > 0 && dpi == as.integer(dpi),
                 "dpi must be a positive integer.")
  list(width_in = 12, height_in = 6.6, dpi = dpi, residue_size_pt = residue_size_pt,
    row_size_pt = row_size_pt, tick_size_pt = tick_size_pt, title_size_pt = title_size_pt,
    legend_size_pt = legend_size_pt, left = 2.50, right = 11.80, xmax = 330,
    row_y = c(5.46, 4.74, 4.02, 3.30, 1.70, 1.085),
    label_offsets = c(.130, -.130, .275, -.275, .420, -.420),
    label_gap_in = .030, dense_marker_offset_in = .030,
    marker_size_pt = 3.5, red = "#FF0000", black = "#000000")
}

choose_panel_a_font <- function(font_family = NULL) {
  panel_a_require("systemfonts")
  available <- unique(systemfonts::system_fonts()$family)
  if (!is.null(font_family)) {
    panel_a_assert(length(font_family) == 1L && font_family %in% available,
                   paste("Requested font family is not installed:", font_family))
    return(font_family)
  }
  preferred <- c("Liberation Sans", "Arial", "DejaVu Sans")
  chosen <- preferred[preferred %in% available]
  panel_a_assert(length(chosen) > 0L,
                 "No preferred sans-serif font found; supply an installed font_family.")
  chosen[[1L]]
}

# Draw a mathematical marker rather than a font glyph. Dimensions are in points.
# R lwd=1 is nominally 1/96 inch, so convert Python's point widths explicitly.
panel_a_marker <- function(x, y, shape, color, size_pt = 3.5) {
  gp <- grid::gpar(fill = color, col = "white", lwd = .20 * 96 / 72, linejoin = "mitre")
  if (shape == "o") {
    grid::grid.circle(x = grid::unit(x, "native"), y = grid::unit(y, "native"),
                       r = grid::unit(size_pt / 2, "bigpts"), gp = gp)
  } else {
    xy <- switch(shape,
      "v" = list(x = c(-.5, .5, 0), y = c(.5, .5, -.5)),
      "^" = list(x = c(-.5, .5, 0), y = c(-.5, -.5, .5)),
      "s" = list(x = c(-.5, .5, .5, -.5), y = c(-.5, -.5, .5, .5)),
      "D" = list(x = c(0, 1, 0, -1) / sqrt(2), y = c(1, 0, -1, 0) / sqrt(2)),
      stop("Unrecognized marker shape.", call. = FALSE))
    grid::grid.polygon(x = grid::unit(x, "native") + grid::unit(xy$x * size_pt, "bigpts"),
                       y = grid::unit(y, "native") + grid::unit(xy$y * size_pt, "bigpts"), gp = gp)
  }
}

# Measurements are made on the active output device using its actual font.
# Marker x coordinates NEVER move. Label x coordinates match the markers;
# labels are staggered in y using the same six tiers as the refined Python panel.
# The small +/-0.030-inch marker y offset only separates adjacent categorical
# markers; it is not a numeric measurement. No leader lines are drawn.
draw_panel_a <- function(result, font_family, style = panel_a_style()) {
  s <- style; width <- s$width_in; height <- s$height_in
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, width), yscale = c(0, height), clip = "off"))
  on.exit(grid::popViewport(), add = TRUE)
  x <- function(p) s$left + p / s$xmax * (s$right - s$left)
  text_gp <- function(size) grid::gpar(fontfamily = font_family, fontface = "bold", fontsize = size, col = s$black)
  metrics <- function(label, size) {
    grob <- grid::textGrob(label, gp = text_gp(size))
    c(width = grid::convertWidth(grid::grobWidth(grob), "inches", valueOnly = TRUE),
      height = grid::convertHeight(grid::grobHeight(grob), "inches", valueOnly = TRUE))
  }
  text_bounds <- list(); label_layout <- list()
  text <- function(xx, yy, label, size, ha = "centre", kind = "annotation") {
    met <- metrics(label, size); w <- unname(met[[1L]]); h <- unname(met[[2L]])
    just_x <- switch(ha, left = 0, centre = .5, right = 1)
    grid::grid.text(label, x = grid::unit(xx, "native"), y = grid::unit(yy, "native"),
                     just = c(ha, "centre"), gp = text_gp(size))
    text_bounds[[length(text_bounds) + 1L]] <<- data.frame(
      text = label, kind = kind, x0 = xx - just_x * w, x1 = xx + (1 - just_x) * w,
      y0 = yy - h / 2, y1 = yy + h / 2, stringsAsFactors = FALSE)
  }
  line <- function(x0, y0, x1, y1, width_pt) grid::grid.segments(
    x0 = grid::unit(x0, "native"), x1 = grid::unit(x1, "native"),
    y0 = grid::unit(y0, "native"), y1 = grid::unit(y1, "native"),
    gp = grid::gpar(col = s$black, lwd = width_pt * 96 / 72, lineend = "butt"))

  text(.09, 6.29, "A", 25, ha = "left")
  panel_a_marker(s$left + .02, 6.36, "o", s$red, 4.1)
  text(s$left + .14, 6.36, "At least one P_strong missense SNV", s$legend_size_pt, ha = "left")
  panel_a_marker(s$left + 5.13, 6.36, "o", s$black, 4.1)
  text(s$left + 5.25, 6.36, "No P_strong missense SNV", s$legend_size_pt, ha = "left")
  text(s$right, 6.10, "Greenhough et al. (2023)", s$legend_size_pt, ha = "right")
  text(s$right, 2.19, "Rawal et al. (2023)", s$legend_size_pt, ha = "right")
  for (bounds in list(c(2.95, 5.98), c(.715, 2.08))) {
    grid::grid.rect(x = grid::unit(s$left, "native"), y = grid::unit(bounds[1], "native"),
      width = grid::unit(s$right - s$left, "native"), height = grid::unit(diff(bounds), "native"),
      just = c("left", "bottom"), gp = grid::gpar(fill = NA, col = s$black, lwd = .9 * 96 / 72))
  }
  for (i in seq_len(6L)) {
    yy <- s$row_y[i]; title <- PANEL_A_ROWS$structural_context[i]
    text(s$left - .095, yy, title, s$row_size_pt, ha = "right", kind = "context")
    items <- result$plotted[result$plotted$row_index == i, , drop = FALSE]
    items <- items[order(items$position), , drop = FALSE]
    lanes <- rep(list(list()), length(s$label_offsets))
    for (j in seq_len(nrow(items))) {
      item <- items[j, , drop = FALSE]; xx <- x(item$position)
      w <- metrics(item$residue, s$residue_size_pt)[[1L]]
      span <- c(xx - w / 2, xx + w / 2); chosen <- NA_integer_
      for (level in seq_along(lanes)) {
        free <- all(vapply(lanes[[level]], function(old) {
          span[2] + s$label_gap_in <= old[1] || old[2] + s$label_gap_in <= span[1]
        }, logical(1)))
        if (free) { chosen <- level; lanes[[level]][[length(lanes[[level]]) + 1L]] <- span; break }
      }
      panel_a_assert(!is.na(chosen), paste("No collision-free label tier for", title, item$residue,
        "-- reduce residue_size_pt or review the layout."))
      label_y <- yy + s$label_offsets[chosen]
      side <- if (j %% 2L == 1L) 1 else -1
      nearby <- (j > 1L && item$position - items$position[j - 1L] <= 1L) ||
        (j < nrow(items) && items$position[j + 1L] - item$position <= 1L)
      point_y <- yy + if (nearby) side * s$dense_marker_offset_in else 0
      panel_a_marker(xx, point_y, item$marker_shape, item$point_color, s$marker_size_pt)
      text(xx, label_y, item$residue, s$residue_size_pt, kind = "residue")
      label_layout[[length(label_layout) + 1L]] <- data.frame(
        row_index = i, context = title, residue = item$residue, position = item$position,
        point_x_inches = xx, point_y_inches = point_y, label_x_inches = xx,
        label_y_inches = label_y, label_font_size_pt = s$residue_size_pt,
        label_tier = chosen, leader_line = FALSE, color = item$point_color,
        marker = item$marker_shape, stringsAsFactors = FALSE)
    }
  }
  for (a in list(c(2.805, 2.605, 2.38), c(.59, .39, .155))) {
    line(s$left, a[1], s$right, a[1], 1.3)
    for (p in seq(0, 330, by = 10)) {
      major <- p %% 50 == 0
      line(x(p), a[1], x(p), a[1] - if (major) .070 else .035, if (major) 1.15 else .8)
      if (major) text(x(p), a[2], as.character(p), s$tick_size_pt, kind = "tick")
    }
    text((s$left + s$right) / 2, a[3], "AA position", s$title_size_pt, kind = "axis_title")
  }
  boxes <- do.call(rbind, text_bounds)
  clipped <- boxes$text[boxes$x0 < 0 | boxes$x1 > width | boxes$y0 < 0 | boxes$y1 > height]
  overlaps <- character()
  if (nrow(boxes) > 1L) for (i in seq_len(nrow(boxes) - 1L)) {
    for (j in seq.int(i + 1L, nrow(boxes))) {
      a <- boxes[i, ]; b <- boxes[j, ]
      if (min(a$x1, b$x1) > max(a$x0, b$x0) && min(a$y1, b$y1) > max(a$y0, b$y0))
        overlaps <- c(overlaps, paste(a$text, "<->", b$text))
    }
  }
  panel_a_assert(!length(clipped) && !length(overlaps), paste(
    "Text layout did not pass. Clipped:", paste(clipped, collapse = ", "),
    "; overlaps:", paste(overlaps, collapse = "; "),
    ". Reduce font sizes or select a metrically compatible sans-serif font."))
  layout <- do.call(rbind, label_layout)
  panel_a_assert(all(abs((layout$point_x_inches - s$left) / (s$right - s$left) * s$xmax - layout$position) < 1e-9),
                 "Marker x coordinate changed.")
  panel_a_assert(all(layout$label_x_inches == layout$point_x_inches), "Label x coordinate shifted.")
  list(layout = layout, text_bounds = boxes, overlap_count = 0L, clipped_count = 0L, leader_line_count = 0L)
}

check_panel_a_marker_reference <- function(layout, reference_dir) {
  path <- file.path(reference_dir, "PanelA_label_layout.csv")
  panel_a_assert(file.exists(path), paste("Missing marker reference:", path))
  ref <- utils::read.csv(path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  panel_a_assert(nrow(layout) == nrow(ref), "Marker reference row count differs.")
  for (name in c("row_index", "context", "residue", "position", "color", "marker"))
    panel_a_assert(identical(as.character(layout[[name]]), as.character(ref[[name]])),
                   paste("Marker reference mismatch:", name))
  for (name in c("point_x_inches", "point_y_inches"))
    panel_a_assert(all(abs(layout[[name]] - ref[[name]]) < 1e-10), paste("Marker coordinate mismatch:", name))
  # Different devices/fonts may assign a different collision-free label tier.
  # The frozen label-tier values are never used to place the R labels.
  sprintf("PASS marker coordinates, identities, shapes and colors: %d markers", nrow(layout))
}

make_panel_a <- function(workbook, output_dir = "output", sheet = NULL,
                         reference_dir = NULL, font_family = NULL,
                         style = panel_a_style(), save_svg = TRUE, overwrite = FALSE) {
  panel_a_require(c("readxl", "ragg", "systemfonts", "digest", if (save_svg) "svglite"))
  result <- analyse_panel_a(workbook, sheet)
  checks <- if (is.null(reference_dir)) "Not requested (reference_dir not supplied)." else
    check_panel_a_reference(result, reference_dir)
  font <- choose_panel_a_font(font_family)
  panel_a_assert(style$width_in == 12 && style$height_in == 6.6,
                 "The refined layout uses a fixed 12 x 6.6-inch canvas.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  panel_a_assert(dir.exists(output_dir), "Could not create output directory.")
  stem <- "RAD51D_SuppFigS3_PanelA_refined_R"
  tiff_name <- sprintf("%s_%ddpi.tiff", stem, style$dpi)
  names_out <- c(tiff_name, paste0(stem, "_preview.png"), if (save_svg) paste0(stem, ".svg"),
    unname(panel_a_table_files), "PanelA_label_layout.csv", "PanelA_text_bounds.csv",
    "PanelA_provenance.txt", "sessionInfo.txt")
  if (!overwrite && any(file.exists(file.path(output_dir, names_out))))
    stop("Output files exist. Choose a new output_dir or set overwrite = TRUE.", call. = FALSE)
  stage <- tempfile("panel_a_stage_", tmpdir = output_dir)
  dir.create(stage); on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  # Write to staging first: no final figure is delivered when an audit fails.
  export <- function(open_device) {
    open_device(); device <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(device), add = TRUE)
    draw_panel_a(result, font, style)
  }
  tiff_qa <- export(function() ragg::agg_tiff(file.path(stage, tiff_name),
    width = style$width_in, height = style$height_in, units = "in", res = style$dpi,
    compression = "lzw", bitsize = 8, background = "white"))
  preview_qa <- export(function() ragg::agg_png(file.path(stage, paste0(stem, "_preview.png")),
    width = style$width_in, height = style$height_in, units = "in", res = 150, background = "white"))
  if (save_svg) svg_qa <- export(function() svglite::svglite(file.path(stage, paste0(stem, ".svg")),
    width = style$width_in, height = style$height_in, bg = "white"))
  if (!is.null(reference_dir)) checks <- c(checks, check_panel_a_marker_reference(tiff_qa$layout, reference_dir))
  for (key in names(panel_a_table_files)) panel_a_write_csv(result[[key]], file.path(stage, panel_a_table_files[[key]]))
  panel_a_write_csv(tiff_qa$layout, file.path(stage, "PanelA_label_layout.csv"))
  panel_a_write_csv(tiff_qa$text_bounds, file.path(stage, "PanelA_text_bounds.csv"))
  provenance <- c(
    paste("Script version:", PANEL_A_VERSION),
    paste("Run UTC:", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
    paste("Input workbook:", basename(workbook)), paste("Input sheet:", result$source_sheet),
    paste("Input SHA-256:", digest::digest(file = workbook, algo = "sha256")),
    paste("Category column:", result$category_column),
    "Structural columns: U ('Greenhouph Nature 2023') and V ('Rawal Nature 2023').",
    "Red: at least one assayed missense SNV at the position has explicit category 'P strong'.",
    "Black: no current P_strong missense hit; NOT a benign classification.",
    "Union nonempty structural annotations across missense SNVs at a reference position.",
    "Distinct SNVs remain distinct; different SNVs with the same protein change are not deduplicated.",
    "No posterior re-thresholding, SpliceAI exclusion, imputation, or VarCall model refit.",
    "Reference CSVs are used only for validation, not to draw the graph.",
    "No new structural/literature verification is performed by this visualization script.",
    "Row titles reproduce the reference categories; they do not establish essentiality or hydrolysis by RAD51D.",
    paste("Missense-labeled SNVs:", nrow(result$snvs)),
    paste("P_strong missense-labeled SNVs:", sum(result$snvs$is_P_strong)),
    paste("All annotated positions:", nrow(result$annotated)),
    paste("Plotted residue/context markers:", nrow(result$plotted)),
    paste("Plotted unique positions:", length(unique(result$plotted$position))),
    paste("Plotted unique P_strong positions:", length(unique(result$plotted$position[result$plotted$P_strong_hit == 1L]))),
    paste("Annotated positions outside displayed categories:", panel_a_join(result$omitted$residue)),
    "Q115 appears in upper ATP and RAD51C rows because both contexts occur in the workbook.",
    "I311 appears in the upper ATP row because the workbook states 'RAD51D (ATP, Non-bonded)'.",
    paste("Font family:", font),
    sprintf("Font points (all bold): residues=%g; row labels=%g; ticks=%g; titles=%g; legend=%g.",
      style$residue_size_pt, style$row_size_pt, style$tick_size_pt, style$title_size_pt, style$legend_size_pt),
    sprintf("TIFF: %g x %g inches, %g dpi, %d x %d pixels, LZW compression.",
      style$width_in, style$height_in, style$dpi,
      round(style$width_in * style$dpi), round(style$height_in * style$dpi)),
    "No leader lines. Label x coordinates equal marker x coordinates.",
    "Adjacent categorical markers have +/-0.030-inch y offsets (no quantitative meaning).",
    "TIFF, preview, and requested SVG passed text overlap and clipping checks.",
    "Label tiers are computed on each output device; minor device/font differences are possible.",
    "", "Reference validation:", checks)
  writeLines(provenance, file.path(stage, "PanelA_provenance.txt"), useBytes = TRUE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(stage, "sessionInfo.txt"))
  copied <- file.copy(file.path(stage, names_out), file.path(output_dir, names_out), overwrite = overwrite)
  panel_a_assert(all(copied), "Could not copy all validated outputs to output_dir.")
  result$rendering <- tiff_qa; result$font_family <- font; result$reference_checks <- checks
  message("Wrote ", file.path(output_dir, tiff_name))
  message(sprintf("%d markers across %d unique positions; %d unique P_strong-hit positions.",
    nrow(result$plotted), length(unique(result$plotted$position)),
    length(unique(result$plotted$position[result$plotted$P_strong_hit == 1L]))))
  invisible(result)
}

panel_a_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Usage: Rscript scripts/panel_a.R --workbook FILE.xlsx [--outdir output]",
    "       [--sheet NAME] [--reference-dir validation/reference] [--font FAMILY]",
    "       [--no-svg] [--overwrite]", sep = "\n")
  if (!length(args) || "--help" %in% args) { cat(usage, "\n"); return(invisible(NULL)) }
  opt <- list(output_dir = "output", save_svg = TRUE, overwrite = FALSE)
  flags <- c("--workbook" = "workbook", "--outdir" = "output_dir", "--sheet" = "sheet",
    "--reference-dir" = "reference_dir", "--font" = "font_family")
  i <- 1L
  while (i <= length(args)) {
    flag <- args[[i]]
    if (flag == "--no-svg") opt$save_svg <- FALSE else if (flag == "--overwrite") opt$overwrite <- TRUE else {
      panel_a_assert(flag %in% names(flags), paste("Unknown option:", flag, "\n", usage))
      panel_a_assert(i < length(args), paste("Missing value for", flag))
      i <- i + 1L
      panel_a_assert(!startsWith(args[[i]], "--"), paste("Missing value for", flag))
      opt[[flags[[flag]]]] <- args[[i]]
    }
    i <- i + 1L
  }
  panel_a_assert(!is.null(opt$workbook), "--workbook is required.")
  do.call(make_panel_a, opt)
}

# source() only defines functions. Direct Rscript execution runs the CLI.
if (sys.nframe() == 0L) panel_a_main()
