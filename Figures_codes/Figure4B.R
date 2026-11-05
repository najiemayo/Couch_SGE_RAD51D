rad51d_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop(paste("Missing R packages:", paste(missing, collapse = ", ")), call. = FALSE)
}

rad51d_check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

rad51d_clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", ".", "NA", "N/A")] <- NA_character_
  x
}

rad51d_read <- function(workbook, sheet = NULL) {
  rad51d_require("readxl")
  rad51d_check(file.exists(workbook), "Workbook not found.")
  header_key <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(as.character(x))))
  sheets <- readxl::excel_sheets(workbook)
  if (!is.null(sheet)) {
    if (is.numeric(sheet)) {
      rad51d_check(length(sheet) == 1L && sheet %in% seq_along(sheets), "Invalid sheet index.")
      sheet <- sheets[[sheet]]
    }
    rad51d_check(length(sheet) == 1L && sheet %in% sheets, "Sheet not found.")
    sheets <- sheet
  }
  candidates <- list()
  for (s in sheets) {
    preview <- as.matrix(readxl::read_excel(workbook, sheet = s, col_names = FALSE,
      col_types = "text", range = readxl::cell_limits(c(1L, 1L), c(20L, NA)),
      .name_repair = "minimal", progress = FALSE))
    hits <- which(apply(preview, 1L, function(row)
      all(c("chrom", "poshg38", "ref", "alt", "variantclass") %in% header_key(row))))
    for (h in hits) candidates[[length(candidates) + 1L]] <- list(sheet = s, row = h)
  }
  rad51d_check(length(candidates) == 1L, "A unique variant-table header is required; specify sheet if necessary.")
  selected <- candidates[[1L]]
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = selected$sheet,
    range = readxl::cell_limits(c(selected$row, 1L), c(NA, NA)), col_types = "text",
    .name_repair = "minimal", progress = FALSE), stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- header_key(names(raw))
  rad51d_check(!anyDuplicated(names(raw)), "Duplicate normalized column names.")
  raw$source_row <- seq_len(nrow(raw)) + selected$row
  chrom <- rad51d_clean(raw$chrom)
  keep <- !is.na(chrom) & grepl("^(chr)?([0-9]+|X|Y|M|MT)$", chrom, ignore.case = TRUE)
  raw <- raw[keep, , drop = FALSE]
  rad51d_check(nrow(raw) > 0L, "No genomic variant rows found.")
  field <- function(name, required = TRUE) {
    if (!name %in% names(raw)) {
      if (required) stop(paste("Missing column:", name), call. = FALSE)
      return(rep(NA_character_, nrow(raw)))
    }
    rad51d_clean(raw[[name]])
  }
  category_column <- grep("^(varcallfunctionalcategory|sgefunctioncategoryofvarcall)[0-9]*$",
    names(raw), value = TRUE)
  rad51d_check(length(category_column) == 1L, "One explicit VarCall functional-category column is required.")
  categories <- c("P strong", "P moderate", "P supporting", "VUS", "B supporting", "B moderate", "B strong")
  category <- tolower(gsub("[_[:space:]]+", " ", field(category_column)))
  category <- categories[match(category, tolower(categories))]
  rad51d_check(!anyNA(category), "Missing or unrecognized VarCall functional categories.")
  pos <- suppressWarnings(as.numeric(field("poshg38")))
  ref <- toupper(field("ref")); alt <- toupper(field("alt"))
  rad51d_check(all(is.finite(pos) & pos > 0 & pos == floor(pos)), "Invalid genomic position.")
  rad51d_check(all(ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") & ref != alt),
    "Each variant must be a single-nucleotide substitution.")
  chrom <- paste0("chr", toupper(sub("^chr", "", field("chrom"), ignore.case = TRUE)))
  out <- data.frame(variant_id = paste(chrom, format(pos, scientific = FALSE, trim = TRUE), ref, alt, sep = "_"),
    cdna = field("c"), protein = field("p"), aa = field("aa"), type = field("variantclass"),
    exon = field("sge"), category = category, source_row = raw$source_row,
    domain = field("proteindomain", FALSE), greenhough = field("greenhouphnature2023", FALSE),
    rawal = field("rawalnature2023", FALSE), stringsAsFactors = FALSE)
  if ("greenhoughnature2023" %in% names(raw)) out$greenhough <- field("greenhoughnature2023")
  rad51d_check(!anyDuplicated(out$variant_id), "Duplicate genomic variants; records were not collapsed.")
  attr(out, "has_domains") <- "proteindomain" %in% names(raw)
  attr(out, "has_structural_annotations") <-
    any(c("greenhouphnature2023", "greenhoughnature2023") %in% names(raw)) && "rawalnature2023" %in% names(raw)
  out
}

rad51d_missense <- function(data) {
  out <- data[!is.na(data$type) & tolower(data$type) == "missense", , drop = FALSE]
  pattern <- "^(Ter|[A-Z*])([0-9]+)(Ter|[A-Z?*])$"
  rad51d_check(nrow(out) > 0L && all(!is.na(out$aa) & grepl(pattern, out$aa)), "Invalid missense amino-acid annotations.")
  out$position <- as.integer(sub(pattern, "\\2", out$aa))
  out$reference_aa <- sub(pattern, "\\1", out$aa)
  rad51d_check(all(out$position >= 1L & out$position <= 329L), "Amino-acid positions must be in 1-329.")
  rad51d_check(all(!is.na(out$cdna) & grepl("^c\\.[0-9]+[ACGT]>[ACGT]$", out$cdna)), "Invalid coding-SNV notation.")
  coding <- as.integer(sub("^c\\.([0-9]+).*", "\\1", out$cdna))
  rad51d_check(all((coding + 2L) %/% 3L == out$position), "cDNA and amino-acid positions disagree.")
  rad51d_check(!anyDuplicated(out$cdna), "Duplicate missense cDNA variants.")
  references <- split(out$reference_aa, out$position)
  rad51d_check(all(vapply(references, function(x) length(unique(x)) == 1L, logical(1))),
    "Conflicting reference amino acids at the same position.")
  out[order(out$position, out$variant_id), , drop = FALSE]
}

rad51d_export <- function(draw, output_dir, stem, width, height, dpi) {
  rad51d_require(c("ragg", "svglite"))
  rad51d_check(length(dpi) == 1L && is.finite(dpi) && dpi >= 72 && dpi == as.integer(dpi), "Invalid DPI.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rad51d_check(dir.exists(output_dir), "Output directory could not be created.")
  tiff <- file.path(output_dir, paste0(stem, "_", dpi, "dpi.tiff"))
  svg <- file.path(output_dir, paste0(stem, ".svg"))
  ragg::agg_tiff(tiff, width = width, height = height, units = "in", res = dpi,
    background = "white", compression = "lzw")
  tryCatch(draw(), finally = grDevices::dev.off())
  svglite::svglite(svg, width = width, height = height, bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(c(tiff = tiff, svg = svg))
}
# SNV-level numerator: explicit P strong calls; denominator: all missense-annotated SNVs.
# Start/stop-loss records are retained as annotated. No SpliceAI filter is applied.
figure4b_data <- function(workbook, sheet = NULL) {
  data <- rad51d_read(workbook, sheet)
  missense <- rad51d_missense(data)
  n <- tabulate(missense$position, nbins = 329L)
  hits <- tabulate(missense$position[missense$category == "P strong"], nbins = 329L)
  positions <- data.frame(position = seq_len(329L), n_missense_SNVs = n, n_P_strong_SNVs = hits,
    P_strong_percent = ifelse(n > 0L, 100 * hits / n, NA_real_))
  rad51d_check(sum(n) == nrow(missense) && sum(hits) == sum(missense$category == "P strong"),
    "Position-level totals do not match variant-level totals.")

  # Data 4 omits domains. These fixed intervals reproduce Data 2, column T.
  domains <- data.frame(domain = c("NTD", "RecA-like", "Walker A", "Walker B"),
    start = c(13L, 101L, 107L, 202L), end = c(61L, 310L, 114L, 206L))
  if (isTRUE(attr(data, "has_domains"))) {
    classes <- list("NTD", c("RecA-like", "RecA-like_Walker A", "RecA-like_Walker B"),
      "RecA-like_Walker A", "RecA-like_Walker B")
    for (i in seq_len(nrow(domains))) {
      p <- missense$position[missense$domain %in% classes[[i]]]
      rad51d_check(length(p) > 0L, paste("No observed positions for domain", domains$domain[i]))
      domains$start[i] <- min(p); domains$end[i] <- max(p)
    }
  }
  exon_names <- unique(data$exon)
  rad51d_check(setequal(exon_names, paste0("E", 1:10)), "Expected assay exons E1-E10.")
  exons <- do.call(rbind, lapply(paste0("E", 1:10), function(exon) {
    cdna <- data$cdna[data$exon == exon & !is.na(data$cdna)]
    coding <- as.integer(sub("^c\\.([0-9]+).*", "\\1", cdna[grepl("^c\\.[0-9]+[ACGT]>[ACGT]$", cdna)]))
    donor <- unique(as.integer(sub("^c\\.([0-9]+).*", "\\1", cdna[grepl("^c\\.[0-9]+\\+[0-9]+", cdna)])))
    acceptor <- unique(as.integer(sub("^c\\.([0-9]+).*", "\\1", cdna[grepl("^c\\.[0-9]+-[0-9]+", cdna)])))
    rad51d_check(length(coding) > 0L && length(donor) <= 1L && length(acceptor) <= 1L,
      paste("Ambiguous coding boundaries for", exon))
    start <- if (length(acceptor)) acceptor else min(coding)
    end <- if (length(donor)) donor else max(coding)
    rad51d_check(all(coding >= start & coding <= end), paste("Coding positions outside exon", exon))
    data.frame(exon = exon, coding_start = start, coding_end = end,
      plot_left = (start - 1) / 3 + 0.5, plot_right = end / 3 + 0.5)
  }))
  rad51d_check(all(head(exons$coding_end, -1L) + 1L == tail(exons$coding_start, -1L)) &&
    exons$coding_start[1L] == 1L && tail(exons$coding_end, 1L) == 987L,
    "Exon boundaries do not cover the expected c.1-c.987 interval.")
  list(positions = positions, variants = missense, domains = domains, exons = exons)
}

figure4b_draw <- function(result, font_family = "sans", panel_letter = "b") {
  rad51d_check(max(result$positions$P_strong_percent, na.rm = TRUE) <= 80,
    "A value exceeds the 0-80% axis; revise the plotting limit before export.")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 7.5), yscale = c(0, 2.2), clip = "off"))
  on.exit(grid::popViewport())
  # Physical x coordinates preserve alignment to the existing Figure 4a artwork.
  x <- function(p) 0.0833333333333333 + (823.300703998283 + 10.7776046145027 * p) * 7.34375 / 4396
  y <- function(v) 0.36 + v / 80 * 1.16
  text <- function(label, xx, yy, size = 7.8, just = "centre", rot = 0, face = "bold")
    grid::grid.text(label, x = xx, y = yy, default.units = "native", just = just, rot = rot,
      gp = grid::gpar(fontfamily = font_family, fontsize = size, fontface = face, lineheight = 1.05))
  line <- function(x0, y0, x1, y1, width = 0.8)
    grid::grid.segments(x0, y0, x1, y1, default.units = "native",
      gp = grid::gpar(lwd = width * 96 / 72, lineend = "butt"))
  p <- result$positions$position; values <- result$positions$P_strong_percent
  finite <- which(is.finite(values))
  groups <- split(finite, cumsum(c(TRUE, diff(finite) != 1L)))
  for (indices in groups) {
    xx <- x(p[indices]); yy <- y(values[indices])
    if (length(indices) == 1L) {
      line(xx, y(0), xx, yy, 0.42)
    } else {
      grid::grid.polygon(c(xx[1L], xx, tail(xx, 1L)), c(y(0), yy, y(0)), default.units = "native",
        gp = grid::gpar(fill = "black", col = NA))
      grid::grid.lines(xx, yy, default.units = "native", gp = grid::gpar(lwd = 0.42 * 96 / 72))
    }
  }
  left <- x(0.5); right <- x(329.5)
  line(left, y(0), right, y(0)); line(left, y(0), left, y(80))
  xticks <- c(1, 50, 100, 150, 200, 250, 300, 329)
  for (p in setdiff(seq(10, 320, 10), xticks)) line(x(p), y(0), x(p), y(0) - 1.3 / 72, 0.6)
  for (p in xticks) {
    line(x(p), y(0), x(p), y(0) - 2.4 / 72, 0.75)
    text(as.character(p), x(p), 0.275, size = 7)
  }
  for (v in seq(0, 80, 20)) {
    line(left, y(v), left - 2.4 / 72, y(v), 0.75)
    text(as.character(v), left - 4.8 / 72, y(v), size = 7, just = "right")
  }
  text("Amino Acid Position", (left + right) / 2, 0.105)
  text("P_strong missense\nSNVs (%)", left - 0.43, y(40), rot = 90)
  for (i in seq_len(nrow(result$exons))) {
    e <- result$exons[i, ]
    grid::grid.rect(x(e$plot_left), 2.01, width = x(e$plot_right) - x(e$plot_left), height = 0.07,
      default.units = "native", just = c("left", "bottom"),
      gp = grid::gpar(fill = if (i %% 2L) "#EFF1F4" else "#C8CDD1", col = "#656A70", lwd = 0.55 * 96 / 72))
  }
  for (i in seq_len(nrow(result$domains))) {
    d <- result$domains[i, ]; motif <- startsWith(d$domain, "Walker")
    yy <- if (motif) 1.665 else 1.785
    line(x(d$start - 0.5), yy, x(d$end + 0.5), yy, 0.75)
    text(d$domain, x((d$start + d$end) / 2), if (motif) 1.64 else 1.825,
      just = c("centre", if (motif) "top" else "bottom"))
  }
  text("RAD51D", 0.87, 2.045, 9, "left", face = "bold.italic")
  if (nzchar(panel_letter)) text(panel_letter, 0.17, 2.04, 18, "left")
}

make_figure4b <- function(workbook, output_dir = "Figure4B_output", sheet = NULL,
                         dpi = 600L, font_family = "sans", panel_letter = "b") {
  result <- figure4b_data(workbook, sheet)
  paths <- rad51d_export(function() figure4b_draw(result, font_family, panel_letter),
    output_dir, "Figure4B", width = 7.5, height = 2.2, dpi = dpi)
  for (name in names(result)) utils::write.csv(result[[name]],
    file.path(output_dir, paste0("Figure4B_", name, ".csv")), row.names = FALSE, na = "")
  invisible(list(files = paths, data = result))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  workbook <- if (length(args)) args[[1L]] else if (interactive()) file.choose() else stop("Workbook path required.")
  output_dir <- if (length(args) >= 2L) args[[2L]] else "Figure4B_output"
  dpi <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 600L
  make_figure4b(workbook, output_dir = output_dir, dpi = dpi)
}
