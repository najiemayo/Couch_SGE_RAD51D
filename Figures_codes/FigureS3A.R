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
# Fixed structural annotations from Supplementary Data 2, columns U/V (092826).
# Used only when the input workbook omits structural annotations; no functional calls are stored here.
figure_s3a_reference <- function() {
  utils::read.delim(text = "position	reference_aa	greenhough	rawal
5	R	RAD51C (Ionic bond)	
20	L	Predicted monomer instability	
25	I	Predicted monomer instability	
31	L	Predicted monomer instability	
39	V	Predicted monomer instability	
53	L	Predicted monomer instability	
62	S	XRCC2 (H-bond)	
66	V	XRCC2 (H-bond)	
67	N	XRCC2 (H-bond)	
68	G	XRCC2 (H-bond)	
89	L	Predicted monomer instability	
93	L	Predicted monomer instability	
102	V	Predicted monomer instability	
110	G	ATP (RAD51D, Non-bonded)	
111	S	ATP (RAD51D, H-bond)	
112	G	ATP (RAD51D, H-bond)	
113	K	ATP (RAD51D, H-bond)	ATP binding
114	T	ATP (RAD51D, H-bond)	Mg ion coordination
115	Q	RAD51C (H-bond), ATP (RAD51D, H-bond)	ATPase site
132	V	Predicted monomer instability	
135	V	Predicted monomer instability	
137	S	RAD51C (H-bond)	
145	R	ATP (RAD51D, Non-bonded)	ATPase site
146	L	Predicted monomer instability, RAD51C (H-bond); RAD51C (H-bond)	
149	L	Predicted monomer instability	
150	L	Predicted monomer instability	
164	L	RAD51C (H-bond)	
167	I	RAD51C (H-bond)	
169	V	RAD51C (H-bond)	
170	V		polymerization motif binding pocket
171	H	RAD51C (H-bond)	
176	F		intermolecular interactions within BCDX2
181	V		polymerization motif binding pocket
188	T		polymerization motif binding pocket
206	D	ATP (RAD51D, H-bond)	
211	V	Predicted monomer instability	
216	L		ssDNA binding
217	G		ssDNA binding
221	R		ssDNA binding
224	L		ssDNA binding
226	L		intermolecular interactions within BCDX2
232	R	XRCC2 (H-bond)	
234	L	Predicted monomer instability	
235	K	XRCC2 (H-bond)	
239	R	XRCC2 (H-bond)	
251	I		ssDNA binding
252	T		ssDNA binding
253	R		ssDNA binding
255	R		ssDNA binding
263	A		ssDNA binding
264	L	Predicted monomer instability	ssDNA binding
265	G		ssDNA binding
266	R		ssDNA binding
267	S	XRCC2 (H-bond)	
270	F	XRCC2 (H-bond)	
278	L	Predicted monomer instability	
297	K	ATP (XRCC2, H-bond)	ATP binding
299	S	XRCC2 (H-bond)	
300	R	ATP (XRCC2, H-bond); XRCC2 (H-bond, ionic), ATP (XRCC2, H-bond)	
301	Q	ATP (XRCC2, H-bond)	
302	P	ATP (XRCC2, Non-bonded)	ATPase site
303	T	ATP (XRCC2, Non-bonded)	
311	I	RAD51D (ATP, Non-bonded)	", sep = "\t", quote = "", na.strings = "", stringsAsFactors = FALSE, check.names = FALSE)
}

# A red marker requires at least one explicit P strong missense SNV at that residue.
# Black denotes no P strong call, not functional neutrality. No SpliceAI filter is applied.
figure_s3a_data <- function(workbook, sheet = NULL) {
  data <- rad51d_read(workbook, sheet)
  missense <- rad51d_missense(data)
  if (isTRUE(attr(data, "has_structural_annotations"))) {
    annotated <- missense[!is.na(missense$greenhough) | !is.na(missense$rawal),
      c("position", "reference_aa", "greenhough", "rawal"), drop = FALSE]
    annotated <- unique(annotated)
  } else {
    annotated <- figure_s3a_reference()
  }
  rad51d_check(nrow(annotated) > 0L, "No structural annotations found.")
  rows <- list()
  contexts <- c("Interaction with ATP", "Interaction with XRCC2", "Interaction with RAD51C",
    "Monomer stability", "ATP binding/hydrolysis", "ssDNA binding")
  for (p in sort(unique(annotated$position))) {
    a <- annotated[annotated$position == p, , drop = FALSE]
    snvs <- missense[missense$position == p, , drop = FALSE]
    rad51d_check(nrow(snvs) > 0L, paste("No assayed missense SNVs at annotated position", p))
    rad51d_check(length(unique(a$reference_aa)) == 1L &&
      unique(a$reference_aa) == unique(snvs$reference_aa), paste("Reference amino-acid mismatch at", p))
    g <- as.character(unique(unlist(strsplit(a$greenhough[!is.na(a$greenhough)], "; ", fixed = TRUE))))
    r <- as.character(unique(unlist(strsplit(a$rawal[!is.na(a$rawal)], "; ", fixed = TRUE))))
    matches <- c(any(grepl("ATP", g, fixed = TRUE)), any(startsWith(g, "XRCC2 (")),
      any(grepl("RAD51C (", g, fixed = TRUE)), any(grepl("Predicted monomer instability", g, fixed = TRUE)),
      any(r %in% c("ATP binding", "Mg ion coordination", "ATPase site")), any(r == "ssDNA binding"))
    hit <- snvs$category == "P strong"
    for (i in which(matches)) rows[[length(rows) + 1L]] <- data.frame(
      context_id = i, context = contexts[i],
      source = if (i <= 4L) "Greenhough et al. (2023)" else "Rawal et al. (2023)",
      position = p, residue = paste0(unique(snvs$reference_aa), p),
      n_missense_SNVs = nrow(snvs), n_P_strong_SNVs = sum(hit), P_strong_hit = any(hit),
      P_strong_variants = paste(snvs$cdna[hit], collapse = "; "), stringsAsFactors = FALSE)
  }
  rad51d_check(length(rows) > 0L, "No residues match the six displayed structural contexts.")
  points <- do.call(rbind, rows)
  points <- points[order(points$context_id, points$position), , drop = FALSE]
  rownames(points) <- NULL
  rad51d_check(!anyDuplicated(paste(points$context_id, points$position)), "Duplicate residue-context markers.")
  list(points = points, variants = missense[missense$position %in% points$position, , drop = FALSE])
}

figure_s3a_draw <- function(result, font_family = "sans", panel_letter = "a") {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 12), yscale = c(0, 6.6), clip = "off"))
  on.exit(grid::popViewport())
  left <- 2.50; right <- 11.80
  x <- function(p) left + p / 330 * (right - left)
  boxes <- list()
  text_gp <- function(size) grid::gpar(fontfamily = font_family, fontface = "bold", fontsize = size, col = "black")
  bounds <- function(label, xx, yy, size, just = "centre") {
    grob <- grid::textGrob(label, gp = text_gp(size))
    w <- grid::convertWidth(grid::grobWidth(grob), "inches", valueOnly = TRUE)
    h <- grid::convertHeight(grid::grobHeight(grob), "inches", valueOnly = TRUE)
    j <- switch(just, left = 0, centre = 0.5, right = 1)
    c(x0 = xx - j * w, x1 = xx + (1 - j) * w, y0 = yy - h / 2, y1 = yy + h / 2)
  }
  overlaps <- function(a, b, gap = 0) {
    min(a[[2L]], b[[2L]]) + gap > max(a[[1L]], b[[1L]]) &&
      min(a[[4L]], b[[4L]]) + gap > max(a[[3L]], b[[3L]])
  }
  text <- function(label, xx, yy, size, just = "centre") {
    b <- bounds(label, xx, yy, size, just)
    rad51d_check(b[[1L]] >= 0 && b[[2L]] <= 12 && b[[3L]] >= 0 && b[[4L]] <= 6.6,
      paste("Text exceeds the canvas:", label))
    rad51d_check(!any(vapply(boxes, function(old) overlaps(b, old), logical(1))),
      paste("Overlapping text:", label))
    grid::grid.text(label, xx, yy, default.units = "native", just = c(just, "centre"), gp = text_gp(size))
    boxes[[length(boxes) + 1L]] <<- b
  }
  line <- function(x0, y0, x1, y1, width = 1)
    grid::grid.segments(x0, y0, x1, y1, default.units = "native",
      gp = grid::gpar(col = "black", lwd = width * 96 / 72, lineend = "butt"))
  marker <- function(xx, yy, context, color, size = 3.5) {
    shapes <- c(25L, 24L, 22L, 21L, 23L, 25L)
    grid::grid.points(xx, yy, default.units = "native", pch = shapes[context],
      size = grid::unit(size, "bigpts"), gp = grid::gpar(col = "white", fill = color, lwd = 0.20 * 96 / 72))
  }
  if (nzchar(panel_letter)) text(panel_letter, 0.09, 6.29, 25, "left")
  marker(left + 0.02, 6.36, 4L, "#FF0000", 4.1)
  text("At least one P_strong missense SNV", left + 0.14, 6.36, 10, "left")
  marker(left + 5.13, 6.36, 4L, "black", 4.1)
  text("No P_strong missense SNV", left + 5.25, 6.36, 10, "left")
  text("Greenhough et al. (2023)", right, 6.10, 10, "right")
  text("Rawal et al. (2023)", right, 2.19, 10, "right")
  for (b in list(c(2.95, 5.98), c(0.715, 2.08)))
    grid::grid.rect(left, b[1L], width = right - left, height = diff(b), default.units = "native",
      just = c("left", "bottom"), gp = grid::gpar(fill = NA, col = "black", lwd = 0.9 * 96 / 72))
  row_y <- c(5.46, 4.74, 4.02, 3.30, 1.70, 1.085)
  context_names <- c("Interaction with ATP", "Interaction with XRCC2", "Interaction with RAD51C",
    "Monomer stability", "ATP binding/hydrolysis", "ssDNA binding")
  for (i in seq_len(6L)) text(context_names[i], left - 0.095, row_y[i], 11.5, "right")
  for (a in list(c(2.805, 2.605, 2.38), c(0.59, 0.39, 0.155))) {
    line(left, a[1L], right, a[1L], 1.3)
    for (p in seq(0, 330, 10)) {
      major <- p %% 50 == 0
      line(x(p), a[1L], x(p), a[1L] - if (major) 0.070 else 0.035, if (major) 1.15 else 0.8)
      if (major) text(as.character(p), x(p), a[2L], 11.5)
    }
    text("AA position", (left + right) / 2, a[3L], 13)
  }
  residue_size <- 9.25
  offsets <- c(0.130, -0.130, 0.275, -0.275, 0.420, -0.420)
  for (i in seq_len(6L)) {
    points <- result$points[result$points$context_id == i, , drop = FALSE]
    lanes <- vector("list", length(offsets))
    for (j in seq_len(nrow(points))) {
      p <- points[j, ]; xx <- x(p$position); chosen <- NA_integer_
      for (level in seq_along(offsets)) {
        b <- bounds(p$residue, xx, row_y[i] + offsets[level], residue_size)
        span_clear <- !any(vapply(lanes[[level]], function(old)
          min(b[[2L]], old[[2L]]) + 0.030 > max(b[[1L]], old[[1L]]), logical(1)))
        clear <- !any(vapply(boxes, function(old) overlaps(b, old, gap = 0.005), logical(1)))
        if (span_clear && clear) { chosen <- level; break }
      }
      rad51d_check(!is.na(chosen), paste("No collision-free label position for", p$residue, "in", p$context))
      lanes[[chosen]][[length(lanes[[chosen]]) + 1L]] <- b
      nearby <- (j > 1L && p$position - points$position[j - 1L] <= 1L) ||
        (j < nrow(points) && points$position[j + 1L] - p$position <= 1L)
      yy <- row_y[i] + if (nearby) { if (j %% 2L == 1L) 0.030 else -0.030 } else 0
      marker(xx, yy, i, if (p$P_strong_hit) "#FF0000" else "black")
      text(p$residue, xx, row_y[i] + offsets[chosen], residue_size)
    }
  }
}

make_figure_s3a <- function(workbook, output_dir = "FigureS3A_output", sheet = NULL,
                            dpi = 600L, font_family = "sans", panel_letter = "a") {
  result <- figure_s3a_data(workbook, sheet)
  paths <- rad51d_export(function() figure_s3a_draw(result, font_family, panel_letter),
    output_dir, "FigureS3A", width = 12, height = 6.6, dpi = dpi)
  for (name in names(result)) utils::write.csv(result[[name]],
    file.path(output_dir, paste0("FigureS3A_", name, ".csv")), row.names = FALSE, na = "")
  invisible(list(files = paths, data = result))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  workbook <- if (length(args)) args[[1L]] else if (interactive()) file.choose() else stop("Workbook path required.")
  output_dir <- if (length(args) >= 2L) args[[2L]] else "FigureS3A_output"
  dpi <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 600L
  make_figure_s3a(workbook, output_dir = output_dir, dpi = dpi)
}
