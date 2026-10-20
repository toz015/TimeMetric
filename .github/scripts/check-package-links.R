# Check local link targets in the built archive, not the development checkout.
# Remote URLs are checked separately with urlchecker::url_check().
check_archive <- function(archive) {
  if (!nzchar(Sys.which("pandoc"))) stop("Pandoc is required; do not skip Markdown links.")
  if (!requireNamespace("xml2", quietly = TRUE)) stop("xml2 is required.")
  unpack <- tempfile("package-links-")
  dir.create(unpack)
  on.exit(unlink(unpack, recursive = TRUE), add = TRUE)
  utils::untar(archive, exdir = unpack)
  descriptions <- list.files(unpack, pattern = "^DESCRIPTION$", recursive = TRUE,
                              full.names = TRUE)
  if (length(descriptions) != 1L) stop("Expected exactly one package DESCRIPTION.")
  root <- dirname(descriptions)
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  links <- data.frame(Parent = character(), URI = character())
  add <- function(parent, uri) {
    if (length(uri)) links <<- rbind(links,
      data.frame(Parent = rep(parent, length(uri)), URI = uri))
  }
  for (file in files[grepl("\\.(md|html?)$", files, ignore.case = TRUE)]) {
    parent <- substring(file, nchar(root) + 2L)
    html <- file
    if (grepl("\\.md$", file, ignore.case = TRUE)) {
      html <- tempfile(fileext = ".html", tmpdir = unpack)
      status <- system2(Sys.which("pandoc"),
        c("--from=gfm", "--to=html", "--output", shQuote(html), shQuote(file)))
      if (status != 0L) stop("Could not parse Markdown: ", parent)
    }
    doc <- xml2::read_html(html)
    for (attribute in c("href", "src")) {
      add(parent, xml2::xml_attr(xml2::xml_find_all(doc,
        paste0("//*[@", attribute, "]")), attribute))
    }
  }
  meta <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
  # DESCRIPTION URL fields allow comma/whitespace-separated targets, optionally
  # followed by parenthesized labels. Include local targets, not just https URLs.
  for (field in intersect(c("URL", "BugReports"), names(meta))) {
    values <- gsub("[[:space:]]+\\([^)]*\\)", "", meta[[field]])
    add("DESCRIPTION", strsplit(trimws(values), "[[:space:],]+")[[1]])
  }
  if ("Description" %in% names(meta)) {
    tokens <- regmatches(meta[["Description"]],
      gregexpr("<[^<>]+>", meta[["Description"]]))[[1]]
    add("DESCRIPTION", sub(">$", "", sub("^<", "", tokens)))
  }
  citation_file <- file.path(root, "inst", "CITATION")
  if (file.exists(citation_file)) {
    citations <- utils::readCitationFile(citation_file, meta = as.list(meta))
    add("inst/CITATION", unlist(citations$url, use.names = FALSE))
  }
  rd <- tools::Rd_db(dir = root)
  walk_rd <- function(node, parent) {
    tag <- attr(node, "Rd_tag")
    if (identical(tag, "\\url") || identical(tag, "\\href")) {
      value <- if (identical(tag, "\\href")) node[[1L]] else node
      add(parent, paste(unlist(value), collapse = ""))
    } else if (is.list(node)) {
      for (child in node) walk_rd(child, parent)
    }
  }
  for (name in names(rd)) walk_rd(rd[[name]], file.path("man", name))
  links <- unique(links)
  bad <- logical(nrow(links))
  for (i in seq_len(nrow(links))) {
    uri <- trimws(links$URI[i])
    if (!nzchar(uri) || startsWith(uri, "#")) next
    # Check drive letters before treating the prefix as a URI scheme.
    # Decode first so encoded separators cannot hide a host-specific path.
    decoded <- utils::URLdecode(uri)
    if (grepl("^[A-Za-z]:", decoded) || startsWith(decoded, "\\") ||
        startsWith(decoded, "/") || startsWith(decoded, "~")) {
      bad[i] <- TRUE
      next
    }
    scheme <- grepl("^[A-Za-z][A-Za-z0-9+.-]*:", uri)
    if (scheme && !grepl("^file:", uri, ignore.case = TRUE)) next
    # Host-specific absolute paths and file: URLs cannot work for recipients.
    if (scheme || startsWith(uri, "/")) {
      bad[i] <- TRUE
      next
    }
    target <- utils::URLdecode(sub("[?#].*$", "", uri))
    resolved <- normalizePath(file.path(root, dirname(links$Parent[i]), target),
                              winslash = "/", mustWork = FALSE)
    canonical_root <- normalizePath(root, winslash = "/")
    bad[i] <- !file.exists(resolved) ||
      !(identical(resolved, canonical_root) || startsWith(resolved, paste0(canonical_root, "/")))
  }
  cat("Inspected", nrow(links), "unique document links in", basename(archive), "\n")
  if (any(bad)) {
    print(links[bad, ], row.names = FALSE)
    stop("Invalid local file URIs in the built package.", call. = FALSE)
  }
  cat("PASS: no invalid local file URIs in the built package.\n")
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L || !file.exists(args[1]))
    stop("Usage: Rscript .github/scripts/check-package-links.R PACKAGE_VERSION.tar.gz")
  check_archive(normalizePath(args[1], mustWork = TRUE))
}
