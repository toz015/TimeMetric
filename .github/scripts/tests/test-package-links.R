library(testthat)
source(file.path('..', 'check-package-links.R'), local = TRUE)

# Build actual archives: a checkout-only test would miss the original defect.
fixture_archive <- function(readme = '[News](NEWS.md)', url = 'https://example.org',
                            bugs = 'https://example.org/issues',
                            citation_url = 'https://example.org',
                            description = 'A link audit fixture.', extra = NULL) {
  base <- tempfile('link-fixture-')
  dir.create(base)
  withr::defer(unlink(base, recursive = TRUE), envir = parent.frame())
  root <- file.path(base, 'fixture')
  dir.create(root)
  dir.create(file.path(root, 'inst'))
  writeLines(c('Package: fixture', 'Version: 1.0.0',
    paste('Description:', description), paste('URL:', url),
    paste('BugReports:', bugs)), file.path(root, 'DESCRIPTION'))
  writeLines(readme, file.path(root, 'README.md'))
  writeLines('News', file.path(root, 'NEWS.md'))
  writeLines('Parentheses', file.path(root, 'doc(v1).md'))
  writeLines('A document', file.path(root, 'with space.md'))
  writeLines(sprintf('bibentry("Manual", title="Fixture", author="A. Tester", year="2026", url=%s)',
    encodeString(citation_url, quote = '"')), file.path(root, 'inst', 'CITATION'))
  if (!is.null(extra)) extra(root)
  archive <- file.path(base, 'fixture.tar.gz')
  withr::with_dir(base, utils::tar(archive, files = 'fixture', compression = 'gzip', tar = 'internal'))
  archive
}

test_that('valid local and external links pass', {
  archive <- fixture_archive(readme = c('[News](NEWS.md)', '[ref][n]', '', '[n]: NEWS.md',
    '[Space](with%20space.md)', '<a href="NEWS.md#release">News</a>',
    '[Remote](https://example.org)', '[Email](mailto:a@example.org)', '[Anchor](#section)'),
    url = 'https://example.org (project), https://example.org/manual, doc(v1).md',
    description = 'Methods <doi:10.1000/182>.')
  expect_output(check_archive(archive), 'PASS')
})

test_that('missing files and host-specific paths fail in Markdown and metadata', {
  targets <- c('docs/development-2026-10-18.md', 'docs/simulation-ici.md',
    'missing%20file.md', '../outside.md', '~/private.md', '/tmp/private.md',
    'file:///tmp/private.md', 'C:/Users/test/private.md',
    'C:%2FUsers/test/private.md', 'C:\\Users\\test\\private.md',
    '\\\\server\\share\\private.md')
  for (target in targets) {
    for (field in c('readme', 'url', 'bugs', 'citation_url')) {
      args <- setNames(list(target), field)
      if (field == 'readme') args[[1]] <- sprintf('<a href="%s">target</a>', target)
      archive <- do.call(fixture_archive, args)
      expect_error(suppressMessages(capture.output(check_archive(archive))),
                   'Invalid local file URIs', info = paste(field, target))
    }
  }
})

test_that('Description angle links, Rd and installed HTML are audited', {
  expect_error(capture.output(check_archive(fixture_archive(
    description = 'See <file:///tmp/private.md>.'))), 'Invalid local file URIs')
  for (format in c('rd', 'html', 'reference')) {
    extra <- function(root) {
      if (format == 'rd') {
        dir.create(file.path(root, 'man'))
        writeLines(c('\\name{fixture}', '\\alias{fixture}', '\\title{Fixture}',
          '\\description{See \\url{missing.md} and \\href{missing2.md}{here}.}'),
          file.path(root, 'man', 'fixture.Rd'))
      } else if (format == 'html') {
        dir.create(file.path(root, 'inst', 'doc'))
        writeLines('<html><body><img src="missing.png"></body></html>',
          file.path(root, 'inst', 'doc', 'fixture.html'))
      } else writeLines(c('[reference][r]', '', '[r]: missing.md'), file.path(root, 'README.md'))
    }
    expect_error(capture.output(check_archive(fixture_archive(extra = extra))),
                 'Invalid local file URIs', info = format)
  }
})

test_that('case mismatches are detected on case-sensitive filesystems', {
  probe <- tempfile(); dir.create(probe)
  withr::defer(unlink(probe, recursive = TRUE))
  writeLines('x', file.path(probe, 'UPPER'))
  expect_output(check_archive(fixture_archive()), 'PASS')
  if (!file.exists(file.path(probe, 'upper'))) {
    expect_error(capture.output(check_archive(fixture_archive(readme = '[News](news.md)'))),
                 'Invalid local file URIs')
  }
})
