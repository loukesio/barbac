#' Report which external command-line tools barbac can find
#'
#' The CLI pipeline shells out to FastQC, MultiQC, PEAR, minimap2 and samtools.
#' They may come from the conda environment \code{configure_environment()}
#' creates, or already be installed system-wide; either works, so both are
#' searched and the result says which one supplied each tool.
#'
#' @param env_name Conda environment to search first. Default: "barbac_env".
#' @param tools Character vector of tool names to look for.
#'
#' @return A data.frame with one row per tool: \code{tool}, \code{available},
#'   \code{source} (the environment name, \code{"PATH"}, or NA when missing),
#'   \code{version} (NA when the tool reports none), and \code{path}.
#' @export
check_barbac_tools <- function(env_name = "barbac_env",
                               tools = c("fastqc", "multiqc", "pear",
                                         "minimap2", "samtools")) {

  # Locate the environment once, if it exists at all. A missing environment is
  # not an error: the tools may still be installed system-wide, and reporting
  # them as absent would send users to reinstall software they already have.
  conda_prefix <- NULL
  conda_exe <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  if (!is.null(conda_exe) && nzchar(conda_exe)) {
    envs <- suppressWarnings(
      system2(conda_exe, c("env", "list"), stdout = TRUE, stderr = FALSE))
    hit <- grep(paste0("[/\\\\]", env_name, "[/\\\\]?\\s*$"), envs, value = TRUE)
    if (length(hit) > 0) {
      conda_prefix <- sub("^.*?\\s+(/.*?)\\s*$", "\\1", hit[[1]])
      if (!dir.exists(conda_prefix)) conda_prefix <- NULL
    }
  }

  # Not every tool implements --version: PEAR exits non-zero on it and prints
  # a usage error. Availability is therefore decided by the executable being
  # present, never by whether a version string could be parsed, so a tool is
  # never reported missing merely because it is quiet about its version.
  tool_version <- function(path) {
    for (flag in c("--version", "-v", "-h")) {
      out <- suppressWarnings(tryCatch(
        system2(path, flag, stdout = TRUE, stderr = TRUE),
        error = function(e) character(0)))
      out <- out[nzchar(out)]
      hit <- grep("[0-9]+\\.[0-9]", out, value = TRUE)
      if (length(hit) > 0) return(trimws(hit[[1]]))
    }
    NA_character_
  }

  rows <- lapply(tools, function(tool) {
    path <- ""
    where <- NA_character_

    if (!is.null(conda_prefix)) {
      candidate <- file.path(conda_prefix, "bin", tool)
      if (file.exists(candidate)) {
        path <- candidate
        where <- env_name
      }
    }
    if (!nzchar(path)) {
      candidate <- unname(Sys.which(tool))
      if (nzchar(candidate)) {
        path <- candidate
        where <- "PATH"
      }
    }

    data.frame(
      tool      = tool,
      available = nzchar(path),
      source    = where,
      version   = if (nzchar(path)) tool_version(path) else NA_character_,
      path      = if (nzchar(path)) path else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# Resolve an external tool to a full path.
#
# configure_environment() installs the CLI tools into a conda environment, which
# is not on PATH. Invoking them by bare name therefore fails with "command not
# found" even though the package installed them itself, so every call site
# resolves the tool first and runs it by path. The search order matches
# check_barbac_tools(), so what that function reports is what actually runs.
#
#' @keywords internal
#' @noRd
.barbac_env_prefix <- local({
  cached <- NULL
  resolved <- FALSE
  function(env_name = "barbac_env") {
    if (resolved) return(cached)
    resolved <<- TRUE
    conda_exe <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
    if (is.null(conda_exe) || !nzchar(conda_exe)) return(cached)
    envs <- suppressWarnings(
      system2(conda_exe, c("env", "list"), stdout = TRUE, stderr = FALSE))
    hit <- grep(paste0("[/\\\\]", env_name, "[/\\\\]?\\s*$"), envs, value = TRUE)
    if (length(hit) > 0) {
      p <- sub("^.*?\\s+(/.*?)\\s*$", "\\1", hit[[1]])
      if (dir.exists(p)) cached <<- p
    }
    cached
  }
})

#' @keywords internal
#' @noRd
.barbac_tool <- function(tool, env_name = "barbac_env") {
  prefix <- .barbac_env_prefix(env_name)
  if (!is.null(prefix)) {
    candidate <- file.path(prefix, "bin", tool)
    if (file.exists(candidate)) return(candidate)
  }
  unname(Sys.which(tool))
}

# Fail before doing any work when a required tool is missing, naming the tool
# and how to get it. A pipeline that starts without its tools reports a
# downstream symptom -- an empty output directory several steps later -- which
# is far harder to act on than the real cause.
#
#' @keywords internal
#' @noRd
.barbac_require_tools <- function(tools, env_name = "barbac_env") {
  paths <- vapply(tools, .barbac_tool, character(1), env_name = env_name)
  missing <- tools[!nzchar(paths)]
  if (length(missing) > 0) {
    stop("Required tool(s) not found: ", paste(missing, collapse = ", "),
         ".\n  Searched the '", env_name, "' conda environment and PATH.",
         "\n  Install them with configure_environment(), or put them on PATH.",
         "\n  check_barbac_tools() shows what barbac can currently find.",
         call. = FALSE)
  }
  paths
}
