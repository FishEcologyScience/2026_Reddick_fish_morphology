## --------------------------------------------------------------#
## Script name: fct01-02_git_functions.R
##
## Purpose of script:
##    Git-aware object comparison utilities. Retrieve R objects
##    from previous commits and compare them to the current session.
##
## Dependencies:
##    - callr  (run R subprocess)
##    - waldo  (human-readable object diffs)
##    - git    (must be on system PATH)
##
## Author: Paul Bzonek
##
## Date Created: 2026-05-05
##
## --------------------------------------------------------------#
## Modification Notes:
##
## --------------------------------------------------------------#
##
## Known limitations:
##   - 01_data/04_gitignore_files/ is excluded from git. Scripts
##     that depend on data in that folder will fail in the subprocess.
##   - Git LFS files require LFS to be installed and authenticated.
##   - All-objects mode runs the full script chain; runtime reflects
##     the full pipeline.
##
## --------------------------------------------------------------#



#-------------------------------------------------------------#
#####fct_git_get_object ##################################----
#-------------------------------------------------------------#
#' Run scripts at a past commit and return R objects
#'
#' @description
#' Creates a temporary git worktree at the specified commit, sources one or
#' more R scripts in an isolated \code{callr} subprocess, and returns the
#' resulting R objects. The worktree is always removed on exit, even if an
#' error occurs.
#'
#' @param commit Character string. Any git ref: commit SHA, \code{"HEAD~1"},
#'   branch name, or tag.
#' @param scripts Character vector of script paths to source in order.
#'   Paths are relative to the project root. Defaults to
#'   \code{"02_scripts/script00-00_user_interface.R"}.
#' @param object_name Character string or \code{NULL}. If a string, return only
#'   that named object. If \code{NULL} (default), return all objects in the
#'   subprocess environment as a named list.
#'
#' @return The named object, or a named list of all objects created by the scripts.
#'
#' @examples
#' # Return a single object from the previous commit
#' # fct_git_get_object("HEAD~1", object_name = "df_walleye")
#'
#' # Return full environment from a specific commit
#' # fct_git_get_object("abc1234")
#'
#' @export

fct_git_get_object <- function(
  commit,
  scripts     = "02_scripts/script00-00_user_interface.R",
  object_name = NULL
) {

  #----------------------------
  # 1. Validate inputs
  #----------------------------
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required. Install with install.packages('callr').",
         call. = FALSE)
  }

  if (!is.character(scripts) || length(scripts) == 0) {
    stop("`scripts` must be a non-empty character vector.", call. = FALSE)
  }

  #----------------------------
  # 2. Create temp worktree
  #----------------------------
  temp_path <- tempfile(pattern = "rcompare_worktree_")

  # Guarantee cleanup even if an error occurs mid-function
  on.exit({
    if (dir.exists(temp_path)) {
      system2("git", c("worktree", "remove", "--force", temp_path),
              stdout = FALSE, stderr = FALSE)
    }
  }, add = TRUE)

  wt_out <- system2(
    "git",
    c("worktree", "add", temp_path, commit),
    stdout = TRUE, stderr = TRUE
  )

  if (!dir.exists(temp_path)) {
    stop(sprintf(
      "Could not create worktree for commit '%s':\n%s",
      commit, paste(wt_out, collapse = "\n")
    ), call. = FALSE)
  }

  #----------------------------
  # 3. Source scripts in subprocess and return objects
  #----------------------------
  tryCatch(
    callr::r(
      function(worktree_path, scripts, object_name) {
        setwd(worktree_path)
        for (script in scripts) {
          source(script)
        }
        if (is.null(object_name)) {
          as.list(.GlobalEnv)
        } else {
          if (!exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
            stop(sprintf(
              "Object '%s' not found in the environment after sourcing scripts.",
              object_name
            ))
          }
          get(object_name, envir = .GlobalEnv, inherits = FALSE)
        }
      },
      args = list(
        worktree_path = temp_path,
        scripts       = scripts,
        object_name   = object_name
      )
    ),
    error = function(e) {
      stop(sprintf(
        "Error running scripts at commit '%s':\n%s",
        commit, conditionMessage(e)
      ), call. = FALSE)
    }
  )
}



#-------------------------------------------------------------#
#####fct_compare_object ##################################----
#-------------------------------------------------------------#
#' Compare R objects between the current session and a past commit
#'
#' @description
#' Compares one or all R objects in the current session against the same
#' objects produced by sourcing scripts at a specified git commit.
#'
#' Uses \code{\link[base]{identical}()} for exact equality. If objects differ,
#' automatically escalates to \code{\link[waldo]{compare}()} for a human-readable
#' diff showing exactly what changed.
#'
#' @param object_name Character string or \code{NULL}. If a string, compare
#'   only that object (\strong{single-object mode}). If \code{NULL} (default),
#'   compare all objects present in the current session and/or the old commit
#'   (\strong{all-objects mode}).
#' @param commit Character string. Git ref for the commit to compare against.
#'   Defaults to \code{"HEAD~1"} (the previous commit).
#' @param scripts Character vector of script paths to source at the old commit,
#'   in order. Paths are relative to the project root. Defaults to
#'   \code{"02_scripts/script00-00_user_interface.R"}.
#' @param object_name_old Character string or \code{NULL}. Override the object
#'   name on the old-commit side (useful if the object was renamed). Only
#'   applies in single-object mode; a warning is issued if supplied in
#'   all-objects mode. Defaults to \code{NULL} (same name as \code{object_name}).
#'
#' @return
#' \strong{Single-object mode:} \code{TRUE} invisibly if identical,
#' \code{FALSE} invisibly if not (waldo diff is printed).
#'
#' \strong{All-objects mode:} a data frame with columns \code{object},
#' \code{identical} (logical or \code{NA}), and \code{note}, returned
#' invisibly. Waldo diffs are printed for any differing objects.
#'
#' @examples
#' # Compare df_walleye in current session vs previous commit
#' # fct_compare_object("df_walleye")
#'
#' # Compare all objects vs two commits back
#' # fct_compare_object(commit = "HEAD~2")
#'
#' # Compare using partial pipeline (load + import only)
#' # fct_compare_object(
#' #   "df_walleye",
#' #   scripts = c("script00-01_load_packages.R", "script01-01_data_import.R")
#' # )
#'
#' # Object was renamed between commits
#' # fct_compare_object("df_fish", object_name_old = "df_walleye")
#'
#' @export

fct_compare_object <- function(
  object_name     = NULL,
  commit          = "HEAD~1",
  scripts         = "02_scripts/script00-00_user_interface.R",
  object_name_old = NULL
) {

  if (!requireNamespace("waldo", quietly = TRUE)) {
    stop("Package 'waldo' is required. Install with install.packages('waldo').",
         call. = FALSE)
  }

  #----------------------------
  # Single-object mode
  #----------------------------
  if (!is.null(object_name)) {

    old_name <- if (!is.null(object_name_old)) object_name_old else object_name

    if (!exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
      stop(sprintf("Object '%s' not found in the current session.", object_name),
           call. = FALSE)
    }
    obj_new <- get(object_name, envir = .GlobalEnv, inherits = FALSE)

    cat(sprintf("Retrieving '%s' from commit %s ...\n", old_name, commit))
    obj_old <- fct_git_get_object(commit = commit, scripts = scripts,
                                   object_name = old_name)

    if (identical(obj_new, obj_old)) {
      cat(sprintf("[MATCH] '%s' is identical between the current session and %s.\n",
                  object_name, commit))
      return(invisible(TRUE))
    }

    cat(sprintf("[DIFFERS] '%s' has changed since %s.\n\n", object_name, commit))
    print(waldo::compare(obj_new, obj_old, x_arg = "current", y_arg = commit))
    return(invisible(FALSE))
  }

  #----------------------------
  # All-objects mode
  #----------------------------
  if (!is.null(object_name_old)) {
    warning("`object_name_old` is ignored in all-objects mode.")
  }

  cat(sprintf("Retrieving full environment from commit %s ...\n", commit))
  old_env  <- fct_git_get_object(commit = commit, scripts = scripts,
                                  object_name = NULL)

  current_names <- ls(envir = .GlobalEnv)
  old_names     <- names(old_env)
  all_names     <- union(current_names, old_names)

  results <- lapply(all_names, function(nm) {
    in_current <- nm %in% current_names
    in_old     <- nm %in% old_names

    if (in_current && in_old) {
      is_identical <- identical(
        get(nm, envir = .GlobalEnv, inherits = FALSE),
        old_env[[nm]]
      )
      note <- if (is_identical) "" else "differs"
    } else if (in_current) {
      is_identical <- NA
      note         <- "only in current session"
    } else {
      is_identical <- NA
      note         <- sprintf("only in commit %s", commit)
    }

    data.frame(object = nm, identical = is_identical, note = note,
               stringsAsFactors = FALSE)
  })

  summary_tbl <- do.call(rbind, results)
  summary_tbl <- summary_tbl[order(
    is.na(summary_tbl$identical),  # non-NA first
    summary_tbl$identical,          # differs (FALSE) before match (TRUE)
    summary_tbl$object              # alphabetical within groups
  ), ]
  row.names(summary_tbl) <- NULL

  cat(sprintf("\n--- Object comparison: current session vs %s ---\n", commit))
  print(summary_tbl, row.names = FALSE)

  differ <- summary_tbl$object[!is.na(summary_tbl$identical) & !summary_tbl$identical]
  if (length(differ) > 0) {
    cat("\n--- waldo diffs for differing objects ---\n")
    for (nm in differ) {
      cat(sprintf("\n[ %s ]\n", nm))
      print(waldo::compare(
        get(nm, envir = .GlobalEnv, inherits = FALSE),
        old_env[[nm]],
        x_arg = "current",
        y_arg = commit
      ))
    }
  }

  invisible(summary_tbl)
}
