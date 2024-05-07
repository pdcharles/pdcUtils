#' @export
img <- function(filename, code, width = 7, height = 7, units = "in", res = 144, env = .GlobalEnv, ...) { # image writing fucntion, jupyter-compatible, render/spin-compatible
  is_jupyter <- "IRkernel" %in% loadedNamespaces()
  is_knitr <- !is.null(getOption("knitr.in.progress"))
  is_png <- grepl("\\.png$", filename)
  is_jpg <- grepl("\\.jpg$", filename) || grepl("\\.jpeg$", filename)
  to_prepend <- ""
  if (!grepl("/", filename)) {
    dir.create("out", showWarnings = F)
    to_prepend <- "out/"
  }
  if (is_jupyter || is_knitr || is_png) {
    fn <- paste0(to_prepend, if (is_png) filename else sub("\\.[^.]+$", ".png", filename))
    png(fn, width = width, height = height, units = units, res = res, ...)
    eval(substitute(code), envir = env)
    dev.off()
    if (is_jupyter) {
      IRdisplay::display_png(file = fn, width = width, height = height, units = units, res = res)
    } else {
      if (is_knitr) knitr::include_graphics(fn, dpi = NA)
    }
  } else if (is_jpg) {
    jpeg(paste0(to_prepend, filename), width = width, height = height, units = units, res = res, ...)
    eval(substitute(code), envir = env)
    dev.off()
  } else {
    if (units == "px") {
      width <- width / 72
      height <- height / 72
    }
    if (!grepl("\\..+$", filename)) {
      filename <- paste0(filename, ".pdf")
    } else {
      filename <- sub("\\.[^.]+$", ".pdf", filename)
    }
    cairo_pdf(paste0(to_prepend, filename), width = width, height = height, ...)
    eval(substitute(code), envir = env)
    dev.off()
  }
}

#' @export
assignPrecomputed <- function(variable, code, env = .GlobalEnv) { # create variable unless exists or file object exists
  dir.create("precomputed", showWarnings = F)
  if (!exists(variable)) {
    fileName <- paste0("precomputed/", variable, ".rds")
    if (!file.exists(fileName)) {
      obj <- eval(substitute(code), envir = env)
      saveRDS(obj, fileName)
      assign(variable, obj, envir = env)
    } else {
      assign(variable, readRDS(fileName), envir = env)
    }
  }
}

#' @export
mean_NAsafe <- function(x) { # NA-safe calculation of means
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = T)
}

#' @export
bp <- function(m, labels = NULL, ...) { # make boxplot from matrix/data.frame
  if (is.null(labels)) labels <- colnames(m)
  m[is.infinite(m)] <- NA
  m <- lapply(as.data.frame(m), na.omit)
  do.call(boxplot, list(m, names = labels, las = 2, pch = 20, ...))
}

#' @export
source <- (function() { # patched version of source, allows specification of a line to start/end on by grep matching
  source_txt <- deparse(source, control = "all")
  fix1 <- grep("function (file, ", source_txt, fixed = T)
  source_txt[fix1] <- sub(
    "function (file, ",
    "function (file, start = NULL, stop = NULL, ",
    source_txt[fix1],
    fixed = T
  )
  fix2 <- grep('.Internal(parse(stdin(), n = -1, lines, "?", ', source_txt, fixed = T)
  source_txt[fix2] <- sub(
    '.Internal(parse(stdin(), n = -1, lines, "?", ',
    '.Internal(parse(stdin(), n = -1, lines[ (function() {
      if (!is.null(start)) {
       from <- grep(start,lines,fixed=T)[1]
      } else from <- 1
      if (!is.null(stop)) {
       to <- (from-1) + grep(stop,lines[from:length(lines)],fixed=T)[1]
      } else to <- length(lines)
      from:to
     })() ], "?", ',
    source_txt[fix2],
    fixed = T
  )
  eval(parse(text = source_txt))
})()

#' @export
whatdidido <- function() { # what did I do?
  objs <- ls(.GlobalEnv)
  histFile <- tempfile("rawhist")
  origHistLen <- Sys.getenv("R_HISTSIZE")
  Sys.setenv("R_HISTSIZE" = "1000000000")
  savehistory(histFile)
  Sys.setenv("R_HISTSIZE" = origHistLen)
  histLines <- readLines(histFile)
  unlink(histFile)
  histExprs <- c()
  lineNumbers <- c()
  b <- 0
  while (b <= length(histLines)) {
    b <- b + 1
    if (!grepl("^library\\(|^require\\(|^data\\(", histLines[b]) && !any(sapply(objs, function(obj) grepl(paste0("^\\s*", obj, "\\s*<-\\s*"), histLines[b])))) next
    # print(histLines[b])
    # cat('starting relevant expr at line ',b,'\n')
    found <- F
    for (e in b:length(histLines)) {
      expr <- tryCatch(as.character(parse(text = histLines[b:e]))[1], error = function(e) e)
      if (!inherits(expr, "error")) {
        histExprs <- c(histExprs, expr)
        lineNumbers <- c(lineNumbers, b)
        # cat('expr closes at line ',e,'\n')
        found <- T
        b <- e
        break
      }
    }
    # if (!found) cat('no valid expr closure found\n')
  }
  indices <- c()
  libI <- grep("^library\\(|^require\\(|^data\\(", histExprs)
  if (length(libI)) {
    indices <- c(indices, sapply(unique(histExprs[libI]), function(cmd) libI[histExprs[libI] %in% cmd][1]))
  }
  if (length(objs)) {
    objIndices <- sapply(objs, function(obj) tail(grep(paste0("^\\s*", obj, "\\s*<-\\s*"), histExprs), n = 1))
    indices <- c(indices, unlist(objIndices[lapply(objIndices, length) > 0]))
  }
  structure(histExprs[indices][order(indices)], lineNumbers = lineNumbers[indices][order(indices)], class = "whatdidido")
}


#' @export
convertID <- function(database, acc, accType, retrieveAccType, removeTextNA = T) { # convert gene/protein ids to another format
  if (class(database) == "UniProt.ws") {
    retrieveAcc <- suppressMessages(select(database, acc[1], keytype = accType, columns = c(retrieveAccType)))
    if (length(acc) > 1) {
      prog <- txtProgressBar(style = 3)
      for (i in seq(2, length(acc), by = 100)) {
        setTxtProgressBar(prog, i / length(acc))
        retrieveAcc <- rbind(retrieveAcc, suppressMessages(select(database, acc[i:min(length(acc), i + 99)], keytype = accType, columns = c(retrieveAccType))))
      }
      setTxtProgressBar(prog, 1)
      close(prog)
    }
  } else {
    retrieveAcc <- select(database, acc, keytype = accType, columns = c(retrieveAccType))
  }
  retrieveAcc <- data.frame(b = tapply(retrieveAcc[[retrieveAccType]], retrieveAcc[[accType]], paste, collapse = ";"))
  retrieveAcc <- data.frame(a = acc, retrieveAcc[acc, , drop = F])
  rownames(retrieveAcc) <- NULL
  retrieveAcc[is.na(retrieveAcc[, "b"]), "b"] <- retrieveAcc[is.na(retrieveAcc[["b"]]), "a"]
  if (removeTextNA) retrieveAcc[retrieveAcc[, "b"] == "NA", "b"] <- retrieveAcc[retrieveAcc[, "b"] == "NA", "a"]
  colnames(retrieveAcc) <- c(accType, retrieveAccType)
  retrieveAcc
}

#' @export
wrap <- Vectorize(function(text, len) {
  words <- strsplit(text, " +")[[1]]
  ret <- words[1]
  if (length(words) > 1) {
    i <- 2
    j <- 1
    while (i <= length(words)) {
      ret[j] <- paste(ret[j], words[i])
      i <- i + 1
      if (nchar(ret[j]) >= len) {
        if (i > length(words)) break
        j <- j + 1
        ret[j] <- words[i]
        i <- i + 1
      }
    }
  }
  paste0(ret, collapse = "\n")
}, vectorize.args = "text")
