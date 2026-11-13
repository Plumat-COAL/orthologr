# The native routines are registered by the code Rcpp::compileAttributes() writes into src/RcppExports.cpp (there is
# no hand-written init file). A registration that does not match the wrappers in R/RcppExports.R makes .Call()
# fail with "Incorrect number of arguments" in code that is not byte-compiled (R CMD INSTALL --no-byte-compile,
# devtools::load_all(), covr) - while a byte-compiled installation hides it.

test_that("every registered native routine has the number of arguments of its R wrapper", {
        reg <- getDLLRegisteredRoutines("orthologr")[[".Call"]]
        expect_gt(length(reg), 0)
        ns <- asNamespace("orthologr")
        for (routine in names(reg)) {
                wrapper <- sub("^_orthologr_", "", routine)
                expect_true(exists(wrapper, envir = ns, inherits = FALSE), info = routine)
                expect_identical(as.integer(reg[[routine]]$numParameters), length(formals(get(wrapper, envir = ns))),
                                 info = routine)
        }
})

# the first .Call(...) inside an expression, or NULL
find_dot_call <- function(e) {
        if (!is.call(e))
                return(NULL)
        if (identical(e[[1]], as.name(".Call")))
                return(e)
        for (part in as.list(e)[-1]) {
                if (missing(part))
                        next
                found <- find_dot_call(part)
                if (!is.null(found))
                        return(found)
        }
        NULL
}

test_that("every .Call() wrapper calls a registered routine with as many arguments as are registered", {
        reg <- getDLLRegisteredRoutines("orthologr")[[".Call"]]
        ns <- asNamespace("orthologr")
        n_checked <- 0
        for (name in ls(ns, all.names = TRUE)) {
                f <- get(name, envir = ns)
                if (!is.function(f) || is.primitive(f))
                        next
                b <- find_dot_call(body(f))
                if (is.null(b))
                        next
                routine <- as.character(b[[2]])
                expect_true(is.element(routine, names(reg)), info = name)
                # .Call(<routine>, arg1, ...): the call holds the routine and one element per argument
                expect_identical(length(b) - 2L, as.integer(reg[[routine]]$numParameters), info = name)
                n_checked <- n_checked + 1
        }
        expect_gte(n_checked, 12)
})

test_that("gestimator() runs from code that is not byte-compiled", {
        # a copy of the wrapper that the byte-code compiler has not seen: .Call() then checks the argument count
        f <- eval(parse(text = deparse(get("gestimator", envir = asNamespace("orthologr")))),
                  envir = asNamespace("orthologr"))
        out <- tempfile(fileext = ".tsv")
        on.exit(unlink(out), add = TRUE)
        expect_error(f(file = system.file("seqs/pal2nal.aln", package = "orthologr"), file_out = out), NA)
        expect_true(file.exists(out) && file.size(out) > 0)
})
