#!/usr/bin/env Rscript

input <- file('stdin', 'r')
lines <- readLines(input)

mat <- as.matrix(read.table(textConnection(lines) , comment.char = "#", header = T))
alphabet <- rownames(mat)
alphabet <- sort(alphabet)
if ("X" %in% alphabet) {
  alphabet <- c(alphabet[!(alphabet %in% c("B", "Z", "X"))], "X")
  mat <- mat[alphabet, alphabet]
  R <- 2**(0.5*mat)
  A <- nrow(R)
  p <- solve(R[-A,-A]) %*% rep(1, A-1)
  avgSeqId <- sum(diag(R[-A,-A])*p*p)
} else {
  alphabet <- alphabet[!(alphabet %in% c("B", "Z"))]
  mat <- mat[alphabet, alphabet]
  R <- 2**(0.5*mat)
  A <- nrow(R)
  p <- solve(R) %*% rep(1, A)
  avgSeqId <- sum(diag(R)*p*p)
}

header <- paste0(lines[startsWith(lines, "#")], collapse = "\n")
p0 <- paste0("# Background (precomputed optional): ", paste(unlist(as.list(round(p, 4))), collapse = " "), " 0.00001")
# FIXME: find out how to compute lambda
lambda <-    "# Lambda     (precomputed optional): 0.34657"
avg <-paste0("# Avg SeqId  (precomputed optional): ", round(avgSeqId, 5))
out <- paste0(capture.output(write.table(mat, quote = F, sep = "\t")), collapse = "\n")
cat(paste0(header, "\n", p0, "\n", lambda, "\n", avg, "\n\t", out, "\n"))

