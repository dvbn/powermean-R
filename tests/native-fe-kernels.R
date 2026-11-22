basis <- getFromNamespace("exact_basis", "powermean")
project <- getFromNamespace("project_fe", "powermean")
groups <- cbind(rep(0:3, each = 4), rep(0:3, times = 4))
levels <- c(4L, 4L)
selected <- basis(groups, levels)
stopifnot(selected$rank == 7L, selected$components == 1L,
          length(selected$columns) == 7L, !selected$multiprecision)
rhs <- cbind(sin(seq_len(16)), cos(seq_len(16)))
weights <- as.double(seq_len(16))
dummy <- cbind(model.matrix(~ 0 + factor(groups[, 1L])),
               model.matrix(~ 0 + factor(groups[, 2L]))[, -1L])
expected <- lm.wfit(dummy, rhs, weights)$residuals
answer <- project(groups, levels, rhs, weights)
stopifnot(answer$status == 0L, max(abs(answer$values - expected)) < 1e-10,
          max(abs(answer$group_residual)) < 1e-10)

weighted_sums <- getFromNamespace("C_pm_group_weighted_sums", "powermean")
aggregate <- .Call(weighted_sums, groups, levels, 2L, rhs, weights)
reference <- rowsum(rhs * weights, groups[, 2L])
stopifnot(identical(dim(aggregate), dim(reference)), max(abs(aggregate - reference)) < 1e-12)

even <- rbind(c(0L, 0L, 0L), c(0L, 1L, 1L), c(1L, 0L, 1L), c(1L, 1L, 0L))
rank_even <- basis(even, rep.int(2L, 3L))
stopifnot(identical(rank_even$columns, 0:3),
          identical(rank_even$position, c(0L, 1L, -1L, 2L, -1L, 3L)))
even[3:4, 3L] <- even[3:4, 2L]
rank_alias <- basis(even, rep.int(2L, 3L))
stopifnot(identical(rank_alias$columns, 0:2),
          identical(rank_alias$position, c(0L, 1L, -1L, 2L, -1L, -1L)))
