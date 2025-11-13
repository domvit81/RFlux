



mode <- function(x) {
	uniq <- unique(x)
	uniq[which.max(tabulate(match(x, uniq)))]
	}