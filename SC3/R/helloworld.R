# .libPaths()
# library()
# install.packages("")

library(SC3)

add <- function(x, y) {
    x + y
}

print(add(1, 2))
print(add(1.0e10, 2e1))
print(paste("hello", "world"))

h <- c(1, 2, 3, 4, 5, 6)
M <- c("a", "b", "c", "d", "e", "f")
barplot(h,
    names.arg = M, xlab = "X", ylab = "Y", col = "#FF0000",
    main = "Chart", border = "#0000FF"
) # can be save into ppt for independent editing
