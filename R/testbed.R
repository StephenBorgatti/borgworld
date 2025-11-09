# testbed

set.seed(42)
x <- seq(-3, 3, length.out = 100)  # Symmetric range around 0
y <- 5 + x^2 + rnorm(100, sd = 2)  # Inverted parabola with modest noise
df <- data.frame(x = x, y = y)

# Plot it
bscatterplot(df, y ~ poly(x, 2))

