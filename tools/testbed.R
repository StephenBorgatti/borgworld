# testbed
#' @exclude


set.seed(42)
x <- seq(-3, 3, length.out = 100)  # Symmetric range around 0
y <- 5 + x^2 + rnorm(100, sd = 2)  # Inverted parabola with modest noise
df <- data.frame(x = x, y = y)

# Plot it
bscatterplot(df, y ~ poly(x, 2))

# Split with names
by_gender <- nhanes |> split(~female)
# regress
for (name in names(by_gender)) {
  cat("Female = ",name)
  by_gender[[name]] |> bregress(leadership ~ homophily)
}


# Create empty list for results
models <- list()

# Loop through each named element
for (name in names(by_gender)) {
  models[[name]] <- by_gender[[name]] |> bregress(leadership ~ homophily)
}

# Now you have models[["male"]] and models[["female"]]

nhanes$homsq = nhanes$homophily^2
bregress(nhanes,leadership ~ homophily + homsq)

center = function(x) {
  (x - mean(x, na.rm=TRUE))/sd(x,na.rm=TRUE)
}
nhanes$chomophily = center(nhanes$homophily)
nhanes$chomsq = nhanes$chomophily^2
bregress(nhanes,leadership ~ chomophily + chomsq)

mod = bregress(nhanes,leadership ~ homophily + homsq)
fit = mod$fitted.values
plot(y=fit,x=nhanes$homophily)

mod = bregress(nhanes,leadership ~ chomophily + chomsq)
fit = mod$fitted.values
plot(y=fit,x=nhanes$chomophily)

nhanes = nhanes |>
  mutate(highbp = 1*(systolic_bp >= 140))

blogisticinteraction(nhanes,"highbp","age","weight")
blogisticinteraction(nhanes,y="highbp",x="age",w="weight")
