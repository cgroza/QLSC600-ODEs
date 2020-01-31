library(tidyverse)
logistic.map <- function(r){
  x_i <- vector('numeric', length=51)
  x_i[1] <- 0.1

  for(i in 2:51 ){
    x_i[i] <- r * x_i[i - 1] * (1 - x_i[i - 1])
  }
  x_i
}
curves <- lapply(c(0.5, 0.9, 2.8, 3.3), logistic.map)
names(curves) <- c("0.5", "0.9", "2.8", "3.3")

theme_set(theme_bw())
as_tibble(curves) %>% mutate(Generation = row_number()) %>% gather(key = r, value = PopSize, -Generation) %>%
  ggplot() + geom_line(aes(Generation, PopSize, color = r)) + ggtitle("Logistic map for values of r") +
  ylab("Population size")

forward.euler <- function(f.x, x.0, t.max, h){
  iter <- ceiling(t.max/h)
  x_i <- vector("numeric", length = iter)
  x_i[1] <- x.0
  t_i <- vector("numeric", length = iter)
  t_i[1] <- 0

  for(i in 2:(iter+1)) {
    x_i[i] <- x_i[i - 1] + h * f.x(x_i[i - 1])
    t_i[i] <- i * h
    }
  tibble(x=x_i, t=t_i, h=as.character(h))
}
solutions <- lapply(c(0.1, 0.01, 0.001), function(h) forward.euler(identity, 1, 6, h))
solutions <- do.call(rbind, solutions)

theme_set(theme_bw())
as_tibble() %>%
  ggplot() + geom_line(data = solutions, mapping = aes(x = t, y = x, color = h)) + ggtitle("Forward Euler") +
  geom_line(aes(x=seq(1,6, by = 0.001), y=exp(seq(1,6, by = 0.001)))) + geom_text(aes(x = 4.5, y = 200, label = "x(t) = exp(t)"))

forward.euler.2d <- function(dv.dt, dw.dt, v.0 = 1, w.0 = 0.1, I = 0.5, a = 0.7,
                             b = 0.8, epsilon = 0.08, t.max = 400, h = 0.01){
  iter <- ceiling(t.max/h)
  v_i <- vector("numeric", length = iter)
  v_i[1] <- v.0

  w_i <- vector("numeric", length = iter)
  w_i[1] <- w.0

  t_i <- vector("numeric", length = iter)
  t_i[1] <- 0

  for(i in 2:(iter+1)) {
    v_i[i] <- v_i[i - 1] + h * dv.dt(v_i[i - 1], w_i[i - 1], I)
    w_i[i] <- w_i[i - 1] + h * dw.dt(v_i[i - 1], w_i[i - 1], a, b, epsilon)
    t_i[i] <- i * h
  }
  tibble(v = v_i, w = w_i, t=t_i, h=as.character(h))
}
simulation = forward.euler.2d(function(v, w, I) {v - v^3/3 - w + I},
                 function(v, w, a, b, epsilon) {epsilon*(v + a - b*w)}

)

library(patchwork)
v.plot <- ggplot() + geom_line(data = simulation, mapping = aes(x = t, y = v))
w.plot <- ggplot() + geom_line(data = simulation, mapping = aes(x = t, y = w))
v.plot + w.plot + plot_annotation(tag_levels = 'a', title  = "FitzHugh-Nagumo") + plot_layout(nrow=2, ncol = 1)

v.coords <- seq(-2, 2, by = 0.01)
w.coords <- seq(-0.8, 1.5, by = 0.01)
ggplot() + geom_point(data = simulation, mapping = aes(x = v, y = w)) +
  geom_line(aes(x = v.coords , v.coords - v.coords ^ 3 /3 + 0.5), color = "red") +
  geom_line(aes(y = w.coords , x = 0.8 * w.coords - 0.7), color = "blue" ) + ggtitle("Phase space of the system")

eigen <- eigen(matrix(c(0.3522326, -1, 0.08, -0.08 * 0.8), nrow = 2, ncol = 2, byrow=T))
data.frame(Eigenvalue = eigen$values)

simulation.newI = forward.euler.2d(function(v, w, I) {v - v^3/3 - w + I},
                 function(v, w, a, b, epsilon) {epsilon*(v + a - b*w)}, I = 0)

library(patchwork)
v.plot.newI <- ggplot() + geom_line(data = simulation.newI, mapping = aes(x = t, y = v)) + ggtitle("Trajectory of v with t")
w.plot.newI <- ggplot() + geom_line(data = simulation.newI, mapping = aes(x = t, y = w)) + ggtitle("Trajectory of w with t")
phasespace.newI <- ggplot() + geom_point(data = simulation.newI, mapping = aes(x = v, y = w)) +
  geom_line(aes(x = v.coords , v.coords - v.coords ^ 3 /3), color = "red") +
  geom_line(aes(y = w.coords , x = 0.8 * w.coords - 0.7), color = "blue" ) + ggtitle("Phase space of the system (I = 0)")
 (v.plot.newI / w.plot.newI | phasespace.newI) + plot_annotation(tag_levels = 'a', title  = "FitzHugh-Nagumo with I = 0")

eigen <- eigen(matrix(c(-0.4385604, -1, 0.08, -0.08 * 0.8), nrow = 2, ncol = 2, byrow=T))
data.frame(Eigenvalue = eigen$values)

v.fixed_i <- vector("numeric", length = 501)
I_i <- vector("numeric", length = 501)
eigen.value_i <- vector("numeric", length = 501)
i <- 1
for(I in seq(0, 0.5, by = 0.001)){
  roots <- polyroot(c(I - 0.7/0.8, 1 - 1/0.8, 0, -1/3))
  v.fixed <- Re(roots[abs(Im(roots)) < 1e-10])

  v.fixed_i[i] <- v.fixed
  I_i[i] <- I

  eigen <- eigen(matrix(c(1 - v.fixed^2, -1, 0.08, -0.08 * 0.8),
                        nrow = 2, ncol = 2, byrow=T))
  eigen.value_i[i] <- Re(eigen$values)[1]
  i <- i + 1
}
stability <- tibble(v.fixed = v.fixed_i, I = I_i,
                    Stable = ifelse(eigen.value_i > 0, "Unstable", "Stable"))

stability %>% ggplot() + geom_line(aes(x=I, y=v.fixed, color = Stable)) + ggtitle("Stability of fixed points versus I") +
  ylab("v*") + theme(legend.position = "top")
