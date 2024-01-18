
e_data <- data.frame()

for (i in 900:920) {
  p <- plot_current(res[[i]]$h, 155)
  es <- p[[1]][[1]]
  vert <- p[[1]][[2]]

  nd <- data.frame(frame = i, x = es$x, y = es$y, group = es$group)

  e_data <- rbind(e_data, nd)
}

fig <- e_data %>%
  plot_ly(
    x = ~x,
    y = ~y,
    frame = ~frame,
    type = 'scatter',
    mode = 'lines',
    transforms = list(
      list(
        type = 'groupby',
        groups = e_data$group
      )
    ),
    showlegend = F
  )

fig

ggplot(data= e_data)+
  geom_path(aes(x=x,y=y,group = group)) +
  transition_states(frame, transition_length = 2,
                    state_length = 1)

ggplot(data= vert)+
  geom_point(aes(x=x,y=y,group = group), size = 0.1) +
  coord_fixed()
