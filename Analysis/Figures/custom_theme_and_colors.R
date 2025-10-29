chosen_colors = c("#6E8B3D", "#8B1A1A", "darkgoldenrod3", "cadetblue4")

chosen_colors_extended = c(
  # A2C, A2G, A2T,
  "darkgoldenrod3", "#6E8B3D",  "mediumpurple4",
  # C2A, C2G, C2T
  "cadetblue4", "rosybrown4", "#8B1A1A")

index_colors = c("darkgoldenrod3", "cadetblue4")

case_control_colors = c("#6E8B3D", "#8B1A1A")

theme_custom <- function(title_size = 14, text_size = 14, legend_position = "bottom") {
  
  # Start with theme_minimal as base
  theme_bw() %+replace%
    
    theme(
      # Text elements
      text = element_text(size = text_size),
      plot.title = element_text(
        hjust = 0.5, 
        vjust = -0.5,
        face = "bold",
        size = title_size,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, 
        vjust = -0.5,
        size = text_size,
        margin = margin(b = 10)
      ),
      
      # # Axis formatting
      # axis.title = element_text(size = text_size),
      # axis.text = element_text(size = text_size * 0.9),
      # axis.title.x = element_text(margin = margin(t = 10)),
      # axis.title.y = element_text(margin = margin(r = 10)),
      # 
      # Legend formatting
      legend.position = legend_position,
      # legend.title = element_text(size = text_size),
      # legend.text = element_text(size = text_size * 0.9)
    )
}

theme_custom_void <- function(title_size = 14, text_size = 14, legend_position = "bottom") {
  
  # Start with theme_minimal as base
  theme_void() %+replace%
    
    theme(
      # white background
      plot.background = element_rect(fill = "white", color = NA),
      
      # Text elements
      text = element_text(size = text_size),
      plot.title = element_text(
        hjust = 0.5, 
        vjust = -0.5,
        face = "bold",
        size = title_size,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, 
        vjust = -0.5,
        size = text_size,
        margin = margin(b = 10)
      ),
      
      # Legend formatting
      legend.position = legend_position,
    )
}
