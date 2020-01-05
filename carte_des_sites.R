# Make a map of the three goudoux
# Sources: https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html and https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html

# load librairies
install.packages("ggspatial");install.packages("rnaturalearth")
library("ggplot2"); theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("ggspatial")
library("RColorBrewer") # coloblnd friendly colors

# Load sites coordinates
sites <- data.frame(name = c('Grande-Plée-Bleue','Grande-Tourbière-de-Villeroy','Mer Bleue Peatland'),
                    lat = c(46.773268,46.380574,45.409408),
                    lon = c(-71.066612,-71.832705,-75.517953))

# Load Canada and USA information
CaUS <- ne_countries(country = c('United states of America', 'Canada'), scale = "medium", returnclass = "sf")

# Plot it
ggplot(data = CaUS) +
  geom_sf(color = "darkgrey", fill = "white") + # Colors of the borders and inland
  coord_sf(xlim = c(-78, -64.5), ylim = c(42, 50), expand = FALSE) + # Zoom in or out
  annotation_scale(location = "br", width_hint = 0.2) + # Scale at the br (bottom right) and size
  annotation_north_arrow(location = "br", which_north = "true",  # Arrow at the br and sizes
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = sites, mapping = aes(x = lon, y = lat), color = "#FB6A4A", size = 3) + # Coordinates of the plots
  geom_text(data = sites, aes(x= lon, y = lat, label=name),hjust=-0.1, vjust=0.1, size = 3.5, fontface = "bold") + # Plot the labels near the coordinates
  annotate("text", -Inf, Inf, label = "Canada", hjust = -2, vjust = 6, fontface = "bold", color = "darkgrey") + # Add Canada on the top left
  annotate("text", -Inf, -Inf, label = "U.S.A.", hjust = -4, vjust = -4, fontface = "bold", color = "darkgrey") + # Add USA on the bottom left
  labs(x = 'Longitude (°)', y = 'Latitude (°)') + # Label x and y axis
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = 'aliceblue')) # No grids, water in blue

# /!\ Note the warning of the inaccurate scale bar: 
# since the map use unprojected data in longitude/latitude (WGS84) on an equidistant cylindrical projection (all meridians being parallel),
# length in (kilo)meters on the map directly depends mathematically on the degree of latitude. Plots of small regions or projected data will often allow for more accurate scale bars.


#ggsave("map.pdf", width = 6, height = 6, )
ggsave("~/Desktop/carte_sites.png", width = 6, height = 6, dpi = 300)


