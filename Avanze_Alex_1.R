# Establecer el directorio de trabajo
setwd("C:/Users/Acer/Documents/carpeta_zi")

# Listar los archivos en el directorio
files <- list.files()
print(files)


trajan <- as.data.frame(readRDS('trajan_recoded.rds'))
View(trajan)
#Estadistica descriptivas de trajan

# Rango
range(trajan$nshoots, na.rm = TRUE)
range(trajan$hormone, na.rm = TRUE)
range(trajan$period, na.rm = TRUE)
range(trajan$period2, na.rm = TRUE)
range(trajan$hormone2, na.rm = TRUE)

# Coeficiente de Variación
cv <- function(x) { sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) }
cv(trajan$nshoots)
cv(trajan$hormone)
cv(trajan$period)
cv(trajan$period2)
cv(trajan$hormone2)

#Distribución y Frecuencia
# Tablas de frecuencia
table(trajan$nshoots)
table(trajan$hormone)
table(trajan$period)
table(trajan$period2)
table(trajan$hormone2)

#Visualización de Datos:
# Histograma
hist(trajan$nshoots, main = "Histograma de nshoots", xlab = "nshoots")
hist(trajan$hormone, main = "Histograma de hormone", xlab = "hormone")
# Crear el histograma
# Crear el histograma con plotly
histograma <- plot_ly(
  data = trajan, 
  x = ~nshoots, 
  type = 'histogram',
  marker = list(color = 'rgba(100, 200, 102, 0.7)', line = list(color = 'rgba(0, 0, 0, 1.0)', width = 1))
) %>%
  layout(
    title = list(text = "Histograma de nshoots", font = list(size = 24, color = "darkblue")),
    xaxis = list(title = "Número de shoots", titlefont = list(size = 18, color = "darkblue")),
    yaxis = list(title = "Frecuencia", titlefont = list(size = 18, color = "darkblue")),
    bargap = 0.2
)

histograma
# Boxplot
boxplot(trajan$nshoots, main = "Boxplot de nshoots", ylab = "nshoots")
boxplot(trajan$hormone, main = "Boxplot de hormone", ylab = "hormone")

# Gráficos de dispersión
plot(trajan$period, trajan$nshoots, main = "Scatterplot de period vs nshoots", xlab = "period", ylab = "nshoots")
plot(trajan$hormone, trajan$nshoots, main = "Scatterplot de hormone vs nshoots", xlab = "hormone", ylab = "nshoots")

  
###################################################################
source('ZI.R')

# Standard Poisson
trajan_po <- glm(nshoots ~ 0 + period3 : hormone3, data = trajan, family = poisson)
# NB-lin (called NBII by gamlss)
trajan_nblin <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBII)
# NB-quad (called NBI by gamlss)
trajan_nbquad <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBI)
# ZI types A to D
trajan_po_typeA <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "A", family = "poisson")