# Establecer el directorio de trabajo
setwd("C:/Users/Acer/Documents/carpeta_zi")

# Listar los archivos en el directorio
files <- list.files()
print(files)


trajan <- as.data.frame(readRDS('trajan_recoded.rds'))
View(trajan)
#Estadistica descriptivas de trajan


###################################################################
source('ZI.R')

# Standard Poisson
trajan_po <- glm(nshoots ~ 0 + period3 : hormone3, data = trajan, family = poisson)
# NB-lin (called NBII by gamlss)
trajan_nblin <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBII)

summarise()

# NB-quad (called NBI by gamlss)
trajan_nbquad <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBI)
# ZI types A to D
trajan_po_typeA <- ZI(nshoots ~ 0 + period3 : hormone3, data = trajan, ZI_type = "A", family = "poisson")

trejan_experiment <- gamlss(nshoots ~ 0 + period3 : hormone3, data = trajan, family = NBII)

############################################################################################################################################
# Cargar los paquetes
library(pscl)
library(MASS)
# Verifica si trajan ya es un data.frame o tibble
class(trajan)

# Verifica los nombres de las columnas en trajan
names(trajan)

# Ajustar el modelo Zero-Inflated Negative Binomial
model <- zeroinfl(nshoots ~ period3 + hormone3 | period3 + hormone3, 
                  data = trajan, 
                  dist = "negbin", 
                  link = "log")

# Ver resultados del modelo Zero-Inflated Negative Binomial
summary(model)

############################################################################################################################################
library(pscl)

# Ajustar un modelo inicial con un conjunto amplio de variables
initial_model <- zeroinfl(nshoots ~ period3 * hormone3 + period3:hormone3 | period3 + hormone3, data = trajan, dist = "negbin", link = "log")

# Usar `step` para ajustar el modelo y seleccionar las mejores variables basadas en AIC
best_model <- step(initial_model, direction = "both", trace = FALSE, k = 2)

# Ver el mejor modelo basado en AIC
summary(best_model)


#################################################################################################################

# Ejemplo de un modelo predictivo ajustado con variables clave
best_predictive_model <- zeroinfl(nshoots ~ period3 * hormone3 + period3:hormone3 | period3 + hormone3,
                                  data = trajan, dist = "negbin", link = "log")

# Ver resumen del mejor modelo predictivo
summary(best_predictive_model)


#####

# Ajustar el mejor modelo predictivo con variables clave
best_predictive_model <- zeroinfl(nshoots ~ period3 * hormone3 + period3:hormone3 | period3 + hormone3,
                                  data = trajan, dist = "negbin", link = "log")

# Ver resumen del mejor modelo predictivo
summary(best_predictive_model)


####

# Ajustar un modelo completo
full_model <- zeroinfl(nshoots ~ period3 * hormone3 + period3:hormone3 | period3 + hormone3,
                       data = trajan, dist = "negbin", link = "log")

# Realizar la selección hacia atrás basada en AIC
backward_model <- step(full_model, direction = "backward")
summary(backward_model)










# Crear una función para ajustar el modelo y calcular AIC
evaluate_model <- function(formula_count, formula_zero, data) {
  model <- zeroinfl(formula_count | formula_zero, data = data, dist = "negbin", link = "log")
  aic <- AIC(model)
  return(list(model = model, aic = aic))
}

# Generar diferentes combinaciones de variables
library(MASS)
combinations <- list(
  list(count_formula = nshoots ~ period3 + hormone3, zero_formula = period3 + hormone3),
  list(count_formula = nshoots ~ period3 * hormone3, zero_formula = period3 + hormone3),
  list(count_formula = nshoots ~ period3 + hormone3 + period3:hormone3, zero_formula = period3 + hormone3),
  list(count_formula = nshoots ~ period3 + hormone3 + period3:hormone3 + hormone3:period3, zero_formula = period3 + hormone3)
)

# Ajustar los modelos y almacenar los resultados
results <- list()
for (i in 1:length(combinations)) {
  formulas <- combinations[[i]]
  count_formula <- formulas$count_formula
  zero_formula <- formulas$zero_formula
  model_results <- evaluate_model(count_formula, zero_formula, trajan)
  results[[i]] <- list(count_formula = count_formula, zero_formula = zero_formula, model = model_results$model, aic = model_results$aic)
}

# Seleccionar el mejor modelo basado en AIC
best_model_index <- which.min(sapply(results, function(x) x$aic))
best_model <- results[[best_model_index]]$model
summary(best_model)



# Ajustar un modelo reducido eliminando variables no significativas
reduced_model <- zeroinfl(nshoots ~ period3 + hormone3 + period3:hormone3 | period3, 
                          data = trajan, dist = "negbin", link = "log")
summary(reduced_model)









library(caret)



# Comparar AIC del modelo reducido con el modelo completo
full_model_aic <- AIC(full_model)
reduced_model_aic <- AIC(reduced_model)

# Validación cruzada del modelo reducido
library(caret)
train_control <- trainControl(method = "cv", number = 10)
cv_model_reduced <- train(nshoots ~ period3 + hormone3 + period3:hormone3,
                          data = trajan, method = "zeroinfl", trControl = train_control, dist = "negbin", link = "log")

print(cv_model_reduced)













library(pscl)

# Definir la función para calcular AIC usando validación cruzada
cv_aic <- function(data, formula, n_splits = 10) {
  set.seed(123)  # Para reproducibilidad
  folds <- createFolds(data$nshoots, k = n_splits, list = TRUE, returnTrain = TRUE)
  
  aic_values <- c()
  
  for (i in 1:length(folds)) {
    train_indices <- folds[[i]]
    test_indices <- setdiff(1:nrow(data), train_indices)
    
    train_data <- data[train_indices, ]
    test_data <- data[test_indices, ]
    
    model <- zeroinfl(formula, data = train_data, dist = "negbin", link = "log")
    test_model <- zeroinfl(formula, data = test_data, dist = "negbin", link = "log")
    
    aic_values <- c(aic_values, AIC(test_model))
  }
  
  return(mean(aic_values))
}

# Fórmulas del modelo completo y reducido
full_formula <- nshoots ~ period3 * hormone3 + period3:hormone3
reduced_formula <- nshoots ~ period3 + hormone3 + period3:hormone3

# Calcular AIC promedio usando validación cruzada
full_model_aic_cv <- cv_aic(trajan, full_formula)
reduced_model_aic_cv <- cv_aic(trajan, reduced_formula)

print(paste("AIC promedio del modelo completo:", full_model_aic_cv))
print(paste("AIC promedio del modelo reducido:", reduced_model_aic_cv))






library(pscl)
library(caret)

# Función para calcular AIC con validación cruzada
cv_aic <- function(data, formula, n_splits = 10) {
  set.seed(123)  # Para reproducibilidad
  folds <- createFolds(data$nshoots, k = n_splits, list = TRUE, returnTrain = FALSE)
  
  aic_values <- c()
  
  for (i in 1:length(folds)) {
    train_indices <- setdiff(1:nrow(data), folds[[i]])
    train_data <- data[train_indices, ]
    test_data <- data[folds[[i]], ]
    
    model <- zeroinfl(formula, data = train_data, dist = "negbin", link = "log")
    
    # Predicción en el conjunto de prueba
    pred <- predict(model, newdata = test_data, type = "response")
    
    # Calcular el AIC en el conjunto de prueba
    loglik <- sum(dnbinom(test_data$nshoots, size = model$theta, mu = pred, log = TRUE))
    aic <- -2 * loglik + 2 * length(coef(model))
    aic_values <- c(aic_values, aic)
  }
  
  return(mean(aic_values))
}

# Fórmulas del modelo completo y reducido
full_formula <- nshoots ~ period3 * hormone3 + period3:hormone3
reduced_formula <- nshoots ~ period3 + hormone3 + period3:hormone3

# Calcular AIC promedio usando validación cruzada
full_model_aic_cv <- cv_aic(trajan, full_formula)
reduced_model_aic_cv <- cv_aic(trajan, reduced_formula)

print(paste("AIC promedio del modelo completo:", full_model_aic_cv))
print(paste("AIC promedio del modelo reducido:", reduced_model_aic_cv))

