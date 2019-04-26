setwd('C:/Users/1035419453/Documents/DUVAN TORRES/Seguimiento a Rentas/Info Repositorio')

library(forecast)
library(lmtest)
library(xtable)
library(TSA)
library(FitAR)

#MEDIDAS

#----------------funcion medidas especial 
#                para loess y holt-winters

medidas.1 = function(yest,y,k){
  T = length(y)
  sse = sum((yest-y)^2)
  ssr = sum((y-mean(y))^2) 
  mse = sse/(T-k)
  R2 = 1 - sse/ssr
  Ra2 = 1 - (T-1)*(1-R2)/(T-k)
  aic = log((T-k)*exp(2*k/T)*mse/T)
  bic = log(T^(k/T)*(T-k)*mse/T)
  
  M = c(sqrt(mse),Ra2,  aic, bic)
  names(M) = c("mse","R2-ad","log.aic","log.bic")
  return(M)
}
#----------------funcion medidas
# source("medidas.r")

medidas = function(m,y,k){
  # m = objeto producido con lm()
  # y = variable dependiente
  # k = número de coeficientes beta
  T = length(y)
  yest = fitted(m)
  sse = sum((yest-y)^2)
  ssr = sum((y-mean(y))^2) 
  mse = sse/(T-k)
  R2 = 1 - sse/ssr
  Ra2 = 1 - (T-1)*(1-R2)/(T-k)
  aic = log((T-k)*exp(2*k/T)*mse/T)
  bic = log(T^(k/T)*(T-k)*mse/T)
  M = c(sqrt(mse),Ra2,  aic, bic)
  names(M) = c("rmse","R2-ad","log.aic","log.bic")
  return(M)
}

Archivo <- "https://raw.githubusercontent.com/duvantorres88/Work/master/Rentas_marzo.txt?token=AJMN2GYOHCB5HXUFBFGYUWS4YM2UM"

#Leyendo el archivo fuente. Se actualiza cada mes con la ejecución mensual

A <- read.table(file = Archivo, header = T, sep = "\t")

#Generador de fechas. desde enero de 2005 hasta la última ejecución mensual
fechas  =  seq(as.Date("2005/1/1"),
               length.out = dim(A)[1],  by  =  "month")

#Dataframe con rentas asignadas.

Rentas <- data.frame(Fecha=fechas, Predial=A$Impuesto.Predial...Impuesto.Predial.D.E.,
                     Rec.Predial= A$Recuperación.Cartera.Impuesto.Predial, 
                     IndyComercio= A$Impuesto.Industria.y.Comercio, 
                     Rec.IndyComercio=A$Recuperación.Cartera.Ind..Y.Ccio,
                     AyT=A$Avisos.y.tableros, Telefonos=A$Impuesto.de.teléfonos,
                     SobretasaGas=A$Sobretasa.a.la.gasolina...Sobretasa.a.la.gasolina.D.E.,
                     Espectaculos=A$Espectáculos.públicos, Vallas=A$Registro.de.Vallas,
                     Alumbrado=A$Alumbrado.público, 
                     IntMoraPredial=A$Intereses.de.mora.predial,
                     IntMoraIyC=A$Intereses.de.mora.Ind..Y.Ccio., 
                     IntMoraEspec=A$Intereses.de.mora.espectáculos.públicos,
                     IntMoraAyT=A$Intereses.de.mora.avisos.y.tableros,
                     IntMoraTelefonos=A$Intereses.de.mora.impuesto.de.teléfono,
                     IntMoraAlumbr=A$Intereses.de.mora.alumbrado.público,
                     IntMoraValla=A$Intereses.de.Mora.Registro.de.Vallas..Ley.140.94.,
                     SancionIyC=A$Sanción.de.industria.y.Comercio,
                     SancionSobreG=A$Sanción.Sobretasa.Gasolina)

#Gráfico de la serie
with(Rentas, plot(Fecha, Predial, type = 'l'))

PredialT <- apply(matrix(Rentas$Predial, nrow =3), 2, sum)
fechasT  =  seq(as.Date("2005/1/1"),
               length.out = dim(A)[1]/3,  by  =  "quarter")

plot(fechasT, PredialT, type = 'l')

length(fechasT)

#Series para estimación y validación cruzada

y <- Rentas$Predial
m = 18
n = length(y)
yi = ts(y[1:(n-m)],frequency=12)
lyi = log(yi)
yf = ts(y[(n-m+1):n],frequency=12)
tf = seq(T+1,T+m,1)
tf2 = tf*tf
tf3 = tf2*tf
T = length(yi)
t = seq(1,T)
t2 = t*t
t3 = t2*t
library(forecast)

It = seasonaldummy(yi)


#Descomposición de la serie
M1 <- stl(yi, s.window = 'per', t.jump = 1)
plot(M1)

#Modelo de Descomposición

X1t = data.frame(rep(1, T), t, It)
DS1 <- data.frame(lyi,X1t)

mod1 = lm(lyi ~ t + It, data = DS1)
summary(mod1) 

M1 <- medidas(mod1, lyi, 12)
# El modelo exponencial-lineal-estacional
# genere una matriz con las variables 
# explicativas, mas columna de 1's

yhat1 = exp(fitted(mod1))

plot(t, yi, type='b')
lines(t, yhat1, col='red') #La serie tiene un efecto multiplicativo. La amplitud del efecto estacional es diferente en el trayecto
legend('topleft', col = c("black", "red"),
       lwd=2, bty='n', legend=c('Yi(PIB)', 'Yhat(PIB ajustado)'))
title(main='Serie real versus ajustada mediante el Modelo 1')

#Pronósticos del Modelo 1

t = seq(T+1,T+m,1)

# Nota: forecast version 7.1 cambiÃ³
# la funcion seasonaldummyf por
# seasonaldummy(y,m)

It = seasonaldummy(yi,m)

Xt = cbind(rep(1,m),t,It)

pron1 <- exp(predict.lm(mod1, newdata = data.frame(Xt = I(Xt))))

par(mfrow=c(1,1))
plot(t,yf, type = 'o', main = "Validación cruzada modelo 1")
lines(t,pron1, type = 'b', pch = 3,col='red' )

legend("topleft", 
       c("Obs","logcub+indic"), 
       pch = c(1, 3),
       col = c("black","red"))

pron1 = ts(pron1,frequency=12)
(A1 = accuracy(pron1,yf))

sum(tail(y, 18)-tail(y, 3))
sum(yf[1:15])

#Para el modelo 2
T = length(yi)
t = seq(1,T)
It = seasonaldummy(yi)

Xt = cbind(rep(1,T),t,It)
Ds = data.frame(yi,Xt) #arreglo con la variable dependiente y una matiz donde se guarda toda la informaciÃ³n de las covariables 
#(intercepto, t, t2, t3, IT)
theta.0 = coef(mod1) #vector donde se guardan los parÃ¡metros estimados

mod2 = nls(yi~exp(Xt%*%theta),
           data=Ds, start= list(theta=theta.0)) 
#Se utiliza el método de mínimos cuadrados no lineales para el modelo exponencial
summary(mod2) #Modelo cÃºbico exponencial ajusta y significacia de los estimadores (El modelo es competitivo porque hizo el trabajo de 
#capturar tendencia y estacionalidad. El objetivo es encontrar un modelo mejor que este)

yhat2 = fitted(mod2)

plot(t,yi,type='b')
lines(t, yhat2,col='blue') #La serie tiene un efecto multiplicativo. La amplitud del efecto estacional es diferente en el trayecto
legend('topleft', col = c("black", "blue"),
       lwd=2, bty='n', legend=c('Yi (PIB)', 'Yhat (PIB ajustado)'))
title(main='Serie real versus ajustada mediante el Modelo 2')

# al aumentar el tiempo, aumenta la amplitud de la componente estacional (en contraposiciÃ³n con un modelo con efecto aditivo)
#con esto se deduce que el modelo tiene otra ventaja y es su capacidad para capturar ese efecto multiplicativo de la tendencia y 
#estacionalidad identificado en la serie

k = 2+11
(M.exp.cub = medidas(mod2,yi,k))

#----------------- m pronosticos

tf = seq(T+1,T+m,1)

# Nota: forecast version 7.1 cambiÃ³
# la funcion seasonaldummyf por
# seasonaldummy(y,m)

Itf = seasonaldummy(yi,m)

Xtf = cbind(rep(1,m),tf,Itf)

pron2 = predict(mod2,
                data.frame(Xt = I(Xtf)))

par(mfrow=c(1,1))
plot(tf,yf, type = 'o')
lines(tf,pron2, type = 'b', pch = 3,col='red' )

legend("topleft", 
       c("Obs","exp cub+indic"), 
       pch = c(1, 3),
       col = c("black","red"))
title(main = "Pronóstico con modelo Exponencial cúbico estacional")

pron2 = ts(pron2,frequency=12)
(A2 = accuracy(pron2,yf))


sum(tail(y, 18)-tail(y, 3))
sum(pron2[1:15])


#Holt - Winters

mhw = HoltWinters(yi)
(c(mhw$alpha,mhw$beta,mhw$gamma))

yhat.hw = mhw$fitted[,"xhat"]
length(yhat.hw)
length(yi)

plot(t[-seq(1,12)],yi[-seq(1,12)],type='l')
lines(t[-seq(1,12)], yhat.hw,col='blue')

(M.hw = medidas.1(yhat.hw,yi,3))

pron.hw = predict(mhw, m, 
                  prediction.interval = FALSE)
#Pronósticos H-W

pron.hw = ts(pron.hw,frequency=12)
(A4=accuracy(pron.hw[,1],yf))

#Gráfico Validación Cruzada

plot(tf[1:15],yf[1:15], type = 'o')
lines(tf[1:15], pron.hw[1:15], type = 'b', pch = 5,col='orange')

sum(as.vector(pron.hw)[1:15])
sum(yf[1:15])
sum(pron1[1:15])
sum(pron2[1:15])
sum(pron.hw[1:15])


#LOESS

#SIMULACIÓN LOESS

It12 = cbind(It,(rep(1,nrow(It)) - apply(It,1,sum)))
Itf12 = cbind(seasonaldummy(yi,m),(rep(1,m) - apply(Itf12,1,sum)))

span <- seq(0.16, 1, 0.01)
Modelo <- list()
pron <- list()
Tt.loess <- list()
St.loess <- list()
ModEst <- list()
yhat4 <- list()
pron4 <- list()
pron34 <- list()
Utheil <- numeric()
Mape <- numeric()
Mpe <- numeric()

for (i in 1:length(span)) {
Modelo[[i]] <- loess(yi ~ t,span = span[i], modelo = TRUE,
                     control = loess.control(surface = "direct"))
pron[[i]] = predict(Modelo[[i]], data.frame(t = tf))
Tt.loess[[i]] = fitted(Modelo[[i]])
St.loess[[i]] = yi - Tt.loess[[i]]
ModEst[[i]] = lm(St.loess[[i]] ~ -1+ It12)
yhat4[[i]] = fitted(ModEst[[i]]) + Tt.loess[[i]]
pron4[[i]] = predict(ModEst[[i]],data.frame(It12=I(Itf12)))
pron34[[i]] = pron[[i]]+pron4[[i]]
pron34[[i]] = ts(pron34[[i]],frequency=12)
Utheil[i] <- accuracy(pron34[[i]],yf)[7]
Mape[i] <- accuracy(pron34[[i]],yf)[5]
Mpe[i] <- accuracy(pron34[[i]],yf)[4]
}

Results <- data.frame(span, Utheil, Mape, Mpe)
par(mfrow=c(2,1))
with(Results, plot(span, Utheil))
with(Results, plot(span, Mpe))

Minimo1 <- subset(Results, Utheil==min(Results$Utheil))
Minimo2 <- subset(Results, Mpe==min(Results$Mpe))

#----------Ejemplo 3: loess + indicadoras

#---------- ajuste Loess
mod3 = loess(yi ~ t,span = Minimo2[1] , modelo = TRUE,
             control = loess.control(surface = "direct"))

#---------------------pronosticos loess

pron3 = predict(mod3, data.frame(t = tf))

#---------- ajuste Loess + estacional


# modelo estacional sin la tendencia

Tt.loess = fitted(mod3)

St.loess = yi - Tt.loess

mod4 = lm(St.loess ~ -1+ It12)
summary(mod4)

yhat4 = fitted(mod4) + Tt.loess

plot(t,yi,type='b')
lines(t, yhat4,col='blue')

k = 1+12
(M.loess = medidas.1(yhat4,yi,k))


#---------------------pronosticos loess+estacional

pron4 = predict(mod4,data.frame(It12=I(Itf12)))

pron34 = pron3+pron4

par(mfrow=c(1,1))
plot(tf,yf, type = 'o')
lines(tf,pron34, type = 'b', pch = 3,col='red' )
legend("topleft", 
       c("Obs","loess+estacional"), 
       pch = c(1, 3),
       col = c("black","red"))

#-----------------errores de pronostico
pron34 = ts(pron34,frequency=12)
A=accuracy(pron34,yf)

predict.interval3 <- predict(mod4, data.frame(It12=I(Itf12)), interval = "prediction")
predict.interval4 <- predict(mod3, data.frame(t = tf), interval = "prediction")

#Pronósticos Fuera de la Muestra
predict.interval3+predict.interval4
colSums(predict.interval3+predict.interval4)

sum(tail(y, m))

#Pronóstico fuera de la muestra:
tt = 1:n
mod5 = loess(y ~ tt,span = Minimo2[1] , modelo = TRUE,
             control = loess.control(surface = "direct"))
tff = seq(n+1, (n+m), 1)
pron5 = predict(mod5, newdata = tff)

Tt.loess1 = fitted(mod5)

y <- ts(y, frequency = 12)
I_t = seasonaldummy(y)
It_12 = cbind(I_t,(rep(1,nrow(I_t)) - apply(I_t,1,sum)))

St.loess1 = y - Tt.loess1
mod6 = lm(St.loess1 ~ -1+ It_12)

It_f12 = seasonaldummy(y,m)
It_f12 = cbind(It_f12,rep(0,m))
It_f12[,12] = rep(1,m) - apply(It_f12,1,sum)

pron6 = predict(mod6, newdata = data.frame(It_12=I(It_f12)))

pron7 <- pron5 + pron6

sum(pron7)
plot(pron7, type = 'b')

ypron <- c(y, pron7)

sum(tail(ypron, 12))

sum(tail(ypron, 18))-sum(tail(ypron, 6))

plot(ypron, type = 'b')

summary(mod6)