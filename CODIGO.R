#LIBRERIAS NECESARIAS
library(kableExtra)
library(sp)
library(ggplot2)
library(ggrepel)
library("factoextra")
library(ggcorrplot)
library(svDialogs)
library(tidyverse)
library(readxl)
library(dplyr)
library(sf)
library(rgdal)
library(ggspatial)
library(dplyr)
library(sf)

#Cargar el shapefile
setwd("D:/USUARIO/Documents/OneDrive/Trabajo de grado/R")
getwd()
deptos<-read_sf("mpio.shp")

#Filtrar por Santander (dpto=68)
stder=deptos %>% filter(DPTO=="68")

#Cargar base de datos a usar
bf<-read_excel("D:/USUARIO/Documents/OneDrive/Trabajo de grado/BD.xlsx")

#Corregir un dato
stder[75,8]<- "EL PECON"

#Cargar base de datos con indices y geometria de los municipios
mpios<-st_read("D:/USUARIO/Documents/OneDrive/Trabajo de grado/R/icv.shp")

 
#INDICE POR ANALISIS DE COMPONENTES PRINCIPALES
i=0
x<-c(0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1) #vector que asigna la polaridad de la variable

#Funcion del ICV_PCA
ICVPCA<- function(b,x,shp){
    #Funcion que normaliza las variables dependiendo de su polaridad asignada
    normalizar<-function(b){
    
    bd<-cbind(b[,1],matrix(0,nrow=nrow(b), ncol=ncol(b)-1))
    colnames(bd)<-colnames(b)

    for (i in 1:length(x)){
    n<-x[i]
    if (n==0){
      bd[,i+1]<-abs((min(b[,i+1])-b[,i+1]))/(max(b[,i+1])-min(b[,i+1]))}
    else
    {bd[,i+1]<-abs((max(b[,i+1])-b[,i+1]))/(max(b[,i+1])-min(b[,i+1]))
    }
  }
  return(bd)
   }

  bd<-normalizar(b)
  
  #Correlaciones para el data frame normalizado
  corr<-cor(bd[,-1])
  
  #Analisis de Componentes Principales
  acp<-prcomp(bd[,-1], scale=TRUE)
  print(summary(acp))#Resumen acp
  print(biplot (acp, cex = c (0.6, 0.8), repel=TRUE)) #biplot PC1 vs PC2

  #Valores Propios
  eig.val <- get_eigenvalue(acp)
  egv<-eig.val$eigenvalue
  ev<-c()
  
  #Seleccionar las componentes principales
  j=1
 while (egv[j]>1){
      ev[j]<-egv[j]
      j<-j+1
    }
  class(ev)<-"numeric"
  
#CONSTRUCCION DEL INDICE POR PCA
  
  #Denominador
  vp<-sqrt(ev)#raiz valores propios
  vp<-as.matrix(vp)
  vlrp<-matrix(rep(t(vp),length(vp)), nrow = length(vp), ncol = length(vp))#VALORES PROPIOS 32x32
  
  den<-colSums(vp)#denominador
  
  #Numerador
  matrizd<-as.matrix(bd[,-1])
  
  pc<-acp$x[,1:length(vp)]#PUNTUACIONES
  numerador<-pc%*%vlrp#NUMERADOR 32 REPETICIONES
  
  numer<-numerador[,1]#numerador

  #Aplicar f??rmula
  indice<-numerador[,1]/den
  
  #Unir al el indice a los municipios
  indmpio<-cbind.data.frame(NOMBRE_MPI=bd[,1], indice)
  
  #Normalizar el indice
  indmpio$icvpca=(indmpio$indice - min(indmpio$indice))/(max(indmpio$indice)-min(indmpio$indice))
  indmpio<-indmpio[,!colnames(indmpio) %in% "indice"]
  
  #Unir al data frame
  bd<-cbind.data.frame(bd,indmpio[,2])
  
  
  #GRAFICAR RESULTADOS
  
  #Leer el shapefile de Santander
  deptos<-read_sf("mpio.shp")
  stder=deptos %>% filter(DPTO=="68") 
  
  #Unir por nombre de municipios
  municipios = merge(stder,indmpio,by="NOMBRE_MPI")

  #Graficar los resultados del icvpca

  grafico <- ggplot(municipios) + 
    geom_sf(aes(fill = icvpca)) + 
    annotation_scale(style = "ticks") +
    annotation_north_arrow(location = 'tl', height = unit(1, "cm"), width = unit(1, "cm")) + 
    xlab("Longitud") + 
    ylab("Latitud") + 
    ggtitle("Indice por PCA") +
    scale_fill_gradient(name = "Indice", low = "red", high = "green")
  
  print(grafico)
  
  #Tablas top 7
  ind<- indmpio[order(indmpio$icvpca, decreasing = TRUE), ]  # Ordenar los datos de manera descendente
  top_7m <- head(ind, 7)  # Seleccionar los 7 primeros valores mas altos
  
  # Imprimir la tabla de los 7 valores mas altos 
  print(top_7m %>%
    kable() %>%
    kable_styling(full_width = FALSE) %>%
    add_header_above(c("Top 7 Valores Mas Altos" = 3)))
  
  ind1<- indmpio[order(indmpio$icvpca, decreasing = FALSE), ]  # Ordenar los datos de manera descendente segun la columna deseada
  top_7p <- head(ind1, 7)  # Seleccionar los 7 primeros valores mas altos
  
  # Imprimir la tabla de los 7 valores mas bajos con kableExtra
  print(top_7p %>%
          kable() %>%
          kable_styling(full_width = FALSE) %>%
          add_header_above(c("Top 7 Valores Mas Bajos" = 3)))
  
  #La funcion devuleve los valores del icvpca
  return(indmpio)
}

#Aplicaci??n de la funcion ICVPCA
ICVPCA(bf,x,shp)


#INDICE POR DISTANCIA P2
#Funcion del ICV_P2
ICVP2<-function(b,epsilon,iteraciones,shp){
  NOMBRE_MPI<-b[,1]
  bd<-as.matrix(b[,-1])
  class(bd)<-"numeric"
  #Contar filas y columnas
  n<-ncol(bd)
  m<-nrow(bd)
  
  #Funci??n para calcular la distancia de Frechet
  frechet<-function(bd){
    class(bd)<-"numeric"
    bd1<-matrix(0,nrow=87,ncol=32)
    colnames(bd1)<-colnames(bd)
    for(i in 1:87) {
      for(j in 1:32) {
        bd1[i,j]=abs((bd[i,j]-bd[9,j]))/(sd(bd[,j])*(sqrt((nrow(bd)-1)/nrow(bd))))
      }
    }
    
    bdf<-as.matrix(cbind(bd1,rowSums(bd1)))
    return(bdf)}
  
  f<-frechet(bd)
  
  #Funcion para ordenar de forma decreciente el data frame con frechet
  ordenar<-function(b){
    c=cor(f)
    orden<-sort(abs(c[1:nrow(c)-1,ncol(c)]),decreasing = TRUE)
    bdo<-subset(f,select=names(orden))
    bdo1<-as.matrix(bdo)
    class(bdo1)<-"numeric"
    return(bdo1)
  }
  
  bdo1<-ordenar(f)
  
  #Funcion para realizar las regresiones y obtener el R^2
  regresiones<-function(b){
    
    class(bdo1)<-"numeric"
    C<-rep(1,32)
    
    for (i in 2:32) {
      modelo<- lm(bdo1[,i]~bdo1[,1:i-1])
      C[i]<-(1-summary(modelo)$r.squared)
    }
    return(C)
  } 
  
  r<-regresiones(bdo1)
  
  #Funcion para calcular el Z
  zeta<-function(b){
    bdd<-matrix(0,nrow=87, ncol=32)
    for(i in 1:87) {
      for(j in 1:32) {
        bdd[i,j]<-bdo1[i,j]*r[j]
      }
    }
    z<-as.array(rowSums(bdd))
    return(z)
  }
  
  z<-zeta(bdo1)
  
  
  dist<-Inf #asignar una distancia infinita
  iteracion <- 0
  Z<-matrix(f[,33])#asigna a Z los valores de la distancia de Frechet
  
  #Proceso de la DP2
  while ((epsilon <= dist) && (iteraciones<=iteracion)){
    
    f[,33]<-z#asigna el valor de z a los valores iniciales de frechet
    bdo1<-ordenar(f) #ordena
    r<-regresiones(bdo1)#realiza las regresiones
    z<-zeta(bdo1)
    Z<-cbind(Z,z)
    distancia<-(f[,ncol(f)]-z)^2
  
    dist<-sum(distancia)#calculo de la condicion para decidir
    
    iteracion=iteracion+1
     if(iteracion>=iteraciones){
       warning("The itertions stops because the error is stable")
       break
    }
    
  }
  
  icvp2<-Z[,ncol(Z)]#renombrar y separar el icvp2
  
  #Unir y normalizar el icvp2
  indp2<-cbind.data.frame(NOMBRE_MPI,icvp2)
  indp2$icvp2=1-((indp2$icvp2-min(indp2$icvp2))/(max(indp2$icvp2)-min(indp2$icvp2)))
  
  
  #Unir a municipios el icvp2
  municipios = merge(stder,indp2,by="NOMBRE_MPI")
  
  #Graficar el icvp2
  grafico = ggplot(municipios) + geom_sf(aes(fill=icvp2))+ annotation_scale(style="ticks")+
    annotation_north_arrow(location='tl',height = unit(1, "cm"), width = unit(1, "cm"),) +
    xlab("Longitud") + ylab("Latitud") + ggtitle("Indice por DP2")+
    scale_fill_gradient(name = "DP2", low = "red", high = "green")
  
  print(grafico)
  
  #Tablas del top
  ind<- indp2[order(indp2$icvp2, decreasing = TRUE), ]  
  top_7m <- head(ind, 7)
  
  # Imprimir la tabla de los 7 valores mas altos 
  print(top_7m %>%
          kable() %>%
          kable_styling(full_width = FALSE) %>%
          add_header_above(c("Top 7 Valores Mas Altos" = 3)))
  
  ind1<- indp2[order(indp2$icvp2, decreasing = FALSE), ]  
  top_7p <- head(ind1, 7)  
  
  # Imprimir la tabla de los 7 valores mas bajos
  print(top_7p %>%
          kable() %>%
          kable_styling(full_width = FALSE) %>%
          add_header_above(c("Top 7 Valores Mas Bajos" = 3)))
  return(indp2)
}

ICVP2(bf,0.00001,17,stder)




