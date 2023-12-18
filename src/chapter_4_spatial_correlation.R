
# 공간 상관 계수 계산
# 필요한 라이브러리 불러오기
library(geosphere)

# 자료 읽기
bst <- read.csv("../data/KOEM_Busan_202208.csv")

# 위경도 정보 추출
LL <- bst[,6:7]

# Target parameters, 9, 10, 11 - WT, S, DO
xx <- bst$dissolved_oxygen 
xm <- mean(xx)
nst <- nrow(LL)
dmat <- matrix(0, nrow=nst, ncol=nst)

## 거리계산은 경도, 위도 순서로 입력 
for (ii in 1:nst) {
  for (jj in 1:nst) {
    p1 <- c(LL[ii,2], LL[ii,1])
    p2 <- c(LL[jj,2], LL[jj,1])
    dmat[ii,jj] <- distVincentyEllipsoid(p1, p2)
  }
}
## Basic computation
isup <- 0; isub <- 0; csup <- 0; csub <- 0; gsup <- 0; gsub <- 0
for (ii in 1:nst) {
  for (jj in 1:nst) {
    isup1 <- dmat[ii,jj]*(xx[ii]-xm)*(xx[jj]-xm)
    csup1 <- dmat[ii,jj]*(xx[ii]-xx[jj])^2
    gsup1 <- dmat[ii,jj]*xx[ii]*xx[jj]
    ifelse(ii != jj, gsub1 <- xx[ii]*xx[jj], gsub1 <- 0)
    isup <- isup + isup1; csup <- csup + csup1
    gsup <- gsup + gsup1; gsub <- gsub + gsub1
  }
  isub1 <- (xx[ii]-xm)^2; csub1 <- isub1 
  isub <- isub + isub1; csub <- csub + csub1
}

## Moran's I computation
wsum <- sum(dmat); nn <- nst
I_Moran <- nn*isup/(wsum*isub)

## Geary's c computation
c_Geary <- (nn-1)*csup/(2*wsum*csub)

## Getis-Ord G (mean = 1 normalization //
G_Getis_Ord <- nn*(nn-1)*gsup/(wsum*gsub)
c(I_Moran, c_Geary, G_Getis_Ord)

