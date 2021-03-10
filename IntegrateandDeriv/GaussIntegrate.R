#---------------------------------------------高斯求积 by 王文涛-----------------------------------------------------
GaussIntegrate <- function(f, a, b, n, gauss.method = NULL)
{
  #高斯求积公式
  #f -- 被积函数,为function类型
  #a,b -- 积分区间端点
  #n -- 使用n个高斯点求积
  #gauss.method -- 具体的高斯求积方法,为"Legendre","Chebyshev","Laguerre","Hermite"中的一种
  #                   未指定或指定错误将被忽略
  #---------------------------------------------
  if (n < 1)
    stop("至少要有1个高斯点!")
  if (a == b)
  {
    result <- 0
  }
  else if (a > b)
  {
    result <- -GaussIntegrate(f, b, a, n, gauss.method)
  }
  else if (!(is.infinite(a) || is.infinite(b)))
  {
    if (is.null(gauss.method))
      result <- GaussLegendre(f, a, b, n)
    else if (any(gauss.method == c("Laguerre","Hermite")))
      stop("有界区间gauss.method只能使用'Legendre','Chebyshev'")
    else if (gauss.method == "Chebyshev")
      result <- GaussChebyshev(f, a, b, n)
    else
      result <- GaussLegendre(f, a, b, n)
  }
  else if (a == -Inf && b != Inf)
  {
    g <- function(x)f(-x)
    result <- GaussLaguerre(g, 0, Inf, n) + GaussLegendre(g, -b, 0, n)
  }
  else if (a != -Inf && b == Inf)
  {
    result <- GaussLegendre(f, a, 0, n) + GaussLaguerre(f, 0, Inf, n)
  }
  else if (a == -Inf && b == Inf)
  {
    result <- GaussHermite(f, a, b, n)
  }
  return(result)
}

GaussLegendre <- function(f, a, b, n)
{
  #高斯-勒让德求积
  #f -- 被积函数,为function类型
  #a,b -- 积分区间端点
  #n -- 使用n点高斯-勒让德求积公式
  f_trans <- function(a, b, x)f((b - a)/2*x + (a + b)/2)*(b - a)/2
  #---------------------------------------------
  Polylegendre <- function(n)
  {
    #勒让德多项式Pn各项系数
    P1 <- c(0,1)
    P2 <- c(-1/2,0,3/2)
    if (n == 1)
      return(P1)
    else if (n == 2)
      return(P2)
    else
      return((2*n - 1)/n*c(0,Polylegendre(n - 1)) - (n - 1)/n*c(Polylegendre(n - 2),0,0))
  }
  #---------------------------------------------
  GaussLegendreCoef <- function(n)
  {
    #n点的高斯-勒让德求积公式的节点x和系数A
    x <- Re(polyroot(Polylegendre(n)))
    Y <- matrix(NA,2*n,n)
    for (i in 1:(2*n))
    {
      Y[i,] <- x^(i - 1)
    }
    b <- (1 - (-1)^(1:(2*n)))/(1:(2*n))
    A <- qr.solve(Y,b)
    list(x = x,A = A)
  }
  #---------------------------------------------
  coef <- GaussLegendreCoef(n)
  x <- coef$x
  A <- coef$A
  result <- sum(A*f_trans(a,b,x))
  return(result)
}

#----------------------------------------------------------------------
GaussChebyshev <- function(f,a,b,n)
{
  #高斯-切比雪夫求积
  #f -- 被积函数,为function类型
  #a,b -- 积分区间端点
  #n -- 使用n点高斯-切比雪夫求积公式
  #将区间[a,b]变换为区间[-1,1]
  f_trans <- function(a,b,x)f((b - a)/2*x + (a + b)/2)*(b - a)/2
  x <- cos((2*(1:n) - 1)/(2*n)*pi)
  result <- pi/n*sum(f_trans(a,b,x)*sqrt(1 - x^2))
  return(result)
}
#---------------------------------------------
GaussLaguerre <- function(f,a,b,n)
{
  #高斯-拉盖尔求积,计算[0,Inf]区间上的积分
  #f -- 被积函数,为function类型
  #a,b -- 积分区间端点
  #n -- 使用n点高斯-拉盖尔求积公式
  #---------------------------------------------
  PolyLaguerre <- function(n)
  {
    #拉盖尔多项式Ln各项系数
    L1 <- c(1,-1)
    L2 <- c(2,-4,1)
    if (n == 1)
      return(L1)
    else if (n == 2)
      return(L2)
    else
      return( - c(0,PolyLaguerre(n - 1)) + (2*n - 1)*c(PolyLaguerre(n - 1),0) - (n - 1)^2*c(PolyLaguerre(n - 2),0,0))
  }
  #---------------------------------------------
  GaussLaguerreCoef <- function(n)
  {
    #n点的高斯-拉盖尔求积公式的节点x和系数A
    vec <- PolyLaguerre(n)
    x <- Re(polyroot(vec))
    L <- NULL
    for (i in 1:n)
    {
      L[i] <- x[i]^(0:(n - 1))%*%(vec[2:(n + 1)]*(1:n))
    }
    
    A <- factorial(n)^2/(x*L[1:n]^2)
    list(x = x,A = A)
  }
  #---------------------------------------------
  coef <- GaussLaguerreCoef(n)
  x <- coef$x
  A <- coef$A
  result <- sum(A*exp(x)*f(x))
  return(result)
}

#----------------------------------------------------------------------
GaussHermite <- function(f,a,b,n)
{
  #高斯-埃尔米特求积,计算[-Inf,Inf]区间上的积分
  #f -- 被积函数,为function类型
  #a,b -- 积分区间端点
  #n -- 使用n点高斯-埃尔米特求积公式
  #---------------------------------------------
  PolyHermite <- function(n)
  {
    #埃尔米特多项式Hn各项系数
    H1 <- c(0,2)
    H2 <- c(-2,0,4)
    if (n == 1)
      return(H1)
    else if (n == 2)
      return(H2)
    else
      return(2*c(0,PolyHermite(n - 1)) - 2*(n - 1)*c(PolyHermite(n - 2),0,0))
  }
  #---------------------------------------------
  GaussHermiteCoef <- function(n)
  {
    #n点的高斯-埃尔米特求积公式的节点x和系数A
    vec <- PolyHermite(n)
    x <- Re(polyroot(vec))
    H <- NULL
    for (i in 1:n)
    {
      H[i] <- x[i]^(0:(n - 1))%*%(vec[2:(n + 1)]*(1:n))
    }
    A <- 2^(n + 1)*factorial(n)*sqrt(pi)/H[1:n]^2
    list(x = x,A = A)
  }
  #---------------------------------------------
  coef <- GaussHermiteCoef(n)
  x <- coef$x
  A <- coef$A
  result <- sum(A*exp(x^2)*f(x))
  return(result)
}