#-------------------------------------复合辛普森和复合梯形求积 by 耿翔-----------------------------------------------------
#复合梯形公式
Trapezoid<-function(expr,f,a,b,eps=1e-5)
{ 
  if (a > b)
  {
    list <- Trapezoid(expr,f,b,a,eps)
    n <- list$区间个数
    result <- - list$积分结果
    return(list(区间个数 = n,积分结果 = result))
  }
  else
  {
    p<-DerivMaxValue(expr,a,b,2)#高阶导数最大值
    if (p == 0)
    {
        n <- 1
        result <- (f(a) + f(b))*(b - a)/2
    }
    else
    {
        n<-sqrt((b-a)^3*p/(12*eps))#计算n等份
        n<-ceiling(n)#n向上取整
        h<-(b-a)/n
        x<-a+(1:n)*h
        result<-(h/2*(f(a)+2*sum(f(x)[1:n-1])+f(b)))#复合梯形公式计算积分
    }
    return(list(区间个数 = n,积分结果 = result))
  }
}
#---------------------------------------------
#复合辛普森公式
Simpson<-function(expr,f,a,b,eps=1e-5)
{
  if (a > b)
  {
    list <- Simpson(expr,f,b,a,eps)
    n <- list$区间个数
    result <- - list$积分结果
    return(list(区间个数 = n,积分结果 = result))
  }
  else
  {
    q<-DerivMaxValue(expr,a,b,4)#高阶导数最大值
    if (q == 0)
    {
        n <- 1
        result <- (f(a) + 4*f((a + b)/2) + f(b))*(b - a)/6
    }
    else
    {
        n<-((b-a)^5*q/(2880*eps))^(1/4)#计算n等份
        n<-ceiling(n)#n向上取整
        h<-(b-a)/n
        x<-a+(0:n)*h
        x.mid<-x[1:n]+h/2
        result<-(h/6*(f(a)+f(b)+4*sum(f(x.mid[1:n]))+2*sum(f(x[2:n]))))#复合辛普森计算积分
    }
    return(list(区间个数 = n,积分结果 = result))
  }
}
#---------------------------------------------
DerivMaxValue <- function(expr,c,d,m)
{
  #求函数F的m阶导的最大值
  #F -- function类型
  #c,d -- 区间端点
  #m -- 求导阶数
  expr <- DD(expr,"x",m)
  Value <- function(x) abs(eval(expr))
  optimize(Value,c(c,d),maximum = TRUE)$objective
}
#---------------------------------------------
DD <- function(expr, name, order = 1) 
{
  if (order < 1)
    stop("'order' must be >= 1")
  if (order == 1)
    D(expr, name) 
  else 
    DD(D(expr, name), name, order - 1)
}