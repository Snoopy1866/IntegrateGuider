#------------------------------------------牛顿-科特斯求积 by 耿翔-----------------------------------------------------
Cotes<-function(n) #定义柯特斯系数函数
{  
  M<-rep(0,n+1)        #定义一个任意向量
  for(k in 0:n)
  {
    F<-function(t)   #定义积分函数
    {
      x<-1
      for(j in 0:n)
      {
        if(j!=k)
          x<-x*(t-j)
      }
      x
    }
    a<-((-1)^(n-k)/(n*factorial(k)*factorial(n-k)))
    b<-integrate(F,0,n) #求积分
    M[k+1]<-a*b$value
  }
  M
}
#---------------------------------------------

#牛顿柯特斯公式
NewtonCotes<-function(f,a,b,n)
{
    h<-(b-a)/n
    x<-a+(0:n)*h
    NC<-(b-a)*sum((Cotes(n))*f(x))
    NC
}