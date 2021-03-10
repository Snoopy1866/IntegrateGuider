#------------------------------------龙贝格求积 by 黄志恒-----------------------------------------------------
Romberg<-function(f,a,b,eps = 1e-5,it_max=50)
  #f为求积函数，a，b分别为积分上下限，eps为预先给定的精度，it_max为最大二分次数
{
  if(a==b)
  {
    stop('计算结果为0，请检查积分上下限')
  }    
  else if(a>b)
  {
    Romberg(f,b,a,eps,it_max=50)
  }
  else
  {
    Table<-matrix(NA,it_max,it_max)
    #建立一个50*50的T矩阵，存放计算结果
    h<-b-a
    Table[1,1]<-h/2*(f(a)+f(b))
    for(i in 2:it_max)
    {
      n<-2^(i-1)
      k<-(b-a)/2^(i-1)
      Table[i,1]<-Table[i-1,1]/2+k*sum(f(seq(1,n-1,2)*(b-a)/n))
      #从梯形公式出发，将区间逐次二分 
      for(j in 2:i)
      {
        m<-4^(j-1)
        Table[i,j]<-m/(m-1)*Table[i,j-1]-1/(m-1)*Table[i-1,j-1]
        #利用递推公式外推   
      } 
      R<-abs(Table[i,i]-Table[i-1,i-1])
      if(R<eps)  break
      #当达到预先给定的精度，中止计算
    }
    Table<-cbind(0:(i-1),h/2^(0:(i-1)),Table[1:i,1:i])
    #将二分次数，子区间长度，T矩阵合并成一个新矩阵，得到最终的T表
    dimnames(Table)<-list(paste(1:i),c('二分次数','h',paste('T',0:(i-1),sep='')))
    #给T表命名
    return(Table[i,i])
  }
}