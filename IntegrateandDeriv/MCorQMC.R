#------------------------------------------(拟)蒙特卡洛求积 by 袁愈新-----------------------------------------------------
MCorQMC <- function(f,n=1e5,set_seed=1000,a,b,fail1=NULL,
                    fail2=NULL,base=NULL,method=c("MC","QMC")){
  #f为被积函数,n为随机数向量个数,set_seed为蒙特卡洛方法预先设置的种子数
  #a,b分别为自变量所能达到的的下限和上限的极限构成的向量；
  #fail1,fail2分别为自变量的积分上限、下限约束条件组成的向量，默认情况为矩形区域
  #即此时其值分别a,b，故将这两向量的值均赋为空值
  #base为QMC所需基数构成的向量，其中基数必须为质数。当采用MC方法其值为空
  
  if (is.null(a) | is.null(b)) stop("下限和上限必须存在")
  if(length(a)!=length(b)) stop("下限和上限不对应")
  if (all(method==c("MC","QMC"))) method <- "QMC"
  #如果未指定方法，默认采用QMC方法
  if(is.null(base)) method <- "MC"
  #如果未指定数基，则改用MC方法
  if (method!="MC" & method!="QMC" ) stop("'method' should be one of 'MC', 'QMC'")
  
  p <- length(a) #p为自变量个数
  A <- matrix(a,n,p,byrow=TRUE)
  B <- matrix(b-a,n,p,byrow=TRUE)
  #A,B为将标准多维空间阵化为一般多维空间所需的矩阵
  
  if (method=="MC"){
    #采用蒙特卡洛方法求积分
    set.seed(set_seed) #设置种子数
    U <- matrix(runif(p*n),n) 
    #构建一个n*p的呈0-1均匀分布矩阵,U[i,j]为第i个随机数向量的第j个元素
  }
  
  else{
    #采用拟蒙特卡洛方法求积分
    if (length(base)!=length(a))  stop("数基中的元素个数过度")
    Judgingprime <- function(x){
      if(x>2){
        if (any(x %% 2:(x - 1) == 0))
          stop("数基必须为质数")
      }}
    sapply(base,Judgingprime) #确保数基都是质数
    
    U <- matrix(0,n,p)
    NumBits <- 1+ceiling(log(n)/log(base)) 
    #NumBits为这些数基构成从0到n的任意整数的表达式的最高阶数加1组成的向量
    g <- function(NumBits) -c(1:NumBits)
    bj <- lapply(NumBits,g) 
    #生成数基向量中的元素对应构成任意整数的所以可能的阶数
    #并生成一个列表来来接收不同数基对应的阶数
    for(r in 1:p){
      u <- rep(0,times=length(bj[[r]]))
      for(i in 1:n){
        j=1
        ok=0
        while (ok==0){
          u[j]=u[j]+1
          if (u[j]<base[r])  ok=1
          else{
            u[j]=0
            j=j+1}}
        
        VetBase=base[r]^(bj[[r]]) 
        U[i,r]=u%*%VetBase 
        #U[i,r]为第i个符合标准超均匀分布的向量对应的第r个数基构成的值
        #其值为整数对应的根式逆函数对应的值
      }
    }
  }
  
  X <- A+B*U  
  #将其化一般多维序列,其中矩阵X的元素的值为X[,i]=A[i]+(B[i]-A[i])*U[,i]
  choosex <- function(x){
    if(is.null(fail1)) f(x)
    else{
      l <- (fail1(x)<=x & x<=fail2(x))
      if(all(l)) f(x)}
  }
  t <- unlist(apply(X,1,choosex)) #当满足约束条件时，执行被积函
  s <- sum(t) #对满足条件的函数值求和
  result <- (prod(b-a)*s)/n#得出积分估计值
  return(result)
}


#example
if(FALSE){
  #1.数值积分127页例14
  f <- function(x) log(x[1]+2*x[2])
  a <- c(1.4,1);b <- c(2,1.5)
  
  MCorQMC(f,a=a,b=b,method = "MC")
  MCorQMC(f,a=a,b=b,base=c(2,3)) #默认method为QMC
  
  #结果：用MC、QMC求得的近似解分别为0.4295554、0.4295518，真实解：0.4295545265
  
  #2.高数下册167页例11(2)
  f <- function(x) x[3]^3*sin(x[2])
  a <- c(0,0,0);b <- c(2*pi,pi/2,1)
  fail1 <- function(x) c(0,0,0)
  fail2 <- function(x) c(2*pi,pi/2,cos(x[2]))
  
  MCorQMC(f,a=a,b=b,fail1=fail1,fail2=fail2,method="MC")
  MCorQMC(f,a=a,b=b,fail1=fail1,fail2=fail2,base=c(2,3,5),method="QMC")
  
  #结果:用MC、QMC求得的近似解分别为0.3186149、0.3144375143，真实解：pi/10=0.31415926
  
}


#比较MC与QMC抽取的二维单位正方形的散点图
if (FALSE){
  
  opar <- par(no.readonly = T)
  par(mfrow=c(1,2))
  par(pin=c(3,3)) #定义图形宽和高
  n <- 2000;p <- 2
  #抽取点的个数为2000,维数为2
  MCU <- matrix(runif(n*p),n) 
  QMCU <- matrix(0,n,p)
  #MCU、QMCU分别为用MC、QMC方法得到的散点坐标构成的矩阵
  base <- c(2,3) #定义数基为2和3
  NumBits <- 1+ceiling(log(n)/log(base)) 
  #NumBits为这些数基构成从0到n的任意整数的表达式的最高阶数加1组成的向量
  g <- function(NumBits) -c(1:NumBits)
  bj <- lapply(NumBits,g) 
  #生成数基向量中的元素对应构成任意整数的所以可能的阶数
  #并生成一个列表来来接收不同数基对应的阶数
  for(r in 1:2){
    u <- rep(0,times=length(bj[[r]]))
    for(i in 1:n){
      j=1
      ok=0
      while (ok==0){
        u[j]=u[j]+1
        if (u[j]<base[r])  ok=1
        else{
          u[j]=0
          j=j+1}}
      
      VetBase=base[r]^(bj[[r]]) 
      QMCU[i,r]=u%*%VetBase 
      #QMCU[i,r]为第i个符合标准超均匀分布的向量对应的第r个数基构成的值
      #其值为整数对应的根式逆函数对应的值
    }
  }
  plot(MCU[,1],MCU[,2],main="用蒙特卡罗方法抽取的二维单位正方形上的点",
       xlab="x1",ylab="x2")
  plot(QMCU[,1],QMCU[,2],main="用拟蒙特卡罗方法抽取的二维单位正方形上的点",
       xlab="x1",ylab="x2")
  par(opar)
}