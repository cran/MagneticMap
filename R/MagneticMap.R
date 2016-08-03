MagneticMap<-function(W, g, Plot2D=FALSE, Plot3D=FALSE)
{
  
  # Conditions on g
  if(g>=1/2  || g< 0)
  {print("wrong value of g (0 <= g < 1/2)")
  break;}
  
  # Symmetric weights
  Ws<-(W+t(W))/2
  
  # Edge flow
  A<-W-t(W)
  
  # Degree matrix
  D<-matrix(0, nrow=dim(W)[1], ncol=dim(W)[1])
  for(i in 1:dim(D)[1])
  {
  D[i,i]<-sum(Ws[i,])}
  
  # Transporter
  t<-matrix(0, nrow=dim(W)[1], ncol=dim(W)[1])
  for(i in 1:dim(t)[1])
  {
    for(j in 1:dim(t)[1])
    {
     t[i,j]<-exp(i*2*pi*g*A[j,i])  
    }
  }
  
  # Magnetic Laplacian
  Lg<-D-W*t
  
  # Normalized Laplacian
  ND<-D
  diag(ND)<-diag(ND)^(-1/2)
  Lng<-ND %*% Lg %*% ND
  
  # Eigenvectors
  Eigenvector<-eigen(Lng)$vectors
  
  # 2D Plot
  if(Plot2D==TRUE){plot(Eigenvector[,1],Eigenvector[,1],xlab="",ylab="")}
  
  # 3D Plot
  if(Plot3D==TRUE){scatterplot3d(Eigenvector[,1],Eigenvector[,2],Eigenvector[,3],xlab="",ylab="",zlab="")}
  

# Return Normalized Laplacian Matrix and Eigenvectors
object<-list(matrix=Lng, eigenvector=Eigenvector)
object
}