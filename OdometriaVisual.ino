#include <TimerOne.h>

double T[3];     double R[3][3];  int ban=0;
String msj; String aux; double mp1[8][2]; double mp2[8][2];
double pi=3.14159265;

//This function sends the Traslation and Rotation from Serial every 100 miliseconds
void callback()
{
      Serial.println("T : ");
      Serial.print(T[0]); Serial.print(","); Serial.print(T[1]); Serial.print(","); Serial.println(T[2]); 
      Serial.println("R : ");
      DespliegaMatriz(R);
}

void setup() 
{
    Timer1.initialize(100000);  //100 miliseconds.    
    Timer1.pwm(9, 512);
    Timer1.attachInterrupt(callback);
    Serial.begin(9600);
}

void loop() 
{
    double A[8][8];      double x[8];      double E[3][3];      double U[3][3];      double S[3][3];      double V[3][3];   
    double K[3][3] = {{1,0,0},{0,1,0},{0,0,1}};    
    //First, we recibe from Serial the matching points.
    //If is the first time,  we save them in mp1.
    //If not, we save them in mp2 and copy mp2 to mp1.
    if (Serial.available())
    {
        int cont=0; int i=0; int j=0;
        msj = Serial.readStringUntil('\n');
        while (msj.length() > 0)
        {
          int index = msj.indexOf(',');
          if (index == -1) // No space found
          {
            break;
          }
          else
          {       
            i=cont/2; j=cont%2;              
            if(ban==0)
            {
                aux=msj.substring(0, index);
                mp1[i][j] = aux.toDouble();
                mp2[i][j] = 0.0;
            }     
            else
            {
                mp1[i][j] = mp2[i][j];
                aux=msj.substring(0, index);
                mp2[i][j] = aux.toDouble();
            }
            msj = msj.substring(index+1);
          }
        }
        ban=1;
    } 
    /* Test Example.   
      double mp1[8][2]={{0.8763,0.2621},{1.3939,0.7830},{1.7518,1.3459},{1.5114,1.0647}, {1.2113,0.3775},{1.3274,0.6700},{0.1429,0.1016},{5.9979,7.1899}};
      double mp2[8][2]={{0.5324,0.0374},{0.4376,0.0356},{0.5929,0.2377},{0.6549,0.2768},{0.6405,0.0718},{0.6062,0.1430},{0.1860,-0.0786},{0.8200,0.6884}};      
    */  
      double nmp1[8][2];  double nmp2[8][2];
      //Normalization of the matching points.
      Normaliza(mp1,nmp1);
      Normaliza(mp2,nmp2);
      //Construction of the matrix A.
      ConstruyeA(nmp1,nmp2,A);
      //Solve for the Fundamental Matrix elements x.
      Resuelve(A,x);      
      //We calculated the Essential Matrix E using the normalization matrix.
      DesNormaliza(mp1,mp2,x,E);
      //Display the Esential Matrix.
      //DespliegaMatriz(E);      
      //We find the USV Decomposition of the E Matrix.
      Factoriza(E,U,S,V);      
      //Finally we estimate the Traslation vector T and the Rotation Matrix R
      CalculaTraslacion(U,T);
      CalculaRotacion(U,V,R);      
}

//Calculate Rotation Matrix as Hartley & Zisserman’s book proposed.
void CalculaRotacion(double U[3][3], double V[3][3], double R[3][3])
{
      double newE[3][3]; double newU[3][3]; double newS[3][3]; double newV[3][3]; double aux[3][3]; double diag[3][3]; double VT[3][3];
      for(int i=0;i<3; i++)
      {
        for(int j=0;j<3; j++)
        {
            diag[i][j]=0.0;            
            VT[i][j]=V[j][i];
        }          
      }
      diag[0][0]=1.0; diag[1][1]=1.0;
      Mul(U,diag,aux);
      Mul(aux,VT,newE);
      Factoriza(newE,newU,newS,newV);            
      diag[0][0]=0.0;       diag[1][1]=0.0;   diag[0][1]=-1.0; diag[1][0]=1.0; diag[2][2]=1.0;
      Mul(newU,diag,aux);
      for(int i=0;i<3; i++)
      {
        for(int j=0;j<3; j++)
        {
            diag[i][j]=0.0;            
            VT[i][j]=newV[j][i];
        }          
      }
      Mul(aux,VT,R);
      if(R[0][0]*R[1][1]*R[2][2]<0.0)
      {
          for(int i=0;i<3; i++)
          {
            for(int j=0;j<3; j++)
            {
                R[i][j]=-1.0*R[i][j];
            }          
          }
      }
}

//Display a 3x3 Matrix.
void DespliegaMatriz(double A[3][3])
{
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
            Serial.print(A[i][j]);        Serial.print(",");
        }                
        Serial.println(".");
      }      
}

//Calculate Translation Vector as Hartley & Zisserman’s book proposed.
void CalculaTraslacion(double U[3][3], double T[3])
{
      double Z[3][3];      double Tt[3][3];      double aux[3][3];
      for(int i=0;i<3; i++)
      {
        for(int j=0;j<3; j++)
        {
            Z[i][j]=0.0;
        }          
      }
      Z[0][1]=1.0;
      Z[1][0]=-1.0;
      Mul(U,Z,aux);
      for (int i=0;i<3; i++)
      {
              for (int j=0; j<3; j++)
              {
                  Z[i][j]=U[j][i];
              }
      }
      Mul(aux,Z,Tt);
      T[0]=Tt[2][1]; 
      T[1]=Tt[2][0];      
      T[2]=Tt[0][1];
}

//Estimate SVD Decomposition.
void Factoriza(double E[3][3], double U[3][3], double S[3][3], double V[3][3])
{    
    double ET[3][3];
    //Calculate Transpose of E
    for (int i=0;i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            ET[i][j]=E[j][i];
        }
    }
    //Calculate A=E^T E
    double A[3][3];
    Mul(ET,E,A);  
    //Calculate the eigenvalues of A.
    double p1=A[0][1]*A[0][1]+A[0][2]*A[0][2]+A[1][2]*A[1][2];
    double q=(A[0][0]+A[1][1]+A[2][2])/3.0;
    double p2 = (A[0][0] - q)*(A[0][0] - q) + (A[1][1] - q)*(A[1][1] - q) + (A[2][2] - q)*(A[2][2] - q) + 2 * p1;
    double p = sqrt(p2 / 6.0);
    //Solve B = (1 / p) * (A - q * I)
    double B[3][3];  
    double aux;
    for (int i=0;i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            if(i==j)
            {
                aux=q;
            }
            else
            {
                aux=0;
            }
            B[i][j]=(1/p)*(A[i][j]-aux);
        }
    }
    double r=B[0][0]*B[1][1]*B[2][2]+B[0][1]*B[1][2]*B[2][0]+B[0][2]*B[1][0]*B[2][1]-(B[2][0]*B[1][1]*B[0][2]+B[2][1]*B[1][2]*B[0][0]+B[2][2]*B[1][0]*B[0][1]);
    //probably unnecessary but helps to correct errors caused by numerical asymmetry
    if (r <= -1)
    {
        aux = pi / 3.0;
    } 
    else
    {
        if(r >=1)
        {
            aux = 0.0;
        }      
        else
        {
            aux = acos(r) / 3;
        }        
    }  
    double eig[3];
    eig[0] = q + 2 * p * cos(aux);
    eig[2] = q + 2 * p * cos(aux + (2*pi/3));
    eig[1] = 3 * q - eig[0] - eig[2];
    //Calculate the eigenvectors of A 
    double eigv1[3];
    double eigv2[3];
    double eigv3[3];
    //Two columns of  A − λ I
    double aux1[3]={A[0][0]-eig[0],A[1][0],A[2][0]};
    double aux2[3]={A[0][1],A[1][1]-eig[0],A[2][1]};        
//    double aux2[3]={A[0][0]-eig[0],A[1][0],-1*A[2][1]};        
    //We use cross product to calculate the vectors of V
    Cruz(aux1,aux2,eigv1);
    Cruz(eigv1,aux1,eigv2);
    Cruz(eigv1,eigv2,eigv3);
    //Make the eigenvectors unitary.
    aux=sqrt(eigv1[0]*eigv1[0]+eigv1[1]*eigv1[1]+eigv1[2]*eigv1[2]);
    for (int i=0;i<3; i++)
    {
        eigv1[i]=eigv1[i]/aux;
    }
    aux=sqrt(eigv2[0]*eigv2[0]+eigv2[1]*eigv2[1]+eigv2[2]*eigv2[2]);
    for (int i=0;i<3; i++)
    {
        eigv2[i]=eigv2[i]/aux;
    }
    aux=sqrt(eigv3[0]*eigv3[0]+eigv3[1]*eigv3[1]+eigv3[2]*eigv3[2]);
    for (int i=0;i<3; i++)
    {
        eigv3[i]=eigv3[i]/aux;
    }
    for (int i=0;i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            
                S[i][j]=0;            
        }
    }
    S[0][0] = sqrt(eig[0]);
    S[1][1] = sqrt(eig[0]);

    V[0][0] = eigv2[0]; V[1][0]=eigv2[1]; V[2][0]=eigv2[2];
    V[0][1] = eigv1[0]; V[1][1]=eigv1[1]; V[2][1]=eigv1[2];
    V[0][2] = eigv3[0]; V[1][2]=eigv3[1]; V[2][2]=eigv3[2];

    //Lastly, calculate U.
    U[0][0] = (1/eig[0])*(E[0][0]*V[0][0]+E[0][1]*V[1][0]+E[0][2]*V[2][0]);
    U[1][0] = (1/eig[0])*(E[1][0]*V[0][0]+E[1][1]*V[1][0]+E[1][2]*V[2][0]);
    U[2][0] = (1/eig[0])*(E[2][0]*V[0][0]+E[2][1]*V[1][0]+E[2][2]*V[2][0]);

    U[0][1] = (1/eig[0])*(E[0][0]*V[0][1]+E[0][1]*V[1][1]+E[0][2]*V[2][1]);
    U[1][1] = (1/eig[0])*(E[1][0]*V[0][1]+E[1][1]*V[1][1]+E[1][2]*V[2][1]);
    U[2][1] = (1/eig[0])*(E[2][0]*V[0][1]+E[2][1]*V[1][1]+E[2][2]*V[2][1]);

    aux1[0]=U[0][0]; aux1[1]=U[1][0]; aux1[2]=U[2][0];
    aux2[0]=U[0][1]; aux2[1]=U[1][1]; aux2[2]=U[2][1];  
    Cruz(aux1,aux2,eigv1);

    U[0][2]=eigv1[0]; U[1][2]=eigv1[1]; U[2][2]=eigv1[2];    

    //And finally make U vectors unitary.
    aux=sqrt(U[0][0]*U[0][0]+U[1][0]*U[1][0]+U[2][0]*U[2][0]);
    for (int i=0;i<3; i++)
    {
        U[i][0]=U[i][0]/aux;
    }
    aux=sqrt(U[0][1]*U[0][1]+U[1][1]*U[1][1]+U[2][1]*U[2][1]);
    for (int i=0;i<3; i++)
    {
        U[i][1]=U[i][1]/aux;
    }
    aux=sqrt(eigv1[0]*eigv1[0]+eigv1[1]*eigv1[1]+eigv1[2]*eigv1[2]);
    for (int i=0;i<3; i++)
    {
        U[i][2]=U[i][2]/aux;
    }

}
//Cross product.
void Cruz(double A[3], double B[3], double C[3]) 
{
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
}

//Solve the system Ax=[-1 ... -1]^T
void Resuelve(double A[8][8], double x[8])
{
    double E[8][9];
    double pivote;
    //Extended Matrix
    for(int i=0; i<8; i++)
    {
        for (int j=0; j<8; j++)
        {
            E[i][j]=A[i][j];
        } 
        E[i][8]=-1;
    }
    //Gauss principal cycle            
    for(int i=0; i<8; i++)
    {         
        pivote=E[i][i];
        //Divide the row i by the pivot.
        for(int j=0; j<9; j++)
        {
            E[i][j]=E[i][j]/pivote;
        }        
        for(int j=0; j<8; j++)
        {
            if(j!=i)
            {
              double aux=E[j][i];
              for(int k=0; k<9; k++)
              {
                  //Row Operation to make the zeroes.
                  E[j][k]=E[j][k]-E[i][k]*aux;
              }
            }                            
        }
    }
    //Extract the values of the last column of the extended matrix E.
    for(int i=0; i<8; i++)
    {
      x[i]=E[i][8];
    }
}
//Build A matrix.
void ConstruyeA(double mp1[8][2], double mp2[8][2], double A[8][8])
{
    for(int i=0; i<8; i++)
    {
          A[i][0]=mp1[i][0]*mp2[i][0];
          A[i][1]=mp2[i][0]*mp1[i][1];
          A[i][2]=mp2[i][0];
          A[i][3]=mp1[i][0]*mp2[i][1];
          A[i][4]=mp2[i][1]*mp1[i][1];
          A[i][5]=mp2[i][1];
          A[i][6]=mp1[i][0];
          A[i][7]=mp1[i][1];
    }  
}
//Normalization of matched points.
void Normaliza(double mp[8][2] , double Nmp[8][2])
{
    double Cx=0.0;
    double Cy=0.0;
    double dis=0.0;
    for(int i=0; i<8; i++)
    {
        Cx=Cx+mp[i][0];
        Cy=Cy+mp[i][1];
     }
    Cx=Cx/8.0;
    Cy=Cy/8.0;
    for(int i=0; i<8; i++)
    {
        dis=dis+sqrt((mp[i][0]-Cx)*(mp[i][0]-Cx)+(mp[i][1]-Cy)*(mp[i][1]-Cy));
    }
    dis=dis/8.0;
    for(int i=0; i<8; i++)
    {
        Nmp[i][0]=(sqrt(2)/dis)*mp[i][0]-(sqrt(2)/dis)*Cx;
        Nmp[i][1]=(sqrt(2)/dis)*mp[i][1]-(sqrt(2)/dis)*Cy;        
    }
}
//Undo the normalization in the E Matrix.
void DesNormaliza(double mp1[8][2], double mp2[8][2], double x[8], double E[3][3])
{
    double Cx1=0.0;
    double Cy1=0.0;
    double dis1=0.0;
    double Cx2=0.0;
    double Cy2=0.0;
    double dis2=0.0;
    for(int i=0; i<8; i++)
    {
        Cx1=Cx1+mp1[i][0];
        Cy1=Cy1+mp1[i][1];
        Cx2=Cx2+mp2[i][0];
        Cy2=Cy2+mp2[i][1];
     }
    Cx1=Cx1/8.0;
    Cy1=Cy1/8.0;
    Cx2=Cx2/8.0;
    Cy2=Cy2/8.0;
    for(int i=0; i<8; i++)
    {
        dis1=dis1+sqrt((mp1[i][0]-Cx1)*(mp1[i][0]-Cx1)+(mp1[i][1]-Cy1)*(mp1[i][1]-Cy1));
        dis2=dis2+sqrt((mp2[i][0]-Cx2)*(mp2[i][0]-Cx2)+(mp2[i][1]-Cy2)*(mp2[i][1]-Cy2));
    }
    dis1=dis1/8.0;
    dis2=dis2/8.0;
    double N2[3][3] = {{sqrt(2)/dis2,0,0},{0,sqrt(2)/dis2,0},{(-1*sqrt(2)/dis2)*Cx2,(-1*sqrt(2)/dis2)*Cy2,1}};
    double N1[3][3] = {{sqrt(2)/dis1,0,(-1*sqrt(2)/dis1)*Cx1},{0,sqrt(2)/dis1,(-1*sqrt(2)/dis1)*Cy1},{0,0,1}};
    double Et[3][3] = {{x[0],x[1],x[2]},{x[3],x[4],x[5]},{x[6],x[7],1}};
    double aux[3][3];
    Mul(N2,Et,aux);
    Mul(aux,N1,E);    
    for (int i = 0; i < 3; i++) 
    {
      for (int j = 0; j < 3; j++) 
      {
          E[i][j] = E[i][j]/E[2][2];          
      }
    }  
}
//Multiply two 3x3 Matrix
void Mul(double A[3][3], double B[3][3], double C[3][3])
{
    for (int a = 0; a < 3; a++) 
    {
      for (int i = 0; i < 3; i++) 
      {
          double suma = 0.0;
          for (int j = 0; j < 3; j++) 
          {
              suma = suma+ A[i][j] * B[j][a];
          }
          C[i][a] = suma;
      }
    }  
}
