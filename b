#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    int counter=0;
    double n=200;
    double rho_min= 5.0;
    double rho_max=10e-8;
    double h = (rho_max-rho_min)/(n-1);
    vec rho = zeros(n+1);
    double faktor1=2.0/(h*h);
    double faktor2=-1.0/(h*h);
    double omega_r;
    cout << "omega_r=";
    cin >> omega_r;

    mat A= mat(n,n);
    //finn rho-vektor:
    for (int i=0;i<=n;i++){
        rho(i)=rho_min+i*h;
    }
    A(0,0)= faktor1+omega_r*omega_r*rho(1)*rho(1);//A(1,1)= 2/h^2+rho(1)^2
    for (int i = 1;i<n;i++){
        A(i,i)=faktor1 +omega_r*omega_r*rho(i+1)*rho(i+1)+1.0/rho(i+1); //d_i
        A(i-1,i)=faktor2; // e_{i-1}
        A(i,i-1)=faktor2;// e_{i-1}
    }
    //begynn med Jacobi-metode:
    double epsilon =10e-8;
    int l=0;
    int k=0;
    double maks = abs(A(0,1)); // maa vaere stoerre enn epsilon!
    double tau;
    double t;

    while (maks>epsilon){
        //start: finn maks
        maks=0;
        for (int i=0; i<n;i++){
            for(int j=0; j<n;j++){

                if ( (i!=j) && (abs(A(i,j))>abs(maks)) ){
                    maks=abs(A(i,j));
                    if (i>j){l=i; k=j;}
                    if (i<j){l=j; k=i;}
                }
            }
        }
        //finner riktig element
        double c;
        double s;
        if (A(k,l)!=0){
               tau= ( A(l,l)-A(k,k) )/(2*A(k,l));
               if(tau<0){
                    t= -tau-sqrt(1+tau*tau);
                        }
               else{
                    t= -tau+sqrt(1+tau*tau);
                   }
               c=1/sqrt(1+t*t);
               s=t*c;
        }
        else{s=0.0;c=1.0;}

        double Akk= A(k,k);
        double All= A(l,l);
        double cc=c*c;
        double cs=c*s;
        double ss=s*s;
        A(k,k)=Akk*cc-2.0*A(k,l)*cs+All*ss;
        A(l,l)=All*cc+2.0*A(k,l)*cs+Akk*ss;
        A(k,l)=0.0;
        A(l,k)=0.0;
        for (int i=0;i<n;i++){
                if( (i!=k) && (i!=l) ){
                    double Aik= A(i,k);
                    double Ail= A(i,l);
                    A(i,k)=Aik*c-Ail*s;
                    A(k,i)=A(i,k);
                    A(i,l)=Ail*c+Aik*s;
                    A(l,i)=A(i,l);
                    }
        }
        counter++;
    }

    vec lambda=zeros(3);
    lambda(0)=10e300;
    lambda(1)=10e300;
    lambda(2)=10e300;
    for (int i=0;i<n;i++){
        if ((A(i,i)<lambda(2)) && (A(i,i)>lambda(1))){
            lambda(2)=A(i,i);
        }
        if (A(i,i)<lambda(1) && A(i,i)>lambda(0)  ){
            lambda(1)=A(i,i);}

        if (A(i,i)<lambda(0)){
            lambda(0)=A(i,i);
        }
    }
    cout << lambda(0)<<endl;
    return 0;
}
