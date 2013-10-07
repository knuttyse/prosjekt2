#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

int main()
{
    int counter=0;
    double n=200;
    double rho_min= 0.0;
    double rho_max=5.0;
    double h = (rho_max-rho_min)/(n-1);
    vec rho = zeros(n+1);
    double faktor1=2.0/(h*h);
    double faktor2=-1.0/(h*h);
    double omega_r;
    // Brukeren må skrive inn verdien på omega_r:
    cout << "omega_r=";
    cin >> omega_r;

    mat A= mat(n,n);
    //finn rho-vektor:
    for (int i=0;i<=n;i++){
        rho(i)=rho_min+i*h;
    }
    //Regn ut A-matrisen vi starter med.
    //Merk at potensialet blir annerledes for denne oppgaven:
    A(0,0)= faktor1+omega_r*omega_r*rho(1)*rho(1);//A(1,1)= 2/h^2+rho(1)^2
    for (int i = 1;i<n;i++){
        A(i,i)=faktor1 +omega_r*omega_r*rho(i+1)*rho(i+1)+1.0/rho(i+1); //d_i
        A(i-1,i)=faktor2; // e_{i-1}
        A(i,i-1)=faktor2;// e_{i-1}
    }
    //begynn med Jacobi-metode:
    double epsilon =10e-8; // 10^-7
    int l=0;
    int k=0;
    double maks = abs(A(0,1)); // maa vaere stoerre enn epsilon!
    double tau;
    double t;

    while (maks>epsilon){
        //Finn maks element i A:
        maks=0;
        for ( int i = 0; i < n; i++ ) {
            // Siden transformasjonen er symmetrisk,
            // sjekker vi bare elementene på den ene siden av diagonalen.
            // Her har vi valgt den øvre slik at j>i
            for ( int j = i + 1; j < n; j++ ) {
                if ( abs(A(i,j)) > maks ) {
                    maks = abs(A(i,j));
                    l = j;
                    k = i;}
            }
        }
        double c;
        double s;
        if (A(k,l)!=0){
               // De trigonometriske identitetene. regner ut tau, t,c,s
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
    //Finn den laveste egenverdien
    double lambda0=1000;
    for (int i=0;i<n;i++){
        if (A(i,i)<lambda0){
            lambda0=A(i,i);
        }
    }
    cout << lambda0<<endl;
    return 0;
}
