#include <iostream>
#include <armadillo>
#include <fstream>
using namespace std;
using namespace arma;

int main()
{
    double n=1000;
    double omega_r;
    cout << "omega_r=";
    cin >> omega_r;
    double rho_min= 0.0;
    double rho_max =5.0;

    double h = (rho_max-rho_min)/(n-1);
    vec rho = zeros(n+1);
    double faktor1=2.0/(h*h);
    double faktor2=-1.0/(h*h);


    //finn rho-vektor:
    for (int i=0;i<=n;i++){
        rho(i)=rho_min+i*h;
    }

    mat A= mat(n,n);
    //Finn A matrisen:
    A(0,0)= faktor1+omega_r*omega_r*rho(1)*rho(1);//A(1,1)= 2/h^2+rho(1)^2
    for (int i = 1;i<n;i++){
        A(i,i)=faktor1 +omega_r*omega_r*rho(i+1)*rho(i+1)+1.0/rho(i+1); //d_i
        A(i-1,i)=faktor2; // e_{i-1}
        A(i,i-1)=faktor2; }

    vec eigval = zeros(n);
    mat eigvec = mat(n,n);

    eig_sym(eigval, eigvec, A);
    //Armadillo har funnet egenverdier og egenvektorer.
    //Naa skal resultatene skrives til fil:
    string filnavn;
    cout << "filnavn (eksempel 05.txt):";
    cin >> filnavn;
    //skriver ut den foerste egenverdien (egentlig oppgave c)
    cout << endl << eigval(0);
    ofstream myfile;
    myfile.open (filnavn.c_str());
    myfile << rho(0)<<"  " << 0.0 <<"  "<< 0.0<<"  " <<0.0<<endl;
    for (int i=0;i<n;i++){
        myfile << rho(i+1)<<"  " << eigvec(i,0) <<"  "<< eigvec(i,1)<<"  " <<eigvec(i,2)<<endl;
    }
    myfile << rho(n)+h<<"  " << 0.0 <<"  "<< 0.0<<"  " <<0.0<<endl;
    myfile.close();
    return 0;
}
