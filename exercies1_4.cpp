#include<iostream>
#include <cmath>
using namespace std;

int sgn(double b){
    if(b>0){
        return 1;
    }else{
        return -1;
    }
}

double x1_1(double a, double b, double c){
    return (-b+sqrt(b*b-4*a*c))/(2*a);
}
double x2_1(double a, double b, double c){
    return (-b-sqrt(b*b-4*a*c))/(2*a);
}

double x1_2(double a, double b, double c){
    return (2*c/(-b-sqrt(b*b-4*a*c)));
}

double x2_2(double a, double b, double c){
    return (2*c/(-b+sqrt(b*b-4*a*c)));
}

double x1_3(double a, double b, double c){
    double q = -0.5*(b+sgn(b)*sqrt(b*b-4*a*c));
    return c/q;
}

double x2_3(double a, double b, double c){
    double q = -0.5*(b+sgn(b)*sqrt(b*b-4*a*c));
    return q/a;
}



int main(){
    double a;
    double b;
    double c;
    cout<<"to terminate the programm set all parameters to 0\n";
    while(1){
        cout<<"Enter parameter a: ";
        cin>>a;
        cout<<"Enter parameter b: ";
        cin>>b;
        cout<<"Enter parameter c: ";
        cin>>c;
        if(a==0 && b==0 && c==0){
            break;
        }
        // cout<<"the term is: f(x)="<<a<<"x^2+"<<b<<"x+"<<c<<"\n";
        cout<<"values with formula 1: "<<x1_1(a,b,c)<<" and "<<x2_1(a,b,c)<<"\n";
        cout<<"values with formula 2: "<<x1_2(a,b,c)<<" and "<<x2_2(a,b,c)<<"\n";
        cout<<"values with formula 3: "<<x1_3(a,b,c)<<" and "<<x2_3(a,b,c)<<"\n";
    }
    return 0;
}