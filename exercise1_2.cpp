#include<iostream>
#include <cmath>
using namespace std;

double round_to(double value, int decimals)
{
    return std::round(value * std::pow(10,decimals)) / std::pow(10,decimals);
}

double integral(int n){
    if(n==0){
        return 1-cos(1);
    }else{
        return 1-n*(n-1)*integral(n-2);
    }
}

double integral_r(int n, int decimals){
    if(n==0){
        return 1-round_to(cos(1), decimals);
    }else{
        return 1-n*(n-1)*integral_r(n-2, decimals);
    }
}

int main(){
    cout<<"\nInput value for I0:  "<<cos(1);
    cout<<"\nIntegral 8 mit genauen Eingangswerten:            "<<integral(8);
    // printf("\nIntegral 8 = %f\n", integral(8));
    cout<<"\nInput value for I0:  "<<round_to(cos(1), 5);
    cout<<"\nIntegral 8 mit Eingangswerten gerundet auf 5nkst: "<<integral_r(8, 5);
    cout<<"\nInput value for I0:  "<<round_to(cos(1), 3);
    cout<<"\nIntegral 8 mit Eingangswerten gerundet auf 3nkst: "<<integral_r(8, 3);
    return 0;
}