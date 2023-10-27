#include<iostream>
#include <cmath>
using namespace std;


double cheated_recursion(int N, double y1){
    int n = 1;
    double ym = 1;
    double yn = y1;
    while(n<N){
        y1 = yn;
        yn = ym-yn;
        ym = y1;
        n++;
    }
    return yn;
}

double recursion(int N, double y1){
    if(N==0){
        return 1;
    }else if(N == 1){
        return y1;
    }else{
        return recursion(N-2, y1)-recursion(N-1,y1);
    }
}

double potentiation(int N, double y1){
    return pow(y1, N);
}

int main(){
    int cycles;
    cout<<"Enter 1 for the positive y and !1 for negative y: ";
    int variation;
    cin>>variation;
    double y1 = (-1-sqrt(5))/2;
    double ypos = (-1+sqrt(5))/2;
    if(variation==1){
        y1 = ypos;
    }
    cout<<"y1 = "<<y1;
    cout<<"\n";
    cout<<"to terminate the programm type 0 after \"number of recursions\"";
    cout<<"\n";
    while(1){
    cout<<"Enter number of recursions: ";
    cin>>cycles;
    if(cycles == 0){
        break;
    }
    cout<<"with cheated recursionmethod: "<<cheated_recursion(cycles, y1);
    cout<<"\n";  
    if(cycles >= 40){
        cout<<"with recursionmethod:          no real recursion, because it takes to long";
        cout<<"\n";  
    }else{
        cout<<"with recursionmethod:         "<<recursion(cycles, y1);
        cout<<"\n"; 
    }   
    cout<<"with potentiationmethod:      "<<potentiation(cycles, y1);    
    cout<<"\n \n"; 
    }
    return 0;
}