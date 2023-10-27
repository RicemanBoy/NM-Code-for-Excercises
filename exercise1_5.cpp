#include<iostream>
#include <cmath>
using namespace std;


float cheated_recursion(int N, float y1){
    int n = 1;
    float ym = 1;
    float yn = y1;
    while(n<N){
        y1 = yn;
        yn = ym-yn;
        ym = y1;
        n++;
    }
    return yn;
}

float recursion(int N, float y1){
    if(N==0){
        return 1;
    }else if(N == 1){
        return y1;
    }else{
        return recursion(N-2, y1)-recursion(N-1,y1);
    }
}

float potentiation(int N, float y1){
    return pow(y1, N);
}

int main(){
    int cycles;
    cout<<"Enter 1 for the positive y and !1 for negative y: ";
    int variation;
    cin>>variation;
    float y1 = (-1-sqrt(5))/2;
    float ypos = (-1+sqrt(5))/2;
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
    if(cycles > 40){
        cout<<"with recursionmethod:          no calculation,because it takes to long";
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


