#include<iostream>
#include <limits>
using namespace std;
int main(){
    float testf = std::numeric_limits<float>::min();
    double testd = std::numeric_limits<double>::min();
    cout<<"single precision test value: "<<testf<<"\n";
    int countf = 0;
    int countd = 0;
    while(testf != 0){
        testf = testf/2;
        cout<<testf<<"\n";
        countf += 1;
    }
    cout<<"Number of divisions before 0: "<<countf<<"\n";
    cout<<"\ndouble precision test value: "<<testd<<"\n";
    while(testd != 0){
        testd = testd/2;
        cout<<testd<<"\n";
        countd += 1;
    }
    cout<<"Number of divisions before 0: "<<countd<<"\n";
    cout<<"press something to close";
    string hi = "ööö";
    while(1){
        cin>>hi;
        if(hi!="ööö"){
            break;
        }
    }
    return 0;
}
//der Abstand ist immer ne verdopplung der zahl - logisch, teile ja durch 2
