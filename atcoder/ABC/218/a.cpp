#include<stdio.h>
#include<atcoder/all>

int main(){
    int a,b,c;
    // scanf("%d %d %d",&a,&b,&c);
    a= 1;
    b= 2;
    c= 3;
    if(a==b && b==c){
        printf("No\n");
    }else if(a==b || b==c || c==a){
        printf("Yes\n");
    }else{
        printf("No\n");
    }
    return 0;
}