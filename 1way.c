#include<stdio.h>
#include<math.h>
#include<gsl/gsl_cdf.h>
#include<stdlib.h>
#define TEST 1

//prototype
typedef struct onewayANOVA_Experiments Experiments;
Experiments* cst_Exper(int l,int m);
void dst_Exper(Experiments* expr);
Experiments* readCSV(char* filename);
int oneWayANOVA(const Experiments* Expr,const double alpha);

struct onewayANOVA_Experiments{
    //levelArr[level][repetition],levelArr is a pointer of a pointer,not a 2Darray
    int level;
    int repetition;
    double**levelArr;
};

Experiments* cst_Exper(int l,int m){
    Experiments* temp=(Experiments*)malloc(sizeof(Experiments));
    temp->level=l;
    temp->repetition=m;
    temp->levelArr=(double**)malloc(sizeof(double*)*l);
    for(int i=0;i<l;i++)temp->levelArr[i]=(double*)malloc(sizeof(double)*m);
    return temp;
}
void dst_Exper(Experiments* expr){
    for(int i=0;i<expr->level;i++)free(expr->levelArr[i]);
    free(expr);
}

Experiments* readCSV(char* filename){
    Experiments* temp;
    int lcount=0,mcount=0,n=0,level=0,rpt=0;
    char c,record[10];
    double d;
    FILE* file=fopen(filename,"r");
    while((c=fgetc(file))!=EOF){
        if(c==','){
            lcount++;
        }
        else if(c=='\n'){
            mcount++;
        }
    }
    lcount=lcount/mcount+1;
    fseek(file,0l,0);

    if(TEST)printf("lcount:%d,mcount:%d\n",lcount,mcount);
    
    temp=cst_Exper(lcount,mcount);
    while((c=fgetc(file))!=EOF){
        if(c==','){
            record[n]='\0';
            d=atof(record);
            temp->levelArr[level++][rpt]=d;
            n=0;
        }
        else if(c=='\n'){
            record[n]='\0';
            d=atof(record);
            temp->levelArr[level][rpt++]=d;
            n=0;
            level=0;
        }
        else{
            record[n++]=c;
        }
    }
    fclose(file);
    return temp;
}



int main(void) {
    Experiments* data = readCSV("experiments.csv");
    double alpha = 0.05;
    int result;
    printf("significance level: %f\n",alpha);

    result = oneWayANOVA(data,alpha);
    if(result) printf("reject null hypothesis\n");
    else printf("rejection is failed\n");
    dst_Exper(data);
    return 0;
}

int oneWayANOVA(const Experiments* Expr,const double alpha){
    int l = Expr->level;
    int m = Expr->repetition;
    double* sumArr = (double*)malloc(sizeof(double)*l);
    double* meanArr = (double*)malloc(sizeof(double)*l);
    double totalSum, pMean, CT, St, Se, Sa, Va, Ve, F, p_value;
    totalSum = 0; St = 0; Sa = 0;
    for(int i = 0; i < l; i++) sumArr[i] = 0;
    //sumArr
    for(int i = 0; i < l; i++){
        for(int j = 0; j < m; j++) sumArr[i] += Expr->levelArr[i][j];
        if(TEST) printf("sumArr[%d]:%f ",i,sumArr[i]);
    }
    if(TEST) printf("\n");
    //meanArr
    for(int i = 0; i < l; i++){ 
        meanArr[i] = sumArr[i] / m;
        if(TEST) printf("meanArr[%d]:%f ",i,meanArr[i]);   
    }
    if(TEST)printf("\n");
    //totalSum
    for(int i=0;i<l;i++)totalSum+=sumArr[i];
    
    if(TEST)printf("totalSum:%f\n",totalSum);
    //pMean
    pMean=totalSum/l;
    
    if(TEST)printf("pMean:%f\n",pMean);
    //CT
    CT=pow(totalSum,2)/(l*m);

    if(TEST) printf("CT:%f\n",CT);
    //St
    for(int i=0;i<l;i++) {
        for(int j=0;j<m;j++) {
            St+=pow(Expr->levelArr[i][j],2);
        }
    }
    St-=CT;
    if(TEST) printf("St:%f\n",St);
    //Sa
    for(int i=0;i<l;i++) {
        Sa += pow(sumArr[i],2);
    }
    Sa = Sa/m-CT;
    if(TEST) printf("Sa:%f\n",Sa);
    //Se
    Se = St-Sa;
    if(TEST) printf("Se:%f\n",Se);
    //Va(DFa = l-1)
    Va = Sa/(l-1);
    if(TEST) printf("Va:%f\n",Va);
    //Ve(DFe = l(m-1))
    Ve = Se/(l*(m-1));
    if(TEST) printf("Ve:%f\n",Ve);
    //F
    F = Va/Ve;
    if(TEST) printf("F:%f\n",F);
    //p_value
    p_value = gsl_cdf_fdist_Q(F,l-1,l*(m-1));
    if(TEST) printf("p-value:%f\n",p_value);
    //free memory
    free(sumArr);
    free(meanArr);
    //return testing result
    if(p_value<alpha) return 1;	// reject null hypothesis
    else return 0; 		//  the rejection is failed 
}