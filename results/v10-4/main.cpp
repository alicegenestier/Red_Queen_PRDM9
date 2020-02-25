#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iterator>
#include <algorithm>
#include <fstream>

#include "model.h"

using namespace std;


int main(int argc, char* argv[])    {
	
	srand(time(NULL));

	cout<< "Hello World !\0"<< endl;
	cout<< " "<< endl;
	
	float temps1, temps2;
	clock_t t1, t2, t3;
	
	t1=clock();
	
	cout<< "Test for Model constructor" <<endl;
	Model model1;
	
	t2=clock();
	temps1=(float)(t2-t1)/CLOCKS_PER_SEC;
	printf("temps1 = %f\n", temps1);
	
	//cout<< "N : "<< model1.N() << endl;
	//cout<< "L : "<< model1.L() << endl;
	/*cout<< "Position of the Prdm9 gene in the genome : "<<model1.indPrdm9()<<endl;
	cout<<'\n';
	
	cout<< "populations : " << endl;
	model1.printpop(0);
	
	
	cout<< "genotypes : "<<endl;
	model1.printgen(0);
	
	cout<< "site positions for each allele"<<endl;
	model1.printposallele();
	
	cout<< "Allele for each position : "<<endl;
	model1.printallelepos();
	
	cout<<"Affinity : "<<endl;
	model1.printaffinity();*/
	
	//cout<<model1.choose(model1.L())<<endl;
	/*cout<<"test vectfreesites :"<<endl;
	vector<int> vect = {1,0,0,1,-1,1,0,-1};
	vector<int> freevect = model1.vectfreesites(vect,1);
	for (auto i : freevect){
		cout<<' '<<i;
	}
	cout<<endl;
	cout<<'\n';*/

	/*cout<<"test choosemany :"<<endl;
	vector<int> vect1 = {1,5,8,9};
	vector<int> vectsite = model1.choosemany(3, vect1);
	for (auto i : vectsite){
		cout<<' '<<i;
	}
	cout<<'\n';*/
	
	//fonction retournant nombre flottant entre 2 bornes
	//double frand_a_b(double a, double b){
		//return(rand()/(double)RAND_MAX)*(b-a)+a;
	//}
	//double c = (rand()/(double)RAND_MAX);
	//cout<<"C : "<<c<<endl
	
	// test occupiedsites function
	/*vector<int> a{0,0,-3,-3,0,-2,-3,-1,1,1};
	vector<int> b=model1.occupiedsites(a);
	for (auto i : b){
		cout<<" "<< i;
	}
	cout<<endl;*/
	
	/*cout<< "Test for sitemutation function : "<<endl;
	for (int l=1; l<11; l++){
		cout<<"Generation "<<l<<endl;
		model1.sitemutation();
		model1.printpop(2);
	}*/
	
	/*cout<< "Test for allelemutation function : "<<endl;
	for (int l=1; l<5; l++){
		cout<<"Generation "<<l<<endl;
		model1.allelemutation();
		model1.printpop(0);
	}*/
	
	/*cout<<"Test sitemutation + allelemutation + updatemissingallele"<<endl;
	for (int l=1; l<5; l++){
		cout<<"Generation "<<l<<endl;
		model1.sitemutation();
		model1.allelemutation();
		model1.updatemissingallele();
		cout<<"pop :"<<endl;
		model1.printpop(0);
		cout<<"map :"<<endl;
		model1.printposallele();
		cout<< "genotypes : "<<endl;
		model1.printgen(0);
		cout<< "Allele for each position : "<<endl;
		model1.printallelepos();
		//model1.Meiosis();
	}
	cout<<endl;*/
	
	//Test meiosis function
	/*cout<< "Test for meiosis function" <<endl;
	model1.sitemutation();
	model1.allelemutation();
	model1.updatemissingallele();
	cout<<"pop :"<<endl;
	model1.printpop(0);
	cout<<"map :"<<endl;
	model1.printposallele();
	cout<< "genotypes : "<<endl;
	model1.printgen(0);
	cout<< "Allele for each position : "<<endl;
	model1.printallelepos();
	int v = model1.Meiosis(3,1);
	cout<< "populations : " << endl;
	model1.printpop(2);
	cout<<v<<endl;
	cout<<endl;*/
	
	/*cout<<"Test fillnewpop function : "<<endl;
	model1.sitemutation();
	model1.allelemutation();
	model1.updatemissingallele();
	cout<<"pop :"<<endl;
	model1.printpop(0);
	cout<<"map :"<<endl;
	model1.printposallele();
	cout<< "genotypes : "<<endl;
	model1.printgen(0);
	cout<< "Allele for each position : "<<endl;
	model1.printallelepos();
	model1.fillnewpop();
	cout<< "populations : " << endl;
	model1.printpop(2);
	cout<<endl;
	cout<<"nb of failed meiosis : "<<model1.nbfailedmeiosis()<<endl;*/
	
	model1.manygenerations();
	
	t3=clock();
	temps2=(float)(t3-t2)/CLOCKS_PER_SEC;
	printf("temps2 = %f\n", temps2);



	/*ofstream os((name + ".trace").c_str());

    for (int gen=0; gen<max_ngen; gen++)    {

        // mutation
        //
        //
        //

        if (! gen % every)  {
            os << gen << '\t' << get_allele_number() << '\t' << get_current_diversity() << '\t'  << get_current_activity() << '\n';
            os.flush();
        }
    }*/

return 0;
}

