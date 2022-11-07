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
#include <string>
#include <cstring>
#include <omp.h>
//#include <boost/math/distributions/beta.hpp>

#include "model.h"

using namespace std;


int main(int argc, char* argv[])    {
	srand(time(NULL));

	cout<< "Hello World !\0"<< endl;
	cout<< " "<< endl;


	// default values
	string name = "simu";
	int N = 1000;
	int L = 100000;
	int nbsite = 400;
	int indPrdm9 = 5;
	int nballele = 1;
	int parityIndex = 0;
    double u = 1e-4;
    double v = 1e-4;
	double w = 1e-3;
	double meanaff = 0.6;
	double varaff = 1;
	int nbDSB = 6;
	int nbGenerations = 10000;
	bool ismigration = false;
	bool zygosity = false;
	bool withDSB = false;
	int everygen = 10;
	double m = 0.1;
	double alpha = 0.5;
	double beta = 0.5;
	int nbgenmig = -1;
	int popsamesize = true;
	int nbloop = 1000;
	int nbcore = 4;
	bool isallele = true;
	bool issampling = true;
	bool isanalytic = true;
	double ctot = 1600;
	bool targetcomp = false;
	int quantilenb = 0;
	int nbmeiperind = 1;
	double cfreethreshold = 0.001;
	bool affinityUniform = 0;
	

	int i=1;
	//cout<<"argc"<<argc<<endl;
    while (i < argc)    {
    	//cout<<argv[i]<<endl;
    	string s = argv[i];

		if (s == "-N") {
        	i++;
            N = atoi(argv[i]);
        }
        else if (s == "-L") {
        	i++;
            L = atoi(argv[i]);
        }
        else if (s == "-nbsite") {
        	i++;
            nbsite = atoi(argv[i]);
        }
        else if (s == "-indPrdm9") {
        	i++;
            indPrdm9 = atoi(argv[i]);
        }
        else if (s == "-nballele") {
        	i++;
            nballele = atoi(argv[i]);
        }
        else if (s == "-parityIndex") {
        	i++;
            parityIndex = atoi(argv[i]);
        }
        else if (s == "-u")  {
        	i++;
            u = atof(argv[i]);
        }
        else if (s == "-v")  {
        	i++;
            v = atof(argv[i]);
        }
        else if (s == "-w") {
        	i++;
            w = atof(argv[i]);
        }
        else if (s == "-meanaff") {
        	i++;
            meanaff = atof(argv[i]);
        }
        else if (s == "-varaff") {
        	i++;
            varaff = atof(argv[i]);
        }
        else if (s == "-nbDSB") {
        	i++;
            nbDSB = atoi(argv[i]);
        }
        else if (s == "-nbGenerations") {
        	i++;
            nbGenerations = atoi(argv[i]);
        }
        else if (s == "-ismigration") {
        	i++;
            ismigration = atoi(argv[i]);
        }
        else if (s == "-zygosity") {
        	i++;
            zygosity = atoi(argv[i]);
        }
        else if (s == "-withDSB") {
        	i++;
            withDSB = atoi(argv[i]);
        }
        else if (s == "-everygen") {
        	i++;
            everygen = atoi(argv[i]);
        }
		else if (s == "-m") {
        	i++;
            m = atof(argv[i]);
        }
		else if (s == "-alpha") {
        	i++;
            alpha = atof(argv[i]);
        }
		else if (s == "-beta") {
        	i++;
            beta = atof(argv[i]);
        }
		else if (s == "-nbgenmig") {
        	i++;
            nbgenmig = atoi(argv[i]);
        }
		else if (s == "-popsamesize") {
        	i++;
            popsamesize = atoi(argv[i]);
        }
        else if (s == "-nbloop"){
			i++;
			nbloop = atoi(argv[i]);
		}
		else if (s == "-nbcore"){
			i++;
			nbcore = atoi(argv[i]);
		}
		else if (s == "-isallele"){
			i++;
			isallele = atoi(argv[i]);
		}
		else if (s == "-issampling"){
			i++;
			issampling = atoi(argv[i]);
		}
		else if (s == "-isanalytic"){
			i++;
			isanalytic = atoi(argv[i]);
		}
		else if (s == "-ctot"){
			i++;
			ctot = atoi(argv[i]);
		}
		else if (s == "-targetcomp"){
			i++;
			targetcomp = atoi(argv[i]);
		}
		else if (s == "-quantilenb"){
			i++;
			quantilenb = atoi(argv[i]);
		}
		else if (s == "-nbmeiperind"){
			i++;
			nbmeiperind = atoi(argv[i]);
		}
		else if (s == "-cfreethreshold"){
			i++;
			cfreethreshold = atoi(argv[i]);
		}
		else if (s == "-affinityUniform"){
			i++;
			affinityUniform = atoi(argv[i]);
		}
        else {
        	// name of the run (name will be followed by extensions: <name>.general…);
            name = argv[i];
        }
        i++;
	}
	
	
	float temps1, temps2;
	clock_t t1, t2, t3;
	
	t1=clock();
	
	cout<< "Test for Model constructor" <<endl;
	Model model1(N, L, nbsite, indPrdm9, nballele, parityIndex, v, u, w, meanaff, varaff, nbDSB, nbGenerations, ismigration, zygosity, withDSB, everygen, m, alpha, beta, nbgenmig, popsamesize, nbloop, nbcore, isallele, issampling, isanalytic, ctot, targetcomp, quantilenb, nbmeiperind, cfreethreshold, affinityUniform, name);
	
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
		model1.sitemutation((model1.populations()),(model1.genotypes()));
		model1.printpop(2,model1.populations());
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
	model1.printpop(0,model1.populations());
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
	
	model1.manygenerations();/////////////////////////////////////////////////////////////////////////////////////
	
	/*cout<<"pop1 :"<<endl;
	model1.printpop(0,model1.populations1());
	cout<<"pop2 :"<<endl;
	model1.printpop(0,model1.populations2());
	model1.migration();
	cout<<"pop1 :"<<endl;
	model1.printpop(0,model1.populations1());
	cout<<"pop2 :"<<endl;
	model1.printpop(0,model1.populations2());*/

	//model1.cfree();
	
	//test quantile loi exponentielle
	/*int nbquantile=10;
	double quantile_length=1/double(nbquantile);
	double quantile=0;
	double meangamma=1/0.44;
	for (int i=0; i<nbquantile; i++){
		quantile=quantile+quantile_length;
		cout<<"quantile "<<quantile<<" : "<<-log(1-quantile)/meangamma<<endl;
	}
	*/
	
	t3=clock();
	temps2=(float)(t3-t2)/CLOCKS_PER_SEC;
	printf("temps2 = %f\n", temps2);
	
	/*cout<<"test choosemanymigration :"<<endl;
	vector<int> vectsite = model1.choosemanymigration(15);
	for (auto i : vectsite){
		cout<<' '<<i;
	}
	cout<<'\n';*/
	

return 0;
}

