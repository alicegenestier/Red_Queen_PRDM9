//============================
//           Includes
//============================
#include <iostream>
#include <fstream>
#include "model.h"
#include <vector>
#include <map>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <string>
#include <cstring>
//#include <boost/math/distributions/beta.hpp>
using namespace std;
//using namespace boost::math;
//============================
//         Constructors
//============================
Model::Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen, double m, double alpha, double beta, int nbgenmig, bool popsamesize,string name): N_(N), L_(L), nbsite_(nbsite), indPrdm9_(indPrdm9), nballele_(nballele), parityIndex_(parityIndex), v_(v), u_(u), w_(w), meanaff_(meanaff), varaff_(varaff), nbDSB_(nbDSB), nbGenerations_(nbGenerations), ismigration_(ismigration), zygosity_(zygosity), withDSB_(withDSB), everygen_(everygen),m_(m),alpha_(alpha),beta_(beta),nbgenmig_(nbgenmig),popsamesize_(popsamesize),name_(name) {
	
	//vector counting the number of failed meiosis per generation
	nbfailedmeiosis_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
	
	//popluations matrix
	populations_ = vector<vector<vector<int>>>(2,vector<vector<int>>(2*N_,vector<int> (L_,1)));
	// matrix initialized with 1 because at the begining all the sites are activated in the whole genome for each individual
	for (auto &i : populations_){
		for (int j=0; j<2*N_; j++)
		i[j][indPrdm9_]=0; 
	}
	
	//affinity vector
	Affinity_=vector<double>(L_,0);
	for(int i=0; i<L_; i++){
		Affinity_[i]=choosegamma(meanaff_, varaff_);
	}
	
	//genotypes vector
	genotypes_=vector<vector<int>>(2,vector<int>(2*N_,0));
	
	// Prdm9 position in the genome
	Alleleforeachpos_ = vector<int>(L_,-1);
	Alleleforeachpos_[indPrdm9_]=-2;
	
	Siteforeacheallele_[-2]={indPrdm9_};
	
	// neutral positions
	vector<int> freeposneutral = vectfreesites(Alleleforeachpos_, -1);
	vector<int> firstposneutral = choosemany(nbsite_, freeposneutral); 
	
	for(auto i : firstposneutral){
		Alleleforeachpos_[i]=-3;
	}
	
	Siteforeacheallele_[-3]=firstposneutral;
	
	//distribution beta for neutral sites
	for(int site=0; site<nbsite_; site++){
		double freqactivsite = choosebeta(alpha_, beta_);
		for(int ind=0; ind<2*N_; ind++){
			double p=bernoulli_draw(freqactivsite);
			//cout<<"p : "<<p<<endl;
			if(not p){
				populations_[parityIndex_][ind][Siteforeacheallele_[-3][site]]=0;
			}
		}
	}
	
	
	// positions of sites for first allele (0)
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> firstpos = choosemany(nbsite_, freepos); 
	
	for(auto i : firstpos){
		Alleleforeachpos_[i]=0;
	}
	
	// positions map
	for (int i = 0; i < nballele_; i++){
		Siteforeacheallele_[i]=firstpos;
		Ageallele_[i]=freqall(i, &genotypes_, &populations_);
		infoperallele_[i]={0,0,0,0,0,0};
	}
	
	/*if(ismigration_==true and nbgenmig_==0){
		populations1_=populations_;
		populations2_=populations_;
		genotypes1_=genotypes_;
		genotypes2_=genotypes_;
	}*/
}

//============================
//        Destructors
//============================
	
//============================
//           Getters
//============================
vector<vector<vector<int>>> Model::populations(){
	return populations_;
}
vector<vector<vector<int>>> Model::populations1(){
	return populations1_;
}
vector<vector<vector<int>>> Model::populations2(){
	return populations2_;
}
vector<vector<int>> Model::genotypes(){
	return genotypes_;
}
vector<vector<int>> Model::genotypes1(){
	return genotypes1_;
}
vector<vector<int>> Model::genotypes2(){
	return genotypes2_;
}
map<int,vector<int>> Model::Siteforeacheallele(){
	return Siteforeacheallele_;
}
vector<int> Model::Alleleforeachpos(){
	return Alleleforeachpos_;
}
vector<double> Model::Affinity(){
	return Affinity_;
}
int Model::parityIndex(){
	return parityIndex_;
}
int Model::N(){
	return N_;
}
int Model::L(){
	return L_;
}
int Model::indPrdm9(){
	return indPrdm9_;
}
int Model::nballele(){
	return nballele_;
}
int Model::nbsite(){
	return nbsite_;
}
double Model::u(){
	return u_;
}
double Model::v(){
	return v_;
}
double Model::m(){
	return m_;
}
double Model::meanaff(){
	return meanaff_;
} 
double Model::varaff(){
	return varaff_;
}
int Model::nbDSB(){
	return nbDSB_;
}
int Model::nbGenerations(){
	return nbGenerations_;
}
vector<vector<int>> Model::nbfailedmeiosis(){
	return nbfailedmeiosis_;
}
vector<vector<int>> Model::nbfailedmeiosis1(){
	return nbfailedmeiosis1_;
}
vector<vector<int>> Model::nbfailedmeiosis2(){
	return nbfailedmeiosis2_;
}
bool Model::zygosity(){
	return zygosity_;
}
int Model::everygen(){
	return everygen_;
}
bool Model::ismigration(){
	return ismigration_;
}
double Model::q(){
	return q_;
}
double Model::q1(){
	return q1_;
}
double Model::q2(){
	return q2_;
}
bool Model::withDSB(){
	return withDSB_;
}
double Model::w(){
	return w_;
}
string Model::name(){
	return name_;
}
map<int,double> Model::Ageallele(){
	return Ageallele_;
}
map<int,double> Model::Ageallele1(){
	return Ageallele1_;
}
map<int,double> Model::Ageallele2(){
	return Ageallele2_;
}
map<int,vector<double>> Model::infoperallele(){
	return infoperallele_;
}
map<int,vector<double>> Model::infoperallele1(){
	return infoperallele1_;
}
map<int,vector<double>> Model::infoperallele2(){
	return infoperallele2_;
}
double Model::alpha(){
	return alpha_;
}
double Model::beta(){
	return beta_;
}
int Model::nbgenmig(){
	return nbgenmig_;
}
bool Model::popsamesize(){
	return popsamesize_;
}
//============================
//           Setters
//============================

	
//============================
//           Methods
//============================
// choose method

static random_device rdev{};
static default_random_engine e{rdev()};

//uniforme
int Model::choose(int n)   {
   uniform_int_distribution<int> d(0,n-1);
   return d(e);
}

//bernoulli
int Model::bernoulli_draw(double p){
	bernoulli_distribution distribution(p);
	return distribution(e);
}

//binomiale
int Model::binomial_draw(int n, double p)  {
    binomial_distribution<int> b(n,p);
    return b(e);
}

//gamma
double Model::choosegamma(double meanaff, double varaff){
	gamma_distribution<double> g(varaff,meanaff);
	return g(e);
}

//beta
double Model::choosebeta(double alpha, double beta){
	//beta_distribution<double> mybeta(alpha, beta);
	//return mybeta(e);
	double X = choosegamma(alpha, 1);
	double Y = choosegamma(beta, 1);
	return X/(X+Y);
}

vector<int> Model::choosemany(int k, vector<int> vect){ //choose k positions among all the free positions (vect)
//return the vector of index of the chosen positions
	vector<int> newsites;
	int n = vect.size();
	try{
		if (k>n){
			throw string("Not enough free positions in the genome for new allele");
		}
	} // assert
	catch(string const& chaine){
		cerr << chaine << endl;
	}
	for (int i=0; i<k; i++){
		int index = choose(n-i);
		//boucle: tant que le find dans le vector est true, on incremnte l'indice
		bool found=true;
		while(found==true){
			vector<int>::iterator it = find(newsites.begin(),newsites.end(),vect[index]);
			if(it != newsites.end()){
				//found=true;
				index+=1;
				if(index==n){
					index=0;
				}
			}
			else{
				found=false;
			}
		}
		newsites.push_back(vect[index]);
	}
	sort(newsites.begin(), newsites.end());
	return(newsites);
}

vector<int> Model::vectfreesites(vector<int> vect, int nb){//return the index of all positions equal to nb
	vector<int> freesites;
	for(int i=0; i<vect.size(); i++){
		if(vect[i]==nb){
			freesites.push_back(i);
		}
	}
	return (freesites);	
}

vector<vector<int>> Model::occupiedsites(vector<int> vect){//return the index of all positions occupied
	vector<int> occupiedsites;
	vector<int> occupiedsitesneutral;
	for(int i=0; i<vect.size(); i++){
		if(vect[i]!= -2 and vect[i]!= -1 and vect[i]!= -3){
			occupiedsites.push_back(i);
		}else if(vect[i]==-3){
			occupiedsitesneutral.push_back(i);
		}
	}
	return (vector<vector<int>>{occupiedsites,occupiedsitesneutral});	
}

void Model::sitemutation(vector<vector<vector<int>>>* population){      ///////////// change with 2 pop : maybe give as argument a vector of pop containing either 1 or 2 pop (and vector of 1 or 2 genotype) => boucle sur les elements du vector
	//vector<vector<vector<int>>>* popul(0);
	//popul=population;
	//vector<vector<vector<int>>>copypop=population;
	// for each position with allele, prob v to mutate and if mutation, choose randomly 1 chrom to mutate
	vector<int> occupied = occupiedsites(Alleleforeachpos_)[0];
	vector<int> occupiedneutral = occupiedsites(Alleleforeachpos_)[1];
	// affichage
	/*for (auto i : occupied){
		cout<<" "<< i;
	}
	cout<<endl;*/
	vector<int> mutsites;
	vector<int> mutsitesneutral;
	for (auto i : occupied){
		if (bernoulli_draw(2*N_*v_)){
			mutsites.push_back(i);
			/*for (auto j : mutsites){
				cout<<" "<< j;
			}
			cout<<endl;*/
		}
	}
	for (auto i : occupiedneutral){
		if (bernoulli_draw(2*N_*w_)){
			mutsitesneutral.push_back(i);
			/*for (auto j : mutsites){
				cout<<" "<< j;
			}
			cout<<endl;*/
		}
	}
	for (auto j : mutsites){
		int mutchrom = choose(2*N_); //only one indiv will preform a mutation at this site
		//vector<int>::iterator itv = find(Siteforeacheallele_[-3].begin(),Siteforeacheallele_[-3].end(),j);
		//if(itv == Siteforeacheallele_[-3].end()){
		if ((*population)[parityIndex_][mutchrom][j]==1){
			/*cout<<"mutchrom : "<<mutchrom<<endl;
			cout<<"j : "<<j<<endl;*/
			(*population)[parityIndex_][mutchrom][j]=0;
		}
		//}else{
		//	popul[parityIndex_][mutchrom][j]=(popul[parityIndex_][mutchrom][j]+1)%2;
		//}
	}
	for (auto j : mutsitesneutral){
		int mutchrom = choose(2*N_);
		//vector<int>::iterator itv = find(Siteforeacheallele_[-3].begin(),Siteforeacheallele_[-3].end(),j);
		//if(itv == Siteforeacheallele_[-3].end()){
		//	if (popul[parityIndex_][mutchrom][j]==1){ 
		//		popul[parityIndex_][mutchrom][j]=0;
		//	}
		//}else{
		(*population)[parityIndex_][mutchrom][j]=((*population)[parityIndex_][mutchrom][j]+1)%2;
		//}
	}
}

void Model::allelemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,double>* Ageallele, map<int,vector<double>>* infoperallele){ //Est ce que lorsque je fais une mutation dans les 2 pop c'est le meme allele (pour le moment je dit que ce sont 2 alleles differents) + est ce que 2 alleles differents dans les pop differentes sont susceptibles d'avoir des sites en communs (je fais comme si c'etait pas le cas)
	// for each chromosome, mutation of allele with proba u and if at least one mutation update populations (change prdm9 allele at its position), genotypes(same), map (site pos associated to the new allele) and allele for each pos (change the allel corresponding to each pos : -1 -> new allele)
	vector<int> mutchromallele;
	for (int i=0; i<2*N_; i++){
		if(bernoulli_draw(u_)){
			mutchromallele.push_back(i);
		}
	}
	for (auto i : mutchromallele){
		nballele_+=1;
		vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
		/*for(auto k : freepos){//affichage freepos
			cout<<" "<<k;
		}*/
		//cout<<endl;
		vector<int> newpos = choosemany(nbsite_, freepos); 
		int newallele = nballele_-1;
		for(auto j : newpos){
			Alleleforeachpos_[j]=newallele;//update allele for each pos
		}
		Siteforeacheallele_[newallele]=newpos; //update map
		(*genotype)[parityIndex_][i]=newallele;//update genotype
		(*population)[parityIndex_][i][indPrdm9_]=newallele;//update pop
		(*Ageallele)[newallele]=freqall(newallele, genotype, population);
		(*infoperallele)[newallele]={0,0,0,0,0,0};
	}
}

void Model::updatemissingallele(){
//fonction mise a jours map si allele disparu
//comparaison entre clefs map et genotype : si un ou plusieurs alleles sont dans la map mais plus dans genotype ou pop -> on les supprime de la map et on remet toutes les pos correspondantes a -1 dans allelforeachpos et on remet toutes les activations des pos a 1
	vector<int> alleletoerase;
	for(auto &it : Siteforeacheallele_){
		//cout<<it.first<<endl;
		if(it.first!=-3 and it.first!=-2){
			//cout<<it.first<<endl;
			//cout<<&it<<endl;
			//cout<<ismigration_<<endl;
			if(not ismigration_){
				vector<int>::iterator itvect = find(genotypes_[parityIndex_].begin(),genotypes_[parityIndex_].end(),it.first);
				if(itvect == genotypes_[parityIndex_].end()){
					for(auto &i : it.second){
						//cout<<"i"<<i<<endl;
						Alleleforeachpos_[i]=-1;//free all sites
						for(int j=0; j<2*N_; j++){
							populations_[parityIndex_][j][i]=1;//reactive all sites
						}
					}
					alleletoerase.push_back(it.first);
					//Siteforeacheallele_.erase(it.first);//erase in the map
				}
			}else{
				vector<int>::iterator itvect1 = find(genotypes1_[parityIndex_].begin(),genotypes1_[parityIndex_].end(),it.first);
				vector<int>::iterator itvect2 = find(genotypes2_[parityIndex_].begin(),genotypes2_[parityIndex_].end(),it.first);
				if(itvect1 == genotypes1_[parityIndex_].end() and itvect2 == genotypes2_[parityIndex_].end()){
					for(auto &i : it.second){
						//cout<<"i"<<i<<endl;
						Alleleforeachpos_[i]=-1;//free all sites
						for(int j=0; j<2*N_; j++){
							populations1_[parityIndex_][j][i]=1;//reactive all sites
							populations2_[parityIndex_][j][i]=1;//reactive all sites
						}
					}
					alleletoerase.push_back(it.first);
					//Siteforeacheallele_.erase(it.first);//erase in the map
				}
			}
		}
	}
	for (auto a : alleletoerase){
		if(not ismigration_){
			Siteforeacheallele_.erase(a);
			Ageallele_.erase(a);
			infoperallele_.erase(a);
		}else{
			Siteforeacheallele_.erase(a);
			Ageallele1_.erase(a);
			infoperallele1_.erase(a);
			Ageallele2_.erase(a);
			infoperallele2_.erase(a);
		}
	}
}

//Print functions

//print populations
void Model::printpop(int n, vector<vector<vector<int>>>population){//n = parityIndexu (0 ou 1) ou 2 si plot 2 pop
	if(n==2){
		for (auto i : population){ // a changer
			for (auto j : i){
				for (auto k : j){
					cout<<' '<<k;
				}
				cout<<'\n';
			}	
			cout<<'\n';
		}
		cout<<'\n';
	}
	else if(n==0 or n==1){
		vector<vector<int>> pop = population[n]; // a changer
		for (auto i : pop){
			for (auto j : i){
				cout<<' '<<j;
			}
			cout<<'\n';
		}
		cout<<'\n';
	}
}

//print genotypes
void Model::printgen(int n, vector<vector<int>>genotype){
	if(n==2){
		for (auto i : genotype){
			for (auto j : i){
				cout<<' '<<j;
			}
			cout<<'\n';
		}
		cout<<'\n';
	}
	else if(n==0 or n==1){
		vector<int> gen = genotype[n];
		for (auto i : gen){
			cout<<' '<<i;
		}
		cout<<'\n';
	}
}

//print site positions for each allele
void Model::printposallele(){
	for (auto const &it : Siteforeacheallele_){
		cout<<it.first<<"=>";
		for(auto const &i : it.second){
			cout<<" "<<i;
		}
		cout<<endl;
	}
	cout<<'\n';
}

//print Allele for each position
void Model::printallelepos(){
		for (auto i : Alleleforeachpos_){
		cout<<' '<<i;
	}
	cout<<'\n'<<endl;
	
}

//print Affinity
void Model::printaffinity(){
		for (auto i : Affinity_){
		cout<<' '<<i;
	}
	cout<<'\n'<<endl;
}

//print age allele
void Model::printageallele(map<int,double>* Ageallele){
	for (auto const &it : (*Ageallele)){
		cout<<it.first<<" => "<<it.second<<endl;
	}
	cout<<'\n';
}

//print info per allele
void Model::printinfoallele(map<int,vector<double>>* infoperallele){
	for (auto const &it : (*infoperallele)){
		cout<<it.first<<" => ";
		for(auto const &i : it.second){
			cout<<" "<<i;
		}
		cout<<endl;
	}
	cout<<'\n';
}


// Meiosis
int Model::Meiosis(int no_chrom_ind, int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, vector<vector<int>>* nbfailedmeiosis, double* q){ //
	// choose one ind
	int indiv = 2*choose(N_);
	//cout<<"indiv "<<indiv<<endl;
	//homozygote or heterozygote
	int ind_gen=2;
	vector<int> zygote{(*genotype)[parityIndex_][indiv]};
	/*for(auto i : genotypes_){
		for (auto j : i){cout<<"j"<<j;}
	}*/
	if((*genotype)[parityIndex_][indiv]!=(*genotype)[parityIndex_][indiv+1]){
		if(zygosity_){
			ind_gen=1;
		}
		zygote.push_back((*genotype)[parityIndex_][indiv+1]);
	}
	//cout<<"ind_gen : "<<ind_gen<<endl;
	//cout<<"zygote : ";
	/*for (auto i : zygote){
		cout<<" "<<i;
	}
	cout<<endl;*/
	
	//PRDM9 binding
	vector<vector<int>>summarysites;
	vector<int>Z;
	vector<int> vectsites;
	int nblinksite = 0;
	for(auto z : zygote){
		/*cout<<"z : "<<z<<endl;
		cout<<"sites and affinity"<<endl;
		for(auto sitez : Siteforeacheallele_[z]){
			cout<<sitez<<" affinity : "<<Affinity_[sitez]<<endl;
		}
		cout<<"zygot z :"<<z<<endl;*/
		for(auto i : Siteforeacheallele_[z]){
			//cout<<"i : "<<i<<endl;
			vector<int>linkedsites=vector<int>(5,0);
			linkedsites[0]=i;
			double p_occup=ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]);
			//cout<<"proba occup : "<<p_occup<<endl;
			//cout<<"p_occup : "<<p_occup<<endl;
			bool islinked=false;
			for(int j=0; j<4; j++){
				//cout<<"+";
				int chrom=indiv+j/2;
				if((*population)[parityIndex_][chrom][i]==1){ // a changer
					if(bernoulli_draw(p_occup)){
						linkedsites[j+1]=1;
						islinked=true;
						nblinksite+=1;
					}
				}
			}
			if(islinked==true){
				vectsites.push_back(i);
				summarysites.push_back(linkedsites);
				Z.push_back(z);
			}
		}
	}
	//nblinksite=summarysites.size();
	//cout<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
	/*cout<<"summary sites"<<endl;
	for(auto i : summarysites){
		for(auto j : i){
			cout<<' '<<j;
		}
		cout<<"\n";
	}
	cout<<endl;*/
	//cout<<"nblinksite "<<nblinksite<<endl;
	//DSB
	double pDSB=double(nbDSB_)/nblinksite;
	//cout<<"pDSB : "<<pDSB<<endl;
	//cout<<"nblinksite : "<<nblinksite<<endl;
	//cout<<"nb linked sites : "<<nblinksite<<endl;
	//cout<<"pDSB : "<<pDSB<<endl;
	vector<vector<int>>vectsitedsb;
	vector<double>alleleDSB; ////peut servir pour q per allele 
	int checknblink =0;
	vector<vector<int>> vect_CO;
	vector<double>alleleCO; ////peut servir pour q per allele 
	for(int i=0; i<summarysites.size(); i++){
		//proba DSB
		vector<int> link = vectfreesites(vector<int>(summarysites[i].begin()+1, summarysites[i].end()), 1);
		int nbdsbpersite=0;
		/*cout<<"linked sites : ";
		for(auto l : link){
			cout<<" "<<l;
		}
		cout<<endl;*/
		bool dsb=false;
		for(auto j : link){
			checknblink++;
			//cout<<"pDSB : "<<pDSB_<<endl;
			if(bernoulli_draw(pDSB)){
				dsb=true;
				summarysites[i][j+1]=2;// DSB -> 2
				nbdsbpersite+=1;
				vectsitedsb.push_back({i,j});
				alleleDSB.push_back(Z[i]); ////peut servir pour q per allele 
				/*
				cout<<"summarysites : "<<endl;
				for(auto i : summarysites){
					for(auto j : i){
						cout<<' '<<j;
					}
					cout<<"\n";
				}
				cout<<endl;
				cout<<"a1 "<<nbdsbpersite<<endl;
				*/
				try{
					if (nbdsbpersite>1 and withDSB_){
						(*nbfailedmeiosis)[nb_gen][0]+=1; 
						(*infoperallele)[Z[i]][0]+=1;
						(*infoperallele)[Z[i]][1]+=1;
						cout<<"case 2 DSB"<<endl;
						throw int(0);
					}
				} // assert
				catch(int const& error_nb){
					if(error_nb==0){
						//cerr << "2 DSB on one site" << endl;
					}
					return-1;
				}
			}
			//cout<<"a2"<<nbdsbpersite<<endl;
		}
		/*cout<<" vectsitedsb : "<<endl;
		for(auto i : vectsitedsb){
			for(auto j : i){
				cout<<' '<<j;
			}
			cout<<"\n";
		}
		cout<<endl;*/
		vector<int> vco;
		if (dsb){ /// pose pb puisqu'on peut avoir plusieurs dsb sur meme site donc .back marche pas...
			for(int indexnbdsb = 0; indexnbdsb<nbdsbpersite; indexnbdsb++){
				vco = {summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]};
				if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==0 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==1){
					if((summarysites[i][3]==1 or summarysites[i][3]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=2){/// in these cases, I suppose that 2 DSB can perform a CO
						//cout<<"cas1"<<endl;
						//vect_CO.push_back({summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1],2});
						vco.push_back(2);
					}
					if((summarysites[i][4]==1 or summarysites[i][4]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=3){
						//cout<<"cas2"<<endl;
						//vect_CO.push_back({summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1],3});
						vco.push_back(3);
					}
				}else if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==2 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==3 or summarysites[i][1]==2){
					if((summarysites[i][1]==1 or summarysites[i][1]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=0 ){
						//cout<<"cas3"<<endl;
						//vect_CO.push_back({summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1],0});
						vco.push_back(0);
					}
					if((summarysites[i][2]==1 or summarysites[i][2]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=1){
						//cout<<"cas4"<<endl;
						//vect_CO.push_back({summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1],1});
						vco.push_back(1);
					}
				}
				if(vco.size()>2){
					vect_CO.push_back(vco);
					alleleCO.push_back(Z[i]); ////peut servir pour q per allele 
				}
			}
		}
		/*cout<<"summarysites : "<<endl;
		for(auto i : summarysites){
			for(auto j : i){
				cout<<' '<<j;
			}
			cout<<"\n";
		}
		cout<<endl;
		cout<<"vect_CO : "<<endl;
		for(auto i : vect_CO){
			for(auto j : i){
				cout<<' '<<j;
			}
			cout<<"\n";
		}
		cout<<endl;*/
	}
	try{
		if(vectsitedsb.size()==0){ // soit je fais ca mais dans ces cs la je suppose directement qu'il y a eu au moins une liaison sur chaque chromosome
//sinon, je peux stocker les noms des alleles qui ont ete lies et je n'incremente dans info allele que pour ceux qui ont eu une liaison et je peux aussi rajouter une erreure qui serai pas de liaison pour cet allele...
			(*nbfailedmeiosis)[nb_gen][1]+=1;
			if(zygote.size()==1){
				(*infoperallele)[zygote[0]][0]+=2;
				(*infoperallele)[zygote[0]][2]+=2;
			}else if(zygote.size()==2){
				(*infoperallele)[zygote[0]][0]+=1;
				(*infoperallele)[zygote[0]][2]+=1;
				(*infoperallele)[zygote[1]][0]+=1;
				(*infoperallele)[zygote[1]][2]+=1;
			}
			//cout<<"case No DSB"<<endl;
			throw int(1);
		}else if(not vect_CO.size()){
			(*nbfailedmeiosis)[nb_gen][2]+=1;
			if(zygote.size()==1){
				(*infoperallele)[zygote[0]][0]+=2;
				(*infoperallele)[zygote[0]][3]+=2;
			}else if(zygote.size()==2){
				(*infoperallele)[zygote[0]][0]+=1;
				(*infoperallele)[zygote[0]][3]+=1;
				(*infoperallele)[zygote[1]][0]+=1;
				(*infoperallele)[zygote[1]][3]+=1;
			}
			//cout<<"case No sym"<<endl;
			throw int(2);
		}
	} // assert
	catch(int const& error_nb){
		if(error_nb==1){
			//cerr << "No DSB" << endl;
		}else if(error_nb==2){
			//cerr << "No symmetrical sites (binding + DSB)" << endl;
		}
		return-1;
	}
	//cout<<"checknblink : "<<checknblink<<endl;
	/*cout<<"summarysites : "<<endl;
	for(auto i : summarysites){
		for(auto j : i){
			cout<<' '<<j;
		}
		cout<<"\n";
	}
	cout<<endl;
	cout<<"vect_CO : "<<endl;
		for(auto i : vect_CO){
			for(auto j : i){
				cout<<' '<<j;
			}
			cout<<"\n";
		}
	cout<<endl;
	cout<<"vectsitedsb : "<<endl;
	for(auto i : vectsitedsb){
		for(auto j : i){
			cout<<' '<<j;
		}
		cout<<"\n";
	}
	cout<<endl;*/
	//cout<<vect_CO.size()<<endl;
	//cout<<nblinksite<<endl;
	*q=*q+double(vect_CO.size())/(vectsitedsb.size()); // suppose qu'il n'y a pas eu d'erreur avant
	if(zygote.size()==1){
		(*infoperallele)[zygote[0]][4]+=2*(double(alleleCO.size())/alleleDSB.size());
		//cout<< "zygote : " << zygote[0] << "; alleleCO.size() : " << alleleCO.size() << "; alleleDSB.size() : "<< alleleDSB.size()<<"; infoperallele_[zygote[0]][4] : " << infoperallele_[zygote[0]][4] << endl;
		(*infoperallele)[zygote[0]][5]+=2;
	}else if(zygote.size()==2){
		(*infoperallele)[zygote[0]][4]+=(double(alleleCO.size())/alleleDSB.size());
		(*infoperallele)[zygote[1]][4]+=(double(alleleCO.size())/alleleDSB.size());
		(*infoperallele)[zygote[0]][5]+=1;
		(*infoperallele)[zygote[1]][5]+=1;
		//cout<< "zygote : " << zygote[0] << "; alleleCO.size() : " << alleleCO.size() << "; alleleDSB.size() : "<< alleleDSB.size()<<"; infoperallele_[zygote[0]][4] : " << infoperallele_[zygote[0]][4] << endl;
		//cout<< "zygote : " << zygote[1] << "; alleleCO.size() : " << alleleCO.size() << "; alleleDSB.size() : "<< alleleDSB.size()<<"; infoperallele_[zygote[1]][4] : " << infoperallele_[zygote[1]][4] << endl;
	}
	//cout<<q_<<endl;
	vector<int> index_CO;
	int choosevect=choose(vect_CO.size());
	if(vect_CO[choosevect].size()==4){
		if(bernoulli_draw(0.5)){
			index_CO={vect_CO[choosevect][0],vect_CO[choosevect][1],vect_CO[choosevect][2]};
		}else{
			index_CO={vect_CO[choosevect][0],vect_CO[choosevect][1],vect_CO[choosevect][3]};
		}
	}else{
		index_CO=vect_CO[choosevect];	
	}
	/*cout<<"index_CO : "<<endl;
	for(auto i : index_CO){
		cout<<' '<<i;
	}
	cout<<'\n';*/
	
	
	// si besoin chercher dans anciencode.cpp
	
	//toute cette partie avec copy
	
	//puis si meiose se fait :
	//tirage 1 parmi 4 : choisi le gamete a mettre dans la nouvelle pop
	int no_chromatide=choose(4);
	//cout<<"no chromatide"<<no_chromatide<<endl;
	int no_current_chrom = indiv;
	int no_homologue_chrom = indiv+1;
	if(no_chromatide==2 or no_chromatide==3){
		no_current_chrom=indiv+1;
		no_homologue_chrom=indiv;
	}
	if(no_chromatide==index_CO[1] or no_chromatide==index_CO[2]){
		//cout<<"chromatide with CO"<<endl;
		copy( (*population)[parityIndex_][no_current_chrom].begin(), (*population)[parityIndex_][no_current_chrom].begin()+index_CO[0], (*population)[(parityIndex_+1)%2][no_chrom_ind].begin() );
		copy( (*population)[parityIndex_][no_homologue_chrom].begin()+index_CO[0], (*population)[parityIndex_][no_homologue_chrom].end(), (*population)[(parityIndex_+1)%2][no_chrom_ind].begin()+index_CO[0] );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide and summarysites[index_dsb[0]][0]<index_CO[0]){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}else if((index_dsb[1]==index_CO[1] or index_dsb[1]==index_CO[2]) and summarysites[index_dsb[0]][0]>index_CO[0]){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_current_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}else{
		//cout<<"chromatide without CO"<<endl;
		copy( (*population)[parityIndex_][no_current_chrom].begin(), (*population)[parityIndex_][no_current_chrom].end(), (*population)[(parityIndex_+1)%2][no_chrom_ind].begin() );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}
	
	//return nb meiose echouee?? Pour le moment, return 0 si eurreur qqpart et 1 si meiose reussie.
	return 1;
}


//methode remplissage nouvelle pop
void Model::fillnewpop(int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, vector<vector<int>>* nbfailedmeiosis, double* q){
	//for(int indnewpop=0; indnewpop<2*N_; indnewpop++){
	for(int indnewpop=0; indnewpop<2*N_; indnewpop++){
		int meiosisState = Meiosis(indnewpop, nb_gen, population, genotype, infoperallele, nbfailedmeiosis, q);//
		while (meiosisState==-1){
			meiosisState = Meiosis(indnewpop, nb_gen, population, genotype, infoperallele, nbfailedmeiosis, q);//
			(*nbfailedmeiosis)[nb_gen][3]+=1;
		}
		//cout<<"-----------------------------------------------------------------------------------------------------"<<endl;
	}
	for(int indpop=0; indpop<2*N_;indpop++){
		(*genotype)[(parityIndex_+1)%2][indpop]=(*population)[(parityIndex_+1)%2][indpop][indPrdm9_];
	}
	parityIndex_=(parityIndex_+1)%2;
}

//methode qui repete tout ce au'on vient de faire pendant X generations
void Model::manygenerations(){
	ofstream generalfile ((name_+".trace").c_str());
	generalfile << "Generation_number" << '\t' << "Total_number_of_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" <<'\n';
    	generalfile.flush();
    	ofstream allelefile ((name_+".allele").c_str());
	allelefile << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" <<'\n';
    	allelefile.flush();
	ofstream paramsfile ((name_+".params").c_str());
	paramsfile << "N" << '\t' << N_ << '\n' << "L" << '\t' << L_ << '\n' << "nbsite" << '\t' << nbsite_ << '\n' << "indPrdm9" << '\t' << indPrdm9_ << '\n' << "nballele" << '\t' << nballele_ << '\n' << "parityIndex" << '\t' << parityIndex_ << '\n' << "u" << '\t' << u_ << '\n' << "v" << '\t' << v_ << '\n' << "w" << '\t' << w_ << '\n' << "meanaff" << '\t' << meanaff_ << '\n' << "varaff" << '\t' << varaff_ << '\n' << "nbDSB" << '\t' << nbDSB_ << '\n' << "nbGenerations" << '\t' << nbGenerations_ << '\n' << "ismigration" << '\t' << ismigration_ << '\n' << "zygosity" << '\t' << zygosity_ << '\n' << "withDSB" << '\t' << withDSB_ << '\n' << "everygen" << '\t' << everygen_ << '\n';
    	paramsfile.flush();
	
	for(int indgeneration=0; indgeneration<nbGenerations_; indgeneration++){
		clock_t t1, t2;
		t1=clock();
		//cout<<"parity index : "<<parityIndex_<<endl;
		cout<<"generation "<<indgeneration<<endl;
		// remise a 0 vector de chaque allele de infoperallele
		for (auto const &it : infoperallele_){//////////////////////// for both pop
			infoperallele_[it.first]=vector<double>{0,0,0,0,0,0};
		}	
		sitemutation(&populations_);///////////////////////for both pop
		allelemutation(&populations_, &genotypes_, &Ageallele_, &infoperallele_);///////////////////////////for both pop
		updatemissingallele();
		/*cout<<"pop :"<<endl;
		printpop(parityIndex_, &populations_);
		cout<<"map :"<<endl;
		printposallele();*/
		/*cout<<"age allele :"<<endl;
		printageallele(&Ageallele);*/
		/*cout<< "genotypes : "<<endl;
		printgen(parityIndex_,&genotypes_);*/
		/*cout<< "Allele for each position : "<<endl;
		printallelepos();*/
		/*cout<<"info per allele : "<<endl;
		printinfoallele(&infoperallele_);*/
		/////////////////////////////////////// when we want to split the pop : set ismigration to true, copy populations_ in populations1_ and populations2_, and then fillnwpop for both pop and then do all the function for both pop for the following generations
		/*if(indgeneration==nbgenmig_){
			ismigration_=true;
			populations1_=populations_;
			populations2_=populations_;
			genotypes1_=genotypes_;
			genotypes2_=genotypes_;
			Ageallele1_=Ageallele_;
			Ageallele2_=Ageallele_;
			infoperallele1_=infoperallele_;
			infoperallele2_=infoperallele_;
			//initialiser fichier generalfile et allelefile 1 et 2
			ofstream generalfile1 ((name_+"_1.trace").c_str());
			generalfile1 << "Generation_number" << '\t' << "Total_number_of_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" <<'\n';
    			generalfile1.flush();
    			ofstream allelefile1 ((name_+"_1.allele").c_str());
			allelefile1 << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" <<'\n';
    			allelefile1.flush();
			ofstream generalfile2 ((name_+"_2.trace").c_str());
			generalfile2 << "Generation_number" << '\t' << "Total_number_of_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" <<'\n';
    			generalfile2.flush();
    			ofstream allelefile2 ((name_+"_2.allele").c_str());
			allelefile2 << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" <<'\n';
    			allelefile2.flush();
		}*/
		q_=0;//////////////////////for both pop
		fillnewpop(indgeneration, &populations_, &genotypes_, &infoperallele_, &nbfailedmeiosis_, &q_);//////////////////////for both pop
		q_=q_/(2*N_);//////////////////////for both pop
		//cout<<q_<<endl;
		updatemissingallele();
		for (auto &it : Ageallele_){//////////////////////for both pop
			it.second=it.second+freqall(it.first, &genotypes_, &populations_);
		}
		for (auto &it : infoperallele_){//////////////////////for both pop
			it.second[4]=it.second[4]/it.second[5];
		}
		/*cout<<"---------"<<endl;
		cout<<"pop :"<<endl;
		printpop(parityIndex_, &popuations_);
		cout<<"map :"<<endl;
		printposallele();*/
		/*cout<<"age allele :"<<endl;
		printageallele(&Ageallele_);*/
		/*cout<< "genotypes : "<<endl;
		printgen(parityIndex_,&genotypes_);*/
		/*cout<< "Allele for each position : "<<endl;
		printallelepos();*/
		/*cout<<"info per allele : "<<endl;
		printinfoallele(&infoperallele_);*/
		//cout<<"a "<<indgeneration % everygen_<<endl;
		t2=clock();
		if (indgeneration % everygen_ ==0)  {//////////////////////for both pop
			//cout<<"b "<<indgeneration % everygen_<<endl;
        		generalfile << indgeneration << '\t' << get_allele_number() << '\t' << get_current_diversity(&genotypes_) << '\t'  << get_current_activity(&genotypes_, &populations_) << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << 1-(double(nbfailedmeiosis_[indgeneration][3])/(2*N_+nbfailedmeiosis_[indgeneration][3]))<< '\t' << double(nbfailedmeiosis_[indgeneration][0])/(2*N_+nbfailedmeiosis_[indgeneration][3]) << '\t'<< double(nbfailedmeiosis_[indgeneration][1])/(2*N_+nbfailedmeiosis_[indgeneration][3]) << '\t' << double(nbfailedmeiosis_[indgeneration][2])/(2*N_+nbfailedmeiosis_[indgeneration][3]) << '\t' << q_ <<'\n'; //ajouter taux meiose echouee du a DSB, no sym...
            		generalfile.flush();
            		for (auto const &it : Siteforeacheallele_){
	    			if(it.first!=-2){
            				allelefile << indgeneration << '\t' << it.first << '\t' << freqall(it.first, &genotypes_, &populations_) << '\t'  << actall(it.first, &populations_) << '\t' << get_age_allele(it.first, &Ageallele_) << '\t' << get_info_allele(it.first, &infoperallele_) << '\n';
            				allelefile.flush();
            			}
            		}
        	}
	}
	generalfile.close();
	allelefile.close();
	paramsfile.close();
}

int Model::get_allele_number(){
	return Siteforeacheallele_.size()-2;
}

double Model::get_age_allele(int allname, map<int,double>* Ageallele){
	if(allname==-3){
		return 0;
	}else{
		return (*Ageallele)[allname];
	}
}

double Model::get_info_allele(int allname, map<int,vector<double>>* infoperallele){
	if(allname==-3){
		return 0;
	}else{
		return (*infoperallele)[allname][4];
	}
}

double Model::freqallele(int allelename, vector<vector<int>>* genotype){
	/*if(allelename==-3){
		return double(0);
	}*/
	//if(allelename==0){
		//cout<<"count : "<<count(genotypes_[parityIndex_].begin(), genotypes_[parityIndex_].end(), allelename)-1<<endl;
		//cout<<"size genotype : "<<genotypes_[parityIndex_].size()-1<<endl;
		//return double(count(genotypes_[parityIndex_].begin(), genotypes_[parityIndex_].end(), allelename)-1)/(genotypes_[parityIndex_].size()-1); /// ????
	//}else{
		//cout<<allelename<<endl;
		//cout<<double (count(genotypes_[parityIndex_].begin(), genotypes_[parityIndex_].end(), allelename))/(genotypes_[parityIndex_].size())<<endl;
		return double (count((*genotype)[parityIndex_].begin(), (*genotype)[parityIndex_].end(), allelename))/((*genotype)[parityIndex_].size());
	//}
}

double Model::get_current_diversity(vector<vector<int>>* genotype){
	double sumfreq=0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			//cout<<"freqallele : "<<freqallele(it.first)<<endl;
			sumfreq+=(freqallele(it.first, genotype))*(freqallele(it.first, genotype));
		}
	}
	return 1/sumfreq;
}

double Model::activitymoyallele(int allele,  vector<vector<vector<int>>>* population){
	double moyact=0;
	//cout<<"allele : "<<allele<<endl;
	for(auto all : Siteforeacheallele_[allele]){
		double moyactsite=0;
		for(int i=0; i<2*N_; i++){
			if((*population)[parityIndex_][i][all]==1){
				moyactsite+=1;			
			}
		}
		moyact+=moyactsite/(2*N_);
		//cout<<allele<<"site "<<all<<"moyactsite "<<moyactsite<<"moyactsite/2N "<<moyactsite/(2*N_)<<endl;
		//moyact+=Affinity_[all];
		//cout<<"all : "<<all<<endl;
		//cout<<"Affinity_[all] : "<<Affinity_[all]<<endl;
		/*if(allele==-3){
			res=res+2*moyact*(1-moyact);
		}*/
	}
	//cout<<"moyact : "<<moyact<<endl;
	//cout<<"moyact/nbsite_ : "<<moyact/nbsite_<<endl;
	moyact=moyact/nbsite_;
	return (moyact);
}

vector<double> Model::freqneutral(vector<vector<vector<int>>>* population){
	double moyfreq=0;
	double moy2f=0;
	for(auto all : Siteforeacheallele_[-3]){
		double moyactsite=0;
		for(int i=0; i<2*N_; i++){
			if((*population)[parityIndex_][i][all]==1){
				moyactsite+=1;			
			}
		}
		double freq=moyactsite/(2*N_);
		moyfreq+=freq;
		double twof=2*freq*(1-freq);
		moy2f+=twof;
	}
	moyfreq=moyfreq/nbsite_;
	moy2f=moy2f/nbsite_;
	vector<double> vectneutral {moyfreq, moy2f};
	return vectneutral;
}

double Model::freqall(int allele, vector<vector<int>>* genotype, vector<vector<vector<int>>>* population){
	if (allele==-3){
		return freqneutral(population)[0];
	}else{
		return freqallele(allele, genotype);
	}
}

double Model::actall(int allele,  vector<vector<vector<int>>>* population){
	if (allele==-3){
		return freqneutral(population)[1];
	}else{
		return activitymoyallele(allele, population);
	}
}

double Model::get_current_activity(vector<vector<int>>* genotype, vector<vector<vector<int>>>* population){
	double moytotact=0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			//cout<<"it.first : "<<it.first<<endl;
			//cout<<"activitymoyallele(it.first) : "<<activitymoyallele(it.first)<<endl;
			moytotact+=activitymoyallele(it.first, population)*freqallele(it.first, genotype);
			//cout<<"d"<<endl;
		}
	}
	//cout<<"activity moy"<<moytotact/meanaff_<<endl;
	return moytotact;
}

vector<int> Model::choosemanymigration(int k){ //choose k individuals in the pop
//return the vector of index of the chosen positions
	vector<int> newsites;
	try{
		if (k>N_){
			throw string("To much migrants");
		}
	} // assert
	catch(string const& chaine){
		cerr << chaine << endl;
		return vector<int>{-1};
	}
	for (int i=0; i<k; i++){
		int index = choose(N_-i);
		//boucle: tant que le find dans le vector est true, on incremnte l'indice
		bool found=true;
		while(found==true){
			vector<int>::iterator it = find(newsites.begin(),newsites.end(),index);
			if(it != newsites.end()){
				//found=true;
				index+=1;
				if(index==N_){
					index=0;
				}
			}
			else{
				found=false;
			}
		}
		newsites.push_back(index);
	}
	sort(newsites.begin(), newsites.end());
	return(newsites);
}

double Model::q_two_hap(vector<int> haplotype1, vector<int> haplotype2){
	//give two haplotype (2 ind of the pop ??)
}

double Model::q_hom(int allele, vector<int> haplotype1, vector<int> haplotype2){
	//
}

double Model::q_hete(int allele1, int allele2, vector<int> haplotype1, vector<int> haplotype2){
	//
}

void Model::migration(){
	//cout<<"m_*N_ :"<< m_*N_ <<endl;
	vector<int> migrated_indiv_pop1 = choosemanymigration(int(m_*N_));
	/*cout<<"migrated_indiv_pop1 : "<<endl;
	for (auto i : migrated_indiv_pop1){
		cout<<i<<" ";
	}
	cout<<endl;*/
	vector<int> migrated_indiv_pop2 = choosemanymigration(int(m_*N_));
	/*cout<<"migrated_indiv_pop2 : "<<endl;
	for (auto i : migrated_indiv_pop2){
		cout<<i<<" ";
	}
	cout<<endl;*/
	if(migrated_indiv_pop1!=vector<int>{-1} and migrated_indiv_pop2!=vector<int>{-1}){
		for (int indiv=0; indiv<migrated_indiv_pop1.size(); indiv++){
			vector<vector<int>> ind_pop1_to_pop2 = vector<vector<int>> {populations1_[parityIndex_][2*migrated_indiv_pop1[indiv]],populations1_[parityIndex_][2*migrated_indiv_pop1[indiv]+1]};
			vector<int> ind_gen1_to_gen2 = vector<int> {genotypes1_[parityIndex_][2*migrated_indiv_pop1[indiv]],genotypes1_[parityIndex_][2*migrated_indiv_pop1[indiv]+1]};
			populations1_[parityIndex_][2*migrated_indiv_pop1[indiv]]=populations2_[parityIndex_][2*migrated_indiv_pop2[indiv]];
			populations1_[parityIndex_][2*migrated_indiv_pop1[indiv]+1]=populations2_[parityIndex_][2*migrated_indiv_pop2[indiv]+1];
			genotypes1_[parityIndex_][2*migrated_indiv_pop1[indiv]]=genotypes2_[parityIndex_][2*migrated_indiv_pop2[indiv]];
			genotypes1_[parityIndex_][2*migrated_indiv_pop1[indiv]+1]=genotypes2_[parityIndex_][2*migrated_indiv_pop2[indiv]+1];
			populations2_[parityIndex_][2*migrated_indiv_pop2[indiv]]=ind_pop1_to_pop2[0];
			populations2_[parityIndex_][2*migrated_indiv_pop2[indiv]+1]=ind_pop1_to_pop2[1];
			genotypes2_[parityIndex_][2*migrated_indiv_pop2[indiv]]=ind_gen1_to_gen2[0];
			genotypes2_[parityIndex_][2*migrated_indiv_pop2[indiv]+1]=ind_gen1_to_gen2[1];
		}
	}
}

