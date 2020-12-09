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
Model::Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen, double m, double alpha, double beta, int nbgenmig, bool popsamesize, int nbloop, string name): N_(N), L_(L), nbsite_(nbsite), indPrdm9_(indPrdm9), nballele_(nballele), parityIndex_(parityIndex), v_(v), u_(u), w_(w), meanaff_(meanaff), varaff_(varaff), nbDSB_(nbDSB), nbGenerations_(nbGenerations), ismigration_(ismigration), zygosity_(zygosity), withDSB_(withDSB), everygen_(everygen),m_(m),alpha_(alpha),beta_(beta),nbgenmig_(nbgenmig),popsamesize_(popsamesize),nbloop_(nbloop),name_(name) {
	
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
		infoperallele_[i]={0,0,0,0,0,0,get_mean_affinity(i,&populations_)}; // infoperallele{nb failed meiosis, 2 DSB on 1 site, no DSB, no symetrical site, q, nb meiosis with at least 1 DSB, mean affinity of the allele}
	}
	q_=0;
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
int Model::nbloop(){
	return nbloop_;
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
		bool found=true;
		while(found==true){
			vector<int>::iterator it = find(newsites.begin(),newsites.end(),vect[index]);
			if(it != newsites.end()){
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

vector<vector<int>> Model::occupiedsites(vector<int> vect, vector<vector<int>>* genotype){//return the index of all positions occupied
	vector<int> occupsites;
	vector<int> occupsitesneutral;
	for(int i=0; i<vect.size(); i++){
		if(vect[i]!= -2 and vect[i]!= -1 and vect[i]!= -3){
			vector<int>::iterator itv = find((*genotype)[parityIndex_].begin(),(*genotype)[parityIndex_].end(),vect[i]);
			if(itv!=(*genotype)[parityIndex_].end()){
				occupsites.push_back(i);
			}
		}else if(vect[i]==-3){
			occupsitesneutral.push_back(i);
		}
	}
	return (vector<vector<int>>{occupsites,occupsitesneutral});	
}

void Model::sitemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype){
	// for each position with allele, prob v to mutate and if mutation, choose randomly 1 chrom to mutate
	vector<int> occupied = occupiedsites(Alleleforeachpos_, genotype)[0];
	vector<int> occupiedneutral = occupiedsites(Alleleforeachpos_, genotype)[1];
	for (auto i : occupied){ //for all the sites recognized by an allele
		if (bernoulli_draw(2*N_*v_)){ //if mutation at this site
			int mutchrom = choose(2*N_); //only one indiv will preform a mutation at this site
			if ((*population)[parityIndex_][mutchrom][i]==1){ //if the site at this position is 
				(*population)[parityIndex_][mutchrom][i]=0; //perform the mutation
			}
		}
	}
	for (auto i : occupiedneutral){ //for all the neutral sites
		if (bernoulli_draw(2*N_*w_)){ //if mutation at this site
			int mutchrom = choose(2*N_); //only one indiv will preform a mutation at this site
			(*population)[parityIndex_][mutchrom][i]=((*population)[parityIndex_][mutchrom][i]+1)%2; //perform the mutation active->inacive or inactive->active
		}
	}
}

void Model::allelemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,double>* Ageallele, map<int,vector<double>>* infoperallele){
	// for each chromosome, mutation of allele with proba u and if at least one mutation update populations (change prdm9 allele at its position), genotypes(same), map (site pos associated to the new allele) and allele for each pos (change the allele corresponding to each pos : -1 -> new allele)
	vector<int> mutchromallele;
	for (int i=0; i<2*N_; i++){
		if(bernoulli_draw(u_)){
			nballele_+=1;
			vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
			vector<int> newpos = choosemany(nbsite_, freepos); 
			int newallele = nballele_-1;
			for(auto j : newpos){
				Alleleforeachpos_[j]=newallele;//update allele for each pos
			}
			Siteforeacheallele_[newallele]=newpos; //update map
			(*genotype)[parityIndex_][i]=newallele;//update genotype
			(*population)[parityIndex_][i][indPrdm9_]=newallele;//update pop
			(*Ageallele)[newallele]=freqall(newallele, genotype, population);//update age of the new allele
			if(ismigration_){
				infoperallele1_[newallele]={0,0,0,0,0,0,get_mean_affinity(newallele,population)};
				infoperallele2_[newallele]={0,0,0,0,0,0,get_mean_affinity(newallele,population)};
			}else{
				infoperallele_[newallele]={0,0,0,0,0,0,get_mean_affinity(newallele,population)}; ///////////////////a ajouter au fichiers
			}
		}
	}
}

void Model::updatemissingallele(){
//update the map if one allele has gone extinct : allele suppressed in map, all positions of this allele are set to -1 (free position) and all the activations are set to one
	vector<int> alleletoerase;
	for(auto &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			if(not ismigration_){
				vector<int>::iterator itvect = find(genotypes_[parityIndex_].begin(),genotypes_[parityIndex_].end(),it.first);
				if(itvect == genotypes_[parityIndex_].end()){
					for(auto &i : it.second){
						Alleleforeachpos_[i]=-1;
						for(int j=0; j<2*N_; j++){
							populations_[parityIndex_][j][i]=1;
						}
					}
					alleletoerase.push_back(it.first);
				}
			}else{
				vector<int>::iterator itvect1 = find(genotypes1_[parityIndex_].begin(),genotypes1_[parityIndex_].end(),it.first);
				vector<int>::iterator itvect2 = find(genotypes2_[parityIndex_].begin(),genotypes2_[parityIndex_].end(),it.first);
				if(itvect1 == genotypes1_[parityIndex_].end() and itvect2 == genotypes2_[parityIndex_].end()){
					for(auto &i : it.second){
						Alleleforeachpos_[i]=-1;
						for(int j=0; j<2*N_; j++){
							populations1_[parityIndex_][j][i]=1;
							populations2_[parityIndex_][j][i]=1;
						}
					}
					alleletoerase.push_back(it.first);
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
void Model::printpop(int n, vector<vector<vector<int>>>population){
	if(n==2){
		for (auto i : population){
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
		vector<vector<int>> pop = population[n];
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
int Model::Meiosis(int no_chrom_ind, int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, vector<vector<int>>* nbfailedmeiosis, double* q){ //Perform the meiosis of one individual
	//Initialisation time meiosis
	double tps_binding, tps_DSB, tps_vCO, tps_symq, tps_CONCO, tps_tot;
	clock_t t_1, t_2, t_3, t_4, t_5, t_6;
	/////////////////////////
	t_1=clock();
	/////////////////////////
	// choose one ind
	int indiv = 2*choose(N_);
	//homozygote or heterozygote
	int ind_gen=2;
	vector<int> zygote{(*genotype)[parityIndex_][indiv]};
	if((*genotype)[parityIndex_][indiv]!=(*genotype)[parityIndex_][indiv+1]){
		if(zygosity_){
			ind_gen=1;
		}
		zygote.push_back((*genotype)[parityIndex_][indiv+1]);
	}
	//PRDM9 binding
	vector<vector<int>>summarysites;
	vector<int>Z;
	vector<int> vectsites;
	int nblinksite = 0;
	for(auto z : zygote){
		for(auto i : Siteforeacheallele_[z]){
			vector<int>linkedsites=vector<int>(5,0);
			linkedsites[0]=i;
			double p_occup=ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]);
			bool islinked=false;
			for(int j=0; j<4; j++){
				int chrom=indiv+j/2;
				if((*population)[parityIndex_][chrom][i]==1){
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
	////////////////////////////
	t_2=clock();
	tps_binding=(float)(t_2-t_1)/CLOCKS_PER_SEC;
	//printf("tps_binding = %f\n", tps_binding);
	////////////////////////////
	//DSB
	double pDSB=double(nbDSB_)/nblinksite;
	vector<vector<int>>vectsitedsb;
	vector<double> alleleDSB;
	int checknblink =0;
	vector<vector<int>> vect_CO;
	vector<double>alleleCO;
	for(int i=0; i<summarysites.size(); i++){
		//proba DSB
		vector<int> link = vectfreesites(vector<int>(summarysites[i].begin()+1, summarysites[i].end()), 1);
		int nbdsbpersite=0;
		bool dsb=false;
		for(auto j : link){
			checknblink++;
			if(bernoulli_draw(pDSB)){
				dsb=true;
				summarysites[i][j+1]=2;// DSB -> 2
				nbdsbpersite+=1;
				vectsitedsb.push_back({i,j});
				alleleDSB.push_back(Z[i]);
				try{
					if (nbdsbpersite>1 and withDSB_){
						(*nbfailedmeiosis)[nb_gen][0]+=1; 
						(*infoperallele)[Z[i]][0]+=1;
						(*infoperallele)[Z[i]][1]+=1;
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
		}
		////////////////////////////
		//t_3=clock();
		//tps_DSB=(float)(t_3-t_2)/CLOCKS_PER_SEC;
		//printf("tps_DSB = %f\n", tps_DSB);
		////////////////////////////
		//Crossing over
		vector<int> vco;
		if (dsb){
			for(int indexnbdsb = 0; indexnbdsb<nbdsbpersite; indexnbdsb++){
				vco = {summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]};
				if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==0 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==1){
					if((summarysites[i][3]==1 or summarysites[i][3]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=2){/// in these cases, I suppose that 2 DSB can perform a CO
						vco.push_back(2);
					}
					if((summarysites[i][4]==1 or summarysites[i][4]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=3){
						vco.push_back(3);
					}
				}else if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==2 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==3 or summarysites[i][1]==2){
					if((summarysites[i][1]==1 or summarysites[i][1]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=0 ){
						vco.push_back(0);
					}
					if((summarysites[i][2]==1 or summarysites[i][2]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=1){
						vco.push_back(1);
					}
				}
				if(vco.size()>2){
					vect_CO.push_back(vco);
					alleleCO.push_back(Z[i]);
				}
			}
		}
	}
	////////////////////////////
	t_4=clock();
	tps_vCO=(float)(t_4-t_2)/CLOCKS_PER_SEC;
	//printf("tps_vCO = %f\n", tps_vCO);
	////////////////////////////
	try{
		if(vectsitedsb.size()==0){
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
			throw int(1);
		}else if(not vect_CO.size()){
			(*nbfailedmeiosis)[nb_gen][2]+=1;
			if(zygote.size()==1){
				(*infoperallele)[zygote[0]][0]+=2;
				(*infoperallele)[zygote[0]][3]+=2;
				(*infoperallele)[zygote[0]][5]+=2;
			}else if(zygote.size()==2){
				(*infoperallele)[zygote[0]][0]+=1;
				(*infoperallele)[zygote[0]][3]+=1;
				(*infoperallele)[zygote[1]][0]+=1;
				(*infoperallele)[zygote[1]][3]+=1;
				(*infoperallele)[zygote[0]][5]+=1;
				(*infoperallele)[zygote[1]][5]+=1;
			}
			throw int(2);
		}
	}
	catch(int const& error_nb){
		if(error_nb==1){
			//cerr << "No DSB" << endl;
		}else if(error_nb==2){
			//cerr << "No symmetrical sites (binding + DSB)" << endl;
		}
		return-1;
	}
	//symmetrical binding
	*q=*q+double(vect_CO.size())/(vectsitedsb.size());
	if(zygote.size()==1){
		(*infoperallele)[zygote[0]][4]+=2*(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0])); 
		(*infoperallele)[zygote[0]][5]+=2;
	}else if(zygote.size()==2){
		if (count(alleleDSB.begin(),alleleDSB.end(),zygote[0]) != 0){
			(*infoperallele)[zygote[0]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0]));
			(*infoperallele)[zygote[0]][5]+=1;
			if(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))==0){////////////////////added
				(*infoperallele)[zygote[0]][3]+=1;
			}	
		}else{
			(*infoperallele)[zygote[0]][2]+=1;
		}
		if (count(alleleDSB.begin(),alleleDSB.end(),zygote[1]) != 0){
			(*infoperallele)[zygote[1]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[1]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[1]));
			(*infoperallele)[zygote[1]][5]+=1;
			if(double(count(alleleCO.begin(),alleleCO.end(),zygote[1]))==0){/////////////////////////added
				(*infoperallele)[zygote[1]][3]+=1;
			}
		}else{
			(*infoperallele)[zygote[1]][2]+=1;
		}
	}
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
	////////////////////////////
	t_5=clock();
	tps_symq=(float)(t_5-t_4)/CLOCKS_PER_SEC;
	//printf("tps_symq = %f\n", tps_symq);
	////////////////////////////
	//Then if meiosis happen
	//choose one gamete
	int no_chromatide=choose(4);
	int no_current_chrom = indiv;
	int no_homologue_chrom = indiv+1;
	if(no_chromatide==2 or no_chromatide==3){
		no_current_chrom=indiv+1;
		no_homologue_chrom=indiv;
	}
	//performe recombination (CO and NCO)
	if(no_chromatide==index_CO[1] or no_chromatide==index_CO[2]){
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
		copy( (*population)[parityIndex_][no_current_chrom].begin(), (*population)[parityIndex_][no_current_chrom].end(), (*population)[(parityIndex_+1)%2][no_chrom_ind].begin() );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}
	//return 0 if fail somewhere and 1 if meiosis successfull
	////////////////////////////
	t_6=clock();
	tps_CONCO=(float)(t_6-t_5)/CLOCKS_PER_SEC;
	//printf("tps_CONCO = %f\n", tps_CONCO);
	tps_tot=(float)(t_6-t_1)/CLOCKS_PER_SEC;
	//printf("tps_tot = %f\n", tps_tot);
	////////////////////////////
	return 1;
}


//methode fill new population (make as many meiosis as it is necessary to fill the entire new population)
void Model::fillnewpop(int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, vector<vector<int>>* nbfailedmeiosis, double* q){
	////////////////////////////////////
	double tps_fill, tps_meiosis;
	clock_t tps1, tps2, tps3, tps4, tps5,tps6;
	tps1=clock();
	double temps_moy=0;
	int nbmei=0;
	*q=0;
	for(int indnewpop=0; indnewpop<2*N_; indnewpop++){
		tps3=clock();
		int meiosisState = Meiosis(indnewpop, nb_gen, population, genotype, infoperallele, nbfailedmeiosis, q);
		////////////////////////////////////
		tps4=clock();
		tps_meiosis=(float)(tps4-tps3)/CLOCKS_PER_SEC;
		//printf("tps_meiosis = %f\n", tps_meiosis);
		temps_moy+=tps_meiosis;
		nbmei+=1;
		////////////////////////////////////
		while (meiosisState==-1){
			tps5=clock();
			meiosisState = Meiosis(indnewpop, nb_gen, population, genotype, infoperallele, nbfailedmeiosis, q);
			(*nbfailedmeiosis)[nb_gen][3]+=1;
			////////////////////////////////////
			tps6=clock();
			tps_meiosis=(float)(tps6-tps5)/CLOCKS_PER_SEC;
			//printf("tps_meiosis = %f\n", tps_meiosis);
			temps_moy+=tps_meiosis;
			nbmei+=1;
			////////////////////////////////////
		}
	}
	////////////////////////////////////
	temps_moy=temps_moy/nbmei;
	cout<< "temps_moy = "<<temps_moy<<endl;
	////////////////////////////////////
	for(int indpop=0; indpop<2*N_;indpop++){
		(*genotype)[(parityIndex_+1)%2][indpop]=(*population)[(parityIndex_+1)%2][indpop][indPrdm9_];
	}
	*q=*q/(2*N_+(*nbfailedmeiosis)[nb_gen][2]);
	////////////////////////////////////
	tps2=clock();
	tps_fill=(float)(tps2-tps1)/CLOCKS_PER_SEC;
	//printf("tps_fill = %f\n", tps_fill);
	///////////////////////////////////
}

//methode which mix together and runs all the previous functions during X generations 
void Model::manygenerations(){
	//Initialisation of files
	//Trace file : contains all the descriptive statistics averaged in the population for each generation
	ofstream generalfile ((name_+".trace").c_str());
	generalfile << "Generation_number" << '\t' << "Total_number_of_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' << "Mean_Age" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" << '\t' << "q_intra" << '\t' << "fertility_intra" << '\t' << "q_analytic" <<  '\t' << "fertility_analytic" << '\n';
    generalfile.flush();
	//Allele file : contains all the mean descriptive statistics for each allele segregating in the population for each generation
    ofstream allelefile ((name_+".allele").c_str());
	allelefile << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" << '\t' << "Fertility_allele" << '\t' << "mean_affinity" << '\t' << "q_analytic" <<  '\t' << "fertility_analytic" << '\n';
    allelefile.flush();
	//Params file : contains all the parameters values for this simulation
	ofstream paramsfile ((name_+".params").c_str());
	paramsfile << "N" << '\t' << N_ << '\n' << "L" << '\t' << L_ << '\n' << "nbsite" << '\t' << nbsite_ << '\n' << "indPrdm9" << '\t' << indPrdm9_ << '\n' << "nballele" << '\t' << nballele_ << '\n' << "parityIndex" << '\t' << parityIndex_ << '\n' << "u" << '\t' << u_ << '\n' << "v" << '\t' << v_ << '\n' << "w" << '\t' << w_ << '\n' << "meanaff" << '\t' << meanaff_ << '\n' << "varaff" << '\t' << varaff_ << '\n' << "nbDSB" << '\t' << nbDSB_ << '\n' << "nbGenerations" << '\t' << nbGenerations_ << '\n' << "ismigration" << '\t' << ismigration_ << '\n' << "zygosity" << '\t' << zygosity_ << '\n' << "withDSB" << '\t' << withDSB_ << '\n' << "everygen" << '\t' << everygen_ << '\n' << "m" << '\t' << m_ << '\n' << "nbgenmig" << '\t' << nbgenmig_ << '\n';
    paramsfile.flush();
	ofstream generalfile1;
	ofstream allelefile1;
	ofstream generalfile2;
	ofstream allelefile2;
	vector<vector<vector<vector<int>>>*> vectpop = vector<vector<vector<vector<int>>>*> {&populations_};
	vector<vector<vector<int>>*> vectgen = vector<vector<vector<int>>*> {&genotypes_};
	vector<map<int,double>*> vectage = vector<map<int,double>*> {&Ageallele_};
	vector<map<int,vector<double>>*> vectinfo = vector<map<int,vector<double>>*> {&infoperallele_};
	double q1_=0;
	double q2_=0;
	vector<double*> vectq = vector<double*> {&q_};
	vector<vector<vector<int>>*> vectfailed = vector<vector<vector<int>>*> {&nbfailedmeiosis_};
	int lenvect = vectpop.size();
	for(int indgeneration=0; indgeneration<nbGenerations_; indgeneration++){
		//Initialisation time
		double temps_1, temps_printfile, temps_fillnewpop, temps_afterfillpop, temps_calculstat, temps_calcunit, temps_calcloop, temps_mutations_1;
		clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
		t1=clock();
		cout<<"generation "<<indgeneration<<endl;
		t9=clock();
		//set to 0 all the vectors containing the informatons from an allele during the meiosis
		for(int i_vect=0; i_vect<lenvect; i_vect++){
			for (auto const &it : *vectinfo[i_vect]){
				(*vectinfo[i_vect])[it.first]=vector<double>{0,0,0,0,0,0,get_mean_affinity(it.first,vectpop[i_vect])};
			}
			//printpop(2, (*vectpop[i_vect]));
			//printgen(2, (*vectgen[i_vect]));
			//printposallele();
			//printallelepos();
			sitemutation(vectpop[i_vect], vectgen[i_vect]);
			allelemutation(vectpop[i_vect], vectgen[i_vect], vectage[i_vect], vectinfo[i_vect]);
			//printpop(2, (*vectpop[i_vect]));
			//printgen(2, (*vectgen[i_vect]));
			//printposallele();
			//printallelepos();
		}
		updatemissingallele();
		t10=clock();
		temps_mutations_1=(float)(t10-t9)/CLOCKS_PER_SEC;
		printf("temps_mutations_1 = %f\n", temps_mutations_1);
		// If we want to split the population in two and perform migration (or not)
		if(indgeneration==nbgenmig_){
			cout<<"MIGRATION _________________________________________________________________________________"<<endl;
			//initialisation of containers and local variables
			ismigration_=true;
			populations1_=populations_;
			populations2_=populations_;
			genotypes1_=genotypes_;
			genotypes2_=genotypes_;
			Ageallele1_=Ageallele_;
			Ageallele2_=Ageallele_;
			infoperallele1_=infoperallele_;
			infoperallele2_=infoperallele_;
			nbfailedmeiosis1_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
			nbfailedmeiosis2_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
			vectpop = vector<vector<vector<vector<int>>>*> {&populations1_,&populations2_};
			vectgen = vector<vector<vector<int>>*> {&genotypes1_,&genotypes2_};
			vectage = vector<map<int,double>*> {&Ageallele1_,&Ageallele2_};
			vectinfo = vector<map<int,vector<double>>*> {&infoperallele1_,&infoperallele2_};
			vectq = vector<double*> {&q1_,&q2_};
			vectfailed = vector<vector<vector<int>>*> {&nbfailedmeiosis1_,&nbfailedmeiosis2_};
			lenvect = vectpop.size();
			//Initialise file generalfile et allelefile
			/////////////////////////////////////////////////add analytical q and fertility
			generalfile1.open ((name_+"_1.trace").c_str());
			generalfile1 << "Generation_number" << '\t' << "Number_of_allele" << '\t' << "Total_nb_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' << "Mean_Age" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" << '\t' << "q_intra" << '\t' << "q_inter" << '\t' << "fertility_intra" << '\t' << "fertility_inter" << '\t' << "FST_neutral" << '\t' << "FST_PRDM9" << '\n';
    			generalfile1.flush();
    			allelefile1.open ((name_+"_1.allele").c_str());
			allelefile1 << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" << '\t' << "Fertility_allele" << '\t' << "mean_affinity" << '\n';
    			allelefile1.flush();
			generalfile2.open ((name_+"_2.trace").c_str());
			generalfile2 << "Generation_number" << '\t' << "Number_of_allele" << '\t' << "Total_nb_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' << "Mean_Age" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" << '\t' << "q_intra" << '\t' << "q_inter" << '\t' << "fertility_intra" << '\t' << "fertility_inter" << '\t' << "FST_neutral" << '\t' << "FST_PRDM9" << '\n';
    			generalfile2.flush();
    			allelefile2.open ((name_+"_2.allele").c_str());
			allelefile2 << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" << '\t' << "Fertility_allele" << '\t' << "mean_affinity" << '\n';
    			allelefile2.flush();
		}
		if(ismigration_){
			migration();
		}
		//////////////////////////
		t4=clock();
		//////////////////////////
		for(int i_vect=0; i_vect<lenvect; i_vect++){
			fillnewpop(indgeneration, vectpop[i_vect], vectgen[i_vect], vectinfo[i_vect], vectfailed[i_vect], vectq[i_vect]);
		}
		////////////////////////
		t5=clock();
		temps_fillnewpop=(float)(t5-t4)/CLOCKS_PER_SEC;
		printf("temps_fillnewpop = %f\n", temps_fillnewpop);
		////////////////////////
		parityIndex_=(parityIndex_+1)%2;
		updatemissingallele();
		for(int i_vect=0; i_vect<lenvect; i_vect++){
			for (auto &it : *vectage[i_vect]){
				it.second=it.second+freqall(it.first, vectgen[i_vect], vectpop[i_vect]);
			}
			for (auto &it : *vectinfo[i_vect]){
				if(it.second[5]!=0){
					it.second[4]=it.second[4]/it.second[5];
				}
			}
		}
		vector<int> tot_allele_nb = vector<int>(2,0);
		vector<int> allele_nb_trace=vector<int>(2,0);
		vector<double> current_div=vector<double>(2,0);
		vector<double> current_act=vector<double>(2,0);
		vector<double> Fertility_rate=vector<double>(2,0);
		vector<double> twoDSB=vector<double>(2,0);
		vector<double> noDSB=vector<double>(2,0);
		vector<double> nosym=vector<double>(2,0);
		vector<double> current_q=vector<double>(2,0);
		vector<double> mean_age=vector<double>(2,0);
		vector<double> q_fertility;
		double FST_neutral;
		double FST_PRDM9;
		vector<int> allele_nb;
		vector<double> q_fertility_hybrid;
		vector<double> q_fertility_1;
		vector<double> q_fertility_2;
		vector<double> freqallele;
		vector<double> freqallele1;
		vector<double> freqallele2;
		vector<double> activityall;
		vector<double> activityall1;
		vector<double> activityall2;
		vector<double> ageallele;
		vector<double> ageallele1;
		vector<double> ageallele2;
		vector<double> qallele;
		vector<double> qallele1;
		vector<double> qallele2;
		vector<double> fertilityallele;
		vector<double> fertilityallele1;
		vector<double> fertilityallele2;
		vector<double> meanaffinity;
		vector<double> meanaffinity1;
		vector<double> meanaffinity2;
		map<int,vector<double>> analytic_allele;///////////////////////////////////////////////////////////
		map<int,vector<double>> analytic_indiv;///////////////////////////////////////////////////////////////
		////////////////////////
		t6=clock();
		temps_afterfillpop=(float)(t6-t5)/CLOCKS_PER_SEC;
		printf("temps_afterfillpop = %f\n", temps_afterfillpop);
		////////////////////////
		if (indgeneration % everygen_ ==0)  {
			if(ismigration_==false){
				// generalfile
				allele_nb_trace[0]=get_allele_number(vectgen,true)[0];
				current_div[0] = get_current_diversity(&genotypes_);
				current_act[0] = get_current_activity(&genotypes_, &populations_);
				Fertility_rate[0] = 1-(double(nbfailedmeiosis_[indgeneration][3])/(2*N_+nbfailedmeiosis_[indgeneration][3]));
				twoDSB[0] = double(nbfailedmeiosis_[indgeneration][0])/(2*N_+nbfailedmeiosis_[indgeneration][3]);
				noDSB[0] = double(nbfailedmeiosis_[indgeneration][1])/(2*N_+nbfailedmeiosis_[indgeneration][3]);
				nosym[0] = double(nbfailedmeiosis_[indgeneration][2])/(2*N_+nbfailedmeiosis_[indgeneration][3]);
				current_q[0] = q_;
				mean_age[0] = get_mean_age(vectgen[0], vectage[0]);
				////////////////////////
				t7=clock();
				temps_calcunit=(float)(t7-t6)/CLOCKS_PER_SEC;
				("temps_calcunit = %f\n", temps_calcunit);
				////////////////////////
				for (auto const &it : Siteforeacheallele_){
	    				if(it.first!=-2){
						allele_nb.push_back(it.first);
						freqallele.push_back(freqall(it.first, &genotypes_, &populations_));
						activityall.push_back(actall(it.first, &populations_));
						ageallele.push_back(get_age_allele(it.first, &Ageallele_));
						qallele.push_back(get_info_allele(it.first, &infoperallele_)[0]);
						fertilityallele.push_back(get_info_allele(it.first, &infoperallele_)[1]);
						if(it.first==-3){
							meanaffinity.push_back(get_mean_affinity(it.first,&populations_));
						}else{
							meanaffinity.push_back(infoperallele_[it.first][6]);
						}
            				}
            			}
				q_fertility = get_q_fertility_indep(vectpop, vectgen, nbloop_, 1);
				/////////////////////////////////////////////////////////////////////////////////////////////////////
				/*analytic_allele=q_fert_allele_analytique(&genotypes_);
				for (auto const &it : analytic_allele){
					cout<<it.first<<"=>";
					for(auto const &i : it.second){
						cout<<" "<<i;
					}
					cout<<endl;
				}
				cout<<'\n';*/
				analytic_indiv=q_fert_individual_analytique(&genotypes_,&populations_);
				/*for (auto const &it : analytic_indiv){
					cout<<it.first<<"=>";
					for(auto const &i : it.second){
						cout<<" "<<i;
					}
					cout<<endl;
				}
				cout<<'\n';*/
				/////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////
				t8=clock();
				temps_calcloop=(float)(t8-t7)/CLOCKS_PER_SEC;
				printf("temps_calcloop = %f\n", temps_calcloop);
				////////////////////////
			}else if(ismigration_==true){
				for(int indfile=0; indfile<2; indfile++){
					allele_nb_trace[indfile] = get_allele_number(vectgen,false)[indfile];
					current_div[indfile] = get_current_diversity(vectgen[indfile]);
					current_act[indfile] = get_current_activity(vectgen[indfile], vectpop[indfile]);
					Fertility_rate[indfile] = 1-(double((*vectfailed[indfile])[indgeneration][3])/(2*N_+(*vectfailed[indfile])[indgeneration][3]));
					twoDSB[indfile] = double((*vectfailed[indfile])[indgeneration][0])/(2*N_+(*vectfailed[indfile])[indgeneration][3]);
					noDSB[indfile] = double((*vectfailed[indfile])[indgeneration][1])/(2*N_+(*vectfailed[indfile])[indgeneration][3]);
					nosym[indfile] = double((*vectfailed[indfile])[indgeneration][2])/(2*N_+(*vectfailed[indfile])[indgeneration][3]);
					current_q[indfile] = *vectq[indfile];
					mean_age[indfile] = get_mean_age(vectgen[indfile], vectage[indfile]);
				}
				tot_allele_nb = get_allele_number(vectgen,true);
				FST_neutral = get_FST_neutral(vectpop);
				FST_PRDM9 = get_FST_PRDM9(vectgen);
				q_fertility_hybrid = get_q_fertility_indep(vectpop, vectgen, nbloop_, 0);
				q_fertility_1 = get_q_fertility_indep(vectpop, vectgen, nbloop_, 1);
				q_fertility_2 = get_q_fertility_indep(vectpop, vectgen, nbloop_, 2);
				for (auto const &it : Siteforeacheallele_){
	    				if(it.first!=-2){
							allele_nb.push_back(it.first);
							freqallele1.push_back(freqall(it.first, &genotypes1_, &populations1_));
							freqallele2.push_back(freqall(it.first, &genotypes2_, &populations2_));
							activityall1.push_back(actall(it.first, &populations1_));
							activityall2.push_back(actall(it.first, &populations2_));
							ageallele1.push_back(get_age_allele(it.first, &Ageallele1_));
							ageallele2.push_back(get_age_allele(it.first, &Ageallele2_));
							qallele1.push_back(get_info_allele(it.first, &infoperallele1_)[0]);
							qallele2.push_back(get_info_allele(it.first, &infoperallele2_)[0]);
							fertilityallele1.push_back(get_info_allele(it.first, &infoperallele1_)[1]);
							fertilityallele2.push_back(get_info_allele(it.first, &infoperallele2_)[1]);
							if(it.first==-3){
								meanaffinity1.push_back(get_mean_affinity(it.first,&populations1_));
								meanaffinity2.push_back(get_mean_affinity(it.first,&populations2_));
							}else{
								meanaffinity1.push_back(infoperallele1_[it.first][6]);
								meanaffinity2.push_back(infoperallele2_[it.first][6]);
							}
            				}
            			}
			}
			///////////////////////////
			t2=clock();
			temps_calculstat=(float)(t2-t6)/CLOCKS_PER_SEC;
			printf("temps_calculstat = %f\n", temps_calculstat);
			temps_1=(float)(t2-t1)/CLOCKS_PER_SEC;
			printf("temps_1 = %f\n", temps_1);
			//////////////////////////
			if(ismigration_==false){
				generalfile << indgeneration << '\t' << allele_nb_trace[0] << '\t' << current_div[0] << '\t'  << current_act[0] << '\t' << mean_age[0] << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << Fertility_rate[0] << '\t' << twoDSB[0] << '\t'<< noDSB[0] << '\t' << nosym[0] << '\t' << current_q[0] << '\t' << q_fertility[0] << '\t' << q_fertility[1] << '\t' << analytic_indiv[-1][0] <<  '\t' << analytic_indiv[-1][1] << '\n';
            			generalfile.flush();
            			for (int indall = 0; indall<allele_nb.size(); indall++){
            				allelefile << indgeneration << '\t' << allele_nb[indall] << '\t' << freqallele[indall] << '\t'  << activityall[indall] << '\t' << ageallele[indall] << '\t' << qallele[indall] << '\t' << fertilityallele[indall] << '\t' << meanaffinity[indall] << '\t' <<analytic_indiv[int(allele_nb[indall])][0] <<  '\t' << analytic_indiv[int(allele_nb[indall])][1] << '\n';
            				allelefile.flush();
            			}
			}else if(ismigration_==true){
				generalfile1 << indgeneration << '\t' << allele_nb_trace[0] << '\t' << tot_allele_nb[0] << '\t' << current_div[0] << '\t'  << current_act[0] << '\t' << mean_age[0] << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << Fertility_rate[0] << '\t' << twoDSB[0] << '\t'<< noDSB[0] << '\t' << nosym[0] << '\t' << current_q[0] << '\t' << q_fertility_1[0] << '\t' << q_fertility_hybrid[0] << '\t' << q_fertility_1[1] << '\t' << q_fertility_hybrid[1] << '\t' << FST_neutral << '\t' << FST_PRDM9 << '\n';
            			generalfile1.flush();
            			for (int indall = 0; indall<allele_nb.size(); indall++){
            				allelefile1 << indgeneration << '\t' << allele_nb[indall] << '\t' << freqallele1[indall] << '\t'  << activityall1[indall] << '\t' << ageallele1[indall] << '\t' << qallele1[indall] << '\t' << fertilityallele1[indall] << '\t' << meanaffinity1[indall] << '\n';
            				allelefile1.flush();
            			}
            			generalfile2 << indgeneration << '\t' << allele_nb_trace[1] << '\t' << tot_allele_nb[0] << '\t' << current_div[1] << '\t'  << current_act[1] << '\t' << mean_age[1] << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << Fertility_rate[1] << '\t' << twoDSB[0] << '\t'<< noDSB[1] << '\t' << nosym[1] << '\t' << current_q[1] << '\t' << q_fertility_2[0] << '\t' << q_fertility_hybrid[0] << '\t' << q_fertility_2[1] << '\t' << q_fertility_hybrid[1] << '\t' << FST_neutral << '\t' << FST_PRDM9 << '\n';
            			generalfile2.flush();
            			for (int indall = 0; indall<allele_nb.size(); indall++){
            				allelefile2 << indgeneration << '\t' << allele_nb[indall] << '\t' << freqallele2[indall] << '\t'  << activityall2[indall] << '\t' << ageallele2[indall] << '\t' << qallele2[indall] << '\t' << fertilityallele2[indall] << '\t' << meanaffinity2[indall] << '\n';
            				allelefile2.flush();
            			}
			}
			///////////////////////////
			t3=clock();
			temps_printfile=(float)(t3-t2)/CLOCKS_PER_SEC;
			printf("temps_printfile = %f\n", temps_printfile);
			//////////////////////////
        	}
	}
	generalfile.close();
	allelefile.close();
	paramsfile.close();
	generalfile1.close();
	allelefile1.close();
	generalfile2.close();
	allelefile2.close();
}

//get functions for descriptive statistics

//get number of different alleles in the population
vector<int> Model::get_allele_number(vector<vector<vector<int>>*> vectgen, bool nbtot){
	if(vectgen.size()==1 or nbtot){
		return vector<int>{int(Siteforeacheallele_.size()-2)};
	}else if(vectgen.size()==2){
		vector<int> allnb=vector<int>(3,0);
		for(int i=0; i<2; i++){
			for(auto &it : Siteforeacheallele_){
				if(it.first!=-3 and it.first!=-2){
					vector<int>::iterator itvect = find((*vectgen[i])[parityIndex_].begin(),(*vectgen[i])[parityIndex_].end(),it.first);
					if(itvect != (*vectgen[i])[parityIndex_].end()){
						allnb[i]+=1;
					}
				}
			}
		}
		allnb[2]=Siteforeacheallele_.size()-2;
		return allnb;
	}
}

//get the age of each allele segregating in the population
double Model::get_age_allele(int allname, map<int,double>* Ageallele){
	if(allname==-3){
		return 0;
	}else{
		return (*Ageallele)[allname];
	}
}

//get the mean age of the allele segregating in the population
double Model::get_mean_age(vector<vector<int>>* genotype, map<int,double>* Ageallele){
	double meanage=0;
	for(auto const allele : (*genotype)[parityIndex_]){
		meanage+=(*Ageallele)[allele];
	}
	return meanage/((*genotype)[parityIndex_].size());
}

//get the containers for informations for each allele (number of successfull meiosis, number of failed meiosis due to a lack of DSB or symmetrical binding or even a lack of linked sie)
vector<double> Model::get_info_allele(int allname, map<int,vector<double>>* infoperallele){
	if(allname==-3){
		return {0,0};
	}else{
		typedef map<int,vector<double>>::iterator mi;
		if ( (*infoperallele).find(allname) != (*infoperallele).end() ) {
			if((*infoperallele)[allname][0]+(*infoperallele)[allname][5]-(*infoperallele)[allname][3]!=0){
				return vector<double>{(*infoperallele)[allname][4],double((*infoperallele)[allname][5]-(*infoperallele)[allname][3])/((*infoperallele)[allname][0]+(*infoperallele)[allname][5]-(*infoperallele)[allname][3])};
			}else{
				return vector<double>{(*infoperallele)[allname][4],0};
			}
		}
	}
}

//give the frequence of each allele in the population
double Model::freqallele(int allelename, vector<vector<int>>* genotype){
		return double (count((*genotype)[parityIndex_].begin(), (*genotype)[parityIndex_].end(), allelename))/((*genotype)[parityIndex_].size());
}

//get the diversity at one given generation
double Model::get_current_diversity(vector<vector<int>>* genotype){
	double sumfreq=0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			sumfreq+=(freqallele(it.first, genotype))*(freqallele(it.first, genotype));
		}
	}
	return 1/sumfreq;
}

//give the mean activity of a given allele
double Model::activitymoyallele(int allele,  vector<vector<vector<int>>>* population){
	double moyact=0;
	for(auto all : Siteforeacheallele_[allele]){
		double moyactsite=0;
		for(int i=0; i<2*N_; i++){
			if((*population)[parityIndex_][i][all]==1){
				moyactsite+=1;			
			}
		}
		moyact+=moyactsite/(2*N_);
	}
	moyact=moyact/nbsite_;
	return (moyact);
}

//???????????
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

//give the frequency of an allele or neural site
double Model::freqall(int allele, vector<vector<int>>* genotype, vector<vector<vector<int>>>* population){
	if (allele==-3){
		return freqneutral(population)[0];
	}else{
		return freqallele(allele, genotype);
	}
}

//give the activity of an allele or neutral site
double Model::actall(int allele,  vector<vector<vector<int>>>* population){
	if (allele==-3){
		return freqneutral(population)[1];
	}else{
		return activitymoyallele(allele, population);
	}
}

//give the average activity of allele the allele in the population
double Model::get_current_activity(vector<vector<int>>* genotype, vector<vector<vector<int>>>* population){
	double moytotact=0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			moytotact+=activitymoyallele(it.first, population)*freqallele(it.first, genotype);
		}
	}
	return moytotact;
}

//give the mean affinity of all the active sites of a given allele
double Model::get_mean_affinity(int allele, vector<vector<vector<int>>>* pop){
	double meanaffinity=0;
	int nbactivesite=0;
	for(auto const site : Siteforeacheallele_[allele]){
		int activitysite=0;
		for(int siteind=0; siteind<2*N_; siteind++){
			if((*pop)[parityIndex_][siteind][site]==1){
				activitysite+=1;
				nbactivesite+=1;
			}
		}
		meanaffinity+=Affinity_[site]*activitysite;
	}
	return meanaffinity/nbactivesite;
}

//give the FST of neutral site
double Model::get_FST_neutral(vector<vector<vector<vector<int>>>*> vectpop){
	vector<double> p1p2 = vector<double>(2,0);
	double p=0;
	vector<double> H1H2 = vector<double>(2,0);
	double Hinter=0;
	for(auto all : Siteforeacheallele_[-3]){
		for(int j=0; j<2; j++){
			for(int i=0; i<2*N_; i++){
				if((*vectpop[j])[parityIndex_][i][all]==1){
					p1p2[j]+=1;		
				}
			}
			p1p2[j]=p1p2[j]/(2*N_);
			H1H2[j]+=2*p1p2[j]*(1-p1p2[j]);
		}
		p=0.5*(p1p2[0]+p1p2[1]);
		Hinter+=2*p*(1-p);
	}
	H1H2[0]=H1H2[0]/nbsite_;
	H1H2[1]=H1H2[1]/nbsite_;
	Hinter=Hinter/nbsite_;
	double Hintra=0.5*(H1H2[0]+H1H2[1]);
	double FST=1-(Hintra)/(Hinter);
	return FST;
}

//give the PRDM9 FST
double Model::get_FST_PRDM9(vector<vector<vector<int>>*> vectgen){
	vector<double> H1H2=vector<double>(2,0);
	for(int i=0; i<2; i++){
		for (auto const &it : Siteforeacheallele_){
			if(it.first!=-3 and it.first!=-2){
				H1H2[i]+=((freqallele(it.first, vectgen[i]))*(freqallele(it.first, vectgen[i])));
			}
		}
		H1H2[i]=1-H1H2[i];
	}
	double Hintra=0.5*(H1H2[0]+H1H2[1]);
	vector<double> p1p2=vector<double>(2,0);
	double p =0;
	double Hinter =0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			for(int i=0; i<2; i++){
				p1p2[i]=freqallele(it.first, vectgen[i]);
			}
			p=0.5*(p1p2[0]+p1p2[1]);
			Hinter+=(p*p);
		}
	}
	Hinter=1-Hinter;
	if (Hinter==0){
		return 0;
	}else{
		double FST=1-Hintra/Hinter;
		return FST;
	}
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
		bool found=true;
		while(found==true){
			vector<int>::iterator it = find(newsites.begin(),newsites.end(),index);
			if(it != newsites.end()){
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

//function that performs the migration between two populations
void Model::migration(){
	int nb_mig = binomial_draw(N_, m_);
	vector<int> migrated_indiv_pop1 = choosemanymigration(int(nb_mig));
	vector<int> migrated_indiv_pop2 = choosemanymigration(int(nb_mig));
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

//give the mean probability of symetrical site and the mean of fertility in the population
vector<double> Model::get_q_fertility_indep(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen, int nbloop_, int nopop){
	int nb_good_meiosis = 0;
	int nb_no_sym = 0;
	double moy_q = 0;
	double moy_fertility;
	double res;
	for(int i=0; i<nbloop_; i++){
		if(vectpop.size()==1){
			res = get_q(vectpop[0], vectgen[0]);
		}else if(vectpop.size()==2){
			if(nopop==1){
				res = get_q(vectpop[0], vectgen[0]);
			}else if(nopop==2){
				res = get_q(vectpop[1], vectgen[1]);
			}else if(nopop==0){
				res = get_q_hybrid(vectpop,vectgen);
			}
		}
		if(res != -1 and res!= -2){
			nb_good_meiosis++;
			moy_q+=res;
		}
		else if(res == -2){
			nb_no_sym++;
		}
	}
	moy_q=moy_q/double(nb_good_meiosis+nb_no_sym);
	moy_fertility = double(nb_good_meiosis)/nbloop_;
	return vector<double> {moy_q,moy_fertility};
}






// calculation of q and fertility indepedently from the system (sampling)
//give q
double Model::get_q(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype){
	int indiv = 2*choose(N_);
	vector<int> genotype_indiv = vector<int>{(*genotype)[parityIndex_][indiv],(*genotype)[parityIndex_][indiv+1]};
	vector<vector<int>> indiv_chrom = vector<vector<int>>{(*population)[parityIndex_][indiv],(*population)[parityIndex_][indiv+1]};
	double q_indiv = q_two_hap(genotype_indiv, indiv_chrom);
	return q_indiv;
}

//give q of hybride
double Model::get_q_hybrid(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen){
	vector<int> genotype_indiv;
	vector<vector<int>> indiv_chrom;
	for(int i=0; i<2; i++){
		int indiv = 2*choose(N_);
		vector<int> genotype_parent = vector<int>{(*vectgen[i])[parityIndex_][indiv],(*vectgen[i])[parityIndex_][indiv+1]};
		vector<vector<int>> parent_chrom = vector<vector<int>>{(*vectpop[i])[parityIndex_][indiv],(*vectpop[i])[parityIndex_][indiv+1]};
		vector<int> gamete = get_one_gamete(genotype_parent,parent_chrom);
		while(gamete == vector<int>{-1}){
			indiv = 2*choose(N_);
			genotype_parent = vector<int> {(*vectgen[i])[parityIndex_][indiv],(*vectgen[i])[parityIndex_][indiv+1]};
			parent_chrom = vector<vector<int>> {(*vectpop[i])[parityIndex_][indiv],(*vectpop[i])[parityIndex_][indiv+1]};
			gamete = get_one_gamete(genotype_parent,parent_chrom);
		}
		indiv_chrom.push_back(gamete);
		genotype_indiv.push_back(gamete[indPrdm9_]);
	}
	double q_indiv = q_two_hap(genotype_indiv, indiv_chrom);
	return q_indiv;
}

//give 
double Model::q_two_hap(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom){ //return a vector of vector
	double q;
	int ind_gen=2;
	vector<int> zygote{genotype_indiv[0]};
	if(genotype_indiv[0]!=genotype_indiv[1]){ // if we give only one indiv don't need gen
		if(zygosity_){
			ind_gen=1;
		}
		zygote.push_back(genotype_indiv[1]);
	}
	vector<vector<int>>summarysites;
	vector<int>Z;
	vector<int> vectsites;
	int nblinksite = 0;
	for(auto z : zygote){
		for(auto i : Siteforeacheallele_[z]){
			vector<int>linkedsites=vector<int>(5,0);
			linkedsites[0]=i;
			double p_occup=ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]);
			bool islinked=false;
			for(int j=0; j<4; j++){
				int chrom=j/2;
				if(indiv_chrom[chrom][i]==1){ // if we give only one indiv don't need pop
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
	double pDSB=double(nbDSB_)/nblinksite;
	vector<vector<int>>vectsitedsb;
	vector<vector<int>> vect_CO;
	//vector<double>alleleDSB;
	//vector<double>alleleCO;
	for(int i=0; i<summarysites.size(); i++){
		vector<int> link = vectfreesites(vector<int>(summarysites[i].begin()+1, summarysites[i].end()), 1);
		int nbdsbpersite=0;
		bool dsb=false;
		for(auto j : link){
			if(bernoulli_draw(pDSB)){
				dsb=true;
				summarysites[i][j+1]=2;
				nbdsbpersite+=1;
				vectsitedsb.push_back({i,j});
				//alleleDSB.push_back(Z[i]);
				try{
					if (nbdsbpersite>1 and withDSB_){
						throw int(0);
					}
				}
				catch(int const& error_nb){
					if(error_nb==0){
						//cerr << "2 DSB on one site" << endl;
					}
					return -1;
				}
			}
		}
		vector<int> vco;
		if (dsb){
			for(int indexnbdsb = 0; indexnbdsb<nbdsbpersite; indexnbdsb++){
				vco = {summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]};
				if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==0 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==1){
					if((summarysites[i][3]==1 or summarysites[i][3]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=2){
						vco.push_back(2);
					}
					if((summarysites[i][4]==1 or summarysites[i][4]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=3){
						vco.push_back(3);
					}
				}else if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==2 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==3 or summarysites[i][1]==2){
					if((summarysites[i][1]==1 or summarysites[i][1]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=0 ){
						vco.push_back(0);
					}
					if((summarysites[i][2]==1 or summarysites[i][2]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=1){
						vco.push_back(1);
					}
				}
				if(vco.size()>2){
					vect_CO.push_back(vco);
					//alleleCO.push_back(Z[i]);
				}
			}
		}
	}
	try{
		if(vectsitedsb.size()==0){
			throw int(1);
		}else if(not vect_CO.size()){
			throw int(2);
		}
	}
	catch(int const& error_nb){
		if(error_nb==1){
			//cerr << "No DSB" << endl;
			return -1;
		}else if(error_nb==2){
			//cerr << "No symmetrical sites (binding + DSB)" << endl;
			return -2;
		}
	}
	q=double(vect_CO.size())/(vectsitedsb.size());
	return q;
}

//
vector<int> Model::get_one_gamete(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom){ //return un vector de vector
	double q;
	int ind_gen=2;
	vector<int> zygote{genotype_indiv[0]};
	if(genotype_indiv[0]!=genotype_indiv[1]){ // if we give only one indiv don't need gen
		if(zygosity_){
			ind_gen=1;
		}
		zygote.push_back(genotype_indiv[1]);
	}
	vector<vector<int>>summarysites;
	vector<int>Z;
	vector<int> vectsites;
	int nblinksite = 0;
	for(auto z : zygote){
		for(auto i : Siteforeacheallele_[z]){
			vector<int>linkedsites=vector<int>(5,0);
			linkedsites[0]=i;
			double p_occup=ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]);
			bool islinked=false;
			for(int j=0; j<4; j++){
				int chrom=j/2;
				if(indiv_chrom[chrom][i]==1){ // if we give only one indiv don't need pop
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
	double pDSB=double(nbDSB_)/nblinksite;
	vector<vector<int>>vectsitedsb;
	vector<vector<int>> vect_CO;
	//vector<double>alleleDSB;
	//vector<double>alleleCO;
	for(int i=0; i<summarysites.size(); i++){
		vector<int> link = vectfreesites(vector<int>(summarysites[i].begin()+1, summarysites[i].end()), 1);
		int nbdsbpersite=0;
		bool dsb=false;
		for(auto j : link){
			if(bernoulli_draw(pDSB)){
				dsb=true;
				summarysites[i][j+1]=2;
				nbdsbpersite+=1;
				vectsitedsb.push_back({i,j});
				//alleleDSB.push_back(Z[i]);
				try{
					if (nbdsbpersite>1 and withDSB_){
						throw int(0);
					}
				}
				catch(int const& error_nb){
					if(error_nb==0){
						//cerr << "2 DSB on one site" << endl;
					}
					return {-1};
				}
			}
		}
		vector<int> vco;
		if (dsb){
			for(int indexnbdsb = 0; indexnbdsb<nbdsbpersite; indexnbdsb++){
				vco = {summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]};
				if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==0 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==1){
					if((summarysites[i][3]==1 or summarysites[i][3]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=2){
						vco.push_back(2);
					}
					if((summarysites[i][4]==1 or summarysites[i][4]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=3){
						vco.push_back(3);
					}
				}else if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==2 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==3 or summarysites[i][1]==2){
					if((summarysites[i][1]==1 or summarysites[i][1]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=0 ){
						vco.push_back(0);
					}
					if((summarysites[i][2]==1 or summarysites[i][2]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=1){
						vco.push_back(1);
					}
				}
				if(vco.size()>2){
					vect_CO.push_back(vco);
					//alleleCO.push_back(Z[i]);
				}
			}
		}
	}
	try{
		if(vectsitedsb.size()==0){
			throw int(1);
		}else if(not vect_CO.size()){
			throw int(2);
		}
	}
	catch(int const& error_nb){
		if(error_nb==1){
			//cerr << "No DSB" << endl;
		}else if(error_nb==2){
			//cerr << "No symmetrical sites (binding + DSB)" << endl;
		}
		return {-1};
	}
	q=double(vect_CO.size())/(vectsitedsb.size());
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
	int no_chromatide=choose(4);
	int no_current_chrom = 0;
	int no_homologue_chrom = 1;
	if(no_chromatide==2 or no_chromatide==3){
		no_current_chrom= 1;
		no_homologue_chrom= 0;
	}
	vector<int> new_indiv = vector<int>(L_,0);
	if(no_chromatide==index_CO[1] or no_chromatide==index_CO[2]){
		copy( indiv_chrom[no_current_chrom].begin(), indiv_chrom[no_current_chrom].begin()+index_CO[0], new_indiv.begin() );
		copy( indiv_chrom[no_homologue_chrom].begin()+index_CO[0], indiv_chrom[no_homologue_chrom].end(), new_indiv.begin()+index_CO[0] );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide and summarysites[index_dsb[0]][0]<index_CO[0]){
				new_indiv[summarysites[index_dsb[0]][0]]=indiv_chrom[no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}else if((index_dsb[1]==index_CO[1] or index_dsb[1]==index_CO[2]) and summarysites[index_dsb[0]][0]>index_CO[0]){
				new_indiv[summarysites[index_dsb[0]][0]]=indiv_chrom[no_current_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}else{
		copy( indiv_chrom[no_current_chrom].begin(), indiv_chrom[no_current_chrom].end(), new_indiv.begin() );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide){
				new_indiv[summarysites[index_dsb[0]][0]]=indiv_chrom[no_homologue_chrom][summarysites[index_dsb[0]][0]];				
			}
		}
	}
	return new_indiv;
}



//Fonction calcul du q et fitness de chaque individu de la population (non analytique)




map<int,vector<double>> Model:: q_fert_individual_analytique(vector<vector<int>>* genotype, vector<vector<vector<int>>>* pop){
//Fonction calcul du q et fitness de chaque individu de la population (analytique)
/*
On donne en entree l'individu, son genotype PRDM9 et les sites associes (leur distribution daffinite) aux alleles de l'indiv
Pour chaque site on calcul le x=cy/1+cy qu'on moyenne sur tous les sites de l'allele.
On fait de meme pour x^2=(cy/1+cy)^2 et x^3=(cy/1+cy)^3
Puis on fait le calcul suivant :
- si l'individu etudie est homozygote : q_hom=(2<x^2>-<x^3>)/<x>
- si l'individu etudie est heterozygote : q_het=(2<x^2>_z1-<x^3>_z1+2<x^2>_z2-<x^3>_z2)/(<x>_z1+<x>_z2)
Ensuite, on calcul la fitness de l'individus :
- si l'individu etudie est homozygote : w_hom=1-exp(-dq_hom)
- si l'individu etudie est heterozygote : w_het=1-exp(-dq_het)
*/	
	//vector<double> res;
	map<int,vector<double>> res;
	map<int,vector<double>> qmoy;
	map<int,vector<double>> fertmoy;
	double qindmoy=0;
	double fertindmoy=0;
	for (int indiv=0; indiv<2*N_; indiv=indiv+2){ // for each individual in the population
		bool homoz=true;
		vector<int>genotype_indiv={(*genotype)[parityIndex_][indiv]}; // suppose that the individual is homozygous : add the first allele of the individual's genotype
		if((*genotype)[parityIndex_][indiv]!=(*genotype)[parityIndex_][indiv+1]){ //if it is heterozygous : add the second allele
			genotype_indiv.push_back((*genotype)[parityIndex_][indiv+1]);
			homoz=false; // it is not homozygous
		}
		map<int,vector<double>> infoallele;
		double qal=0;
		for(int i=0; i<genotype_indiv.size(); i++){
			int ind_gen=2;
			double moyprobalink=0; //cy/(1+cy)
			double moyprobalink2=0; //(cy/(1+cy))^2
			double moyprobalink3=0; //(cy/(1+cy))^3
			double nbactivesite=0;
			for(auto const &it : Siteforeacheallele_[genotype_indiv[i]]){ // for each site
				if((*pop)[parityIndex_][indiv][it]==1 or (*pop)[parityIndex_][indiv+1][it]==1){
					nbactivesite+=1;
					if(zygosity_){
						if(homoz){
							moyprobalink+=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
							moyprobalink2+=puissance_double(2, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
							moyprobalink3+=puissance_double(3, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
						}else{
							ind_gen=1;
							moyprobalink+=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
							moyprobalink2+=puissance_double(2, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
							moyprobalink3+=puissance_double(3, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
						}
					}else{
						moyprobalink+=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
						moyprobalink2+=puissance_double(2, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
						moyprobalink3+=puissance_double(3, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
					}
				}
			}
			//cout<<nbactivesite<<endl;
			infoallele[genotype_indiv[i]].push_back(moyprobalink/nbactivesite);/////////////////// ???????????
			infoallele[genotype_indiv[i]].push_back(moyprobalink2/nbactivesite);/////////////////// ???????????
			infoallele[genotype_indiv[i]].push_back(moyprobalink3/nbactivesite);/////////////////// ???????????
				/*if(zygosity_){
					if(homoz){
						moyprobalink+=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
						moyprobalink2+=puissance_double(2, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
						moyprobalink3+=puissance_double(3, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
					}else{
						ind_gen=1;
						moyprobalink+=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
						moyprobalink2+=puissance_double(2, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
						moyprobalink3+=puissance_double(3, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
					}
				}else{
					moyprobalink+=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
					moyprobalink2+=puissance_double(2, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
					moyprobalink3+=puissance_double(3, ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]));
				}
				
			}
			// for this allele add cy/(1+cy), (cy/(1+cy))^2, and (cy/(1+cy))^3
			infoallele[genotype_indiv[i]].push_back(moyprobalink/nbsite_);
			infoallele[genotype_indiv[i]].push_back(moyprobalink2/nbsite_);
			infoallele[genotype_indiv[i]].push_back(moyprobalink3/nbsite_);*/
		}
		if(zygosity_){
			if(homoz){
				qal=(2*infoallele[genotype_indiv[0]][1]-infoallele[genotype_indiv[0]][2])/(infoallele[genotype_indiv[0]][0]);
			}else{
				qal=(2*infoallele[genotype_indiv[0]][1]+2*infoallele[genotype_indiv[1]][1]-infoallele[genotype_indiv[0]][2]-infoallele[genotype_indiv[1]][2])/(infoallele[genotype_indiv[0]][0]+infoallele[genotype_indiv[1]][0]);
			}
		}else{
			if(homoz){
				qal=(2*infoallele[genotype_indiv[0]][1]-infoallele[genotype_indiv[0]][2])/(infoallele[genotype_indiv[0]][0]);
			}else{
				qal=(2*infoallele[genotype_indiv[0]][1]+2*infoallele[genotype_indiv[1]][1]-infoallele[genotype_indiv[0]][2]-infoallele[genotype_indiv[1]][2])/(infoallele[genotype_indiv[0]][0]+infoallele[genotype_indiv[1]][0]);
			}
		}
		if(zygosity_){
			if(homoz){
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
			}else{
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy[genotype_indiv[1]].push_back(qal);
				fertmoy[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal));
			}
		}else{
			if(homoz){
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
			}else{
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy[genotype_indiv[1]].push_back(qal);
				fertmoy[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal));
			}
		}
		qindmoy+=qal;
		fertindmoy+=(1-exp(-nbDSB_*qal));
	}
	//res.push_back(qmoy/N_); //pop
	//res.push_back(fertmoy/N_); //pop
	res[-1].push_back(qindmoy/N_);
	res[-1].push_back(fertindmoy/N_);
	double qmoyallpop=0;
	double fertmoyallpop=0;
	for(auto const &it : qmoy){
		//res[it.first].push_back(qmoy[it.first].mean());
		//res[it.first].push_back(fertmoy[it.first].mean());
		double qmoyall=0;
		double fertmoyall=0;
		for(int i=0; i<qmoy[it.first].size(); i++){
			qmoyall+=qmoy[it.first][i];
			fertmoyall+=fertmoy[it.first][i];
		}
		res[it.first].push_back(qmoyall/ qmoy[it.first].size());
		res[it.first].push_back(fertmoyall/ fertmoy[it.first].size());
		qmoyallpop+=qmoyall;
		fertmoyallpop+=fertmoyall;
	}
	res[-2].push_back(qmoyallpop/(2*N_));
	res[-2].push_back(fertmoyallpop/(2*N_));
	res[-3].push_back(0);
	res[-3].push_back(0);
	return res;
}


//Fonction calcul du q et fitness de chaque allele de la population (non analytique)


//Fonction puissance double
double Model::puissance_double(int puiss, double x){
	double res=1;
	for(int i=0; i<puiss; i++){
		res=res*x;
	}
	return(res);
}

int Model::puissance_int(int puiss, int x){
	int res=1;
	for(int i=0; i<puiss; i++){
		res=res*x;
	}
	return(res);
}


map<int,vector<double>>Model:: q_fert_allele_analytique(vector<vector<int>>* genotype){//???????????????????????????????? utile ???????????????????????????????????????
//Fonction calcul du q et fitness de chaque allele de la population (analytique)
/*
On donne en entree un allele, et ses sites associes (leur distribution daffinite)
Pour chaque site on calcul le x=cy/1+cy qu'on moyenne sur tous les sites de l'allele.
On fait de meme pour x^2=(cy/1+cy)^2 et x^3=(cy/1+cy)^3
En realite, il faut claculer prealablement tous les x, x^2 et x^3 de chaque allele de la population
Puis on fait les calcul suivant :
- pour prendre en compte le cas homozygote, on calcul : q_hom=(2<x^2>-<x^3>)/<x>
- pour prendre en compte tous les cas heterozygote (contenant l'allele en question), on calcul : q_het=(2<x^2>_z1-<x^3>_z1+2<x^2>_z2-<x^3>_z2)/(<x>_z1+<x>_z2) -> ou z_1 est l'allele etudie et z_2 est l'un des autres alleles de la popuation
=> le q total de l'allele est donc la moyenne de tous les q associes a cet allele pondere par leur frequence sous Hardy-Weinberg : f_hom=fallele^2 et f_het=fallele*fallele2
Ensuite, on calcul la fitness de l'allele :
- pour prendre en compte le cas homozygote, on calcul : w_hom=1-exp(-dq_hom)
- pour prendre en compte tous les cas heterozygote (contenant l'allele en question), on calcul : w_het=1-exp(-dq_het)
=> le w total de l'allele est donc la moyenne de tous les q associes a cet allele pondere par leur frequence sous Hardy-Weinberg : f_hom=fallele^2 et f_het=fallele*fallele2
*/
	int ind_gen=2;
	map<int,vector<double>> infoallele;
	double qal=0;
	double qmoy=0;
	double fertmoy=0;
	for (auto const &it : Siteforeacheallele_){ // for each allele in the population
		if (it.first!=-3){
			double moyprobalink_c2=0;
			double moyprobalink2_c2=0;
			double moyprobalink3_c2=0;
			double moyprobalink_c1=0;
			double moyprobalink2_c1=0;
			double moyprobalink3_c1=0;
			double frequency_allele=0;
			frequency_allele=freqallele(it.first, genotype);
			for(auto const &i : it.second){ // for each site
				moyprobalink_c2+=ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]);
				moyprobalink2_c2+=puissance_double(2, ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]));
				moyprobalink3_c2+=puissance_double(3, ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]));
				if(zygosity_){
					ind_gen=1;
					for(auto const &i : it.second){ // for each site
						moyprobalink_c1+=ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]);
						moyprobalink2_c1+=puissance_double(2, ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]));
						moyprobalink3_c1+=puissance_double(3, ind_gen*Affinity_[i]/(1+ind_gen*Affinity_[i]));
					}
				}
			}
			infoallele[it.first].push_back(frequency_allele);
			infoallele[it.first].push_back(moyprobalink_c2/nbsite_);
			infoallele[it.first].push_back(moyprobalink2_c2/nbsite_);
			infoallele[it.first].push_back(moyprobalink3_c2/nbsite_);
			infoallele[it.first].push_back(moyprobalink_c1/nbsite_);
			infoallele[it.first].push_back(moyprobalink2_c1/nbsite_);
			infoallele[it.first].push_back(moyprobalink3_c1/nbsite_);
		}
	}
	map<int,vector<double>> results;
	for (auto const &it1 : Siteforeacheallele_){ // for each allele in the population
		double qmoy=0;
		double fertmoy=0;
		if (it1.first!=-3){
			for (auto const &it2 : Siteforeacheallele_){ // for each allele in the population
				double qal=0;
				if (it2.first!=-3){
					if(zygosity_){
						if(it1.first==it2.first){
							qal=(2*infoallele[it1.first][2]-infoallele[it1.first][3])/(infoallele[it1.first][1]);
						}else if(it1.first!=it2.first){
							qal=(2*infoallele[it1.first][5]+2*infoallele[it2.first][5]-infoallele[it1.first][6]-infoallele[it2.first][6])/(infoallele[it1.first][4]+infoallele[it2.first][4]);
						}
					}else{
						qal=(2*infoallele[it1.first][2]+2*infoallele[it2.first][2]-infoallele[it1.first][3]-infoallele[it2.first][3])/(infoallele[it1.first][1]+infoallele[it2.first][1]);
					}
					//qmoy+=qal*infoallele[it1.first][0]*infoallele[it2.first][0];
					//fertmoy+=(1-exp(-nbDSB_*qal))*infoallele[it1.first][0]*infoallele[it2.first][0];
					qmoy+=qal;
					fertmoy+=(1-exp(-nbDSB_*qal));
				}
			}
			//results[it1.first].push_back(qmoy);
			//results[it1.first].push_back(fertmoy);
			results[it1.first].push_back(qmoy/infoallele.size());
			results[it1.first].push_back(fertmoy/infoallele.size());
		}
	}
	return (results);
}





