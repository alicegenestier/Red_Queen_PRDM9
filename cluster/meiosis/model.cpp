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
#include <omp.h>
//#include <boost/math/distributions/beta.hpp>
using namespace std;
//using namespace boost::math;

// pour paralleliser une boucle, on ajoute   #pragma omp parallel for num_threads(nb)   avant la boucle for.
// Ajouter dans le makefile -fopenmp

//============================//
//         Constructors	      //
//============================//
Model::Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen, double m, double alpha, double beta, int nbgenmig, bool popsamesize, int nbloop, int nbcore, bool isallele, bool issampling, bool isanalytic, double ctot, bool targetcomp, int quantilenb, int nbmeiperind, double cfreethreshold, string name): N_(N), L_(L), nbsite_(nbsite), indPrdm9_(indPrdm9), nballele_(nballele), parityIndex_(parityIndex), v_(v), u_(u), w_(w), meanaff_(meanaff), varaff_(varaff), nbDSB_(nbDSB), nbGenerations_(nbGenerations), ismigration_(ismigration), zygosity_(zygosity), withDSB_(withDSB), everygen_(everygen), m_(m), alpha_(alpha), beta_(beta), nbgenmig_(nbgenmig), popsamesize_(popsamesize), nbloop_(nbloop), nbcore_(nbcore), isallele_(isallele), issampling_(issampling), isanalytic_(isanalytic), ctot_(ctot), targetcomp_(targetcomp), quantilenb_(quantilenb), nbmeiperind_(nbmeiperind), cfreethreshold_(cfreethreshold), name_(name) {

	//Create the model : the population of N individuals with 2 homologous chromosomes per individual, on each chromosome there is a Prdm9 locus (at position indPrdm9_) to witch is associated an allele (number >=0). Each allele recognizes some target sites along the chromosome. Each site is either active (1) or inactive (0) and posess its own binding affinity.

		if(targetcomp_){
			zygosity_=1;
		}
	
	//vector counting the number of failed meiosis per generation
	nbfailedmeiosis_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
	
	//popluations matrix
	populations_ = vector<vector<vector<int>>>(2,vector<vector<int>>(2*N_,vector<int> (L_,1)));
	// matrix initialized with 1 because at the begining all the sites are activated in the whole genome for each individual
	for (auto &i : populations_){
		#pragma omp parallel for num_threads(nbcore_)
		for (int j=0; j<2*N_; j++)
		i[j][indPrdm9_]=0; //site for all individuals the PRDM9 allele 0 at the PRDM9 locus
	}
	
	//affinity vector
	Affinity_=vector<double>(L_,0);
	#pragma omp parallel for num_threads(nbcore_)
	for(int i=0; i<L_; i++){
		Affinity_[i]=choosegamma(meanaff_, varaff_);
		//affinity chosen in gamma law
	}
	
	// Affinity categories : only if targetcomp_=1 and quantilenb>0
	if(quantilenb_!=0){
		double quantile_length=1/double(quantilenb_);
		double quantile=0;
		double meanquantile=0;
		for (int i=0; i<quantilenb_; i++){
			if(quantile+quantile_length<1){ // to compute the mean affinity between these 2 quantiles
				meanquantile=double((-log(1-quantile)/meanaff_)+(-log(1-(quantile+quantile_length))/meanaff_))/double(2);
			}else{
				meanquantile=double((-log(1-quantile)/meanaff_)+(-log(1-(0.99999999999))/meanaff_))/double(2);
			}
			quantile=quantile+quantile_length;
			if(quantile+quantile_length<1){ // to set all the quantiles
				nbsitesperquantile_[-log(1-quantile)/meanaff_]={meanquantile,0};
			}else{
				nbsitesperquantile_[-log(1-0.99999999999)/meanaff_]={meanquantile,0};
			}
		}
		//Print for each quantile the mean affinity in each category and the number of site attributed to this category (here 0 because initiation)
		/*for (auto const it : nbsitesperquantile_){
			cout<<it.first<<" => ";
			for(auto const &i : it.second){
				cout<<" "<<i;
			}
			cout<<endl;
		}
		cout<<'\n';*/
	}
	
	// genotypes vector
	genotypes_=vector<vector<int>>(2,vector<int>(2*N_,0));
	
	// Prdm9 position in the genome
	Alleleforeachpos_ = vector<int>(L_,-1);
	Alleleforeachpos_[indPrdm9_]=-2;
	
	
	// neutral positions
	vector<int> freeposneutral = vectfreesites(Alleleforeachpos_, -1);
	vector<int> firstposneutral = choosemany(nbsite_, freeposneutral); 
	
	for(auto i : firstposneutral){
		Alleleforeachpos_[i]=-3;
	}
	
	
	// positions of sites for first allele (0)
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> firstpos = choosemany(nbsite_, freepos); 
	
	for(auto i : firstpos){
		Alleleforeachpos_[i]=0;
	}
	
	// positions map
	Siteforeacheallele_[-2]={indPrdm9_};
	#pragma omp parallel for num_threads(nbcore_)
	for (int i = 0; i < nballele_; i++){
		Siteforeacheallele_[i]=firstpos;
		Ageallele_[i]=freqall(i, &genotypes_, &populations_)*(N_*v_*nbDSB_)/(2*nbsite_);
		double meanaffall = get_mean_affinity(i,&populations_);
		infoperallele_[i]={0,0,0,0,0,0,meanaffall,0,0,0,0}; // infoperallele{nb failed meiosis, 2 DSB on 1 site, no DSB, no symetrical site, q, nb meiosis with at least 1 DSB, mean affinity of the allele, relative_nb_meiosis, nblinksiteall, nblinkposall, absolute_nb_meiosis}
		
		infoperallele_hom_[i]={0,0,0,0,0,0,meanaffall,0,0,0,0,0,0,0}; // infoperallele_hom{nb failed meiosis, 2 DSB on 1 site, no DSB, no symetrical site, q, nb meiosis with at least 1 DSB, mean affinity of the allele, relative_nb_meiosis, nblinksiteall, nblinkposall, absolute_nb_meiosis, cfree_moy, cfree_tot, c_occup}
		infoperallele_het_[i]={0,0,0,0,0,0,meanaffall,0,0,0,0,0,0,0}; // infoperallele_het{nb failed meiosis, 2 DSB on 1 site, no DSB, no symetrical site, q, nb meiosis with at least 1 DSB, mean affinity of the allele, relative_nb_meiosis, nblinksiteall, nblinkposall, absolute_nb_meiosis, cfree_moy, cfree_tot, c_occup}
	}
	
	//distribution beta for neutral sites
	Siteforeacheallele_[-3]=firstposneutral;
	for(int site=0; site<nbsite_; site++){
		double freqactivsite = choosebeta(alpha_, beta_);
		for(int ind=0; ind<2*N_; ind++){
			double p=bernoulli_draw(freqactivsite);
			if(not p){
				populations_[parityIndex_][ind][Siteforeacheallele_[-3][site]]=0;
			}
		}
	}
	
	q_=0;
	qsym_=0;
	qnum_=0;
	qdenom_=0;
	
}

//============================//
//        Destructors         //
//============================//
	
//============================//
//           Getters          //
//============================//
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

map<int,vector<double>> Model::infoperallele_hom(){ 
	return infoperallele_hom_;
}
/*
map<int,vector<double>> Model::infoperallele1_hom(){
	return infoperallele1_hom_;
}
map<int,vector<double>> Model::infoperallele2_hom(){
	return infoperallele2_hom_;
}
*/

map<int,vector<double>> Model::infoperallele_het(){ 
	return infoperallele_het_;
}
/*
map<int,vector<double>> Model::infoperallele1_het(){
	return infoperallele1_het_;
}
map<int,vector<double>> Model::infoperallele2_het(){
	return infoperallele2_het_;
}
*/
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
int Model::nbcore(){
	return nbcore_;
}
bool Model::isallele(){
	return isallele_;
}
bool Model::issampling(){
	return issampling_;
}
bool Model::isanalytic(){
	return isanalytic_;
}
bool Model::targetcomp(){
	return targetcomp_;
}
double Model::qnum(){
	return qnum_;
}
double Model::qdenom(){
	return qdenom_;
}
double Model::ctot(){
	return ctot_;
}

int Model::quantilenb(){
	return quantilenb_;
}
int Model::nbmeiperind(){
	return nbmeiperind_;
}

//============================//
//           Setters          //
//============================//


//============================//
//           Methods          //
//============================//

//----------------------------//
// 	choose methods        //
//----------------------------//

static random_device rdev{};
static default_random_engine e{rdev()};

//1) uniforme
int Model::choose(int n)   {
   uniform_int_distribution<int> d(0,n-1);
   return d(e);
}

//2) bernoulli
int Model::bernoulli_draw(double p){
	bernoulli_distribution distribution(p);
	return distribution(e);
}

//3) binomiale
int Model::binomial_draw(int n, double p)  {
    binomial_distribution<int> b(n,p);
    return b(e);
}

//4) gamma
double Model::choosegamma(double meanaff, double varaff){
	gamma_distribution<double> g(varaff,meanaff);
	return g(e);
}

//5) beta
double Model::choosebeta(double alpha, double beta){
	double X = choosegamma(alpha, 1);
	double Y = choosegamma(beta, 1);
	return X/(X+Y);
}


vector<int> Model::choosemany(int k, vector<int> vect){
//choose k positions among all the free positions (vect)
//return the vector of index of the chosen positions
	vector<int> newsites;
	int n = vect.size(); //number of free positions
	try{
		if (k>n){
			throw string("Not enough free positions in the genome for new allele");
		}
	} // assert
	catch(string const& chaine){
		cerr << chaine << endl;
	}
	for (int i=0; i<k; i++){ //for i in k (number of elements we need)
		int j=choose(n-i);//pick randomly an index in n minus the i last elements
		newsites.push_back(vect[j]);//add vect[j] in newsites
		vect[j]=vect[n-i-1];//the value of index j (just chosen) have the value of the n-i-1 index => the last element that has not been chosen replace the element chosen
	}
	sort(newsites.begin(), newsites.end());
	return(newsites);
}

//----------------------------//
// 	other methods         //
//----------------------------//
vector<int> Model::vectfreesites(vector<int> vect, int nb){
	//return the index of all positions equal to nb
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
	//int num1=0;
	//int num2=0;
	for(int i=0; i<vect.size(); i++){
		if(vect[i]!= -2 and vect[i]!= -1 and vect[i]!= -3){
			//num1+=1;
			//vector<int>::iterator itv = find((*genotype)[parityIndex_].begin(),(*genotype)[parityIndex_].end(),vect[i]);/// ??????????
			//if(itv!=(*genotype)[parityIndex_].end()){
				//num2+=1;
				occupsites.push_back(i);
			//}
		}else if(vect[i]==-3){
			occupsitesneutral.push_back(i);
		}
	}
	//if(num1!=num2){
		//cout<<"FALSE"<<endl;
	//}
	return (vector<vector<int>>{occupsites,occupsitesneutral});	
}

void Model::sitemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype){
	// for each position with allele, prob v to mutate and if mutation, choose randomly 1 chrom to mutate
	vector<vector<int>> occupsite = occupiedsites(Alleleforeachpos_, genotype);
	#pragma omp parallel for num_threads(nbcore_)
	for (int ind=0; ind<occupsite[0].size(); ind++){ //for all the sites recognized by an allele
		int i=occupsite[0][ind];
		if (bernoulli_draw(2*N_*v_)){ //if mutation at this site
			int mutchrom = choose(2*N_); //only one indiv will preform a mutation at this site
			//if ((*population)[parityIndex_][mutchrom][i]==1){ //inutile
			(*population)[parityIndex_][mutchrom][i]=0; //perform the mutation
			//}
		}
	}
	#pragma omp parallel for num_threads(nbcore_)
	for (int ind=0; ind<occupsite[1].size(); ind++){ //for all the neutral sites
		int i=occupsite[1][ind];
		if (bernoulli_draw(2*N_*w_)){ //if mutation at this site
			int mutchrom = choose(2*N_); //only one indiv will preform a mutation at this site
			(*population)[parityIndex_][mutchrom][i]=((*population)[parityIndex_][mutchrom][i]+1)%2; //perform the mutation active->inacive or inactive->active
		}
	}
}

void Model::allelemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,double>* Ageallele, vector<int> vectnbind){
	// for each chromosome, mutation of allele with proba u and if at least one mutation update populations (change prdm9 allele at its position), genotypes(same), map (site pos associated to the new allele) and allele for each pos (change the allele corresponding to each pos : -1 -> new allele)
	int k = binomial_draw(2*N_, u_);
	vector<int> alltomutate = choosemany(k, vectnbind);// vector de 0 a 2N-1
	for(int i=0; i<k;i++){
		nballele_+=1;
		vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
		vector<int> newpos = choosemany(nbsite_, freepos); 
		int newallele = nballele_-1;
		for(auto j : newpos){
			Alleleforeachpos_[j]=newallele;//update allele for each pos
		}
		Siteforeacheallele_[newallele]=newpos; //update map
		(*genotype)[parityIndex_][alltomutate[i]]=newallele;//update genotype
		(*population)[parityIndex_][alltomutate[i]][indPrdm9_]=newallele;//update pop
		(*Ageallele)[newallele]=freqall(newallele, genotype, population)*(N_*v_*nbDSB_)/(2*nbsite_);//update age of the new allele
		if(ismigration_){
			infoperallele1_[newallele]={0,0,0,0,0,0,0,0,0,0,0};
			infoperallele2_[newallele]={0,0,0,0,0,0,0,0,0,0,0};
			/*
			infoperallele1_hom_[newallele]={0,0,0,0,0,0,0,0,0,0};
			infoperallele2_hom_[newallele]={0,0,0,0,0,0,0,0,0,0};
			infoperallele1_het_[newallele]={0,0,0,0,0,0,0,0,0,0};
			infoperallele2_het_[newallele]={0,0,0,0,0,0,0,0,0,0};
			*/
		}else{
			infoperallele_[newallele]={0,0,0,0,0,0,0,0,0,0,0};
			
			infoperallele_hom_[newallele]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			infoperallele_het_[newallele]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			
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
			infoperallele_hom_.erase(a);
			infoperallele_het_.erase(a);
		}else{
			Siteforeacheallele_.erase(a);
			Ageallele1_.erase(a);
			infoperallele1_.erase(a);
			//infoperallele1_hom_.erase(a);
			//infoperallele1_het_.erase(a);
			Ageallele2_.erase(a);
			infoperallele2_.erase(a);
			//infoperallele2_hom_.erase(a);
			//infoperallele2_het_.erase(a);
		}
	}
}

//----------------------------//
//	Print functions       //
//----------------------------//

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

//----------------------------//
// 	Meiosis method        //
//----------------------------//

// Meiosis
int Model::Meiosis(int no_chrom_ind, int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het, vector<vector<int>>* nbfailedmeiosis, double* q, double* qsym, double* qnum, double* qdenom, int indiv/*, int nb_meiosis*/){
//Perform the meiosis of one individual and return 0 if the meiosis fails somewhere or 1 if the meiosis is successfull

	//----------------------------//
	//Initialisation time meiosis //
	//----------------------------//
	double tps_binding, tps_DSB, tps_vCO, tps_symq, tps_CONCO, tps_tot;
	clock_t t_1, t_2, t_3, t_4, t_5, t_6;
	t_1=clock();
	
	//----------------------------//
	//Initialisation of zygot vect//
	//----------------------------//
	//creation of the vector zygote containing the 1 or 2 allele(s) of the individual
	double ind_gen=1; //if no target competition

	vector<int> zygote{(*genotype)[parityIndex_][indiv]};
	
	/*for (auto zyg : zygote){
		cout<<zyg<<" ";
	}
	cout<<"\n"<<endl;*/

	
	if((*genotype)[parityIndex_][indiv]==(*genotype)[parityIndex_][indiv+1]){//for homozygous
		//cout<<"homozygous"<<endl;
		if(zygosity_){ // if genetic dosage
			ind_gen=2; //if no target competition
		}
		(*infoperallele)[zygote[0]][7]+=2; //relative nb meiosis realised per allele (*2 if hom)
		(*infoperallele)[zygote[0]][10]+=1; //absolute nb meiosis realised per allele
		
		(*infoperallele_hom)[zygote[0]][7]+=2; //relative nb de meiose realise per allele
		(*infoperallele_hom)[zygote[0]][10]+=1; //absolute nb meiosis realised per allele
		
	}
	if((*genotype)[parityIndex_][indiv]!=(*genotype)[parityIndex_][indiv+1]){//for heterozygous
		//cout<<"heterozygous"<<endl;
		zygote.push_back((*genotype)[parityIndex_][indiv+1]);
		(*infoperallele)[zygote[0]][7]+=1; //relative nb meiosis realised per allele 1
		(*infoperallele)[zygote[1]][7]+=1; //relative nb meiosis realised per allele 2
		(*infoperallele)[zygote[0]][10]+=1; //relative nb meiosis realised per allele 1
		(*infoperallele)[zygote[1]][10]+=1; //relative nb meiosis realised per allele 2
		
		(*infoperallele_het)[zygote[1]][7]+=1; //relative nb de meiose realise par l'allele 2
		(*infoperallele_het)[zygote[0]][7]+=1; //relative nb de meiose realise par l'allele 1
		(*infoperallele_het)[zygote[0]][10]+=1; //absolute nb meiosis realised per allele 1
		(*infoperallele_het)[zygote[1]][10]+=1; //absolute nb meiosis realised per allele 2
		
	}
	
	
	//----------------------------//
	//    If target competition   //
	//----------------------------//
	/*if (zygote.size()==1){
		cout<<"infoperallele_hom_before"<<endl;
		printinfoallele(infoperallele_hom);
	}else{
		cout<<"infoperallele_het_before"<<endl;
		printinfoallele(infoperallele_het);
	}*/
	if(targetcomp_){//option param if we want to take into account target competition
	
		vector<map<int,vector<double>>*> vectinfo_hom_het = vector<map<int,vector<double>>*> {infoperallele_hom, infoperallele_het};
		for(auto z : zygote){
			double is_diff=1;
			double ctot=ctot_;
			double cfree;
			double threshold=cfreethreshold_;
			double index_vectinfo_hom_het=1;
			double sum_p_occup=0;
			if (zygote.size()==1){
				//cout<<"homozygote"<<endl;////////////////////////
				if(zygosity_){
					ctot=2*ctot_;
				}
				index_vectinfo_hom_het=0;
			}
			//else{//////////////////////
				//cout<<"heterozygote"<<endl;//////////////////
			//}/////////////////////
			if((*vectinfo_hom_het[index_vectinfo_hom_het])[z][11]==0){
				cfree=1/float(1000)*ctot;
			}else{
				cfree=(*vectinfo_hom_het[index_vectinfo_hom_het])[z][11];
			}
			//Two cases possible :
			//1) If we want to take the affinity categories map (quantilenb_!=0),if the site is active, the category in which the affinity of this site enters gains one site
			for(auto quantile_index : nbsitesperquantile_){
				nbsitesperquantile_[quantile_index.first][1]=0;
			}
			if(quantilenb_!=0){
				for(auto i : Siteforeacheallele_[z]){
					for(int j=0; j<4; j++){
						int chrom=indiv+j/2;
						if((*population)[parityIndex_][chrom][i]==1){
							bool is_quantile=0;
							for(auto quantile_index : nbsitesperquantile_){
								if(Affinity_[i]<=quantile_index.first and is_quantile==0){
									nbsitesperquantile_[quantile_index.first][1]+=1;//if this site is active, the category in which the affinity of this site enters gains one site
									is_quantile=1;
								}
							}
						}
					}
				}
			}
			double cfree_min=0;
			double cfree_max=ctot;
			while(is_diff==1){
				is_diff=0;
				sum_p_occup=0;//the mean sum of occupied site
				//---------------
				//sum_p_occup
				//---------------
				//In the case of affinity categories
				if(quantilenb_!=0){
					for(auto quantile_index : nbsitesperquantile_){
						sum_p_occup+=double(nbsitesperquantile_[quantile_index.first][1]*(cfree*nbsitesperquantile_[quantile_index.first][0])/(1+cfree*nbsitesperquantile_[quantile_index.first][0]));//The probability of binding in each category multiplied by the number of sites in this category is added to the mean sum of occupied site
					}
				}

				//In the other case
				else{
					for(auto i : Siteforeacheallele_[z]){//for each position
						double p_occup=cfree*Affinity_[i]/(1+cfree*Affinity_[i]);//the probability of binding is calculated
						for(int j=0; j<4; j++){//for each site of this position
							int chrom=indiv+j/2;
							if((*population)[parityIndex_][chrom][i]==1){//if the site is active, its probability of binding is added to the mean sum of occupied site
								sum_p_occup+=p_occup;
							}
						}
					}
				}
				if(cfree+sum_p_occup>ctot){
					//if cfree is to high
					is_diff=1;
					cfree_max=cfree;
					cfree=double(cfree_max+cfree_min)/2;
				}else{
					//if cfree is to low or the right one
					double new_cfree=ctot-sum_p_occup;
					//cout<<"cfree/ctot : "<<double(cfree)/double(ctot)<<endl;
					//cout<<"cfree/ctot_diff : "<<abs(double(cfree)/double(ctot)-double(new_cfree)/double(ctot))<<endl;
					if(abs(double(cfree)/double(ctot)-double(new_cfree)/double(ctot))<=threshold){
						//the difference between the current and the next cfree/ctot is below or equal to the thershold
						is_diff=0;
					}else{
						//the difference between the current and the next cfree/ctot is above the thershold
						is_diff=1;
						cfree_min=cfree;
						cfree=double(cfree_max+cfree_min)/2;
					}
				}
			}
			(*vectinfo_hom_het[index_vectinfo_hom_het])[z][11]=cfree;
			(*vectinfo_hom_het[index_vectinfo_hom_het])[z][12]+=cfree;
			(*vectinfo_hom_het[index_vectinfo_hom_het])[z][13]+=sum_p_occup;
		}
	}
	/*if (zygote.size()==1){
		cout<<"infoperallele_hom_after"<<endl;
		printinfoallele(infoperallele_hom);
	}else{
		cout<<"infoperallele_het_after"<<endl;
		printinfoallele(infoperallele_het);
	}*/
	
	

	
	//----------------------------//
	//	PRDM9 binding         //
	//----------------------------//
	vector<vector<int>>summarysites;
	vector<int>Z;
	vector<int> vectsites;
	int nblinksite = 0;
	int nblinksitesym = 0;
	
	for(auto z : zygote){ // for each allele
		//If we want to study target competition
		if(zygote.size()==1 and targetcomp_){
			ind_gen=(*infoperallele_hom)[z][11];
			//cout<<"ind_gen_hom : "<<ind_gen<<endl;
		}
		else if(zygote.size()==2 and targetcomp_){	
			ind_gen=(*infoperallele_het)[z][11];
			//cout<<"ind_gen_het : "<<ind_gen<<endl;
		}
		//cout <<ind_gen<<endl;
		double nblinksiteall = 0; //
		double nblinkposall = 0; //
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
						nblinksiteall+=1; //
					}
				}
			}
			if(islinked==true){
				vectsites.push_back(i);
				nblinkposall+=1; //
				summarysites.push_back(linkedsites);
				Z.push_back(z);
				if (linkedsites[1]==1 and linkedsites[3]==1){
					nblinksitesym++;
				}
				if (linkedsites[1]==1 and linkedsites[4]==1){
					nblinksitesym++;
				}
				if (linkedsites[2]==1 and linkedsites[3]==1){
					nblinksitesym++;
				}
				if (linkedsites[2]==1 and linkedsites[4]==1){
					nblinksitesym++;
				}
			}
		}
		(*infoperallele)[z][8]+=nblinksiteall; //
		(*infoperallele)[z][9]+=nblinkposall; //
		//cout<<"nblinksiteall : "<<nblinksiteall<<endl;
	}
	*qsym=*qsym+double(nblinksitesym)/double(nblinksite);
	t_2=clock();
	tps_binding=(float)(t_2-t_1)/CLOCKS_PER_SEC;
	//printf("tps_binding = %f\n", tps_binding);

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
				vectsitedsb.push_back({i,j});//vectsitedsb.push_back({line in summarysites of the corresponding position, site undergoing the DSB (1,2,3 or 4)})
				alleleDSB.push_back(Z[i]);
				try{
					if (nbdsbpersite>1 and withDSB_){
						//if (nb_meiosis==nbmeiperind_-1){
						(*nbfailedmeiosis)[nb_gen][0]+=1; 
						(*infoperallele)[Z[i]][0]+=1;
						(*infoperallele)[Z[i]][1]+=1;
						
						if (zygote.size()==1){
							(*infoperallele_hom)[Z[i]][0]+=1;
							(*infoperallele_hom)[Z[i]][1]+=1;
						}else if(zygote.size()==2){
							(*infoperallele_het)[Z[i]][0]+=1;
							(*infoperallele_het)[Z[i]][1]+=1;
						}
						//}
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
		//t_3=clock();
		//tps_DSB=(float)(t_3-t_2)/CLOCKS_PER_SEC;
		//printf("tps_DSB = %f\n", tps_DSB);
		
		//Crossing over
		vector<int> vco;
		if (dsb){
			for(int indexnbdsb = 0; indexnbdsb<nbdsbpersite; indexnbdsb++){// for all the dsb performed at this position (0 to 4)
				vco = {summarysites[i][0],vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]};//vco={position,last indexnbdsb DSB added to vectsitedsb=all DSB performed at this site}
				if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==0 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==1){//if there is a DSB on the first or second chromatide
					//if((summarysites[i][3]==1 or summarysites[i][3]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=2){
					if(summarysites[i][3]==1 or summarysites[i][3]==2){//if there is a PRDM9 linked or even a DSB on the third chromatide
						// in these cases, I suppose that 2 DSB can perform a CO
						vco.push_back(2);//add third chromatide to the potential site for CO
					}
					//if((summarysites[i][4]==1 or summarysites[i][4]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=3){
					if(summarysites[i][4]==1 or summarysites[i][4]==2){//if there is a PRDM9 linked or even a DSB on the fourth chromatide
						vco.push_back(3);//add fourth chromatide to the potential site for CO
					}
					//that means that 1 DSB = 1 potential CO (even if there is 2 linked sites on the other side)
				//}else if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==2 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==3 or summarysites[i][1]==2){
				}else if(vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==2 or vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]==3){//same thing but with DSB on third or fourth chromatide
					//if((summarysites[i][1]==1 or summarysites[i][1]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=0){
					if(summarysites[i][1]==1 or summarysites[i][1]==2){
						vco.push_back(0);
					}
					//if((summarysites[i][2]==1 or summarysites[i][2]==2) and vectsitedsb[vectsitedsb.size()-indexnbdsb-1][1]!=1){
					if(summarysites[i][2]==1 or summarysites[i][2]==2){
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
	t_4=clock();
	tps_vCO=(float)(t_4-t_2)/CLOCKS_PER_SEC;
	//printf("tps_vCO = %f\n", tps_vCO);
	
	//Failed meiosis
	try{
		//If no DSB
		if(vectsitedsb.size()==0){
			//if (nb_meiosis==nbmeiperind_-1){
			(*nbfailedmeiosis)[nb_gen][1]+=1;
			if(zygote.size()==1){
				(*infoperallele)[zygote[0]][0]+=2;//////////////////1 si zygosity_?
				(*infoperallele)[zygote[0]][2]+=2;//////////////////1 si zygosity_?
				
				(*infoperallele_hom)[zygote[0]][0]+=2;//////////////////1 si zygosity_?
				(*infoperallele_hom)[zygote[0]][2]+=2;//////////////////1 si zygosity_?
				
			}else if(zygote.size()==2){
				(*infoperallele)[zygote[0]][0]+=1;
				(*infoperallele)[zygote[0]][2]+=1;
				(*infoperallele)[zygote[1]][0]+=1;
				(*infoperallele)[zygote[1]][2]+=1;
				
				(*infoperallele_het)[zygote[0]][0]+=1;
				(*infoperallele_het)[zygote[0]][2]+=1;
				(*infoperallele_het)[zygote[1]][0]+=1;
				(*infoperallele_het)[zygote[1]][2]+=1;
				
			}
			//}
			throw int(1);
		//If no potential CO
		}else if(not vect_CO.size()){
			//if (nb_meiosis==nbmeiperind_-1){
			(*nbfailedmeiosis)[nb_gen][2]+=1;
			if(zygote.size()==1){
				(*infoperallele)[zygote[0]][0]+=2;//////////////////1 si zygosity_?
				(*infoperallele)[zygote[0]][3]+=2;//////////////////1 si zygosity_?
				(*infoperallele)[zygote[0]][5]+=2;//////////////////1 si zygosity_?
				
				(*infoperallele_hom)[zygote[0]][0]+=2;//////////////////1 si zygosity_?
				(*infoperallele_hom)[zygote[0]][3]+=2;//////////////////1 si zygosity_?
				(*infoperallele_hom)[zygote[0]][5]+=2;//////////////////1 si zygosity_?
				
			}else if(zygote.size()==2){
				(*infoperallele)[zygote[0]][0]+=1;
				(*infoperallele)[zygote[0]][3]+=1;
				(*infoperallele)[zygote[1]][0]+=1;
				(*infoperallele)[zygote[1]][3]+=1;
				(*infoperallele)[zygote[0]][5]+=1;
				(*infoperallele)[zygote[1]][5]+=1;
				
				(*infoperallele_het)[zygote[0]][0]+=1;
				(*infoperallele_het)[zygote[0]][3]+=1;
				(*infoperallele_het)[zygote[1]][0]+=1;
				(*infoperallele_het)[zygote[1]][3]+=1;
				(*infoperallele_het)[zygote[0]][5]+=1;
				(*infoperallele_het)[zygote[1]][5]+=1;
				
			}
			//}
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
	
	//Statistics computation
	
	//*qnum=*qnum+double(vect_CO.size())
	//*qdenom=*qdenom+double((vectsitedsb.size()))
	*q=*q+double(vect_CO.size())/(vectsitedsb.size());
	if(zygote.size()==1){
		if(zygosity_){
			(*infoperallele)[zygote[0]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0]));
			(*infoperallele)[zygote[0]][5]+=2;//////////////////1 si zygosity_?
			
			(*infoperallele_hom)[zygote[0]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0]));
			(*infoperallele_hom)[zygote[0]][5]+=2;//////////////////1 si zygosity_?
			
		}else{
			(*infoperallele)[zygote[0]][4]+=2*(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0])); 
			(*infoperallele)[zygote[0]][5]+=2;
			
			(*infoperallele_hom)[zygote[0]][4]+=2*(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0])); 
			(*infoperallele_hom)[zygote[0]][5]+=2;
			
		}
	}else if(zygote.size()==2){
		if (count(alleleDSB.begin(),alleleDSB.end(),zygote[0]) != 0){
			(*infoperallele)[zygote[0]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0]));
			(*infoperallele)[zygote[0]][5]+=1;
			
			(*infoperallele_het)[zygote[0]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[0]));
			(*infoperallele_het)[zygote[0]][5]+=1;
			
			if(double(count(alleleCO.begin(),alleleCO.end(),zygote[0]))==0){
				(*infoperallele)[zygote[0]][3]+=1;
				(*infoperallele)[zygote[0]][0]+=1;
				
				(*infoperallele_het)[zygote[0]][3]+=1;
				
			}	
		}else{
			(*infoperallele)[zygote[0]][2]+=1;
			(*infoperallele)[zygote[0]][0]+=1;
			
			(*infoperallele_het)[zygote[0]][2]+=1;
			
		}
		if (count(alleleDSB.begin(),alleleDSB.end(),zygote[1]) != 0){
			(*infoperallele)[zygote[1]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[1]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[1]));
			(*infoperallele)[zygote[1]][5]+=1;
			
			(*infoperallele_het)[zygote[1]][4]+=(double(count(alleleCO.begin(),alleleCO.end(),zygote[1]))/count(alleleDSB.begin(),alleleDSB.end(),zygote[1]));
			(*infoperallele_het)[zygote[1]][5]+=1;
			
			if(double(count(alleleCO.begin(),alleleCO.end(),zygote[1]))==0){
				(*infoperallele)[zygote[1]][3]+=1;
				(*infoperallele)[zygote[1]][0]+=1;
				
				(*infoperallele_het)[zygote[1]][3]+=1;
				
			}
		}else{
			(*infoperallele)[zygote[1]][2]+=1;
			(*infoperallele)[zygote[1]][0]+=1;
			
			(*infoperallele_het)[zygote[1]][2]+=1;
			
		}
	}
	
	//Performing CO
	//choose 1 position and the 2 chromatides for CO
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
	t_5=clock();
	tps_symq=(float)(t_5-t_4)/CLOCKS_PER_SEC;
	//printf("tps_symq = %f\n", tps_symq);
	
	//choose one gamete
	int no_chromatide=choose(4);
	int no_current_chrom = indiv;
	int no_homologue_chrom = indiv+1;
	if(no_chromatide==2 or no_chromatide==3){
		no_current_chrom=indiv+1;
		no_homologue_chrom=indiv;
	}
	//performe recombination (CO and NCO) depending on the chosen gamete
	if(no_chromatide==index_CO[1] or no_chromatide==index_CO[2]){ // If the chosen gamete is the one performing the CO
		copy( (*population)[parityIndex_][no_current_chrom].begin(), (*population)[parityIndex_][no_current_chrom].begin()+index_CO[0], (*population)[(parityIndex_+1)%2][no_chrom_ind].begin() );
		copy( (*population)[parityIndex_][no_homologue_chrom].begin()+index_CO[0], (*population)[parityIndex_][no_homologue_chrom].end(), (*population)[(parityIndex_+1)%2][no_chrom_ind].begin()+index_CO[0] );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide and summarysites[index_dsb[0]][0]<index_CO[0]){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}else if((index_dsb[1]==index_CO[1] or index_dsb[1]==index_CO[2]) and summarysites[index_dsb[0]][0]>index_CO[0]){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_current_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}else{ // If the chosen gamete is not performing the CO (only NCO)
		copy( (*population)[parityIndex_][no_current_chrom].begin(), (*population)[parityIndex_][no_current_chrom].end(), (*population)[(parityIndex_+1)%2][no_chrom_ind].begin() );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide){
				(*population)[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=(*population)[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}
	t_6=clock();
	tps_CONCO=(float)(t_6-t_5)/CLOCKS_PER_SEC;
	//printf("tps_CONCO = %f\n", tps_CONCO);
	tps_tot=(float)(t_6-t_1)/CLOCKS_PER_SEC;
	//printf("tps_tot = %f\n", tps_tot);

	return 1;
}


//methode fill new population (make as many meiosis as it is necessary to fill the entire new population)
void Model::fillnewpop(int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele,map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het, vector<vector<int>>* nbfailedmeiosis, double* q, double* qsym, double* qnum, double* qdenom){
	//times settings for this function
	double tps_fill, tps_meiosis;
	clock_t tps1, tps2, tps3, tps4, tps5,tps6;
	tps1=clock();
	
	//Declaration and initialisation of parameters
	double temps_moy=0;
	int nbmei=0;
	int nbmei1=0;
	int nbmei2=0;
	*q=0;
	*qsym=0;
	*qnum=0;
	*qdenom=0;
	//int nb_meiosis=0;
	#pragma omp paralell for num_threads(nbcore_)
	for(int indnewpop=0; indnewpop<2*N_; indnewpop++){
		tps3=clock();
		int indiv = 2*choose(N_);
		int nb_meiosis=0;
		int meiosisState=-1;
		while(meiosisState==-1 and nb_meiosis<nbmeiperind_){ //allows an ind to try performing many meiosis (nbmeiperind_) before stopping and choosing another one
			meiosisState = Meiosis(indnewpop, nb_gen, population, genotype, infoperallele,infoperallele_hom,infoperallele_het, nbfailedmeiosis, q, qsym, qnum, qdenom, indiv/*, nb_meiosis*/);
			nb_meiosis+=1;
		}
		tps4=clock();
		tps_meiosis=(float)(tps4-tps3)/CLOCKS_PER_SEC;
		//printf("tps_meiosis = %f\n", tps_meiosis);
		temps_moy+=tps_meiosis;
		
		nbmei+=1; //=> count only if one of the meiosis performed for this individual successed, it doesn't count how many meiosis the individual performed before successing

		while (meiosisState==-1){
			nbmei2+=1;
			tps5=clock();
			int indiv = 2*choose(N_);
			nb_meiosis=0;
			while(meiosisState==-1 and nb_meiosis<nbmeiperind_){ //allows an ind to try performing many meiosis (nbmeiperind_) before stopping and choosing another one
				meiosisState = Meiosis(indnewpop, nb_gen, population, genotype, infoperallele,infoperallele_hom,infoperallele_het, nbfailedmeiosis, q, qsym, qnum, qdenom, indiv/*, nb_meiosis*/);
				nb_meiosis+=1;
			}
			(*nbfailedmeiosis)[nb_gen][3]+=1;
			tps6=clock();
			tps_meiosis=(float)(tps6-tps5)/CLOCKS_PER_SEC;
			//printf("tps_meiosis = %f\n", tps_meiosis);
			temps_moy+=tps_meiosis;
			
			nbmei+=1;

		}
		nbmei1+=1;
		//cout<<nbmei1<<endl;
	}
	
	temps_moy=temps_moy/nbmei;
	cout<< "temps_moy = "<<temps_moy<<endl;
	
	for(int indpop=0; indpop<2*N_;indpop++){
		(*genotype)[(parityIndex_+1)%2][indpop]=(*population)[(parityIndex_+1)%2][indpop][indPrdm9_];
	}
	//*qnum=*qnum/(2*N_+(*nbfailedmeiosis)[nb_gen][2]);
	//*qdenom=*qdenom/(2*N_+(*nbfailedmeiosis)[nb_gen][2]);
	//*q=(*qnum)/(*qdenom);
	*q=*q/(2*N_+(*nbfailedmeiosis)[nb_gen][2]);
	*qsym=*qsym/(2*N_+(*nbfailedmeiosis)[nb_gen][2]);
	
	tps2=clock();
	tps_fill=(float)(tps2-tps1)/CLOCKS_PER_SEC;
	//printf("tps_fill = %f\n", tps_fill);
	
} //////////////////////////////////////////////////////////////////////////////////////////// changer le info per allele si plusieurs meioses possibles par indiv !!!!!

//methode which mix together and runs all the previous functions during X generations 
void Model::manygenerations(){
	//Initialisation of files
	//Trace file : contains all the descriptive statistics averaged in the population for each generation
	ofstream generalfile ((name_+".trace").c_str());
	generalfile << "Generation_number" << '\t' << "Total_number_of_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' << "Mean_Age" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" << '\t' << "q_intra" << '\t' << "fertility_intra" << '\t' << "q_analytic" <<  '\t' << "fertility_analytic" << '\t' << "q_sym" << '\t' << "Fert_newall" << '\t' << "Mean_c_occup_hom" << '\t' << "Mean_c_occup_het" << '\t' << "Mean_c_free_hom" << '\t' << "Mean_c_free_het" << '\t' << "Mean_sigma" << '\n';
    generalfile.flush();
	//Allele file : contains all the mean descriptive statistics for each allele segregating in the population for each generation
    ofstream allelefile ((name_+".allele").c_str());
	allelefile << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" << '\t' << "Fertility_allele" << '\t' << "mean_affinity" << '\t' << "q_analytic" <<  '\t' << "fertility_analytic" << '\t' << "relative_nb_meiosis" <<  '\t' << "nb_linked_site" <<  '\t' << "nb_linked_pos" << '\t' << "absolute_nb_meiosis" << '\t' << "cfree_ctot_hom" << '\t' << "cfree_ctot_het" << '\t' << "q_hom" << '\t' << "q_het" << '\t' << "fertility_hom" << '\t' << "fertility_het" << '\t' << "q_hom_analytic" << '\t' << "q_het_analytique" << '\t' << "fertility_hom_analytic" << '\t' << "fertility_het_analytic" << '\t' << "q_hemi_analytic" << '\t' << "fertility_hemi_analytic" << '\t' << "sigma" <<'\n';//
    allelefile.flush();
	//Params file : contains all the parameters values for this simulation
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	double tps_cfree;
	//clock_t tps_1, tps_2;
	//tps_1=clock();
	//cout<<float(tps_1)/CLOCKS_PER_SEC<<endl;
	//double sig=0;
	double qhet_0=0;
	double qhom_0=0;
	double cfree_ctot_het_0=0;
	double cfree_ctot_hom_0=0;
	double sigma_0_c=0;
	vector<double> cf;
	#pragma omp parallel for num_threads(nbcore_)
	for(int nbloop=0; nbloop<100; nbloop++){
		//cout<<nbloop<<endl;
		cf=sigma_q_w_0();
		//sig+=sigma_0();
		qhet_0+=cf[2];
		qhom_0+=cf[3];
		cfree_ctot_het_0+=cf[0];
		cfree_ctot_hom_0+=cf[1];
		sigma_0_c+=cf[6];
	}
	//sig=sig/100;
	qhet_0=qhet_0/100;
	qhom_0=qhom_0/100;
	cfree_ctot_het_0=cfree_ctot_het_0/100;
	cfree_ctot_hom_0=cfree_ctot_hom_0/100;
	sigma_0_c=sigma_0_c/100;
	
	//tps_2=clock();
	//cout<<float(tps_2)/CLOCKS_PER_SEC<<endl;
	//tps_cfree=(float)(tps_2-tps_1)/CLOCKS_PER_SEC;
	//cout<<"tps_cfree = "<<tps_cfree<<endl;
	//printf("tps_cfree = %f\n", tps_cfree);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	ofstream paramsfile ((name_+".params").c_str());
	paramsfile << "N" << '\t' << N_ << '\n' << "L" << '\t' << L_ << '\n' << "nbsite" << '\t' << nbsite_ << '\n' << "indPrdm9" << '\t' << indPrdm9_ << '\n' << "nballele" << '\t' << nballele_ << '\n' << "parityIndex" << '\t' << parityIndex_ << '\n' << "u" << '\t' << u_ << '\n' << "v" << '\t' << v_ << '\n' << "w" << '\t' << w_ << '\n' << "meanaff" << '\t' << meanaff_ << '\n' << "varaff" << '\t' << varaff_ << '\n' << "nbDSB" << '\t' << nbDSB_ << '\n' << "nbGenerations" << '\t' << nbGenerations_ << '\n' << "ismigration" << '\t' << ismigration_ << '\n' << "zygosity" << '\t' << zygosity_ << '\n' << "withDSB" << '\t' << withDSB_ << '\n' << "everygen" << '\t' << everygen_ << '\n' << "m" << '\t' << m_ << '\n' << "nbgenmig" << '\t' << nbgenmig_ << '\n' << "nbcore" << '\t' << nbcore_ << '\n' << "isallele" << '\t' << isallele_ << '\n' << "issampling" << '\t' << issampling_ << '\n' << "isanalytic" << '\t' << isanalytic_ << '\n' << "targetcomp" << '\t' << targetcomp_ << '\n' << "quantilenb" << '\t' << quantilenb_ << '\n' << "nbmeiperind" << '\t' << nbmeiperind_ << '\n' << "ctot" << '\t' << ctot_ << '\n' << "qhet_0" << '\t' << qhet_0 << '\n' << "qhom_0" << '\t' << qhom_0 << '\n' << "cfree_ctot_het_0" << '\t' << cfree_ctot_het_0 << '\n' << "cfree_ctot_hom_0" << '\t' << cfree_ctot_hom_0 << '\n' << "sigma_0" << '\t' << sigma_0_c << '\n';
    paramsfile.flush();
    	ofstream timefile ((name_+".time").c_str());
    	timefile << "Generation_number" << '\t' << "Mutations" << '\t' << "Fillnewpop" << '\t' << "afterfillpop" << '\t' << "calculstat" << '\t' << "Totaltime" <<'\n';
    	timefile.flush();
	ofstream generalfile1;
	ofstream allelefile1;
	ofstream generalfile2;
	ofstream allelefile2;
	vector<vector<vector<vector<int>>>*> vectpop = vector<vector<vector<vector<int>>>*> {&populations_};
	vector<vector<vector<int>>*> vectgen = vector<vector<vector<int>>*> {&genotypes_};
	vector<map<int,double>*> vectage = vector<map<int,double>*> {&Ageallele_};
	vector<map<int,vector<double>>*> vectinfo = vector<map<int,vector<double>>*> {&infoperallele_};
	vector<map<int,vector<double>>*> vectinfo_hom = vector<map<int,vector<double>>*> {&infoperallele_hom_};
	vector<map<int,vector<double>>*> vectinfo_het = vector<map<int,vector<double>>*> {&infoperallele_het_};
	double q1_=0;
	double q2_=0;
	vector<double*> vectq = vector<double*> {&q_};
	vector<double*> vectqsym = vector<double*> {&qsym_};//////////////////////////////
	vector<double*> vectqnum = vector<double*> {&qnum_};
	vector<double*> vectqdenom = vector<double*> {&qdenom_};
	vector<vector<vector<int>>*> vectfailed = vector<vector<vector<int>>*> {&nbfailedmeiosis_};
	int lenvect = vectpop.size();
	vector<int> one_to_twoN=vector<int>(2*N_,0);
	for(int i=0; i<2*N_; i++){
		one_to_twoN[i]=i;
	}
	for(int indgeneration=0; indgeneration<nbGenerations_; indgeneration++){
		//Initialisation time
		double temps_1, temps_printfile, temps_fillnewpop, temps_afterfillpop, temps_calculstat, temps_calcunit, temps_calcloop, temps_mutations_1, temps_mutations_site, temps_mutations_all, temps_update, temps_settozero;
		double t_nball, t_freqall, t_actall, t_ageall, t_qfertall, t_affall;
		clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22;//ajout t15 a t21 pour temps calcul stat par allele
		t1=clock();
		cout<<"generation "<<indgeneration<<endl;
		t9=clock();
		//set to 0 all the vectors containing the informatons from an allele during the meiosis
		for(int i_vect=0; i_vect<lenvect; i_vect++){
			t14=clock();
			for (auto const &it : *vectinfo[i_vect]){
				//(*vectinfo[i_vect])[it.first]=vector<double>{0,0,0,0,0,0,get_mean_affinity(it.first,vectpop[i_vect])}; ////////////////////////prends beaucoup de temps
				(*vectinfo[i_vect])[it.first]=vector<double>{0,0,0,0,0,0,0,0,0,0,0};
				(*vectinfo_hom[i_vect])[it.first]=vector<double>{0,0,0,0,0,0,0,0,0,0,0,(*vectinfo_hom[i_vect])[it.first][11],0,0};
				(*vectinfo_het[i_vect])[it.first]=vector<double>{0,0,0,0,0,0,0,0,0,0,0,(*vectinfo_hom[i_vect])[it.first][11],0,0};
			}
			t11=clock();
			temps_settozero=(float)(t11-t14)/CLOCKS_PER_SEC;
			//printf("temps_settozero = %f\n", temps_settozero);
			sitemutation(vectpop[i_vect], vectgen[i_vect]);
			t12=clock();
			temps_mutations_site=(float)(t12-t11)/CLOCKS_PER_SEC;
			//printf("temps_mutations_site = %f\n", temps_mutations_site);
			allelemutation(vectpop[i_vect], vectgen[i_vect], vectage[i_vect], one_to_twoN);
			t13=clock();
			temps_mutations_all=(float)(t13-t12)/CLOCKS_PER_SEC;
			//printf("temps_mutations_all = %f\n", temps_mutations_all);
		}
		updatemissingallele();
		t10=clock();
		temps_update=(float)(t10-t13)/CLOCKS_PER_SEC;
		//printf("temps_update = %f\n", temps_update);
		temps_mutations_1=(float)(t10-t9)/CLOCKS_PER_SEC;
		//printf("temps_mutations_1 = %f\n", temps_mutations_1);
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
			/*
			infoperallele1_hom_=infoperallele_;
			infoperallele2_hom_=infoperallele_;
			infoperallele1_het_=infoperallele_;
			infoperallele2_het_=infoperallele_;
			*/
			nbfailedmeiosis1_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
			nbfailedmeiosis2_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
			vectpop = vector<vector<vector<vector<int>>>*> {&populations1_,&populations2_};
			vectgen = vector<vector<vector<int>>*> {&genotypes1_,&genotypes2_};
			vectage = vector<map<int,double>*> {&Ageallele1_,&Ageallele2_};
			vectinfo = vector<map<int,vector<double>>*> {&infoperallele1_,&infoperallele2_};
			/*
			vectinfo_hom = vector<map<int,vector<double>>*> {&infoperallele1_hom_,&infoperallele2_hom_};
			vectinfo_het = vector<map<int,vector<double>>*> {&infoperallele1_het_,&infoperallele2_het_};
			*/
			vectq = vector<double*> {&q1_,&q2_};
			vectfailed = vector<vector<vector<int>>*> {&nbfailedmeiosis1_,&nbfailedmeiosis2_};
			lenvect = vectpop.size();
			//Initialise file generalfile et allelefile
			/////////////////////////////////////////////////add analytical q and fertility
			generalfile1.open ((name_+"_1.trace").c_str());
			generalfile1 << "Generation_number" << '\t' << "Number_of_allele" << '\t' << "Total_nb_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' << "Mean_Age" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" << '\t' << "q_intra" << '\t' << "q_inter" << '\t' << "fertility_intra" << '\t' << "fertility_inter" << '\t' << "FST_neutral" << '\t' << "FST_PRDM9" << '\t' << "Fert_newall" << '\t' << "Mean_c_occup_hom" << '\t' << "Mean_c_occup_het" << '\t' << "Mean_c_free_hom" << '\t' << "Mean_c_free_het" << '\n';
    			generalfile1.flush();
    			allelefile1.open ((name_+"_1.allele").c_str());
			allelefile1 << "Generation_number" << '\t' << "Allele_number" << '\t' << "Frequency" << '\t'  << "Activity" << '\t' << "Age" << '\t' << "q_allele" << '\t' << "Fertility_allele" << '\t' << "mean_affinity" << '\n';
    			allelefile1.flush();
			generalfile2.open ((name_+"_2.trace").c_str());
			generalfile2 << "Generation_number" << '\t' << "Number_of_allele" << '\t' << "Total_nb_allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' << "Mean_Age" << '\t' <<"Time" << '\t' << "Fertility_rate" << '\t' << "2_DSB_on_one_site_rate" << '\t' << "No_DSB_rate" << '\t' << "No_symmetrical_sites_rate" << '\t' << "q" << '\t' << "q_intra" << '\t' << "q_inter" << '\t' << "fertility_intra" << '\t' << "fertility_inter" << '\t' << "FST_neutral" << '\t' << "FST_PRDM9" << '\t' << "Fert_newall" << '\t' << "Mean_c_occup_hom" << '\t' << "Mean_c_occup_het" << '\n';
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
			fillnewpop(indgeneration, vectpop[i_vect], vectgen[i_vect], vectinfo[i_vect],vectinfo_hom[i_vect],vectinfo_het[i_vect], vectfailed[i_vect], vectq[i_vect], vectqsym[i_vect], vectqnum[i_vect], vectqdenom[i_vect]);
		}
		
		
		////////////////////////
		t5=clock();
		temps_fillnewpop=(float)(t5-t4)/CLOCKS_PER_SEC;
		printf("temps_fillnewpop = %f\n", temps_fillnewpop);
		////////////////////////
		parityIndex_=(parityIndex_+1)%2;
		updatemissingallele();
		////////////////////////////////////////////////////////////////////////////
		for(int i_vect=0; i_vect<lenvect; i_vect++){
			for (auto &it : *vectinfo_hom[i_vect]){
				if(it.second[10]!=0){
					it.second[11]=double(it.second[12])/it.second[10];
					it.second[13]=double(it.second[13])/it.second[10];
				}
			}
			for (auto &it : *vectinfo_het[i_vect]){
				if(it.second[10]!=0){
					it.second[11]=double(it.second[12])/it.second[10];
					it.second[13]=double(it.second[13])/it.second[10];
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////
		
		t22=clock();
		temps_update=(float)(t22-t5)/CLOCKS_PER_SEC;
		//printf("temps_update = %f\n", temps_update);
		for(int i_vect=0; i_vect<lenvect; i_vect++){
			for (auto &it : *vectage[i_vect]){
				it.second=it.second+freqall(it.first, vectgen[i_vect], vectpop[i_vect])*(N_*v_*nbDSB_)/(2*nbsite_);
			}
			for (auto &it : *vectinfo[i_vect]){
				if(it.second[5]!=0){
					it.second[4]=it.second[4]/it.second[5];
				}
			}
			
			for (auto &it : *vectinfo_hom[i_vect]){
				if(it.second[5]!=0){
					it.second[4]=it.second[4]/it.second[5];
				}
			}
			for (auto &it : *vectinfo_het[i_vect]){
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
		vector<double> current_qsym=vector<double>(2,0);/////////////////////////////////
		vector<double> mean_age=vector<double>(2,0);
		vector<double> mean_sigma=vector<double>(2,0);
		vector<vector<double>> mean_c_occup=vector<vector<double>>(2,vector<double>(4,0));
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
		vector<double> qfertall;
		vector<double> qfertall1;
		vector<double> qfertall2;
		vector<double> qallele;
		vector<double> qallele1;
		vector<double> qallele2;
		vector<double> fertilityallele;
		vector<double> fertilityallele1;
		vector<double> fertilityallele2;
		vector<double> meanaffinity;
		vector<double> meanaffinity1;
		vector<double> meanaffinity2;
		////////////////////////////////////////////////////////////////////////
		vector<double> relativenbmeiosisperall;
		vector<double> relativenbmeiosisperall1;
		vector<double> relativenbmeiosisperall2;
		vector<double> nblinkedsiteall;
		vector<double> nblinkedsiteall1;
		vector<double> nblinkedsiteall2;
		vector<double> nblinkedposall;
		vector<double> nblinkedposall1;
		vector<double> nblinkedposall2;
		vector<double> absolutenbmeiosisperall;
		vector<double> absolutenbmeiosisperall1;
		vector<double> absolutenbmeiosisperall2;
		vector<double> cfreectothom;
		vector<double> cfreectothom1;
		vector<double> cfreectothom2;
		vector<double> cfreectothet;
		vector<double> cfreectothet1;
		vector<double> cfreectothet2;
		vector<double> coccuphet;
		vector<double> coccuphet1;
		vector<double> coccuphet2;
		vector<double> coccuphom;
		vector<double> coccuphom1;
		vector<double> coccuphom2;
		vector<double> qfertallhom;
		vector<double> qfertallhom1;
		vector<double> qfertallhom2;
		vector<double> qfertallhet;
		vector<double> qfertallhet1;
		vector<double> qfertallhet2;
		vector<double> qallelehom;
		vector<double> qallelehom1;
		vector<double> qallelehom2;
		vector<double> qallelehet;
		vector<double> qallelehet1;
		vector<double> qallelehet2;
		vector<double> fertilityallelehom;
		vector<double> fertilityallelehom1;
		vector<double> fertilityallelehom2;
		vector<double> fertilityallelehet;
		vector<double> fertilityallelehet1;
		vector<double> fertilityallelehet2;
		////////////////////////////////////////////////////////////////////////
		//map<int,vector<double>> analytic_allele;
		vector<map<int,vector<double>>> analytic_indiv;///////////////////////////////////////////////////////////////
		double mean_fertnewall;
		////////////////////////
		t6=clock();
		temps_afterfillpop=(float)(t6-t5)/CLOCKS_PER_SEC;
		printf("temps_afterfillpop = %f\n", temps_afterfillpop);
		////////////////////////
		if (indgeneration % everygen_ ==0)  {
			if(ismigration_==false){
				// generalfile
				if(issampling_==true){
					allele_nb_trace[0]=get_allele_number(vectgen,true)[0];
					current_div[0] = get_current_diversity(&genotypes_);
					current_act[0] = get_current_activity(&genotypes_, &populations_);
					Fertility_rate[0] = 1-(double(nbfailedmeiosis_[indgeneration][3])/(2*N_+nbfailedmeiosis_[indgeneration][3]));
					twoDSB[0] = double(nbfailedmeiosis_[indgeneration][0])/(2*N_+nbfailedmeiosis_[indgeneration][3]);
					noDSB[0] = double(nbfailedmeiosis_[indgeneration][1])/(2*N_+nbfailedmeiosis_[indgeneration][3]);
					nosym[0] = double(nbfailedmeiosis_[indgeneration][2])/(2*N_+nbfailedmeiosis_[indgeneration][3]);
					current_q[0] = q_;
					current_qsym[0] = qsym_;
					mean_age[0] = get_mean_age(vectgen[0], vectage[0]);
					mean_c_occup[0] = get_mean_coccup_cfree(vectgen[0], vectinfo_hom[0], vectinfo_het[0]);// to add in migration
					q_fertility = vector<double>{0,0};
				}else{
					allele_nb_trace[0]= 0;
					current_div[0] = 0;
					current_act[0] = 0;
					Fertility_rate[0] = 0;
					twoDSB[0] = 0;
					noDSB[0] = 0;
					nosym[0] = 0;
					current_q[0] = 0;
					current_qsym[0] = 0;
					mean_age[0] = 0;
					mean_c_occup[0] = vector<double>{0,0,0,0};
					q_fertility = vector<double>{0,0};
				}
				//if(isanalytic_==true){}else{}
				analytic_indiv=q_fert_individual_analytique(&genotypes_, &populations_, &infoperallele_hom_, &infoperallele_het_);
				mean_sigma[0]=get_mean_sigma(&genotypes_, analytic_indiv);
				//mean_fertnewall = Mean_fert_new_allele(&genotypes_,&populations_);
				mean_fertnewall = log(Mean_fert_new_allele(&genotypes_,&populations_, &infoperallele_het_))-log(analytic_indiv[0][-1][1]);////////////
				//q_fertility = get_q_fertility_indep(vectpop, vectgen, nbloop_, 1);//utile pour migration
				////////////////////////
				t7=clock();
				temps_calcunit=(float)(t7-t6)/CLOCKS_PER_SEC;
				printf("temps_calcunit = %f\n", temps_calcunit);
				////////////////////////
				if(isallele_==true){
					for (auto const &it : Siteforeacheallele_){//////////////////////////////////////////////////////////////// =========> openMP peut marcher si on enleve auto (for cassique) et push_back (declaration des vecteurs de longueur la taille de Siteforeachallele)
    					if(it.first!=-2){
    						t15=clock();
							allele_nb.push_back(it.first);
							t16=clock();
							t_nball=t_nball+(float)(t16-t15)/CLOCKS_PER_SEC;
							freqallele.push_back(freqall(it.first, &genotypes_, &populations_));
							t17=clock();
							t_freqall=t_freqall+(float)(t17-t16)/CLOCKS_PER_SEC;
							activityall.push_back(actall(it.first, &populations_));
							t18=clock();
							t_actall=t_actall+(float)(t18-t17)/CLOCKS_PER_SEC;
							ageallele.push_back(get_age_allele(it.first, &Ageallele_));
							t19=clock();
							t_ageall=t_ageall+(float)(t19-t18)/CLOCKS_PER_SEC;
							qfertall=get_info_allele(it.first, &infoperallele_);
							qallele.push_back(qfertall[0]);
							fertilityallele.push_back(qfertall[1]);
							///////////////////////////
							qfertallhom=get_info_allele(it.first, &infoperallele_hom_);
							qallelehom.push_back(qfertallhom[0]);
							fertilityallelehom.push_back(qfertallhom[1]);
							qfertallhet=get_info_allele(it.first, &infoperallele_het_);
							qallelehet.push_back(qfertallhet[0]);
							fertilityallelehet.push_back(qfertallhet[1]);
							///////////////////////////
							t20=clock();
							t_qfertall=t_qfertall+(float)(t20-t19)/CLOCKS_PER_SEC;
							if(it.first==-3){
								meanaffinity.push_back(get_mean_affinity(it.first,&populations_));
								relativenbmeiosisperall.push_back(0); //to add in migration
								absolutenbmeiosisperall.push_back(0); //to add in migration
								nblinkedsiteall.push_back(0); //to add in migration
								nblinkedposall.push_back(0); //to add in migration
								cfreectothom.push_back(0); //to add in migration
								cfreectothet.push_back(0); //to add in migration
							}else{
								//meanaffinity.push_back(infoperallele_[it.first][6]);
								meanaffinity.push_back(get_mean_affinity(it.first,&populations_));
								relativenbmeiosisperall.push_back(infoperallele_[it.first][7]); //to add in migration
								absolutenbmeiosisperall.push_back(infoperallele_[it.first][10]); //to add in migration
								nblinkedsiteall.push_back(infoperallele_[it.first][8]/infoperallele_[it.first][10]); //to add in migration
								nblinkedposall.push_back(infoperallele_[it.first][9]/infoperallele_[it.first][10]); //to add in migration
								cfreectothom.push_back(infoperallele_hom_[it.first][11]/(2*ctot_)); //to add in migration
								cfreectothet.push_back(infoperallele_het_[it.first][11]/ctot_); //to add in migration
								coccuphet.push_back(infoperallele_het_[it.first][13]); //to add in migration
								coccuphom.push_back(infoperallele_hom_[it.first][13]); //to add in migration
							}
							t21=clock();
							t_affall=t_affall+(float)(t21-t20)/CLOCKS_PER_SEC;
        				}
        			}
    				t_nball=t_nball/double((Siteforeacheallele_.size())-1);
    				//printf("t_nball = %f\n", t_nball);
    				t_freqall=t_freqall/double((Siteforeacheallele_.size())-1);
    				//printf("t_freqall = %f\n", t_freqall);
    				t_actall=t_actall/double((Siteforeacheallele_.size())-1);
    				//printf("t_actall = %f\n", t_actall);
    				t_ageall=t_ageall/double((Siteforeacheallele_.size())-1);
    				//printf("t_ageall = %f\n", t_ageall);
    				t_qfertall=t_qfertall/double((Siteforeacheallele_.size())-1);
    				//printf("t_qfertall = %f\n", t_qfertall);
    				t_affall=t_affall/double((Siteforeacheallele_.size())-1);
    				//printf("t_affall = %f\n", t_affall);
        				
				}else{
					allele_nb.push_back(0);
					freqallele.push_back(0);
					activityall.push_back(0);
					ageallele.push_back(0);
					qfertall.push_back(0);
					qallele.push_back(0);
					fertilityallele.push_back(0);
					meanaffinity.push_back(0);
				}
				////////////////////////
				t8=clock();
				temps_calcloop=(float)(t8-t7)/CLOCKS_PER_SEC;
				printf("temps_calcloop = %f\n", temps_calcloop);
				////////////////////////
			}else if(ismigration_==true){////////////////////////////////to add: analytical calcul, isallele, isanalytic, issampling 
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
						qfertall1=get_info_allele(it.first, &infoperallele1_);
						qfertall2=get_info_allele(it.first, &infoperallele2_);
						qallele1.push_back(qfertall1[0]);
						qallele2.push_back(qfertall2[0]);
						fertilityallele1.push_back(qfertall1[1]);
						fertilityallele2.push_back(qfertall2[1]);
						if(it.first==-3){
							meanaffinity1.push_back(get_mean_affinity(it.first,&populations1_));
							meanaffinity2.push_back(get_mean_affinity(it.first,&populations2_));
						}else{
							//meanaffinity1.push_back(infoperallele1_[it.first][6]);//////////////// changer avec get_mean_affinity
							meanaffinity1.push_back(get_mean_affinity(it.first,&populations1_));
							//meanaffinity2.push_back(infoperallele2_[it.first][6]);//////////////// changer avec get_mean_affinity
							meanaffinity2.push_back(get_mean_affinity(it.first,&populations2_));
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
				generalfile << indgeneration << '\t' << allele_nb_trace[0] << '\t' << current_div[0] << '\t'  << current_act[0] << '\t' << mean_age[0] << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << Fertility_rate[0] << '\t' << twoDSB[0] << '\t'<< noDSB[0] << '\t' << nosym[0] << '\t' << current_q[0] << '\t' << q_fertility[0] << '\t' << q_fertility[1] << '\t' << analytic_indiv[0][-1][0] <<  '\t' << analytic_indiv[0][-1][1] << '\t' << current_qsym[0] << '\t' << mean_fertnewall << '\t' << mean_c_occup[0][0] << '\t' << mean_c_occup[0][1] << '\t' << mean_c_occup[0][2] << '\t' << mean_c_occup[0][3] << '\t' << mean_sigma[0] << '\n';
            	generalfile.flush();
            	timefile << indgeneration << '\t' << temps_mutations_1 << '\t' << temps_fillnewpop << '\t' << temps_afterfillpop << '\t' << temps_calculstat << '\t' << temps_1 << '\n';
            	timefile.flush();
            	for (int indall = 0; indall<allele_nb.size(); indall++){
            		allelefile << indgeneration << '\t' << allele_nb[indall] << '\t' << freqallele[indall] << '\t'  << activityall[indall] << '\t' << ageallele[indall] << '\t' << qallele[indall] << '\t' << fertilityallele[indall] << '\t' << meanaffinity[indall] << '\t' <<analytic_indiv[0][int(allele_nb[indall])][0] <<  '\t' << analytic_indiv[0][int(allele_nb[indall])][1] << '\t' << relativenbmeiosisperall[indall] << '\t' << nblinkedsiteall[indall] << '\t' << nblinkedposall[indall] << '\t' << absolutenbmeiosisperall[indall] << '\t' << cfreectothom[indall] << '\t' << cfreectothet[indall] << '\t' << qallelehom[indall] << '\t' << qallelehet[indall] << '\t' << fertilityallelehom[indall] << '\t' << fertilityallelehet[indall] << '\t' << analytic_indiv[1][int(allele_nb[indall])][0] << '\t' << analytic_indiv[2][int(allele_nb[indall])][0] << '\t' << analytic_indiv[1][int(allele_nb[indall])][1] << '\t' << analytic_indiv[2][int(allele_nb[indall])][1] << '\t' << analytic_indiv[3][int(allele_nb[indall])][0] << '\t' << analytic_indiv[3][int(allele_nb[indall])][1] << '\t' << analytic_indiv[0][int(allele_nb[indall])][2] <<'\n';//
            		allelefile.flush();
            	}
			}else if(ismigration_==true){
				generalfile1 << indgeneration << '\t' << allele_nb_trace[0] << '\t' << tot_allele_nb[0] << '\t' << current_div[0] << '\t'  << current_act[0] << '\t' << mean_age[0] << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << Fertility_rate[0] << '\t' << twoDSB[0] << '\t'<< noDSB[0] << '\t' << nosym[0] << '\t' << current_q[0] << '\t' << q_fertility_1[0] << '\t' << q_fertility_hybrid[0] << '\t' << q_fertility_1[1] << '\t' << q_fertility_hybrid[1] << '\t' << FST_neutral << '\t' << FST_PRDM9 << '\t' << mean_fertnewall << '\n';
            	generalfile1.flush();
            	for (int indall = 0; indall<allele_nb.size(); indall++){
            		allelefile1 << indgeneration << '\t' << allele_nb[indall] << '\t' << freqallele1[indall] << '\t'  << activityall1[indall] << '\t' << ageallele1[indall] << '\t' << qallele1[indall] << '\t' << fertilityallele1[indall] << '\t' << meanaffinity1[indall] << '\n';
            		allelefile1.flush();
            	}
            	generalfile2 << indgeneration << '\t' << allele_nb_trace[1] << '\t' << tot_allele_nb[0] << '\t' << current_div[1] << '\t'  << current_act[1] << '\t' << mean_age[1] << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << Fertility_rate[1] << '\t' << twoDSB[0] << '\t'<< noDSB[1] << '\t' << nosym[1] << '\t' << current_q[1] << '\t' << q_fertility_2[0] << '\t' << q_fertility_hybrid[0] << '\t' << q_fertility_2[1] << '\t' << q_fertility_hybrid[1] << '\t' << FST_neutral << '\t' << FST_PRDM9 << '\t' << mean_fertnewall << '\n';
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
	timefile.close();
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
	}else{
		vector<int> allnb=vector<int>(3,0);
		for(int i=0; i<2; i++){
			for(auto &it : Siteforeacheallele_){
				if(it.first!=-3 and it.first!=-2){
					vector<int>::iterator itvect = find((*vectgen[i])[parityIndex_].begin(),(*vectgen[i])[parityIndex_].end(),it.first);
					if(itvect != (*vectgen[i])[parityIndex_].end()){
					//vector<int>::iterator itvect = find((*vectgen[i])[(parityIndex_+1)%2].begin(),(*vectgen[i])[(parityIndex_+1)%2].end(),it.first);
					//if(itvect != (*vectgen[i])[(parityIndex_+1)%2].end()){
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
	return (meanage)/((*genotype)[parityIndex_].size());
}

//get the mean sigma of the allele segregating in the population
double Model::get_mean_sigma(vector<vector<int>>* genotype, vector<map<int,vector<double>>> analytic_indiv){
	double meansigma=0;
	for(auto const allele : (*genotype)[parityIndex_]){
		meansigma+=analytic_indiv[0][allele][2];
	}
	return (meansigma)/((*genotype)[parityIndex_].size());
}

////////////////////////////////////////////////////////////////////////////////////////////
//get the mean number of linked protein of the allele segregating in the population
vector<double> Model::get_mean_coccup_cfree(vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het){
	double mean_coccup_hom=0;
	double mean_coccup_het=0;
	double mean_cfree_hom=0;
	double mean_cfree_het=0;
	int denom_coccup_het=0;
	int denom_coccup_hom=0;
	int denom_cfree_het=0;
	int denom_cfree_hom=0;
	for(auto const allele : (*genotype)[parityIndex_]){
		if((*infoperallele_hom)[allele][13]>0){
			mean_coccup_hom+=(*infoperallele_hom)[allele][13];
			denom_coccup_hom+=1;
		}
		if((*infoperallele_het)[allele][13]>0){
			mean_coccup_het+=(*infoperallele_het)[allele][13];
			denom_coccup_het+=1;
		}
		if((*infoperallele_hom)[allele][11]>0){
			mean_cfree_hom+=(*infoperallele_hom)[allele][11];
			denom_cfree_hom+=1;
		}
		if((*infoperallele_het)[allele][11]>0){
			mean_cfree_het+=(*infoperallele_het)[allele][11];
			denom_cfree_het+=1;
		}
	}
	return {double(mean_coccup_hom)/double(denom_coccup_hom), double(mean_coccup_het)/double(denom_coccup_het), double(mean_cfree_hom)/double(denom_cfree_hom), double(mean_cfree_het)/double(denom_cfree_het)};
}
////////////////////////////////////////////////////////////////////////////////////////////


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
		else{
			return vector<double>{0};
		}
	}
}

//give the frequence of each allele in the population
double Model::freqallele(int allelename, vector<vector<int>>* genotype){
		return double (count((*genotype)[parityIndex_].begin(), (*genotype)[parityIndex_].end(), allelename))/((*genotype)[parityIndex_].size());
		//return double (count((*genotype)[(parityIndex_+1)%2].begin(), (*genotype)[(parityIndex_+1)%2].end(), allelename))/((*genotype)[(parityIndex_+1)%2].size());
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
			//if((*population)[(parityIndex_+1)%2][i][all]==1){
				moyactsite+=1;			
			}
		}
		moyact+=moyactsite/(2*N_);
	}
	moyact=moyact/nbsite_;
	return (moyact);
}

//freq, activity of neutral sites
vector<double> Model::freqneutral(vector<vector<vector<int>>>* population){
	double moyfreq=0;
	double moy2f=0;
	for(auto all : Siteforeacheallele_[-3]){
		double moyactsite=0;
		for(int i=0; i<2*N_; i++){
			if((*population)[parityIndex_][i][all]==1){
			//if((*population)[(parityIndex_+1)%2][i][all]==1){
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
	int lengthsite=Siteforeacheallele_[allele].size();
	vector<int> nbactsiteperind = vector<int>(lengthsite,0);
	for(auto const site : Siteforeacheallele_[allele]){//auto replaced bu basi for
	//#pragma omp parallel for num_threads(nbcore_)// ne marche pas car les boucles ne sont pas indep => pour nbactivesite je peux creer un tableau de longueur lengthsite et pour chaque boucle inserer la some au bon indice puis en sortant de la boucle faire la somme des elements de tableau
	//int lengthsite=Siteforeacheallele_[allele].size();
	//for(int s=0; s<lengthsite; s++){
		//int site=Siteforeacheallele_[allele][s];
		int activitysite=0;
		for(int siteind=0; siteind<2*N_; siteind++){
			if((*pop)[parityIndex_][siteind][site]==1){
			//if((*pop)[(parityIndex_+1)%2][siteind][site]==1){
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
				//if((*vectpop[j])[(parityIndex_+1)%2][i][all]==1){
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






// calculation of q and fertility indepedently from the system (sampling) for the next generation !!!!
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
double Model::q_two_hap(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom){ //return q///////////////////////////////////////////////// CHANGER IND_GENE (prise en compte concentration)
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
vector<int> Model::get_one_gamete(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom){ //return un vector de vector///////////////////////////////////////////////////////CHANGER IND_GENE
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





vector<map<int,vector<double>>> Model:: q_fert_individual_analytique(vector<vector<int>>* genotype, vector<vector<vector<int>>>* pop , map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het){ 
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

	
	map<int,vector<double>> res;
	map<int,vector<double>> qmoy;
	map<int,vector<double>> fertmoy;
	////////////////////////////////////////////////////////
	map<int,vector<double>> res_hom;
	map<int,vector<double>> qmoy_hom;
	map<int,vector<double>> fertmoy_hom;
	map<int,vector<double>> res_het;
	map<int,vector<double>> qmoy_het;
	map<int,vector<double>> fertmoy_het;
	map<int,vector<double>> qmoy_hemi;/////
	map<int,vector<double>> fertmoy_hemi;/////
	map<int,vector<double>> res_hemi;/////
	
	map<int,double> fertmoy_sig_het;
	map<int,double> fertmoy_sig_hom;

	////////////////////////////////////////////////////////
	double qindmoy=0;
	double fertindmoy=0;
	//initialisation of countainers for new allele when an allele is not in a heterozygous state in the population
	double cfree_newall=0;
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> newpos = choosemany(nbsite_, freepos);
	vector<double> infoallele_sig; 
	
	for (int indiv=0; indiv<2*N_; indiv=indiv+2){ // for each individual in the population
		bool homoz=true;
		vector<int>genotype_indiv={(*genotype)[parityIndex_][indiv]};
		if((*genotype)[parityIndex_][indiv]!=(*genotype)[parityIndex_][indiv+1]){
			genotype_indiv.push_back((*genotype)[parityIndex_][indiv+1]);
			homoz=false; // it is not homozygous
		}
		map<int,vector<double>> infoallele;
		map<int,vector<double>> infoallele_sig_hom;
		map<int,vector<double>> infoallele_sig_hemi;
	
		double qal=0;
		double qal_hemi1=0;
		double qal_hemi2=0;
		
		double qal_sig_hom_0=0;
		double qal_sig_hom_1=0;
		double qal_sig_het=0;
		
		
		int nbloop;
		if (indiv==0){
			nbloop=genotype_indiv.size()+1;
		}else{
			nbloop=genotype_indiv.size();
		}
		
		for(int i=0; i<nbloop; i++){
			for(int ind_hom_het=0; ind_hom_het<2; ind_hom_het++){ //2 loops : the first one is what happen normaly (if the individual is homozygous (heterozygous), the cfree and moyprobalinks will be calculated for a homozygous (heterozygous)) , the second loop will be for the contrary (if the indvidual is homozygous, all will be calculated as a heterozygous and vice et versa)
				double z;
				if(i!=genotype_indiv.size()){
					z=genotype_indiv[i];;
				}else{
					z=0;
				}
				double ind_gen=1;
				if((zygosity_ and homoz and ind_hom_het==0) or (zygosity_ and homoz==0 and ind_hom_het==1) ){ //////////////////////////////////////
					ind_gen=2;
				}
				////////////////////////////
				if(targetcomp_){
					vector<map<int,vector<double>>*> vectinfo_hom_het = vector<map<int,vector<double>>*> {infoperallele_hom, infoperallele_het};
					double index_vectinfo_hom_het=1;
					double is_diff=1;
					double ctot=ctot_;
					double cfree;
					double threshold=cfreethreshold_;
					if(indiv==0 and i==genotype_indiv.size()){ //only time we calculate cfree_newall
						cfree=1/float(1000)*ctot;
					}else{
						if ((homoz and ind_hom_het==0) or (homoz==0 and ind_hom_het==1)){
							if(zygosity_){
								ctot=2*ctot_;
							}
							index_vectinfo_hom_het=0;
						}
						if((*vectinfo_hom_het[index_vectinfo_hom_het])[z][11]==0){
							cfree=1/float(1000)*ctot;
						}else{
							cfree=(*vectinfo_hom_het[index_vectinfo_hom_het])[z][11];
						}
					}
					
					//Two cases possible :
					//1) If we want to take the affinity categories map (quantilenb_!=0),if the site is active, the category in which the affinity of this site enters gains one site
					for(auto quantile_index : nbsitesperquantile_){
						nbsitesperquantile_[quantile_index.first][1]=0;
					}
					if(quantilenb_!=0){
						if(indiv==0 and i==genotype_indiv.size()){////////////////////////////////////1
							for(auto const &it : newpos){
								bool is_quantile=0;
								for(auto quantile_index : nbsitesperquantile_){
									if(Affinity_[it]<=quantile_index.first and is_quantile==0){
										nbsitesperquantile_[quantile_index.first][1]+=4;//if this site is active, the category in which the affinity of this site enters gains one site
										is_quantile=1;
									}
								}
							}
						}else{///////////////////////////////////////////////////////1
							for(auto const &it : Siteforeacheallele_[z]){
								for(int j=0; j<4; j++){
									int chrom=indiv+j/2;
									if((*pop)[parityIndex_][chrom][it]==1){
										bool is_quantile=0;
										for(auto quantile_index : nbsitesperquantile_){
											if(Affinity_[it]<=quantile_index.first and is_quantile==0){
												nbsitesperquantile_[quantile_index.first][1]+=1;//if this site is active, the category in which the affinity of this site enters gains one site
												is_quantile=1;
											}
										}
									}
								}
							}
						}//////////////////////////////////////////////////////////1
					}
					//2) If we want don't want affinity categories we don't need previous setup
					double cfree_min=0;
					double cfree_max=ctot;
					double sum_p_occup=0;
					while(is_diff==1){
						is_diff=0;
						sum_p_occup=0;//the mean sum of occupied site
						//---------------
						//sum_p_occup
						//---------------
						//In the case of affinity categories
						if(quantilenb_!=0){
							for(auto quantile_index : nbsitesperquantile_){
								sum_p_occup+=double(nbsitesperquantile_[quantile_index.first][1]*(cfree*nbsitesperquantile_[quantile_index.first][0])/(1+cfree*nbsitesperquantile_[quantile_index.first][0]));//The probability of binding in each category multiplied by the number of sites in this category is added to the mean sum of occupied site
							}
						}
						//In the other case
						else{
							if(indiv==0 and i==genotype_indiv.size()){////////////////////////////////////1
								for(auto const &it : newpos){//for each position
									double p_occup=cfree*Affinity_[it]/(1+cfree*Affinity_[it]);//the probability of binding is calculated
									sum_p_occup+=4*p_occup; //coeff 4 because no erosion in new alleles
								}
							}else{///////////////////////////////////////////////////////////1
								for(auto const &it : Siteforeacheallele_[z]){//for each position
									double p_occup=cfree*Affinity_[it]/(1+cfree*Affinity_[it]);//the probability of binding is calculated
									for(int j=0; j<4; j++){//for each site of this position
										int chrom=indiv+j/2;
										if((*pop)[parityIndex_][chrom][it]==1){//if the site is active, its probability of binding is added to the mean sum of occupied site
											sum_p_occup+=p_occup;
										}
									}
								}
							}//////////////////////////////////////////////////////////1
						}
						if(cfree+sum_p_occup>ctot){
							//if cfree is to high
							is_diff=1;
							cfree_max=cfree;
							cfree=double(cfree_max+cfree_min)/2;
						}else{
							//if cfree is to low or the right one
							double new_cfree=ctot-sum_p_occup;
							if(abs(double(cfree)/double(ctot)-double(new_cfree)/double(ctot))<=threshold){
								//the difference between the current and the next cfree/ctot is below or equal to the thershold
								is_diff=0;
							}else{
								//the difference between the current and the next cfree/ctot is above the thershold
								is_diff=1;
								cfree_min=cfree;
								cfree=double(cfree_max+cfree_min)/2;
							}
						}
					}
					if(indiv==0 and i==genotype_indiv.size()){
						cfree_newall=cfree;
					}else{
						ind_gen=cfree;
					}
				}
				////////////////////////////
				if(i<genotype_indiv.size()){
					double moyprobalinkab=0;
					double moyprobalinka2b=0;
					double moyprobalinkb2a=0;
					double moyprobalinka=0;
					double moyprobalinkb=0;
					for(auto const &it : Siteforeacheallele_[genotype_indiv[i]]){
						double xa=0;
						double xb=0;
						if((*pop)[parityIndex_][indiv][it]==1){
							xa=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
						}else{
							xa=0;
						}
						if((*pop)[parityIndex_][indiv+1][it]==1){
							xb=ind_gen*Affinity_[it]/(1+ind_gen*Affinity_[it]);
						}else{
							xb=0;
						}
						moyprobalinkab+=xa*xb;
						moyprobalinka2b+=puissance_double(2,xa)*xb;
						moyprobalinkb2a+=puissance_double(2,xb)*xa;
						moyprobalinka+=xa;
						moyprobalinkb+=xb;
					}
					if(ind_hom_het==0){
						infoallele[genotype_indiv[i]].push_back(moyprobalinkab/nbsite_);
						infoallele[genotype_indiv[i]].push_back(moyprobalinka2b/nbsite_);
						infoallele[genotype_indiv[i]].push_back(moyprobalinkb2a/nbsite_);
						infoallele[genotype_indiv[i]].push_back(moyprobalinka/nbsite_);
						infoallele[genotype_indiv[i]].push_back(moyprobalinkb/nbsite_);
					}
					
					if((homoz and ind_hom_het==0) or (homoz==0 and ind_hom_het==1)){
						infoallele_sig_hom[genotype_indiv[i]].push_back(moyprobalinkab/nbsite_);
						infoallele_sig_hom[genotype_indiv[i]].push_back(moyprobalinka2b/nbsite_);
						infoallele_sig_hom[genotype_indiv[i]].push_back(moyprobalinkb2a/nbsite_);
						infoallele_sig_hom[genotype_indiv[i]].push_back(moyprobalinka/nbsite_);
						infoallele_sig_hom[genotype_indiv[i]].push_back(moyprobalinkb/nbsite_);
					}else{
						infoallele_sig_hemi[genotype_indiv[i]].push_back(moyprobalinkab/nbsite_);
						infoallele_sig_hemi[genotype_indiv[i]].push_back(moyprobalinka2b/nbsite_);
						infoallele_sig_hemi[genotype_indiv[i]].push_back(moyprobalinkb2a/nbsite_);
						infoallele_sig_hemi[genotype_indiv[i]].push_back(moyprobalinka/nbsite_);
						infoallele_sig_hemi[genotype_indiv[i]].push_back(moyprobalinkb/nbsite_);
					}
					
				}

				if(indiv==0 and i==genotype_indiv.size()){
					double moyprobalinkab=0;
					double moyprobalinka2b=0;
					double moyprobalinkb2a=0;
					double moyprobalinka=0;
					double moyprobalinkb=0;
					for(auto const &it : newpos){
						double xa=0;
						double xb=0;
						xa=cfree_newall*Affinity_[it]/(1+cfree_newall*Affinity_[it]);
						xb=cfree_newall*Affinity_[it]/(1+cfree_newall*Affinity_[it]);
						moyprobalinkab+=xa*xb;
						moyprobalinka2b+=puissance_double(2,xa)*xb;
						moyprobalinkb2a+=puissance_double(2,xb)*xa;
						moyprobalinka+=xa;
						moyprobalinkb+=xb;
					}
					infoallele_sig.push_back(moyprobalinkab/nbsite_);
					infoallele_sig.push_back(moyprobalinka2b/nbsite_);
					infoallele_sig.push_back(moyprobalinkb2a/nbsite_);
					infoallele_sig.push_back(moyprobalinka/nbsite_);
					infoallele_sig.push_back(moyprobalinkb/nbsite_);
				}
			}
		}
		
		if(homoz){
			qal=(4*infoallele[genotype_indiv[0]][0]-infoallele[genotype_indiv[0]][1]-infoallele[genotype_indiv[0]][2])/(infoallele[genotype_indiv[0]][3]+infoallele[genotype_indiv[0]][4]);
			
			//qal calculate in homozygous and heterozygous context for each allele to compute the sigma(z)
			qal_sig_hom_0=qal;
			//qal_sig_het=(4*infoallele[genotype_indiv[0]][0]-infoallele[genotype_indiv[0]][1]-infoallele[genotype_indiv[0]][2]+4*infoallele_sig[0]-infoallele_sig[1]-infoallele_sig[2])/(infoallele[genotype_indiv[0]][3]+infoallele[genotype_indiv[0]][4]+infoallele_sig[3]+infoallele_sig[4]);
			qal_sig_het=(4*infoallele_sig_hemi[genotype_indiv[0]][0]-infoallele_sig_hemi[genotype_indiv[0]][1]-infoallele_sig_hemi[genotype_indiv[0]][2])/(infoallele_sig_hemi[genotype_indiv[0]][3]+infoallele_sig_hemi[genotype_indiv[0]][4]);
		}else{
			//qal here correspond to the q of the heterozygot and not the hemizygot
			qal=(4*infoallele[genotype_indiv[0]][0]-infoallele[genotype_indiv[0]][1]-infoallele[genotype_indiv[0]][2]+4*infoallele[genotype_indiv[1]][0]-infoallele[genotype_indiv[1]][1]-infoallele[genotype_indiv[1]][2])/(infoallele[genotype_indiv[0]][3]+infoallele[genotype_indiv[0]][4]+infoallele[genotype_indiv[1]][3]+infoallele[genotype_indiv[1]][4]);
			
			//qal calculate in homozygous and heterozygous context for each allele to compute the sigma(z)
			qal_sig_het=qal;
			qal_sig_hom_0=(4*infoallele_sig_hom[genotype_indiv[0]][0]-infoallele_sig_hom[genotype_indiv[0]][1]-infoallele_sig_hom[genotype_indiv[0]][2])/(infoallele_sig_hom[genotype_indiv[0]][3]+infoallele_sig_hom[genotype_indiv[0]][4]);
			qal_sig_hom_1=(4*infoallele_sig_hom[genotype_indiv[1]][0]-infoallele_sig_hom[genotype_indiv[1]][1]-infoallele_sig_hom[genotype_indiv[1]][2])/(infoallele_sig_hom[genotype_indiv[1]][3]+infoallele_sig_hom[genotype_indiv[1]][4]);

			
			//for hemizygots : homozygot but with same concentration as heterozygot
			qal_hemi1=(4*infoallele[genotype_indiv[0]][0]-infoallele[genotype_indiv[0]][1]-infoallele[genotype_indiv[0]][2])/(infoallele[genotype_indiv[0]][3]+infoallele[genotype_indiv[0]][4]);/////
			qal_hemi2=(4*infoallele[genotype_indiv[1]][0]-infoallele[genotype_indiv[1]][1]-infoallele[genotype_indiv[1]][2])/(infoallele[genotype_indiv[1]][3]+infoallele[genotype_indiv[1]][4]);/////
		}
		if(zygosity_){
			// if genetic dosage : homozygous count for 1 allele
			if(homoz){
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				qmoy_hom[genotype_indiv[0]].push_back(qal);
				fertmoy_hom[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				
				if(fertmoy_sig_het[genotype_indiv[0]]==0){
					fertmoy_sig_het[genotype_indiv[0]]=1-exp(-nbDSB_*qal_sig_het);
				}
				if(fertmoy_sig_hom[genotype_indiv[0]]==0){
					fertmoy_sig_hom[genotype_indiv[0]]=1-exp(-nbDSB_*qal);
				}
				
				
			}else{
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy[genotype_indiv[1]].push_back(qal);
				fertmoy[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				qmoy_het[genotype_indiv[0]].push_back(qal);
				fertmoy_het[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy_het[genotype_indiv[1]].push_back(qal);
				fertmoy_het[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				qmoy_hemi[genotype_indiv[0]].push_back(qal_hemi1);/////
				fertmoy_hemi[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal_hemi1));/////
				qmoy_hemi[genotype_indiv[1]].push_back(qal_hemi2);/////
				fertmoy_hemi[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal_hemi2));/////
				////////////////////////////////////////////////////////
				
				if(fertmoy_sig_het[genotype_indiv[0]]==0){
					fertmoy_sig_het[genotype_indiv[0]]=1-exp(-nbDSB_*qal);
				}
				if(fertmoy_sig_het[genotype_indiv[1]]==0){
					fertmoy_sig_het[genotype_indiv[1]]=1-exp(-nbDSB_*qal);
				}
				if(fertmoy_sig_hom[genotype_indiv[0]]==0){
					fertmoy_sig_hom[genotype_indiv[0]]=1-exp(-nbDSB_*qal_sig_hom_0);
				}
				if(fertmoy_sig_hom[genotype_indiv[1]]==0){
					fertmoy_sig_hom[genotype_indiv[1]]=1-exp(-nbDSB_*qal_sig_hom_1);
				}
				
			}
		}else{// if no genetic dosage : homozygous count for 2 allele
			if(homoz){ // if homozygous for the allele PRDM9 : the q and the fertility count for 2
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				qmoy_hom[genotype_indiv[0]].push_back(qal);
				fertmoy_hom[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy_hom[genotype_indiv[0]].push_back(qal);
				fertmoy_hom[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				
				if(fertmoy_sig_het[genotype_indiv[0]]==0){
					fertmoy_sig_het[genotype_indiv[0]]=1-exp(-nbDSB_*qal_sig_het);
				}
				if(fertmoy_sig_hom[genotype_indiv[0]]==0){
					fertmoy_sig_hom[genotype_indiv[0]]=1-exp(-nbDSB_*qal);
				}
				
			}else{
				qmoy[genotype_indiv[0]].push_back(qal);
				fertmoy[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy[genotype_indiv[1]].push_back(qal);
				fertmoy[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				qmoy_het[genotype_indiv[0]].push_back(qal);
				fertmoy_het[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal));
				qmoy_het[genotype_indiv[1]].push_back(qal);
				fertmoy_het[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal));
				////////////////////////////////////////////////////////
				qmoy_hemi[genotype_indiv[0]].push_back(qal_hemi1);/////
				fertmoy_hemi[genotype_indiv[0]].push_back(1-exp(-nbDSB_*qal_hemi1));/////
				qmoy_hemi[genotype_indiv[1]].push_back(qal_hemi2);/////
				fertmoy_hemi[genotype_indiv[1]].push_back(1-exp(-nbDSB_*qal_hemi2));/////
				////////////////////////////////////////////////////////
				
				if(fertmoy_sig_het[genotype_indiv[0]]==0){
					fertmoy_sig_het[genotype_indiv[0]]=1-exp(-nbDSB_*qal);
				}
				if(fertmoy_sig_het[genotype_indiv[1]]==0){
					fertmoy_sig_het[genotype_indiv[1]]=1-exp(-nbDSB_*qal);
				}
				if(fertmoy_sig_hom[genotype_indiv[0]]==0){
					fertmoy_sig_hom[genotype_indiv[0]]=1-exp(-nbDSB_*qal_sig_hom_0);
				}
				if(fertmoy_sig_hom[genotype_indiv[1]]==0){
					fertmoy_sig_hom[genotype_indiv[1]]=1-exp(-nbDSB_*qal_sig_hom_1);
				}
				
			}
		}
		qindmoy+=qal;
		fertindmoy+=(1-exp(-nbDSB_*qal));
	}
	res[-1].push_back(qindmoy/N_);
	res[-1].push_back(fertindmoy/N_);
	double qmoyallpop=0;
	double fertmoyallpop=0;
	for(auto const &it : qmoy){
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
	////////////////////////////////////////////////////////
	for(auto const &it : qmoy){
		if (qmoy_hom.find(it.first) != qmoy_hom.end()) {
			double qmoyall_hom=0;
			double fertmoyall_hom=0;
			for(int i=0; i<qmoy_hom[it.first].size(); i++){
				qmoyall_hom+=qmoy_hom[it.first][i];
				fertmoyall_hom+=fertmoy_hom[it.first][i];
			}
			res_hom[it.first].push_back(qmoyall_hom/ qmoy_hom[it.first].size());
			res_hom[it.first].push_back(fertmoyall_hom/ fertmoy_hom[it.first].size());
		}else{
			res_hom[it.first].push_back(0);
			res_hom[it.first].push_back(0);
		}
	}
	for(auto const &it : qmoy){
		if (qmoy_het.find(it.first) != qmoy_het.end()) {
			double qmoyall_het=0;
			double fertmoyall_het=0;
			for(int i=0; i<qmoy_het[it.first].size(); i++){
				qmoyall_het+=qmoy_het[it.first][i];
				fertmoyall_het+=fertmoy_het[it.first][i];
			}
			res_het[it.first].push_back(qmoyall_het/ qmoy_het[it.first].size());
			res_het[it.first].push_back(fertmoyall_het/ fertmoy_het[it.first].size());
		}else{
			res_het[it.first].push_back(0);
			res_het[it.first].push_back(0);
		}
	}
	
	/////////////////////////////////////////
	for(auto const &it : qmoy){
		if (qmoy_hemi.find(it.first) != qmoy_hemi.end()) {
			double qmoyall_hemi=0;
			double fertmoyall_hemi=0;
			for(int i=0; i<qmoy_hemi[it.first].size(); i++){
				qmoyall_hemi+=qmoy_hemi[it.first][i];
				fertmoyall_hemi+=fertmoy_hemi[it.first][i];
			}
			res_hemi[it.first].push_back(qmoyall_hemi/ qmoy_hemi[it.first].size());
			res_hemi[it.first].push_back(fertmoyall_hemi/ fertmoy_hemi[it.first].size());
		}else{
			res_hemi[it.first].push_back(0);
			res_hemi[it.first].push_back(0);
		}
	}
	//////////////////////////////////////////
	
	
	for(auto const &it : qmoy){
		double w_hom=0;
		double w_hemi=0;
		if(res_hom[it.first][1]==0){
			w_hom=fertmoy_sig_hom[it.first];
		}else{
			w_hom=res_hom[it.first][1];
		}
		if(res_hemi[it.first][1]==0){
			w_hemi=fertmoy_sig_het[it.first];
		}else{
			w_hemi=res_hemi[it.first][1];
		}
		if(w_hemi!=0){
			double sigma=(w_hom-w_hemi)/w_hemi;
			res[it.first].push_back(sigma);
		}else{
			res[it.first].push_back(0);
		}
	}
	
	
	////////////////////////
	//Print functions
	///////////////////////
	/*
	cout<<"res :"<<endl;
	for (auto const &it : res){
		cout<<it.first<<"=>";
		for(auto const &i : it.second){
			cout<<" "<<i;
		}
		cout<<endl;
	}
	cout<<'\n'<<endl;
	*/
	/*
	cout<<"res_hom :"<<endl;
	for (auto const &it : res_hom){
		cout<<it.first<<"=>";
		for(auto const &i : it.second){
			cout<<" "<<i;
		}
		cout<<endl;
	}
	cout<<'\n'<<endl;
	
	cout<<"res_het :"<<endl;
	for (auto const &it : res_het){
		cout<<it.first<<"=>";
		for(auto const &i : it.second){
			cout<<" "<<i;
		}
		cout<<endl;
	}
	cout<<'\n'<<endl;
	
	cout<<"res_hemi :"<<endl;
	for (auto const &it : res_hemi){
		cout<<it.first<<"=>";
		for(auto const &i : it.second){
			cout<<" "<<i;
		}
		cout<<endl;
	}
	cout<<'\n'<<endl;
	*/
	
	////////////////////////////////////////////////////////
	res[-2].push_back(qmoyallpop/(2*N_));
	res[-2].push_back(fertmoyallpop/(2*N_));
	res[-3].push_back(0);
	res[-3].push_back(0);
	res[-3].push_back(0);
	////////////////////////////////////////////////////////
	res_hom[-3].push_back(0);
	res_het[-3].push_back(0);
	res_hom[-3].push_back(0);
	res_het[-3].push_back(0);
	
	res_hemi[-3].push_back(0);/////
	res_hemi[-3].push_back(0);/////
	
	////////////////////////////////////////////////////////
	return vector<map<int,vector<double>>>{res,res_hom,res_het,res_hemi};
}


//Fonction calcul du q, fitness, nblinksite et nblinkpos de chaque allele de la population dans tous les contextes het possibles (alleles presents dans la pop + 1 nouvel allele) et en contexte homoaygote.


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

double Model::Mean_fert_new_allele(vector<vector<int>>* genotype, vector<vector<vector<int>>>* pop, map<int,vector<double>>* infoperallele_het){
	//give sites to the new allele
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> newpos = choosemany(nbsite_, freepos); 
	int newallele = nballele_;
	map<int,double> fertperall;
	map<int,double> qperall;
	int ind_gen=1;
	/////////////////////////////////////////////////////////////////////////////////
	if(targetcomp_){
		double is_diff=1;
		double ctot=ctot_;
		double threshold=cfreethreshold_;
		double cfree=1/float(1000)*ctot;
		for(auto quantile_index : nbsitesperquantile_){
			nbsitesperquantile_[quantile_index.first][1]=0;
		}
		if(quantilenb_!=0){
			for(auto const &it : newpos){
				bool is_quantile=0;
				for(auto quantile_index : nbsitesperquantile_){
					if(Affinity_[it]<=quantile_index.first and is_quantile==0){
						nbsitesperquantile_[quantile_index.first][1]+=4;//if this site is active, the category in which the affinity of this site enters gains one site
						is_quantile=1;
					}
				}
			}
		}
		double cfree_min=0;
		double cfree_max=ctot_;
		double sum_p_occup=0;
		while(is_diff==1){
			is_diff=0;
			sum_p_occup=0;
			//In the case of affinity categories
			if(quantilenb_!=0){
				for(auto quantile_index : nbsitesperquantile_){
					sum_p_occup+=double(nbsitesperquantile_[quantile_index.first][1]*(cfree*nbsitesperquantile_[quantile_index.first][0])/(1+cfree*nbsitesperquantile_[quantile_index.first][0]));//The probability of binding in each category multiplied by the number of sites in this category is added to the mean sum of occupied site
				}
			}else{
				for(auto const &it : newpos){//for each position
					double p_occup=cfree*Affinity_[it]/(1+cfree*Affinity_[it]);//the probability of binding is calculated
					sum_p_occup+=4*p_occup;
				}
			}
			if(cfree+sum_p_occup>ctot){
				//if cfree is to high
				is_diff=1;
				cfree_max=cfree;
				cfree=double(cfree_max+cfree_min)/2;
			}else{
				//if cfree is to low or the right one
				double new_cfree=ctot-sum_p_occup;
				if(abs(double(cfree)/double(ctot)-double(new_cfree)/double(ctot))<=threshold){
					//the difference between the current and the next cfree/ctot is below or equal to the thershold
					is_diff=0;
				}else{
					//the difference between the current and the next cfree/ctot is above the thershold
					is_diff=1;
					cfree_min=cfree;
					cfree=double(cfree_max+cfree_min)/2;
				}
			}
		}
		ind_gen=cfree;
	}
	/////////////////////////////////////////////////////////////////////////////////
	double qal;
	double mean_fert=0;
	//give the coefficient for the new allele (its a new allele so none of its site has been eroded in any individual, so its the same for everyone)
	vector<double> infonewallele;
	double xa=0;
	double xb=0;
	double moyprobalinkab=0;
	double moyprobalinka2b=0;
	double moyprobalinkb2a=0;
	double moyprobalinka=0;
	double moyprobalinkb=0;
	for(auto const &it1 : newpos){ //for the allele that does not change
		xa=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
		xb=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
		moyprobalinkab+=xa*xb;
		moyprobalinka2b+=puissance_double(2,xa)*xb;
		moyprobalinkb2a+=puissance_double(2,xb)*xa;
		moyprobalinka+=xa;
		moyprobalinkb+=xb;
	}
	infonewallele.push_back(moyprobalinkab/nbsite_);
	infonewallele.push_back(moyprobalinka2b/nbsite_);
	infonewallele.push_back(moyprobalinkb2a/nbsite_);
	infonewallele.push_back(moyprobalinka/nbsite_);
	infonewallele.push_back(moyprobalinkb/nbsite_);
	//for each allele, we take the first individual that possess this allele and calculate the coefficients for this allele
	for(auto const &it : (*infoperallele_het)){
		if(targetcomp_){
			ind_gen=(*infoperallele_het)[it.first][11];
		}
		int simplechrom;
		int chromtochange;
		if(it.first!=-3){
			vector<int>::iterator itv = find((*genotype)[parityIndex_].begin(),(*genotype)[parityIndex_].end(),it.first);
			if(itv!=(*genotype)[parityIndex_].end()){
				for(int ind=0; ind<2*N_; ind+=1){
					if((*genotype)[parityIndex_][ind]==it.first and ind%2==0){
						simplechrom=ind;
						chromtochange=ind+1;
					}else if((*genotype)[parityIndex_][ind]==it.first and ind%2==1){
						simplechrom=ind-1;
						chromtochange=ind;
					}
				}
				vector<double> infooldallele;
				for(auto const &it1 : Siteforeacheallele_[it.first]){ //for the allele that does not change
					if((*pop)[parityIndex_][simplechrom][it1]==1){
						xa=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
					}else{
						xa=0;
					}
					if((*pop)[parityIndex_][chromtochange][it1]==1){
						xb=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
					}else{
						xb=0;
					}
					moyprobalinkab+=xa*xb;
					moyprobalinka2b+=puissance_double(2,xa)*xb;
					moyprobalinkb2a+=puissance_double(2,xb)*xa;
					moyprobalinka+=xa;
					moyprobalinkb+=xb;
				}
				infooldallele.push_back(moyprobalinkab/nbsite_);
				infooldallele.push_back(moyprobalinka2b/nbsite_);
				infooldallele.push_back(moyprobalinkb2a/nbsite_);
				infooldallele.push_back(moyprobalinka/nbsite_);
				infooldallele.push_back(moyprobalinkb/nbsite_);
				qal=(4*infooldallele[0]-infooldallele[1]-infooldallele[2]+4*infonewallele[0]-infonewallele[1]-infonewallele[2])/(infooldallele[3]+infooldallele[4]+infonewallele[3]+infonewallele[4]);
				qperall[it.first]=qal;
				fertperall[it.first]=1-exp(-nbDSB_*qal);
				//mean_fert+=1-exp(-nbDSB_*qal);
				mean_fert+=(1-exp(-nbDSB_*qal))*freqallele(it.first, genotype);
			}
		}
	}
	//mean_fert=mean_fert/infoperallele_.size();
	//cout<<"mean_fert = "<<mean_fert<<endl;
	return mean_fert;
}




vector<double> Model::sigma_q_w_0(){
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> newpos = choosemany(2*nbsite_, freepos); 
	vector<double> res = vector<double>(7,0);
	double cfree_het=1;
	double cfree_hom=1;
	if(zygosity_){
		cfree_hom=2;
	}
	double cfree_ctot_het_final;
	double cfree_ctot_hom_final;
	bool even=0;
	//if(targetcomp_)
	if(targetcomp_){
		for(int allele=0; allele<2; allele++){
			double is_diff=1;
			double ctot_het=ctot_;
			double ctot_hom=ctot_*2;
			double threshold=cfreethreshold_;
			cfree_hom=1/float(1000)*ctot_hom;
			cfree_het=1/float(1000)*ctot_het;
			for(auto quantile_index : nbsitesperquantile_){
				nbsitesperquantile_[quantile_index.first][1]=0;
			}
			if(quantilenb_!=0){
				for(auto const &it : newpos){
					bool is_quantile=0;
					for(auto quantile_index : nbsitesperquantile_){
						if(Affinity_[it]<=quantile_index.first and is_quantile==0){
							nbsitesperquantile_[quantile_index.first][1]+=4;//if this site is active, the category in which the affinity of this site enters gains one site
							is_quantile=1;
						}
					}
				}
			}
			double cfree_min_het=0;
			double cfree_max_het=ctot_het;
			double cfree_min_hom=0;
			double cfree_max_hom=ctot_hom;
			double sum_p_occup_het=0;
			double sum_p_occup_hom=0;
			while(is_diff==1){
				is_diff=0;
				sum_p_occup_het=0;
				sum_p_occup_hom=0;
				//In the case of affinity categories
				if(quantilenb_!=0){
					for(auto quantile_index : nbsitesperquantile_){
						sum_p_occup_het+=double(nbsitesperquantile_[quantile_index.first][1]*(cfree_het*nbsitesperquantile_[quantile_index.first][0])/(1+cfree_het*nbsitesperquantile_[quantile_index.first][0]));//The probability of binding in each category multiplied by the number of sites in this category is added to the mean sum of occupied site
						sum_p_occup_hom+=double(nbsitesperquantile_[quantile_index.first][1]*(cfree_hom*nbsitesperquantile_[quantile_index.first][0])/(1+cfree_hom*nbsitesperquantile_[quantile_index.first][0]));//The probability of binding in each category multiplied by the number of sites in this category is added to the mean sum of occupied site
					}
				}
				else{
					for(auto const &it : newpos){//for each position
						double p_occup_het=cfree_het*Affinity_[it]/(1+cfree_het*Affinity_[it]);//the probability of binding is calculated
						double p_occup_hom=cfree_hom*Affinity_[it]/(1+cfree_hom*Affinity_[it]);//the probability of binding is calculated
						if(even==0){
							sum_p_occup_hom+=4*p_occup_hom;
							sum_p_occup_het+=4*p_occup_het;
						}
						even=(even+1)%2;
					}
				}
				if(cfree_hom+sum_p_occup_hom>ctot_hom){
					//if cfree is to high
					is_diff=1;
					cfree_max_hom=cfree_hom;
					cfree_hom=double(cfree_max_hom+cfree_min_hom)/2;
				}else{
					//if cfree is to low or the right one
					double new_cfree_hom=ctot_hom-sum_p_occup_hom;
					cfree_ctot_hom_final=double(new_cfree_hom)/double(ctot_hom);
					if(abs(double(cfree_hom)/double(ctot_hom)-double(new_cfree_hom)/double(ctot_hom))<=threshold){
						//the difference between the current and the next cfree/ctot is below or equal to the thershold
						is_diff=0;
					}else{
						//the difference between the current and the next cfree/ctot is above the thershold
						is_diff=1;
						cfree_min_hom=cfree_hom;
						cfree_hom=double(cfree_max_hom+cfree_min_hom)/2;
					}
				}
				if(cfree_het+sum_p_occup_het>ctot_het){
					//if cfree is to high
					is_diff=1;
					cfree_max_het=cfree_het;
					cfree_het=double(cfree_max_het+cfree_min_het)/2;
				}else{
					//if cfree is to low or the right one
					double new_cfree_het=ctot_het-sum_p_occup_het;
					cfree_ctot_het_final=double(new_cfree_het)/double(ctot_het);
					if(abs(double(cfree_het)/double(ctot_het)-double(new_cfree_het)/double(ctot_het))<=threshold){
						//the difference between the current and the next cfree/ctot is below or equal to the thershold
						is_diff=0;
					}else{
						//the difference between the current and the next cfree/ctot is above the thershold
						is_diff=1;
						cfree_min_het=cfree_het;
						cfree_het=double(cfree_max_het+cfree_min_het)/2;
					}
				}
			}
		}
	}
	/////////////////////////////////////////////////////////////
	if(targetcomp_){
		res[0]=cfree_ctot_het_final;
		res[1]=cfree_ctot_hom_final;
	}else{
		res[0]=0;
		res[1]=0;
	}
	
	//q and w
	
	//qhet0 and whet0
	vector<double> infonewallele_1;
	vector<double> infonewallele_2;
	vector<double> infonewallele;
	even=0;
	for(int allele=0; allele<3; allele++){
		double ind_gen;
		if(allele<2){
			ind_gen=cfree_het;
		}else{
			ind_gen=cfree_hom;
		}
		double xa=0;
		double xb=0;
		double moyprobalinkab=0;
		double moyprobalinka2b=0;
		double moyprobalinkb2a=0;
		double moyprobalinka=0;
		double moyprobalinkb=0;
		for(auto const &it1 : newpos){
			if(even==0){
				xa=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
				xb=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
				moyprobalinkab+=xa*xb;
				moyprobalinka2b+=puissance_double(2,xa)*xb;
				moyprobalinkb2a+=puissance_double(2,xb)*xa;
				moyprobalinka+=xa;
				moyprobalinkb+=xb;
			}
			even=(even+1)%2;
		}
		if(allele==0){
			infonewallele_1.push_back(moyprobalinkab/nbsite_);
			infonewallele_1.push_back(moyprobalinka2b/nbsite_);
			infonewallele_1.push_back(moyprobalinkb2a/nbsite_);
			infonewallele_1.push_back(moyprobalinka/nbsite_);
			infonewallele_1.push_back(moyprobalinkb/nbsite_);
		}else if (allele==1){
			infonewallele_2.push_back(moyprobalinkab/nbsite_);
			infonewallele_2.push_back(moyprobalinka2b/nbsite_);
			infonewallele_2.push_back(moyprobalinkb2a/nbsite_);
			infonewallele_2.push_back(moyprobalinka/nbsite_);
			infonewallele_2.push_back(moyprobalinkb/nbsite_);
		}else{
			infonewallele.push_back(moyprobalinkab/nbsite_);
			infonewallele.push_back(moyprobalinka2b/nbsite_);
			infonewallele.push_back(moyprobalinkb2a/nbsite_);
			infonewallele.push_back(moyprobalinka/nbsite_);
			infonewallele.push_back(moyprobalinkb/nbsite_);
		}
		if(allele==0){
			even=1;
		}else{
			even=0;
		}
	}
	double qhet=(4*infonewallele_1[0]-infonewallele_1[1]-infonewallele_1[2]+4*infonewallele_2[0]-infonewallele_2[1]-infonewallele_2[2])/(infonewallele_1[3]+infonewallele_1[4]+infonewallele_2[3]+infonewallele_2[4]);
	double whet=1-exp(-nbDSB_*qhet);
	res[2]=qhet;
	res[4]=whet;
	double qhom=(4*infonewallele[0]-infonewallele[1]-infonewallele[2])/(infonewallele[3]+infonewallele[4]);
	double whom=1-exp(-nbDSB_*qhom);
	res[3]=qhom;
	res[5]=whom;
	
	//sigma0
	double sigma=(whom-whet)/whet;
	res[6]=sigma;

	return res;
}



/*double Model::sigma_0(){
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> newpos = choosemany(2*nbsite_, freepos); 
	int ind_gen=1;
	//////////////////////////////////////////if(targetcomp_)
	//w_het(0,0)
	vector<double> infonewallele_1;
	vector<double> infonewallele_2;
	bool even=0;
	for(int allele=0; allele<2; allele++){
		double xa=0;
		double xb=0;
		double moyprobalinkab=0;
		double moyprobalinka2b=0;
		double moyprobalinkb2a=0;
		double moyprobalinka=0;
		double moyprobalinkb=0;
		for(auto const &it1 : newpos){
			if(even==0){
				xa=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
				xb=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
				moyprobalinkab+=xa*xb;
				moyprobalinka2b+=puissance_double(2,xa)*xb;
				moyprobalinkb2a+=puissance_double(2,xb)*xa;
				moyprobalinka+=xa;
				moyprobalinkb+=xb;
			}
			even=(even+1)%2;
		}
		if(allele==0){
			infonewallele_1.push_back(moyprobalinkab/nbsite_);
			infonewallele_1.push_back(moyprobalinka2b/nbsite_);
			infonewallele_1.push_back(moyprobalinkb2a/nbsite_);
			infonewallele_1.push_back(moyprobalinka/nbsite_);
			infonewallele_1.push_back(moyprobalinkb/nbsite_);
		}else{
			infonewallele_2.push_back(moyprobalinkab/nbsite_);
			infonewallele_2.push_back(moyprobalinka2b/nbsite_);
			infonewallele_2.push_back(moyprobalinkb2a/nbsite_);
			infonewallele_2.push_back(moyprobalinka/nbsite_);
			infonewallele_2.push_back(moyprobalinkb/nbsite_);
		}
		even=1;
	}
	double qhet=(4*infonewallele_1[0]-infonewallele_1[1]-infonewallele_1[2]+4*infonewallele_2[0]-infonewallele_2[1]-infonewallele_2[2])/(infonewallele_1[3]+infonewallele_1[4]+infonewallele_2[3]+infonewallele_2[4]);
	double whet=1-exp(-nbDSB_*qhet);
	
	//whom(0)
	ind_gen=2;
	vector<double> infonewallele;
	even=0;
	double xa=0;
	double xb=0;
	double moyprobalinkab=0;
	double moyprobalinka2b=0;
	double moyprobalinkb2a=0;
	double moyprobalinka=0;
	double moyprobalinkb=0;
	for(auto const &it1 : newpos){
		if(even==0){
			xa=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
			xb=ind_gen*Affinity_[it1]/(1+ind_gen*Affinity_[it1]);
			moyprobalinkab+=xa*xb;
			moyprobalinka2b+=puissance_double(2,xa)*xb;
			moyprobalinkb2a+=puissance_double(2,xb)*xa;
			moyprobalinka+=xa;
			moyprobalinkb+=xb;
		}
		even=(even+1)%2;	
	}
	infonewallele.push_back(moyprobalinkab/nbsite_);
	infonewallele.push_back(moyprobalinka2b/nbsite_);
	infonewallele.push_back(moyprobalinkb2a/nbsite_);
	infonewallele.push_back(moyprobalinka/nbsite_);
	infonewallele.push_back(moyprobalinkb/nbsite_);
	double qhom=(4*infonewallele[0]-infonewallele[1]-infonewallele[2])/(infonewallele[3]+infonewallele[4]);
	double whom=1-exp(-nbDSB_*qhom);
	
	double sigma=(whom-whet)/whet;
	return sigma;
}*/


