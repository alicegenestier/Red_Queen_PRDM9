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

using namespace std;

//============================
//         Constructors
//============================
Model::Model(): N_(1000), L_(100000), nbsite_(400), indPrdm9_(5), nballele_(1), parityIndex_(0), v_(1e-7), u_(1e-6), meanaff_(0.6), varaff_(1), nbDSB_(6), nbGenerations_(2), ismigration_(false) {
	
	//vector counting the number of failed meiosis per generation
	nbfailedmeiosis_=vector<vector<int>>(nbGenerations_,vector<int>(4,0));
	
	//popluations matrix
	populations_ = vector<vector<vector<int>>>(2,vector<vector<int>>(2*N_,vector<int> (L_,1)));
	// matrix initialized with 1 because at the begining all the sites are activated in the whole genme for each individual
	for (auto &i : populations_){
		for (int j=0; j<2*N_; j++)
		i[j][indPrdm9_]=0; 
	}

	//parity = 0 ou 1 => pour 2 pop : 
	//int current = parity
	//int next = 1-parity
	
	//affinity vector
	Affinity_=vector<double>(L_,0);
	for(int i=0; i<L_; i++){
		Affinity_[i]=chosegamma(meanaff_, varaff_);
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
	
	// positions of sites for first allele (0)
	vector<int> freepos = vectfreesites(Alleleforeachpos_, -1);
	vector<int> firstpos = choosemany(nbsite_, freepos); 
	
	for(auto i : firstpos){
		Alleleforeachpos_[i]=0;
	}
	
	// positions map
	for (int i = 0; i < nballele_; i++){
		Siteforeacheallele_[i]=firstpos;
	}
	
	if(ismigration_){
		populations1_ = vector<vector<vector<int>>>(2,vector<vector<int>>(2*N_,vector<int> (L_,1)));
		for (auto &i : populations1_){
			for (int j=0; j<2*N_; j++)
			i[j][indPrdm9_]=0; 
		}
		genotypes1_=vector<vector<int>>(2,vector<int>(2*N_,0));
	}

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
vector<vector<int>> Model::genotypes(){
	return genotypes_;
}
vector<vector<int>> Model::genotypes1(){
	return genotypes1_;
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
float Model::u(){
	return u_;
}
float Model::v(){
	return v_;
}
float Model::m(){
	return m_;
}
double Model::meanaff(){
	return meanaff_;
} 
double Model::varaff(){
	return varaff_;
}
float Model::nbDSB(){
	return nbDSB_;
}
int Model::nbGenerations(){
	return nbGenerations_;
}
vector<vector<int>> Model::nbfailedmeiosis(){
	return nbfailedmeiosis_;
}
bool Model::ismigration(){
	return ismigration_;
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
double Model::chosegamma(double meanaff, double varaff){
	gamma_distribution<double> g(varaff,meanaff);
	return g(e);
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

vector<int> Model::occupiedsites(vector<int> vect){//return the index of all positions occupied
	vector<int> occupiedsites;
	for(int i=0; i<vect.size(); i++){
		if(vect[i]!= -2 and vect[i]!= -1){
			occupiedsites.push_back(i);
		}
	}
	return (occupiedsites);	
}

void Model::sitemutation(){
	// for each position with allele, prob v to mutate and if mutation, choose randomly 1 chrom to mutate
	vector<int> occupied = occupiedsites(Alleleforeachpos_);
	
	// affichage
	/*for (auto i : occupied){
		cout<<" "<< i;
	}
	cout<<endl;*/
	
	vector<int> mutsites;
	for (auto i : occupied){
		if (bernoulli_draw(2*N_*v_)){
			mutsites.push_back(i);
			/*for (auto j : mutsites){
				cout<<" "<< j;
			}
			cout<<endl;*/
		}
	}
	for (auto j : mutsites){
		int mutchrom = choose(2*N_);
		if (populations_[parityIndex_][mutchrom][j]==1){ 
			populations_[parityIndex_][mutchrom][j]=0;
		}
	}
}

void Model::allelemutation(){
	// for each chromosome, mutation of allele with proba u and if at least one mutation (if more than one : same allele ?) update populations (change prdm9 allele at its position), genotypes(same), map (site pos associated to the new allele) and allele for each pos (change the allel corresponding to each pos : -1 -> new allele)
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
		genotypes_[parityIndex_][i]=newallele;//update genotype
		populations_[parityIndex_][i][indPrdm9_]=newallele;//update pop
	}
}

void Model::updatemissingallele(){
//fonction mise a jours map si allele disparu
//comparaison entre clefs map et genotype : si un ou plusieurs alleles sont dans la map mais plus dans genotype ou pop -> on les supprime de la map et on remet toutes les pos correspondantes a -1 dans allelforeachpos et on remet toutes les activations des pos a 1
	for(auto const &it : Siteforeacheallele_){
		//cout<<it.first<<endl;
		if(it.first!=-3 and it.first!=-2){
			//cout<<it.first<<endl;
			vector<int>::iterator itvect = find(genotypes_[parityIndex_].begin(),genotypes_[parityIndex_].end(),it.first);
			if(itvect == genotypes_[parityIndex_].end()){
				for(auto const &i : it.second){
					//cout<<"i"<<i<<endl;
					Alleleforeachpos_[i]=-1;//free all sites
					for(int j=0; j<2*N_; j++)
						populations_[parityIndex_][j][i]=1;//reactive all sites
				}
				Siteforeacheallele_.erase(it.first);//erase in the map
			}
		}
	}
}

//Print functions

//print populations
void Model::printpop(int n){//n = parityIndexu (0 ou 1) ou 2 si plot 2 pop
	if(n==2){
		for (auto i : populations_){ // a changer
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
		vector<vector<int>> pop = populations_[n]; // a changer
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
void Model::printgen(int n){
	if(n==2){
		for (auto i : genotypes_){
			for (auto j : i){
				cout<<' '<<j;
			}
			cout<<'\n';
		}
		cout<<'\n';
	}
	else if(n==0 or n==1){
		vector<int> gen = genotypes_[n];
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


// Meiosis
int Model::Meiosis(int no_chrom_ind, int nb_gen){
	// choose one ind
	int indiv = 2*choose(N_);
	//cout<<"indiv "<<indiv<<endl;
	//homozygote or heterozygote
	int ind_gen=2;
	vector<int> zygote{genotypes_[parityIndex_][indiv]};
	/*for(auto i : genotypes_){
		for (auto j : i){cout<<"j"<<j;}
	}*/
	if(genotypes_[parityIndex_][indiv]!=genotypes_[parityIndex_][indiv+1]){
		ind_gen=1;
		zygote.push_back(genotypes_[parityIndex_][indiv+1]);
	}
	//cout<<"ind_gen : "<<ind_gen<<endl;
	//cout<<"zygote : ";
	/*for (auto i : zygote){
		cout<<" "<<i;
	}
	cout<<endl;*/
	
	//PRDM9 binding
	vector<vector<int>>summarysites;
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
				if(populations_[parityIndex_][chrom][i]==1){ // a changer
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
			}
		}
	}
	//nblinksite=summarysites.size();
	//cout<<"+++++++++++++++++++++++++++++++++++++++++"<<endl;
	/*for(auto i : summarysites){
		for(auto j : i){
			cout<<' '<<j;
		}
		cout<<"\n";
	}
	cout<<endl;*/
	cout<<"nblinksite "<<nblinksite<<endl;
	//DSB
	double pDSB=double(nbDSB_)/nblinksite;
	//cout<<"pDSB : "<<pDSB<<endl;
	//cout<<"nblinksite : "<<nblinksite<<endl;
	//cout<<"nb linked sites : "<<nblinksite<<endl;
	//cout<<"pDSB : "<<pDSB<<endl;
	vector<vector<int>>vectsitedsb;
	int checknblink =0;
	vector<vector<int>> vect_CO;
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
				/*cout<<"summarysites : "<<endl;
				for(auto i : summarysites){
					for(auto j : i){
						cout<<' '<<j;
					}
					cout<<"\n";
				}
				cout<<endl;
				cout<<"a1 "<<nbdsbpersite<<endl;*/
				try{
					if (nbdsbpersite>1){
						nbfailedmeiosis_[nb_gen][0]+=1;
						throw int(0);
					}
				} // assert
				catch(int const& error_nb){
					if(error_nb==0){
						cerr << "2 DSB on one site" << endl;
					}
					return-1;
				}
			}
			//cout<<"a2"<<nbdsbpersite<<endl;
		}
		if (dsb){
			if(vectsitedsb.back()[1]==0 or vectsitedsb.back()[1]==1){
				if(summarysites[i][3]==1){
					vect_CO.push_back({summarysites[i][0],vectsitedsb.back()[1],2});
				}
				if(summarysites[i][4]==1){
					vect_CO.push_back({summarysites[i][0],vectsitedsb.back()[1],3});
				}
			}else if(vectsitedsb.back()[1]==2 or vectsitedsb.back()[1]==3){
				if(summarysites[i][1]==1){
					vect_CO.push_back({summarysites[i][0],vectsitedsb.back()[1],0});
				}
				if(summarysites[i][2]==1){
					vect_CO.push_back({summarysites[i][0],vectsitedsb.back()[1],1});
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
		cout<<endl;*/
		/*cout<<"vect_CO : "<<endl;
		for(auto i : vect_CO){
			for(auto j : i){
				cout<<' '<<j;
			}
			cout<<"\n";
		}
		cout<<endl;*/
	}
	try{
		if(vectsitedsb.size()==0){
			nbfailedmeiosis_[nb_gen][1]+=1;
			throw int(1);
		}else if(not vect_CO.size()){
			nbfailedmeiosis_[nb_gen][2]+=1;
			throw int(2);
		}
	} // assert
	catch(int const& error_nb){
		if(error_nb==1){
			cerr << "No DSB" << endl;
		}else if(error_nb==2){
			cerr << "No symmetrical sites (binding + DSB)" << endl;
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
	cout<<endl;*/
	/*cout<<"vect_CO : "<<endl;
		for(auto i : vect_CO){
			for(auto j : i){
				cout<<' '<<j;
			}
			cout<<"\n";
		}
	cout<<endl;*/
	/*cout<<"vectsitedsb : "<<endl;
	for(auto i : vectsitedsb){
		for(auto j : i){
			cout<<' '<<j;
		}
		cout<<"\n";
	}
	cout<<endl;*/
	
	vector<int> index_CO=vect_CO[choose(vect_CO.size())];
	//cout<<"index_CO : "<<endl;
	/*for(auto i : index_CO){
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
		//remplacement et crossing over
		int indvectsite=0;
		//faire d'abord CO
		for (int ind_site_pop=0; ind_site_pop<index_CO[0]; ind_site_pop++){
			if(indvectsite!=vectsitedsb.size() and ind_site_pop==summarysites[vectsitedsb[indvectsite][0]][0]){
				if(vectsitedsb[indvectsite][1]==no_chromatide){
					populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_homologue_chrom][ind_site_pop];
					//cout<<"pop1 "<<populations_[parityIndex_][no_homologue_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
				}else{
					populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_current_chrom][ind_site_pop];
				//cout<<"pop1 "<<populations_[parityIndex_][no_current_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
				}
				indvectsite++;
			}else{
				populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_current_chrom][ind_site_pop];
				//cout<<"pop1 "<<populations_[parityIndex_][no_current_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
			}
		}
		if(no_chromatide==index_CO[1]){
			no_chromatide==index_CO[2];
		}else{
			no_chromatide==index_CO[1];
		}
		for (int ind_site_pop=index_CO[0]; ind_site_pop<L_; ind_site_pop++){
			if(indvectsite!=vectsitedsb.size() and ind_site_pop==summarysites[vectsitedsb[indvectsite][0]][0]){
				if(vectsitedsb[indvectsite][1]==no_chromatide){
					populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_current_chrom][ind_site_pop];
					//cout<<"pop1 "<<populations_[parityIndex_][no_current_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
				}else{
					populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_homologue_chrom][ind_site_pop];
					//cout<<"pop1 "<<populations_[parityIndex_][no_homologue_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
				}
				indvectsite++;
			}else{
				populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_homologue_chrom][ind_site_pop];
				//cout<<"pop1 "<<populations_[parityIndex_][no_homologue_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
			}
		}
	}else{
		//remplacement
		//cout<<"chromatide without CO"<<endl;
		int indvectsite=0;
		for (int ind_site_pop=0; ind_site_pop<L_; ind_site_pop++){
			//cout<<"indvectsite "<<indvectsite<<endl;
			//cout<<"ind_site_pop "<<ind_site_pop<<endl;
			if(indvectsite!=vectsitedsb.size() and ind_site_pop==summarysites[vectsitedsb[indvectsite][0]][0]){
				if(vectsitedsb[indvectsite][1]==no_chromatide){
					populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_homologue_chrom][ind_site_pop];
					//cout<<"pop1 "<<populations_[parityIndex_][no_homologue_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
				}else{
					populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_current_chrom][ind_site_pop];
					//cout<<"pop1 "<<populations_[parityIndex_][no_current_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
				}
				indvectsite++;
			}else{
				populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]=populations_[parityIndex_][no_current_chrom][ind_site_pop];
				//cout<<"pop1 "<<populations_[parityIndex_][no_current_chrom][ind_site_pop]<<" pop2 "<<populations_[(parityIndex_+1)%2][no_chrom_ind][ind_site_pop]<<endl;
			}
		}
	}
	
	//return nb meiose echouee?? Pour le moment, return 0 si eurreur qqpart et 1 si meiose reussie.
	return 1;
}


//methode remplissage nouvelle pop
void Model::fillnewpop(int nb_gen){
	//for(int indnewpop=0; indnewpop<2*N_; indnewpop++){
	for(int indnewpop=0; indnewpop<2*N_; indnewpop++){
		int meiosisState = Meiosis(indnewpop, nb_gen);
		while (meiosisState==-1){
			meiosisState = Meiosis(indnewpop, nb_gen);
			nbfailedmeiosis_[nb_gen][3]+=1;
		}
		//cout<<"-----------------------------------------------------------------------------------------------------"<<endl;
	}
	for(int indpop=0; indpop<2*N_;indpop++){
		genotypes_[(parityIndex_+1)%2][indpop]=populations_[(parityIndex_+1)%2][indpop][indPrdm9_];
	}
	parityIndex_=(parityIndex_+1)%2;
}

//methode qui repete tout ce au'on vient de faire pendant X generations
void Model::manygenerations(){
	for(int indgeneration=0; indgeneration<nbGenerations_; indgeneration++){
		//cout<<"parity index : "<<parityIndex_<<endl;
		cout<<"generation "<<indgeneration<<endl;
		sitemutation();
		allelemutation();
		updatemissingallele();
		/*cout<<"pop :"<<endl;
		printpop(parityIndex_);
		cout<<"map :"<<endl;
		printposallele();
		cout<< "genotypes : "<<endl;
		printgen(parityIndex_);
		cout<< "Allele for each position : "<<endl;
		printallelepos();*/
		fillnewpop(indgeneration);
		updatemissingallele();
		/*cout<<"---------"<<endl;
		cout<<"pop :"<<endl;
		printpop(parityIndex_);
		cout<<"map :"<<endl;
		printposallele();
		cout<< "genotypes : "<<endl;
		printgen(parityIndex_);
		cout<< "Allele for each position : "<<endl;
		printallelepos();*/
		cout<<"nb Failed meiosis : 2 DSB on one site, No DSB, No symmetrical sites (binding + DSB), Total"<<endl;
		for(auto indgen : nbfailedmeiosis_){
			for(auto indexfailed : indgen){
				cout<<' '<<indexfailed;
			}
			cout<<'\n';
		}
		cout<<endl;
	}
}

void migration(){
	//choisir le nombre de migrant voulu dans la pop (m_*N_ migrants parmi N_ indiv et independants)
	//copier les deux chromosomes de l'indiv de la pop1, remplacer cet indiv dans la pop1 par l'indiv choisi dans la pop2 puis remplacer l'indiv de la pop2 par l'indiv copier de la pop1.
	//faire ce meme echange dans les genotypes
}

//vecteur : num generation ; nbtotalallele par generation ; diversite ; activite moy ; 
