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

using namespace std;

//============================
//         Constructors
//============================
Model::Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen,string name): N_(N), L_(L), nbsite_(nbsite), indPrdm9_(indPrdm9), nballele_(nballele), parityIndex_(parityIndex), v_(v), u_(u), w_(w), meanaff_(meanaff), varaff_(varaff), nbDSB_(nbDSB), nbGenerations_(nbGenerations), ismigration_(ismigration), zygosity_(zygosity), withDSB_(withDSB), everygen_(everygen),name_(name) {
	
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
bool Model::withDSB(){
	return withDSB_;
}
double Model::w(){
	return w_;
}
string Model::name(){
	return name_;
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

void Model::sitemutation(){
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
		int mutchrom = choose(2*N_);
		//vector<int>::iterator itv = find(Siteforeacheallele_[-3].begin(),Siteforeacheallele_[-3].end(),j);
		//if(itv == Siteforeacheallele_[-3].end()){
			if (populations_[parityIndex_][mutchrom][j]==1){ 
				populations_[parityIndex_][mutchrom][j]=0;
			}
		//}else{
		//	populations_[parityIndex_][mutchrom][j]=(populations_[parityIndex_][mutchrom][j]+1)%2;
		//}
	}
	for (auto j : mutsitesneutral){
		int mutchrom = choose(2*N_);
		//vector<int>::iterator itv = find(Siteforeacheallele_[-3].begin(),Siteforeacheallele_[-3].end(),j);
		//if(itv == Siteforeacheallele_[-3].end()){
		//	if (populations_[parityIndex_][mutchrom][j]==1){ 
		//		populations_[parityIndex_][mutchrom][j]=0;
		//	}
		//}else{
			populations_[parityIndex_][mutchrom][j]=(populations_[parityIndex_][mutchrom][j]+1)%2;
		//}
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
	vector<int> alleletoerase;
	for(auto &it : Siteforeacheallele_){
		//cout<<it.first<<endl;
		if(it.first!=-3 and it.first!=-2){
			//cout<<it.first<<endl;
			//cout<<&it<<endl;
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
		}
	}
	for (auto a : alleletoerase){
		Siteforeacheallele_.erase(a);
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
		if(zygosity_){
			ind_gen=1;
		}
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
						nbfailedmeiosis_[nb_gen][0]+=1;
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
	q_=q_+double(vect_CO.size())/(vectsitedsb.size());
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
		copy( populations_[parityIndex_][no_current_chrom].begin(), populations_[parityIndex_][no_current_chrom].begin()+index_CO[0], populations_[(parityIndex_+1)%2][no_chrom_ind].begin() );
		copy( populations_[parityIndex_][no_homologue_chrom].begin()+index_CO[0], populations_[parityIndex_][no_homologue_chrom].end(), populations_[(parityIndex_+1)%2][no_chrom_ind].begin()+index_CO[0] );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide and summarysites[index_dsb[0]][0]<index_CO[0]){
				populations_[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=populations_[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
			}else if((index_dsb[1]==index_CO[1] or index_dsb[1]==index_CO[2]) and summarysites[index_dsb[0]][0]>index_CO[0]){
				populations_[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=populations_[parityIndex_][no_current_chrom][summarysites[index_dsb[0]][0]];
			}
		}
	}else{
		//cout<<"chromatide without CO"<<endl;
		copy( populations_[parityIndex_][no_current_chrom].begin(), populations_[parityIndex_][no_current_chrom].end(), populations_[(parityIndex_+1)%2][no_chrom_ind].begin() );
		for(auto index_dsb : vectsitedsb){
			if(index_dsb[1]==no_chromatide){
				populations_[(parityIndex_+1)%2][no_chrom_ind][summarysites[index_dsb[0]][0]]=populations_[parityIndex_][no_homologue_chrom][summarysites[index_dsb[0]][0]];
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
	ofstream generalfile ((name_+".trace").c_str());
	generalfile << "Generation number" << '\t' << "Total number of allele" << '\t' << "Diversity" << '\t'  << "Activity" << '\t' <<"Time" << '\t' << "Fertility rate" << '\t' << "2 DSB on one site rate" << '\t' << "No DSB rate" << '\t' << "No symmetrical sites (binding + DSB) rate" << '\t' << "q" <<'\n';
    	generalfile.flush();
    	ofstream allelefile ((name_+".allele").c_str());
	allelefile << "Generation number" << '\t' << "Allele number" << '\t' << "Frequency" << '\t'  << "Activity" << '\n';
    	allelefile.flush();
	ofstream paramsfile ((name_+".params").c_str());
	paramsfile << "N" << '\t' << N_ << '\n' << "L" << '\t' << L_ << '\n' << "nbsite" << '\t' << nbsite_ << '\n' << "indPrdm9" << '\t' << indPrdm9_ << '\n' << "nballele" << '\t' << nballele_ << '\n' << "parityIndex" << '\t' << parityIndex_ << '\n' << "u" << '\t' << u_ << '\n' << "v" << '\t' << v_ << '\n' << "w" << '\t' << w_ << '\n' << "meanaff" << '\t' << meanaff_ << '\n' << "varaff" << '\t' << varaff_ << '\n' << "nbDSB" << '\t' << nbDSB_ << '\n' << "nbGenerations" << '\t' << nbGenerations_ << '\n' << "ismigration" << '\t' << ismigration_ << '\n' << "zygosity" << '\t' << zygosity_ << '\n' << "withDSB" << '\t' << withDSB_ << '\n' << "everygen" << '\t' << everygen_ << '\n';
    	paramsfile.flush();
	
	for(int indgeneration=0; indgeneration<nbGenerations_; indgeneration++){
		clock_t t1, t2;
		t1=clock();
		//cout<<"parity index : "<<parityIndex_<<endl;
		cout<<"generation "<<indgeneration<<endl;
		sitemutation();
		allelemutation();
		updatemissingallele();
		/*cout<<"pop :"<<endl;
		printpop(parityIndex_);
		cout<<"map :"<<endl;
		printposallele();*/
		/*cout<< "genotypes : "<<endl;
		printgen(parityIndex_);
		cout<< "Allele for each position : "<<endl;
		printallelepos();*/
		q_=0;
		fillnewpop(indgeneration);
		q_=q_/(2*N_);
		//cout<<q_<<endl;
		updatemissingallele();
		/*cout<<"---------"<<endl;
		cout<<"pop :"<<endl;
		printpop(parityIndex_);
		cout<<"map :"<<endl;
		printposallele();*/
		/*cout<< "genotypes : "<<endl;
		printgen(parityIndex_);*/
		/*cout<< "Allele for each position : "<<endl;
		printallelepos();*/
		//cout<<"a "<<indgeneration % everygen_<<endl;
		t2=clock();
		if (indgeneration % everygen_ ==0)  {
			//cout<<"b "<<indgeneration % everygen_<<endl;
        		generalfile << indgeneration << '\t' << get_allele_number() << '\t' << get_current_diversity() << '\t'  << get_current_activity() << '\t' << (float)(t2-t1)/CLOCKS_PER_SEC << '\t' << 1-(double(nbfailedmeiosis_[indgeneration][3])/(2*N_+nbfailedmeiosis_[indgeneration][3]))<< '\t' << double(nbfailedmeiosis_[indgeneration][0])/(2*N_+nbfailedmeiosis_[indgeneration][3]) << '\t'<< double(nbfailedmeiosis_[indgeneration][1])/(2*N_+nbfailedmeiosis_[indgeneration][3]) << '\t' << double(nbfailedmeiosis_[indgeneration][2])/(2*N_+nbfailedmeiosis_[indgeneration][3]) << '\t' << q_ <<'\n'; //ajouter taux meiose echouee du a DSB, no sym...
            		generalfile.flush();
            		for (auto const &it : Siteforeacheallele_){
	    				if(it.first!=-2){
            				allelefile << indgeneration << '\t' << it.first << '\t' << freqall(it.first) << '\t'  << actall(it.first) << '\n';
            				allelefile.flush();
            			}
            		}
        	}
	}
	generalfile.close();
	allelefile.close();
	paramsfile.close();
	/*cout<<"nb Failed meiosis : 2 DSB on one site, No DSB, No symmetrical sites (binding + DSB), Total"<<endl;
	for(auto indgen : nbfailedmeiosis_){
		for(auto indexfailed : indgen){
			cout<<' '<<indexfailed;
		}
		cout<<'\n';
	}
	cout<<endl;*/
}

int Model::get_allele_number(){
	return Siteforeacheallele_.size()-2;
}

double Model::freqallele(int allelename){
	/*if(allelename==-3){
		return double(0);
	}*/
	if(allelename==0){
		//cout<<"count : "<<count(genotypes_[parityIndex_].begin(), genotypes_[parityIndex_].end(), allelename)-1<<endl;
		//cout<<"size genotype : "<<genotypes_[parityIndex_].size()-1<<endl;
		return double(count(genotypes_[parityIndex_].begin(), genotypes_[parityIndex_].end(), allelename)-1)/(genotypes_[parityIndex_].size()-1);
	}else{
		return double (count(genotypes_[parityIndex_].begin(), genotypes_[parityIndex_].end(), allelename))/(genotypes_[parityIndex_].size()-1);
	}
}

double Model::get_current_diversity(){
	double sumfreq=0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			//cout<<"freqallele : "<<freqallele(it.first)<<endl;
			sumfreq+=(freqallele(it.first))*(freqallele(it.first));
		}
	}
	return 1/sumfreq;
}

double Model::activitymoyallele(int allele){
	double moyact=0;
	//cout<<"allele : "<<allele<<endl;
	for(auto all : Siteforeacheallele_[allele]){
		double moyactsite=0;
		for(int i=0; i<2*N_; i++){
			if(populations_[parityIndex_][i][all]==1){
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

vector<double> Model::freqneutral(){
	double moyfreq=0;
	double moy2f=0;
	for(auto all : Siteforeacheallele_[-3]){
		double moyactsite=0;
		for(int i=0; i<2*N_; i++){
			if(populations_[parityIndex_][i][all]==1){
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

double Model::freqall(int allele){
	if (allele==-3){
		return freqneutral()[0];
	}else{
		return freqallele(allele);
	}
}

double Model::actall(int allele){
	if (allele==-3){
		return freqneutral()[1];
	}else{
		return activitymoyallele(allele);
	}
}

double Model::get_current_activity(){
	double moytotact=0;
	for (auto const &it : Siteforeacheallele_){
		if(it.first!=-3 and it.first!=-2){
			//cout<<"it.first : "<<it.first<<endl;
			//cout<<"activitymoyallele(it.first) : "<<activitymoyallele(it.first)<<endl;
			moytotact+=activitymoyallele(it.first)*freqallele(it.first);
			//cout<<"d"<<endl;
		}
	}
	//cout<<"activity moy"<<moytotact/meanaff_<<endl;
	return moytotact;
}

void migration(){
	//choisir le nombre de migrant voulu dans la pop (m_*N_ migrants parmi N_ indiv et independants)
	//copier les deux chromosomes de l'indiv de la pop1, remplacer cet indiv dans la pop1 par l'indiv choisi dans la pop2 puis remplacer l'indiv de la pop2 par l'indiv copier de la pop1.
	//faire ce meme echange dans les genotypes
}

//vecteur : num generation ; nbtotalallele par generation ; diversite ; activite moy 
//tableau : generation ; no allele ; freq ; activite moy 
