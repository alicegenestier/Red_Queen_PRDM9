
//============================
//          Includes
//============================
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

class Model {
	public:
	//============================
	//        Constructors
	//============================
	// Default	
	Model();
	//with arg
	Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen,string name);
	
	//============================
	//        Destructors
	//============================
	
	//============================
	//           Getters
	//============================
	vector<vector<vector<int>>> populations(); 
	vector<vector<vector<int>>> populations1();
	vector<vector<int>> genotypes(); 
	vector<vector<int>> genotypes1();
	map<int,vector<int>> Siteforeacheallele(); 
	vector<int> Alleleforeachpos();
	vector<double> Affinity();
	int parityIndex();
	int N(); 
	int L(); 
	int indPrdm9(); 
	int nballele(); 
	int nbsite(); 
	double u();
	double v();
	double m();
	double meanaff(); 
	double varaff();
	int nbDSB();
	int nbGenerations();
	vector<vector<int>> nbfailedmeiosis();
	bool zygosity();
	int everygen();
	bool ismigration();
	double q();
	bool withDSB();
	double w();
	string name();
	map<int,double> Ageallele();
	map<int,vector<double>> infoperallele();
	//============================
	//           Setters
	//============================

	//============================
	//           Methods
	//============================
	int choose(int n);
	int bernoulli_draw(double p);
	int binomial_draw(int n, double p);  
	double chosegamma(double meanaff, double varaff);
	vector<int> choosemany(int k, vector<int> vect);
	vector<int> vectfreesites(vector<int> vect, int nb);
	vector<vector<int>> occupiedsites(vector<int> vect);
	void sitemutation();
	void allelemutation();
	void updatemissingallele();
	void printpop(int n);
	void printgen(int n);
	void printposallele();
	void printallelepos();
	void printaffinity();
	int Meiosis(int no_chrom_ind, int nb_gen);
	void fillnewpop(int nb_gen);
	void manygenerations();
	int get_allele_number();
	double freqallele(int allelename);
	double get_current_diversity();
	double activitymoyallele(int allele);
	double get_current_activity();
	void migration();
	vector<double> freqneutral();
	double freqall(int allele);
	double actall(int allele);
	void printageallele();
	double get_age_allele(int allname);
	void printinfoallele();
	double get_info_allele(int allname);
	
	protected:
	//============================
	//       Data members
	//============================
	vector<vector<vector<int>>> populations_; // marix of the populations current and next: columns => number of chromosomes (2*N), rows => all sites in the genome (L)
	vector<vector<vector<int>>> populations1_; //second population for the migration
	vector<vector<int>> genotypes_; // prdm9 allele for each chromosom
	vector<vector<int>> genotypes1_; //genotype of the second pop if migration
	map<int,vector<int>> Siteforeacheallele_; // all sites positions for each allele
	vector<int> Alleleforeachpos_; // prdm9 allele corresponding to each position
	vector<double> Affinity_; //sites affinity
	int parityIndex_; // index corresponding to the current population
	int N_; // effective population size
	int L_; // size of the genome
	int indPrdm9_; // index of the site of the prdm9 gene
	int nballele_; // number of allele which possess active motifs in the popluation (nb of alleles in the pop since the begining) -> nb current allele???
	int nbsite_; // number of sites for each new allele
	double u_; //prdm9 mutation rate
	double v_; // target mutation rate
	double m_; // migration rate
	double meanaff_; //mean gamma law 
	double varaff_; // variance gamma law
	int nbDSB_; //nb of DSB
	int nbGenerations_; //nb of generations
	vector<vector<int>> nbfailedmeiosis_;
	bool zygosity_;
	int everygen_;
	bool ismigration_;
	double q_;
	bool withDSB_;
	double w_;
	string name_;
	map<int,double> Ageallele_;
	map<int,vector<double>> infoperallele_;
};
