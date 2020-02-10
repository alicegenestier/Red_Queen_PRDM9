
//============================
//          Includes
//============================
#include <vector>
#include <map>
#include <random>
#include <fstream>

using namespace std;

class Model {
	public:
	//============================
	//        Constructors
	//============================
	Model();
	
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
	float u();
	float v();
	float m();
	double meanaff(); 
	double varaff();
	float nbDSB();
	int nbGenerations();
	vector<vector<int>> nbfailedmeiosis();
	bool zygosity();
	int everygen();
	bool ismigration();
	double q();
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
	vector<int> occupiedsites(vector<int> vect);
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
	float u_; //prdm9 mutation rate
	float v_; // target mutation rate
	float m_; // migration rate
	double meanaff_; //mean gamma law 
	double varaff_; // variance gamma law
	float nbDSB_; //nb of DSB
	int nbGenerations_; //nb of generations
	vector<vector<int>> nbfailedmeiosis_;
	bool zygosity_;
	int everygen_;
	bool ismigration_;
	double q_;
};
