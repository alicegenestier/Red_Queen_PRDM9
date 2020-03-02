
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
	Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen,double m,double alpha,double beta,int nbgenmig,bool popsamesize,string name);
	
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
	double alpha();
	double beta();
	int nbgenmig();
	bool popsamesize();
	//============================
	//           Setters
	//============================

	//============================
	//           Methods
	//============================
	int choose(int n);
	int bernoulli_draw(double p);
	int binomial_draw(int n, double p);  
	double choosegamma(double meanaff, double varaff);
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
	vector<int> choosemanymigration(int k);
	double choosebeta(double alpha, double beta);
	double q_two_hap(vector<int> haplotype1, vector<int> haplotype2);
	double q_hom(int allele, vector<int> haplotype1, vector<int> haplotype2);
	double q_hete(int allele1, int allele2, vector<int> haplotype1, vector<int> haplotype2);
	
	protected:
	//============================
	//       Data members
	//============================
	vector<vector<vector<int>>> populations_; // marix of the populations current and next: columns => number of chromosomes (2*N), rows => all sites in the genome (L)
	vector<vector<vector<int>>> populations1_; //first population for the migration
	vector<vector<vector<int>>> populations2_; //second population for the migration
	vector<vector<int>> genotypes_; // prdm9 allele for each chromosom
	vector<vector<int>> genotypes1_; //first genotype for the migration
	vector<vector<int>> genotypes2_; //second genotype for the migration
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
	vector<vector<int>> nbfailedmeiosis_; //store the number of failed meiosis of each type
	bool zygosity_; //make difference between heterozygots and homozygots
	int everygen_; //nb of generation at which we want to print the reults in the files
	bool ismigration_; //is there migration
	double q_; //
	bool withDSB_; //do we take into account the 2 DSB at one site as a cause of failed meiosis
	double w_; //neutral site mutation rate
	string name_; //name of the files
	map<int,double> Ageallele_; //store the age of each allele
	map<int,vector<double>> infoperallele_; //store information for each allele such as the numer of symetrical binding or the nb of failed meiosis per allele
	double alpha_; //first param of the beta distribution
	double beta_; //second param of the beta distribution
	int nbgenmig_; //nb of the generation at which we want to split de pop for migration (if = 0 => begin directly with 2 pop)
	bool popsamesize_; //two pop for migration has the same size of the initial pop or devided by two
};
