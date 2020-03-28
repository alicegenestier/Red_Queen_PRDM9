
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
	Model(int N,int L,int nbsite,int indPrdm9,int nballele,int parityIndex,double v,double u,double w,double meanaff,double varaff,int nbDSB,int nbGenerations,bool ismigration,bool zygosity,bool withDSB,int everygen,double m,double alpha,double beta,int nbgenmig,bool popsamesize,int nbloop,string name);
	
	//============================
	//        Destructors
	//============================
	
	//============================
	//           Getters
	//============================
	vector<vector<vector<int>>> populations(); 
	vector<vector<vector<int>>> populations1();
	vector<vector<vector<int>>> populations2();
	vector<vector<int>> genotypes(); 
	vector<vector<int>> genotypes1();
	vector<vector<int>> genotypes2();
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
	vector<vector<int>> nbfailedmeiosis1();
	vector<vector<int>> nbfailedmeiosis2();
	bool zygosity();
	int everygen();
	bool ismigration();
	double q();
	double q1();
	double q2();
	bool withDSB();
	double w();
	string name();
	map<int,double> Ageallele();
	map<int,double> Ageallele1();
	map<int,double> Ageallele2();
	map<int,vector<double>> infoperallele();
	map<int,vector<double>> infoperallele1();
	map<int,vector<double>> infoperallele2();
	double alpha();
	double beta();
	int nbgenmig();
	bool popsamesize();
	int nbloop();
	//============================
	//           Setters
	//============================

	//============================
	//           Methods
	//============================
	int choose(int n); //choose a nb between 0 and n-1 with a uniform law
	int bernoulli_draw(double p); //give bernoulli distrib with prob p
	int binomial_draw(int n, double p); // 
	double choosegamma(double meanaff, double varaff); //
	vector<int> choosemany(int k, vector<int> vect); //
	vector<int> vectfreesites(vector<int> vect, int nb); //
	vector<vector<int>> occupiedsites(vector<int> vect, vector<vector<int>>* genotype); //
	void sitemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype);
	void allelemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,double>* Ageallele, map<int,vector<double>>* infoperallele); //
	void updatemissingallele(); //
	void printpop(int n, vector<vector<vector<int>>> population); //
	void printgen(int n, vector<vector<int>> genotype); //
	void printposallele(); //
	void printallelepos(); //
	void printaffinity(); //
	int Meiosis(int no_chrom_ind, int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, vector<vector<int>>* nbfailedmeiosis, double* q); //
	void fillnewpop(int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, vector<vector<int>>* nbfailedmeiosis, double* q); //
	void manygenerations(); //
	vector<int> get_allele_number(vector<vector<vector<int>>*> vectgen); //
	double freqallele(int allelename, vector<vector<int>>* genotype); //
	double get_current_diversity(vector<vector<int>>* genotype); //
	double activitymoyallele(int allele,  vector<vector<vector<int>>>* population); //
	double get_current_activity(vector<vector<int>>* genotype, vector<vector<vector<int>>>* population); //
	void migration(); //
	vector<double> freqneutral(vector<vector<vector<int>>>* population); //
	double freqall(int allele, vector<vector<int>>* genotype, vector<vector<vector<int>>>* population); //
	double actall(int allele, vector<vector<vector<int>>>* population); //
	void printageallele(map<int,double>* Ageallele); //
	double get_age_allele(int allname, map<int,double>* Ageallele); //
	void printinfoallele(map<int,vector<double>>* infoperallele); //
	vector<double> get_info_allele(int allname, map<int,vector<double>>* infoperallele); //
	vector<int> choosemanymigration(int k); //
	double choosebeta(double alpha, double beta); //
	double q_two_hap(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom); //
	//double q_hom(int allele, vector<int> haplotype1, vector<int> haplotype2); //
	//double q_hete(int allele1, int allele2, vector<int> haplotype1, vector<int> haplotype2); //
	double get_q(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype);
	double get_q_hybrid(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen);
	vector<int> get_one_gamete(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom);
	vector<double> get_q_fertility_indep(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen, int nbloop_, int nopop);
	double get_FST_neutral(vector<vector<vector<vector<int>>>*> vectpop);
	double get_FST_PRDM9(vector<vector<vector<int>>*> vectgen);
	
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
	int nballele_; // nb of alleles in the pop since the begining
	int nbsite_; // number of sites for each new allele
	double u_; //prdm9 mutation rate
	double v_; // target mutation rate
	double m_; // migration rate
	double meanaff_; //mean gamma law 
	double varaff_; // variance gamma law
	int nbDSB_; //nb of DSB
	int nbGenerations_; //nb of generations
	vector<vector<int>> nbfailedmeiosis_; //store the number of failed meiosis of each type
	vector<vector<int>> nbfailedmeiosis1_; //
	vector<vector<int>> nbfailedmeiosis2_; //
	bool zygosity_; //make difference between heterozygots and homozygots
	int everygen_; //nb of generation at which we want to print the results in the files
	bool ismigration_; //is there migration
	double q_; //
	double q1_;//
	double q2_;//
	bool withDSB_; //do we take into account the 2 DSB at one site as a cause of failed meiosis
	double w_; //neutral site mutation rate
	string name_; //name of the files
	map<int,double> Ageallele_; //store the age of each allele
	map<int,double> Ageallele1_; //
	map<int,double> Ageallele2_; //
	map<int,vector<double>> infoperallele_; //store information for each allele such as the numer of symetrical binding or the nb of failed meiosis per allele
	map<int,vector<double>> infoperallele1_; // totnbfail, 2dsb, nodsb, nosym, q, totnbok
	map<int,vector<double>> infoperallele2_; //
	double alpha_; //first param of the beta distribution
	double beta_; //second param of the beta distribution
	int nbgenmig_; //nb of the generation at which we want to split de pop for migration (if = 0 => begin directly with 2 pop)
	bool popsamesize_; //two pop for migration has the same size of the initial pop or devided by two
	int nbloop_; //nb loop for indep q and fertility
};
