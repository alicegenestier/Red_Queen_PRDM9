
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
	Model(int N, int L, int nbsite, int indPrdm9, int nballele, int parityIndex, double v, double u, double w, double meanaff, double varaff, int nbDSB, int nbGenerations, bool ismigration, bool zygosity, bool withDSB, int everygen, double m, double alpha, double beta, int nbgenmig, bool popsamesize, int nbhyb, int nbcore, bool isallele, bool issampling, bool isanalytic, double ctot, bool targetcomp, int quantilenb, int nbmeiperind, double cfreethreshold, bool affinityUniform, double dosagecoeff, string name);
	
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
	bool targetcomp();
	int everygen();
	bool ismigration();
	double q();
	double q1();
	double q2();
	double qsym();
	double qsym1();
	double qsym2();
	bool withDSB();
	double w();
	string name();
	map<int,double> Ageallele();
	map<int,double> Ageallele1();
	map<int,double> Ageallele2();
	map<int,vector<double>> infoperallele();
	map<int,vector<double>> infoperallele1();
	map<int,vector<double>> infoperallele2();
	map<int,vector<double>> infoperallele_hom();
	map<int,vector<double>> infoperallele1_hom();//////////////////
	map<int,vector<double>> infoperallele2_hom();//////////////////
	map<int,vector<double>> infoperallele_het();
	map<int,vector<double>> infoperallele1_het();//////////////////
	map<int,vector<double>> infoperallele2_het();//////////////////
	double alpha();
	double beta();
	int nbgenmig();
	bool popsamesize();
	int nbhyb();
	int nbcore();
	bool isallele();
	bool issampling();
	bool isanalytic();
	double qnum();
	double qdenom();
	double qnum1();
	double qnum2();
	double qdenom1();
	double qdenom2();
	double ctot();
	int quantilenb();
	int nbmeiperind(); 
	double dosagecoeff();
	//============================
	//           Setters
	//============================

	//============================
	//           Methods
	//============================
	int choose(int n); //choose a nb between 0 and n-1 with a uniform law
	int bernoulli_draw(double p); //give bernoulli distrib with prob p
	int binomial_draw(int n, double p); // binomial law
	double choosegamma(double meanaff, double varaff); // gamma law
	vector<int> choosemany(int k, vector<int> vect); //return the vector of index of the chosen positions
	vector<int> vectfreesites(vector<int> vect, int nb); //return the index, in the vector vect, of all positions equal to nb
	vector<vector<int>> occupiedsites(vector<int> vect, vector<vector<int>>* genotype); //return the index of all occupied positions
	void sitemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype); // for each position with allele, probability v to mutate and if mutation, choose randomly 1 chrom to mutate
	void allelemutation(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,double>* Ageallele, vector<int>vectnbind, double num_pop); // for each chromosome, mutation of allele with probability u
	void updatemissingallele(); //update the map if one allele has gone extinct : allele suppressed in map, all positions of this allele are set to -1 (free position) and all the activations are set to one
	void printpop(int n, vector<vector<vector<int>>> population); //print popuation
	void printgen(int n, vector<vector<int>> genotype); //print genotype vector
	void printposallele(); //print site positions for each allele
	void printallelepos(); //print Allele for each position
	void printaffinity(); //print Affinity
	int Meiosis(int no_chrom_ind, int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele,map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het, vector<vector<int>>* nbfailedmeiosis, double* q, double* qsym, double* qnum, double* qdenom, int indiv, int nb_meiosis); //Perform the meiosis of one individual and return 0 if the meiosis fails somewhere or 1 if the meiosis is successfull
	void fillnewpop(int nb_gen, vector<vector<vector<int>>>* population, vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele, map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het, vector<vector<int>>* nbfailedmeiosis, double* q, double* qsym, double* qnum, double* qdenom); //methode that fill new population (make as many meiosis as it is necessary to fill the entire new population)
	void manygenerations(); //methode which mix together and runs all the previous functions during X generations 
	vector<int> get_allele_number(vector<vector<vector<int>>*> vectgen, bool nbtot); //get number of different alleles in the population
	double freqallele(int allelename, vector<vector<int>>* genotype); //give the frequence of each allele in the population
	double get_current_diversity(vector<vector<int>>* genotype); //get the PRDM9 diversity
	double activitymoyallele(int allele,  vector<vector<vector<int>>>* population); //give the mean activity of a given allele
	double get_current_activity(vector<vector<int>>* genotype, vector<vector<vector<int>>>* population); //give the average activity of an allele in the population
	void migration(); //give the list of k index individual in the population that will migrate from one generation to the other
	vector<double> freqneutral(vector<vector<vector<int>>>* population); //give the frequency and the activity of neutral sites
	double freqall(int allele, vector<vector<int>>* genotype, vector<vector<vector<int>>>* population); //give the frequency of an allele or neural site
	double actall(int allele, vector<vector<vector<int>>>* population); //give the activity of an allele or neutral site
	void printageallele(map<int,double>* Ageallele); //print age allele
	double get_age_allele(int allname, map<int,double>* Ageallele); //get the age of each allele segregating in the population
	void printinfoallele(map<int,vector<double>>* infoperallele); //print info per allele
	vector<double> get_info_allele(int allname, map<int,vector<double>>* infoperallele, int index_meiosis_success); //get the symmetrical binding rate and the fertility rate per allele and individual
	vector<int> choosemanymigration(int k); //give the list of k index individual in the population that will migrate from one generation to the other
	double choosebeta(double alpha, double beta); //beta law
	double q_two_hap(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom); //give q of a specific individual
	double get_q(vector<vector<vector<int>>>* population, vector<vector<int>>* genotype);//give q indepedently from the system (sampling) for the next generation
	double get_q_hybrid(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen);//give q of hybride
	vector<int> get_one_gamete(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom);//give one gamete of an individual
	vector<double> get_q_fertility_indep(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen, int hyb_, int nopop);//give the mean probability of symetrical site and the mean of fertility in the population
	double get_FST_neutral(vector<vector<vector<vector<int>>>*> vectpop);//give the FST of neutral site
	double get_FST_PRDM9(vector<vector<vector<int>>*> vectgen);//give the PRDM9 FST
	double get_mean_age(vector<vector<int>>* genotype, map<int,double>* Ageallele);//get the mean age of the alleles segregating in the population
	double get_mean_sigma(vector<vector<int>>* genotype, vector<map<int,vector<double>>> analytic_indiv);//get the mean sigma of the alleles segregating in the population
	vector<double> get_mean_coccup_cfree(vector<vector<int>>* genotype, map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het);//get the mean number of bound protein of the allele segregating in the population
	double get_mean_affinity(int allele, vector<vector<vector<int>>>* pop);//give the mean affinity of all the active sites of a given allele
	double puissance_double(int puiss, double x);//Power function for double
	int puissance_int(int puiss, int x);//Power function for int
	vector<map<int,vector<double>>> q_fert_individual_analytique(vector<vector<int>>* genotype, vector<vector<vector<int>>>* pop, map<int,vector<double>>* infoperallele_hom, map<int,vector<double>>* infoperallele_het);//Give the mean q, fertility and sigma for each allele present in the population
	double Mean_fert_new_allele(vector<vector<int>>* genotype, vector<vector<vector<int>>>* pop, map<int,vector<double>>* infoperallele_het);//Give the mean fertility of a new allele in all possible heterozygot context in the population
	vector<double> sigma_q_w_0();//Give the mean q, fertility, sigma and cfree/ctot for a new allele : in parameters file
	double if_allele_print_else_nan(map<int,vector<double>> map_allele, int allele_nb, int third_int);/////////////////////////
	///////////////////////////////////////////
	vector<map<int,vector<double>>> q_fert_hybrid_analytic_general(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen, vector<map<int,vector<double>>*> vectinfo_hom, vector<map<int,vector<double>>*> vectinfo_het);
	void get_q_fert_hybrid_analytic(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen, vector<map<int,vector<double>>*> vectinfo_hom, vector<map<int,vector<double>>*> vectinfo_het, map<int,vector<double>>* res_q, map<int,vector<double>>* res_fert, map<int,vector<double>>* res_number_meiosis);
	void q_fert_two_hap_analytic(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom, vector<map<int,vector<double>>*> vectinfo_hom, vector<map<int,vector<double>>*> vectinfo_het, map<int,vector<double>>* res_q, map<int,vector<double>>* res_fert, map<int,vector<double>>* res_number_meiosis);
	///////////////////////////////////////////
	
	////////////////////////////////////////////
	/*vector<vector<double>> get_q_fertility_indep_bis(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen, int nbhyb_);
	vector<double> get_q_hybrid_bis(vector<vector<vector<vector<int>>>*> vectpop, vector<vector<vector<int>>*> vectgen);
	vector<vector<double>> get_one_gamete_bis(vector<int> genotype_indiv, vector<vector<int>> indiv_chrom);*/
	//////////////////////////////////////////
	
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
	vector<vector<int>> nbfailedmeiosis1_; //store the number of failed meiosis of each type for pop 1
	vector<vector<int>> nbfailedmeiosis2_; //store the number of failed meiosis of each type for pop 2
	bool zygosity_; //make difference between heterozygots and homozygots
	bool targetcomp_; //If we want to take into account target competition
	int everygen_; //nb of generation at which we want to print the results in the files
	bool ismigration_; //is there migration
	double q_; //probability of symmetrical binding + DSB
	double q1_;//probability of symmetrical binding + DSB for pop 1
	double q2_;//probability of symmetrical binding + DSB for pop 2
	double qsym_;//
	double qsym1_;//
	double qsym2_;//
	bool withDSB_; //do we take into account the 2 DSB at one site as a cause of failed meiosis
	double w_; //neutral site mutation rate
	string name_; //name of the files
	map<int,double> Ageallele_; //store the age of each allele
	map<int,double> Ageallele1_; //store the age of each allele for pop 1
	map<int,double> Ageallele2_; //store the age of each allele for pop 2
	map<int,vector<double>> infoperallele_; //store information for each allele such as the numer of symetrical binding or the nb of failed meiosis per allele
	map<int,vector<double>> infoperallele_hom_;//store information for each allele in homozygous state such as the numer of symetrical binding or the nb of failed meiosis per allele
	map<int,vector<double>> infoperallele_het_;//store information for each allele in heterozygous state such as the numer of symetrical binding or the nb of failed meiosis per allele
	map<int,vector<double>> infoperallele1_; // store information for each allele such as the numer of symetrical binding or the nb of failed meiosis per allele for pop 1
	map<int,vector<double>> infoperallele2_; // store information for each allele such as the numer of symetrical binding or the nb of failed meiosis per allele for pop 2
	map<int,vector<double>> infoperallele1_hom_; //////////////////
	map<int,vector<double>> infoperallele2_hom_; //////////////////
	map<int,vector<double>> infoperallele1_het_; //////////////////
	map<int,vector<double>> infoperallele2_het_; //////////////////
	double alpha_; //first param of the beta distribution
	double beta_; //second param of the beta distribution
	int nbgenmig_; //nb of the generation at which we want to split de pop for migration (if = 0 => begin directly with 2 pop)
	bool popsamesize_; //two pop for migration has the same size of the initial pop or devided by two
	int nbhyb_; //nb loop for indep q and fertility
	int nbcore_;//nb core for open mp
	bool isallele_;//do we want the result for each allele
	bool issampling_;//do we want the result from the sampling
	bool isanalytic_;//do we want the analytical results
	double qdenom_; //denominator of q
	double qdenom1_; //denominator of q
	double qdenom2_; //denominator of q
	double qnum_; //numerator of q
	double qnum1_; //numerator of q
	double qnum2_; //numerator of q
	double ctot_;//total concentration for 1 PRDM9 allele in heterozygot
	int quantilenb_;//number of categories for the affinity distribution
	map<double,vector<double>> nbsitesperquantile_;//[quantile category]{mean of the quantile category, number of active site in each category of quantile}
	int nbmeiperind_;//number of meiosis that an individual can perform before being caracterizes as steril
	double cfreethreshold_;//threshold for the determination of cfree
	bool affinityUniform_;//affinity distribution follows a uniform law
	double dosagecoeff_;
};
