#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*definition(structures)*/
typedef struct {
	double *coef;
	int *power;
	int size;
} polynomial;

typedef struct {
	polynomial poly1;
	polynomial poly2;
	polynomial poly3;
	polynomial poly4;
} polynomial_individual;

typedef struct {
	polynomial poly1;
	polynomial poly2;
	polynomial poly3;
	polynomial poly4;
	polynomial poly5;
	polynomial poly6;
	polynomial poly7;
	polynomial poly8;
	polynomial poly9;
} polynomial_sibling;


/*definition(variables)*/
static int **pedigree=NULL;
static int **nuclear=NULL;
static int NIndividual;
static int N_NUCLEAR;
static int MAX_SIBLING;
static int *ID_branch_candidate;
static double fA;
static int MAX_POWER;
static polynomial_individual *likelihood;
static polynomial_individual *likelihood_conditional;
static polynomial_sibling *likelihood_sibling;
static int ID;
static int ROW;


/*definition(functions)*/
void allocate_memory_polynomial(polynomial *poly);
void allocate_memory_polynomial_individual(polynomial_individual *poly);
void allocate_memory_polynomial_individual_all(polynomial_individual *poly);
void allocate_memory_polynomial_sibling(polynomial_sibling *poly);
void free_likelihood(int *ids);
void free_likelihood_conditioanl(void);
void free_likelihood_sibling(int *ids);
void free_polynomial(polynomial *poly);
void free_polynomial_individual(polynomial_individual *poly);
void free_pedigree_nuclear(void);
void make_nuclear(void);
int is_terminal(int id);
int is_founder(int id);
int get_id_ith(int id);

void search_terminal(int id, int *terminals);
void sibling(int id, int *siblingID);
void descendant(int id, int *resultID);
void search_spouse(int id, int *spouseID);
void search_children(int id, int *childrenID);
void search_ancesters(int id, int *ancestersID);
void search_parent(int id, int *parentID);
void search_descendants(int id, int *descendantsID);
void search_branch2(int id, int *branchID);
void search_branch_in(int id, int *branchID);
void search_branch(int id, int *branchID);
void branch(int *terminalID, int *resultID);

void remove_0terms(polynomial *poly);
void eval_min(int power1, int power2, int power3, int power4, int *min, int *order);
void sum(polynomial_individual *poly, polynomial *result);
void prod_scalar(double a, polynomial *poly, polynomial *result);
void prod_polynomial(polynomial *poly1, polynomial *poly2, polynomial *res);
void set_likelihood_conditional(polynomial_individual *poly, int row);
void initialise_likelihood0(polynomial *poly);
double InheritP(int genotype_parent, int origin_id, int allele_id);
void penetrance(int Phenotype, int row, polynomial *poly);
void Xi_cal_founderP_3paras(polynomial_individual *poly, int Phenotype);

void AncestralLikelihoodDiplo_C(int id_ith, int *siblingID);
void SumOfSiblingLikelihood(int row_f, int row_m, int id, polynomial *poly);
void substitute_polynomial(polynomial *to, polynomial *from);
void calculate_likelihood_sibling(polynomial_sibling *poly, int *siblingID);
void AncestralLikelihood_C(int id);
void ConditionalLikelihoodDiplo(int id, int row);
void clkh(int *ids);
void AncestralLikelihoodDiplo(int id_ith, int *siblingID);
void AncestralLikelihood(int id);




void mlep(int *N_ind, int *ped, double *freqA, int *maxpower, double *coef, int *power, int *size){


	NIndividual = *N_ind;
	fA = *freqA;
	MAX_POWER = *maxpower;

	int i, j;
	
	pedigree = (int **)malloc(sizeof(int *)*NIndividual);

//	if( pedigree == NULL ){
//		printf("can not allocate memory");
//		exit(1);
//	}

	for( i=0; i<NIndividual; i++ ){
		pedigree[i] = (int *)malloc(sizeof(int *)*5);
//		if( pedigree[i] == NULL ){
//			printf("can not allocate memory");
//			exit(1);
//		}
		for( j=0; j<5; j++ ){
			pedigree[i][j] = *(ped+NIndividual*j+i);
		}
	}

	make_nuclear();

	int *terminals = (int *)calloc(NIndividual, sizeof(int));
	search_terminal(nuclear[1][0], terminals);

	ID_branch_candidate = (int *)calloc(NIndividual, sizeof(int));

	int *branchesID;
	branchesID = (int *)calloc(NIndividual, sizeof(int));
	branch(terminals, branchesID);

	likelihood = (polynomial_individual *)malloc(sizeof(polynomial_individual)*NIndividual);
	allocate_memory_polynomial_individual_all(likelihood);
	likelihood_sibling = (polynomial_sibling *)malloc(sizeof(polynomial_sibling)*NIndividual);
	allocate_memory_polynomial_sibling(likelihood_sibling);

	clkh(branchesID);
	AncestralLikelihood(get_id_ith(*terminals));

	polynomial *lkh;
	lkh = (polynomial *)malloc(sizeof(polynomial));
	allocate_memory_polynomial(lkh);

	sum(&likelihood[get_id_ith(*terminals)],lkh);

	for( i=0; i<(*lkh).size; i++ ){
		*(coef+i) = (*lkh).coef[i];
		*(power+i) = (*lkh).power[i];
	}
	*size = (*lkh).size;

	free(terminals);
	free(ID_branch_candidate);
	free(branchesID);
	free(likelihood_sibling);
	free_pedigree_nuclear();
}


void allocate_memory_polynomial(polynomial *poly){
	
	(*poly).coef = (double *)malloc(sizeof(double));
	(*poly).power = (int *)malloc(sizeof(int));
	(*poly).size = 0;

}


void allocate_memory_polynomial_individual(polynomial_individual *poly){
	
	(*poly).poly1.coef = (double *)malloc(sizeof(double));
	(*poly).poly1.power = (int *)malloc(sizeof(int));
	(*poly).poly1.size = 0;
	(*poly).poly2.coef = (double *)malloc(sizeof(double));
	(*poly).poly2.power = (int *)malloc(sizeof(int));
	(*poly).poly2.size = 0;
	(*poly).poly3.coef = (double *)malloc(sizeof(double));
	(*poly).poly3.power = (int *)malloc(sizeof(int));
	(*poly).poly3.size = 0;
	(*poly).poly4.coef = (double *)malloc(sizeof(double));
	(*poly).poly4.power = (int *)malloc(sizeof(int));
	(*poly).poly4.size = 0;

}


void allocate_memory_polynomial_individual_all(polynomial_individual *poly){
	
	int i;
	for( i=0; i<NIndividual; i++ ){
		(*(poly+i)).poly1.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly1.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly1.size = 0;
		(*(poly+i)).poly2.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly2.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly2.size = 0;
		(*(poly+i)).poly3.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly3.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly3.size = 0;
		(*(poly+i)).poly4.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly4.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly4.size = 0;
	}
}


void allocate_memory_polynomial_sibling(polynomial_sibling *poly){
	
	int i;
	for( i=0; i<NIndividual; i++ ){
		(*(poly+i)).poly1.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly1.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly1.size = 0;
		(*(poly+i)).poly2.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly2.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly2.size = 0;
		(*(poly+i)).poly3.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly3.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly3.size = 0;
		(*(poly+i)).poly4.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly4.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly4.size = 0;
		(*(poly+i)).poly5.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly5.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly5.size = 0;
		(*(poly+i)).poly6.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly6.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly6.size = 0;
		(*(poly+i)).poly7.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly7.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly7.size = 0;
		(*(poly+i)).poly8.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly8.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly8.size = 0;
		(*(poly+i)).poly9.coef = (double *)malloc(sizeof(double));
		(*(poly+i)).poly9.power = (int *)malloc(sizeof(int));
		(*(poly+i)).poly9.size = 0;
	}
}


void free_likelihood(int *ids){
	
	int i=0;
	int id_ith;

	while( *(ids+i) != 0 ){
		id_ith = get_id_ith(*(ids+i));
		if( likelihood[id_ith].poly1.size != 0 ){
			free(likelihood[id_ith].poly1.coef);
			free(likelihood[id_ith].poly1.power);
			free(likelihood[id_ith].poly2.coef);
			free(likelihood[id_ith].poly2.power);
			free(likelihood[id_ith].poly4.coef);
			free(likelihood[id_ith].poly4.power);
			likelihood[id_ith].poly1.size=0;
			likelihood[id_ith].poly2.size=0;
			likelihood[id_ith].poly4.size=0;
		}
		//For branchID
		if( likelihood[id_ith].poly3.size != 0 ){
			free(likelihood[id_ith].poly3.coef);
			free(likelihood[id_ith].poly3.power);
			likelihood[id_ith].poly3.size=0;
		}
		i++;
	}
}


void free_likelihood_conditional(void){
	
	int i;
	for( i=0; i<NIndividual; i++ ){
		free(likelihood_conditional[i].poly1.coef);
		free(likelihood_conditional[i].poly1.power);
		free(likelihood_conditional[i].poly2.coef);
		free(likelihood_conditional[i].poly2.power);
		free(likelihood_conditional[i].poly3.coef);
		free(likelihood_conditional[i].poly3.power);
		free(likelihood_conditional[i].poly4.coef);
		free(likelihood_conditional[i].poly4.power);
	}
	free(likelihood_conditional);
}


void free_likelihood_sibling(int *ids){
	
	int i=0;
	int id_ith;

	while( *(ids+i) != 0 ){
		id_ith = get_id_ith(*(ids+i));
		if( likelihood_sibling[id_ith].poly1.size != 0 ){
			free(likelihood_sibling[id_ith].poly1.coef);
			free(likelihood_sibling[id_ith].poly1.power);
			free(likelihood_sibling[id_ith].poly2.coef);
			free(likelihood_sibling[id_ith].poly2.power);
			free(likelihood_sibling[id_ith].poly3.coef);
			free(likelihood_sibling[id_ith].poly3.power);
			free(likelihood_sibling[id_ith].poly4.coef);
			free(likelihood_sibling[id_ith].poly4.power);
			free(likelihood_sibling[id_ith].poly5.coef);
			free(likelihood_sibling[id_ith].poly5.power);	
			free(likelihood_sibling[id_ith].poly6.coef);
			free(likelihood_sibling[id_ith].poly6.power);
			free(likelihood_sibling[id_ith].poly7.coef);
			free(likelihood_sibling[id_ith].poly7.power);
			free(likelihood_sibling[id_ith].poly8.coef);
			free(likelihood_sibling[id_ith].poly8.power);
			free(likelihood_sibling[id_ith].poly9.coef);
			free(likelihood_sibling[id_ith].poly9.power);
			likelihood_sibling[id_ith].poly1.size=0;
			likelihood_sibling[id_ith].poly2.size=0;
			likelihood_sibling[id_ith].poly3.size=0;
			likelihood_sibling[id_ith].poly4.size=0;
			likelihood_sibling[id_ith].poly5.size=0;
			likelihood_sibling[id_ith].poly6.size=0;
			likelihood_sibling[id_ith].poly7.size=0;
			likelihood_sibling[id_ith].poly8.size=0;
			likelihood_sibling[id_ith].poly9.size=0;
		}
		i++;
	}
}


void free_polynomial(polynomial *poly){
	
	free((*poly).coef);
	free((*poly).power);
	free(poly);
}


void free_polynomial_individual(polynomial_individual *poly){
	
	free((*poly).poly1.coef);
	free((*poly).poly1.power);
	free((*poly).poly2.coef);
	free((*poly).poly2.power);
	free((*poly).poly3.coef);
	free((*poly).poly3.power);
	free((*poly).poly4.coef);
	free((*poly).poly4.power);
	free(poly);
}


void free_pedigree_nuclear(void){
	
	int i;
	
	for( i=0; i<NIndividual; i++ ){
		free(pedigree[i]);
	}
	for( i=0; i<N_NUCLEAR; i++ ){
		free(nuclear[i]);
	}
	free(pedigree);
	free(nuclear);
}


void make_nuclear(void){

	N_NUCLEAR = 1;
	MAX_SIBLING = 1;
	
	int i, j, k;
	int fid, mid;
	int nfounder;
	int flag;
	int **tmp;

	tmp = (int **)malloc(sizeof(int *) * NIndividual);
//	if( tmp == NULL ){
//		printf("can not allocate memory");
//		exit(1);
//	}

	tmp[0] = (int *)calloc(NIndividual, sizeof(int *));
//	if( tmp[0] == NULL ){
//		printf("can not allocate memory");
//		exit(1);
//	}
	
	nfounder = 2;
	for( i=0; i<NIndividual; i++ ){
		/*case: founder*/
		if( pedigree[i][1] == 0 ){
			tmp[0][nfounder] = pedigree[i][0];
			if( MAX_SIBLING < nfounder-1 ){
				MAX_SIBLING = nfounder-1;
			}
			nfounder++;
		}
		/*case non-founder*/
		else{
			fid = pedigree[i][1];
			mid = pedigree[i][2];
			flag = 0;
			if( N_NUCLEAR > 1 ){
				for( j=1; j<N_NUCLEAR; j++ ){
					if( tmp[j][0]==fid && tmp[j][1]==mid ){
						k = 2;
						while( tmp[j][k]!=0 ){
							k++;
						}
						tmp[j][k] = pedigree[i][0];
						if(MAX_SIBLING < k-1){
							MAX_SIBLING = k-1;
						}
						flag = 1;
						break;
					}
				}
			}
			if( flag == 0 ){
				tmp[N_NUCLEAR] = (int *)calloc(NIndividual, sizeof(int *));
//				if( tmp[N_NUCLEAR] == NULL ){
//					printf("can not allocate memory");
//					exit(1);
//				}
				tmp[N_NUCLEAR][0] = fid;
				tmp[N_NUCLEAR][1] = mid;
				tmp[N_NUCLEAR][2] = pedigree[i][0];
				N_NUCLEAR++;
			}
		}
	}
	
	nuclear = (int **)malloc(sizeof(int *) * N_NUCLEAR);
//	if( nuclear == NULL ){
//		printf("can not allocate memory");
//		exit(1);
//	}
	for( i=0; i<N_NUCLEAR; i++ ){
		nuclear[i] = (int *)calloc((MAX_SIBLING+2), sizeof(int *));
//		if( nuclear[i] == NULL ){
//			printf("can not allocate memory");
//			exit(1);
//		}
		for( j=0; j<MAX_SIBLING+2; j++ ){
			nuclear[i][j] = tmp[i][j];
		}
		free(tmp[i]);
	}
	free(tmp);
}


int is_terminal(int id){

	int i;
	for( i=0; i<N_NUCLEAR; i++ ){
		if( nuclear[i][0] == id || nuclear[i][1] == id ){
			return 0;
		}
	}
	return 1;
}


int is_founder(int id){

	int i;
	for( i=2; i<MAX_SIBLING+2; i++ ){
		if( nuclear[0][i] == id ){
			return 1;
		}
	}
	return 0;
}


int get_id_ith(int id){

	int i;
	for ( i=0; i<NIndividual; i++ ){
		if( pedigree[i][0] == id ){
			return i;
		}
	}
	return -1;
}


void search_terminal(int id, int *terminals){

	int i,j,k,l;
	int flag;
	int terminals_index = 0;

	for( i=1; i<N_NUCLEAR; i++ ){
		if( nuclear[i][0] == id || nuclear[i][1] == id ){
			j = 2;
			while( nuclear[i][j] != 0 ){
				j++;
			}

			flag = 0;
			for( k=2; k<j; k++ ){
				if( is_terminal(nuclear[i][k]) == 0 ){
					int *tmp = (int *)calloc(NIndividual, sizeof(int));
					search_terminal(nuclear[i][k], tmp);
					l = 0;
					while( *(tmp+l) != 0 ){
						*(terminals+terminals_index) = *(tmp+l);
						terminals_index++;
						l++;
					}
					free(tmp);
					flag = 1;
				}
			}
			if( flag == 0 ){
				*terminals = nuclear[i][2];
			}
		}
	}
}


void sibling(int id, int *siblingID){

	int i,j;
	int flag;
	int index;

	for( i=1; i<N_NUCLEAR; i++ ){
		j = 2;
		flag = 0;
		while( nuclear[i][j] != 0 ){
			if( nuclear[i][j]==id ){
				flag = 1;
				break;
			}
			j++;
		}
		if( flag == 1 ){
			j = 2;
			index = 0;
			while( nuclear[i][j] != 0 ){
				if( nuclear[i][j] != id ){
					*(siblingID+index) = nuclear[i][j];
					index++;
				}
				j++;
			}
		}
	}
}


void descendant(int id, int *resultID){

	int i,j;
	int Nresult = 0;

	int *spouseID = (int *)calloc(1, sizeof(int));;
	search_spouse(id, spouseID);

	int *childrenID = (int *)calloc(MAX_SIBLING, sizeof(int));
	search_children(id, childrenID);

	if( *spouseID != 0 ){
		*(resultID+Nresult) = *spouseID;
		Nresult++;

		int *ancestersID = (int *)calloc(NIndividual, sizeof(int));
		search_ancesters(*spouseID, ancestersID);
		i = 0;
		while( *(ancestersID+i) != 0 ){
			*(resultID+Nresult) = *(ancestersID+i);
			Nresult++;
			i++;
		}
		free(ancestersID);
	}

	if( *childrenID != 0 ){
		i = 0;
		while( *(childrenID+i) != 0 ){
			*(resultID+Nresult) = *(childrenID+i);
			Nresult++;
			int *tmp = (int *)calloc(NIndividual, sizeof(int));
			search_descendants(*(childrenID+i), tmp);
			j = 0;
			while( *(tmp+j) != 0 ){
				*(resultID+Nresult) = *(tmp+j);
				Nresult++;
				j++;
			}
			free(tmp);
			i++;
		}
	}
	free(spouseID);
	free(childrenID);
}


void search_spouse(int id, int *spouseID){

	int i;

	for( i=1; i<N_NUCLEAR; i++ ){
		if( nuclear[i][0] == id ){
			*spouseID = nuclear[i][1];
			return;
		}
		else if( nuclear[i][1] == id ){
			*spouseID = nuclear[i][0];
			return;
		}
	}
}


void search_children(int id, int *childrenID){

	int i,j;

	for( i=1; i<N_NUCLEAR; i++ ){
		if( nuclear[i][0] == id || nuclear[i][1] == id ){
			j = 0;
			while( nuclear[i][j+2] != 0 ){
				*(childrenID+j) = nuclear[i][j+2];
				j++;
			}
			return;
		}
	}
}


void search_ancesters(int id, int *ancestersID){

	if( id == 0 ){
		return;
	}

	int i,j; 
	int NAncesters = 0;
	int *tmp;
	int *parentID = (int *)calloc(2, sizeof(int));
	search_parent(id, parentID);

	if( *parentID != 0 ){
		*ancestersID = *parentID;
		*(ancestersID+1) = *(parentID+1);
		NAncesters = 2;

		tmp = (int *)calloc(NIndividual, sizeof(int));
		search_ancesters(*parentID, tmp);
		j = 0;
		while( *(tmp+j) != 0 ){
			*(ancestersID+NAncesters) = *(tmp+j);
			NAncesters++;
			j++;
		}
		free(tmp);
		tmp = (int *)calloc(NIndividual, sizeof(int));
		search_ancesters(*(parentID+1), tmp);
		j = 0;
		while( *(tmp+j) != 0 ){
			*(ancestersID+NAncesters) = *(tmp+j);
			NAncesters++;
			j++;
		}
		free(tmp);
	}

	int *siblingID = (int *)calloc(MAX_SIBLING, sizeof(int));
	sibling(id, siblingID);

	if( *siblingID != 0 ){
		i = 0;
		while( *(siblingID+i) != 0 ){
			*(ancestersID+NAncesters) = *(siblingID+i);
			NAncesters++;
			tmp = (int *)calloc(NIndividual, sizeof(int));
			search_descendants(*(siblingID+i), tmp);
			j = 0;
			while( *(tmp+j) != 0 ){
				*(ancestersID+NAncesters) = *(tmp+j);
				NAncesters++;
				j++;
			}
			free(tmp);
			i++;
		}
	}

	free(parentID);
	free(siblingID);
}


void search_parent(int id, int *parentID){

	if( is_founder(id) == 1 ){
		return;
	}

	int i;
	for( i=0; i<NIndividual; i++ ){
		if( pedigree[i][0] == id ){
			*parentID = pedigree[i][1];
			*(parentID+1) = pedigree[i][2];
			return;
		}
	}
}


void search_descendants(int id, int *descendantsID){

	if( id == 0 ){
		return;
	}
	int NDescendants = 0;
	int i,j;
	
	int *spouseID = (int *)calloc(1, sizeof(int));
	search_spouse(id, spouseID);

	int *childrenID = (int *)calloc(MAX_SIBLING, sizeof(int));
	search_children(id, childrenID);

	if( *spouseID != 0 ){
		*(descendantsID+NDescendants) = *spouseID;
		NDescendants++;

		int *ancestersID = (int *)calloc(NIndividual, sizeof(int));
		search_ancesters(*spouseID, ancestersID);
		i = 0;
		while( *(ancestersID+i) != 0 ){
			*(descendantsID+NDescendants) = *(ancestersID+i);
			NDescendants++;
			i++;
		}
		free(ancestersID);
	}

	if( *childrenID != 0 ){
		i = 0;
		int *tmp;
		while( *(childrenID+i) != 0 ){
			*(descendantsID+NDescendants) = *(childrenID+i);
			NDescendants++;
			tmp = (int *)calloc(NIndividual, sizeof(int));
			search_descendants(*(childrenID+i), tmp);
			j = 0;
			while( *(tmp+j) != 0 ){
				*(descendantsID+NDescendants) = *(tmp+j);
				NDescendants++;
				j++;
			}
			free(tmp);
			i++;
		}
	}
	free(spouseID);
	free(childrenID);
}


void search_branch2(int id, int *branchID){

	int i;
	for( i=0; i<NIndividual; i++ ){
		if( *(ID_branch_candidate+i) == 0 ){
			*(ID_branch_candidate+i) = id;
			break;
		}
	}

	int *terminalID = (int *)calloc(NIndividual, sizeof(int));
	search_terminal(id, terminalID);
	int *tmp = (int *)calloc(NIndividual, sizeof(int));
	search_branch_in(*terminalID, tmp);

	i = 0;
	while( *(tmp+i)!=0 ){
		*(branchID+i) = *(tmp+i);
		i++;
	}
	free(terminalID);
	free(tmp);
}


void search_branch_in(int id, int *branchID){

	if( is_founder(id) == 1 ){
		return;
	}

	int i;
	for( i=0; i<NIndividual; i++ ){
		if( *(ID_branch_candidate+i) == 0 ){
			break;
		}
		if( *(ID_branch_candidate+i) == id ){
			return;
		}
	}

	int fid=0;
	int mid=0;
	for( i=0; i<NIndividual; i++ ){
		if( pedigree[i][0] == id ){
			fid = pedigree[i][1];
			mid = pedigree[i][2];
			break;
		}
	}

	int *siblingID = (int *)calloc(MAX_SIBLING, sizeof(int));
	sibling(id, siblingID);

	int *tmp_fid = (int *)calloc(NIndividual, sizeof(int));
	int *tmp_mid = (int *)calloc(NIndividual, sizeof(int));

	search_branch_in(fid, tmp_fid);
	search_branch_in(mid, tmp_mid);

	int NBranch = 0;
	i = 0;
	while( *(tmp_fid+i) != 0 ){
		*(branchID+NBranch) = *(tmp_fid+i);
		NBranch++;
		i++;
	}
	i = 0;
	while( *(tmp_mid+i) != 0 ){
		*(branchID+NBranch) = *(tmp_mid+i);
		NBranch++;
		i++;
	}
	free(tmp_fid);
	free(tmp_mid);

	if( *siblingID != 0 ){
		int j;
		i = 0;
		int *tmp;
		while( *(siblingID+i)!=0 ){
			if( is_terminal(*(siblingID+i)) == 0 ){
				tmp = (int *)calloc(NIndividual, sizeof(int));
				search_branch2(*(siblingID+i), tmp);
				*(branchID+NBranch) = *(siblingID+i);
				NBranch++;
				j = 0;
				while( *(tmp+j) != 0 ){
					*(branchID+NBranch) = *(tmp+j);
					NBranch++;
					j++;
				}
				free(tmp);
			}
			i++;
		}
	}
	free(siblingID);
}


void search_branch(int id, int *branchID){

	int i;
	int fid=0;
	int mid=0;
	int NBranch;

	for( i=0; i<NIndividual; i++ ){
		if( pedigree[i][0] == id ){
			fid = pedigree[i][1];
			mid = pedigree[i][2];
			break;
		}
	}

	int *tmp_fid = (int *)calloc(NIndividual, sizeof(int));
	int *tmp_mid = (int *)calloc(NIndividual, sizeof(int));

	search_branch_in(fid, tmp_fid);
	search_branch_in(mid, tmp_mid);

	NBranch = 0;
	i = 0;
	while( *(tmp_fid+i)!=0 ){
		*(branchID+NBranch) = *(tmp_fid+i);
		NBranch++;
		i++;
	}
	i = 0;
	while( *(tmp_mid+i)!=0 ){
		*(branchID+NBranch) = *(tmp_mid+i);
		NBranch++;
		i++;
	}
	free(tmp_fid);
	free(tmp_mid);
}




void branch(int *terminalID, int *resultID){

	int i,j,k,l;
	int N;
	int flag;
	int *order;
	int *branchID = (int *)calloc(NIndividual, sizeof(int));
	search_branch(*terminalID, branchID);

	int *branchesID = (int *)calloc(NIndividual, sizeof(int));
	int NBranches = 0;
	i = 0;
	while( *(branchID+i) != 0 ){
		j = 0;
		flag = 0;
		while( *(terminalID+j) != 0 ){
			if( *(branchID+i) == *(terminalID+j) ){
				flag = 1;
				break;
			}
			j++;
		}
		if( flag == 0 ){
			*(branchesID+NBranches) = *(branchID+i);
			NBranches++;
		}
		i++;
	}

	int *SizeBranchesID = (int *)calloc(NIndividual, sizeof(int));
	i = 0;
	int *tmp;
	while( *(branchesID+i) != 0 ){
		tmp = (int *)calloc(NIndividual, sizeof(int));
		descendant(*(branchesID+i), tmp);
		N = 0;
		while( *(tmp+N) != 0 ){
			N++;
		}
		*(SizeBranchesID+i) = N;
		free(tmp);
		i++;
	}

	order = (int *)malloc(i*sizeof(int));
	for( k=0; k<i; k++ ){
		*(order+k) = 0;
		for( l=0; l<i; l++ ){
			if(*(SizeBranchesID+k) > *(SizeBranchesID+l)){
				*(order+k) = *(order+k)+1;
			}
			else if(*(SizeBranchesID+k) == *(SizeBranchesID+l) && k>l){
				*(order+k) = *(order+k)+1;
			}
		}
	}

	for( k=0; k<i; k++ ){
		for( l=0; l<i; l++ ){
			if( *(order+l) == k ){
				*(resultID+k) = *(branchesID+l);
				break;
			}
		}
	}

	free(SizeBranchesID);
	free(order);
	free(branchID);
	free(branchesID);
}



void remove_0terms(polynomial *poly){

	int i,n;
	double *result_coef;
	int *result_power;
	result_coef = (double *)malloc((*poly).size*sizeof(double));
	result_power = (int *)malloc((*poly).size*sizeof(int));

	n=0;
	for( i=0; i<(*poly).size; i++ ){
		if ( (*poly).coef[i] != 0 ){
			*(result_coef+n) = (*poly).coef[i];
			*(result_power+n) = (*poly).power[i];
			n++;
		}
	}

	if( n == 0 ){
		(*poly).coef = (double *)realloc((*poly).coef, sizeof(double));
		(*poly).power = (int *)realloc((*poly).power, sizeof(int));
		*((*poly).coef)=0;
		*((*poly).power)=0;
		(*poly).size = 1;
	}
	else {
		(*poly).coef = (double *)realloc((*poly).coef, n*sizeof(double));
		(*poly).power = (int *)realloc((*poly).power, n*sizeof(int));
		(*poly).size = n;
		for(i=0; i<n; i++){
			(*poly).coef[i] = *(result_coef+i);
			(*poly).power[i] = *(result_power+i);
		}
	}
	free(result_coef);
	free(result_power);
}


void eval_min(int power1, int power2, int power3, int power4, int *min, int *order){

	int tmin=-1;
	int torder=0;

	if( power1!=-1 ){
		tmin = power1;
		torder = 1;
	}
	if( power2!=-1 && (tmin==-1||power2<tmin) ){
		tmin = power2;
		torder = 2;
	}
	if( power3!=-1 && (tmin==-1||power3 < tmin) ){
		tmin = power3;
		torder = 3;
	}
	if( power4!=-1 && (tmin==-1||power4 < tmin) ){
		tmin = power4;
		torder = 4;
	}

	*min = tmin;
	*order = torder;
}


void sum(polynomial_individual *poly, polynomial *result){

	int i;
	int n=0;
	int n1=0;
	int n2=0;
	int n3=0;
	int n4=0;
	int order;
	int mini;
	int power1;
	int power2;
	int power3;
	int power4;

	int size = (*poly).poly1.size + (*poly).poly2.size + (*poly).poly3.size + (*poly).poly4.size;

	double *tmp_coef;
	int *tmp_power;
	tmp_coef = (double *)malloc(size * sizeof(double));
	tmp_power = (int *)malloc(size * sizeof(int));

	if( (*poly).poly1.size == 0 ){
		power1 = -1;
	}
	else{
		power1 = (*poly).poly1.power[n1];
	}
	if( (*poly).poly2.size == 0 ){
		power2 = -1;
	}
	else{
		power2 = (*poly).poly2.power[n2];
	}
	if( (*poly).poly3.size == 0 ){
		power3 = -1;
	}
	else{
		power3 = (*poly).poly3.power[n3];
	}
	if( (*poly).poly4.size == 0 ){
		power4 = -1;
	}
	else{
		power4 = (*poly).poly4.power[n4];
	}

	eval_min(power1,power2,power3,power4,&mini,&order);

	while(1){
		if( order == 1 ){
			tmp_coef[n] = (*poly).poly1.coef[n1];
			tmp_power[n] = (*poly).poly1.power[n1];
			n1++;
			if( n2<(*poly).poly2.size && mini==power2){
				tmp_coef[n] += (*poly).poly2.coef[n2];
				n2++;
			}
			if( n3<(*poly).poly3.size && mini==power3){
				tmp_coef[n] += (*poly).poly3.coef[n3];
				n3++;
			}
			if( n4<(*poly).poly4.size && mini==power4){
				tmp_coef[n] += (*poly).poly4.coef[n4];
				n4++;
			}
		}
		else if( order == 2 ){
			tmp_coef[n] = (*poly).poly2.coef[n2];
			tmp_power[n] = (*poly).poly2.power[n2];
			n2++;
			if( n3<(*poly).poly3.size && mini==power3){
				tmp_coef[n] += (*poly).poly3.coef[n3];
				n3++;
			}
			if( n4<(*poly).poly4.size && mini==power4){
				tmp_coef[n] += (*poly).poly4.coef[n4];
				n4++;
			}
		}
		else if( order == 3 ){
			tmp_coef[n] = (*poly).poly3.coef[n3];
			tmp_power[n] = (*poly).poly3.power[n3];
			n3++;
			if( n4<(*poly).poly4.size && mini==power4){
				tmp_coef[n] += (*poly).poly4.coef[n4];
				n4++;
			}
		}
		else if( order == 4 ){
			tmp_coef[n] = (*poly).poly4.coef[n4];
			tmp_power[n] = (*poly).poly4.power[n4];
			n4++;
		}

		if( tmp_coef[n] != 0 ){
			n++;
		}

		if( power1==-1 || n1 == (*poly).poly1.size ){
			power1 = -1;
		}
		else{
			power1 = (*poly).poly1.power[n1];
		}
		if( power2==-1 || n2 == (*poly).poly2.size ){
			power2 = -1;
		}
		else{
			power2 = (*poly).poly2.power[n2];
		}
		if( power3==-1 || n3 == (*poly).poly3.size ){
			power3 = -1;
		}
		else{
			power3 = (*poly).poly3.power[n3];
		}
		if( power4==-1 || n4 == (*poly).poly4.size ){
			power4 = -1;
		}
		else{
			power4 = (*poly).poly4.power[n4];
		}
		if(((power1==-1 && power2==-1) && power3==-1) && power4==-1){
			break;
		}
		eval_min(power1,power2,power3,power4,&mini,&order);
	}

	if( n == 0 ){
		(*result).coef = (double *)realloc((*result).coef, sizeof(double));
		(*result).power = (int *)realloc((*result).power, sizeof(int));
		*((*result).coef)=0;
		*((*result).power)=0;
		(*result).size = 1;
		return;
	}
	(*result).coef = (double *)realloc((*result).coef, n*sizeof(double));
	(*result).power = (int *)realloc((*result).power, n*sizeof(int));
	(*result).size = n;

	for(i=0; i<n; i++){
		(*result).coef[i] = *(tmp_coef+i);
		(*result).power[i] = *(tmp_power+i);
	}
	free(tmp_coef);
	free(tmp_power);

}


void prod_scalar(double a, polynomial *poly, polynomial *result){

	int i;

	if( a == 0 ){
		(*result).coef = (double *)realloc((*result).coef,sizeof(double));
		(*result).power = (int *)realloc((*result).power,sizeof(int));
		*((*result).coef) = 0;
		*((*result).power) = 0;
		(*result).size = 1;
		return;
	}

	substitute_polynomial(result, poly);
	for(i=0; i<(*poly).size; i++){
		(*result).coef[i] = (*poly).coef[i] * a;
	}
}


void prod_polynomial(polynomial *poly1, polynomial *poly2, polynomial *result){

	int i,j;
	int n=0;
	int size = (*poly1).size * (*poly2).size;

	if( (*poly1).size==1 || (*poly2).size==1 ){
		(*result).coef = (double *)realloc((*result).coef, size*sizeof(double));
		(*result).power = (int *)realloc((*result).power, size*sizeof(int));
		(*result).size = size;
		for(i=0; i<(*poly1).size; i++){
			for(j=0; j<(*poly2).size; j++){
				(*result).coef[n] = (*poly1).coef[i] * (*poly2).coef[j];
				(*result).power[n] = (*poly1).power[i] + (*poly2).power[j];
				n++;
			}
		}
		return;	
	}

	double *tmp_coef;
	int *tmp_power;
	tmp_coef = (double *)calloc(size, sizeof(double));
	tmp_power = (int *)malloc(size * sizeof(int));

	int tmp,tmpi;
	int current_position=0;
	int *xposition;
	int *xvalue;

	xposition = (int *)calloc((*poly2).size, sizeof(int));
	*xposition = 1;
	xvalue = (int *)malloc((*poly2).size*sizeof(int));
	*xvalue = (*poly1).power[1] + (*poly2).power[0];
	for(i=1; i<(*poly2).size; i++){
		*(xvalue+i) = (*poly1).power[0] + (*poly2).power[i];
	}

	*tmp_coef = (*poly1).coef[0] * (*poly2).coef[0];
	*tmp_power = (*poly1).power[0] + (*poly2).power[0];
	n=1;

	while(1){
		if(current_position == (*poly2).size-1){
			for(i=*(xposition+current_position); i<(*poly1).size; i++){
				*(tmp_coef+n) = (*poly1).coef[i] * (*poly2).coef[(*poly2).size-1];
				*(tmp_power+n) = (*poly1).power[i] + (*poly2).power[(*poly2).size-1];
				n++;
			}
			break;
		}

		tmp = *(xvalue+current_position);
		tmpi = current_position;
		for(i=current_position+1; i<(*poly2).size; i++){
			if( tmp > *(xvalue+i) ){
				tmp = *(xvalue+i);
				tmpi = i;
			}
		}
		
		*(tmp_power+n) = tmp;
		for(i=tmpi; i<(*poly2).size; i++){
			if( *(xvalue+i) == tmp ){
				*(tmp_coef+n) += (*poly1).coef[*(xposition+i)] * (*poly2).coef[i];
				*(xposition+i) = *(xposition+i)+1;
				if(*(xposition+i)==(*poly1).size){
					current_position++;
				}
				else{
					*(xvalue+i) = (*poly1).power[*(xposition+i)] + (*poly2).power[i];
				}
			}
		}
		if(current_position==(*poly2).size && *(xposition+current_position)==(*poly1).size){
			break;
		}
		if( *(tmp_coef+n) != 0 ){
			n++;
		}
	}

	(*result).coef = (double *)realloc((*result).coef, n*sizeof(double));
	(*result).power = (int *)realloc((*result).power, n*sizeof(int));
	(*result).size = n;

	for(i=0; i<n; i++){
		(*result).coef[i] = *(tmp_coef+i);
		(*result).power[i] = *(tmp_power+i);
	}

	free(tmp_coef);
	free(tmp_power);
	free(xvalue);
	free(xposition);
}


void set_likelihood_conditional(polynomial_individual *poly, int row){

	if( row == 1 ){
		*((*poly).poly1.coef) = 1;
		*((*poly).poly1.power) = 0;
		*((*poly).poly2.coef) = 0;
		*((*poly).poly2.power) = 0;
		*((*poly).poly3.coef) = 0;
		*((*poly).poly3.power) = 0;
		*((*poly).poly4.coef) = 0;
		*((*poly).poly4.power) = 0;
	}
	else if( row == 2 ){
		*((*poly).poly1.coef) = 0;
		*((*poly).poly1.power) = 0;
		*((*poly).poly2.coef) = 1;
		*((*poly).poly2.power) = 0;
		*((*poly).poly3.coef) = 0;
		*((*poly).poly3.power) = 0;
		*((*poly).poly4.coef) = 0;
		*((*poly).poly4.power) = 0;
	}
	else if( row == 3 ){
		*((*poly).poly1.coef) = 0;
		*((*poly).poly1.power) = 0;
		*((*poly).poly2.coef) = 0;
		*((*poly).poly2.power) = 0;
		*((*poly).poly3.coef) = 1;
		*((*poly).poly3.power) = 0;
		*((*poly).poly4.coef) = 0;
		*((*poly).poly4.power) = 0;
	}
	else if( row == 4 ){
		*((*poly).poly1.coef) = 0;
		*((*poly).poly1.power) = 0;
		*((*poly).poly2.coef) = 0;
		*((*poly).poly2.power) = 0;
		*((*poly).poly3.coef) = 0;
		*((*poly).poly3.power) = 0;
		*((*poly).poly4.coef) = 1;
		*((*poly).poly4.power) = 0;
	}
	(*poly).poly1.size = 1;
	(*poly).poly2.size = 1;
	(*poly).poly3.size = 1;
	(*poly).poly4.size = 1;
}


void initialise_likelihood0(polynomial *poly){

	(*poly).coef = (double *)realloc((*poly).coef, sizeof(double));
	(*poly).power = (int *)realloc((*poly).power, sizeof(int));
	*((*poly).coef) = 0;
	*((*poly).power) = 0;
	(*poly).size = 1;

}


double InheritP(int genotype_parent, int origin_id, int allele_id){

	int allele1=0;
	int allele2=0;
	int allele3=0;

	if( genotype_parent == 1){
		allele1 = allele2 = 1;
	}
	else if( genotype_parent == 2 || genotype_parent == 3 ){
		allele1 = 1;
		allele2 = 2;
	}
	else if( genotype_parent == 4){
		allele1 = allele2 = 2;
	}
	if( origin_id == 1 ){
		if( allele_id == 1 || allele_id ==  2){
			allele3 = 1;
		}
		else if( allele_id == 3 || allele_id == 4 ){
			allele3 = 2;
		}
	}
	else if( origin_id == 2 ){
		if( allele_id == 1 || allele_id ==  3){
			allele3 = 1;
		}
		else if( allele_id == 2 || allele_id == 4 ){
			allele3 = 2;
		}
	}

	if( allele1 == allele2 ){
		if( allele1 == allele3 ){
			return 1;
		}
		else{
			return 0;
		}
	}
	else{
		return 0.5;
	}
}


void penetrance(int Phenotype, int row, polynomial *poly){

	if( Phenotype == 0 ){
		return;
	}

	int i,j,n;
	int flag = 0;
	for( i=0; i<(*poly).size; i++){
		if( (*poly).coef != 0 ){
			flag = 1;
			break;
		}
	}

	if( flag == 0 ){
		initialise_likelihood0(&(*poly));
		return;
	}

	if( Phenotype == 2 ){
		if( row == 1 ){
			for( i=0; i<(*poly).size; i++ ){
				(*poly).power[i]++;
			}
		}
		else if( row == 2 || row == 3 ){
			for( i=0; i<(*poly).size; i++ ){
				(*poly).power[i] += MAX_POWER;
			}
		}
		else if( row == 4 ){
			for( i=0; i<(*poly).size; i++ ){
				(*poly).power[i] += MAX_POWER*MAX_POWER;
			}
		}
	}
	else if( Phenotype == 1 ){

		double *coef;
		int *power;
		coef = (double *)malloc((*poly).size*2*sizeof(double));
		power = (int *)malloc((*poly).size*2*sizeof(int));

		for( i=0; i<(*poly).size; i++){
			*(coef+i) = (*poly).coef[i];
			*(power+i) = (*poly).power[i];
		}
		n = (*poly).size;
		for( i=0; i<(*poly).size; i++){
			*(coef+n) = -(*poly).coef[i];
			n++;
		}

		n = (*poly).size;
		if( row == 1 ){
			for( i=0; i<(*poly).size; i++ ){
				*(power+n) = (*poly).power[i]+1;
				n++;
			}
		}
		else if( row == 2 || row == 3 ){
			for( i=0; i<(*poly).size; i++ ){
				*(power+n) = (*poly).power[i]+MAX_POWER;
				n++;
			}
		}
		else if( row == 4 ){
			for( i=0; i<(*poly).size; i++ ){
				*(power+n) = (*poly).power[i]+MAX_POWER*MAX_POWER;
				n++;
			}
		}

		double *tcoef;
		int *tpower;
		tcoef = (double *)malloc((*poly).size*2*sizeof(double));
		tpower = (int *)malloc((*poly).size*2*sizeof(int));

		int n1 = 0;
		int n2 = (*poly).size;
		n = 0;
		while(1){
			if( n1 == (*poly).size ){
				for( j=n2; j<(*poly).size*2; j++ ){
					*(tcoef+n) = *(coef+j);
					*(tpower+n) = *(power+j);
					n++;
				}
				break;
			}
			if( n2 == 2*(*poly).size ){
				for( j=n1; j<(*poly).size; j++ ){
					*(tcoef+n) = *(coef+j);
					*(tpower+n) = *(power+j);
					n++;
				}
				break;
			}
			if( *(power+n1) < *(power+n2) ){
				*(tpower+n) = *(power+n1);
				*(tcoef+n) = *(coef+n1);
				n++;
				n1++;
			}
			else if( *(power+n1) > *(power+n2) ){
				*(tpower+n) = *(power+n2);
				*(tcoef+n) = *(coef+n2);
				n++;
				n2++;
			}
			else if( *(power+n1) == *(power+n2) ){
				*(tcoef+n) = *(coef+n1) + *(coef+n2);
				if( *(tcoef+n) != 0 ){
					*(tpower+n) = *(power+n1);
					n++;
				}
				n1++;
				n2++;
			}
		}

		if( n == 0 ){
			(*poly).coef = (double *)realloc((*poly).coef, sizeof(double));
			(*poly).power = (int *)realloc((*poly).power, sizeof(int));
			*((*poly).coef)=0;
			*((*poly).power)=0;
			(*poly).size = 1;
		}
		else{
			(*poly).coef = (double *)realloc((*poly).coef, n*sizeof(double));
			(*poly).power = (int *)realloc((*poly).power, n*sizeof(int));
			(*poly).size = n;
			for(i=0; i<n; i++){
				(*poly).coef[i] = *(tcoef+i);
				(*poly).power[i] = *(tpower+i);
			}
		}
		free(coef);
		free(power);
		free(tcoef);
		free(tpower);
	}
}


void Xi_cal_founderP_3paras(polynomial_individual *poly, int Phenotype){

	if( Phenotype == 0 ){
		*((*poly).poly1.coef) = fA*fA;
		*((*poly).poly1.power) = 0;
		(*poly).poly1.size = 1;
		*((*poly).poly2.coef) = fA*(1-fA);
		*((*poly).poly2.power) = 0;
		(*poly).poly2.size = 1;
		*((*poly).poly3.coef) = fA*(1-fA);
		*((*poly).poly3.power) = 0;
		(*poly).poly3.size = 1;
		*((*poly).poly4.coef) = (1-fA)*(1-fA);
		*((*poly).poly4.power) = 0;
		(*poly).poly4.size = 1;
	}
	else if( Phenotype == 1 ){
		(*poly).poly1.coef = (double *)realloc((*poly).poly1.coef, 2 * sizeof(double));
		(*poly).poly1.coef[0] = fA*fA;
		(*poly).poly1.coef[1] = -(fA*fA);
		(*poly).poly1.power = (int *)realloc((*poly).poly1.power, 2 * sizeof(int));
		(*poly).poly1.power[0] = 0;
		(*poly).poly1.power[1] = 1;
		(*poly).poly1.size = 2;
		(*poly).poly2.coef = (double *)realloc((*poly).poly2.coef, 2 * sizeof(double));
		(*poly).poly2.coef[0] = fA*(1-fA);
		(*poly).poly2.coef[1] = -(fA*(1-fA));
		(*poly).poly2.power = (int *)realloc((*poly).poly2.power, 2 * sizeof(int));
		(*poly).poly2.power[0] = 0;
		(*poly).poly2.power[1] = MAX_POWER;
		(*poly).poly2.size = 2;
		(*poly).poly3.coef = (double *)realloc((*poly).poly3.coef, 2 * sizeof(double));
		(*poly).poly3.coef[0] = fA*(1-fA);
		(*poly).poly3.coef[1] = -(fA*(1-fA));
		(*poly).poly3.power = (int *)realloc((*poly).poly3.power, 2 * sizeof(int));
		(*poly).poly3.power[0] = 0;
		(*poly).poly3.power[1] = MAX_POWER;
		(*poly).poly3.size = 2;
		(*poly).poly4.coef = (double *)realloc((*poly).poly4.coef, 2 * sizeof(double));
		(*poly).poly4.coef[0] = (1-fA)*(1-fA);
		(*poly).poly4.coef[1] = -((1-fA)*(1-fA));
		(*poly).poly4.power = (int *)realloc((*poly).poly4.power, 2 * sizeof(int));
		(*poly).poly4.power[0] = 0;
		(*poly).poly4.power[1] = MAX_POWER*MAX_POWER;
		(*poly).poly4.size = 2;
	}
	else if( Phenotype == 2 ){
		*((*poly).poly1.coef) = fA*fA;
		*((*poly).poly1.power) = 1;
		(*poly).poly1.size = 1;
		*((*poly).poly2.coef) = fA*(1-fA);
		*((*poly).poly2.power) = MAX_POWER;
		(*poly).poly2.size = 1;
		*((*poly).poly3.coef) = fA*(1-fA);
		*((*poly).poly3.power) = MAX_POWER;
		(*poly).poly3.size = 1;
		*((*poly).poly4.coef) = (1-fA)*(1-fA);
		*((*poly).poly4.power) = MAX_POWER*MAX_POWER;
		(*poly).poly4.size = 1;
	}
}


void AncestralLikelihoodDiplo_C(int id_ith, int *siblingID){

	int Phenotype = pedigree[id_ith][4];

	polynomial *polyf,*polym,*polym_1,*polym_23,*polym_4;
	polynomial_individual *tmp,*tmp2_1,*tmp2_23,*tmp2_4,*tmp3;
	int fid_ith,mid_ith;
	int i;

	fid_ith = get_id_ith(pedigree[id_ith][1]);
	mid_ith = get_id_ith(pedigree[id_ith][2]);

	AncestralLikelihood_C(fid_ith);
	AncestralLikelihood_C(mid_ith);

	if( *siblingID == 0 ){
		for( i=1; i<5; i++ ){
			tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp);
			prod_scalar(InheritP(1,1,i), &likelihood_conditional[fid_ith].poly1, &((*tmp).poly1));
			prod_scalar(InheritP(2,1,i), &likelihood_conditional[fid_ith].poly2, &((*tmp).poly2));
			prod_scalar(InheritP(3,1,i), &likelihood_conditional[fid_ith].poly3, &((*tmp).poly3));
			prod_scalar(InheritP(4,1,i), &likelihood_conditional[fid_ith].poly4, &((*tmp).poly4));
			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			sum(tmp, polyf);
			free_polynomial_individual(tmp);

			tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp);
			prod_scalar(InheritP(1,2,i), &likelihood_conditional[mid_ith].poly1, &((*tmp).poly1));
			prod_scalar(InheritP(2,2,i), &likelihood_conditional[mid_ith].poly2, &((*tmp).poly2));
			prod_scalar(InheritP(3,2,i), &likelihood_conditional[mid_ith].poly3, &((*tmp).poly3));
			prod_scalar(InheritP(4,2,i), &likelihood_conditional[mid_ith].poly4, &((*tmp).poly4));
			polym = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polym);
			sum(tmp, polym);
			free_polynomial_individual(tmp);

			if( i==1 ){
				prod_polynomial(polyf, polym, &likelihood_conditional[id_ith].poly1);
				penetrance(Phenotype, 1, &likelihood_conditional[id_ith].poly1);
			}
			else if( i==2 ){
				prod_polynomial(polyf, polym, &likelihood_conditional[id_ith].poly2);
				penetrance(Phenotype, 2, &likelihood_conditional[id_ith].poly2);
			}
			else if( i==3 ){
				prod_polynomial(polyf, polym, &likelihood_conditional[id_ith].poly3);
				penetrance(Phenotype, 3, &likelihood_conditional[id_ith].poly3);
			}
			else if( i==4 ){
				prod_polynomial(polyf, polym, &likelihood_conditional[id_ith].poly4);
				penetrance(Phenotype, 4, &likelihood_conditional[id_ith].poly4);
			}
			free_polynomial(polyf);
			free_polynomial(polym);
		}
	}
	else{
		for( i=1; i<5; i++ ){
			tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp);
			prod_scalar(InheritP(1,2,i), &likelihood_conditional[mid_ith].poly1, &((*tmp).poly1));
			prod_scalar(InheritP(2,2,i), &likelihood_conditional[mid_ith].poly2, &((*tmp).poly2));
			prod_scalar(InheritP(3,2,i), &likelihood_conditional[mid_ith].poly3, &((*tmp).poly3));
			prod_scalar(InheritP(4,2,i), &likelihood_conditional[mid_ith].poly4, &((*tmp).poly4));

			tmp2_1 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			tmp2_23 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			tmp2_4 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp2_1);
			allocate_memory_polynomial_individual(tmp2_23);
			allocate_memory_polynomial_individual(tmp2_4);

			prod_polynomial(&((*tmp).poly1), &likelihood_sibling[id_ith].poly1, &((*tmp2_1).poly1));
			prod_polynomial(&((*tmp).poly2), &likelihood_sibling[id_ith].poly2, &((*tmp2_1).poly2));
			prod_polynomial(&((*tmp).poly3), &likelihood_sibling[id_ith].poly2, &((*tmp2_1).poly3));
			prod_polynomial(&((*tmp).poly4), &likelihood_sibling[id_ith].poly3, &((*tmp2_1).poly4));

			prod_polynomial(&((*tmp).poly1), &likelihood_sibling[id_ith].poly4, &((*tmp2_23).poly1));
			prod_polynomial(&((*tmp).poly2), &likelihood_sibling[id_ith].poly5, &((*tmp2_23).poly2));
			prod_polynomial(&((*tmp).poly3), &likelihood_sibling[id_ith].poly5, &((*tmp2_23).poly3));
			prod_polynomial(&((*tmp).poly4), &likelihood_sibling[id_ith].poly6, &((*tmp2_23).poly4));

			prod_polynomial(&((*tmp).poly1), &likelihood_sibling[id_ith].poly7, &((*tmp2_4).poly1));
			prod_polynomial(&((*tmp).poly2), &likelihood_sibling[id_ith].poly8, &((*tmp2_4).poly2));
			prod_polynomial(&((*tmp).poly3), &likelihood_sibling[id_ith].poly8, &((*tmp2_4).poly3));
			prod_polynomial(&((*tmp).poly4), &likelihood_sibling[id_ith].poly9, &((*tmp2_4).poly4));

			polym_1 = (polynomial *)malloc(sizeof(polynomial));
			polym_23 = (polynomial *)malloc(sizeof(polynomial));
			polym_4 = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polym_1);
			allocate_memory_polynomial(polym_23);
			allocate_memory_polynomial(polym_4);
			sum(tmp2_1,polym_1);
			sum(tmp2_23,polym_23);
			sum(tmp2_4,polym_4);

			free_polynomial_individual(tmp2_1);
			free_polynomial_individual(tmp2_23);
			free_polynomial_individual(tmp2_4);

			tmp3 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp3);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(1,1,i),&likelihood_conditional[fid_ith].poly1,polyf);

			prod_polynomial(polyf,polym_1,&((*tmp3).poly1));
			free_polynomial(polyf);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(2,1,i),&likelihood_conditional[fid_ith].poly2,polyf);
			prod_polynomial(polyf,polym_23,&((*tmp3).poly2));
			free_polynomial(polyf);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(3,1,i),&likelihood_conditional[fid_ith].poly3,polyf);
			prod_polynomial(polyf,polym_23,&((*tmp3).poly3));
			free_polynomial(polyf);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(4,1,i),&likelihood_conditional[fid_ith].poly4,polyf);
			prod_polynomial(polyf,polym_4,&((*tmp3).poly4));
			free_polynomial(polyf);

			if( i==1 ){
				sum(tmp3,&likelihood_conditional[id_ith].poly1);
				penetrance(Phenotype,1,&likelihood_conditional[id_ith].poly1);
			}
			else if( i==2 ){
				sum(tmp3,&likelihood_conditional[id_ith].poly2);
				penetrance(Phenotype,2,&likelihood_conditional[id_ith].poly2);
			}
			else if( i==3 ){
				sum(tmp3,&likelihood_conditional[id_ith].poly3);
				penetrance(Phenotype,3,&likelihood_conditional[id_ith].poly3);
			}
			else if( i==4 ){
				sum(tmp3,&likelihood_conditional[id_ith].poly4);
				penetrance(Phenotype,4,&likelihood_conditional[id_ith].poly4);
			}

			free_polynomial_individual(tmp);
			free_polynomial_individual(tmp3);
			free_polynomial(polym_1);
			free_polynomial(polym_23);
			free_polynomial(polym_4);
		}
	}
}


void SumOfSiblingLikelihood(int row_f, int row_m, int id, polynomial *poly){

	int j;
	int id_ith = get_id_ith(id);
	int Phenotype = pedigree[id_ith][4];

	double *prob = (double *)malloc(4*sizeof(double));
	for( j=1; j<5; j++ ){
		*(prob+j-1) = InheritP(row_f, 1, j) * InheritP(row_m, 2, j);
	}
	if( is_terminal(id) == 1 ){
		if( Phenotype == 0 ){
			*((*poly).coef) = *prob + *(prob+1) + *(prob+2) + *(prob+3);
			*((*poly).power) = 0;
			(*poly).size = 1;
		}
		else if( Phenotype == 2 ){
			(*poly).coef = (double *)realloc((*poly).coef, 3*sizeof(double));
			(*poly).power = (int *)realloc((*poly).power, 3*sizeof(int));
			(*poly).size = 3;
			*((*poly).coef) = *prob;
			*((*poly).power) = 1;
			*((*poly).coef+1) = *(prob+1) + *(prob+2);
			*((*poly).power+1) = MAX_POWER;
			*((*poly).coef+2) = *(prob+3);
			*((*poly).power+2) = MAX_POWER * MAX_POWER;
			remove_0terms(poly);
		}
		else if( Phenotype == 1 ){
			(*poly).coef = (double *)realloc((*poly).coef, 4*sizeof(double));
			(*poly).power = (int *)realloc((*poly).power, 4*sizeof(int));
			(*poly).size = 4;
			*((*poly).coef) = *prob + *(prob+1) + *(prob+2) + *(prob+3);
			*((*poly).power) = 0;
			*((*poly).coef+1) = -(*prob);
			*((*poly).power+1) = 1;
			*((*poly).coef+2) = -(*(prob+1) + *(prob+2));
			*((*poly).power+2) = MAX_POWER;
			*((*poly).coef+3) = -*(prob+3);
			*((*poly).power+3) = MAX_POWER * MAX_POWER;
			remove_0terms(poly);
		}
	}
	else{
		polynomial_individual *tmp;
		tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
		allocate_memory_polynomial_individual(tmp);
		substitute_polynomial(&(*tmp).poly1, &likelihood[id_ith].poly1);
		substitute_polynomial(&(*tmp).poly2, &likelihood[id_ith].poly2);
		substitute_polynomial(&(*tmp).poly3, &likelihood[id_ith].poly2);
		substitute_polynomial(&(*tmp).poly4, &likelihood[id_ith].poly4);

		for( j=0; j<(*tmp).poly1.size; j++ ){
			*((*tmp).poly1.coef+j) *= *prob;
		}
		for( j=0; j<(*tmp).poly2.size; j++ ){
			*((*tmp).poly2.coef+j) *= *(prob+1);
		}
		for( j=0; j<(*tmp).poly3.size; j++ ){
			*((*tmp).poly3.coef+j) *= *(prob+2);
		}
		for( j=0; j<(*tmp).poly4.size; j++ ){
			*((*tmp).poly4.coef+j) *= *(prob+3);
		}
		penetrance(Phenotype, 1, &((*tmp).poly1));
		penetrance(Phenotype, 2, &((*tmp).poly2));
		penetrance(Phenotype, 3, &((*tmp).poly3));
		penetrance(Phenotype, 4, &((*tmp).poly4));
		sum(tmp, poly);
		free_polynomial_individual(tmp);
	}
}


void substitute_polynomial(polynomial *to, polynomial *from){

	int i;
	(*to).coef = (double *)realloc((*to).coef, (*from).size*sizeof(double));
	(*to).power = (int *)realloc((*to).power, (*from).size*sizeof(int));

	for( i=0; i<(*from).size; i++){
		*((*to).coef+i) = *((*from).coef+i);
		*((*to).power+i) = *((*from).power+i);
	}
	(*to).size = (*from).size;
}

void calculate_likelihood_sibling(polynomial_sibling *poly, int *siblingID){

	if( (*poly).poly1.size == 0 ){
		int i,j,l;
		polynomial *tmp1,*tmp2,*result;
		for( i=1; i<5; i++ ){
			if( i != 3 ){
				for( j=1; j<5; j++ ){
					if( j != 3 ){
						result = (polynomial *)malloc(sizeof(polynomial));
						allocate_memory_polynomial(result);
						SumOfSiblingLikelihood(i, j, *siblingID, result);
						l = 1;
						while(*(siblingID+l)!=0){
							tmp1 = (polynomial *)malloc(sizeof(polynomial));
							allocate_memory_polynomial(tmp1);
							SumOfSiblingLikelihood(i, j, *(siblingID+l), tmp1);
							tmp2 = (polynomial *)malloc(sizeof(polynomial));
							allocate_memory_polynomial(tmp2);
							prod_polynomial(result, tmp1, tmp2);
							substitute_polynomial(result,tmp2);
							free_polynomial(tmp1);
							free_polynomial(tmp2);
							l++;
						}
						if( i == 1 && j == 1 ){
							substitute_polynomial(&((*poly).poly1), result);
						}
						else if( i == 1 && j == 2 ){
							substitute_polynomial(&((*poly).poly2), result);
						}
						else if( i == 1 && j == 4 ){
							substitute_polynomial(&((*poly).poly3), result);
						}
						else if( i == 2 && j == 1 ){
							substitute_polynomial(&((*poly).poly4), result);
						}
						else if( i == 2 && j == 2 ){
							substitute_polynomial(&((*poly).poly5), result);
						}
						else if( i == 2 && j == 4 ){
							substitute_polynomial(&((*poly).poly6), result);
						}
						else if( i == 4 && j == 1 ){
							substitute_polynomial(&((*poly).poly7), result);
						}
						else if( i == 4 && j == 2 ){
							substitute_polynomial(&((*poly).poly8), result);
						}
						else if( i == 4 && j == 4 ){
							substitute_polynomial(&((*poly).poly9), result);
						}
						free_polynomial(result);
					}
				}
			}
		}
	}
}


void AncestralLikelihood_C(int id_ith){

	int id = pedigree[id_ith][0];

	if( likelihood_conditional[id_ith].poly1.size == 0 ){
		if( id == ID ){
			set_likelihood_conditional(&likelihood_conditional[id_ith], ROW);
		}
		else{
			if( is_founder(id) == 1 ){
				int Phenotype = pedigree[id_ith][4];
				Xi_cal_founderP_3paras(&likelihood_conditional[id_ith], Phenotype);
			}
			else{
				int *siblingID = (int *)calloc(NIndividual, sizeof(int));
				sibling(id, siblingID);
				if( *siblingID != 0 ){
					calculate_likelihood_sibling(&likelihood_sibling[id_ith], siblingID);
				}
				AncestralLikelihoodDiplo_C(id_ith, siblingID);
			}
		}
	}
}


void ConditionalLikelihoodDiplo(int id_ith, int row){
	
	int terminal_ith;

	likelihood_conditional = (polynomial_individual *)malloc(sizeof(polynomial_individual)*NIndividual);
	allocate_memory_polynomial_individual_all(likelihood_conditional);

	ID = pedigree[id_ith][0];
	ROW = row;

	int *terminals = (int *)calloc(NIndividual, sizeof(int));
	search_terminal(pedigree[id_ith][0], terminals);

	terminal_ith = get_id_ith(*terminals);
	AncestralLikelihood_C(terminal_ith);

	if( row == 1 ){
		sum(&likelihood_conditional[terminal_ith], &likelihood[id_ith].poly1);
	}
	else if( row == 2 ){
		sum(&likelihood_conditional[terminal_ith], &likelihood[id_ith].poly2);
	}
	else if( row == 4 ){
		sum(&likelihood_conditional[terminal_ith], &likelihood[id_ith].poly4);
	}

	free(terminals);
	free_likelihood_conditional();
}



void clkh(int *ids){

	int id_ith;
	int i = 0;
	int *descendantID;

	while( *(ids+i) != 0 ){
		id_ith = get_id_ith(*(ids+i));
		ConditionalLikelihoodDiplo(id_ith, 1);
		ConditionalLikelihoodDiplo(id_ith, 2);
		ConditionalLikelihoodDiplo(id_ith, 4);

		descendantID = (int *)calloc(NIndividual, sizeof(int));
		descendant(*(ids+i), descendantID);
		free_likelihood(descendantID);
		free_likelihood_sibling(descendantID);
		free(descendantID);
		i++;
	}
}


void AncestralLikelihoodDiplo(int id_ith, int *siblingID){

	int Phenotype = pedigree[id_ith][4];

	polynomial *polyf,*polym,*polym_1,*polym_23,*polym_4;
	polynomial_individual *tmp,*tmp2_1,*tmp2_23,*tmp2_4,*tmp3;
	int fid_ith,mid_ith;
	int i;

	fid_ith = get_id_ith(pedigree[id_ith][1]);
	mid_ith = get_id_ith(pedigree[id_ith][2]);

	AncestralLikelihood(fid_ith);
	AncestralLikelihood(mid_ith);

	if( *siblingID == 0 ){
		for( i=1; i<5; i++ ){
			tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp);
			prod_scalar(InheritP(1,1,i),&likelihood[fid_ith].poly1,&((*tmp).poly1));
			prod_scalar(InheritP(2,1,i),&likelihood[fid_ith].poly2,&((*tmp).poly2));
			prod_scalar(InheritP(3,1,i),&likelihood[fid_ith].poly3,&((*tmp).poly3));
			prod_scalar(InheritP(4,1,i),&likelihood[fid_ith].poly4,&((*tmp).poly4));
			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			sum(tmp, polyf);
			free_polynomial_individual(tmp);

			tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp);
			prod_scalar(InheritP(1,2,i),&likelihood[mid_ith].poly1,&((*tmp).poly1));
			prod_scalar(InheritP(2,2,i),&likelihood[mid_ith].poly2,&((*tmp).poly2));
			prod_scalar(InheritP(3,2,i),&likelihood[mid_ith].poly3,&((*tmp).poly3));
			prod_scalar(InheritP(4,2,i),&likelihood[mid_ith].poly4,&((*tmp).poly4));
			polym = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polym);
			sum(tmp, polym);
			free_polynomial_individual(tmp);

			if( i==1 ){
				prod_polynomial(polyf,polym,&likelihood[id_ith].poly1);
				penetrance(Phenotype,1,&likelihood[id_ith].poly1);
			}
			else if( i==2 ){
				prod_polynomial(polyf,polym,&likelihood[id_ith].poly2);
				penetrance(Phenotype,2,&likelihood[id_ith].poly2);
			}
			else if( i==3 ){
				prod_polynomial(polyf,polym,&likelihood[id_ith].poly3);
				penetrance(Phenotype,3,&likelihood[id_ith].poly3);
			}
			else if( i==4 ){
				prod_polynomial(polyf,polym,&likelihood[id_ith].poly4);
				penetrance(Phenotype,4,&likelihood[id_ith].poly4);
			}
			free_polynomial(polyf);
			free_polynomial(polym);
		}
	}
	else{
		for( i=1; i<5; i++ ){

			tmp = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp);
			prod_scalar(InheritP(1,2,i),&likelihood[mid_ith].poly1,&((*tmp).poly1));
			prod_scalar(InheritP(2,2,i),&likelihood[mid_ith].poly2,&((*tmp).poly2));
			prod_scalar(InheritP(3,2,i),&likelihood[mid_ith].poly3,&((*tmp).poly3));
			prod_scalar(InheritP(4,2,i),&likelihood[mid_ith].poly4,&((*tmp).poly4));

			tmp2_1 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			tmp2_23 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			tmp2_4 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp2_1);
			allocate_memory_polynomial_individual(tmp2_23);
			allocate_memory_polynomial_individual(tmp2_4);

			prod_polynomial(&((*tmp).poly1),&likelihood_sibling[id_ith].poly1,&((*tmp2_1).poly1));
			prod_polynomial(&((*tmp).poly2),&likelihood_sibling[id_ith].poly2,&((*tmp2_1).poly2));
			prod_polynomial(&((*tmp).poly3),&likelihood_sibling[id_ith].poly2,&((*tmp2_1).poly3));
			prod_polynomial(&((*tmp).poly4),&likelihood_sibling[id_ith].poly3,&((*tmp2_1).poly4));

			prod_polynomial(&((*tmp).poly1),&likelihood_sibling[id_ith].poly4,&((*tmp2_23).poly1));
			prod_polynomial(&((*tmp).poly2),&likelihood_sibling[id_ith].poly5,&((*tmp2_23).poly2));
			prod_polynomial(&((*tmp).poly3),&likelihood_sibling[id_ith].poly5,&((*tmp2_23).poly3));
			prod_polynomial(&((*tmp).poly4),&likelihood_sibling[id_ith].poly6,&((*tmp2_23).poly4));

			prod_polynomial(&((*tmp).poly1),&likelihood_sibling[id_ith].poly7,&((*tmp2_4).poly1));
			prod_polynomial(&((*tmp).poly2),&likelihood_sibling[id_ith].poly8,&((*tmp2_4).poly2));
			prod_polynomial(&((*tmp).poly3),&likelihood_sibling[id_ith].poly8,&((*tmp2_4).poly3));
			prod_polynomial(&((*tmp).poly4),&likelihood_sibling[id_ith].poly9,&((*tmp2_4).poly4));

			polym_1 = (polynomial *)malloc(sizeof(polynomial));
			polym_23 = (polynomial *)malloc(sizeof(polynomial));
			polym_4 = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polym_1);
			allocate_memory_polynomial(polym_23);
			allocate_memory_polynomial(polym_4);
			sum(tmp2_1,polym_1);
			sum(tmp2_23,polym_23);
			sum(tmp2_4,polym_4);

			free_polynomial_individual(tmp2_1);
			free_polynomial_individual(tmp2_23);
			free_polynomial_individual(tmp2_4);

			tmp3 = (polynomial_individual *)malloc(sizeof(polynomial_individual));
			allocate_memory_polynomial_individual(tmp3);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(1,1,i),&likelihood[fid_ith].poly1,polyf);

			prod_polynomial(polyf,polym_1,&((*tmp3).poly1));
			free_polynomial(polyf);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(2,1,i),&likelihood[fid_ith].poly2,polyf);
			prod_polynomial(polyf,polym_23,&((*tmp3).poly2));
			free_polynomial(polyf);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(3,1,i),&likelihood[fid_ith].poly3,polyf);
			prod_polynomial(polyf,polym_23,&((*tmp3).poly3));
			free_polynomial(polyf);

			polyf = (polynomial *)malloc(sizeof(polynomial));
			allocate_memory_polynomial(polyf);
			prod_scalar(InheritP(4,1,i),&likelihood[fid_ith].poly4,polyf);
			prod_polynomial(polyf,polym_4,&((*tmp3).poly4));
			free_polynomial(polyf);

			if( i==1 ){
				sum(tmp3,&likelihood[id_ith].poly1);
				penetrance(Phenotype,1,&likelihood[id_ith].poly1);
			}
			else if( i==2 ){
				sum(tmp3,&likelihood[id_ith].poly2);
				penetrance(Phenotype,2,&likelihood[id_ith].poly2);
			}
			else if( i==3 ){
				sum(tmp3,&likelihood[id_ith].poly3);
				penetrance(Phenotype,3,&likelihood[id_ith].poly3);
			}
			else if( i==4 ){
				sum(tmp3,&likelihood[id_ith].poly4);
				penetrance(Phenotype,4,&likelihood[id_ith].poly4);
			}

			free_polynomial_individual(tmp);
			free_polynomial_individual(tmp3);
			free_polynomial(polym_1);
			free_polynomial(polym_23);
			free_polynomial(polym_4);
		}
	}
}


void AncestralLikelihood(int id_ith){

	int id = pedigree[id_ith][0];

	if( likelihood[id_ith].poly1.size == 0 ){
		if( is_founder(id) == 1 ){
			int Phenotype = pedigree[id_ith][4];
			Xi_cal_founderP_3paras(&likelihood[id_ith], Phenotype);
		}
		else{
			int *siblingID = (int *)calloc(NIndividual, sizeof(int));
			sibling(id, siblingID);

			if( *siblingID != 0 ){
				calculate_likelihood_sibling(&likelihood_sibling[id_ith], siblingID);
			}

			AncestralLikelihoodDiplo(id_ith, siblingID);

		}
	}
}







/*evaluation of MLE*/
void eval_grr(double *penetrance, double *poly, int *powers, int *max_power, int *length, double *result)
{

	int i,j;
	int matrix[3][*length];
	int tmp;
	int a;
	for(i=0;i<*length;i++){
		tmp = *(powers+i);
		for(j=2;j>=1;j--){
//			matrix[j][i] = tmp / (int)pow(*max_power, j);
			a = tmp % (int)pow(*max_power, j);
//			if( a == 0 ){
//				matrix[j][i] = (int)(tmp / pow(*max_power, j) + 0.1);
//			}
//			else{
				matrix[j][i] = tmp / (int)pow(*max_power, j);
//			}
			tmp = a;
		}
		matrix[0][i] = tmp;
	}
	double value=0.0;
	double tmp1,tmp2;
	int k;
	int x = 0;

	while( x < *length ){
		k = matrix[2][x];
		tmp1 = 0.0;
		while( (x < *length) && k == matrix[2][x] ){
			j = matrix[1][x];
			tmp2 = 0.0;
			while( ((x < *length) && j == matrix[1][x]) && k == matrix[2][x] ){
				tmp2 += *(poly+x) * pow(*penetrance, matrix[0][x]);
				x++;
			}
			tmp1 += tmp2 * pow(*(penetrance+1), j);
		}
		value += tmp1 * pow(*(penetrance+2), k);
	}
//	*result = log(value);



/* gradient */

	double coef_alpha[*length],coef_beta[*length],coef_gamma[*length];
	int dalpha[*length],dbeta[*length],dgamma[*length];

	for(i=0;i<*length;i++){
		coef_alpha[i] = *(poly+i) * matrix[0][i];
		coef_beta[i] = *(poly+i) * matrix[1][i];
		coef_gamma[i] = *(poly+i) * matrix[2][i];
		if(matrix[0][i]==0){
			dalpha[i]=0;
		}
		else{
			dalpha[i]=matrix[0][i]-1;
		}
		if(matrix[1][i]==0){
			dbeta[i]=0;
		}
		else{
			dbeta[i]=matrix[1][i]-1;
		}
		if(matrix[2][i]==0){
			dgamma[i]=0;
		}
		else{
			dgamma[i]=matrix[2][i]-1;
		}
	}

	double valuealpha=0.0;
	double valuebeta=0.0;
	double valuegamma=0.0;

/* gradient(alpha) */
	x=0;
	while( x < *length ){
		k = matrix[2][x];
		tmp1 = 0.0;
		while( (x < *length) && k == matrix[2][x] ){
			j = matrix[1][x];
			tmp2 = 0.0;
			while( ((x < *length) && j == matrix[1][x]) && k == matrix[2][x] ){
				tmp2 += *(coef_alpha+x) * pow(*penetrance, *(dalpha+x));
				x++;
			}
			tmp1 += tmp2 * pow(*(penetrance+1), j);
		}
		valuealpha += tmp1 * pow(*(penetrance+2), k);
	}
/* gradient(beta) */
	x=0;
	while( x < *length ){
		k = matrix[2][x];
		tmp1 = 0.0;
		while( (x < *length) && k == matrix[2][x] ){
			j = *(dbeta+x);
			tmp2 = 0.0;
			while( ((x < *length) && j == *(dbeta+x)) && k == matrix[2][x] ){
				tmp2 += *(coef_beta+x) * pow(*penetrance, matrix[0][x]);
				x++;
			}
			tmp1 += tmp2 * pow(*(penetrance+1), j);
		}
		valuebeta += tmp1 * pow(*(penetrance+2), k);
	}

/* gradient(gamma) */
	x=0;
	while( x < *length ){
		k = *(dgamma+x);
		tmp1 = 0.0;
		while( (x < *length) && k == *(dgamma+x) ){
			j = matrix[1][x];
			tmp2 = 0.0;
			while( ((x < *length) && j == matrix[1][x]) && k == *(dgamma+x) ){
				tmp2 += *(coef_gamma+x) * pow(*penetrance, matrix[0][x]);
				x++;
			}
			tmp1 += tmp2 * pow(*(penetrance+1), j);
		}
		valuegamma += tmp1 * pow(*(penetrance+2), k);
	}

	*(result) = valuealpha / value;
	*(result+1) = valuebeta / value;
	*(result+2) = valuegamma / value;


}


void eval_fr(double *penetrance, double *poly, int *powers, int *max_power, int *length, double *result)
{

	int i,j;
	int matrix[3][*length];
	int tmp;
	int a;
	for(i=0;i<*length;i++){
		tmp = *(powers+i);
		for(j=2;j>=1;j--){
//			matrix[j][i] = tmp / (int)pow(*max_power, j);
			a = tmp % (int)pow(*max_power, j);
//			if( a == 0 ){
//				matrix[j][i] = (int)(tmp / pow(*max_power, j) + 0.1);
//			}
//			else{
				matrix[j][i] = tmp / (int)pow(*max_power, j);
//			}
			tmp = a;
		}
		matrix[0][i] = tmp;
	}

	double value=0.0;
	double tmp1,tmp2;
	int k;
	int x = 0;

	while( x < *length ){
		k = matrix[2][x];
		tmp1 = 0.0;
		while( (x < *length) && k == matrix[2][x] ){
			j = matrix[1][x];
			tmp2 = 0.0;
			while( ((x < *length) && j == matrix[1][x]) && k == matrix[2][x] ){
				tmp2 += *(poly+x) * pow(*penetrance, matrix[0][x]);
				x++;
			}
			tmp1 += tmp2 * pow(*(penetrance+1), j);
		}
		value += tmp1 * pow(*(penetrance+2), k);
	}

	*result = log(value);

}




