#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include "Ngram.h"
#include "File.h"
#include "Vocab.h"
#include <limits>
#include <unordered_map>

#define MAX_T 1500
#define MAX_Q 15000


using namespace std;
unsigned cand_cnts[MAX_Q]={0};
char wdlookup[MAX_Q][MAX_Q][3]={0};

double inf = std::numeric_limits<double>::infinity();
unordered_map<int, int> lmap;

void disambig(VocabString* outseq, VocabString* sentence, int numWords, Vocab& voc, Ngram& lm);
void disambig3(VocabString* outseq, VocabString* sentence, int numWords, Vocab& voc, Ngram& lm);

int main(int argc, char** argv){
	if(argc != 9){
		cerr << "Input format: ./disambig -text [text] -map [map] -lm [lm] -order [ord]" << endl;
		return 1;
	}
	char* textfile = argv[2];
	char* mapfile = argv[4];
	char* lmfile = argv[6];
	int ngram_order = argv[8][0]-'0';

	// initialize language model
	Vocab voc;
	Ngram lm( voc, ngram_order );
	File lmFile( lmfile, "r" );
    lm.read(lmFile);
    lmFile.close();


	VocabString sentence[MAX_T];
    File file(mapfile, "r");
    char *line;    
	unsigned lookupcnt = 0;
    while ((line = file.getline())) {	
		unsigned numWords = Vocab::parseWords(line, sentence, MAX_T);

		// set counts of candidates
		cand_cnts[lookupcnt] = numWords - 1;

		// add all candidates
		sentence[numWords] = 0;
		for(int i = 0; i < numWords; i++){
			strncpy(wdlookup[lookupcnt][i], sentence[i], 2);
		}

		// create mapping
		int val = (-sentence[0][0])*256 + sentence[0][1] + 128;
		lmap[val] = lookupcnt;

		// increment entries
		lookupcnt++;
	}
	file.close();


    
	File file2(textfile, "r");
	File file3("result.txt", "w");
    while ((line = file2.getline())) {
		unsigned numWords = Vocab::parseWords(line, sentence, MAX_T);
		// run viterbi here
		
		VocabString out[MAX_T];
		if(ngram_order==2)
			disambig(out, sentence, numWords, voc, lm);
		else
			disambig3(out, sentence, numWords, voc, lm);

		// print best sentence to results.
		cout << "<s> ";
		for(int t = 0; t < numWords; t++)
			cout << out[t] << " ";
		cout << "</s>" << endl;
	}
	file2.close();
}


int getMyIndex(VocabString c){
	return lmap.find((-c[0])*256 + c[1] + 128)->second;
}


// Get P(W1) -- unigram
double getUnigramProb(const char *w1, Vocab& voc, Ngram& lm)
{
    VocabIndex wid1 = voc.getIndex(w1);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);

    VocabIndex emptyContext[] = { Vocab_None };
    return lm.wordProb(wid1, emptyContext);
}

// Get P(W2 | W1) -- bigram
double getBigramProb(const char *w1, const char *w2, Vocab& voc, Ngram& lm)
{
    VocabIndex wid1 = voc.getIndex(w1);
    VocabIndex wid2 = voc.getIndex(w2);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);
    if(wid2 == Vocab_None)  //OOV
        wid2 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { wid1, Vocab_None };
    return lm.wordProb( wid2, context);
}

// Get P(W3 | W1, W2) -- trigram
double getTrigramProb(const char *w1, const char *w2, const char *w3, Vocab& voc, Ngram& lm)
{
    VocabIndex wid1 = voc.getIndex(w1);
    VocabIndex wid2 = voc.getIndex(w2);
    VocabIndex wid3 = voc.getIndex(w3);

    if(wid1 == Vocab_None)  //OOV
        wid1 = voc.getIndex(Vocab_Unknown);
    if(wid2 == Vocab_None)  //OOV
        wid2 = voc.getIndex(Vocab_Unknown);
    if(wid3 == Vocab_None)  //OOV
        wid3 = voc.getIndex(Vocab_Unknown);

    VocabIndex context[] = { wid2, wid1, Vocab_None };
    return lm.wordProb( wid3, context );
}


typedef struct v{
	VocabString wd;
	double delta;
	struct v* parent;
} VNode;

void disambig(VocabString* outseq, VocabString* sentence, int numWords, Vocab& voc, Ngram& lm){
	
	int nodes = 0, max_q = 0;
	int c_count;
	int index;
	int indeces[numWords];
	for(int t = 0; t < numWords; t++){
		index = getMyIndex(sentence[t]);
		indeces[t] = index;
		c_count = cand_cnts[index];
		nodes += c_count;
		if(c_count > max_q)
			max_q = c_count;
	}

	int usage = 0, flip = 1;
	VNode* pool = new VNode[nodes];
	VNode** A = new VNode*[max_q];
	VNode** B = new VNode*[max_q];



	// t = 1: P(W1=qi)
	index = indeces[0];
	c_count = cand_cnts[index];
	for(int j = 0; j < c_count; j++){
		VocabString jth = wdlookup[index][j+1];

		VNode* ptr = &(pool[usage++]);
		ptr->wd = jth;
		ptr->delta = getUnigramProb(jth, voc, lm);
		ptr->parent = NULL;
		A[j] = ptr;
	}
	
	// t > 1: P(Wt=qj | W(t-1)=qi)
	VNode **context, **tail;
	int index0, c_count0;
	for(int t = 1; t < numWords; t++){
		context = flip ? A : B;
		tail = flip ? B : A;

		index = indeces[t];
		c_count = cand_cnts[index];
		index0 = indeces[t-1];
		c_count0 = cand_cnts[index0];

		for(int j = 0; j < c_count; j++){
			VocabString jth = wdlookup[index][j+1];


			int bestind = -1;
			double bestprob = 0.;
			for(int i = 0; i < c_count0; i++){
				VocabString prev = context[i]->wd;
				double thisprob = context[i]->delta+getBigramProb(prev, jth, voc, lm);
				if(bestind==-1 || thisprob>bestprob){
					bestind = i;
					bestprob = thisprob;
				}
			}

			VNode* ptr = &(pool[usage++]);
			ptr->wd = jth;
			ptr->delta = bestprob;
			ptr->parent = context[bestind];
			tail[j] = ptr;
		}
		flip = !flip;
	}

	int bestend = -1;
	double bestendprob = 0.;
	for(int i = 0; i < c_count; i++){
		double thisprob = tail[i]->delta;
		if(bestend==-1 || thisprob>bestendprob){
			bestend = i;
			bestendprob = thisprob;
		}
	}

	VNode* ptr2 = tail[bestend];
	for(int t = numWords-1; t>=0; t--){
		outseq[t] = ptr2->wd;
		ptr2 = ptr2->parent;
	}

	delete[] pool;
	delete[] A;
	delete[] B;
}

void disambig3(VocabString* outseq, VocabString* sentence, int numWords, Vocab& voc, Ngram& lm){

	unsigned v_index[numWords];

	unsigned width = 0;
	for(int t = 0; t < numWords; t++){
		int idx = getMyIndex(sentence[t]);
		if (cand_cnts[idx] > width)
			width = cand_cnts[idx];
		v_index[t] = idx;
	}

 	double ***delta = new double **[2];
    for (int i = 0; i < 2; i++)
    {
        delta[i] = new double *[width];
        for (int j = 0; j < width; j++)
            delta[i][j] = new double [width];
    }

    int ***phi_k = new int **[numWords];
    for (int i = 0; i < numWords; i++)
    {
        phi_k[i] = new int *[width];
        for (int j = 0; j < width; j++)
            phi_k[i][j] = new int [width];
    }
	
	//assume k >> j >> i in sentence

	// t = 1
	unsigned i_index = v_index[0];
	unsigned i_count = cand_cnts[i_index];

	for(int i = 0; i < i_count; i++){
		VocabString qi = wdlookup[i_index][i+1];
		double uni = getUnigramProb(qi, voc, lm);

		for(int j = 0; j < width; j++)
			delta[0][i][j] = uni;
	}

	// t = 2
	int j_index = v_index[0];
	int j_count = cand_cnts[j_index];

	i_index = v_index[1];
	i_count = cand_cnts[i_index];

	for(int i = 0; i < i_count; i++){
		VocabString qi = wdlookup[i_index][i+1];
		
		for(int j = 0; j < j_count; j++){
			VocabString qj = wdlookup[j_index][j+1];
			
			double bi = getBigramProb(qj, qi, voc, lm);
			
			delta[1][i][j] = bi + delta[0][j][0];
		}
	}
	// t > 2
	int k_index, k_count;
	int fl = 0;
	for(int t = 2; t < numWords; t++){
		k_index = v_index[t-2];
		k_count = cand_cnts[k_index];

		j_index = v_index[t-1];
		j_count = cand_cnts[j_index];

		i_index = v_index[t];
		i_count = cand_cnts[i_index];
		
		for(int i = 0; i < i_count; i++){
			VocabString qi = wdlookup[i_index][i+1];
			for(int j = 0; j < j_count; j++){
				VocabString qj = wdlookup[j_index][j+1];

				int best_k = -1;
				double bestprob = -inf;
				for(int k = 0; k < k_count; k++){
					VocabString qk = wdlookup[k_index][k+1];
					double trig = getTrigramProb(qk, qj, qi, voc, lm);
					double thisprob = trig + delta[!fl][j][k];
					if(best_k==-1 || thisprob>bestprob){
						best_k = k;
						bestprob = thisprob;
					}
				}
				delta[fl][i][j] = bestprob;
				phi_k[t][i][j] = best_k;
			}
		}
		fl = !fl;
	}


	j_index = v_index[numWords-2];
	j_count = cand_cnts[j_index];

	i_index = v_index[numWords-1];
	i_count = cand_cnts[i_index];

	int best_i = -1, best_j, best_k = -1;
	double bestprob=-inf; 
	for(int i = 0; i < i_count; i++){		
		for(int j = 0; j < j_count; j++){
			double thisprob = delta[!fl][i][j];
			if(best_i==-1 || thisprob>bestprob){
				best_i = i;
				best_j = j;
				best_k = phi_k[numWords-1][i][j];
				bestprob = thisprob;
			}
		}
	}



	k_index = v_index[numWords-3];

	outseq[numWords-1] = wdlookup[i_index][best_i+1];
	outseq[numWords-2] = wdlookup[j_index][best_j+1];
	outseq[numWords-3] = wdlookup[k_index][best_k+1];

	for(int t = numWords-4; t >= 0; t--){
		i_index = v_index[t+2];
		j_index = v_index[t+1];
		k_index = v_index[t];

		best_i = best_j;
		best_j = best_k;
		best_k = phi_k[t+2][best_i][best_j];

		outseq[t] = wdlookup[k_index][best_k+1];
	}

	// k_index = v_index[numWords-3];

	// outseq[numWords-1] = wdlookup[i_index][best_i+1];
	// outseq[numWords-2] = wdlookup[j_index][best_j+1];
	// outseq[numWords-3] = wdlookup[k_index][best_k+1];

	// for(int t = numWords-4; t >= 0; t--){
	// 	i_index = v_index[t+2];
	// 	j_index = v_index[t+1];
	// 	k_index = v_index[t];

	// 	best_i = best_j;
	// 	best_j = best_k;
	// 	best_k = phi_k[t+2][best_i][best_j];

	// 	outseq[t] = wdlookup[k_index][best_k+1];
	// }


	// dealloc
	for (int i = 0; i < 2; i++)
	{
	    for (int j = 0; j < width; j++)
	        delete[] delta[i][j];
	    delete[] delta[i];
	}
	delete[] delta;
	for (int i = 0; i < numWords; i++)
	{
	    for (int j = 0; j < width; j++)
	        delete[] phi_k[i][j];
	    delete[] phi_k[i];
	}
	delete[] phi_k;
}