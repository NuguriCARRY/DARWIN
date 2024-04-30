#include "afl_darwin.h"
#include "rand.h"
#include <stdio.h>

unsigned kAflNumMutationOperators;

#ifdef REAL_VALUED
	#define VALUE_TYPE double
#else
	#define VALUE_TYPE bool
#endif

#ifndef _nParents
#define _nParents 5u
#endif

#ifndef _lambda
#define _lambda 4u
#endif

unsigned int *_nextToEvaluate;  // next child to evaluate
unsigned int *_currentParent;   // current parent being mutated
unsigned int *_bestSoFar;       // best so far parent
unsigned int *_bestSoFarChild;  // best so far child (for a single parent)
int **_parentsFitness;
int **_childrenFitness;
unsigned int variation_intensity_best = 0;


VALUE_TYPE ***_parents;
VALUE_TYPE ***_children;
VALUE_TYPE **_currentRelativeProb;

/*
 * Initialize data structures for DARWIN
 * @intensity: number of different initial seeds
 * @nr_mutations: number of mutation operators
*/
void DARWIN_init(uint32_t intensity, unsigned nr_mutations) {
	// init RNG
	rand_init();

	kAflNumMutationOperators = nr_mutations;

	// initialize opt alg data structures
	_nextToEvaluate = (unsigned int *) malloc(intensity * sizeof(unsigned int));
	_currentParent = (unsigned int *) malloc(intensity * sizeof(unsigned int));
	_bestSoFar = (unsigned int *) malloc(intensity * sizeof(unsigned int));         // best so far parent
	_bestSoFarChild = (unsigned int *) malloc(intensity * sizeof(unsigned int));    // best so far child (for a single parent)
	_parentsFitness = (int **) malloc(intensity * sizeof(int *));
	_childrenFitness = (int **) malloc(intensity * sizeof(int *));

	_parents = (VALUE_TYPE ***) malloc(intensity * sizeof(VALUE_TYPE **));
	_children = (VALUE_TYPE ***) malloc(intensity * sizeof(VALUE_TYPE **));
	
	_currentRelativeProb = (VALUE_TYPE **) malloc(intensity * sizeof(VALUE_TYPE *));

	printf("nParents: %u, lambda %u\n", _nParents, _lambda);

	for (int cur_intensity = 0; cur_intensity < intensity; cur_intensity++) {
		_nextToEvaluate[cur_intensity] = 0;
		_bestSoFar[cur_intensity] = 0;
		_bestSoFarChild[cur_intensity] = 0;
		_currentParent[cur_intensity] = 0;

		_parentsFitness[cur_intensity] = malloc(_nParents * sizeof(int));
		_childrenFitness[cur_intensity] = malloc(_lambda * sizeof(int));

		_children[cur_intensity] = (VALUE_TYPE **) malloc(_lambda * sizeof(VALUE_TYPE *));
		for (int i = 0; i < _lambda; i++)
			_children[cur_intensity][i] = malloc(kAflNumMutationOperators * sizeof(VALUE_TYPE));

		_parents[cur_intensity] = (VALUE_TYPE **) malloc(_nParents * sizeof(VALUE_TYPE *));
		for (int i = 0; i < _nParents; i++)
			_parents[cur_intensity][i] = malloc(kAflNumMutationOperators * sizeof(VALUE_TYPE));

		memset(_parentsFitness[cur_intensity], 0, _nParents);
		memset(_childrenFitness[cur_intensity], 0, _lambda);

		// initial random values for the parents and first child individual
		for(uint i = 0; i < _nParents; i++) {
			for(uint j = 0; j < kAflNumMutationOperators; j++) {
				// TODO: optionally replace with an alternate randomizer
				_parents[cur_intensity][i][j] = rand() > (RAND_MAX / 2);
			}
		}

		for(uint i = 0; i < kAflNumMutationOperators; i++) {
			_children[cur_intensity][_nextToEvaluate[cur_intensity]][i] = rand() > (RAND_MAX / 2);
		}

		_currentRelativeProb[cur_intensity] = _children[cur_intensity][_nextToEvaluate[cur_intensity]];
	}
}

/*
 * Choose an AFL mutation operator
 * @cur_intensity: cur_intensity to select per-cur_intensity vector
*/
int DARWIN_SelectOperator(uint32_t cur_intensity)
{
#ifdef REAL_VALUED
	double *v = _currentRelativeProb[cur_intensity];
	uint v_size = kAflNumMutationOperators;
	double cumulative[kA flNumMutationOperators];
	for (uint i = 0; i < kAflNumMutationOperators; i++)	
		cumulative[i] = 0.0;

	double minVal = 0.0;
	double maxVal = 1.0;
	for (uint i = 0; i < kAflNumMutationOperators; i++) {
		if (v[kAflNumMutationOperators] < minVal) {
			minVal = v[kAflNumMutationOperators];
		} else {
			if (v[kAflNumMutationOperators] > maxVal) {
				maxVal = v[kAflNumMutationOperators];
			}
		}
	}

	cumulative[0] = 1 + (kAflNumMutationOperators - 1) * (v[0] - minVal) / (maxVal - minVal); // selection pressure is kAflNumMutationOperators
	for(uint i = 1; i < kAflNumMutationOperators; i++) {
		cumulative[i] = 1 + (kAflNumMutationOperators - 1) * (v[i] - minVal)/(maxVal - minVal);
		cumulative[i] += cumulative[i - 1];
	}

	double rngNormVal = rand_32_double();

	double randVal = rngNormVal * cumulative[kAflNumMutationOperators - 1];

	int chosen = 0;
	// didn't bother with binary search since negligible compared to other modules
	while(cumulative[chosen] < randVal)
		chosen++;
	return chosen;
#else
	// baseline:
	bool *v = (bool *) _currentRelativeProb[cur_intensity];

	// select a random mutation operator with flag == true
	// revert to random if all false
	uint operatorId = rand_32_int(kAflNumMutationOperators);
	uint nTries = 0;
	while (v[operatorId] == false && nTries < kAflNumMutationOperators) {
		nTries++;
		operatorId = (operatorId + 1) % kAflNumMutationOperators;
	}
    	return operatorId;
#endif
}


/*
 * Report feedback to DARWIN
 * @cur_intensity: cur_intensity to attribute this to
 * @numPaths: number of new paths
*/
uint32_t DARWIN_NotifyFeedback(uint32_t cur_intensity, unsigned numPaths)
{
	// update this candidate solution fitness
	_childrenFitness[cur_intensity][_nextToEvaluate[cur_intensity]] = numPaths;

	// update if new best child found
	if(_childrenFitness[cur_intensity][_nextToEvaluate[cur_intensity]] > _childrenFitness[cur_intensity][_bestSoFarChild[cur_intensity]]) {
		_bestSoFarChild[cur_intensity] = _nextToEvaluate[cur_intensity];
	}
	// move to next child candidate
	_nextToEvaluate[cur_intensity]++;
	
	// if all children evaluated
	if(_nextToEvaluate[cur_intensity] == _lambda) {
		// set best child as future parent (only if better than the parent)
		if(_childrenFitness[cur_intensity][_bestSoFarChild[cur_intensity]] > _parentsFitness[cur_intensity][_currentParent[cur_intensity]]) {
			for(uint i = 0; i < kAflNumMutationOperators; i++) {
				_parents[cur_intensity][_currentParent[cur_intensity]][i] = _children[cur_intensity][_bestSoFarChild[cur_intensity]][i];
			}

			_parentsFitness[cur_intensity][_currentParent[cur_intensity]] = _childrenFitness[cur_intensity][_bestSoFarChild[cur_intensity]];

			variation_intensity_best = cur_intensity;
		}
		// update best parent solution if needed
		if(_parentsFitness[cur_intensity][_currentParent[cur_intensity]] > _parentsFitness[cur_intensity][_bestSoFar[cur_intensity]]) {
			_bestSoFar[cur_intensity] = _currentParent[cur_intensity];
		}

		// move to next parent (or return to first)
		_currentParent[cur_intensity] = (_currentParent[cur_intensity] + 1) % _nParents;

		// reset indices
		_bestSoFarChild[cur_intensity] = 0;
		_nextToEvaluate[cur_intensity] = 0;

		// reset children scores
		memset(_childrenFitness[cur_intensity], 0, _lambda); 
	}
	
	// if there are children to evaluate, generate new candidate and return
	if(_nextToEvaluate[cur_intensity] < _lambda) {
		_currentRelativeProb[cur_intensity] = &(_children[cur_intensity][_nextToEvaluate[cur_intensity]][0]);

		for(uint i = 1; i < kAflNumMutationOperators; i++) {
			_currentRelativeProb[cur_intensity][i] = _parents[cur_intensity][_currentParent[cur_intensity]][i];
		}
		
		// select a single random gene and invert
		int randomGene = rand_32_int(kAflNumMutationOperators);
#ifdef REAL_VALUED
		double mutation = (rand_32_double_gauss() - 0.5) / 4;
		_currentRelativeProb[cur_intensity][randomGene] += mutation;
		if(_currentRelativeProb[cur_intensity][randomGene] < 0)
			_currentRelativeProb[cur_intensity][randomGene] = 0;
#else		
		_currentRelativeProb[cur_intensity][randomGene] = !_currentRelativeProb[cur_intensity][randomGene];
#endif
	}
	return variation_intensity_best;
}

/*
 * Get the best parent solution so far as a vector (hardcoded to max mutation operators in AFL)
 * @cur_intensity: cur_intensity to attribute this to
*/
uint32_t DARWIN_get_parent_repr(uint32_t cur_intensity) {
	uint32_t value;
	bool *_currentParentBool = (bool *) (_parents[cur_intensity][_bestSoFar[cur_intensity]]); 
	for (int i = 0; i < 15; i++) {
		printf("%d: %s\n", i,_currentParentBool[i] ? "true" : "false");
	}
	/* encode parent in a number */
	for (int i = 14; i >= 0; i--) {
		value |= (_parents[cur_intensity][_bestSoFar[cur_intensity]][i]) ? (1 << i) : 0;
	}
	printf("val: %d\n", value);
	return value;
}
