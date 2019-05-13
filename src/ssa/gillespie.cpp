//
//  main.cpp
//  gillespie
//

// simulates trajectories of a ctmc based on gillespies direct method

// compile with: mex gillespie.cpp

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <ctime>
#include <mex.h>

struct system_model {
    size_t num_species;
    size_t num_reactions;
    size_t num_elements;
    double *Pre;
    double *Post;
    double *S;
    double *state;
    double *rates;
    double t_min;
    double t_max;
    int seed;
};

void gillespie (std::vector<double> &time, std::vector<double> &state_history, std::vector<double> &events, std::vector<double> &total_propensity, system_model &sys) ;
void update_propensity(std::vector<double> &propensity, std::vector<double> &state, system_model &sys);
inline double comb_factor (int n, int k);
bool next_reaction (double *delta_t, int *index, system_model &sys, std::vector<double> &propensity, double rand_1, double rand_2);
void update_state(std::vector<double> &state,int index,system_model &sys);
void update_history(std::vector<double> &state_history,std::vector<double> &state, system_model &sys);
void update_total_propensity(std::vector<double> &total_propensity, const std::vector<double> &propensity, system_model &sys, double delta_t );

// int main(int argc, const char * argv[]) {
//     
//     // set up model
//     size_t size_estimate = 1000;
//     std::vector<double> Pre = {0,1};
//     std::vector<double> Post = {1,0};
//     std::vector<double> S = {1,-1};
//     std::vector<double> state = {0,0};
//     std::vector<double> rates = {1.0,0.01};
//     system_model sys;
//     sys.num_species = 1;
//     sys.num_reactions = 2;
//     sys.num_elements = 2;
//     sys.Pre = Pre.data();
//     sys.Post = Post.data();
//     sys.S = S.data();
//     sys.state = state.data();
//     sys.rates = rates.data();
//     sys.t_min = 0.0;
//     sys.t_max = 1e4;
//     sys.seed = std::time(NULL);
//     
//     // set up output
//     std::vector<double> time;
//     time.reserve(size_estimate);
//     std::vector<double> state_history;
//     time.reserve(size_estimate*sys.num_species);
//     std::vector<double> events(sys.num_reactions);
//     std::vector<double> total_propensity(sys.num_reactions);
//     
//     // perform calculation
//     gillespie(time,state_history,events,total_propensity,sys);
//     
//     // test output
//     for (size_t i = 0; i < events.size(); i++) {
//         std::cout << i << " " << events[i] << std::endl;
//     }
//     
// }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* This function expects one input argument: A structure defining the system
     * containing the fields num_species, num_reactions, Pre, Post, S, rates, t_min, t_max, 
     * initial, seed
     */
    // check number of output arguments
    if (nrhs != 1)
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin","Expected one input argument!");
    if (nlhs > 4)
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout","Too many output arguments!");
    
    // The input argument is the system struct. This must be transformed into something that c++ can work with
    mxArray *mxArr;
    
    // Extract number of species and reactions
    mxArr = mxGetField(prhs[0],0,"num_species");
    unsigned num_species = (unsigned)(*mxGetPr(mxArr));
    mxArr = mxGetField(prhs[0],0,"num_reactions");
    unsigned num_reactions = (unsigned)(*mxGetPr(mxArr));
    
    // Extract Pre matrix
    mxArr = mxGetField(prhs[0],0,"Pre");
    double *Pre = mxGetPr(mxArr);
    
    // Extract Post matrix
    mxArr = mxGetField(prhs[0],0,"Post");
    double *Post = mxGetPr(mxArr);
    
    // Extract S matrix
    mxArr = mxGetField(prhs[0],0,"S");
    double *S = mxGetPr(mxArr);
    
    // Extract start state
    mxArr = mxGetField(prhs[0],0,"initial");
    double *state = mxGetPr(mxArr);
    
    // Extract rate constant
    mxArr = mxGetField(prhs[0],0,"rates");
    double *rates = mxGetPr(mxArr);
    
    // Extract start and beginning times
    mxArr = mxGetField(prhs[0],0,"t_min");
    double t_min = *mxGetPr(mxArr);
    mxArr = mxGetField(prhs[0],0,"t_max");
    double t_max = *mxGetPr(mxArr);
    
    // Extract rng seed
    mxArr = mxGetField(prhs[0],0,"seed");
    unsigned seed = (unsigned)(*mxGetPr(mxArr));
    
    // set up objects required for the calculation
    system_model sys;
    sys.num_species = num_species;
    sys.num_reactions = num_reactions;
    sys.num_elements = num_reactions*num_species;
    sys.Pre = Pre;
    sys.Post = Post;
    sys.S = S;
    sys.state = state;
    sys.rates = rates;
    sys.t_min = t_min;
    sys.t_max = t_max;
    sys.seed = seed;
    
    // set up output
    size_t size_estimate = 1000;
    std::vector<double> time;
    time.reserve(size_estimate);
    std::vector<double> state_history;
    state_history.reserve(size_estimate*sys.num_species);
    std::vector<double> events(sys.num_reactions);
    std::vector<double> total_propensity(sys.num_reactions);
    
    // perform calculation
    gillespie(time,state_history,events,total_propensity,sys);
    
    // create output matrix for time
    plhs[0] = mxCreateDoubleMatrix(1,time.size(),mxREAL);
    double *time_out = mxGetPr(plhs[0]);
    for (int i = 0; i < time.size(); i++) {
        time_out[i] = time[i];
    }
    
    // create output matrix for the state
    plhs[1] = mxCreateDoubleMatrix(num_species,state_history.size()/num_species,mxREAL);
    double *state_out = mxGetPr(plhs[1]);
    for (int i = 0; i < state_history.size(); i++) {
        state_out[i] = state_history[i];
    }
    
    
}

void gillespie (std::vector<double> &time, std::vector<double> &state_history, std::vector<double> &events, std::vector<double> &total_propensity, system_model &sys) {
    // preparations
    double t = sys.t_min;
    double t_max = sys.t_max;
    double delta_t = 0.0;
    int index = 0;
    std::vector<double> state(sys.state,sys.state+sys.num_species);
    std::vector<double> propensity(sys.num_reactions);
    std::mt19937 rng(sys.seed);
    std::uniform_real_distribution<double> U;
    // store initals in history vector
    time.push_back(t);
    update_history(state_history,state,sys);
    // sample path
    while (t < t_max) {
        // determine next event
        update_propensity(propensity,state,sys);
        bool success = next_reaction(&delta_t,&index,sys,propensity,U(rng),U(rng));
        if (not success)
            break;
        // update system
        t += delta_t;
        update_state(state,index,sys);
        // update output statistics
        events[index] += 1.0;
        time.push_back(t);
        update_history(state_history,state,sys);
        update_total_propensity(total_propensity,propensity,sys,delta_t);
    }
    return;
}

void update_propensity(std::vector<double> &propensity, std::vector<double> &state, system_model &sys) {
    // initialise to one
    for (int i = 0; i < sys.num_reactions; i++) {
        propensity[i] = 1.0;
    }
    // calculate the stoichiometric factors
    for (int i = 0; i < sys.num_elements; i++) {
        if ( sys.Pre[i] > 0 ) {
            size_t species = i/sys.num_reactions;
            size_t reaction = i%sys.num_reactions;
            propensity[reaction] *= comb_factor(state[species],sys.Pre[i]);
        }
    }
    return;
}

inline double comb_factor (int n, int k) {
    double res;
    if (n < k)
        res = 0.0;
    else {
        res = 1.0;
        for (int i = 0; i < k; i++) {
            res *= n-i;
        }
    }
    return(res);
}

bool next_reaction (double *delta_t, int *index, system_model &sys, std::vector<double> &propensity, double rand_1, double rand_2) {
    /* Calculates reaction times for all channels. The mimimum time and the corresponding index are saved in delta_t and index. */
    // calculate the reaction hazards
    std::vector<double> hazard(sys.num_reactions);
    double total_hazard = 0.0;
    for (size_t i = 0; i < sys.num_reactions; i++ ) {
        total_hazard += sys.rates[i]*propensity[i];
        hazard[i] = total_hazard;
    }
    // calculate reaction time
    *delta_t = -std::log(rand_1)/total_hazard;
    // sample random event from the individual hazards
    rand_2 *= total_hazard;
    *index = 0;
    while ( hazard[*index] < rand_2) {
        (*index)++;
    }
    return( total_hazard>0.0);
}

void update_state(std::vector<double> &state,int index,system_model &sys) {
    for (size_t i = 0; i < sys.num_species; i++ ) {
        size_t ind = i*sys.num_reactions+index;
        state[i] += sys.S[ind];
    }
    return;
}

void update_history(std::vector<double> &state_history,std::vector<double> &state, system_model &sys) {
    for (size_t i = 0; i < sys.num_species; i++) {
        state_history.push_back(state[i]);
    }
    return;
}

void update_total_propensity(std::vector<double> &total_propensity, const std::vector<double> &propensity, system_model &sys, double delta_t ) {
    for (size_t i = 0; i < sys.num_reactions; i++) {
        total_propensity[i] += delta_t*propensity[i];
    }
    return;
}



