#include <complex>
#include <vector>


#define OMEGA_FUN(name) __ ## opr_ww_i1 ## _MOD_ ## name


extern "C" {
    void OMEGA_FUN(new_event)(double*);
    void OMEGA_FUN(print_momenta)(double*);
    double _Complex OMEGA_FUN(get_amplitude)(int& flv, int& hel, int& col);
    double OMEGA_FUN(color_sum)(int& flv, int& hel);
    // void __parameters_sm_MOD_init_parameters(void);
    void OMEGA_FUN(init)(double*, int& scheme);
}

void new_event(double* data)
{
    OMEGA_FUN(new_event)(data);
}

void init_parameters(std::vector<double> pars)
{
    int foo = 0;
    return OMEGA_FUN(init)(pars.data(), foo);
}


std::complex<double> get_amplitude(int flv, int hel, int col)
{
    return OMEGA_FUN(get_amplitude)(flv, hel, col);
}

double color_sum(int flv, int hel)
{
    return OMEGA_FUN(color_sum)(flv, hel);
}

// TODO: make this take a beam polarisation
double sqme_sum(int flv)
{
    double res = 0.;
    for (int i = 1; i <= 64; i++) {
        res += color_sum(flv, i);
    }
    return res;
}

double sqme(int flv)
{
    // to average over input helicities instead of summing
    return sqme_sum(flv) / 4.0;
}

// FIXME: not very efficient because the amplitude will be recalculate like this :(
// need to find a better interface to match with RDataFrame...
double sqme(double* p, int flv)
{
    new_event(p);
    return sqme(flv);
}

// FIXME: urgh maybe make this a functor struct instead...
// or make this even properly object oriented?
// where each instance of the "sqme calculator" owns a set of parameters and takes care of all the set up before
// but will I be able to feed that into RDF from python?
// should work if I can get the operator() to have the correct signature
// object might be the correct thing, it can also keep track of the number of helicities etc.
double sqme_flv_sum(double* p, std::vector<double> pars)
{
    init_parameters(pars);
    new_event(p);
    return sqme(1) + sqme(2);
}

std::complex<double> hel_sum(int flv, int col)
{
    std::complex<double> res;
    for (int i = 1; i <= 64; i++) {
        res += get_amplitude(flv, i, col);
    }
    return res;
}