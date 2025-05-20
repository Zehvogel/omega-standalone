#include <complex>
#include <vector>


#define AMP_PREFIX cc20_ac_amplitude
#define OMEGA_FUN(name) __ ## cc20_ac_amplitude ## _MOD_ ## name


extern "C" {
    void OMEGA_FUN(new_event)(double*);
    void OMEGA_FUN(print_momenta)(double*);
    double _Complex OMEGA_FUN(get_amplitude)(int& flv, int& hel, int& col);
    double OMEGA_FUN(color_sum)(int& flv, int& hel);
    // void __parameters_sm_MOD_init_parameters(void);
    void OMEGA_FUN(init)(double*, int& scheme);
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

std::complex<double> hel_sum(int flv, int col)
{
    std::complex<double> res;
    for (int i = 1; i <= 64; i++) {
        res += get_amplitude(flv, i, col);
    }
    return res;
}