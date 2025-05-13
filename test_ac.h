#include <complex>
#include <vector>

extern "C" {
    void __cc20_ac_amplitude_MOD_new_event(double*);
    void __cc20_ac_amplitude_MOD_print_momenta(double*);
    double _Complex __cc20_ac_amplitude_MOD_get_amplitude(int& flv, int& hel, int& col);
    double __cc20_ac_amplitude_MOD_color_sum(int& flv, int& hel);
    // void __parameters_sm_MOD_init_parameters(void);
    void __cc20_ac_amplitude_MOD_init(double*, int& scheme);
}


void init_parameters(std::vector<double> pars)
{
    int foo = 0;
    return __cc20_ac_amplitude_MOD_init(pars.data(), foo);
}


std::complex<double> get_amplitude(int flv, int hel, int col)
{
    return __cc20_ac_amplitude_MOD_get_amplitude(flv, hel, col);
}

double color_sum(int flv, int hel)
{
    return __cc20_ac_amplitude_MOD_color_sum(flv, hel);
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