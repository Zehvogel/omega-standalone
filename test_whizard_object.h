#include <complex>
#include <vector>
#include <iostream>


#define OMEGA_FUN(name) __ ## opr_ww_i1 ## _MOD_ ## name


extern "C" {
    void OMEGA_FUN(new_event)(const double*);
    double OMEGA_FUN(color_sum)(int& flv, int& hel);
    void OMEGA_FUN(init)(const double*, int& scheme);
}

class OmegaWrapper {
    public:
    OmegaWrapper(const std::vector<double>& parameters) : m_parameters(parameters) {}

    double sqme(double* lvec, int flv)
    {
        init_parameters();
        new_event(lvec);

        double res = 0;
        for (int i = 1; i <= 64; i++) {
            res += OMEGA_FUN(color_sum)(flv, i);
        }
        //average over input helicity combinations
        return res / 4.0;
    }


    double sqme1(std::vector<double> lvec) {
        return sqme(lvec.data(), 1);
    }

    double sqme2(std::vector<double> lvec) {
        return sqme(lvec.data(), 2);
    }

    private:
    const std::vector<double> m_parameters;
    void init_parameters() {
        int scheme = 0;
        OMEGA_FUN(init)(m_parameters.data(), scheme);
    }

    void new_event(const double* data) {OMEGA_FUN(new_event)(data);}
};

struct OmegaWrapperFunctor {
    const int m_flv;
    OmegaWrapper m_omw;

    OmegaWrapperFunctor(const int flv, OmegaWrapper omw) : m_flv(flv), m_omw(omw) {}

    double operator()(std::vector<double> lvec) {
        return m_omw.sqme(lvec.data(), m_flv);
    } 
};
