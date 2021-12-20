#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>
#include <vector>
#include <ctime>
#include <cmath>
const double TwoPi = 6.283185307179586;

template <class T>
class isignal
{
    private:
        std::int16_t id;
        vector <T>  data;
        vector <T> spectrum;

    public:

        isignal(int id_s)
        {
            id = id_s;
        }
        ~isignal()
        {

        }

        void save_vector (vector<T> sign)
        {
            data = sign; //optimize copy procedure
        }   

        void print_vector ()
        {
            for (vector<int>::iterator it = data.begin(); it != data.end(); ++it)
            std::cout<<*it<<std::endl;
        }

        void FFTAnalysis ()
        {
            int i, j, n, m, Mmax, Istp;
            double Tmpr, Tmpi, Wtmp, Theta;
            double Wpr, Wpi, Wr, Wi;
            double *Tmvl;

            n = data.size() * 2; 
            Tmvl = new double[n];

            for (i = 0; i < n; i+=2) {
                Tmvl[i] = 0;
                Tmvl[i+1] = data[i/2];
            }

            i = 1; j = 1;
            while (i < n) {
                if (j > i) {
                Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
                Tmpr = Tmvl[i+1]; 
                Tmvl[i+1] = Tmvl[j+1]; 
                Tmvl[j+1] = Tmpr;
                }
                i = i + 2; m = data.size();
                while ((m >= 2) && (j > m)) {
                j = j - m; m = m >> 1;
                }
                j = j + m;
            }

             Mmax = 2;
                while (n > Mmax) {
                Theta = -TwoPi / Mmax; 
                Wpi = sin(Theta);
                Wtmp = sin(Theta / 2); 
                Wpr = Wtmp * Wtmp * 2;
                Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

                    while (m < Mmax) {
                        i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
                        Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
                        Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

                        while (i < n) {
                            j = i + Mmax;
                            Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
                            Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

                             Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
                             Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
                             i = i + Istp;
                        }
                    }

                    Mmax = Istp;
                }
            
            for (i = 0; i <data.size(); i++) {
                j = i * 2; 
                spectrum.push_back(2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/data.size());
            }

            delete [] Tmvl;

        }

        void print_spectrum ()
        {
             for (vector<int>::iterator it = spectrum.begin(); it != spectrum.end(); ++it)
                std::cout<<*it<<std::endl;
        }

        

};


#endif