#ifndef __M_KURAMOTO__
#define __M_KURAMOTO__

#include <vector>
#include <cmath>

// TODO: increase readability by also considering matrices
namespace nonlinear_systems{
  template<typename state_type = std::vector<double> >
    class MKuramotoODE {
      public:
        std::vector<double> omega, K, alpha, ensemble_size, group_number;

        std::vector<std::vector<double> > CalculateMeanField(
            const state_type& y) {
          std::vector<std::vector<double> > mean_field(ensemble_size.size());
          for (size_t i = 0; i < ensemble_size.size(); ++i) {
            mean_field[i] = CalculateMeanFieldGroup(i, y);
          }
          return mean_field;
        } 


        std::vector<double> CalculateMeanFieldGroup(int group_number, 
            const state_type& y) {
          double real_part = 0, imag_part = 0;
          for (unsigned int i = group_index[group_number]; 
              i < group_index[group_number+1]; ++i) {
            real_part += cos(y[i]);
            imag_part += sin(y[i]); 
          }
          double group_size = group_index[group_number+1] 
            - group_index[group_number];
          std::vector<double> mean_field(2);
          mean_field[0] = 1./group_size*sqrt(real_part*real_part 
              + imag_part*imag_part);
          mean_field[1] = atan2(imag_part, real_part);
          return mean();
        }


        void operator() (const state_type& y, state_type& dy, 
            const double dt) {
          std::vector<std::vector<double> > mean_field = 
            CalculateMeanField(y); 
          for (int i = 0; i < ensemble_size.size(); ++i) {
            for (int j = group_index[i]; j < group_index[i+1]; ++j) {
              double sum = 0;
              for (int k = 0; k < ensemble_size.size(); ++k) {
                int coupling_index = i*ensemble_size.size()+k;
                sum += 
                  coupling[coupling_index]*ensemble_size[k]/N
                  *mean_ensemble[k][0]*sin(mean_ensemble[k][1] - x[j] + 
                      phase_shift[coupling_index]);
              }
              dx[j] = frequency[j] + sum;
            }
          }

        }
    };
} // nonlinear_systems
#endif
