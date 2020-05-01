#pragma once
#include <iostream>
#include <vector>

void wwcFit(const std::vector<double>&x, 
         const std::vector<double>&y,
         double (*) (double, const std::vector<double>& params),
         std::vector<double>& params,
         int max_time = 1000, bool any = true, const char* write_to_file = "./fit.dat", const char* psr_name = "j0218");

