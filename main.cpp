#pragma once
#define _USE_MATH_DEFINES
#include <chrono>
#include <array>
#include <string>

#include "WaveFunction.h"


int main()
{
	const auto t1 = std::chrono::high_resolution_clock::now();
	{
		constexpr float omega = 0.5 * 3 * M_PI*M_PI;
		constexpr std::array<float, 10> omegaMultipliers = { 0.9, 0.92, 0.94, 0.98, 1, 1.02, 1.04, 1.06, 1.08, 1.1 };
		for (int ii=0; ii< omegaMultipliers.size(); ++ii)
		{
			printf("omega id: %d\n", ii);
			WaveFunction wave;
			wave.Simulation(omega*omegaMultipliers[ii], "results/omega" + std::to_string(ii) + ".txt");
		}
	}
	const auto t2 = std::chrono::high_resolution_clock::now();
	printf("took %lld  ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());

	return 0;
}