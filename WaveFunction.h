#pragma once
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <limits>


class WaveFunction
{
public:
	WaveFunction(int n = 101)noexcept
		: m_X(n)
		, m_RealPsi(n)
		, m_ImaginaryPsi(n)
		, m_RealHamiltonian(n)
		, m_ImaginaryHamiltonian(n)
		, m_Densities(n)
		, m_N(n)
	{}

	void Simulation(const float omega, std::string const& fileName, const float kappa = 1.f, const int n = 1, const float  tau = 500000, const float deltaTau = 0.0001)noexcept
	{
		float maxEnergy = std::numeric_limits<float>::min();
		std::ofstream normsFiles(fileName);
		normsFiles << "Nr:\t N:\t X:\t E:\n";
		ResetX(kappa);
		ResetRealPart(n, kappa, omega);
		ResetImaginaryPart();

		const int numberOfSteps = tau * deltaTau;
		for (float ii = 0; ii < numberOfSteps; ii += deltaTau)
		{
			UpdateRealPsi(deltaTau);
			UpdateHamiltonian(m_RealHamiltonian, m_RealPsi, kappa, omega, ii);

			UpdateImaginaryPsi(deltaTau);
			UpdateHamiltonian(m_ImaginaryHamiltonian, m_ImaginaryPsi, kappa, omega, ii);

			UpdateRealPsi(deltaTau);

			if ((int)ii % 100 == 0)
			{
				UpdateDensities();
				const float energy = CalcEnergy();
				normsFiles << ii * deltaTau << "\t" << CalcN() << "\t" << CalcX() << "\t" << energy << "\n";
				if (maxEnergy < energy) maxEnergy = energy;
			}
		}
		printf("max Energy :%f\n", maxEnergy);

	}

protected:

	inline void ResetX(const float kappa)noexcept
	{
		const float tempA = kappa / m_N;
		for (int ii = 0; ii < m_N; ++ii)
		{
			m_X[ii] = tempA * ii;
		}
	}

	inline void ResetRealPart(const int n, const float kappa, const float omega)noexcept
	{
		const float sqrtOfTwo = std::sqrt(2);
		const float tempA = 1 * M_PI;
		for (int ii = 0; ii < m_N; ++ii)
		{
			m_RealPsi[ii] = sqrtOfTwo * std::sin(tempA * m_X[ii]);
		}
		UpdateHamiltonian(m_RealHamiltonian, m_RealPsi, kappa, omega, 0);
	}
	inline void ResetImaginaryPart()noexcept
	{
		std::fill(m_ImaginaryPsi.begin(), m_ImaginaryPsi.end(), 0);
		std::fill(m_ImaginaryHamiltonian.begin(), m_ImaginaryHamiltonian.end(), 0);
	}
	inline void UpdateHamiltonian(std::vector<double>& hamiltonian, std::vector<double>const& psi, const float kappa, const float omega, const float tau)const noexcept
	{
		const float deltaX = 1.0f / (m_N - 1);
		const float tempA = -0.5 / (deltaX * deltaX);
		const float tempB = kappa * std::sin(omega * tau);

		hamiltonian[0] = 0;
		for (int ii = 1; ii < m_N - 1; ++ii)
		{
			hamiltonian[ii] = tempA * (psi[ii + 1] + psi[ii - 1] - 2 * psi[ii])
				+ tempB * (m_X[ii] - 0.5) * psi[ii];
		}
		hamiltonian[m_N - 1] = 0;
	}

	inline void UpdateRealPsi(const float deltaTau)noexcept
	{
		const float halfDeltaTau = deltaTau * 0.5;
		for (int ii = 0; ii < m_N; ++ii)
		{
			m_RealPsi[ii] += m_ImaginaryHamiltonian[ii] * halfDeltaTau;
		}
	}

	inline void UpdateImaginaryPsi(const float deltaTau)noexcept
	{
		for (int ii = 0; ii < m_N; ++ii)
		{
			m_ImaginaryPsi[ii] -= m_RealHamiltonian[ii] * deltaTau;
		}
	}

	inline void UpdateDensities() noexcept
	{
		for (int ii = 0; ii < m_N; ++ii)
		{
			m_Densities[ii] = m_RealPsi[ii] * m_RealPsi[ii] + m_ImaginaryPsi[ii] * m_ImaginaryPsi[ii];
		}
	}

	inline float CalcN()const noexcept
	{
		return (1.0f / (m_N - 1)) * std::accumulate(m_Densities.begin(), m_Densities.end(), 0.f);
	}

	inline float CalcX()const noexcept
	{
		return (1.0f / (m_N - 1)) * std::inner_product(m_X.begin(), m_X.end(), m_Densities.begin(), 0.f);
	}

	inline float CalcEnergy()const noexcept
	{
		float sum = 0;
		sum += std::inner_product(m_RealPsi.begin(), m_RealPsi.end(), m_RealHamiltonian.begin(), 0.f);
		sum += std::inner_product(m_ImaginaryPsi.begin(), m_ImaginaryHamiltonian.end(), m_ImaginaryPsi.begin(), 0.f);
		return (1.0f / (m_N - 1)) * sum;
	}

protected:
	std::vector<double> m_X;

	std::vector<double> m_RealPsi;
	std::vector<double> m_RealHamiltonian;

	std::vector<double> m_ImaginaryPsi;
	std::vector<double> m_ImaginaryHamiltonian;

	std::vector<double> m_Densities;

	int m_N;

};