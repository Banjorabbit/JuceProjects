#pragma once
#include "../BaseClasses/PureCRTP.h"

namespace IIR2ndOrderBaseClass
{
	// define general IIR2ndOrder class that other classes can be specialized from
	template<typename Tinput = I::Real2D, typename Toutput = O::Real2D, typename Tpersistent = I::Real>
	class IIR2ndOrder : public AsynchronousBase<IIR2ndOrder<Tinput, Toutput, Tpersistent>, Tinput, Toutput, Tpersistent>
	{
		friend Base<IIR2ndOrder<Tinput, Toutput, Tpersistent>, Tinput, Toutput, Tpersistent>;

	public:
		struct Coefficients
		{
			int NChannelsIn = 2;
			int BufferSize = 128;
			float SampleRate = 16e3f;
			AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
		} C;

		struct Parameters
		{
			ConstrainedType<float> Frequency = { 1000.f, 0.f, 20000.f };
			ConstrainedType<float> Q = { static_cast<float>(BUTTERWORTH_Q), 0.f, 1000.f };
			ConstrainedType<float> GaindB = { 0.f, -100.f, 100.f };
			enum FilterTypes { Flat, Peaking, Lowpass, Highpass, HighPassVariableQ, LowShelving, HighShelving, Notch, Bandpass, Custom };
			FilterTypes FilterType = Flat;
		} P;

		struct Data
		{
			float A1, A2, B0, B1, B2, Gain;
			typename I::GetType<Tinput>::type State1, State2; // set equal to input type
			void Reset() { State1.setZero(); State2.setZero(); }
			bool InitializeMemory(const Coefficients& c)
			{
				A1 = 0.f; A2 = 0.f; B0 = 1.f; A1 = 0.f; A1 = 0.f;
				State1.resize(c.NChannelsIn, 1);
				State2.resize(c.NChannelsIn, 1);
				return true;
			}
			size_t GetAllocatedMemorySize() const { return State1.GetAllocatedMemorySize() + State2.GetAllocatedMemorySize(); }

			void OnParameterChange(const Parameters& p, const Coefficients& c)
			{
				double Beta, dB0, dB1, dB2, dA0, dA1, dA2, t0, Omega, CosOmega, Alpha, K, Qinv, KK, GainKK, Denominator, SqrtGain, G;

				Gain = p.GaindB == 0 ? 1.f : powf(10.f, p.GaindB*0.05f);

				switch (p.FilterType)
				{
				case Parameters::Flat:
					B0 = Gain;
					B1 = 0.0f;
					B2 = 0.0f;
					A1 = 0.0f;
					A2 = 0.0f;
					break;
				case Parameters::Peaking:
					if (Gain == 1.0f)
					{
						B0 = 1.0f;
						B1 = 0.0f;
						B2 = 0.0f;
						A1 = 0.0f;
						A2 = 0.0f;
						break;
					}
					t0 = 2.0 * PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate);

					if (Gain >= 1.0f)
					{
						Beta = t0 / (2.0 * static_cast<double>(p.Q));
					}
					else
					{
						Beta = t0 / (2.0 * static_cast<double>(Gain) * static_cast<double>(p.Q));
					}
					dA2 = -0.5 * (1.0 - Beta) / (1.0 + Beta);
					dA1 = (0.5 - dA2) * cos(t0);
					dB0 = (static_cast<double>(Gain) - 1.0) * (0.25 + 0.5 * dA2) + 0.5;
					dB1 = -dA1;
					dB2 = -dB0 + 0.5 - dA2;

					B0 = static_cast<float>(dB0 * 2.0);
					B1 = static_cast<float>(dB1 * 2.0);
					B2 = static_cast<float>(dB2 * 2.0);
					A1 = static_cast<float>(dA1 * -2.0);
					A2 = static_cast<float>(dA2 * -2.0);
					break;
				case Parameters::Lowpass:
					Omega = 2.0 * PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate);
					CosOmega = cos(Omega);
					Alpha = sin(Omega) / (2.0* BUTTERWORTH_Q);

					dB0 = (1.0 - CosOmega) / 2.0;
					dB1 = 1.0 - CosOmega;
					dB2 = dB0;
					dA0 = 1.0 + Alpha;
					dA1 = -2.0 * CosOmega;
					dA2 = 1.0 - Alpha;

					B0 = static_cast<float>(dB0 / dA0);
					B1 = static_cast<float>(dB1 / dA0);
					B2 = static_cast<float>(dB2 / dA0);
					A1 = static_cast<float>(dA1 / dA0);
					A2 = static_cast<float>(dA2 / dA0);
					break;
				case Parameters::Highpass:
					Omega = 2.0 * PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate);
					CosOmega = cos(Omega);
					Alpha = sin(Omega) / (2.0 * BUTTERWORTH_Q);

					dB0 = (1.0 + CosOmega) / 2.0;
					dB1 = -1.0 - CosOmega;
					dB2 = dB0;
					dA0 = 1.0 + Alpha;
					dA1 = -2.0 * CosOmega;
					dA2 = 1.0 - Alpha;

					B0 = static_cast<float>(dB0 / dA0);
					B1 = static_cast<float>(dB1 / dA0);
					B2 = static_cast<float>(dB2 / dA0);
					A1 = static_cast<float>(dA1 / dA0);
					A2 = static_cast<float>(dA2 / dA0);
					break;
				case Parameters::HighPassVariableQ:						// High pass filter with variable Q factor
					Omega = 2.0 * PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate);
					CosOmega = cos(Omega);
					Alpha = sin(Omega) / (2.0 * p.Q);

					dB0 = (1.0 + CosOmega) / 2.0;
					dB1 = -1.0 - CosOmega;
					dB2 = dB0;
					dA0 = 1.0 + Alpha;
					dA1 = -2.0 * CosOmega;
					dA2 = 1.0 - Alpha;

					B0 = static_cast<float>(dB0 / dA0);
					B1 = static_cast<float>(dB1 / dA0);
					B2 = static_cast<float>(dB2 / dA0);
					A1 = static_cast<float>(dA1 / dA0);
					A2 = static_cast<float>(dA2 / dA0);
					break;
				case Parameters::LowShelving:						// Low pass shelving filter
					K = tan((PI * static_cast<double>(p.Frequency)) / static_cast<double>(c.SampleRate));
					Qinv = 1.0 / static_cast<double>(p.Q);
					KK = K * K;

					if (Gain > 1.0f)
					{
						GainKK = static_cast<double>(Gain) * KK;
						Denominator = (1.0 + Qinv * K + KK);
						SqrtGain = sqrt(static_cast<double>(Gain));
						dB0 = (1.0 + SqrtGain * Qinv * K + GainKK) / Denominator;
						dB1 = (2.0 * (GainKK - 1.0)) / Denominator;
						dB2 = (1.0 - SqrtGain * Qinv * K + GainKK) / Denominator;
						dA1 = (2.0 * (KK - 1.0)) / Denominator;
						dA2 = (1.0 - Qinv * K + KK) / Denominator;
					}
					else if (Gain < 1.0f)
					{
						Gain = 1.0f / Gain;
						GainKK = static_cast<double>(Gain) * KK;
						SqrtGain = sqrt(static_cast<double>(Gain));
						Denominator = (1.0 + Qinv * SqrtGain * K + GainKK);
						dB0 = (1.0 + Qinv * K + KK) / Denominator;
						dB1 = (2.0 * (KK - 1.0)) / Denominator;
						dB2 = (1.0 - Qinv * K + KK) / Denominator;
						dA1 = (2.0 * (GainKK - 1.0)) / Denominator;
						dA2 = (1.0 - Qinv * SqrtGain * K + GainKK) / Denominator;
					}
					else
					{
						dB0 = 1.0;
						dB1 = 0.0;
						dB2 = 0.0;
						dA1 = 0.0;
						dA2 = 0.0;
					}
					B0 = static_cast<float>(dB0);						// Output from intermediary coefficients
					B1 = static_cast<float>(dB1);
					B2 = static_cast<float>(dB2);
					A1 = static_cast<float>(dA1);
					A2 = static_cast<float>(dA2);
					break;
				case Parameters::HighShelving:						// High pass shelving filter
					K = tan(PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate));
					Qinv = 1.0 / static_cast<double>(p.Q);
					KK = K * K;

					if (Gain > 1.0f)
					{
						Denominator = (1.0 + Qinv * K + KK);
						SqrtGain = sqrt(static_cast<double>(Gain));
						dB0 = (static_cast<double>(Gain) + Qinv * SqrtGain * K + KK) / Denominator;
						dB1 = (2.0 * (KK - static_cast<double>(Gain))) / Denominator;
						dB2 = (static_cast<double>(Gain) - Qinv * SqrtGain * K + KK) / Denominator;
						dA1 = (2.0 * (KK - 1.0)) / Denominator;
						dA2 = (1.0 - Qinv * K + KK) / Denominator;
					}
					else if (Gain < 1.0f)
					{
						Gain = 1.0f / Gain;
						SqrtGain = sqrt(static_cast<double>(Gain));
						Denominator = (static_cast<double>(Gain) + Qinv * SqrtGain * K + KK);
						dB0 = (1.0 + Qinv * K + KK) / Denominator;
						dB1 = (2.0 * (KK - 1.0)) / Denominator;
						dB2 = (1.0 - Qinv * K + KK) / Denominator;
						dA1 = (2.0 * ((KK) / static_cast<double>(Gain) - 1.0)) / (1.0 + Qinv / SqrtGain * K + (KK) / static_cast<double>(Gain));
						dA2 = (1.0 - Qinv / SqrtGain * K + (KK) / static_cast<double>(Gain)) / (1.0 + Qinv / SqrtGain * K + (KK) / static_cast<double>(Gain));
					}
					else
					{
						dB0 = 1.0;
						dB1 = 0.0;
						dB2 = 0.0;
						dA1 = 0.0;
						dA2 = 0.0;
					}
					B0 = static_cast<float>(dB0);						// Output from intermediary coefficients
					B1 = static_cast<float>(dB1);
					B2 = static_cast<float>(dB2);
					A1 = static_cast<float>(dA1);
					A2 = static_cast<float>(dA2);
					break;
				case Parameters::Notch:
					Omega = 2.0 * PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate);

					Beta = tan(Omega / (static_cast<double>(p.Q) + static_cast<double>(p.Q)));
					G = 1.0 / (1.0 + Beta);

					dB0 = G;
					dB1 = (-G - G) * cos(Omega);
					dB2 = G;
					dA1 = dB1;
					dA2 = G + G - 1.0;

					B0 = static_cast<float>(dB0);						// Output from intermediary coefficients
					B1 = static_cast<float>(dB1);
					B2 = static_cast<float>(dB2);
					A1 = static_cast<float>(dA1);
					A2 = static_cast<float>(dA2);
					break;
				case Parameters::Bandpass:
					Omega = 2.0 * PI * static_cast<double>(p.Frequency) / static_cast<double>(c.SampleRate);
					Alpha = sin(Omega) / (static_cast<double>(p.Q) + static_cast<double>(p.Q));
					dB0 = Alpha;
					dB1 = 0.0;
					dB2 = -Alpha;
					dA0 = 1.0 + Alpha;
					dA1 = -2.0 * cos(Omega);
					dA2 = 1.0 - Alpha;

					B0 = static_cast<float>(dB0 / dA0);
					B1 = static_cast<float>(dB1 / dA0);
					B2 = static_cast<float>(dB2 / dA0);
					A1 = static_cast<float>(dA1 / dA0);
					A2 = static_cast<float>(dA2 / dA0);
					break;
				default:
					break;
				}
			}
		} D;

		void ProcessOn(Input xTime, Output yTime) { yTime = xTime; }
		void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
	};
} // end namespace IIR2ndOrderBaseClass

// IIR 2nd order Direct-Form 2 filter with support for any input/output types
template<typename Tinput = I::Real2D, typename Toutput = O::Real2D>
class IIR2ndDF2AnyType : public AsynchronousBase<IIR2ndDF2AnyType<Tinput, Toutput>, Tinput, Toutput>
{
	friend Base<IIR2ndDF2AnyType<Tinput, Toutput>, Tinput, Toutput>;

	typedef typename IIR2ndOrderBaseClass::IIR2ndOrder<Tinput, Toutput>::Coefficients Coefficients;
	typedef typename IIR2ndOrderBaseClass::IIR2ndOrder<Tinput, Toutput>::Parameters Parameters;
	typedef typename IIR2ndOrderBaseClass::IIR2ndOrder<Tinput, Toutput>::Data Data;

	Coefficients C;
	Parameters P;
	Data D;

	void ProcessOn(Input xTime, Output yTime)
	{
		for (auto sample = 0; sample < C.BufferSize; sample++)
		{
			for (auto channel = 0; channel < C.NChannelsIn; channel++) // Channel in inner loop is faster according to profiling
			{
				I::GetScalarType<Tinput>::type sub = xTime(sample, channel) - D.A1 * D.State1(channel) - D.A2 * D.State2(channel);
				yTime(sample, channel) = D.B0 * sub + D.B1 * D.State1(channel) + D.B2 * D.State2(channel);

				// Update state variables
				D.State2(channel) = D.State1(channel);
				D.State1(channel) = sub;
			}

		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

// IIR 2nd order Direct-Form 2 Transposed filter with support for any input/output types
template<typename Tinput = I::Real2D, typename Toutput = O::Real2D>
class IIR2ndDF2TransposedAnyType : public AsynchronousBase<IIR2ndDF2TransposedAnyType<Tinput, Toutput>, Tinput, Toutput, I::Real>
{
	friend Base<IIR2ndDF2TransposedAnyType<Tinput, Toutput>, Tinput, Toutput, I::Real>;

	typedef typename IIR2ndOrderBaseClass::IIR2ndOrder<Tinput, Toutput>::Coefficients Coefficients;
	typedef typename IIR2ndOrderBaseClass::IIR2ndOrder<Tinput, Toutput>::Parameters Parameters;
	typedef typename IIR2ndOrderBaseClass::IIR2ndOrder<Tinput, Toutput>::Data Data;

	Coefficients C;
	Parameters P;
	Data D;

	void ProcessPersistentInput(InputPersistent filterCoefficients)
	{
		if (P.FilterType == P.Custom)
		{
			if (filterCoefficients(3) != 0.f)
			{
				D.B0 = filterCoefficients(0) / filterCoefficients(3);
				D.B1 = filterCoefficients(1) / filterCoefficients(3);
				D.B2 = filterCoefficients(2) / filterCoefficients(3);
				D.A1 = filterCoefficients(4) / filterCoefficients(3);
				D.A2 = filterCoefficients(5) / filterCoefficients(3);
			}
			else // this is not a valid filter, so set all coefficients to zero
			{
				D.B0 = 0.f;
				D.B1 = 0.f;
				D.B1 = 0.f;
				D.A1 = 0.f;
				D.A2 = 0.f;
			}
		}
	}

	void ProcessOn(Input xTime, Output yTime)
	{
		for (auto sample = 0; sample < C.BufferSize; sample++)
		{
			for (auto channel = 0; channel < C.NChannelsIn; channel++) // Channel in inner loop is faster according to profiling
			{
				O::GetScalarType<Toutput>::type out = xTime(sample, channel) * D.B0 + D.State1(channel); // can not write directly to yTime if yTime and xTime are same memory
				// Update state variables
				D.State1(channel) = xTime(sample, channel) * D.B1 - out * D.A1 + D.State2(channel);
				D.State2(channel) = xTime(sample, channel) * D.B2 - out * D.A2;
				yTime(sample, channel) = out;
			}
		}
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }
};

// Cascaded IIR 2nd order Direct-Form 2 Transposed filter with support for any input/output types
template<typename Tinput = I::Real2D, typename Toutput = O::Real2D>
class IIR2ndCascadedAnyType : public AsynchronousBase<IIR2ndCascadedAnyType<Tinput, Toutput>, Tinput, Toutput>
{
	friend Base<IIR2ndCascadedAnyType<Tinput, Toutput>, Tinput, Toutput>;

public:
	VectorAlgo<IIR2ndDF2TransposedAnyType<Tinput, Toutput>> Filters;

private:

	struct Coefficients 
	{
		int Nsos = 4;
		int BufferSize = 128;
		int NChannelsIn = 2;
		AsynchronousBufferType AsynchronousBuffer = VARIABLE_SIZE;
	} C;

	struct Parameters { float Gain = 1;  } P;

	struct Data 
	{
		void Reset() { }
		bool InitializeMemory(const Coefficients& c) { return true; }
		size_t GetAllocatedMemorySize() const { return 0; }
		void OnParameterChange(const Parameters& p, const Coefficients& c) {}
	} D;

	DEFINEMEMBERALGORITHMS(1, Filters);

	auto InitializeMembers()
	{
		Filters.resize(C.Nsos);
		auto s = Filters[0].GetSetup();
		s.Coefficients.NChannelsIn = C.NChannelsIn;
		s.Coefficients.BufferSize = C.BufferSize;
		s.Coefficients.SampleRate = 0.f; // not used since P.FilterType == P.Custom
		s.Parameters.FilterType = s.Parameters.Custom;
		return Filters.Initialize(s);
	}

	void ProcessPersistentInput(InputPersistent filterCoefficients)
	{
		for (auto i = 0; i < C.Nsos; i++)
		{
			Filters[i].SetPersistentInput(filterCoefficients.row(i).transpose());
		}
	}

	void ProcessOn(Input xTime, Output yTime) 
	{
		yTime = xTime * P.Gain;
		for (auto i = 0; i < C.Nsos; i++) { Filters[i].Process(yTime, yTime); }
	}

	void ProcessOff(Input xTime, Output yTime) { yTime = xTime; }

};

// Alias for specific types
using IIR2ndDF2 = IIR2ndDF2AnyType<I::Real2D, O::Real2D>;
using IIR2ndDF2Complex = IIR2ndDF2AnyType<I::Complex2D, O::Complex2D>;
using IIR2ndDF2Transposed = IIR2ndDF2TransposedAnyType<I::Real2D, O::Real2D>;
using IIR2ndDF2TransposedComplex = IIR2ndDF2TransposedAnyType<I::Complex2D, O::Complex2D>;
using IIR2ndCascaded = IIR2ndCascadedAnyType<I::Real2D, O::Real2D>;
using IIR2ndCascadedComplex = IIR2ndCascadedAnyType<I::Complex2D, O::Complex2D>;