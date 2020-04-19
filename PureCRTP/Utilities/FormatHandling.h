#pragma once
#include <typeinfo> 
#include <Eigen/Dense>

// FormatHandling
//
// Converts between interleaved audio and Eigen::ArrayXXf type.
//
// Example with 8 channel interleaved audio, where we only want 160 samples from the first 2 channels:
// 
// std::unique_ptr<FormantHandlingBase> FormatHandler;
// if ( x.Format == '16') // pseudo code that checks whether audio is 16 or 32 bit
//		{ FormatHandler = std::make_unique<FormatHandling<int16_t>(); }
// else { FormatHandler = std::make_unique<FormatHandling<int32_t>(); }
// int channels = 2;
// int samples = 160;
// Eigen::ArrayXXf XBuffer(samples,channels);
// FormatHandling.ReadIn(8, x, XBuffer); // number of samples/channels is deduced from XBuffer size
//
// author: Kristian Timm Andersen

class FormatHandlingBase
{
public:
	virtual ~FormatHandlingBase() = default;

	virtual void ReadIn(const int stride, const void* input, Eigen::Ref<Eigen::ArrayXXf> buffer) = 0;
	virtual void ReadIn(const int stride, const void* input, int index, Eigen::Ref<Eigen::ArrayXXf> buffer) = 0;
	virtual void ReadIn(const void* input, Eigen::Ref<Eigen::ArrayXXf> buffer) = 0;
	virtual void ReadIn(const void* input, int index, Eigen::Ref<Eigen::ArrayXXf> buffer) = 0;
	virtual void WriteOut(const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output) = 0;
	virtual void WriteOut(const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output, int index) = 0;
	virtual void WriteOut(const int stride, const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output) = 0;
	virtual void WriteOut(const int stride, const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output, int index) = 0;
};

template<typename T>
class FormatHandling final : public FormatHandlingBase
{
public:
	FormatHandling()
	{
		if (typeid(T) == typeid(int32_t))
		{
			BittoFloat = 4.6566129e-10f;			// 2^-31
			FloattoBit = 2147483648.f;				// 2^31
		}
		else if (typeid(T) == typeid(int16_t))
		{
			BittoFloat = 3.0517578125e-05f;			// 2^-15
			FloattoBit = 32768.f;					// 2^15
		}
		else
		{
			BittoFloat = 1.f;
			FloattoBit = 1.f;
		}
	}
	
	void ReadIn(const void* input, Eigen::Ref<Eigen::ArrayXXf> buffer) override
	{
		ReadIn( static_cast<int>(buffer.cols()), input, 0, buffer);
	}

	void ReadIn(const int stride, const void* input, Eigen::Ref<Eigen::ArrayXXf> buffer) override
	{
		ReadIn(stride, input, 0, buffer);
	}

	void ReadIn(const void* input, const int index, Eigen::Ref<Eigen::ArrayXXf> buffer) override
	{
		ReadIn(static_cast<int>(buffer.cols()), input, index, buffer);
	}

	void ReadIn(const int stride, const void* input, const int index, Eigen::Ref<Eigen::ArrayXXf> buffer) override
	{
		const Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::OuterStride<>> mapInput(static_cast<const T*>(input) + index * buffer.cols(), buffer.cols(), buffer.rows(), Eigen::OuterStride<>(stride));
		buffer = mapInput.transpose().template cast<float>() * BittoFloat;
	}

	void WriteOut(const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output) override
	{
		WriteOut( static_cast<int>(buffer.cols()), buffer, output, 0);
	}

	void WriteOut(const int stride, const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output) override
	{
		WriteOut(stride, buffer, output, 0);
	}

	void WriteOut(const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output, const int index) override
	{
		WriteOut( static_cast<int>(buffer.cols()), buffer, output, 0);
	}

	void WriteOut(const int stride, const Eigen::Ref<const Eigen::ArrayXXf>& buffer, void* output, const int index) override
	{
		Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, 0, Eigen::OuterStride<>> mapOutput(static_cast<T*>(output) + index * buffer.cols(), buffer.cols(), buffer.rows(), Eigen::OuterStride<>(stride));
		mapOutput = (buffer * FloattoBit).transpose().cast<T>();
	}
private:
	float BittoFloat, FloattoBit;
};

// if implementation is moved to .cpp file then the templates have to be put at the end of .cpp files, for instance:
//template class FormatHandling<int32_t>;
//template class FormatHandling<int16_t>;
//template class FormatHandling<float>;