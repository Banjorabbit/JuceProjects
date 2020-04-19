#pragma once
#include <tuple>
//#include <utility>
//#include <string>

//template<typename... Tmembers>
//class NamedTuple
//{
//public:
//	NamedTuple(std::string tnames...) { Names = std::make_tuple(tnames...); }
//private:
//	std::tuple<Tmembers...> M;
//	std::tuple<std::string...> Names;
//};

template<typename... Tmembers>
class AlgorithmTuple
{
public:

	void Initialize()
	{
		auto c = GetCoefficients();
		Initialize(c);
	}

	template<std::size_t I = 0, typename... Tcs>
	typename std::enable_if<I == sizeof...(Tmembers), void>::type Initialize(const std::tuple<Tcs...>& c) { }

	template<std::size_t I = 0, typename... Tcs>
	typename std::enable_if < I < sizeof...(Tmembers), void>::type Initialize(const std::tuple<Tcs...>& c)
	{
		Initialize(std::get<I>(c), std::get<I>(M));
		Initialize<I + 1>(c);
	}

	template<std::size_t I = 0>
	typename std::enable_if<I == sizeof...(Tmembers), void>::type Reset() { }

	template<std::size_t I = 0>
	typename std::enable_if < I < sizeof...(Tmembers), void>::type Reset()
	{
		Reset(std::get<I>(M));
		Reset<I + 1>();
	}

	template<std::size_t I = 0>
	typename std::enable_if<I == sizeof...(Tmembers), size_t>::type	GetAllocatedMemorySize() const { return 0; }

	template<std::size_t I = 0>
	typename std::enable_if < I < sizeof...(Tmembers), size_t>::type GetAllocatedMemorySize() const
	{
		return GetAllocatedMemorySize(std::get<I>(M)) + GetAllocatedMemorySize<I + 1>();
	}

	template<int I> auto& Get() { return std::get<I>(M); } // this is used to modify returned element so it must return reference (c++14 decltype(auto) should also do that)

	auto GetCoefficients() const { return GetCoefficientsImpl(std::index_sequence_for<Tmembers...>{}); }

	template<std::size_t I = 0, typename... Tcs>
	typename std::enable_if<I == sizeof...(Tmembers), void>::type SetCoefficients(const std::tuple<Tcs...>& c) {}

	template<std::size_t I = 0, typename... Tcs>
	typename std::enable_if < I < sizeof...(Tmembers), void>::type SetCoefficients(const std::tuple<Tcs...>& c)
	{
		SetCoefficients(std::get<I>(c), std::get<I>(M));
		SetCoefficients<I + 1>(c);
	}

	auto GetParameters() const { return GetParametersImpl(std::index_sequence_for<Tmembers...>{}); }

	template<std::size_t I = 0, typename... Tps>
	typename std::enable_if<I == sizeof...(Tmembers), void>::type SetParameters(const std::tuple<Tps...>& p) {}

	template<std::size_t I = 0, typename... Tps>
	typename std::enable_if<I < sizeof...(Tmembers), void>::type SetParameters(const std::tuple<Tps...>& p) 
	{
		SetParameters(std::get<I>(p), std::get<I>(M));
		SetParameters<I + 1>(p);
	}

	auto GetCoefficientsAll() const { return GetCoefficientsAllImpl(std::index_sequence_for<Tmembers...>{}); }
	auto GetParametersAll() const { return GetParametersAllImpl(std::index_sequence_for<Tmembers...>{}); }
	auto GetSetupAll() const { return GetSetupAllImpl(std::index_sequence_for<Tmembers...>{}); }

	//template<std::size_t I = 0, typename Tcoefficient, typename... Tcoefficients>
	//void Initialize(const std::tuple<std::vector<Tcoefficient>, Tcoefficients...>& c)
	//{
	//	for (auto i = 0; i < std::get<0>(c).size(); i++)
	//	{
	//		std::get<I>(M)[i].Initialize(std::get<0>(c)[i]);
	//	}
	//	Initialize<I + 1>(Tail(c));
	//}

	//template<std::size_t I = 0, typename Tcoefficient, typename... Tcoefficients>
	//void Initialize(const std::tuple<Tcoefficient, Tcoefficients...>& c)
	//{
	//	std::get<I>(M).Initialize(std::get<0>(c));
	//	Initialize<I + 1>(Tail(c));
	//}

	//template<std::size_t I = 0, typename... Tcoefficients>
	//void Initialize(const std::tuple<Tcoefficients...>& c) { }

	//template<std::size_t I = 0>	typename std::enable_if < I < sizeof...(Tmembers), void>::type
	//	Reset() const
	//{
	//	for (auto i = 0; i < std::get<I>(M).size(); i++)
	//	{
	//		std::get<I>(M)[i].Reset();
	//	}
	//	Reset<I + 1>();
	//}

	//template<std::size_t I = 0>	typename std::enable_if<I == sizeof...(Tmembers), void>::type
	//Reset() const { }

	//template<std::size_t I = 0> typename std::enable_if<I == sizeof...(Tmembers), void>::type
	//	Initialize(const std::tuple<typename Tmembers::Coefficients...>& c) { }

	//template<std::size_t I = 0> typename std::enable_if<I < sizeof...(Tmembers), void>::type
	//	Initialize(const std::tuple<typename Tmembers::Coefficients...>& c) { std::get<I>(M).Initialize(std::get<I>(c)); Initialize<I + 1>(c); }

	//template<std::size_t I = 0>	typename std::enable_if<I == sizeof...(Tmembers), void>::type
	//	Reset() const { }

	//template<std::size_t I = 0>	typename std::enable_if<I < sizeof...(Tmembers), void>::type
	//	Reset() const { std::get<I>(M).Reset(); Reset<I + 1>(); }

	//template<std::size_t I = 0>	typename std::enable_if<I == sizeof...(Tmembers), size_t>::type
	//	GetAllocatedMemorySize() const { return 0; }

	//template<std::size_t I = 0>	typename std::enable_if < I < sizeof...(Tmembers), size_t>::type
	//	GetAllocatedMemorySize() const { return (std::get<I>(M)).GetAllocatedMemorySize() + GetAllocatedMemorySize<I + 1>(); }

	/*template<typename... Tcoefficients>
	void Initialize(const std::tuple<Tcoefficients...>& c) { Initialize(c, M); }*/

	//template<typename... Tcs>
	//void Initialize(const std::tuple<Tcs...>&c) { Initialize(c, M); }

	//void Reset() const { Reset(M); }

	//size_t GetAllocatedMemorySize() const { return GetAllocatedMemorySize(M); }

	//decltype(auto) GetCoefficients() const { return GetCoefficients(M); }

	//decltype(auto) GetParameters() const { return GetParameters(M); }

private:
	std::tuple<Tmembers...> M;

	template<typename Tc, typename Tm>
	void Initialize(const Tc& c, Tm& m) { m.Initialize(c); }

	template<typename Tc, typename Tm>
	void Initialize(const std::vector<Tc>& c, std::vector<Tm>& m)
	{
		for (auto i = 0; i < m.size(); i++) { m[i].Initialize(c[i]); }
	}

	template<typename Tm>
	void Reset(Tm& m) { m.Reset(); }

	template<typename Tm>
	void Reset(std::vector<Tm>& m) { for (auto& member : m) { member.Reset(); } }

	template<typename Tm>
	size_t GetAllocatedMemorySize(const Tm& m) const { return m.GetAllocatedMemorySize(); }

	template<typename Tm>
	size_t GetAllocatedMemorySize(const std::vector<Tm>& m) const
	{
		size_t size = 0;
		for (auto& member : m)
		{
			size += member.GetAllocatedMemorySize() + sizeof(member); // sizeof(member) is not counted as part of static memory since it's inside a vector 
		}
		return size;
	}

	template<size_t... Is>
	auto GetCoefficientsImpl(const std::index_sequence<Is...>) const
	{
		return std::make_tuple(GetCoefficients(std::get<Is>(M))...);
	}

	template<typename Tm>
	auto GetCoefficients(const Tm& m) const { return m.GetCoefficients(); }

	template<typename Tm>
	auto GetCoefficients(const std::vector<Tm>& m) const
	{
		std::vector<decltype(m[0].GetCoefficients())> cVec(m.size());
		for (auto i = 0; i < m.size(); i++) { cVec[i] = m[i].GetCoefficients(); }
		return cVec;
	}

	template<size_t... Is>
	auto GetParametersImpl(const std::index_sequence<Is...>) const
	{
		return std::make_tuple(GetParameters(std::get<Is>(M))...);
	}

	template<typename Tm>
	auto GetParameters(const Tm& m) const { return m.GetParameters(); }

	template<typename Tm>
	auto GetParameters(const std::vector<Tm>& m) const
	{
		std::vector<decltype(m[0].GetParameters())> pVec(m.size());
		for (auto i = 0; i < m.size(); i++) { pVec[i] = m[i].GetParameters(); }
		return pVec;
	}

	template<typename Tp, typename Tm>
	void SetParameters(const Tp& p, Tm& m) { m.SetParameters(p); }

	template<typename Tp, typename Tm>
	void SetParameters(const std::vector<Tp>& p, std::vector<Tm>& m)
	{
		for (auto i = 0; i < m.size(); i++) { m[i].SetParameters(p[i]); }
	}
	
	template<size_t... Is>
	auto GetCoefficientsAllImpl(const std::index_sequence<Is...>) const
	{
		return std::make_tuple(GetCoefficientsAll(std::get<Is>(M))...);
	}

	template<typename Tm>
	auto GetCoefficientsAll(const Tm& m) const { return m.GetCoefficientsAll(); }

	template<typename Tm>
	auto GetCoefficientsAll(const std::vector<Tm>& m) const
	{
		std::vector<decltype(m[0].GetCoefficientsAll())> cVec(m.size());
		for (auto i = 0; i < m.size(); i++) { cVec[i] = m[i].GetCoefficientsAll(); }
		return cVec;
	}

	template<size_t... Is>
	auto GetSetupAllImpl(const std::index_sequence<Is...>) const
	{
		return std::make_tuple(GetSetupAll(std::get<Is>(M))...);
	}

	template<typename Tm>
	auto GetSetupAll(const Tm& m) const { return m.GetSetupAll(); }

	template<typename Tm>
	auto GetSetupAll(const std::vector<Tm>& m) const
	{
		std::vector<decltype(m[0].GetSetupAll())> cVec(m.size());
		for (auto i = 0; i < m.size(); i++) { cVec[i] = m[i].GetSetupAll(); }
		return cVec;
	}

	template<typename Tc, typename Tm>
	void SetCoefficients(const Tc& c, Tm& m) { m.SetCoefficients(c); }

	template<typename Tc, typename Tm>
	void SetCoefficients(const std::vector<Tc>& c, std::vector<Tm>& m)
	{
		for (auto i = 0; i < m.size(); i++) { m[i].SetCoefficients(c[i]); }
	}

	template<size_t... Is>
	auto GetParametersAllImpl(const std::index_sequence<Is...>) const
	{
		return std::make_tuple(GetParametersAll(std::get<Is>(M))...);
	}

	template<typename Tm>
	auto GetParametersAll(const Tm& m) const { return m.GetParametersAll(); }

	template<typename Tm>
	auto GetParametersAll(const std::vector<Tm>& m) const
	{
		std::vector<decltype(m[0].GetParametersAll())> pVec(m.size());
		for (auto i = 0; i < m.size(); i++) { pVec[i] = m[i].GetParametersAll(); }
		return pVec;
	}

	//template<typename... Tcs, typename... Tms>
	//void Initialize(const std::tuple<Tcs...>& c, std::tuple<Tms...>& m) 
	//{
	//	int dummy = 1;
	//}

	//template<typename Tc, typename Tm, typename... Tcs, typename... Tms>
	//void Initialize(const std::tuple<Tc, Tcs...>& c, std::tuple<Tm, Tms...>& m) 
	//{
	//	std::get<0>(m).Initialize(std::get<0>(c));
	//	Initialize(Tail(c), Tail(m));
	//}

	//template<typename Tc, typename Tm, typename... Tcs, typename... Tms>
	//void Initialize(const std::tuple<std::vector<Tc>, Tcs...>& c, std::tuple<std::vector<Tm>, Tms...>& m)
	//{
	//	for (auto i = 0; i < std::get<0>(c).size(); i++)
	//	{
	//		std::get<0>(m)[i].Initialize(std::get<0>(c)[i]);
	//	}
	//	Initialize(Tail(c), Tail(m));
	//}
	

	//template<typename... Tms>
	//auto GetParameters(std::tuple<Tms...> m) const { return std::make_tuple(); }

	//template<typename Tm, typename... Tms>
	//auto GetParameters(std::tuple<Tm, Tms...> m) const
	//{
	//	return std::tuple_cat(std::make_tuple(std::get<0>(m).GetParameters()), GetParameters(Tail(m)));
	//}

	//template<typename Tm, typename... Tms>
	//auto GetParameters(std::tuple<std::vector<Tm>, Tms...> m) const
	//{
	//	auto memberVec = std::get<0>(m);
	//	std::vector<decltype(memberVec[0].GetParameters())> cVec(memberVec.size());
	//	for (auto i = 0; i < memberVec.size(); i++) { cVec[i] = memberVec[i].GetParameters(); }
	//	return std::tuple_cat(std::make_tuple(cVec), GetParameters(Tail(m)));
	//}

	//template<typename... Tms>
	//auto GetCoefficients(std::tuple<Tms...> m) const { return std::make_tuple(); }

	//template<typename Tm, typename... Tms>
	//auto GetCoefficients(std::tuple<Tm, Tms...> m) const
	//{ 
	//	return std::tuple_cat(std::make_tuple(std::get<0>(m).GetCoefficients()), GetCoefficients(Tail(m))); 
	//}

	//template<typename Tm, typename... Tms>
	//auto GetCoefficients(std::tuple<std::vector<Tm>, Tms...> m) const
	//{
	//	auto memberVec = std::get<0>(m);
	//	std::vector<decltype(memberVec[0].GetCoefficients())> cVec(memberVec.size());
	//	for (auto i = 0; i < memberVec.size(); i++) { cVec[i] = memberVec[i].GetCoefficients(); }
	//	return std::tuple_cat(std::make_tuple(cVec), GetCoefficients(Tail(m))); 
	//}

	//template<typename... Tms>
	//void Reset(std::tuple<Tms...> m) const {}

	//template<typename Tm, typename... Tms>
	//void Reset(std::tuple<Tm, Tms...> m) const
	//{
	//	std::get<0>(m).Reset();
	//	Reset(Tail(m));
	//}

	//template<typename Tm, typename... Tms>
	//void Reset(std::tuple<std::vector<Tm>, Tms...> m) const
	//{
	//	auto memberVec = std::get<0>(m);
	//	for (auto& member : memberVec) { member.Reset(); }
	//	Reset(Tail(m));
	//}

	//template<typename... Tms>
	//size_t GetAllocatedMemorySize(std::tuple<Tms...> m) const { return 0; }

	//template<typename Tm, typename... Tms>
	//size_t GetAllocatedMemorySize(std::tuple<Tm, Tms...> m) const 
	//{
	//	return std::get<0>(m).GetAllocatedMemorySize() + GetAllocatedMemorySize(Tail(m));
	//}

	//template<typename Tm, typename... Tms>
	//size_t GetAllocatedMemorySize(std::tuple<std::vector<Tm>, Tms...> m) const
	//{
	//	size_t size = 0;
	//	auto memberVec = std::get<0>(m);
	//	for (auto& member : memberVec) { size += member.GetAllocatedMemorySize() + sizeof(member); } // sizeof(member) is not counted as part of static memory since it's inside of a vector
	//	return size + GetAllocatedMemorySize(Tail(m));
	//}

	//template <class Tm, class... Tms >
	//std::tuple<Tms...> Tail(const std::tuple<Tm, Tms...>& t) const
	//{
	//	return  TailImpl(std::index_sequence_for<Tms...>(), t);
	//}

	//template <std::size_t... Ns, typename... Ts >
	//auto TailImpl(const std::index_sequence<Ns...>, std::tuple<Ts...> t) const
	//{
	//	return  std::make_tuple(std::get<Ns + 1u>(t)...);
	//}

};
