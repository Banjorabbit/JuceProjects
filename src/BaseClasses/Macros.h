// author: Kristian Timm Andersen

#define EVAL(...) __VA_ARGS__
#define TESTMEMBER Member
#define TESTMEMBERS members

#define DEFINEMEMBER10(method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER9)(method, __VA_ARGS__))
#define DEFINEMEMBER9( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER8)(method, __VA_ARGS__))
#define DEFINEMEMBER8( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER7)(method, __VA_ARGS__))
#define DEFINEMEMBER7( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER6)(method, __VA_ARGS__))
#define DEFINEMEMBER6( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER5)(method, __VA_ARGS__))
#define DEFINEMEMBER5( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER4)(method, __VA_ARGS__))
#define DEFINEMEMBER4( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER3)(method, __VA_ARGS__))
#define DEFINEMEMBER3( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER2)(method, __VA_ARGS__))
#define DEFINEMEMBER2( method, member, ...) decltype(member.method()) member; EVAL(EVAL(DEFINEMEMBER1)(method, __VA_ARGS__))
#define DEFINEMEMBER1( method, member)      decltype(member.method()) member; 

#define SETMEMBERMETHOD10(method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD9)(method, __VA_ARGS__))
#define SETMEMBERMETHOD9( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD8)(method, __VA_ARGS__))
#define SETMEMBERMETHOD8( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD7)(method, __VA_ARGS__))
#define SETMEMBERMETHOD7( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD6)(method, __VA_ARGS__))
#define SETMEMBERMETHOD6( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD5)(method, __VA_ARGS__))
#define SETMEMBERMETHOD5( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD4)(method, __VA_ARGS__))
#define SETMEMBERMETHOD4( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD3)(method, __VA_ARGS__))
#define SETMEMBERMETHOD3( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD2)(method, __VA_ARGS__))
#define SETMEMBERMETHOD2( method, member, ...) member.method(TESTMEMBERS.member); EVAL(EVAL(SETMEMBERMETHOD1)(method, __VA_ARGS__))
#define SETMEMBERMETHOD1( method, member)      member.method(TESTMEMBERS.member); 

#define APPLYMEMBERMETHOD10(method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD9)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD9( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD8)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD8( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD7)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD7( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD6)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD6( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD5)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD5( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD4)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD4( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD3)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD3( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD2)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD2( method, member, ...) member.method EVAL(EVAL(APPLYMEMBERMETHOD1)( method, __VA_ARGS__))
#define APPLYMEMBERMETHOD1( method, member     ) member.method 

#define DEFINEMEMBERSETGETFUNCTIONS(N, type, ...)  \
	struct EVAL(TESTMEMBER)type { \
		type type; \
		EVAL(EVAL(DEFINEMEMBER##N)(EVAL(Get)type, __VA_ARGS__)) \
	}; \
	template<typename...Ts> \
	EVAL(TESTMEMBER)type EVAL(EVAL(Get)type)EVAL(TESTMEMBER)(const Ts&... ts) const \
	{ \
		return EVAL(EVAL(TESTMEMBER)type){ EVAL(Get)type(), ts.EVAL(Get)type()... }; \
	} \
	EVAL(TESTMEMBER)type EVAL(EVAL(Get)type)AllImpl() const { return EVAL(EVAL(Get)type)EVAL(TESTMEMBER)(__VA_ARGS__ ); } \
	void EVAL(EVAL(Set)type)AllImpl(const EVAL(TESTMEMBER)type& TESTMEMBERS) { \
		EVAL(Set)type(TESTMEMBERS.type); \
		EVAL(EVAL(SETMEMBERMETHOD##N)(EVAL(Set)type, __VA_ARGS__ )) \
	} 


#define DEFINEMEMBERALGORITHMS(N,...) \
	DEFINEMEMBERSETGETFUNCTIONS(N, Parameters, __VA_ARGS__) \
	DEFINEMEMBERSETGETFUNCTIONS(N, Coefficients, __VA_ARGS__) \
	DEFINEMEMBERSETGETFUNCTIONS(N, Setup, __VA_ARGS__) \
	void ResetMembers()      { EVAL(EVAL(APPLYMEMBERMETHOD##N)( EVAL(Reset();),      __VA_ARGS__ )) } \
	size_t GetMemberAllocatedMemorySize() const { return EVAL(EVAL(APPLYMEMBERMETHOD##N)(EVAL(GetAllocatedMemorySize() + ), __VA_ARGS__ ))0; }