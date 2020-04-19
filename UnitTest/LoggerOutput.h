// OutputLog prints to output log using the << operator
//
// author: Kristian Timm Andersen

#pragma once

class LoggerOutput {
public:
	template<class T>
	LoggerOutput& operator<<(const T& output) 
	{ 
		S << output; 
		return *this; 
	}
	LoggerOutput& operator<<(std::ostream& (*fn)(std::ostream&)) // this is used for printing
	{ 
		fn(S); 
		Print(); 
		S = std::stringstream();  
		return *this; 
	} 
private:
	std::stringstream S;
	void Print() const { 
		Microsoft::VisualStudio::CppUnitTestFramework::Logger::WriteMessage(S.str().c_str()); 
	}
};
