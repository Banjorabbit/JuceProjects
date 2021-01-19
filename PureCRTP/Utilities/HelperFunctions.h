
// Assign values to array with an offset
//
// author: Kristian Timm Andersen
template<typename T, size_t m, size_t n>
void AssignArray(T(&array)[m], const T(&values)[n], const int offset = 0) { for (size_t i = 0; i < std::min(n, m - offset); i++) { array[i + offset] = values[i]; } }