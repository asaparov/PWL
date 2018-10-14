#ifndef ARRAY_VIEW_H_
#define ARRAY_VIEW_H_

template<typename T>
struct array_view {
	T* array;
	unsigned int length;

	array_view(T* array, unsigned int length) : array(array), length(length) { }

	inline T& operator[] (size_t index) {
		return array[index];
	}
};

template<typename T>
array_view<T> make_array_view(T* array, unsigned int length) {
	return array_view<T>(array, length);
}

template<typename T>
struct indexed_array_view {
	T* array;
	unsigned int* indices;
	unsigned int length;

	indexed_array_view(T* array, unsigned int* indices, unsigned int length) : array(array), indices(indices), length(length) { }

	inline T& operator[] (size_t index) {
		return array[indices[index]];
	}
};

template<typename T>
indexed_array_view<T> make_indexed_array_view(T* array, unsigned int* indices, unsigned int length) {
	return indexed_array_view<T>(array, indices, length);
}

#endif /* ARRAY_VIEW_H_ */
