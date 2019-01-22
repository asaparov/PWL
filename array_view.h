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

	inline const T& operator[] (size_t index) const {
		return array[index];
	}
};

template<typename T>
array_view<T> make_array_view(T* array, unsigned int length) {
	return array_view<T>(array, length);
}

template<typename T>
inline unsigned int size(const array_view<T>& view) {
	return view.length;
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

	inline const T& operator[] (size_t index) const {
		return array[indices[index]];
	}
};

template<typename T>
indexed_array_view<T> make_indexed_array_view(T* array, unsigned int* indices, unsigned int length) {
	return indexed_array_view<T>(array, indices, length);
}

template<typename T>
inline unsigned int size(const indexed_array_view<T>& view) {
	return view.length;
}

template<typename T, template<typename> class Array>
struct composed_array_view {
	T& first;
	Array<T>& second;

	composed_array_view(T& first, Array<T>& second) : first(first), second(second) { }

	inline T& operator[] (size_t index) {
		if (index == 0) return first;
		else return second[index - 1];
	}

	inline const T& operator[] (size_t index) const {
		if (index == 0) return first;
		else return second[index - 1];
	}
};

template<typename T, template<typename> class Array>
composed_array_view<T, Array> make_composed_array_view(T& first, Array<T>& second) {
	return composed_array_view<T, Array>(first, second);
}

template<typename T, template<typename> class Array>
inline unsigned int size(const composed_array_view<T, Array>& view) {
	return 1 + size(view.second);
}

#endif /* ARRAY_VIEW_H_ */
